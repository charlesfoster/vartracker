"""Run provenance capture for vartracker."""

from __future__ import annotations

import datetime as dt
import json
import os
import platform
import shutil
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from ._version import __version__
from .core import FILE_COLUMNS
from .schemas import RESULTS_SCHEMA_VERSION

MANIFEST_SCHEMA_VERSION = "1.0"

_TOOL_VERSION_COMMANDS: dict[str, list[list[str]]] = {
    "snakemake": [["snakemake", "--version"], ["snakemake", "-v"]],
    "bcftools": [["bcftools", "--version"], ["bcftools", "-v"]],
    "samtools": [["samtools", "--version"], ["samtools", "version"]],
    "tabix": [["tabix", "--version"], ["tabix", "-V"]],
    "bgzip": [["bgzip", "--version"], ["bgzip", "-V"]],
    "bwa": [["bwa"], ["bwa", "version"]],
    "fastp": [["fastp", "--version"], ["fastp", "-V"]],
    "lofreq": [["lofreq", "--version"], ["lofreq", "-V"]],
}

_DEFAULT_TOOLS = tuple(_TOOL_VERSION_COMMANDS.keys())


def _utc_now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def _first_nonempty_line(text: str) -> str | None:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return None


def compute_sha256(path: Path) -> tuple[str | None, str | None]:
    """Compute sha256 for a file; return (digest, error_message)."""
    try:
        import hashlib

        hasher = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                hasher.update(chunk)
        return hasher.hexdigest(), None
    except Exception as exc:  # pragma: no cover - non-critical metadata
        return None, str(exc)


def _file_entry(path: Path, *, include_hash: bool) -> dict[str, Any]:
    record: dict[str, Any] = {"path": str(path)}
    try:
        record["size_bytes"] = path.stat().st_size
    except Exception as exc:  # pragma: no cover - filesystem edge cases
        record["size_bytes_error"] = str(exc)
    if include_hash:
        digest, error = compute_sha256(path)
        record["sha256"] = digest
        if error:
            record["sha256_error"] = error
    return record


def _resolve_csv_path(value: Any, csv_dir: Path) -> Path | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = (csv_dir / path).resolve()
    else:
        path = path.resolve()
    return path


def collect_input_files(
    input_rows: Iterable[Mapping[str, Any]],
    *,
    csv_path: Path,
    include_hashes: bool,
) -> list[dict[str, Any]]:
    """Collect referenced input files with optional hashes."""
    csv_dir = csv_path.parent.resolve()
    seen: set[Path] = set()
    records: list[dict[str, Any]] = []

    for row in input_rows:
        for column in FILE_COLUMNS:
            if column not in row:
                continue
            path = _resolve_csv_path(row[column], csv_dir)
            if path is None:
                continue
            if path in seen:
                continue
            seen.add(path)
            record = _file_entry(path, include_hash=include_hashes)
            record["column"] = column
            records.append(record)

    return records


def _run_version_command(cmd: Sequence[str]) -> str | None:
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None

    output = result.stdout or ""
    if not output.strip():
        output = result.stderr or ""
    return _first_nonempty_line(output)


def get_tool_version(tool: str) -> dict[str, Any]:
    path = shutil.which(tool)
    record: dict[str, Any] = {"path": path, "version": None}
    if path is None:
        return record

    for cmd in _TOOL_VERSION_COMMANDS.get(tool, [[tool, "--version"], [tool, "-V"]]):
        line = _run_version_command(cmd)
        if line:
            record["version"] = line
            record["command"] = " ".join(cmd)
            break

    return record


def collect_tool_versions(tools: Iterable[str] = _DEFAULT_TOOLS) -> dict[str, Any]:
    return {tool: get_tool_version(tool) for tool in tools}


def get_git_info(repo_dir: Path) -> dict[str, Any]:
    if shutil.which("git") is None:
        return {"available": False}

    try:
        inside = subprocess.run(
            ["git", "rev-parse", "--is-inside-work-tree"],
            cwd=repo_dir,
            capture_output=True,
            text=True,
            check=False,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return {"available": False}

    if inside.returncode != 0 or inside.stdout.strip().lower() != "true":
        return {"available": False}

    def _git(cmd: Sequence[str]) -> str | None:
        try:
            result = subprocess.run(
                cmd,
                cwd=repo_dir,
                capture_output=True,
                text=True,
                check=False,
                timeout=5,
            )
        except (OSError, subprocess.TimeoutExpired):
            return None
        return _first_nonempty_line(result.stdout or result.stderr or "")

    commit = _git(["git", "rev-parse", "HEAD"])
    describe = _git(["git", "describe", "--tags", "--dirty", "--always"])
    status = _git(["git", "status", "--porcelain"])

    return {
        "available": True,
        "commit": commit,
        "describe": describe,
        "dirty": bool(status),
    }


class RunManifest:
    """Manage writing of run_metadata.json provenance manifests."""

    def __init__(
        self,
        *,
        outdir: str | Path,
        command: str | None,
        cli_args: Sequence[str] | None,
        invocation: str | None,
        manifest_level: str = "light",
    ) -> None:
        self.outdir = Path(outdir).expanduser().resolve()
        self.run_id = str(uuid.uuid4())
        self.started_at = _utc_now_iso()
        self.manifest_level = manifest_level
        self.payload: dict[str, Any] = {
            "manifest_schema_version": MANIFEST_SCHEMA_VERSION,
            "run_id": self.run_id,
            "status": "running",
            "started_at": self.started_at,
            "command": command,
            "cli_invocation": invocation,
            "cli_args": list(cli_args or []),
            "outdir": str(self.outdir),
            "manifest_level": manifest_level,
            "vartracker_version": __version__,
            "results_schema_version": RESULTS_SCHEMA_VERSION,
            "python": {
                "version": platform.python_version(),
                "executable": sys.executable,
            },
            "platform": {
                "system": platform.system(),
                "release": platform.release(),
                "machine": platform.machine(),
                "platform": platform.platform(),
            },
            "git": get_git_info(Path(__file__).resolve().parents[1]),
        }

        container_image = os.environ.get("VARTRACKER_CONTAINER_IMAGE")
        container_digest = os.environ.get("VARTRACKER_CONTAINER_DIGEST")
        if container_image or container_digest:
            self.payload["container"] = {
                "image": container_image,
                "digest": container_digest,
            }

    def start(self) -> None:
        self._write_payload()

    def update(self, **fields: Any) -> None:
        self._deep_update(self.payload, fields)
        self._write_payload()

    def finish(
        self,
        *,
        status: str,
        outputs: Mapping[str, Any] | None = None,
        error: str | None = None,
        extra: Mapping[str, Any] | None = None,
    ) -> None:
        ended_at = _utc_now_iso()
        self.payload["status"] = status
        self.payload["ended_at"] = ended_at
        try:
            started = dt.datetime.fromisoformat(self.started_at)
            ended = dt.datetime.fromisoformat(ended_at)
            self.payload["duration_seconds"] = (ended - started).total_seconds()
        except ValueError:
            self.payload["duration_seconds"] = None

        if outputs is not None:
            merged_outputs = dict(self.payload.get("outputs", {}))
            merged_outputs.update(outputs)
            self.payload["outputs"] = merged_outputs
        if error:
            self.payload["error"] = error
        if extra:
            self._deep_update(self.payload, extra)

        self._write_payload()

    def _write_payload(self) -> None:
        try:
            self.outdir.mkdir(parents=True, exist_ok=True)
            runs_dir = self.outdir / "runs"
            runs_dir.mkdir(parents=True, exist_ok=True)
            latest_path = self.outdir / "run_metadata.json"
            run_path = runs_dir / f"run_metadata.{self.run_id}.json"
            payload_json = json.dumps(self.payload, indent=2, sort_keys=True)
            latest_path.write_text(payload_json, encoding="utf-8")
            run_path.write_text(payload_json, encoding="utf-8")
        except Exception as exc:  # pragma: no cover - avoid breaking runs
            print(f"Warning: failed to write run manifest: {exc}")

    @staticmethod
    def _deep_update(target: dict[str, Any], updates: Mapping[str, Any]) -> None:
        for key, value in updates.items():
            if (
                key in target
                and isinstance(target[key], dict)
                and isinstance(value, Mapping)
            ):
                RunManifest._deep_update(target[key], value)
            else:
                target[key] = value


__all__ = [
    "MANIFEST_SCHEMA_VERSION",
    "collect_input_files",
    "collect_tool_versions",
    "compute_sha256",
    "RunManifest",
]
