"""Version information for the vartracker package."""

from __future__ import annotations

from pathlib import Path

try:  # Python 3.8+
    from importlib import metadata
except ImportError:  # pragma: no cover -- fallback for very old Pythons
    import importlib_metadata as metadata  # type: ignore

_LAST_RESORT_VERSION = "0.0.0+unknown"


def _read_pyproject_version() -> str | None:
    """Read `project.version` straight from pyproject.toml.

    Only reachable when vartracker isn't pip-installed (no dist-info to read
    metadata from), e.g. a git checkout run in place - so pyproject.toml is
    still sitting one directory above this file. Keeps the version fallback
    in sync with pyproject.toml automatically, rather than a second
    hand-maintained copy of the version string.
    """
    pyproject_path = Path(__file__).resolve().parent.parent / "pyproject.toml"
    if not pyproject_path.exists():
        return None
    try:
        import tomllib
    except ImportError:  # pragma: no cover -- Python <3.11
        return None
    with pyproject_path.open("rb") as handle:
        data = tomllib.load(handle)
    return data.get("project", {}).get("version")


try:
    __version__ = metadata.version("vartracker")
except metadata.PackageNotFoundError:  # pragma: no cover - during local dev
    __version__ = _read_pyproject_version() or _LAST_RESORT_VERSION

__all__ = ["__version__"]
