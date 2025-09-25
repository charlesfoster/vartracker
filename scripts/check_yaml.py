#!/usr/bin/env python3
"""Basic YAML validator for pre-commit."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from typing import Any, Protocol, cast


class _YamlModule(Protocol):
    class YAMLError(Exception): ...

    def safe_load(self, stream: Any) -> Any: ...


def _import_yaml_module() -> _YamlModule:
    try:
        module = importlib.import_module("yaml")
    except ModuleNotFoundError as exc:  # pragma: no cover - dependency missing
        raise SystemExit("PyYAML is required to run check_yaml") from exc
    return cast(_YamlModule, module)


yaml = _import_yaml_module()
YamlError = getattr(yaml, "YAMLError", Exception)


def validate_yaml(path: Path) -> None:
    if path.is_dir():
        return
    with path.open("r", encoding="utf-8") as fh:
        yaml.safe_load(fh)


def main(argv: list[str]) -> int:
    exit_code = 0
    for arg in argv:
        path = Path(arg)
        try:
            validate_yaml(path)
        except YamlError as exc:
            print(f"YAML error in {path}: {exc}", file=sys.stderr)
            exit_code = 1
        except Exception as exc:  # pragma: no cover - unexpected issues
            print(f"Failed to validate {path}: {exc}", file=sys.stderr)
            exit_code = 1
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
