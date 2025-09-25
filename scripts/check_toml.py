#!/usr/bin/env python3
"""Basic TOML validator for pre-commit."""

from __future__ import annotations

import importlib
import importlib.util
import sys
from pathlib import Path
from typing import Any, Protocol, cast


class _TomlModule(Protocol):
    class TOMLDecodeError(Exception): ...

    def load(self, fp: Any) -> Any: ...


def _import_toml_module() -> _TomlModule:
    if importlib.util.find_spec("tomllib") is not None:  # Python 3.11+
        module = importlib.import_module("tomllib")
    else:  # pragma: no cover - compatibility with older interpreters
        module = importlib.import_module("tomli")
    return cast(_TomlModule, module)


tomllib = _import_toml_module()


def validate_toml(path: Path) -> None:
    if path.is_dir():
        return
    with path.open("rb") as fh:
        tomllib.load(fh)


def main(argv: list[str]) -> int:
    exit_code = 0
    for arg in argv:
        path = Path(arg)
        try:
            validate_toml(path)
        except tomllib.TOMLDecodeError as exc:
            print(f"TOML error in {path}: {exc}", file=sys.stderr)
            exit_code = 1
        except Exception as exc:  # pragma: no cover - unexpected issues
            print(f"Failed to validate {path}: {exc}", file=sys.stderr)
            exit_code = 1
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
