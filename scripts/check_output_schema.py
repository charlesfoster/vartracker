#!/usr/bin/env python3
"""Check that docs/OUTPUT_SCHEMA.md matches the canonical schema."""

from __future__ import annotations

import difflib
import importlib.util
from pathlib import Path

SCHEMA_PATH = Path("vartracker") / "schemas.py"
DOC_PATH = Path("docs") / "OUTPUT_SCHEMA.md"


def load_schema_module():
    spec = importlib.util.spec_from_file_location("schemas", SCHEMA_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    if not DOC_PATH.exists():
        print(f"Missing {DOC_PATH}")
        return 1

    module = load_schema_module()
    expected = module.render_output_schema_markdown().strip()
    actual = DOC_PATH.read_text(encoding="utf-8").strip()

    if expected == actual:
        print("Output schema documentation is up to date.")
        return 0

    diff = "\n".join(
        difflib.unified_diff(
            actual.splitlines(),
            expected.splitlines(),
            fromfile=str(DOC_PATH),
            tofile="generated",
            lineterm="",
        )
    )
    print("Output schema documentation is out of date.\n")
    print(diff)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
