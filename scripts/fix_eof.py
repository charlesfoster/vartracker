#!/usr/bin/env python3
"""Ensure text files end with a single newline."""

from __future__ import annotations

import sys
from pathlib import Path


def fix_file(path: Path) -> bool:
    if path.is_dir():
        return False
    try:
        data = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return False  # skip binary files
    if not data:
        return False
    if data.endswith("\n") and not data.endswith("\n\n"):
        return False
    # Trim excess blank lines at end and ensure single newline
    stripped = data.rstrip("\n") + "\n"
    path.write_text(stripped, encoding="utf-8")
    return True


def main(argv: list[str]) -> int:
    changed = False
    for arg in argv:
        if fix_file(Path(arg)):
            print(f"Fixed end-of-file newline: {arg}")
            changed = True
    return 1 if changed else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
