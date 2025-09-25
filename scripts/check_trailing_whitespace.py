#!/usr/bin/env python3
"""Remove trailing whitespace from text files."""

from __future__ import annotations

import sys
from pathlib import Path

TRAILING = " \t"


def strip_trailing_whitespace(path: Path) -> bool:
    if path.is_dir():
        return False
    try:
        original = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return False  # skip binary files
    lines = original.splitlines(keepends=True)
    stripped_lines = []
    changed = False
    for line in lines:
        if line.endswith("\n"):
            core = line[:-1].rstrip(TRAILING)
            new_line = core + "\n"
        else:
            new_line = line.rstrip(TRAILING)
        if new_line != line:
            changed = True
        stripped_lines.append(new_line)
    if changed:
        path.write_text("".join(stripped_lines), encoding="utf-8")
    return changed


def main(argv: list[str]) -> int:
    changed = False
    for arg in argv:
        if strip_trailing_whitespace(Path(arg)):
            print(f"Trimmed trailing whitespace: {arg}")
            changed = True
    return 1 if changed else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
