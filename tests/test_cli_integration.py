"""Integration tests for the vartracker CLI."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys

import pytest


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not available")
@pytest.mark.skipif(shutil.which("tabix") is None, reason="tabix not available")
def test_cli_test_flag_runs_successfully(tmp_path):
    outdir = tmp_path / "results"
    mpl_dir = tmp_path / "mpl"
    mpl_dir.mkdir()

    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", str(mpl_dir))
    if not env.get("VARTRACKER_REQUIRE_BCFTOOLS"):
        env["VARTRACKER_SKIP_BCFTOOLS"] = "1"

    cmd = [
        sys.executable,
        "-m",
        "vartracker.main",
        "vcf",
        "--test",
        "--outdir",
        str(outdir),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert result.returncode == 0, result.stderr

    results_csv = outdir / "results.csv"
    assert results_csv.exists()

    content = results_csv.read_text(encoding="utf-8")
    assert "variant" in content
