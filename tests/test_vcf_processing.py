"""Tests for VCF processing helpers."""

from __future__ import annotations

import os

from vartracker.vcf_processing import _derive_vcf_output_paths


def test_derive_vcf_output_paths_handles_gz(tmp_path):
    out, csq, log = _derive_vcf_output_paths(
        "/data/sample1.vcf.gz", tmp_path, "sample1"
    )

    assert out.endswith("sample1.cyvcf2.vcf.gz")
    assert csq.endswith("sample1.csq.vcf.gz")
    assert log.endswith("sample1.log")
    assert os.path.dirname(out) == str(tmp_path)


def test_derive_vcf_output_paths_handles_plain_vcf(tmp_path):
    out, csq, log = _derive_vcf_output_paths("/data/alpha.vcf", tmp_path, "alpha")

    assert out.endswith("alpha.cyvcf2.vcf.gz")
    assert csq.endswith("alpha.csq.vcf.gz")
    assert log.endswith("alpha.log")
    assert os.path.dirname(out) == str(tmp_path)
