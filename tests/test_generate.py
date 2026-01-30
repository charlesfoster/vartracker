"""Tests for the CSV generation helpers."""

from __future__ import annotations

from vartracker.generate import (
    generate_csv,
    generate_from_bams,
    generate_from_reads,
    generate_from_vcfs,
)


def test_generate_from_reads_parses_numbers(tmp_path):
    r1 = tmp_path / "sample_5_R1.fastq.gz"
    r2 = tmp_path / "sample_5_R2.fastq.gz"
    r1.write_text("", encoding="utf-8")
    r2.write_text("", encoding="utf-8")

    records = generate_from_reads(tmp_path)
    assert len(records) == 1
    record = records[0]
    assert record.sample_name == "sample_5"
    assert record.sample_number == 5
    assert record.reads1.endswith("sample_5_R1.fastq.gz")
    assert record.reads2.endswith("sample_5_R2.fastq.gz")


def test_generate_from_bams_requires_bam(tmp_path):
    bam = tmp_path / "test_psg_3.bam"
    bam.write_text("", encoding="utf-8")

    records = generate_from_bams(tmp_path)
    assert len(records) == 1
    record = records[0]
    assert record.sample_name == "test_psg_3"
    assert record.sample_number == 3
    assert record.bam.endswith("test_psg_3.bam")


def test_generate_from_vcfs_handles_suffixes(tmp_path):
    vcf = tmp_path / "sample_8.lofreq.vcf.gz"
    vcf.write_text("", encoding="utf-8")

    records = generate_from_vcfs(tmp_path)
    assert len(records) == 1
    record = records[0]
    assert record.sample_name == "sample_8"
    assert record.sample_number == 8
    assert record.vcf.endswith("sample_8.lofreq.vcf.gz")


def test_generate_csv_writes_output(tmp_path):
    vcf = tmp_path / "sample_1.vcf"
    vcf.write_text("", encoding="utf-8")
    output = tmp_path / "out.csv"

    records, warnings = generate_csv("vcf", tmp_path, output, dry_run=False)

    assert output.exists()
    assert isinstance(records, list)
    assert isinstance(warnings, list)
    data = output.read_text(encoding="utf-8")
    assert "sample_1" in data


def test_generate_csv_dry_run_returns_records(tmp_path):
    vcf = tmp_path / "sample_2.vcf"
    vcf.write_text("", encoding="utf-8")
    records, warnings = generate_csv(
        "vcf", tmp_path, tmp_path / "out.csv", dry_run=True
    )

    assert len(records) == 1
    assert records[0].sample_name == "sample_2"
    assert warnings == []
