"""Tests for VCF processing helpers."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest
from cyvcf2 import VCF

from vartracker.analysis import process_joint_variants
from vartracker.vcf_processing import (
    _derive_vcf_output_paths,
    annotate_vcf,
    format_vcf,
    merge_consequences,
    process_vcf,
)


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


def _write_depth_file(path: Path, *, seqid: str = "chr1", length: int = 9) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for pos in range(1, length + 1):
            handle.write(f"{seqid}\t{pos}\t100\n")


def _write_minimal_reference_bundle(tmp_path: Path) -> tuple[Path, Path]:
    ref = tmp_path / "ref.fa"
    gff = tmp_path / "ref.gff3"
    ref.write_text(">chr1\nATGAAATAA\n", encoding="utf-8")
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tvartracker_test\tgene\t1\t9\t.\t+\t.\tID=gene:chr1_1;Name=GENE1;biotype=protein_coding\n"
        "chr1\tvartracker_test\tmRNA\t1\t9\t.\t+\t.\tID=transcript:chr1_1;Parent=gene:chr1_1;Name=GENE1;biotype=protein_coding\n"
        "chr1\tvartracker_test\tCDS\t1\t9\t.\t+\t0\tID=cds:chr1_1;Parent=transcript:chr1_1;gene=GENE1\n",
        encoding="utf-8",
    )
    return ref, gff


def _prepare_merged_multiallelic_vcf(
    tmp_path: Path, record_line: str, sample_name: str = "s1"
) -> tuple[Path, Path, Path]:
    ref, gff = _write_minimal_reference_bundle(tmp_path)
    sample_vcf = tmp_path / f"{sample_name}.vcf"
    sample_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=9>\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">\n'
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
        f"{record_line}\n",
        encoding="utf-8",
    )

    formatted_vcf, _ = format_vcf(
        str(sample_vcf),
        sample_name,
        str(tmp_path),
        0.0,
        0.0,
        str(ref),
        str(gff),
        False,
    )

    vcf_list = tmp_path / "vcf_list.txt"
    vcf_list.write_text(f"{formatted_vcf}\n", encoding="utf-8")
    sample_names = tmp_path / "sample_names.txt"
    sample_names.write_text(f"{sample_name}\n", encoding="utf-8")

    merged_vcf = tmp_path / "merged.vcf"
    merge_consequences(str(tmp_path), str(merged_vcf), str(sample_names), debug=False)
    return ref, gff, merged_vcf


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not available")
def test_format_vcf_preserves_distinct_alt_records_at_same_position(tmp_path):
    ref, gff = _write_minimal_reference_bundle(tmp_path)
    vcf_path = tmp_path / "same_pos_two_alts.vcf"
    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=9>\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">\n'
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr1\t5\t.\tA\tC\t.\tPASS\tAF=0.20;DP=100\tGT:DP:AF\t1:100:0.20\n"
        "chr1\t5\t.\tA\tG\t.\tPASS\tAF=0.30;DP=100\tGT:DP:AF\t1:100:0.30\n",
        encoding="utf-8",
    )

    formatted_vcf, _ = format_vcf(
        str(vcf_path),
        "s1",
        str(tmp_path),
        0.0,
        0.0,
        str(ref),
        str(gff),
        False,
    )

    records = list(VCF(formatted_vcf))

    assert len(records) == 2
    assert {record.ALT[0] for record in records} == {"C", "G"}


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not available")
def test_multiallelic_input_survives_merge_and_csq_annotation(tmp_path):
    ref, gff = _write_minimal_reference_bundle(tmp_path)
    sample_vcf = tmp_path / "multiallelic.vcf"
    sample_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=9>\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">\n'
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr1\t5\t.\tA\tC,G\t.\tPASS\tAF=0.20,0.30;DP=100\tGT:DP:AF\t1/2:100:0.20,0.30\n",
        encoding="utf-8",
    )

    formatted_vcf, _ = format_vcf(
        str(sample_vcf),
        "s1",
        str(tmp_path),
        0.0,
        0.0,
        str(ref),
        str(gff),
        False,
    )

    vcf_list = tmp_path / "vcf_list.txt"
    vcf_list.write_text(f"{formatted_vcf}\n", encoding="utf-8")
    sample_names = tmp_path / "sample_names.txt"
    sample_names.write_text("s1\n", encoding="utf-8")

    merged_vcf = tmp_path / "merged.vcf"
    merge_consequences(str(tmp_path), str(merged_vcf), str(sample_names), debug=False)

    annotated_vcf = tmp_path / "annotated.vcf"
    annotate_vcf(
        str(merged_vcf),
        str(annotated_vcf),
        str(ref),
        str(gff),
        debug=False,
    )

    depth = tmp_path / "s1.depth.txt"
    _write_depth_file(depth)

    table = process_vcf(str(annotated_vcf), [str(depth)], 10, ["s1"])

    observed = (
        table[["variant", "amino_acid_consequence", "bcsq_aa_notation", "alt_freq"]]
        .sort_values("variant")
        .reset_index(drop=True)
    )

    expected = pd.DataFrame(
        [
            {
                "variant": "A5C",
                "amino_acid_consequence": "K2T",
                "bcsq_aa_notation": "2K>2T",
                "alt_freq": "0.200",
            },
            {
                "variant": "A5G",
                "amino_acid_consequence": "K2R",
                "bcsq_aa_notation": "2K>2R",
                "alt_freq": "0.300",
            },
        ]
    )

    pd.testing.assert_frame_equal(observed, expected)


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not available")
def test_annotate_vcf_reports_informative_error_for_multiallelic_overflow(tmp_path):
    ref, gff, merged_vcf = _prepare_merged_multiallelic_vcf(
        tmp_path,
        "chr1\t5\t.\tA\tG,T,C\t.\tPASS\tAF=0.04,0.34,0.12;DP=100\tGT:DP:AF\t1/2:100:0.04,0.34,0.12",
    )

    with pytest.raises(RuntimeError) as excinfo:
        annotate_vcf(
            str(merged_vcf),
            str(tmp_path / "annotated.vcf"),
            str(ref),
            str(gff),
            debug=False,
            multiallelic_overflow="error",
        )

    message = str(excinfo.value)
    assert "chr1:5" in message
    assert "Sample s1:" in message
    assert "A>G AF=0.040" in message
    assert "A>T AF=0.340" in message
    assert "A>C AF=0.120" in message
    assert "--min-snv-freq above 0.040" in message


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not available")
def test_annotate_vcf_drop_lowest_af_continues_with_warning(tmp_path, capsys):
    ref, gff, merged_vcf = _prepare_merged_multiallelic_vcf(
        tmp_path,
        "chr1\t5\t.\tA\tG,T,C\t.\tPASS\tAF=0.04,0.34,0.12;DP=100\tGT:DP:AF\t1/2:100:0.04,0.34,0.12",
    )

    annotated_vcf = tmp_path / "annotated.vcf"
    annotate_vcf(
        str(merged_vcf),
        str(annotated_vcf),
        str(ref),
        str(gff),
        debug=False,
        multiallelic_overflow="drop-lowest-af",
    )

    stdout = capsys.readouterr().out
    assert "Warning:" in stdout
    assert "Dropping the lowest-frequency ALT allele(s) for this sample" in stdout
    assert "A>G AF=0.040" in stdout
    assert "--min-snv-freq above 0.040" in stdout

    depth = tmp_path / "s1.depth.txt"
    _write_depth_file(depth)
    table = process_vcf(str(annotated_vcf), [str(depth)], 10, ["s1"])

    observed = (
        table[["variant", "amino_acid_consequence", "alt_freq"]]
        .sort_values("variant")
        .reset_index(drop=True)
    )

    expected = pd.DataFrame(
        [
            {
                "variant": "A5C",
                "amino_acid_consequence": "K2T",
                "alt_freq": "0.120",
            },
            {
                "variant": "A5T",
                "amino_acid_consequence": "K2I",
                "alt_freq": "0.340",
            },
        ]
    )

    pd.testing.assert_frame_equal(observed, expected)


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not available")
def test_annotate_vcf_skip_site_keeps_unannotated_rows(tmp_path, capsys):
    ref, gff, merged_vcf = _prepare_merged_multiallelic_vcf(
        tmp_path,
        "chr1\t5\t.\tA\tG,T,C\t.\tPASS\tAF=0.04,0.34,0.12;DP=100\tGT:DP:AF\t1/2:100:0.04,0.34,0.12",
    )

    annotated_vcf = tmp_path / "annotated.vcf"
    annotate_vcf(
        str(merged_vcf),
        str(annotated_vcf),
        str(ref),
        str(gff),
        debug=False,
        multiallelic_overflow="skip-site",
    )

    stdout = capsys.readouterr().out
    assert "Warning:" in stdout
    assert "Skipping consequence calling for this site" in stdout
    assert "--min-snv-freq above 0.040" in stdout

    depth = tmp_path / "s1.depth.txt"
    _write_depth_file(depth)
    table = process_vcf(str(annotated_vcf), [str(depth)], 10, ["s1"])

    observed = (
        table[["variant", "amino_acid_consequence", "alt_freq"]]
        .sort_values("variant")
        .reset_index(drop=True)
    )

    expected = pd.DataFrame(
        [
            {
                "variant": "A5C",
                "amino_acid_consequence": "None",
                "alt_freq": "0.120",
            },
            {
                "variant": "A5G",
                "amino_acid_consequence": "None",
                "alt_freq": "0.040",
            },
            {
                "variant": "A5T",
                "amino_acid_consequence": "None",
                "alt_freq": "0.340",
            },
        ]
    )

    pd.testing.assert_frame_equal(observed, expected)


def test_process_vcf_splits_sample_specific_bcsq_annotations(tmp_path):
    vcf_path = tmp_path / "annotated.vcf"
    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=9>\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##INFO=<ID=BCSQ,Number=.,Type=String,Description="Consequence">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">\n'
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        '##FORMAT=<ID=BCSQ,Number=.,Type=Integer,Description="Bitmask of indexes to INFO/BCSQ">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        "chr1\t4\t.\tA\tG\t.\tPASS\tDP=100;AF=0.1;BCSQ=missense|GENE1|tx|protein_coding|+|2K>2A|4A>G+5A>C\tGT:DP:AF:BCSQ\t1:100:0.1:1\t0:.:.:0\n"
        "chr1\t5\t.\tA\tC\t.\tPASS\tDP=200;AF=0.1;BCSQ=@4,missense|GENE1|tx|protein_coding|+|2K>2T|5A>C\tGT:DP:AF:BCSQ\t1:100:0.1:1\t1:100:0.1:4\n",
        encoding="utf-8",
    )

    cov1 = tmp_path / "s1.depth.txt"
    cov2 = tmp_path / "s2.depth.txt"
    _write_depth_file(cov1)
    _write_depth_file(cov2)

    table = process_vcf(str(vcf_path), [str(cov1), str(cov2)], 10, ["s1", "s2"])
    observed = (
        table[
            [
                "start",
                "amino_acid_consequence",
                "bcsq_aa_notation",
                "presence_absence",
                "alt_freq",
            ]
        ]
        .sort_values(["start", "amino_acid_consequence"])
        .reset_index(drop=True)
    )

    expected = pd.DataFrame(
        [
            {
                "start": 4,
                "amino_acid_consequence": "K2A",
                "bcsq_aa_notation": "2K>2A",
                "presence_absence": "Y / N",
                "alt_freq": "0.100 / .",
            },
            {
                "start": 5,
                "amino_acid_consequence": "@4",
                "bcsq_aa_notation": "@4",
                "presence_absence": "Y / N",
                "alt_freq": "0.100 / .",
            },
            {
                "start": 5,
                "amino_acid_consequence": "K2T",
                "bcsq_aa_notation": "2K>2T",
                "presence_absence": "N / Y",
                "alt_freq": ". / 0.100",
            },
        ]
    )

    pd.testing.assert_frame_equal(observed, expected)


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not available")
def test_merge_then_annotate_preserves_joint_annotations_across_samples(tmp_path):
    ref, gff = _write_minimal_reference_bundle(tmp_path)

    sample1_vcf = tmp_path / "s1.vcf"
    sample2_vcf = tmp_path / "s2.vcf"
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=9>\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">\n'
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    )
    sample1_vcf.write_text(
        header.format(sample="s1")
        + "chr1\t4\t.\tA\tG\t.\tPASS\tAF=0.10;DP=100\tGT:DP:AF\t1:100:0.10\n"
        + "chr1\t5\t.\tA\tC\t.\tPASS\tAF=0.10;DP=100\tGT:DP:AF\t1:100:0.10\n",
        encoding="utf-8",
    )
    sample2_vcf.write_text(
        header.format(sample="s2")
        + "chr1\t5\t.\tA\tC\t.\tPASS\tAF=0.10;DP=100\tGT:DP:AF\t1:100:0.10\n",
        encoding="utf-8",
    )

    for path in (sample1_vcf, sample2_vcf):
        gz_path = Path(f"{path}.gz")
        subprocess.run(
            ["bcftools", "view", "-Oz", "-o", str(gz_path), str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
        subprocess.run(
            ["bcftools", "index", "-f", str(gz_path)],
            check=True,
            capture_output=True,
            text=True,
        )

    vcf_list = tmp_path / "vcf_list.txt"
    vcf_list.write_text(
        f"{sample1_vcf}.gz\n{sample2_vcf}.gz\n",
        encoding="utf-8",
    )
    sample_names = tmp_path / "sample_names.txt"
    sample_names.write_text("s1\ns2\n", encoding="utf-8")

    merged_vcf = tmp_path / "merged.vcf"
    merge_consequences(str(tmp_path), str(merged_vcf), str(sample_names), debug=False)

    annotated_vcf = tmp_path / "annotated.vcf"
    annotate_vcf(
        str(merged_vcf),
        str(annotated_vcf),
        str(ref),
        str(gff),
        debug=False,
    )

    cov1 = tmp_path / "s1.depth.txt"
    cov2 = tmp_path / "s2.depth.txt"
    _write_depth_file(cov1)
    _write_depth_file(cov2)

    table = process_vcf(
        str(annotated_vcf),
        [str(cov1), str(cov2)],
        10,
        ["s1", "s2"],
    )
    results_csv = tmp_path / "results.csv"
    table.to_csv(results_csv, index=False)
    processed = process_joint_variants(str(results_csv))

    observed = (
        processed[
            [
                "start",
                "variant",
                "amino_acid_consequence",
                "presence_absence",
                "type_of_change",
                "joint_variant",
            ]
        ]
        .sort_values(["start", "amino_acid_consequence", "presence_absence"])
        .reset_index(drop=True)
    )

    expected = pd.DataFrame(
        [
            {
                "start": 4,
                "variant": "A4G",
                "amino_acid_consequence": "K2A",
                "presence_absence": "Y / N",
                "type_of_change": "joint_missense",
                "joint_variant": True,
            },
            {
                "start": 5,
                "variant": "A5C",
                "amino_acid_consequence": "K2A",
                "presence_absence": "Y / N",
                "type_of_change": "joint_missense",
                "joint_variant": True,
            },
            {
                "start": 5,
                "variant": "A5C",
                "amino_acid_consequence": "K2T",
                "presence_absence": "N / Y",
                "type_of_change": "missense",
                "joint_variant": False,
            },
        ]
    )

    pd.testing.assert_frame_equal(observed, expected)
