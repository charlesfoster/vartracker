"""Unit tests for vartracker.main helpers."""

from argparse import Namespace
from unittest import mock

import importlib
import pandas as pd
import pytest

from vartracker.core import DependencyError

main_module = importlib.import_module("vartracker.main")


def test_setup_default_paths_uses_package_defaults():
    args = Namespace(reference=None, gff3=None)

    updated = main_module.setup_default_paths(args)

    assert updated.reference.endswith("NC_045512.fasta")
    assert updated.gff3.endswith("NC_045512.gff3")


def test_setup_default_paths_preserves_explicit_values():
    args = Namespace(reference="/tmp/custom.fasta", gff3="/tmp/custom.gff3")

    updated = main_module.setup_default_paths(args)

    assert updated.reference == "/tmp/custom.fasta"
    assert updated.gff3 == "/tmp/custom.gff3"


@mock.patch.object(
    main_module, "check_dependencies", return_value={"bcftools": True, "tabix": True}
)
def test_validate_dependencies_passes_when_available(mock_check):
    main_module.validate_dependencies()
    mock_check.assert_called_once()


@mock.patch.object(
    main_module,
    "check_dependencies",
    return_value={"bcftools": False, "tabix": True},
)
def test_validate_dependencies_raises_when_missing(mock_check):
    with pytest.raises(DependencyError, match="Missing required tools: bcftools"):
        main_module.validate_dependencies()
    mock_check.assert_called_once()


@pytest.fixture()
def minimal_vcf(tmp_path):
    vcf_path = tmp_path / "sample.vcf"
    vcf_path.write_text(
        """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">
##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Variant AF\">
##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
NC_045512.2\t266\t.\tA\tC\t.\tPASS\tDP=10;VAF=0.5
""",
        encoding="utf-8",
    )
    return vcf_path


def test_main_resolves_relative_paths(tmp_path, monkeypatch, minimal_vcf):
    coverage_path = tmp_path / "sample.cov.txt"
    coverage_path.write_text("NC_045512.2\t266\t100\n", encoding="utf-8")

    csv_path = tmp_path / "inputs.csv"
    csv_path.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        "Sample1,0,,,,sample.vcf,sample.cov.txt\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(main_module, "validate_dependencies", lambda: None)

    def fake_setup(args):
        args.reference = "/tmp/mock_reference.fasta"
        args.gff3 = "/tmp/mock_annotation.gff3"
        return args

    monkeypatch.setattr(main_module, "setup_default_paths", fake_setup)
    monkeypatch.setattr(
        main_module, "validate_reference_and_annotation", lambda *a, **k: None
    )
    monkeypatch.setattr(
        main_module, "generate_cumulative_lineplot", lambda *a, **k: None
    )
    monkeypatch.setattr(main_module, "generate_variant_heatmap", lambda *a, **k: None)
    monkeypatch.setattr(
        main_module, "process_joint_variants", lambda path: pd.read_csv(path)
    )
    monkeypatch.setattr(
        main_module, "generate_gene_table", lambda table, *_a, **_k: table
    )
    monkeypatch.setattr(main_module, "plot_gene_table", lambda *a, **k: None)
    monkeypatch.setattr(main_module, "search_pokay", lambda *a, **k: None)

    formatted_csq = tmp_path / "formatted.csq.vcf.gz"
    monkeypatch.setattr(
        main_module,
        "format_vcf",
        lambda *a, **k: (str(tmp_path / "formatted.vcf.gz"), str(formatted_csq)),
    )
    monkeypatch.setattr(main_module, "merge_consequences", lambda *a, **k: None)
    monkeypatch.setattr(
        main_module,
        "process_vcf",
        lambda *a, **k: pd.DataFrame(
            {
                "gene": ["S"],
                "variant": ["A266C"],
                "amino_acid_consequence": ["S:A1C"],
                "nsp_aa_change": [""],
                "presence_absence": ["Y"],
                "variant_status": ["new"],
                "persistence_status": ["new_persistent"],
                "samples": ["Sample1"],
                "alt_freq": ["0.5"],
            }
        ),
    )

    exit_code = main_module.main(
        [
            "vcf",
            str(csv_path),
            "--name",
            "Example",
        ]
    )

    assert exit_code == 0


def test_e2e_runs_snakemake_then_vcf(monkeypatch, tmp_path):
    updated_csv = tmp_path / "samples_updated.csv"
    vcf_out = tmp_path / "vcf.gz"
    cov_out = tmp_path / "coverage.txt"
    vcf_out.write_text("", encoding="utf-8")
    cov_out.write_text("", encoding="utf-8")
    updated_csv.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        f"Sample1,0,,,,{vcf_out},{cov_out}\n",
        encoding="utf-8",
    )

    recorded = {}

    def fake_workflow(**kwargs):
        recorded["workflow_kwargs"] = kwargs
        return str(updated_csv)

    monkeypatch.setattr(main_module, "run_e2e_workflow", fake_workflow)

    def fake_vcf(args):
        recorded["vcf_input"] = args.input_csv
        return 0

    monkeypatch.setattr(main_module, "_run_vcf_command", fake_vcf)

    reads1 = tmp_path / "reads1.fastq"
    reads2 = tmp_path / "reads2.fastq"
    reads1.write_text("", encoding="utf-8")
    reads2.write_text("", encoding="utf-8")
    samples_csv = tmp_path / "reads.csv"
    samples_csv.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        f"Sample1,0,{reads1},{reads2},,,\n",
        encoding="utf-8",
    )

    exit_code = main_module.main(
        [
            "end-to-end",
            str(samples_csv),
            "--reference",
            "ref.fasta",
            "--snakemake-outdir",
            str(tmp_path / "snakemake"),
            "--outdir",
            str(tmp_path / "results"),
        ]
    )

    assert exit_code == 0
    assert recorded["workflow_kwargs"]["samples_csv"] == str(samples_csv)
    assert recorded["workflow_kwargs"]["force_all"] is False
    assert recorded["workflow_kwargs"]["quiet"] is True
    assert recorded["workflow_kwargs"]["mode"] == "reads"
    assert recorded["vcf_input"] == str(updated_csv)


def test_e2e_dryrun_skips_vcf(monkeypatch, tmp_path_factory):
    reads1 = tmp_path_factory.mktemp("reads") / "reads1.fastq"
    reads2 = reads1.parent / "reads2.fastq"
    reads1.write_text("", encoding="utf-8")
    reads2.write_text("", encoding="utf-8")
    samples_csv = reads1.parent / "reads.csv"
    samples_csv.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        f"Sample1,0,{reads1},{reads2},,,\n",
        encoding="utf-8",
    )

    def fake_workflow(**kwargs):
        assert kwargs["force_all"] is True
        assert kwargs["quiet"] is True
        assert kwargs["mode"] == "reads"
        return None

    monkeypatch.setattr(main_module, "run_e2e_workflow", fake_workflow)

    def fake_vcf(_):
        raise AssertionError("VCF pipeline should not run during dry-run")

    monkeypatch.setattr(main_module, "_run_vcf_command", fake_vcf)

    exit_code = main_module.main(
        [
            "end-to-end",
            str(samples_csv),
            "--reference",
            "ref.fasta",
            "--snakemake-dryrun",
            "--redo",
        ]
    )

    assert exit_code == 0


def test_e2e_verbose_enables_snakemake_commands(monkeypatch, tmp_path):
    updated_csv = tmp_path / "samples_updated.csv"
    vcf_out = tmp_path / "vcf4.gz"
    cov_out = tmp_path / "coverage4.txt"
    vcf_out.write_text("", encoding="utf-8")
    cov_out.write_text("", encoding="utf-8")
    updated_csv.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        f"Sample1,0,,,,{vcf_out},{cov_out}\n",
        encoding="utf-8",
    )

    recorded = {}

    reads1 = tmp_path / "reads1.fastq"
    reads2 = tmp_path / "reads2.fastq"
    reads1.write_text("", encoding="utf-8")
    reads2.write_text("", encoding="utf-8")
    samples_csv = tmp_path / "reads.csv"
    samples_csv.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        f"Sample1,0,{reads1},{reads2},,,\n",
        encoding="utf-8",
    )

    def fake_workflow(**kwargs):
        recorded.update(kwargs)
        return str(updated_csv)

    monkeypatch.setattr(main_module, "run_e2e_workflow", fake_workflow)

    def fake_vcf(args):
        recorded["vcf_input"] = args.input_csv
        recorded["suppress_logo"] = getattr(args, "_suppress_logo", False)
        return 0

    monkeypatch.setattr(main_module, "_run_vcf_command", fake_vcf)

    exit_code = main_module.main(
        [
            "end-to-end",
            str(samples_csv),
            "--reference",
            "ref.fasta",
            "--verbose",
        ]
    )

    assert exit_code == 0
    assert recorded["quiet"] is False
    assert recorded["mode"] == "reads"
    assert recorded["vcf_input"] == str(updated_csv)
    assert recorded["suppress_logo"] is True


def test_bam_runs_snakemake_then_vcf(monkeypatch, tmp_path):
    bam_path = tmp_path / "sample.bam"
    bam_path.write_text("", encoding="utf-8")
    samples_csv = tmp_path / "bam.csv"
    samples_csv.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        f"Sample1,0,,,{bam_path},,\n",
        encoding="utf-8",
    )

    updated_csv = tmp_path / "samples_updated.csv"
    vcf_out = tmp_path / "bam_out.vcf.gz"
    cov_out = tmp_path / "bam_out.coverage.txt"
    vcf_out.write_text("", encoding="utf-8")
    cov_out.write_text("", encoding="utf-8")
    updated_csv.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        f"Sample1,0,,,,{vcf_out},{cov_out}\n",
        encoding="utf-8",
    )

    recorded = {}

    def fake_workflow(**kwargs):
        recorded.update(kwargs)
        return str(updated_csv)

    monkeypatch.setattr(main_module, "run_e2e_workflow", fake_workflow)

    def fake_vcf(args):
        recorded["vcf_input"] = args.input_csv
        recorded["suppress_logo"] = getattr(args, "_suppress_logo", False)
        return 0

    monkeypatch.setattr(main_module, "_run_vcf_command", fake_vcf)

    exit_code = main_module.main(
        [
            "bam",
            str(samples_csv),
            "--reference",
            "ref.fasta",
        ]
    )

    assert exit_code == 0
    assert recorded["samples_csv"] == str(samples_csv)
    assert recorded["force_all"] is False
    assert recorded["quiet"] is True
    assert recorded["mode"] == "bam"
    assert recorded["vcf_input"] == str(updated_csv)
    assert recorded["suppress_logo"] is True


def test_bam_dryrun_skips_vcf(monkeypatch, tmp_path):
    bam_path = tmp_path / "sample.bam"
    bam_path.write_text("", encoding="utf-8")
    samples_csv = tmp_path / "bam.csv"
    samples_csv.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        f"Sample1,0,,,{bam_path},,\n",
        encoding="utf-8",
    )

    def fake_workflow(**kwargs):
        assert kwargs["mode"] == "bam"
        return str(tmp_path / "dummy.csv")

    monkeypatch.setattr(main_module, "run_e2e_workflow", fake_workflow)

    def fake_vcf(_):
        raise AssertionError("VCF stage should not run during dry-run")

    monkeypatch.setattr(main_module, "_run_vcf_command", fake_vcf)

    exit_code = main_module.main(
        [
            "bam",
            str(samples_csv),
            "--reference",
            "ref.fasta",
            "--snakemake-dryrun",
        ]
    )

    assert exit_code == 0
