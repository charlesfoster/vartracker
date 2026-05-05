"""Unit tests for vartracker.main helpers."""

from argparse import Namespace
from pathlib import Path
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


def test_vcf_parser_accepts_multiallelic_overflow_option():
    parser = main_module.create_parser()

    args = parser.parse_args(
        ["vcf", "inputs.csv", "--multiallelic-overflow", "skip-site"]
    )

    assert args.multiallelic_overflow == "skip-site"


def test_bam_parser_accepts_lofreq_primer_rescue_options():
    parser = main_module.create_parser()

    args = parser.parse_args(
        [
            "bam",
            "inputs.csv",
            "--primer-bed",
            "primers.bed",
            "--lofreq-primer-rescue",
            "off",
            "--lofreq-rescue-min-af",
            "0.9",
        ]
    )

    assert args.primer_bed == "primers.bed"
    assert args.lofreq_primer_rescue == "off"
    assert args.lofreq_rescue_min_af == 0.9


def test_drop_exact_duplicate_result_rows_removes_only_exact_duplicates(capsys):
    table = pd.DataFrame(
        [
            {"variant": "A1C", "gene": "S", "amino_acid_consequence": "S:A1C"},
            {"variant": "A1C", "gene": "S", "amino_acid_consequence": "S:A1C"},
            {"variant": "A1C", "gene": "N", "amino_acid_consequence": "N:A1C"},
        ]
    )

    deduped = main_module._drop_exact_duplicate_result_rows(table)

    assert len(deduped) == 2
    assert deduped.to_dict(orient="records") == [
        {"variant": "A1C", "gene": "S", "amino_acid_consequence": "S:A1C"},
        {"variant": "A1C", "gene": "N", "amino_acid_consequence": "N:A1C"},
    ]
    assert capsys.readouterr().out == ""


def test_prepare_reference_command_invokes_bundle(monkeypatch, tmp_path):
    recorded = {}

    def fake_parse_accessions(**kwargs):
        recorded["parse_kwargs"] = kwargs
        return ["CY114381"]

    def fake_prepare_reference_bundle(**kwargs):
        recorded["bundle_kwargs"] = kwargs
        return {
            "outputs": {
                "fasta": str(tmp_path / "reference.fa"),
                "gff3": str(tmp_path / "reference.gff3"),
                "fai": str(tmp_path / "reference.fa.fai"),
                "metadata": str(tmp_path / "prepare_metadata.json"),
            }
        }

    monkeypatch.setattr(main_module, "parse_accessions", fake_parse_accessions)
    monkeypatch.setattr(
        main_module, "prepare_reference_bundle", fake_prepare_reference_bundle
    )

    exit_code = main_module.main(
        [
            "prepare",
            "reference",
            "--accessions",
            "CY114381",
            "--outdir",
            str(tmp_path),
        ]
    )

    assert exit_code == 0
    assert recorded["parse_kwargs"]["accessions"] == "CY114381"
    assert recorded["bundle_kwargs"]["accessions"] == ["CY114381"]
    assert recorded["bundle_kwargs"]["outdir"] == str(tmp_path)


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
    with pytest.raises(DependencyError) as excinfo:
        main_module.validate_dependencies()
    message = str(excinfo.value)
    assert (
        "Dependency error: to run the vcf mode the following tools must be available: bcftools"
        in message
    )
    assert "conda install -c bioconda bcftools" in excinfo.value.tip
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

    monkeypatch.setattr(main_module, "validate_dependencies", lambda mode="vcf": None)

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
    monkeypatch.setattr(main_module, "search_literature", lambda *a, **k: None)

    formatted_csq = tmp_path / "formatted.csq.vcf.gz"
    monkeypatch.setattr(
        main_module,
        "format_vcf",
        lambda *a, **k: (str(tmp_path / "formatted.vcf.gz"), str(formatted_csq)),
    )
    monkeypatch.setattr(main_module, "merge_consequences", lambda *a, **k: None)
    monkeypatch.setattr(main_module, "annotate_vcf", lambda *a, **k: None)
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
            "--outdir",
            str(tmp_path / "results"),
        ]
    )

    assert exit_code == 0


def test_search_pokay_retains_parsed_database_csv(tmp_path, monkeypatch, minimal_vcf):
    coverage_path = tmp_path / "sample.cov.txt"
    coverage_path.write_text("NC_045512.2\t266\t100\n", encoding="utf-8")

    csv_path = tmp_path / "inputs.csv"
    csv_path.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        "Sample1,0,,,,sample.vcf,sample.cov.txt\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(main_module, "validate_dependencies", lambda mode="vcf": None)

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
    monkeypatch.setattr(
        main_module, "search_literature", lambda *a, **k: pd.DataFrame()
    )

    formatted_csq = tmp_path / "formatted.csq.vcf.gz"
    monkeypatch.setattr(
        main_module,
        "format_vcf",
        lambda *a, **k: (str(tmp_path / "formatted.vcf.gz"), str(formatted_csq)),
    )
    monkeypatch.setattr(main_module, "merge_consequences", lambda *a, **k: None)
    monkeypatch.setattr(main_module, "annotate_vcf", lambda *a, **k: None)
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

    def fake_parse_pokay(argv):
        output_path = Path(argv[0])
        output_path.write_text(
            "gene,mutation,information,reference\nS,S:A1C,Mock hit,PMID123\n",
            encoding="utf-8",
        )
        return 0

    monkeypatch.setattr(main_module.parse_pokay_module, "main", fake_parse_pokay)

    outdir = tmp_path / "results"
    exit_code = main_module.main(
        [
            "vcf",
            str(csv_path),
            "--name",
            "Example",
            "--outdir",
            str(outdir),
            "--search-pokay",
        ]
    )

    assert exit_code == 0
    assert (outdir / "literature_database.csv").exists()


def test_vcf_workflow_uses_default_heatmap_options(monkeypatch, tmp_path, minimal_vcf):
    coverage_path = tmp_path / "sample.cov.txt"
    coverage_path.write_text("NC_045512.2\t266\t100\n", encoding="utf-8")

    csv_path = tmp_path / "inputs.csv"
    csv_path.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        "Sample1,0,,,,sample.vcf,sample.cov.txt\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(main_module, "validate_dependencies", lambda mode="vcf": None)

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
    monkeypatch.setattr(
        main_module, "process_joint_variants", lambda path: pd.read_csv(path)
    )
    monkeypatch.setattr(
        main_module, "generate_gene_table", lambda table, *_a, **_k: table
    )
    monkeypatch.setattr(main_module, "plot_gene_table", lambda *a, **k: None)
    monkeypatch.setattr(main_module, "search_literature", lambda *a, **k: None)

    recorded = {}

    def fake_heatmap(*args, **kwargs):
        recorded.update(kwargs)

    monkeypatch.setattr(main_module, "generate_variant_heatmap", fake_heatmap)

    formatted_csq = tmp_path / "formatted.csq.vcf.gz"
    monkeypatch.setattr(
        main_module,
        "format_vcf",
        lambda *a, **k: (str(tmp_path / "formatted.vcf.gz"), str(formatted_csq)),
    )
    monkeypatch.setattr(main_module, "merge_consequences", lambda *a, **k: None)
    monkeypatch.setattr(main_module, "annotate_vcf", lambda *a, **k: None)
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
            "--outdir",
            str(tmp_path / "results"),
        ]
    )

    assert exit_code == 0
    for key in (
        "excluded_consequence_types",
        "included_consequence_types",
        "include_joint",
        "only_persistent",
        "only_new",
        "gene_include",
        "gene_exclude",
        "variant_type_include",
        "qc_include",
        "min_prop_passing_qc",
        "min_persistence",
        "min_max_af",
        "min_sample_af",
        "sample_subset",
        "hide_singletons",
        "min_depth",
    ):
        assert key not in recorded


def test_plot_heatmap_replots_results_csv(monkeypatch, tmp_path):
    results_csv = tmp_path / "results.csv"
    results_csv.write_text(
        "samples,sample_number,name,alt_freq,variant_site_depth,presence_absence,variant_status,persistence_status,type_of_variant,type_of_change,gene,variant,start\n"
        "P0 / P1,0 / 1,Example,0.0 / 0.5,0 / 100,N / Y,new,new_persistent,snp,missense,S,A266C,266\n",
        encoding="utf-8",
    )

    recorded = {}
    literature_csv = tmp_path / "literature_hits.csv"
    literature_csv.write_text(
        "gene,amino_acid_consequence,information,reference\n"
        "S,S:A266C,Mock evidence,PMID123\n",
        encoding="utf-8",
    )

    def fake_heatmap(*args, **kwargs):
        recorded["table"] = args[0]
        recorded["sample_names"] = args[1]
        recorded["sample_numbers"] = args[2]
        recorded["outdir"] = args[3]
        recorded["project_name"] = args[4]
        recorded["min_snv_freq"] = args[5]
        recorded["min_indel_freq"] = args[6]
        recorded["kwargs"] = kwargs

    monkeypatch.setattr(main_module, "generate_variant_heatmap", fake_heatmap)

    exit_code = main_module.main(
        [
            "plot",
            "heatmap",
            str(results_csv),
            "--aa-exclude",
            "*frameshift*",
            "--include-joint",
            "--title",
            "Custom heatmap",
            "--x-labels",
            "sample-number",
            "--literature-csv",
            str(literature_csv),
        ]
    )

    assert exit_code == 0
    assert list(recorded["sample_names"]) == ["P0", "P1"]
    assert list(recorded["sample_numbers"]) == ["0", "1"]
    assert recorded["outdir"] == str(tmp_path)
    assert recorded["project_name"] == ""
    assert recorded["min_snv_freq"] == 0.03
    assert recorded["min_indel_freq"] == 0.1
    assert recorded["kwargs"]["excluded_consequence_types"] == ["*frameshift*"]
    assert recorded["kwargs"]["include_joint"] is True
    assert recorded["kwargs"]["x_tick_labels"] == ["0", "1"]
    assert recorded["kwargs"]["plot_title"] == "Custom heatmap"
    assert recorded["kwargs"]["literature_table_path"] == str(literature_csv.resolve())
    assert recorded["kwargs"]["literature_hits"]["gene"].tolist() == ["S"]


def _write_plot_results_csv(path: Path, n_variants: int = 4) -> None:
    rows = [
        "chrom,start,end,samples,sample_number,name,alt_freq,per_sample_variant_qc,presence_absence,variant_status,persistence_status,type_of_variant,type_of_change,gene,variant,amino_acid_consequence,nsp_aa_change,reference"
    ]
    for idx in range(n_variants):
        gene = "S" if idx % 2 == 0 else "N"
        chrom = "segA" if gene == "S" else "segB"
        aa = f"{gene}:V{idx + 1}A"
        effect = "missense" if idx % 3 else "synonymous"
        persistence = "new_persistent" if idx % 2 == 0 else "new_transient"
        reference = "PMID123" if idx == 1 else ""
        rows.append(
            f"{chrom},"
            f"{100 + idx},"
            f"{100 + idx},"
            "P0 / P1 / P2 / P3,"
            "0 / 1 / 2 / 3,"
            "Example,"
            f"0.0 / 0.{idx + 2} / 0.{idx + 3} / 0.{idx + 4},"
            "P / P / P / P,"
            "N / Y / Y / Y,"
            "new,"
            f"{persistence},"
            "snp,"
            f"{effect},"
            f"{gene},"
            f"A{100 + idx}G,"
            f"{aa},,"
            f"{reference}"
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _write_reference_features_json(path: Path) -> None:
    path.write_text(
        """{
  "contig_lengths": {"segA": 30000, "segB": 15000},
  "contig_order": ["segA", "segB"],
  "features": [
    {"aa_length": 1200, "contig": "segA", "end": 25000, "name": "S", "start": 21000, "strand": "+", "type": "mRNA"},
    {"aa_length": 400, "contig": "segB", "end": 12000, "name": "N", "start": 9000, "strand": "+", "type": "mRNA"}
  ],
  "source_gff3": "/tmp/reference.gff3"
}
""",
        encoding="utf-8",
    )


def test_plot_trajectory_variants_override_auto_selection(monkeypatch, tmp_path):
    results_csv = tmp_path / "results.csv"
    _write_plot_results_csv(results_csv, n_variants=6)

    recorded = {}

    def fake_plot(summary, _long_df, **kwargs):
        selected_ids = kwargs["selected_variants"]
        recorded["selected_labels"] = (
            summary.set_index("variant_id")
            .loc[selected_ids, "variant_label"]
            .str.replace(r" \(.*\)$", "", regex=True)
            .tolist()
        )

    monkeypatch.setattr(main_module, "plot_variant_trajectory", fake_plot)

    exit_code = main_module.main(
        [
            "plot",
            "trajectory",
            str(results_csv),
            "--variants",
            "S:V3A,N:V2A",
        ]
    )

    assert exit_code == 0
    assert recorded["selected_labels"] == ["S:V3A", "N:V2A"]


def test_plot_trajectory_threshold_options_are_forwarded(monkeypatch, tmp_path):
    results_csv = tmp_path / "results.csv"
    _write_plot_results_csv(results_csv, n_variants=6)

    recorded = {}

    def fake_plot(*args, **kwargs):
        recorded.update(kwargs)

    monkeypatch.setattr(main_module, "plot_variant_trajectory", fake_plot)

    exit_code = main_module.main(
        [
            "plot",
            "trajectory",
            str(results_csv),
            "--thresholds",
            "0.5,0.9",
            "--label-threshold-crossers",
            "--crossing-rule",
            "strictly_above",
        ]
    )

    assert exit_code == 0
    assert recorded["thresholds"] == [0.5, 0.9]
    assert recorded["label_threshold_crossers"] is True
    assert recorded["crossing_rule"] == "strictly_above"


def test_plot_trajectory_label_mode_is_forwarded(monkeypatch, tmp_path):
    results_csv = tmp_path / "results.csv"
    _write_plot_results_csv(results_csv, n_variants=6)

    recorded = {}

    def fake_plot(*args, **kwargs):
        recorded.update(kwargs)

    monkeypatch.setattr(main_module, "plot_variant_trajectory", fake_plot)

    exit_code = main_module.main(
        [
            "plot",
            "trajectory",
            str(results_csv),
            "--label-mode",
            "nt",
        ]
    )

    assert exit_code == 0
    assert recorded["label_mode"] == "nt"


def test_plot_sample_axis_mode_is_forwarded(monkeypatch, tmp_path):
    results_csv = tmp_path / "results.csv"
    _write_plot_results_csv(results_csv, n_variants=6)

    recorded = {}

    def fake_turnover(*args, **kwargs):
        recorded["turnover"] = kwargs["sample_axis_mode"]

    def fake_lifespan(*args, **kwargs):
        recorded["lifespan"] = kwargs["sample_axis_mode"]

    monkeypatch.setattr(main_module, "plot_variant_turnover", fake_turnover)
    monkeypatch.setattr(main_module, "plot_variant_lifespan", fake_lifespan)

    assert (
        main_module.main(
            ["plot", "turnover", str(results_csv), "--sample-axis", "name"]
        )
        == 0
    )
    assert (
        main_module.main(
            ["plot", "lifespan", str(results_csv), "--sample-axis", "name"]
        )
        == 0
    )

    assert recorded["turnover"] == "name"
    assert recorded["lifespan"] == "name"


def test_plot_trajectory_crossing_only_requires_thresholds(tmp_path):
    results_csv = tmp_path / "results.csv"
    _write_plot_results_csv(results_csv, n_variants=4)

    exit_code = main_module.main(
        ["plot", "trajectory", str(results_csv), "--crossing-only"]
    )

    assert exit_code == 1


def test_plot_trajectory_crossing_only_filters_variants(monkeypatch, tmp_path):
    results_csv = tmp_path / "results.csv"
    _write_plot_results_csv(results_csv, n_variants=6)

    recorded = {}

    def fake_plot(summary, _long_df, **kwargs):
        selected_ids = kwargs["selected_variants"]
        recorded["selected_labels"] = (
            summary.set_index("variant_id")
            .loc[selected_ids, "variant_label"]
            .str.replace(r" \(.*\)$", "", regex=True)
            .tolist()
        )
        recorded["summary_variant_ids"] = list(summary["variant_id"])

    monkeypatch.setattr(main_module, "plot_variant_trajectory", fake_plot)

    exit_code = main_module.main(
        [
            "plot",
            "trajectory",
            str(results_csv),
            "--thresholds",
            "0.5",
            "--crossing-only",
        ]
    )

    assert exit_code == 0
    assert recorded["selected_labels"]
    assert len(recorded["selected_labels"]) <= main_module.DEFAULT_TRAJECTORY_TOP_N


def test_plot_turnover_uses_all_filtered_variants_by_default(monkeypatch, tmp_path):
    results_csv = tmp_path / "results.csv"
    _write_plot_results_csv(results_csv, n_variants=5)

    recorded = {}

    def fake_turnover(summary, long_df, **kwargs):
        recorded["summary_count"] = len(summary)
        recorded["variant_labels"] = (
            summary["variant_label"].str.replace(r" \(.*\)$", "", regex=True).tolist()
        )

    monkeypatch.setattr(main_module, "plot_variant_turnover", fake_turnover)

    exit_code = main_module.main(["plot", "turnover", str(results_csv), "--gene", "S"])

    assert exit_code == 0
    assert recorded["summary_count"] == 3
    assert all(label.startswith("S:") for label in recorded["variant_labels"])


def test_plot_commands_auto_limit_variants_by_top_n(monkeypatch, tmp_path):
    results_csv = tmp_path / "results.csv"
    _write_plot_results_csv(results_csv, n_variants=20)

    recorded = {}

    def fake_trajectory(*args, **kwargs):
        recorded["trajectory"] = kwargs["selected_variants"]

    def fake_lifespan(*args, **kwargs):
        recorded["lifespan"] = kwargs["selected_variants"]

    monkeypatch.setattr(main_module, "plot_variant_trajectory", fake_trajectory)
    monkeypatch.setattr(main_module, "plot_variant_lifespan", fake_lifespan)

    assert main_module.main(["plot", "trajectory", str(results_csv)]) == 0
    assert main_module.main(["plot", "lifespan", str(results_csv)]) == 0

    assert len(recorded["trajectory"]) == main_module.DEFAULT_TRAJECTORY_TOP_N
    assert len(recorded["lifespan"]) == main_module.DEFAULT_LIFESPAN_TOP_N


def test_plot_genome_uses_sidecar_metadata(monkeypatch, tmp_path):
    results_csv = tmp_path / "results.csv"
    _write_plot_results_csv(results_csv, n_variants=6)
    _write_reference_features_json(tmp_path / "reference_features.json")

    recorded = {}

    def fake_genome(table, metadata, **kwargs):
        recorded["rows"] = len(table)
        recorded["metadata"] = metadata
        recorded.update(kwargs)

    monkeypatch.setattr(main_module, "plot_variant_genome", fake_genome)

    exit_code = main_module.main(
        [
            "plot",
            "genome",
            str(results_csv),
            "--gene",
            "S",
            "--aa-scale",
            "--focus-coords",
            "50-120,180-220",
            "--min-af",
            "0.2",
        ]
    )

    assert exit_code == 0
    assert recorded["rows"] == 6
    assert recorded["gene"] == "S"
    assert recorded["aa_scale"] is True
    assert recorded["focus_ranges"] == [(50.0, 120.0, 0), (180.0, 220.0, 1)]
    assert recorded["focus_labels"] == [None, None]
    assert recorded["min_af"] == 0.2
    assert recorded["include_indels"] is False
    assert recorded["show_intersections"] is False
    assert recorded["metadata"]["contig_order"] == ["segA", "segB"]


def test_plot_genome_named_focus_regions_are_forwarded(monkeypatch, tmp_path):
    results_csv = tmp_path / "results.csv"
    _write_plot_results_csv(results_csv, n_variants=4)
    _write_reference_features_json(tmp_path / "reference_features.json")

    recorded = {}

    def fake_genome(_table, _metadata, **kwargs):
        recorded.update(kwargs)

    monkeypatch.setattr(main_module, "plot_variant_genome", fake_genome)

    exit_code = main_module.main(
        [
            "plot",
            "genome",
            str(results_csv),
            "--focus-coords",
            "I:50-120,180-220;II:300-330",
        ]
    )

    assert exit_code == 0
    assert recorded["focus_ranges"] == [
        (50.0, 120.0, 0),
        (180.0, 220.0, 0),
        (300.0, 330.0, 1),
    ]
    assert recorded["focus_labels"] == ["I", "II"]


def test_plot_genome_cds_scale_is_forwarded(monkeypatch, tmp_path):
    results_csv = tmp_path / "results.csv"
    _write_plot_results_csv(results_csv, n_variants=4)
    _write_reference_features_json(tmp_path / "reference_features.json")

    recorded = {}

    def fake_genome(_table, _metadata, **kwargs):
        recorded.update(kwargs)

    monkeypatch.setattr(main_module, "plot_variant_genome", fake_genome)

    exit_code = main_module.main(
        [
            "plot",
            "genome",
            str(results_csv),
            "--gene",
            "S",
            "--cds-scale",
            "--focus-coords",
            "10-30",
        ]
    )

    assert exit_code == 0
    assert recorded["gene"] == "S"
    assert recorded["cds_scale"] is True
    assert recorded["focus_ranges"] == [(10.0, 30.0, 0)]


def test_plot_genome_show_intersections_is_forwarded(monkeypatch, tmp_path):
    results_csv = tmp_path / "results.csv"
    _write_plot_results_csv(results_csv, n_variants=4)
    _write_reference_features_json(tmp_path / "reference_features.json")

    recorded = {}

    def fake_genome(_table, _metadata, **kwargs):
        recorded.update(kwargs)

    monkeypatch.setattr(main_module, "plot_variant_genome", fake_genome)

    exit_code = main_module.main(
        [
            "plot",
            "genome",
            str(results_csv),
            "--focus-coords",
            "I:50-120,180-220",
            "--show-intersections",
        ]
    )

    assert exit_code == 0
    assert recorded["show_intersections"] is True


def test_plot_genome_include_indels_is_forwarded(monkeypatch, tmp_path):
    results_csv = tmp_path / "results.csv"
    _write_plot_results_csv(results_csv, n_variants=4)
    _write_reference_features_json(tmp_path / "reference_features.json")

    recorded = {}

    def fake_genome(_table, _metadata, **kwargs):
        recorded.update(kwargs)

    monkeypatch.setattr(main_module, "plot_variant_genome", fake_genome)

    exit_code = main_module.main(
        ["plot", "genome", str(results_csv), "--include-indels"]
    )

    assert exit_code == 0
    assert recorded["include_indels"] is True


def test_plot_genome_aa_scale_requires_gene(tmp_path):
    results_csv = tmp_path / "results.csv"
    _write_plot_results_csv(results_csv, n_variants=4)
    _write_reference_features_json(tmp_path / "reference_features.json")

    exit_code = main_module.main(["plot", "genome", str(results_csv), "--aa-scale"])

    assert exit_code == 1


def test_vcf_workflow_writes_default_genome_plot(monkeypatch, tmp_path, minimal_vcf):
    coverage_path = tmp_path / "sample.cov.txt"
    coverage_path.write_text("NC_045512.2\t266\t100\n", encoding="utf-8")

    csv_path = tmp_path / "inputs.csv"
    csv_path.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        "Sample1,0,,,,sample.vcf,sample.cov.txt\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(main_module, "validate_dependencies", lambda mode="vcf": None)

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
    monkeypatch.setattr(main_module, "search_literature", lambda *a, **k: None)
    monkeypatch.setattr(
        main_module, "write_reference_feature_metadata", lambda *a, **k: None
    )

    recorded = {}

    def fake_genome(*args, **kwargs):
        recorded["output_path"] = kwargs["output_path"]

    monkeypatch.setattr(main_module, "plot_variant_genome", fake_genome)
    monkeypatch.setattr(
        main_module,
        "load_reference_feature_metadata",
        lambda *_a, **_k: {
            "contig_order": ["NC_045512.2"],
            "contig_lengths": {"NC_045512.2": 29903},
            "features": [],
        },
    )

    formatted_csq = tmp_path / "formatted.csq.vcf.gz"
    monkeypatch.setattr(
        main_module,
        "format_vcf",
        lambda *a, **k: (str(tmp_path / "formatted.vcf.gz"), str(formatted_csq)),
    )
    monkeypatch.setattr(main_module, "merge_consequences", lambda *a, **k: None)
    monkeypatch.setattr(main_module, "annotate_vcf", lambda *a, **k: None)
    monkeypatch.setattr(
        main_module,
        "process_vcf",
        lambda *a, **k: pd.DataFrame(
            {
                "chrom": ["NC_045512.2"],
                "start": [266],
                "end": [266],
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
            "--outdir",
            str(tmp_path / "results"),
        ]
    )

    assert exit_code == 0
    assert recorded["output_path"].endswith("variant_genome_plot.pdf")


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

    modes_checked = []

    def fake_validate(mode="vcf"):
        modes_checked.append(mode)

    monkeypatch.setattr(main_module, "validate_dependencies", fake_validate)

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
            "--consensus-snp-min-af",
            "0.20",
            "--consensus-snp-thresh",
            "0.80",
            "--consensus-indel-thresh",
            "0.70",
            "--ampliconclip-tolerance",
            "2",
            "--outdir",
            str(tmp_path / "results"),
        ]
    )

    assert exit_code == 0
    assert recorded["workflow_kwargs"]["samples_csv"] == str(samples_csv)
    assert recorded["workflow_kwargs"]["force_all"] is False
    assert recorded["workflow_kwargs"]["quiet"] is True
    assert recorded["workflow_kwargs"]["mode"] == "reads"
    assert recorded["workflow_kwargs"]["min_depth"] == 10
    assert recorded["workflow_kwargs"]["consensus_snp_min_af"] == 0.20
    assert recorded["workflow_kwargs"]["consensus_snp_thresh"] == 0.80
    assert recorded["workflow_kwargs"]["consensus_indel_thresh"] == 0.70
    assert recorded["workflow_kwargs"]["ampliconclip_tolerance"] == 2
    assert recorded["workflow_kwargs"]["lofreq_primer_rescue"] == "auto"
    assert recorded["workflow_kwargs"]["lofreq_rescue_min_af"] == 0.95
    assert recorded["vcf_input"] == str(updated_csv)
    assert modes_checked == ["e2e"]


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

    modes_checked = []

    def fake_validate(mode="vcf"):
        modes_checked.append(mode)

    monkeypatch.setattr(main_module, "validate_dependencies", fake_validate)

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
            "--outdir",
            str(tmp_path_factory.mktemp("e2e_out")),
        ]
    )

    assert exit_code == 0
    assert modes_checked == ["e2e"]


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

    modes_checked = []

    def fake_validate(mode="vcf"):
        modes_checked.append(mode)

    monkeypatch.setattr(main_module, "validate_dependencies", fake_validate)

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
            "--outdir",
            str(tmp_path / "results"),
        ]
    )

    assert exit_code == 0
    assert recorded["quiet"] is False
    assert recorded["mode"] == "reads"
    assert recorded["vcf_input"] == str(updated_csv)
    assert recorded["suppress_logo"] is True
    assert modes_checked == ["e2e"]


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

    modes_checked = []

    def fake_validate(mode="vcf"):
        modes_checked.append(mode)

    monkeypatch.setattr(main_module, "validate_dependencies", fake_validate)

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
            "--min-snv-freq",
            "0.05",
            "--min-depth",
            "12",
            "--outdir",
            str(tmp_path / "results"),
        ]
    )

    assert exit_code == 0
    assert recorded["samples_csv"] == str(samples_csv)
    assert recorded["force_all"] is False
    assert recorded["quiet"] is True
    assert recorded["mode"] == "bam"
    assert recorded["min_depth"] == 12
    assert recorded["consensus_snp_min_af"] == 0.25
    assert recorded["consensus_snp_thresh"] == 0.75
    assert recorded["consensus_indel_thresh"] == 0.75
    assert recorded["primer_bed"] is None
    assert recorded["lofreq_primer_rescue"] == "auto"
    assert recorded["vcf_input"] == str(updated_csv)
    assert recorded["suppress_logo"] is True
    assert modes_checked == ["bam"]


def test_bam_dryrun_skips_vcf(monkeypatch, tmp_path):
    bam_path = tmp_path / "sample.bam"
    bam_path.write_text("", encoding="utf-8")
    samples_csv = tmp_path / "bam.csv"
    samples_csv.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        f"Sample1,0,,,{bam_path},,\n",
        encoding="utf-8",
    )

    modes_checked = []

    def fake_validate(mode="vcf"):
        modes_checked.append(mode)

    monkeypatch.setattr(main_module, "validate_dependencies", fake_validate)

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
            "--outdir",
            str(tmp_path / "results"),
        ]
    )

    assert exit_code == 0
    assert modes_checked == ["bam"]


def test_bam_test_mode_uses_demo_dataset(monkeypatch, tmp_path):
    recorded = {}

    def fake_validate(mode="vcf"):
        recorded.setdefault("modes", []).append(mode)

    monkeypatch.setattr(main_module, "validate_dependencies", fake_validate)

    mock_csv = tmp_path / "demo.csv"
    mock_bam = tmp_path / "demo.bam"
    mock_bam.write_text("", encoding="utf-8")
    mock_literature = tmp_path / "mock.csv"
    mock_literature.write_text("", encoding="utf-8")
    mock_csv.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        f"Demo,0,,,{mock_bam},,\n",
        encoding="utf-8",
    )

    class DummyTemp:
        def __init__(self):
            self.cleaned = False

        def cleanup(self):
            self.cleaned = True
            recorded["cleanup_called"] = True

    temp_holder = DummyTemp()

    def fake_prepare(args, mode):
        recorded["prepare_mode"] = mode
        args.search_pokay = False
        args.literature_csv = str(mock_literature)
        args._test_data_dir = str(tmp_path)
        if args.name is None:
            args.name = "demo"
        return temp_holder, mock_csv

    monkeypatch.setattr(main_module, "_prepare_test_run", fake_prepare)

    updated_csv = tmp_path / "updated.csv"
    updated_csv.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        f"Demo,0,,,,{tmp_path / 'demo.vcf'},{tmp_path / 'demo.cov'}\n",
        encoding="utf-8",
    )

    def fake_workflow(**kwargs):
        recorded["workflow_kwargs"] = kwargs
        return str(updated_csv)

    monkeypatch.setattr(main_module, "run_e2e_workflow", fake_workflow)

    def fake_vcf(args):
        recorded["vcf_args"] = {
            "input_csv": args.input_csv,
            "test_flag": getattr(args, "test", None),
            "suppress": getattr(args, "_suppress_logo", False),
        }
        return 0

    monkeypatch.setattr(main_module, "_run_vcf_command", fake_vcf)

    exit_code = main_module.main(
        [
            "bam",
            "--test",
        ]
    )

    assert exit_code == 0
    assert recorded["prepare_mode"] == "bam"
    assert recorded["workflow_kwargs"]["samples_csv"] == str(mock_csv)
    assert recorded["vcf_args"]["input_csv"] == str(updated_csv)
    assert recorded["vcf_args"]["test_flag"] is False
    assert recorded["vcf_args"]["suppress"] is True
    assert recorded["modes"] == ["bam"]
    assert recorded.get("cleanup_called") is True


def test_e2e_test_mode_uses_demo_dataset(monkeypatch, tmp_path):
    recorded = {}

    def fake_validate(mode="vcf"):
        recorded.setdefault("modes", []).append(mode)

    monkeypatch.setattr(main_module, "validate_dependencies", fake_validate)

    mock_csv = tmp_path / "reads.csv"
    mock_r1 = tmp_path / "reads_1.fq"
    mock_r2 = tmp_path / "reads_2.fq"
    mock_r1.write_text("", encoding="utf-8")
    mock_r2.write_text("", encoding="utf-8")
    mock_literature = tmp_path / "mock.csv"
    mock_literature.write_text("", encoding="utf-8")
    mock_csv.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        f"Demo,0,{mock_r1},{mock_r2},,,\n",
        encoding="utf-8",
    )

    class DummyTemp:
        def __init__(self):
            self.cleaned = False

        def cleanup(self):
            self.cleaned = True
            recorded["cleanup_called"] = True

    temp_holder = DummyTemp()

    def fake_prepare(args, mode):
        recorded["prepare_mode"] = mode
        args.search_pokay = False
        args.literature_csv = str(mock_literature)
        args._test_data_dir = str(tmp_path)
        if args.name is None:
            args.name = "demo"
        return temp_holder, mock_csv

    monkeypatch.setattr(main_module, "_prepare_test_run", fake_prepare)

    updated_csv = tmp_path / "updated.csv"
    updated_csv.write_text(
        "sample_name,sample_number,reads1,reads2,bam,vcf,coverage\n"
        f"Demo,0,,,,{tmp_path / 'demo.vcf'},{tmp_path / 'demo.cov'}\n",
        encoding="utf-8",
    )

    def fake_workflow(**kwargs):
        recorded["workflow_kwargs"] = kwargs
        return str(updated_csv)

    monkeypatch.setattr(main_module, "run_e2e_workflow", fake_workflow)

    def fake_vcf(args):
        recorded["vcf_args"] = {
            "input_csv": args.input_csv,
            "test_flag": getattr(args, "test", None),
            "suppress": getattr(args, "_suppress_logo", False),
        }
        return 0

    monkeypatch.setattr(main_module, "_run_vcf_command", fake_vcf)

    exit_code = main_module.main(
        [
            "end-to-end",
            "--test",
        ]
    )

    assert exit_code == 0
    assert recorded["prepare_mode"] == "e2e"
    assert recorded["workflow_kwargs"]["samples_csv"] == str(mock_csv)
    assert recorded["workflow_kwargs"]["mode"] == "reads"
    assert recorded["vcf_args"]["input_csv"] == str(updated_csv)
    assert recorded["vcf_args"]["test_flag"] is False
    assert recorded["vcf_args"]["suppress"] is True
    assert recorded["modes"] == ["e2e"]
    assert recorded.get("cleanup_called") is True
