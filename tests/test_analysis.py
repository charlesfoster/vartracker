"""Tests for analysis helpers."""

from __future__ import annotations

import pandas as pd
import pytest

from vartracker.analysis import (
    _heatmap_figure_size,
    search_literature,
    _prepare_variant_heatmap_matrix,
    process_joint_variants,
    generate_variant_heatmap,
)


def test_search_literature_handles_nullable_boolean_masks(tmp_path):
    table = pd.DataFrame(
        {
            "gene": ["S"],
            "amino_acid_consequence": ["S:N501Y"],
            "variant_status": ["new"],
        }
    )

    literature = pd.DataFrame(
        {
            "gene": pd.Series(["S", pd.NA], dtype="string"),
            "mutation": pd.Series(["S:N501Y", ""], dtype="string"),
            "category": ["Functional impact", ""],
            "information": ["Info", ""],
            "reference": ["PMID123", ""],
        }
    )

    result = search_literature(table, literature, tmp_path, "sample")

    assert not result.empty
    assert (tmp_path / "sample.literature_database_hits.full.csv").exists()
    assert (tmp_path / "sample.literature_database_hits.concise.csv").exists()


def test_prepare_variant_heatmap_matrix_orders_variants_by_genome():
    table = pd.DataFrame(
        [
            {
                "gene": "ORF1ab",
                "amino_acid_consequence": "809T",
                "nsp_aa_change": "",
                "type_of_change": "synonymous",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.5",
                "samples": "P0 / P1",
                "variant": "C1059T",
                "start": 1059,
            },
            {
                "gene": "S",
                "amino_acid_consequence": "D215G",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 1.0",
                "samples": "P0 / P1",
                "variant": "A22206G",
                "start": 22206,
            },
            {  # duplicate row should be ignored
                "gene": "S",
                "amino_acid_consequence": "D215G",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 1.0",
                "samples": "P0 / P1",
                "variant": "A22206G",
                "start": 22206,
            },
            {  # should be dropped: below SNV threshold
                "gene": "S",
                "amino_acid_consequence": "A10T",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.1",
                "samples": "P0 / P1",
                "variant": "G21500A",
                "start": 21500,
            },
            {  # should be dropped: indel below threshold
                "gene": "S",
                "amino_acid_consequence": "144del",
                "nsp_aa_change": "",
                "type_of_change": "inframe_deletion",
                "type_of_variant": "indel",
                "alt_freq": "0.05 / 0.05",
                "samples": "P0 / P1",
                "variant": "CT21991C",
                "start": 21991,
            },
            {
                "gene": "5' UTR",
                "amino_acid_consequence": "",
                "nsp_aa_change": "",
                "type_of_change": "None",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.9",
                "samples": "P0 / P1",
                "variant": "C200T",
                "start": 200,
            },
            {
                "gene": "ORF1ab",
                "amino_acid_consequence": "",
                "nsp_aa_change": "nsp2:VLQKAAITILDGISQYSLRLIDAMMFTSDLATNNLVVMAYIVEELKAA",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.8",
                "samples": "P0 / P1",
                "variant": "G1946GT",
                "start": 1946,
            },
        ]
    )

    matrix = _prepare_variant_heatmap_matrix(table, ["P0", "P1"], 0.2, 0.3)

    long_aa = "VLQKAAITILDGISQYSLRLIDAMMFTSDLATNNLVVMAYIVEELKAA"
    head = long_aa[:10]
    tail = long_aa[-10:]
    middle = len(long_aa) - 20

    expected_long_label = f"nsp2:{head}+{middle}{tail}\n(G1946GT)"

    expected_index = [
        "nsp2:T809=\n(C1059T)",
        expected_long_label,
        "S:D215G\n(A22206G)",
    ]

    assert list(matrix.index) == expected_index
    assert matrix.loc["nsp2:T809=\n(C1059T)", "P0"] == 0.0
    assert matrix.loc["nsp2:T809=\n(C1059T)", "P1"] == 0.5
    assert matrix.loc[expected_long_label, "P1"] == 0.8
    assert matrix.loc["S:D215G\n(A22206G)", "P1"] == 1.0


def test_prepare_variant_heatmap_matrix_excludes_selected_consequence_types():
    table = pd.DataFrame(
        [
            {
                "gene": "S",
                "amino_acid_consequence": "D215G",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.7",
                "samples": "P0 / P1",
                "variant": "A22206G",
                "start": 22206,
            },
            {
                "gene": "S",
                "amino_acid_consequence": "T716T",
                "nsp_aa_change": "",
                "type_of_change": "synonymous",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.8",
                "samples": "P0 / P1",
                "variant": "C23403T",
                "start": 23403,
            },
            {
                "gene": "S",
                "amino_acid_consequence": "145del",
                "nsp_aa_change": "",
                "type_of_change": "frameshift",
                "type_of_variant": "indel",
                "alt_freq": "0.0 / 0.9",
                "samples": "P0 / P1",
                "variant": "A22029-",
                "start": 22029,
            },
        ]
    )

    matrix = _prepare_variant_heatmap_matrix(
        table,
        ["P0", "P1"],
        0.2,
        0.2,
        excluded_consequence_types=["synonymous", "frameshift"],
    )

    assert list(matrix.index) == ["S:D215G\n(A22206G)"]


def test_prepare_variant_heatmap_matrix_excludes_wildcard_consequence_types():
    table = pd.DataFrame(
        [
            {
                "gene": "S",
                "amino_acid_consequence": "D215G",
                "nsp_aa_change": "",
                "type_of_change": "joint_frameshift",
                "type_of_variant": "indel",
                "alt_freq": "0.0 / 0.7",
                "samples": "P0 / P1",
                "variant": "A22206G",
                "start": 22206,
            },
            {
                "gene": "S",
                "amino_acid_consequence": "N501Y",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.8",
                "samples": "P0 / P1",
                "variant": "A23063T",
                "start": 23063,
            },
        ]
    )

    matrix = _prepare_variant_heatmap_matrix(
        table,
        ["P0", "P1"],
        0.2,
        0.2,
        excluded_consequence_types=["*frameshift*"],
    )

    assert list(matrix.index) == ["S:N501Y\n(A23063T)"]


def test_prepare_variant_heatmap_matrix_applies_extended_filters():
    table = pd.DataFrame(
        [
            {
                "gene": "S",
                "amino_acid_consequence": "D215G",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "variant_status": "new",
                "persistence_status": "new_persistent",
                "all_samples_pass_qc": True,
                "proportion_samples_passing_qc": 1.0,
                "presence_absence": "N / Y / Y",
                "alt_freq": "0.0 / 0.6 / 0.7",
                "variant_site_depth": "0 / 120 / 125",
                "samples": "P0 / P1 / P2",
                "variant": "A22206G",
                "start": 22206,
            },
            {
                "gene": "S",
                "amino_acid_consequence": "145fs",
                "nsp_aa_change": "",
                "type_of_change": "joint_frameshift",
                "type_of_variant": "indel",
                "variant_status": "new",
                "persistence_status": "new_persistent",
                "all_samples_pass_qc": True,
                "proportion_samples_passing_qc": 0.33,
                "presence_absence": "N / N / Y",
                "alt_freq": "0.0 / 0.0 / 0.8",
                "variant_site_depth": "0 / 0 / 130",
                "samples": "P0 / P1 / P2",
                "variant": "A22029-",
                "start": 22029,
            },
            {
                "gene": "N",
                "amino_acid_consequence": "R203K",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "variant_status": "original",
                "persistence_status": "original_retained",
                "all_samples_pass_qc": False,
                "proportion_samples_passing_qc": 0.67,
                "presence_absence": "Y / Y / Y",
                "alt_freq": "0.5 / 0.5 / 0.5",
                "variant_site_depth": "110 / 115 / 120",
                "samples": "P0 / P1 / P2",
                "variant": "G28881A",
                "start": 28881,
            },
        ]
    )

    matrix = _prepare_variant_heatmap_matrix(
        table,
        ["P0", "P1", "P2"],
        0.2,
        0.2,
        only_new=True,
        gene_include=["S"],
        variant_type_include=["snp"],
        qc_include=["true"],
        min_prop_passing_qc=0.9,
        min_persistence=2,
        min_max_af=0.6,
        sample_subset=["P1", "P2"],
        hide_singletons=True,
        min_depth=100,
    )

    assert list(matrix.index) == ["S:D215G\n(A22206G)"]


def test_process_joint_variants_only_adds_single_joint_prefix(tmp_path):
    csv_path = tmp_path / "results.csv"
    pd.DataFrame(
        [
            {
                "start": 100,
                "gene": "S",
                "amino_acid_consequence": "N501Y",
                "nsp_aa_change": "",
                "bcsq_nt_notation": "c.1A>T",
                "bcsq_aa_notation": "p.N501Y",
                "aa1_total_properties": "",
                "aa2_total_properties": "",
                "aa1_unique_properties": "",
                "aa2_unique_properties": "",
                "aa1_weight": "",
                "aa2_weight": "",
                "weight_difference": "",
                "type_of_change": "joint_joint_frameshift",
            },
            {
                "start": 101,
                "gene": "",
                "amino_acid_consequence": "",
                "nsp_aa_change": "",
                "bcsq_nt_notation": "",
                "bcsq_aa_notation": "@100",
                "aa1_total_properties": "",
                "aa2_total_properties": "",
                "aa1_unique_properties": "",
                "aa2_unique_properties": "",
                "aa1_weight": "",
                "aa2_weight": "",
                "weight_difference": "",
                "type_of_change": "frameshift",
            },
        ]
    ).to_csv(csv_path, index=False)

    result = process_joint_variants(str(csv_path))

    assert result.loc[0, "type_of_change"] == "joint_frameshift"
    assert result.loc[1, "type_of_change"] == "joint_frameshift"


def test_process_joint_variants_matches_main_row_by_presence_pattern(tmp_path):
    csv_path = tmp_path / "results.csv"
    pd.DataFrame(
        [
            {
                "start": 100,
                "gene": "S",
                "amino_acid_consequence": "K2A",
                "nsp_aa_change": "",
                "bcsq_nt_notation": "4A>G+5A>C",
                "bcsq_aa_notation": "2K>2A",
                "aa1_total_properties": "",
                "aa2_total_properties": "",
                "aa1_unique_properties": "",
                "aa2_unique_properties": "",
                "aa1_weight": "",
                "aa2_weight": "",
                "weight_difference": "",
                "type_of_change": "missense",
                "presence_absence": "Y / N",
            },
            {
                "start": 100,
                "gene": "S",
                "amino_acid_consequence": "K2V",
                "nsp_aa_change": "",
                "bcsq_nt_notation": "4A>T+5A>G",
                "bcsq_aa_notation": "2K>2V",
                "aa1_total_properties": "",
                "aa2_total_properties": "",
                "aa1_unique_properties": "",
                "aa2_unique_properties": "",
                "aa1_weight": "",
                "aa2_weight": "",
                "weight_difference": "",
                "type_of_change": "missense",
                "presence_absence": "N / Y",
            },
            {
                "start": 101,
                "gene": "",
                "amino_acid_consequence": "",
                "nsp_aa_change": "",
                "bcsq_nt_notation": "",
                "bcsq_aa_notation": "@100",
                "aa1_total_properties": "",
                "aa2_total_properties": "",
                "aa1_unique_properties": "",
                "aa2_unique_properties": "",
                "aa1_weight": "",
                "aa2_weight": "",
                "weight_difference": "",
                "type_of_change": "@100",
                "presence_absence": "N / Y",
            },
        ]
    ).to_csv(csv_path, index=False)

    result = process_joint_variants(str(csv_path))

    assert result.loc[2, "amino_acid_consequence"] == "K2V"
    assert result.loc[2, "type_of_change"] == "joint_missense"


@pytest.mark.xfail(
    strict=True,
    reason=(
        "process_joint_variants currently resolves tied main-row candidates by row "
        "order, which is unsafe for overlapping genes"
    ),
)
def test_process_joint_variants_is_order_invariant_for_overlapping_gene_rows(tmp_path):
    shared_rows = [
        {
            "start": 25470,
            "gene": "ORF3a",
            "amino_acid_consequence": "ORF3a:A10V",
            "nsp_aa_change": "",
            "bcsq_nt_notation": "c.30C>T",
            "bcsq_aa_notation": "p.A10V",
            "aa1_total_properties": "",
            "aa2_total_properties": "",
            "aa1_unique_properties": "",
            "aa2_unique_properties": "",
            "aa1_weight": "",
            "aa2_weight": "",
            "weight_difference": "",
            "type_of_change": "missense",
            "presence_absence": "Y / N / Y",
        },
        {
            "start": 25470,
            "gene": "ORF3c",
            "amino_acid_consequence": "ORF3c:M5I",
            "nsp_aa_change": "",
            "bcsq_nt_notation": "c.15G>A",
            "bcsq_aa_notation": "p.M5I",
            "aa1_total_properties": "",
            "aa2_total_properties": "",
            "aa1_unique_properties": "",
            "aa2_unique_properties": "",
            "aa1_weight": "",
            "aa2_weight": "",
            "weight_difference": "",
            "type_of_change": "missense",
            "presence_absence": "Y / N / Y",
        },
        {
            "start": 25471,
            "gene": "",
            "amino_acid_consequence": "",
            "nsp_aa_change": "",
            "bcsq_nt_notation": "",
            "bcsq_aa_notation": "@25470",
            "aa1_total_properties": "",
            "aa2_total_properties": "",
            "aa1_unique_properties": "",
            "aa2_unique_properties": "",
            "aa1_weight": "",
            "aa2_weight": "",
            "weight_difference": "",
            "type_of_change": "@25470",
            "presence_absence": "Y / N / Y",
        },
    ]

    first_csv = tmp_path / "overlap_first.csv"
    second_csv = tmp_path / "overlap_second.csv"
    pd.DataFrame(shared_rows).to_csv(first_csv, index=False)
    pd.DataFrame([shared_rows[1], shared_rows[0], shared_rows[2]]).to_csv(
        second_csv, index=False
    )

    first_result = process_joint_variants(str(first_csv))
    second_result = process_joint_variants(str(second_csv))

    first_joint = first_result.loc[2, ["gene", "amino_acid_consequence"]].to_dict()
    second_joint = second_result.loc[2, ["gene", "amino_acid_consequence"]].to_dict()

    assert first_joint == second_joint


def test_generate_variant_heatmap_creates_interactive_html(tmp_path, monkeypatch):
    mpl_dir = tmp_path / "mpl"
    mpl_dir.mkdir()
    monkeypatch.setenv("MPLCONFIGDIR", str(mpl_dir))

    table = pd.DataFrame(
        [
            {
                "gene": "S",
                "amino_acid_consequence": "D215G",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.5",
                "per_sample_variant_qc": "P / F",
                "samples": "P0 / P1",
                "variant": "A22206G",
                "start": 22206,
                "variant_status": "new",
                "presence_absence": "N / Y",
            },
            {
                "gene": "ORF1ab",
                "amino_acid_consequence": "",
                "nsp_aa_change": "nsp6_TM:L37F",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.2 / 0.6",
                "per_sample_variant_qc": "P / P",
                "samples": "P0 / P1",
                "variant": "C21575T",
                "start": 21575,
                "variant_status": "new",
                "presence_absence": "N / Y",
            },
        ]
    )

    literature_hits = pd.DataFrame(
        {
            "gene": ["S", "nsp6"],
            "amino_acid_consequence": ["D215G", "L37F"],
            "database_mutation_string": ["S:D215G", "nsp6:L37F"],
            "category": ["Functional impact", "Functional impact"],
            "prior_information": ["Info", "Additional"],
            "reference": ["PMID123; 10.1000/xyz123", "10.2000/abc456"],
        }
    )

    literature_csv = tmp_path / "sample.literature_database_hits.full.csv"
    literature_hits.to_csv(literature_csv, index=False)

    command_str = "vartracker vcf input.csv --outdir results"

    generate_variant_heatmap(
        table,
        ["P0", "P1"],
        [0, 1],
        str(tmp_path),
        "Example",
        0.05,
        0.05,
        literature_hits=literature_hits,
        literature_table_path=str(literature_csv),
        cli_command=command_str,
    )

    pdf_path = tmp_path / "variant_allele_frequency_heatmap.pdf"
    html_path = tmp_path / "variant_allele_frequency_heatmap.html"

    assert pdf_path.exists()
    assert html_path.exists()

    content = html_path.read_text(encoding="utf-8")
    assert "Interactive variant heatmap" in content
    assert "Workflow summary" in content
    assert command_str in content
    assert "heatmap-grid" in content
    assert "heatmap-scroll" in content
    assert "Literature results" in content
    assert "table-scroll" in content
    assert "heatmap-anchor" in content
    assert "cell-qc-fail" in content
    assert "AF=0.50, QC=FAIL" in content
    assert 'data-anchor="s:d215g' in content
    assert 'data-anchor="nsp6:l37f' in content
    assert ">10.1000/xyz123</a>" in content
    assert ".cell:hover .cell-value" in content
    assert "clearActive" in content


def test_heatmap_figure_size_enforces_minimum_row_height():
    width, height = _heatmap_figure_size(6, 2)

    assert width == 4.0
    assert height == 4.0

    width, height = _heatmap_figure_size(8, 2)

    assert width == 4.0
    assert height == pytest.approx(5.2)
