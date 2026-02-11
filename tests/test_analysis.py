"""Tests for analysis helpers."""

from __future__ import annotations

import pandas as pd

from vartracker.analysis import (
    search_literature,
    _prepare_variant_heatmap_matrix,
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
    assert 'data-anchor="s:d215g' in content
    assert 'data-anchor="nsp6:l37f' in content
    assert ">10.1000/xyz123</a>" in content
    assert ".cell:hover .cell-value" in content
    assert "clearActive" in content
