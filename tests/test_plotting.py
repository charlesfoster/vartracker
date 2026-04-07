"""Tests for standalone plotting helpers."""

from __future__ import annotations

import json

import pandas as pd
import pytest

from vartracker.core import InputValidationError
from vartracker.plotting import (
    _collapse_variants_for_genome_plot,
    _genome_axis_label,
    _parse_focus_ranges,
    _subset_collapsed_variants,
    build_reference_feature_metadata,
    apply_shared_plot_filters,
    auto_select_variants,
    get_threshold_crossing_variants,
    plot_variant_genome,
    plot_variant_lifespan,
    plot_variant_trajectory,
    plot_variant_turnover,
    prepare_plot_inputs,
    resolve_focus_regions,
    variant_crosses_thresholds,
)


def _results_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "chrom": "segA",
                "start": 23403,
                "end": 23403,
                "gene": "S",
                "variant": "A23403G",
                "amino_acid_consequence": "S:D614G",
                "nsp_aa_change": "",
                "type_of_variant": "snp",
                "type_of_change": "missense",
                "variant_status": "new",
                "persistence_status": "new_persistent",
                "presence_absence": "N / Y / Y / Y",
                "alt_freq": "0.0 / 0.20 / 0.45 / 0.70",
                "samples": "P0 / P1 / P2 / P3",
                "sample_number": "0 / 1 / 2 / 3",
                "per_sample_variant_qc": "P / P / P / P",
                "reference": "PMID123",
            },
            {
                "chrom": "segA",
                "start": 23012,
                "end": 23012,
                "gene": "S",
                "variant": "G23012A",
                "amino_acid_consequence": "S:E484K",
                "nsp_aa_change": "",
                "type_of_variant": "snp",
                "type_of_change": "missense",
                "variant_status": "new",
                "persistence_status": "new_transient",
                "presence_absence": "N / N / Y / N",
                "alt_freq": "0.0 / 0.0 / 0.40 / 0.0",
                "samples": "P0 / P1 / P2 / P3",
                "sample_number": "0 / 1 / 2 / 3",
                "per_sample_variant_qc": "P / P / P / P",
                "reference": "",
            },
            {
                "chrom": "segA",
                "start": 23063,
                "end": 23063,
                "gene": "S",
                "variant": "A23063T",
                "amino_acid_consequence": "S:N501Y",
                "nsp_aa_change": "",
                "type_of_variant": "snp",
                "type_of_change": "missense",
                "variant_status": "original",
                "persistence_status": "original_retained",
                "presence_absence": "Y / Y / Y / Y",
                "alt_freq": "0.55 / 0.58 / 0.60 / 0.62",
                "samples": "P0 / P1 / P2 / P3",
                "sample_number": "0 / 1 / 2 / 3",
                "per_sample_variant_qc": "P / P / P / P",
                "reference": "",
            },
            {
                "chrom": "segB",
                "start": 28977,
                "end": 28977,
                "gene": "N",
                "variant": "C28977T",
                "amino_acid_consequence": "N:S202=",
                "nsp_aa_change": "",
                "type_of_variant": "snp",
                "type_of_change": "synonymous",
                "variant_status": "new",
                "persistence_status": "new_persistent",
                "presence_absence": "N / Y / Y / Y",
                "alt_freq": "0.0 / 0.30 / 0.35 / 0.40",
                "samples": "P0 / P1 / P2 / P3",
                "sample_number": "0 / 1 / 2 / 3",
                "per_sample_variant_qc": "P / P / P / P",
                "reference": "",
            },
        ]
    )


def test_prepare_plot_inputs_and_filters():
    summary, long_df, sample_names, sample_numbers = prepare_plot_inputs(
        _results_table()
    )

    assert sample_names == ["P0", "P1", "P2", "P3"]
    assert sample_numbers == [0, 1, 2, 3]
    assert {label.split(" (", 1)[0] for label in summary["variant_label"]} >= {
        "S:D614G",
        "S:E484K",
        "S:N501Y",
        "N:N202=",
    }

    filtered_summary, filtered_long = apply_shared_plot_filters(
        summary,
        long_df,
        gene="S",
        effects=["missense"],
        min_af=0.5,
        sample_min=1,
        sample_max=3,
    )

    assert {label.split(" (", 1)[0] for label in filtered_summary["variant_label"]} == {
        "S:D614G",
        "S:N501Y",
    }
    assert filtered_long["sample_number"].min() == 1
    assert filtered_long["sample_number"].max() == 3


def test_auto_select_variants_prefers_literature_and_persistent():
    summary, _, _, _ = prepare_plot_inputs(_results_table())

    selected = auto_select_variants(summary, top_n=2)

    first_label = summary.loc[
        summary["variant_id"] == selected[0], "variant_label"
    ].iloc[0]
    assert first_label.split(" (", 1)[0] == "S:D614G"
    assert len(selected) == 2


def test_genome_axis_label_uses_gene_name_for_cds_scale():
    assert (
        _genome_axis_label(gene="F", aa_scale=False, cds_scale=True)
        == "CDS position in F"
    )


def test_plot_functions_create_output_files(tmp_path, monkeypatch):
    mpl_dir = tmp_path / "mpl"
    mpl_dir.mkdir()
    monkeypatch.setenv("MPLCONFIGDIR", str(mpl_dir))

    summary, long_df, _, _ = prepare_plot_inputs(_results_table())
    selected = summary.loc[
        summary["variant_label"]
        .str.replace(r" \(.*\)$", "", regex=True)
        .isin(["S:D614G", "S:N501Y"]),
        "variant_id",
    ].tolist()

    trajectory = tmp_path / "trajectory.pdf"
    turnover = tmp_path / "turnover.pdf"
    lifespan = tmp_path / "lifespan.pdf"
    trajectory_threshold = tmp_path / "trajectory_threshold.pdf"

    plot_variant_trajectory(
        summary,
        long_df,
        selected_variants=selected,
        output_path=trajectory,
        sample_axis_mode="name",
    )
    plot_variant_turnover(
        summary, long_df, output_path=turnover, sample_axis_mode="name"
    )
    plot_variant_lifespan(
        summary,
        selected_variants=selected,
        output_path=lifespan,
        sample_axis_mode="name",
        sample_names=["P0", "P1", "P2", "P3"],
        sample_numbers=[0, 1, 2, 3],
    )
    plot_variant_trajectory(
        summary,
        long_df,
        selected_variants=selected,
        output_path=trajectory_threshold,
        thresholds=[0.5, 0.9],
    )

    assert trajectory.exists()
    assert turnover.exists()
    assert lifespan.exists()
    assert trajectory_threshold.exists()


def test_threshold_crossing_helpers():
    assert variant_crosses_thresholds([0.1, 0.5], [0.5]) is True
    assert (
        variant_crosses_thresholds([0.1, 0.5], [0.5], crossing_rule="strictly_above")
        is False
    )

    summary, long_df, _, _ = prepare_plot_inputs(_results_table())
    crossed = get_threshold_crossing_variants(summary, long_df, [0.5])
    assert {label.split(" (", 1)[0] for label in crossed["variant_label"]} == {
        "S:D614G",
        "S:N501Y",
    }


def test_genome_plot_collapses_variants_by_max_af():
    table = pd.concat([_results_table(), _results_table().iloc[[0]]], ignore_index=True)

    collapsed = _collapse_variants_for_genome_plot(table)

    assert len(collapsed) == 4
    d614g = collapsed.loc[
        collapsed["variant_label"].str.replace(r" \(.*\)$", "", regex=True) == "S:D614G"
    ].iloc[0]
    assert d614g["summary_af"] == 0.7
    assert d614g["af_values"] == [0.0, 0.2, 0.45, 0.7]
    assert d614g["chrom"] == "segA"
    assert d614g["start"] == 23403


def test_prepare_plot_inputs_uses_unique_internal_variant_keys():
    table = pd.DataFrame(
        [
            {
                "chrom": "segA",
                "start": 100,
                "end": 100,
                "ref": "A",
                "alt": "G",
                "gene": "INTERGENIC",
                "variant": "A100G",
                "amino_acid_consequence": "None",
                "nsp_aa_change": "None",
                "type_of_variant": "snp",
                "type_of_change": "intergenic",
                "variant_status": "new",
                "persistence_status": "new_persistent",
                "presence_absence": "Y / N",
                "alt_freq": "0.8 / 0.0",
                "samples": "P1 / P2",
                "sample_number": "1 / 2",
                "per_sample_variant_qc": "P / P",
            },
            {
                "chrom": "segA",
                "start": 200,
                "end": 200,
                "ref": "C",
                "alt": "T",
                "gene": "INTERGENIC",
                "variant": "C200T",
                "amino_acid_consequence": "None",
                "nsp_aa_change": "None",
                "type_of_variant": "snp",
                "type_of_change": "intergenic",
                "variant_status": "new",
                "persistence_status": "new_persistent",
                "presence_absence": "N / Y",
                "alt_freq": "0.0 / 0.7",
                "samples": "P1 / P2",
                "sample_number": "1 / 2",
                "per_sample_variant_qc": "P / P",
            },
        ]
    )

    summary, long_df, _, _ = prepare_plot_inputs(table)

    assert summary["variant_id"].nunique() == 2
    assert {label.split(" (", 1)[0] for label in summary["variant_label"]} == {
        "INTERGENIC:None"
    }
    assert (
        long_df.loc[
            long_df["variant_label"].str.startswith("INTERGENIC:None"), "variant_id"
        ].nunique()
        == 2
    )


def test_genome_plot_defaults_to_snps_only(capsys):
    table = pd.DataFrame(
        [
            {
                "chrom": "segA",
                "start": 100,
                "end": 100,
                "gene": "S",
                "variant": "A100G",
                "amino_acid_consequence": "S:V1A",
                "nsp_aa_change": "",
                "type_of_variant": "snp",
                "type_of_change": "missense",
                "variant_status": "new",
                "persistence_status": "new_persistent",
                "presence_absence": "Y / Y",
                "alt_freq": "0.1 / 0.4",
            },
            {
                "chrom": "segA",
                "start": 101,
                "end": 101,
                "gene": "S",
                "variant": "A101AG",
                "amino_acid_consequence": "S:V2del",
                "nsp_aa_change": "",
                "type_of_variant": "indel",
                "type_of_change": "frameshift",
                "variant_status": "new",
                "persistence_status": "new_persistent",
                "presence_absence": "Y / Y",
                "alt_freq": "0.2 / 0.5",
            },
        ]
    )

    collapsed = _subset_collapsed_variants(_collapse_variants_for_genome_plot(table))

    assert list(collapsed["type_of_variant"]) == ["snp"]
    assert capsys.readouterr().out == ""


def test_genome_plot_builds_outputs_and_supports_gene_modes(tmp_path, monkeypatch):
    mpl_dir = tmp_path / "mpl"
    mpl_dir.mkdir()
    monkeypatch.setenv("MPLCONFIGDIR", str(mpl_dir))

    gff = tmp_path / "reference.gff3"
    gff.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "##sequence-region segA 1 30000",
                "##sequence-region segB 1 15000",
                "segA\tRefSeq\tgene\t21563\t25384\t.\t+\t.\tID=gene:S;Name=S",
                "segA\tRefSeq\tmRNA\t21563\t25384\t.\t+\t.\tID=transcript:S;Name=S",
                "segB\tRefSeq\tgene\t28274\t29533\t.\t+\t.\tID=gene:N;Name=N",
                "segB\tRefSeq\tmRNA\t28274\t29533\t.\t+\t.\tID=transcript:N;Name=N",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    metadata = build_reference_feature_metadata(gff)
    table = _results_table()

    genome_out = tmp_path / "genome.pdf"
    focus_out = tmp_path / "genome_focus.pdf"
    aa_out = tmp_path / "genome_gene_aa.pdf"
    cds_out = tmp_path / "genome_gene_cds.pdf"

    plot_variant_genome(table, metadata, output_path=genome_out)
    plot_variant_genome(
        table,
        metadata,
        output_path=focus_out,
        focus_ranges=resolve_focus_regions(
            "Fusion peptide:23000-24000;Heptad repeat:50000-51000", ""
        )[0],
        focus_labels=resolve_focus_regions(
            "Fusion peptide:23000-24000;Heptad repeat:50000-51000", ""
        )[1],
        gene="S",
        show_intersections=True,
    )
    plot_variant_genome(
        table,
        metadata,
        output_path=aa_out,
        gene="S",
        aa_scale=True,
        focus_ranges=_parse_focus_ranges("450-650"),
    )
    plot_variant_genome(
        table,
        metadata,
        output_path=cds_out,
        gene="S",
        cds_scale=True,
        focus_ranges=_parse_focus_ranges("440-500"),
    )

    assert genome_out.exists()
    assert focus_out.exists()
    assert aa_out.exists()
    assert cds_out.exists()


def test_parse_focus_ranges_caps_number_of_regions():
    with pytest.raises(InputValidationError):
        _parse_focus_ranges("1-2,3-4,5-6,7-8,9-10,11-12,13-14")


def test_parse_focus_ranges_supports_semicolon_color_groups():
    parsed = _parse_focus_ranges(
        "62-69,196-210;31-42,323-332,379-399;254-277;50-50,305-310;422-438;163-181"
    )

    assert parsed == [
        (62.0, 69.0, 0),
        (196.0, 210.0, 0),
        (31.0, 42.0, 1),
        (323.0, 332.0, 1),
        (379.0, 399.0, 1),
        (254.0, 277.0, 2),
        (50.0, 50.0, 3),
        (305.0, 310.0, 3),
        (422.0, 438.0, 4),
        (163.0, 181.0, 5),
    ]


def test_parse_focus_ranges_caps_number_of_groups():
    with pytest.raises(InputValidationError):
        _parse_focus_ranges("1-2;3-4;5-6;7-8;9-10;11-12;13-14")


def test_resolve_focus_regions_supports_named_inline_groups():
    ranges, labels = resolve_focus_regions(
        "Ø:62-69,196-210;I:31-42,323-332,379-399;II:254-277", ""
    )

    assert labels == ["Ø", "I", "II"]
    assert ranges == [
        (62.0, 69.0, 0),
        (196.0, 210.0, 0),
        (31.0, 42.0, 1),
        (323.0, 332.0, 1),
        (379.0, 399.0, 1),
        (254.0, 277.0, 2),
    ]


def test_resolve_focus_regions_supports_json_and_csv_files(tmp_path):
    json_path = tmp_path / "regions.json"
    json_path.write_text(
        json.dumps({"Ø": ["62-69", "196-210"], "I": ["31-42", "323-332"]}),
        encoding="utf-8",
    )
    csv_path = tmp_path / "regions.csv"
    csv_path.write_text(
        "label,start,end\nII,254,277\nIII,50,50\nIII,305,310\n",
        encoding="utf-8",
    )

    json_ranges, json_labels = resolve_focus_regions("", str(json_path))
    csv_ranges, csv_labels = resolve_focus_regions("", str(csv_path))

    assert json_labels == ["Ø", "I"]
    assert json_ranges == [
        (62.0, 69.0, 0),
        (196.0, 210.0, 0),
        (31.0, 42.0, 1),
        (323.0, 332.0, 1),
    ]
    assert csv_labels == ["II", "III"]
    assert csv_ranges == [
        (254.0, 277.0, 0),
        (50.0, 50.0, 1),
        (305.0, 310.0, 1),
    ]
