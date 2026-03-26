"""Tests for standalone plotting helpers."""

from __future__ import annotations

import pandas as pd

from vartracker.plotting import (
    apply_shared_plot_filters,
    auto_select_variants,
    get_threshold_crossing_variants,
    plot_variant_lifespan,
    plot_variant_trajectory,
    plot_variant_turnover,
    prepare_plot_inputs,
    variant_crosses_thresholds,
)


def _results_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
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
    summary, long_df, sample_names, sample_numbers = prepare_plot_inputs(_results_table())

    assert sample_names == ["P0", "P1", "P2", "P3"]
    assert sample_numbers == [0, 1, 2, 3]
    assert set(summary["variant_id"]) >= {"S:D614G", "S:E484K", "S:N501Y", "N:N202="}

    filtered_summary, filtered_long = apply_shared_plot_filters(
        summary,
        long_df,
        gene="S",
        effects=["missense"],
        min_af=0.5,
        sample_min=1,
        sample_max=3,
    )

    assert set(filtered_summary["variant_id"]) == {"S:D614G", "S:N501Y"}
    assert filtered_long["sample_number"].min() == 1
    assert filtered_long["sample_number"].max() == 3


def test_auto_select_variants_prefers_literature_and_persistent():
    summary, _, _, _ = prepare_plot_inputs(_results_table())

    selected = auto_select_variants(summary, top_n=2)

    assert selected[0] == "S:D614G"
    assert len(selected) == 2


def test_plot_functions_create_output_files(tmp_path, monkeypatch):
    mpl_dir = tmp_path / "mpl"
    mpl_dir.mkdir()
    monkeypatch.setenv("MPLCONFIGDIR", str(mpl_dir))

    summary, long_df, _, _ = prepare_plot_inputs(_results_table())
    selected = ["S:D614G", "S:N501Y"]

    trajectory = tmp_path / "trajectory.pdf"
    turnover = tmp_path / "turnover.pdf"
    lifespan = tmp_path / "lifespan.pdf"
    trajectory_threshold = tmp_path / "trajectory_threshold.pdf"

    plot_variant_trajectory(
        summary,
        long_df,
        selected_variants=selected,
        output_path=trajectory,
    )
    plot_variant_turnover(summary, long_df, output_path=turnover)
    plot_variant_lifespan(
        summary,
        selected_variants=selected,
        output_path=lifespan,
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
    assert variant_crosses_thresholds([0.1, 0.5], [0.5], crossing_rule="strictly_above") is False

    summary, long_df, _, _ = prepare_plot_inputs(_results_table())
    crossed = get_threshold_crossing_variants(summary, long_df, [0.5])
    assert set(crossed["variant_id"]) == {"S:D614G", "S:N501Y"}
