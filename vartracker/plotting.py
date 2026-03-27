"""Standalone plotting helpers for vartracker results tables."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Sequence, cast
from urllib.parse import unquote

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

from .analysis import (
    _coerce_frequency,
    _extract_numeric_position,
    _resolve_variant_labels,
)
from .core import InputValidationError, ProcessingError

DEFAULT_TRAJECTORY_TOP_N = 12
DEFAULT_LIFESPAN_TOP_N = 20
DEFAULT_TRAJECTORY_THRESHOLDS = (0.5, 0.9)
MAX_FOCUS_RANGES = 6
FOCUS_RANGE_COLORS = [
    "#f8bbd0",
    "#c5cae9",
    "#c8e6c9",
    "#ffe0b2",
    "#d1c4e9",
    "#b2dfdb",
]
FeatureRecord = dict[str, object]


def _parse_csv_option_list(value: str | None) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _read_variant_list_file(path: str | None) -> list[str]:
    if not path:
        return []
    variant_path = Path(path).expanduser().resolve()
    if not variant_path.exists():
        raise InputValidationError(f"Variant file not found: {variant_path}")
    lines = [
        line.strip()
        for line in variant_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    return lines


def load_results_table(results_csv: str | Path) -> pd.DataFrame:
    path = Path(results_csv).expanduser().resolve()
    if not path.exists():
        raise InputValidationError(f"Results CSV not found: {path}")
    table = pd.read_csv(path, keep_default_na=False)
    if table.empty:
        raise InputValidationError("Results CSV is empty")
    return table


def _parse_slash_tokens(value: object) -> list[str]:
    if value is None:
        return []
    return [token.strip() for token in str(value).split(" / ") if token.strip()]


def _get_sample_axis(table: pd.DataFrame) -> tuple[list[str], list[int]]:
    if "samples" not in table.columns:
        raise InputValidationError(
            "Results CSV must contain a 'samples' column for standalone plotting"
        )
    if "sample_number" not in table.columns:
        raise InputValidationError(
            "Results CSV must contain a 'sample_number' column for standalone plotting"
        )

    sample_names = _parse_slash_tokens(table.iloc[0]["samples"])
    sample_numbers_raw = _parse_slash_tokens(table.iloc[0]["sample_number"])
    if not sample_names or not sample_numbers_raw:
        raise InputValidationError("Could not determine sample axis from results CSV")
    if len(sample_names) != len(sample_numbers_raw):
        raise InputValidationError(
            "Results CSV has mismatched 'samples' and 'sample_number' columns"
        )
    try:
        sample_numbers = [int(float(value)) for value in sample_numbers_raw]
    except ValueError as exc:
        raise InputValidationError(
            "Results CSV 'sample_number' values must be numeric"
        ) from exc
    return sample_names, sample_numbers


def _project_name_from_results(table: pd.DataFrame, explicit_name: str | None) -> str:
    if explicit_name is not None:
        return explicit_name
    if "name" not in table.columns:
        return ""
    names = [
        str(value).strip() for value in table["name"].unique() if str(value).strip()
    ]
    return names[0] if len(names) == 1 else ""


def _build_variant_key(row: object, fallback_label: str) -> str:
    chrom = str(getattr(row, "chrom", "")).strip()
    start = str(getattr(row, "start", "")).strip()
    ref = str(getattr(row, "ref", "")).strip()
    alt = str(getattr(row, "alt", "")).strip()
    variant_name = str(getattr(row, "variant", "")).strip()

    if chrom and start and ref and alt:
        return f"{chrom}:{start}:{ref}>{alt}"
    if chrom and start and variant_name:
        return f"{chrom}:{start}:{variant_name}"
    if variant_name:
        return variant_name
    return fallback_label


def _display_label_prefix(label: object) -> str:
    return str(label).split(" (", 1)[0].strip()


def _sample_qc_flags(row: object, length: int) -> list[str]:
    flags = _parse_slash_tokens(getattr(row, "per_sample_variant_qc", ""))
    if not flags:
        return [""] * length
    if len(flags) < length:
        flags.extend([""] * (length - len(flags)))
    return [str(flag).strip().upper() for flag in flags[:length]]


def prepare_plot_inputs(
    table: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[int]]:
    sample_names, sample_numbers = _get_sample_axis(table)
    long_records: list[dict[str, object]] = []

    for row in table.drop_duplicates().itertuples(index=False):
        _, display_label, base_label = _resolve_variant_labels(row)
        variant_key = _build_variant_key(row, base_label)
        names = _parse_slash_tokens(getattr(row, "samples", ""))
        numbers = _parse_slash_tokens(getattr(row, "sample_number", ""))
        freqs = [
            _coerce_frequency(token)
            for token in _parse_slash_tokens(getattr(row, "alt_freq", ""))
        ]
        presence = _parse_slash_tokens(getattr(row, "presence_absence", ""))
        qc_flags = _sample_qc_flags(row, len(names))

        length = (
            min(len(names), len(numbers), len(freqs))
            if freqs
            else min(len(names), len(numbers))
        )
        if length == 0:
            continue

        if len(freqs) < length:
            freqs.extend([0.0] * (length - len(freqs)))
        if len(presence) < length:
            presence.extend(["N"] * (length - len(presence)))

        for sample_name, sample_number, af, present, qc_flag in zip(
            names[:length],
            numbers[:length],
            freqs[:length],
            presence[:length],
            qc_flags[:length],
        ):
            long_records.append(
                {
                    "variant_id": variant_key,
                    "variant_label": display_label.replace("\n", " "),
                    "variant_name": str(getattr(row, "variant", "")).strip(),
                    "gene": str(getattr(row, "gene", "")).strip(),
                    "type_of_change": str(getattr(row, "type_of_change", "")).strip(),
                    "type_of_variant": str(getattr(row, "type_of_variant", "")).strip(),
                    "variant_status": str(getattr(row, "variant_status", "")).strip(),
                    "persistence_status": str(
                        getattr(row, "persistence_status", "")
                    ).strip(),
                    "sample_name": sample_name,
                    "sample_number": int(float(sample_number)),
                    "allele_frequency": float(af),
                    "present": str(present).strip().upper() == "Y" or float(af) > 0,
                    "sample_qc": qc_flag,
                    "sample_pass_qc": qc_flag != "F",
                    "has_literature": any(
                        str(getattr(row, column, "")).strip()
                        not in {"", "False", "false", "0", "None"}
                        for column in (
                            "key_mutation",
                            "category",
                            "database_mutation_string",
                            "prior_information",
                            "reference",
                        )
                        if hasattr(row, column)
                    ),
                }
            )

    if not long_records:
        raise ProcessingError("No plottable variant records found in results CSV")

    long_df = pd.DataFrame(long_records).drop_duplicates()
    long_df = long_df.sort_values(
        ["sample_number", "variant_id", "sample_name"]
    ).reset_index(drop=True)

    summary = long_df.groupby("variant_id", as_index=False).agg(
        variant_label=("variant_label", "first"),
        variant_name=("variant_name", "first"),
        gene=("gene", "first"),
        type_of_change=("type_of_change", "first"),
        type_of_variant=("type_of_variant", "first"),
        variant_status=("variant_status", "first"),
        persistence_status=("persistence_status", "first"),
        has_literature=("has_literature", "max"),
        max_af=("allele_frequency", "max"),
    )

    present_df = long_df[long_df["present"]]
    first_seen = (
        present_df.groupby("variant_id")["sample_number"].min().rename("first_seen")
    )
    last_seen = (
        present_df.groupby("variant_id")["sample_number"].max().rename("last_seen")
    )
    summary = summary.merge(first_seen, on="variant_id", how="left")
    summary = summary.merge(last_seen, on="variant_id", how="left")
    summary["first_seen"] = summary["first_seen"].fillna(
        summary["max_af"].map(lambda _x: np.nan)
    )
    summary["last_seen"] = summary["last_seen"].fillna(summary["first_seen"])
    summary["duration"] = summary["last_seen"].fillna(0).astype(float) - summary[
        "first_seen"
    ].fillna(0).astype(float)
    summary["is_persistent_new"] = summary["persistence_status"].eq("new_persistent")
    summary["is_nonsynonymous"] = (
        ~summary["type_of_change"].str.lower().str.contains("synonymous", na=False)
    )

    return summary, long_df, sample_names, sample_numbers


def apply_shared_plot_filters(
    summary: pd.DataFrame,
    long_df: pd.DataFrame,
    *,
    gene: str | None = None,
    effects: Sequence[str] | None = None,
    min_af: float | None = None,
    max_af: float | None = None,
    variants: Sequence[str] | None = None,
    sample_min: int | None = None,
    sample_max: int | None = None,
    include_synonymous: bool = True,
    persistent_only: bool = False,
    new_only: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    filtered_long = long_df.copy()
    if sample_min is not None:
        filtered_long = filtered_long[filtered_long["sample_number"] >= sample_min]
    if sample_max is not None:
        filtered_long = filtered_long[filtered_long["sample_number"] <= sample_max]
    if filtered_long.empty:
        raise ProcessingError(
            "No variant observations remain after sample-range filtering"
        )

    filtered_summary = summary.copy()

    if gene:
        filtered_summary = filtered_summary[
            filtered_summary["gene"].astype(str).str.lower()
            == str(gene).strip().lower()
        ]

    if effects:
        lowered = {
            str(effect).strip().lower() for effect in effects if str(effect).strip()
        }
        filtered_summary = filtered_summary[
            filtered_summary["type_of_change"].astype(str).str.lower().isin(lowered)
        ]

    if not include_synonymous:
        filtered_summary = filtered_summary[
            ~filtered_summary["type_of_change"]
            .astype(str)
            .str.lower()
            .str.contains("synonymous", na=False)
        ]

    if persistent_only:
        filtered_summary = filtered_summary[
            filtered_summary["persistence_status"].eq("new_persistent")
        ]

    if new_only:
        filtered_summary = filtered_summary[
            filtered_summary["variant_status"].eq("new")
        ]

    if variants:
        order_map = {name: idx for idx, name in enumerate(variants)}
        filtered_summary = filtered_summary[
            filtered_summary.apply(
                lambda row: any(
                    candidate in order_map
                    for candidate in (
                        row["variant_id"],
                        row["variant_label"],
                        _display_label_prefix(row["variant_label"]),
                        row["variant_name"],
                    )
                ),
                axis=1,
            )
        ]
        if filtered_summary.empty:
            raise ProcessingError("No variants remain after applying plot filters")
        filtered_summary["selection_order"] = filtered_summary.apply(
            lambda row: min(
                order_map[candidate]
                for candidate in (
                    row["variant_id"],
                    row["variant_label"],
                    _display_label_prefix(row["variant_label"]),
                    row["variant_name"],
                )
                if candidate in order_map
            ),
            axis=1,
        )
        filtered_summary = filtered_summary.sort_values("selection_order").drop(
            columns=["selection_order"]
        )

    if filtered_summary.empty:
        raise ProcessingError("No variants remain after applying plot filters")

    filtered_long = filtered_long[
        filtered_long["variant_id"].isin(filtered_summary["variant_id"])
    ]
    recomputed = filtered_long.groupby("variant_id", as_index=False).agg(
        max_af=("allele_frequency", "max"),
    )
    filtered_summary = filtered_summary.drop(columns=["max_af"], errors="ignore").merge(
        recomputed, on="variant_id", how="left"
    )
    first_seen = (
        filtered_long[filtered_long["present"]]
        .groupby("variant_id")["sample_number"]
        .min()
        .rename("first_seen")
    )
    last_seen = (
        filtered_long[filtered_long["present"]]
        .groupby("variant_id")["sample_number"]
        .max()
        .rename("last_seen")
    )
    filtered_summary = filtered_summary.drop(
        columns=["first_seen", "last_seen", "duration"], errors="ignore"
    )
    filtered_summary = filtered_summary.merge(first_seen, on="variant_id", how="left")
    filtered_summary = filtered_summary.merge(last_seen, on="variant_id", how="left")
    filtered_summary["duration"] = filtered_summary["last_seen"].fillna(0).astype(
        float
    ) - filtered_summary["first_seen"].fillna(0).astype(float)

    if min_af is not None:
        filtered_summary = filtered_summary[filtered_summary["max_af"] >= min_af]
    if max_af is not None:
        filtered_summary = filtered_summary[filtered_summary["max_af"] <= max_af]
    if filtered_summary.empty:
        raise ProcessingError("No variants remain after allele-frequency filtering")

    filtered_long = filtered_long[
        filtered_long["variant_id"].isin(filtered_summary["variant_id"])
    ]
    return filtered_summary.reset_index(drop=True), filtered_long.reset_index(drop=True)


def auto_select_variants(
    summary: pd.DataFrame,
    *,
    top_n: int,
    prefer_crossing: bool = False,
    thresholds: Sequence[float] | None = None,
) -> list[str]:
    if summary.empty:
        return []

    ranked = summary.copy()
    threshold_values = [float(value) for value in (thresholds or [])]
    if prefer_crossing and threshold_values:
        ranked["crosses_threshold"] = ranked["max_af"].apply(
            lambda value: any(
                float(value) >= threshold for threshold in threshold_values
            )
        )
    else:
        ranked["crosses_threshold"] = False

    ranked = ranked.sort_values(
        by=[
            "crosses_threshold",
            "has_literature",
            "is_persistent_new",
            "is_nonsynonymous",
            "max_af",
            "first_seen",
            "variant_id",
        ],
        ascending=[False, False, False, False, False, True, True],
        na_position="last",
    )
    return ranked["variant_id"].head(max(1, int(top_n))).tolist()


def _resolve_variant_selection(
    summary: pd.DataFrame,
    *,
    explicit_variants: Sequence[str],
    top_n: int,
    prefer_crossing: bool = False,
    thresholds: Sequence[float] | None = None,
) -> list[str]:
    if explicit_variants:
        selected: list[str] = []
        seen: set[str] = set()
        for requested in explicit_variants:
            match = summary[
                summary.apply(
                    lambda row: requested
                    in {
                        row["variant_id"],
                        row["variant_label"],
                        _display_label_prefix(row["variant_label"]),
                        row["variant_name"],
                    },
                    axis=1,
                )
            ]
            if match.empty:
                continue
            variant_id = str(match.iloc[0]["variant_id"])
            if variant_id not in seen:
                selected.append(variant_id)
                seen.add(variant_id)
        if not selected:
            raise ProcessingError(
                "None of the requested variants were found in the filtered results"
            )
        return selected
    return auto_select_variants(
        summary,
        top_n=top_n,
        prefer_crossing=prefer_crossing,
        thresholds=thresholds,
    )


def resolve_plot_output_path(
    results_csv: str | Path,
    *,
    out: str | None,
    outdir: str | None,
    fmt: str,
    filename: str,
) -> Path:
    results_path = Path(results_csv).expanduser().resolve()
    if out:
        return Path(out).expanduser().resolve()
    destination_dir = (
        Path(outdir).expanduser().resolve() if outdir else results_path.parent
    )
    destination_dir.mkdir(parents=True, exist_ok=True)
    return destination_dir / f"{filename}.{fmt}"


def _apply_plot_title(ax, title: str | None, default_title: str) -> None:
    ax.set_title(title or default_title, fontweight="bold")


def _sample_axis_ticks_from_long_df(
    long_df: pd.DataFrame,
) -> tuple[list[int], list[str]]:
    sample_axis = (
        long_df[["sample_number", "sample_name"]]
        .drop_duplicates()
        .sort_values(["sample_number", "sample_name"])
    )
    return (
        sample_axis["sample_number"].astype(int).tolist(),
        sample_axis["sample_name"].astype(str).tolist(),
    )


def _apply_sample_axis(
    ax,
    *,
    sample_numbers: Sequence[int],
    sample_names: Sequence[str],
    sample_axis_mode: str,
) -> None:
    ax.set_xticks(list(sample_numbers))
    if sample_axis_mode == "name":
        ax.set_xlabel("Sample Name")
        ax.set_xticklabels(list(sample_names), rotation=45, ha="right")
    else:
        ax.set_xlabel("Sample Number")
        ax.set_xticklabels([str(value) for value in sample_numbers])


def plot_variant_trajectory(
    summary: pd.DataFrame,
    long_df: pd.DataFrame,
    *,
    selected_variants: Sequence[str],
    output_path: str | Path,
    thresholds: Sequence[float] | None = None,
    title: str | None = None,
    width: float = 8.5,
    height: float = 5.5,
    dpi: int = 300,
    label_lines: bool = False,
    label_threshold_crossers: bool = False,
    crossing_rule: str = "at_or_above",
    label_mode: str = "aa",
    sample_axis_mode: str = "number",
) -> None:
    plot_df = long_df[long_df["variant_id"].isin(selected_variants)].copy()
    if plot_df.empty:
        raise ProcessingError("No data available for the trajectory plot")

    fig, ax = plt.subplots(figsize=(width, height))
    threshold_values = [float(value) for value in (thresholds or [])]
    for variant_id in selected_variants:
        subset = plot_df[plot_df["variant_id"] == variant_id].sort_values(
            "sample_number"
        )
        if label_mode == "nt":
            label = (
                str(subset["variant_name"].iloc[0]).strip()
                or str(subset["variant_label"].iloc[0]).split(" (", 1)[0]
            )
        else:
            label = str(subset["variant_label"].iloc[0]).split(" (", 1)[0]
        ax.plot(
            subset["sample_number"],
            subset["allele_frequency"],
            marker="o",
            linewidth=1.8,
            label=label,
        )
        crosses_threshold = variant_crosses_thresholds(
            subset["allele_frequency"].tolist(),
            threshold_values,
            crossing_rule=crossing_rule,
        )
        if (
            label_lines or (label_threshold_crossers and crosses_threshold)
        ) and not subset.empty:
            last_row = subset.iloc[-1]
            ax.text(
                float(last_row["sample_number"]) + 0.1,
                float(last_row["allele_frequency"]),
                label,
                fontsize=8,
                va="center",
            )

    ax.set_ylabel("Allele frequency")
    ax.set_ylim(0, 1.02)
    axis_numbers, axis_names = _sample_axis_ticks_from_long_df(plot_df)
    _apply_sample_axis(
        ax,
        sample_numbers=axis_numbers,
        sample_names=axis_names,
        sample_axis_mode=sample_axis_mode,
    )
    for threshold in threshold_values:
        ax.axhline(
            float(threshold), linestyle="--", linewidth=1.0, color="black", alpha=0.6
        )
    _apply_plot_title(ax, title, "Variant trajectories")
    if not label_lines and not label_threshold_crossers:
        ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_variant_turnover(
    summary: pd.DataFrame,
    long_df: pd.DataFrame,
    *,
    output_path: str | Path,
    title: str | None = None,
    width: float = 8.5,
    height: float = 5.0,
    dpi: int = 300,
    count_mode: str = "count",
    sample_axis_mode: str = "number",
) -> None:
    plot_df = long_df.sort_values(["variant_id", "sample_number"]).copy()
    records: list[dict[str, float | int]] = []

    for _, variant_df in plot_df.groupby("variant_id"):
        variant_df = variant_df.sort_values("sample_number").reset_index(drop=True)
        for idx in range(1, len(variant_df)):
            previous = variant_df.iloc[idx - 1]
            current = variant_df.iloc[idx]
            if not bool(previous["present"]) and bool(current["present"]):
                value = (
                    1.0 if count_mode == "count" else float(current["allele_frequency"])
                )
                records.append(
                    {
                        "sample_number": int(current["sample_number"]),
                        "new": value,
                        "lost": 0.0,
                    }
                )
            if bool(previous["present"]) and not bool(current["present"]):
                value = (
                    1.0
                    if count_mode == "count"
                    else float(previous["allele_frequency"])
                )
                records.append(
                    {
                        "sample_number": int(current["sample_number"]),
                        "new": 0.0,
                        "lost": value,
                    }
                )

    sample_numbers = sorted(plot_df["sample_number"].unique())
    turnover = pd.DataFrame({"sample_number": sample_numbers, "new": 0.0, "lost": 0.0})
    if records:
        events = pd.DataFrame(records).groupby("sample_number", as_index=False).sum()
        turnover = turnover.merge(
            events, on="sample_number", how="left", suffixes=("", "_event")
        )
        turnover["new"] = turnover["new_event"].fillna(turnover["new"])
        turnover["lost"] = turnover["lost_event"].fillna(turnover["lost"])
        turnover = turnover.drop(columns=["new_event", "lost_event"])

    fig, ax = plt.subplots(figsize=(width, height))
    ax.bar(turnover["sample_number"], turnover["new"], label="New")
    ax.bar(turnover["sample_number"], -turnover["lost"], label="Lost")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel(
        "Variant count" if count_mode == "count" else "Summed allele frequency"
    )
    axis_numbers, axis_names = _sample_axis_ticks_from_long_df(plot_df)
    _apply_sample_axis(
        ax,
        sample_numbers=axis_numbers,
        sample_names=axis_names,
        sample_axis_mode=sample_axis_mode,
    )
    _apply_plot_title(ax, title, "Variant turnover")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_variant_lifespan(
    summary: pd.DataFrame,
    *,
    selected_variants: Sequence[str],
    output_path: str | Path,
    title: str | None = None,
    width: float = 8.5,
    height: float = 6.0,
    dpi: int = 300,
    sort_by: str = "duration",
    annotate_class: bool = False,
    sample_axis_mode: str = "number",
    sample_names: Sequence[str] | None = None,
    sample_numbers: Sequence[int] | None = None,
) -> None:
    plot_df = summary[summary["variant_id"].isin(selected_variants)].copy()
    if plot_df.empty:
        raise ProcessingError("No data available for the lifespan plot")

    ascending = {
        "first_seen": True,
        "last_seen": True,
        "duration": False,
        "max_af": False,
    }
    plot_df = plot_df.sort_values(
        sort_by, ascending=ascending.get(sort_by, False), na_position="last"
    )
    plot_df = plot_df.reset_index(drop=True)

    labels = plot_df["variant_id"].astype(str).tolist()
    if annotate_class:
        labels = [
            f"{label} [{status}; {persistence}]"
            for label, status, persistence in zip(
                plot_df["variant_id"],
                plot_df["variant_status"],
                plot_df["persistence_status"],
            )
        ]

    fig, ax = plt.subplots(figsize=(width, height))
    y_positions = np.arange(len(plot_df))
    starts = plot_df["first_seen"].fillna(plot_df["last_seen"]).astype(float)
    ends = plot_df["last_seen"].fillna(plot_df["first_seen"]).astype(float)
    ax.hlines(y_positions, starts, ends, linewidth=2.2)
    ax.scatter(starts, y_positions, s=24, zorder=3)
    ax.scatter(ends, y_positions, s=24, zorder=3)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    if sample_numbers is None:
        inferred_numbers = sorted(
            {
                int(value)
                for value in pd.concat([starts, ends])
                .dropna()
                .astype(float)
                .astype(int)
                .tolist()
            }
        )
    else:
        inferred_numbers = [int(value) for value in sample_numbers]
    inferred_names = (
        [str(value) for value in sample_names]
        if sample_names is not None
        else [str(value) for value in inferred_numbers]
    )
    _apply_sample_axis(
        ax,
        sample_numbers=inferred_numbers,
        sample_names=inferred_names,
        sample_axis_mode=sample_axis_mode,
    )
    _apply_plot_title(ax, title, "Variant lifespan")
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def get_threshold_crossing_variants(
    summary: pd.DataFrame,
    long_df: pd.DataFrame,
    thresholds: Sequence[float],
    *,
    crossing_rule: str = "at_or_above",
) -> pd.DataFrame:
    threshold_values = [float(value) for value in thresholds]
    if not threshold_values:
        return summary.copy()
    crossed_ids: list[str] = []
    for variant_id, subset in long_df.groupby("variant_id"):
        if variant_crosses_thresholds(
            subset["allele_frequency"].tolist(),
            threshold_values,
            crossing_rule=crossing_rule,
        ):
            crossed_ids.append(str(variant_id))
    return summary[summary["variant_id"].isin(crossed_ids)].copy()


def parse_thresholds(value: str | None) -> list[float]:
    thresholds = _parse_csv_option_list(value)
    if not thresholds:
        return []
    try:
        return [float(item) for item in thresholds]
    except ValueError as exc:
        raise InputValidationError(
            "Thresholds must be comma-separated numbers"
        ) from exc


def variant_crosses_thresholds(
    values: Iterable[float],
    thresholds: Sequence[float],
    *,
    crossing_rule: str = "at_or_above",
) -> bool:
    if not thresholds:
        return False
    values_list = [float(value) for value in values]
    if crossing_rule == "strictly_above":
        return any(
            value > threshold for value in values_list for threshold in thresholds
        )
    return any(value >= threshold for value in values_list for threshold in thresholds)


def collect_explicit_variants(
    variants: str | None, variant_file: str | None
) -> list[str]:
    requested = _parse_csv_option_list(variants)
    file_variants = _read_variant_list_file(variant_file)
    if requested and file_variants:
        raise InputValidationError("Use either --variants or --variant-file, not both.")
    return requested or file_variants


def _parse_focus_ranges(value: str | None) -> list[tuple[float, float]]:
    if not value:
        return []
    ranges: list[tuple[float, float]] = []
    for part in _parse_csv_option_list(value):
        if "-" not in part:
            raise InputValidationError(
                "Focus ranges must use start-end syntax, e.g. 150-300,900-1800"
            )
        start_text, end_text = part.split("-", 1)
        try:
            start = float(start_text)
            end = float(end_text)
        except ValueError as exc:
            raise InputValidationError(
                "Focus ranges must use numeric start-end coordinates"
            ) from exc
        if end < start:
            start, end = end, start
        ranges.append((start, end))
    if len(ranges) > MAX_FOCUS_RANGES:
        raise InputValidationError(
            f"Use at most {MAX_FOCUS_RANGES} focus ranges so each can retain a distinct highlight color"
        )
    return ranges


def _parse_gff_attrs(attr_text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in attr_text.split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key] = unquote(value)
    return attrs


def build_reference_feature_metadata(gff_path: str | Path) -> dict[str, object]:
    path = Path(gff_path).expanduser().resolve()
    if not path.exists():
        raise InputValidationError(f"GFF3 not found: {path}")

    contig_lengths: dict[str, int] = {}
    contig_order: list[str] = []
    transcript_features: list[FeatureRecord] = []
    gene_features: list[FeatureRecord] = []
    max_end_by_contig: dict[str, int] = {}

    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    contig = parts[1]
                    seq_start = int(parts[2])
                    seq_end = int(parts[3])
                    if contig not in contig_order:
                        contig_order.append(contig)
                    contig_lengths[contig] = max(
                        contig_lengths.get(contig, 0), seq_end - seq_start + 1
                    )
                continue
            if line.startswith("#"):
                continue

            seqid, _source, feature_type, start, end, _score, strand, _phase, attrs = (
                line.split("\t")
            )
            start_i = int(start)
            end_i = int(end)
            max_end_by_contig[seqid] = max(max_end_by_contig.get(seqid, 0), end_i)
            if seqid not in contig_order:
                contig_order.append(seqid)

            parsed = _parse_gff_attrs(attrs)
            name = (
                parsed.get("Name")
                or parsed.get("gene")
                or parsed.get("product")
                or parsed.get("ID", "").split(":")[-1]
            )
            if not name:
                continue
            record = {
                "contig": seqid,
                "start": start_i,
                "end": end_i,
                "name": str(name),
                "strand": strand,
                "type": feature_type,
                "aa_length": max(1, int((end_i - start_i + 1) / 3)),
            }
            if feature_type in {"mRNA", "transcript"}:
                transcript_features.append(record)
            elif feature_type == "gene":
                gene_features.append(record)

    for contig, max_end in max_end_by_contig.items():
        contig_lengths.setdefault(contig, max_end)

    features = transcript_features or gene_features
    deduped_features: list[FeatureRecord] = []
    seen = set()
    for feature in features:
        key = (
            feature["contig"],
            feature["start"],
            feature["end"],
            feature["name"],
        )
        if key in seen:
            continue
        seen.add(key)
        deduped_features.append(feature)

    return {
        "source_gff3": str(path),
        "contig_order": contig_order,
        "contig_lengths": contig_lengths,
        "features": deduped_features,
    }


def write_reference_feature_metadata(gff_path: str | Path, outdir: str | Path) -> Path:
    outdir_path = Path(outdir).expanduser().resolve()
    outdir_path.mkdir(parents=True, exist_ok=True)
    metadata = build_reference_feature_metadata(gff_path)
    destination = outdir_path / "reference_features.json"
    destination.write_text(
        json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8"
    )
    return destination


def load_reference_feature_metadata(results_csv: str | Path) -> dict[str, object]:
    results_path = Path(results_csv).expanduser().resolve()
    sidecar = results_path.parent / "reference_features.json"
    if not sidecar.exists():
        raise InputValidationError(
            f"Reference feature metadata not found beside results CSV: {sidecar}"
        )
    return cast(dict[str, object], json.loads(sidecar.read_text(encoding="utf-8")))


def _collapse_variants_for_genome_plot(
    table: pd.DataFrame,
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for row in table.drop_duplicates().itertuples(index=False):
        gene_label, display_label, base_label = _resolve_variant_labels(row)
        variant_key = _build_variant_key(row, base_label)
        af_values = [
            _coerce_frequency(token)
            for token in _parse_slash_tokens(getattr(row, "alt_freq", ""))
        ]
        aa_position = _extract_numeric_position(
            base_label.split(":", 1)[1] if ":" in base_label else base_label
        )
        records.append(
            {
                "variant_id": variant_key,
                "variant_label": display_label.replace("\n", " "),
                "chrom": str(getattr(row, "chrom", "")).strip(),
                "start": int(getattr(row, "start", 0)),
                "end": int(getattr(row, "end", getattr(row, "start", 0))),
                "plot_gene": str(gene_label).strip(),
                "raw_gene": str(getattr(row, "gene", "")).strip(),
                "type_of_change": str(getattr(row, "type_of_change", "")).strip(),
                "type_of_variant": str(getattr(row, "type_of_variant", "")).strip(),
                "variant_status": str(getattr(row, "variant_status", "")).strip(),
                "persistence_status": str(
                    getattr(row, "persistence_status", "")
                ).strip(),
                "af_values": af_values,
                "summary_af": max(af_values) if af_values else 0.0,
                "aa_position": aa_position,
            }
        )

    collapsed = pd.DataFrame(records)
    if collapsed.empty:
        raise ProcessingError("No variants available for genome plotting")

    grouped_rows: list[dict[str, object]] = []
    skipped_variants = 0
    for variant_id, group in collapsed.groupby("variant_id", dropna=False):
        chroms = group["chrom"].dropna().unique().tolist()
        starts = group["start"].dropna().unique().tolist()
        ends = group["end"].dropna().unique().tolist()
        if len(chroms) > 1 or len(starts) > 1 or len(ends) > 1:
            skipped_variants += 1
            print(
                f"Warning: skipped variant '{variant_id}' due to inconsistent positional metadata"
            )
            continue
        first = group.iloc[0].to_dict()
        merged_af_values: list[float] = []
        for value in group["af_values"].tolist():
            merged_af_values.extend([float(item) for item in cast(list[float], value)])
        first["af_values"] = merged_af_values
        first["summary_af"] = float(group["summary_af"].max())
        grouped_rows.append(first)

    if skipped_variants:
        print(
            f"Warning: skipped {skipped_variants} variant(s) with inconsistent positional metadata"
        )
    if not grouped_rows:
        raise ProcessingError("No variants available for genome plotting")

    return (
        pd.DataFrame(grouped_rows)
        .sort_values(["chrom", "start", "variant_id"])
        .reset_index(drop=True)
    )


def _subset_collapsed_variants(
    collapsed: pd.DataFrame,
    *,
    gene: str | None = None,
    effects: Sequence[str] | None = None,
    min_af: float | None = None,
    max_af: float | None = None,
    persistent_only: bool = False,
    new_only: bool = False,
    include_indels: bool = False,
) -> pd.DataFrame:
    subset = collapsed.copy()
    if not include_indels and "type_of_variant" in subset.columns:
        subset = subset[subset["type_of_variant"].astype(str).str.lower() == "snp"]
    if gene:
        subset = subset[
            subset["plot_gene"].astype(str).str.lower() == str(gene).strip().lower()
        ]
    if effects:
        lowered = {
            str(effect).strip().lower() for effect in effects if str(effect).strip()
        }
        subset = subset[subset["type_of_change"].astype(str).str.lower().isin(lowered)]
    if min_af is not None:
        subset = subset[subset["summary_af"] >= min_af]
    if max_af is not None:
        subset = subset[subset["summary_af"] <= max_af]
    if persistent_only:
        subset = subset[subset["persistence_status"].eq("new_persistent")]
    if new_only:
        subset = subset[subset["variant_status"].eq("new")]
    if subset.empty:
        raise ProcessingError("No variants remain for the genome plot after filtering")
    return subset


def _find_gene_feature(metadata: dict[str, object], gene_name: str) -> FeatureRecord:
    features = cast(list[FeatureRecord], metadata.get("features", []))
    matches = [
        feature
        for feature in features
        if str(feature.get("name", "")).strip().lower() == gene_name.strip().lower()
    ]
    if not matches:
        raise InputValidationError(f"Gene not found in reference features: {gene_name}")
    unique_regions = {
        (match["contig"], match["start"], match["end"], match["name"])
        for match in matches
    }
    if len(unique_regions) != 1:
        raise InputValidationError(
            f"Gene '{gene_name}' maps to multiple regions; --aa-scale requires a unique gene region"
        )
    return matches[0]


def _clip_focus_ranges(
    ranges: Sequence[tuple[float, float]], xmin: float, xmax: float
) -> tuple[list[tuple[float, float]], int]:
    clipped: list[tuple[float, float]] = []
    ignored = 0
    for start, end in ranges:
        if end < xmin or start > xmax:
            ignored += 1
            continue
        clipped.append((max(start, xmin), min(end, xmax)))
    return clipped, ignored


def plot_variant_genome(
    table: pd.DataFrame,
    reference_metadata: dict[str, object],
    *,
    output_path: str | Path,
    gene: str | None = None,
    aa_scale: bool = False,
    focus_ranges: Sequence[tuple[float, float]] | None = None,
    min_af: float | None = None,
    max_af: float | None = None,
    effects: Sequence[str] | None = None,
    persistent_only: bool = False,
    new_only: bool = False,
    include_indels: bool = False,
    title: str | None = None,
    width: float = 10.0,
    height: float = 6.5,
    dpi: int = 300,
) -> None:
    if aa_scale and not gene:
        raise InputValidationError("--aa-scale requires --gene")
    if include_indels:
        print(
            "Warning: including indels in the genome plot may produce ambiguous or hard-to-interpret positions"
        )

    collapsed = _collapse_variants_for_genome_plot(table)
    collapsed = _subset_collapsed_variants(
        collapsed,
        gene=gene,
        effects=effects,
        min_af=min_af,
        max_af=max_af,
        persistent_only=persistent_only,
        new_only=new_only,
        include_indels=include_indels,
    )

    contig_order = cast(list[str], reference_metadata.get("contig_order", []))
    contig_lengths = {
        str(key): int(cast(int | float | str, value))
        for key, value in cast(
            dict[str, object], reference_metadata.get("contig_lengths", {})
        ).items()
    }

    region_by_contig: dict[str, tuple[float, float]] = {}
    plot_df = collapsed.copy()

    if gene:
        gene_feature = _find_gene_feature(reference_metadata, gene)
        contig = str(gene_feature["contig"])
        start = float(cast(int | float | str, gene_feature["start"]))
        end = float(cast(int | float | str, gene_feature["end"]))
        plot_df = plot_df[plot_df["chrom"] == contig]
        plot_df = plot_df[(plot_df["start"] >= start) & (plot_df["start"] <= end)]
        if plot_df.empty:
            raise ProcessingError(f"No variants remain in gene region '{gene}'")
        if aa_scale:
            omitted = plot_df["aa_position"].isna().sum()
            if omitted:
                print(
                    f"Warning: omitted {omitted} variants without amino-acid coordinates"
                )
            plot_df = plot_df.dropna(subset=["aa_position"]).copy()
            if plot_df.empty:
                raise ProcessingError(
                    "No variants could be mapped to amino-acid coordinates"
                )
            plot_df["x_coord"] = plot_df["aa_position"].astype(float)
            region_by_contig[contig] = (
                1.0,
                float(cast(int | float | str, gene_feature.get("aa_length", 1))),
            )
        else:
            plot_df["x_coord"] = plot_df["start"].astype(float)
            region_by_contig[contig] = (start, end)
    else:
        plot_df["x_coord"] = plot_df["start"].astype(float)
        for contig in plot_df["chrom"].unique():
            region_by_contig[str(contig)] = (
                1.0,
                float(
                    contig_lengths.get(
                        str(contig),
                        cast(
                            int | float,
                            plot_df.loc[plot_df["chrom"] == contig, "end"].max(),
                        ),
                    )
                ),
            )

    if not contig_order:
        contig_order = sorted(plot_df["chrom"].astype(str).unique().tolist())
    selected_contigs = [
        contig
        for contig in contig_order
        if contig in plot_df["chrom"].astype(str).unique().tolist()
    ]
    if not selected_contigs:
        raise ProcessingError("No contigs remain for genome plotting")

    n_panels = len(selected_contigs)
    fig = plt.figure(
        figsize=(width, max(height, 2.1 * n_panels)),
        constrained_layout=True,
    )
    grid = GridSpec(
        n_panels,
        1,
        figure=fig,
        hspace=0.35,
    )

    focus = list(focus_ranges or [])
    for idx, contig in enumerate(selected_contigs):
        ax = fig.add_subplot(grid[idx, 0])
        subset = plot_df[plot_df["chrom"].astype(str) == contig].sort_values("x_coord")
        xmin, xmax = region_by_contig[contig]
        clipped_focus, ignored_focus = _clip_focus_ranges(focus, xmin, xmax)
        for focus_index, (start, end) in enumerate(clipped_focus):
            ax.axvspan(
                start,
                end,
                color=FOCUS_RANGE_COLORS[focus_index % len(FOCUS_RANGE_COLORS)],
                alpha=0.25,
                zorder=0,
            )
        if ignored_focus:
            print(
                f"Warning: ignored {ignored_focus} focus range(s) outside the plotted region"
            )

        for row in subset.itertuples(index=False):
            af_values = [float(value) for value in cast(list[float], row.af_values)]
            if af_values:
                ax.vlines(
                    float(row.x_coord),
                    min(af_values),
                    max(af_values),
                    color="#607d8b",
                    linewidth=0.9,
                    alpha=0.55,
                    zorder=1,
                )
            ax.scatter(
                [float(row.x_coord)] * len(af_values),
                af_values,
                s=18,
                alpha=0.45,
                color="#1f77b4",
                edgecolors="none",
                zorder=2,
            )
        ax.axhline(0.5, color="black", linestyle="--", linewidth=0.8, alpha=0.4)
        ax.set_ylim(0, 1.02)
        ax.set_ylabel("Allele frequency")
        ax.set_xlim(xmin, xmax)
        ax.set_title(str(contig), loc="left", fontsize=10, fontweight="bold")
        ax.set_xlabel("Amino-acid position" if aa_scale else "Genomic position")
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    if title:
        fig.suptitle(title, fontweight="bold")
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
