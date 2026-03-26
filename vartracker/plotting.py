"""Standalone plotting helpers for vartracker results tables."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .analysis import _coerce_frequency, _resolve_variant_labels
from .core import InputValidationError, ProcessingError

DEFAULT_TRAJECTORY_TOP_N = 12
DEFAULT_LIFESPAN_TOP_N = 20
DEFAULT_TRAJECTORY_THRESHOLDS = (0.5, 0.9)


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
    names = [str(value).strip() for value in table["name"].unique() if str(value).strip()]
    return names[0] if len(names) == 1 else ""


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
        names = _parse_slash_tokens(getattr(row, "samples", ""))
        numbers = _parse_slash_tokens(getattr(row, "sample_number", ""))
        freqs = [_coerce_frequency(token) for token in _parse_slash_tokens(getattr(row, "alt_freq", ""))]
        presence = _parse_slash_tokens(getattr(row, "presence_absence", ""))
        qc_flags = _sample_qc_flags(row, len(names))

        length = min(len(names), len(numbers), len(freqs)) if freqs else min(len(names), len(numbers))
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
                    "variant_id": base_label,
                    "variant_label": display_label.replace("\n", " "),
                    "variant_name": str(getattr(row, "variant", "")).strip(),
                    "gene": str(getattr(row, "gene", "")).strip(),
                    "type_of_change": str(getattr(row, "type_of_change", "")).strip(),
                    "type_of_variant": str(getattr(row, "type_of_variant", "")).strip(),
                    "variant_status": str(getattr(row, "variant_status", "")).strip(),
                    "persistence_status": str(getattr(row, "persistence_status", "")).strip(),
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
    long_df = long_df.sort_values(["sample_number", "variant_id", "sample_name"]).reset_index(drop=True)

    summary = (
        long_df.groupby("variant_id", as_index=False)
        .agg(
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
    )

    present_df = long_df[long_df["present"]]
    first_seen = (
        present_df.groupby("variant_id")["sample_number"].min().rename("first_seen")
    )
    last_seen = present_df.groupby("variant_id")["sample_number"].max().rename("last_seen")
    summary = summary.merge(first_seen, on="variant_id", how="left")
    summary = summary.merge(last_seen, on="variant_id", how="left")
    summary["first_seen"] = summary["first_seen"].fillna(summary["max_af"].map(lambda _x: np.nan))
    summary["last_seen"] = summary["last_seen"].fillna(summary["first_seen"])
    summary["duration"] = (
        summary["last_seen"].fillna(0).astype(float)
        - summary["first_seen"].fillna(0).astype(float)
    )
    summary["is_persistent_new"] = summary["persistence_status"].eq("new_persistent")
    summary["is_nonsynonymous"] = ~summary["type_of_change"].str.lower().str.contains(
        "synonymous", na=False
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
        raise ProcessingError("No variant observations remain after sample-range filtering")

    filtered_summary = summary.copy()

    if gene:
        filtered_summary = filtered_summary[
            filtered_summary["gene"].astype(str).str.lower() == str(gene).strip().lower()
        ]

    if effects:
        lowered = {str(effect).strip().lower() for effect in effects if str(effect).strip()}
        filtered_summary = filtered_summary[
            filtered_summary["type_of_change"].astype(str).str.lower().isin(lowered)
        ]

    if not include_synonymous:
        filtered_summary = filtered_summary[
            ~filtered_summary["type_of_change"].astype(str).str.lower().str.contains("synonymous", na=False)
        ]

    if persistent_only:
        filtered_summary = filtered_summary[
            filtered_summary["persistence_status"].eq("new_persistent")
        ]

    if new_only:
        filtered_summary = filtered_summary[filtered_summary["variant_status"].eq("new")]

    if variants:
        order_map = {name: idx for idx, name in enumerate(variants)}
        filtered_summary = filtered_summary[
            filtered_summary.apply(
                lambda row: any(
                    candidate in order_map
                    for candidate in (
                        row["variant_id"],
                        row["variant_label"],
                        row["variant_name"],
                    )
                ),
                axis=1,
            )
        ]
        filtered_summary["selection_order"] = filtered_summary.apply(
            lambda row: min(
                order_map[candidate]
                for candidate in (row["variant_id"], row["variant_label"], row["variant_name"])
                if candidate in order_map
            ),
            axis=1,
        )
        filtered_summary = filtered_summary.sort_values("selection_order").drop(
            columns=["selection_order"]
        )

    if filtered_summary.empty:
        raise ProcessingError("No variants remain after applying plot filters")

    filtered_long = filtered_long[filtered_long["variant_id"].isin(filtered_summary["variant_id"])]
    recomputed = (
        filtered_long.groupby("variant_id", as_index=False)
        .agg(
            max_af=("allele_frequency", "max"),
        )
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
    filtered_summary["duration"] = (
        filtered_summary["last_seen"].fillna(0).astype(float)
        - filtered_summary["first_seen"].fillna(0).astype(float)
    )

    if min_af is not None:
        filtered_summary = filtered_summary[filtered_summary["max_af"] >= min_af]
    if max_af is not None:
        filtered_summary = filtered_summary[filtered_summary["max_af"] <= max_af]
    if filtered_summary.empty:
        raise ProcessingError("No variants remain after allele-frequency filtering")

    filtered_long = filtered_long[filtered_long["variant_id"].isin(filtered_summary["variant_id"])]
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
            lambda value: any(float(value) >= threshold for threshold in threshold_values)
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
                    in {row["variant_id"], row["variant_label"], row["variant_name"]},
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
            raise ProcessingError("None of the requested variants were found in the filtered results")
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
) -> None:
    plot_df = long_df[long_df["variant_id"].isin(selected_variants)].copy()
    if plot_df.empty:
        raise ProcessingError("No data available for the trajectory plot")

    fig, ax = plt.subplots(figsize=(width, height))
    threshold_values = [float(value) for value in (thresholds or [])]
    for variant_id in selected_variants:
        subset = plot_df[plot_df["variant_id"] == variant_id].sort_values("sample_number")
        label = str(subset["variant_id"].iloc[0])
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
        if (label_lines or (label_threshold_crossers and crosses_threshold)) and not subset.empty:
            last_row = subset.iloc[-1]
            ax.text(
                float(last_row["sample_number"]) + 0.1,
                float(last_row["allele_frequency"]),
                label,
                fontsize=8,
                va="center",
            )

    ax.set_xlabel("Sample")
    ax.set_ylabel("Allele frequency")
    ax.set_ylim(0, 1.02)
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    for threshold in threshold_values:
        ax.axhline(float(threshold), linestyle="--", linewidth=1.0, color="black", alpha=0.6)
        xmax = max(plot_df["sample_number"]) if not plot_df.empty else 0
        ax.text(
            float(xmax) + 0.1,
            float(threshold),
            f"{threshold:g}",
            fontsize=8,
            va="center",
            ha="left",
            color="black",
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
) -> None:
    plot_df = long_df.sort_values(["variant_id", "sample_number"]).copy()
    records: list[dict[str, float | int]] = []

    for _, variant_df in plot_df.groupby("variant_id"):
        variant_df = variant_df.sort_values("sample_number").reset_index(drop=True)
        for idx in range(1, len(variant_df)):
            previous = variant_df.iloc[idx - 1]
            current = variant_df.iloc[idx]
            if not bool(previous["present"]) and bool(current["present"]):
                value = 1.0 if count_mode == "count" else float(current["allele_frequency"])
                records.append(
                    {
                        "sample_number": int(current["sample_number"]),
                        "new": value,
                        "lost": 0.0,
                    }
                )
            if bool(previous["present"]) and not bool(current["present"]):
                value = 1.0 if count_mode == "count" else float(previous["allele_frequency"])
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
        turnover = turnover.merge(events, on="sample_number", how="left", suffixes=("", "_event"))
        turnover["new"] = turnover["new_event"].fillna(turnover["new"])
        turnover["lost"] = turnover["lost_event"].fillna(turnover["lost"])
        turnover = turnover.drop(columns=["new_event", "lost_event"])

    fig, ax = plt.subplots(figsize=(width, height))
    ax.bar(turnover["sample_number"], turnover["new"], label="New")
    ax.bar(turnover["sample_number"], -turnover["lost"], label="Lost")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Variant count" if count_mode == "count" else "Summed allele frequency")
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
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
) -> None:
    plot_df = summary[summary["variant_id"].isin(selected_variants)].copy()
    if plot_df.empty:
        raise ProcessingError("No data available for the lifespan plot")

    ascending = {"first_seen": True, "last_seen": True, "duration": False, "max_af": False}
    plot_df = plot_df.sort_values(sort_by, ascending=ascending.get(sort_by, False), na_position="last")
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
    ax.set_xlabel("Sample")
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
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
        raise InputValidationError("Thresholds must be comma-separated numbers") from exc


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
        return any(value > threshold for value in values_list for threshold in thresholds)
    return any(value >= threshold for value in values_list for threshold in thresholds)


def collect_explicit_variants(
    variants: str | None, variant_file: str | None
) -> list[str]:
    requested = _parse_csv_option_list(variants)
    file_variants = _read_variant_list_file(variant_file)
    if requested and file_variants:
        raise InputValidationError("Use either --variants or --variant-file, not both.")
    return requested or file_variants
