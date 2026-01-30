"""Output schemas for vartracker results."""

from __future__ import annotations

import datetime as dt
import json
import shutil
import textwrap
from pathlib import Path
from typing import Any, Iterable, Sequence

RESULTS_SCHEMA_VERSION = "1.0"

RESULTS_SCHEMA: list[dict[str, str]] = [
    {
        "name": "chrom",
        "type": "string",
        "description": "Reference contig name (e.g., NC_045512.2).",
        "units": "",
        "values": "",
    },
    {
        "name": "start",
        "type": "integer",
        "description": "1-based start position of the variant on the reference.",
        "units": "bp",
        "values": "",
    },
    {
        "name": "end",
        "type": "integer",
        "description": "1-based end position (inclusive) of the variant.",
        "units": "bp",
        "values": "",
    },
    {
        "name": "gene",
        "type": "string",
        "description": "Gene or genomic region assigned to the variant.",
        "units": "",
        "values": "",
    },
    {
        "name": "ref",
        "type": "string",
        "description": "Reference allele.",
        "units": "",
        "values": "",
    },
    {
        "name": "alt",
        "type": "string",
        "description": "Alternate allele.",
        "units": "",
        "values": "",
    },
    {
        "name": "variant",
        "type": "string",
        "description": "Concise variant label (ref + position + alt, e.g. A23403G).",
        "units": "",
        "values": "",
    },
    {
        "name": "amino_acid_consequence",
        "type": "string",
        "description": "Amino acid consequence for the gene-level annotation.",
        "units": "",
        "values": "e.g. S:D614G",
    },
    {
        "name": "nsp_aa_change",
        "type": "string",
        "description": "NSP-level amino acid change for ORF1ab annotations.",
        "units": "",
        "values": "",
    },
    {
        "name": "bcsq_nt_notation",
        "type": "string",
        "description": "bcftools csq nucleotide notation.",
        "units": "",
        "values": "",
    },
    {
        "name": "bcsq_aa_notation",
        "type": "string",
        "description": "bcftools csq amino acid notation.",
        "units": "",
        "values": "",
    },
    {
        "name": "type_of_variant",
        "type": "string",
        "description": "Variant type derived from the VCF entry.",
        "units": "",
        "values": "snp, indel",
    },
    {
        "name": "type_of_change",
        "type": "string",
        "description": "Functional change classification from bcftools csq.",
        "units": "",
        "values": "synonymous, missense, frameshift, ...",
    },
    {
        "name": "variant_status",
        "type": "string",
        "description": "Whether the variant is present in the first sample.",
        "units": "",
        "values": "original, new",
    },
    {
        "name": "persistence_status",
        "type": "string",
        "description": "Persistence class based on first and last sample presence.",
        "units": "",
        "values": "original_retained, original_lost, new_persistent, new_transient",
    },
    {
        "name": "presence_absence",
        "type": "string (slash-separated)",
        "description": "Per-sample presence (Y) or absence (N), ordered by input.",
        "units": "",
        "values": "Y/N",
    },
    {
        "name": "first_appearance",
        "type": "string",
        "description": "Sample name where the variant first appears.",
        "units": "",
        "values": "",
    },
    {
        "name": "last_appearance",
        "type": "string",
        "description": "Sample name where the variant last appears.",
        "units": "",
        "values": "",
    },
    {
        "name": "overall_variant_qc",
        "type": "string",
        "description": "Aggregated QC status across samples.",
        "units": "",
        "values": "PASS, FAIL",
    },
    {
        "name": "per_sample_variant_qc",
        "type": "string (slash-separated)",
        "description": "Per-sample QC flags (P/F) ordered by input.",
        "units": "",
        "values": "P, F",
    },
    {
        "name": "aa1_total_properties",
        "type": "string",
        "description": "Physicochemical properties for the reference amino acid.",
        "units": "",
        "values": "semicolon-separated properties",
    },
    {
        "name": "aa2_total_properties",
        "type": "string",
        "description": "Physicochemical properties for the alternate amino acid.",
        "units": "",
        "values": "semicolon-separated properties",
    },
    {
        "name": "aa1_unique_properties",
        "type": "string",
        "description": "Properties unique to the reference amino acid.",
        "units": "",
        "values": "semicolon-separated properties",
    },
    {
        "name": "aa2_unique_properties",
        "type": "string",
        "description": "Properties unique to the alternate amino acid.",
        "units": "",
        "values": "semicolon-separated properties",
    },
    {
        "name": "aa1_weight",
        "type": "number",
        "description": "Molecular weight of the reference amino acid.",
        "units": "Da",
        "values": "",
    },
    {
        "name": "aa2_weight",
        "type": "number",
        "description": "Molecular weight of the alternate amino acid.",
        "units": "Da",
        "values": "",
    },
    {
        "name": "weight_difference",
        "type": "number",
        "description": "aa2_weight - aa1_weight.",
        "units": "Da",
        "values": "",
    },
    {
        "name": "alt_freq",
        "type": "string (slash-separated)",
        "description": "Allele frequency per sample (fraction).",
        "units": "fraction",
        "values": "0-1",
    },
    {
        "name": "variant_depth",
        "type": "string (slash-separated)",
        "description": "Alternate-allele read depth per sample.",
        "units": "reads",
        "values": "",
    },
    {
        "name": "variant_site_depth",
        "type": "string (slash-separated)",
        "description": "Total read depth at the variant site per sample.",
        "units": "reads",
        "values": "",
    },
    {
        "name": "variant_window_depth",
        "type": "string (slash-separated)",
        "description": "Mean read depth in the variant window per sample.",
        "units": "reads",
        "values": "",
    },
    {
        "name": "samples",
        "type": "string (slash-separated)",
        "description": "Sample names corresponding to per-sample fields.",
        "units": "",
        "values": "",
    },
    {
        "name": "total_genome_coverage",
        "type": "string (slash-separated)",
        "description": "Total genome coverage (bases covered) per sample.",
        "units": "bases",
        "values": "",
    },
]


def get_results_schema() -> Sequence[dict[str, str]]:
    return RESULTS_SCHEMA


def render_output_schema_markdown(
    schema: Iterable[dict[str, str]] | None = None,
) -> str:
    schema = list(schema or RESULTS_SCHEMA)
    lines = [
        "# Output schema",
        "",
        f"Schema version: `{RESULTS_SCHEMA_VERSION}`",
        "",
        "Generated from `vartracker.schemas.RESULTS_SCHEMA`.",
        "",
        "This document describes the columns produced in `results.csv`.",
        "",
        "Columns that encode per-sample values are slash-separated and ordered by the input CSV.",
        "",
        "| Column | Type | Description | Units | Values |",
        "| --- | --- | --- | --- | --- |",
    ]

    for entry in schema:
        lines.append(
            "| {name} | {type} | {description} | {units} | {values} |".format(
                name=entry["name"],
                type=entry["type"],
                description=entry["description"],
                units=entry.get("units", ""),
                values=entry.get("values", ""),
            )
        )

    lines.append("")
    return "\n".join(lines)


def _wrap_cell(text: str, width: int) -> list[str]:
    if width <= 0:
        return [text]
    return textwrap.wrap(text, width=width) or [""]


def render_output_schema_pretty(
    schema: Iterable[dict[str, str]] | None = None, *, width: int | None = None
) -> str:
    schema = list(schema or RESULTS_SCHEMA)
    term_width = width or shutil.get_terminal_size((120, 20)).columns

    columns = [
        ("Column", "name", 18),
        ("Type", "type", 16),
        ("Units", "units", 10),
        ("Values", "values", 22),
        ("Description", "description", None),
    ]

    fixed = sum(col_width for *_rest, col_width in columns if col_width is not None)
    separators = 3 * (len(columns) - 1)
    description_width = max(30, term_width - fixed - separators)

    computed_widths: list[int] = []
    for header, _key, col_width in columns:
        if col_width is None:
            computed_widths.append(description_width)
        else:
            computed_widths.append(col_width)

    header_cells = [
        header.ljust(width)
        for (header, _key, _), width in zip(columns, computed_widths)
    ]
    header_line = " | ".join(header_cells)
    divider_line = "-+-".join("-" * width for width in computed_widths)

    lines = [
        "Output schema",
        f"Schema version: {RESULTS_SCHEMA_VERSION}",
        "Generated from vartracker.schemas.RESULTS_SCHEMA.",
        "Columns that encode per-sample values are slash-separated and ordered by the input CSV.",
        "",
        header_line,
        divider_line,
    ]

    for entry in schema:
        wrapped_cells: list[list[str]] = []
        for (_header, key, _), col_width in zip(columns, computed_widths):
            value = str(entry.get(key, "")) if entry.get(key, "") is not None else ""
            wrapped_cells.append(_wrap_cell(value, col_width))

        row_height = max(len(cell) for cell in wrapped_cells)
        for idx in range(row_height):
            row_parts = []
            for cell, col_width in zip(wrapped_cells, computed_widths):
                text = cell[idx] if idx < len(cell) else ""
                row_parts.append(text.ljust(col_width))
            lines.append(" | ".join(row_parts))

    return "\n".join(lines)


def render_output_schema_text(schema: Iterable[dict[str, str]] | None = None) -> str:
    return render_output_schema_pretty(schema)


def write_results_metadata(
    outdir: str | Path,
    *,
    results_filename: str = "results.csv",
    schema_version: str = RESULTS_SCHEMA_VERSION,
) -> Path:
    outdir_path = Path(outdir).expanduser().resolve()
    metadata_path = outdir_path / "results_metadata.json"
    payload: dict[str, Any] = {
        "schema_version": schema_version,
        "results_file": str(outdir_path / results_filename),
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    metadata_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
    )
    return metadata_path


__all__ = [
    "RESULTS_SCHEMA_VERSION",
    "RESULTS_SCHEMA",
    "get_results_schema",
    "render_output_schema_markdown",
    "render_output_schema_pretty",
    "render_output_schema_text",
    "write_results_metadata",
]
