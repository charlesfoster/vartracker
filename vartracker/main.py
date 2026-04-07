"""
Main module for vartracker.

Contains the main function and argument parsing for the vartracker command-line interface.
"""

import argparse
import copy
import csv
import json
import io
import os
import shutil
import sys
import tempfile
from contextlib import ExitStack
from importlib import resources
from pathlib import Path
from typing import Sequence

import pandas as pd

from .analysis_launcher import run_workflow as run_e2e_workflow
from .generate import generate_csv
from .core import (
    get_logo,
    get_package_data_path,
    validate_input_file,
    check_dependencies,
    DependencyError,
    InputValidationError,
    ProcessingError,
    FILE_COLUMNS,
)
from .vcf_processing import format_vcf, merge_consequences, process_vcf
from .analysis import (
    process_joint_variants,
    generate_cumulative_lineplot,
    generate_gene_table,
    plot_gene_table,
    generate_variant_heatmap,
    search_literature,
)
from .plotting import (
    DEFAULT_LIFESPAN_TOP_N,
    DEFAULT_TRAJECTORY_TOP_N,
    apply_shared_plot_filters,
    auto_select_variants,
    load_reference_feature_metadata,
    collect_explicit_variants,
    get_threshold_crossing_variants,
    load_results_table,
    parse_thresholds,
    plot_variant_genome,
    plot_variant_lifespan,
    plot_variant_trajectory,
    plot_variant_turnover,
    prepare_plot_inputs,
    _project_name_from_results,
    resolve_focus_regions,
    resolve_plot_output_path,
    write_reference_feature_metadata,
)
from ._version import __version__
from .provenance import (
    RunManifest,
    collect_input_files,
    collect_tool_versions,
    compute_sha256,
)
from .schemas import (
    get_results_schema,
    render_output_schema_text,
    write_results_metadata,
)
from .data import parse_pokay as parse_pokay_module
from .annotation_processing import (
    gene_lengths_from_gff3,
    validate_reference_and_annotation,
)
from .reference_prepare import parse_accessions, prepare_reference_bundle

_RED = "\033[91m"
_YELLOW = "\033[93m"
_RESET = "\033[0m"
HelpFormatter = argparse.RawDescriptionHelpFormatter

LITERATURE_SCHEMA: list[dict[str, str]] = [
    {
        "name": "gene",
        "type": "string",
        "description": "Gene symbol.",
        "units": "",
        "values": "e.g. S",
    },
    {
        "name": "category",
        "type": "string",
        "description": 'Functional category/grouping criterion (e.g. "homoplasy").',
        "units": "",
        "values": "",
    },
    {
        "name": "mutation",
        "type": "string",
        "description": "Short-hand amino acid consequence notation.",
        "units": "",
        "values": "e.g. D614G",
    },
    {
        "name": "information",
        "type": "string",
        "description": "Information relating to the putative effect of this mutation.",
        "units": "",
        "values": "",
    },
    {
        "name": "reference",
        "type": "string",
        "description": "Semi-colon-delimited list of DOIs/references.",
        "units": "",
        "values": "e.g. 10.1038/s41586-020-2012-7; 10.1016/S0140-6736(20)30154-9",
    },
]


def _parse_csv_option_list(value: str | None) -> list[str]:
    if value is None:
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _add_heatmap_option_arguments(group: argparse._ArgumentGroup) -> None:
    group.add_argument(
        "--heatmap-aa-exclude",
        action="store",
        required=False,
        default="",
        help=(
            "Comma-separated `type_of_change` patterns to exclude from heatmaps "
            "(wildcards supported, e.g. synonymous,*frameshift*)"
        ),
    )
    group.add_argument(
        "--heatmap-aa-include",
        action="store",
        required=False,
        default="",
        help=(
            "Comma-separated `type_of_change` patterns to include in heatmaps "
            "(wildcards supported)"
        ),
    )
    group.add_argument(
        "--heatmap-include-joint",
        action="store_true",
        default=False,
        help="Include joint variants in heatmaps (default: exclude them)",
    )
    group.add_argument(
        "--heatmap-only-persistent",
        action="store_true",
        default=False,
        help="Only include variants with persistence_status == new_persistent",
    )
    group.add_argument(
        "--heatmap-only-new",
        action="store_true",
        default=False,
        help="Only include variants with variant_status == new",
    )
    group.add_argument(
        "--heatmap-gene-include",
        action="store",
        required=False,
        default="",
        help="Comma-separated gene patterns to include in heatmaps",
    )
    group.add_argument(
        "--heatmap-gene-exclude",
        action="store",
        required=False,
        default="",
        help="Comma-separated gene patterns to exclude from heatmaps",
    )
    group.add_argument(
        "--heatmap-variant-type",
        action="store",
        required=False,
        default="",
        help="Comma-separated variant type patterns to include (e.g. snp,indel)",
    )
    group.add_argument(
        "--heatmap-qc",
        action="store",
        required=False,
        default="",
        help=(
            "Comma-separated all-samples QC patterns to include "
            "(e.g. true,false,pass,fail)"
        ),
    )
    group.add_argument(
        "--min-prop-passing-qc",
        action="store",
        type=float,
        default=None,
        help="Minimum proportion of samples that must pass per-sample QC (0-1)",
    )
    group.add_argument(
        "--heatmap-min-persistence",
        action="store",
        type=int,
        default=None,
        help="Minimum number of samples in which a variant must be present",
    )
    group.add_argument(
        "--heatmap-min-max-af",
        action="store",
        type=float,
        default=None,
        help="Minimum maximum allele frequency across included samples",
    )
    group.add_argument(
        "--heatmap-min-sample-af",
        action="store",
        type=float,
        default=None,
        help="Minimum allele frequency that must be reached in at least one included sample",
    )
    group.add_argument(
        "--heatmap-sample-subset",
        action="store",
        required=False,
        default="",
        help="Comma-separated sample name patterns to plot",
    )
    group.add_argument(
        "--heatmap-hide-singletons",
        action="store_true",
        default=False,
        help="Hide variants present in only one included sample",
    )
    group.add_argument(
        "--heatmap-min-depth",
        action="store",
        type=int,
        default=None,
        help="Minimum site depth a variant must reach in at least one included sample",
    )


def _collect_heatmap_kwargs(args) -> dict[str, object]:
    return {
        "excluded_consequence_types": _parse_csv_option_list(
            getattr(args, "heatmap_aa_exclude", "")
        ),
        "included_consequence_types": _parse_csv_option_list(
            getattr(args, "heatmap_aa_include", "")
        ),
        "include_joint": getattr(args, "heatmap_include_joint", False),
        "only_persistent": getattr(args, "heatmap_only_persistent", False),
        "only_new": getattr(args, "heatmap_only_new", False),
        "gene_include": _parse_csv_option_list(
            getattr(args, "heatmap_gene_include", "")
        ),
        "gene_exclude": _parse_csv_option_list(
            getattr(args, "heatmap_gene_exclude", "")
        ),
        "variant_type_include": _parse_csv_option_list(
            getattr(args, "heatmap_variant_type", "")
        ),
        "qc_include": _parse_csv_option_list(getattr(args, "heatmap_qc", "")),
        "min_prop_passing_qc": getattr(args, "min_prop_passing_qc", None),
        "min_persistence": getattr(args, "heatmap_min_persistence", None),
        "min_max_af": getattr(args, "heatmap_min_max_af", None),
        "min_sample_af": getattr(args, "heatmap_min_sample_af", None),
        "sample_subset": _parse_csv_option_list(
            getattr(args, "heatmap_sample_subset", "")
        ),
        "hide_singletons": getattr(args, "heatmap_hide_singletons", False),
        "min_depth": getattr(args, "heatmap_min_depth", None),
    }


def _add_shared_plot_filter_arguments(group: argparse._ArgumentGroup) -> None:
    group.add_argument("--gene", default=None, help="Restrict to a single gene")
    group.add_argument(
        "--effect",
        default="",
        help="Comma-separated effect classes to include (e.g. missense,synonymous)",
    )
    group.add_argument(
        "--min-af",
        type=float,
        default=None,
        help="Minimum max allele frequency across included samples",
    )
    group.add_argument(
        "--max-af",
        type=float,
        default=None,
        help="Maximum max allele frequency across included samples",
    )
    group.add_argument(
        "--variants",
        default="",
        help="Comma-separated variant identifiers to plot, in order",
    )
    group.add_argument(
        "--variant-file",
        default=None,
        help="Path to a file with one variant identifier per line",
    )
    group.add_argument(
        "--sample-min",
        type=int,
        default=None,
        help="Minimum sample_number to include",
    )
    group.add_argument(
        "--sample-max",
        type=int,
        default=None,
        help="Maximum sample_number to include",
    )
    group.add_argument(
        "--include-synonymous",
        action="store_true",
        default=True,
        help="Preserve synonymous variants in standalone plot filtering",
    )
    group.add_argument(
        "--persistent-only",
        action="store_true",
        default=False,
        help="Only include variants with persistence_status == new_persistent",
    )
    group.add_argument(
        "--new-only",
        action="store_true",
        default=False,
        help="Only include variants with variant_status == new",
    )


def _add_plot_output_arguments(parser: argparse.ArgumentParser) -> None:
    output_group = parser.add_argument_group("Output")
    output_group.add_argument(
        "--out", default=None, help="Write the plot to this exact path"
    )
    output_group.add_argument(
        "--outdir",
        default=None,
        help="Output directory for plot files (default: beside results.csv)",
    )
    output_group.add_argument(
        "--format",
        choices=["pdf", "png", "svg"],
        default="pdf",
        help="Output format when using --outdir or default naming (default: pdf)",
    )
    output_group.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output DPI (default: 300)",
    )


def _collect_shared_plot_kwargs(args) -> dict[str, object]:
    return {
        "gene": getattr(args, "gene", None),
        "effects": _parse_csv_option_list(getattr(args, "effect", "")),
        "min_af": getattr(args, "min_af", None),
        "max_af": getattr(args, "max_af", None),
        "variants": collect_explicit_variants(
            getattr(args, "variants", ""), getattr(args, "variant_file", None)
        ),
        "sample_min": getattr(args, "sample_min", None),
        "sample_max": getattr(args, "sample_max", None),
        "include_synonymous": getattr(args, "include_synonymous", True),
        "persistent_only": getattr(args, "persistent_only", False),
        "new_only": getattr(args, "new_only", False),
    }


def _resolve_plot_variants(
    summary: pd.DataFrame,
    *,
    explicit_variants: Sequence[str],
    top_n: int,
    prefer_crossing: bool = False,
    thresholds: Sequence[float] | None = None,
) -> list[str]:
    def _label_prefix(value: object) -> str:
        return str(value).split(" (", 1)[0].strip()

    if explicit_variants:
        selected = []
        seen = set()
        for requested in explicit_variants:
            match = summary[
                summary.apply(
                    lambda row: requested
                    in {
                        row["variant_id"],
                        row["variant_label"],
                        _label_prefix(row["variant_label"]),
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

    selected = auto_select_variants(
        summary,
        top_n=top_n,
        prefer_crossing=prefer_crossing,
        thresholds=thresholds,
    )
    if not selected:
        raise ProcessingError("No variants were available for plotting")
    return selected


def _print_dependency_error(error: DependencyError) -> None:
    """Render a dependency error with optional remediation tips."""
    message = str(error)
    print(f"\n{_RED}{message}{_RESET}")
    tip = getattr(error, "tip", None)
    if tip:
        print(f"{_YELLOW}[Tip]:{_RESET} {tip}")
    print()


def _start_manifest(args) -> RunManifest | None:
    outdir = getattr(args, "outdir", None)
    if not outdir:
        return None

    manifest = getattr(args, "_manifest", None)
    if manifest is not None:
        return manifest

    manifest = RunManifest(
        outdir=outdir,
        command=getattr(args, "command", None),
        cli_args=getattr(args, "_argv", None),
        invocation=getattr(args, "_invocation", None),
        manifest_level=getattr(args, "manifest_level", "light"),
    )
    manifest.start()
    setattr(args, "_manifest", manifest)
    return manifest


def _update_manifest_inputs(
    manifest: RunManifest,
    *,
    input_csv: str | None,
    reference: str | None,
    gff3: str | None,
    input_rows=None,
    input_key: str = "input_csv",
    input_files_key: str = "input_files",
    manifest_level: str = "light",
) -> None:
    inputs: dict[str, object] = {}

    if input_csv:
        csv_path = Path(input_csv).expanduser().resolve()
        csv_entry: dict[str, object] = {"path": str(csv_path)}
        csv_hash, csv_error = compute_sha256(csv_path)
        csv_entry["sha256"] = csv_hash
        if csv_error:
            csv_entry["sha256_error"] = csv_error
        try:
            csv_entry["size_bytes"] = csv_path.stat().st_size
        except OSError as exc:
            csv_entry["size_bytes_error"] = str(exc)
        inputs[input_key] = csv_entry

        if input_rows is not None and manifest_level == "deep":
            input_files = collect_input_files(
                input_rows, csv_path=csv_path, include_hashes=True
            )
            inputs[input_files_key] = input_files
            inputs[f"{input_files_key}_count"] = len(input_files)
            inputs[f"{input_files_key}_total_bytes"] = sum(
                entry.get("size_bytes", 0)
                for entry in input_files
                if isinstance(entry.get("size_bytes"), int)
            )

    def _add_ref(label: str, value: str | None) -> None:
        if not value:
            return
        path = Path(value).expanduser().resolve()
        entry: dict[str, object] = {"path": str(path)}
        digest, error = compute_sha256(path)
        entry["sha256"] = digest
        if error:
            entry["sha256_error"] = error
        inputs[label] = entry

    _add_ref("reference_fasta", reference)
    _add_ref("annotation_gff3", gff3)

    if inputs:
        manifest.update(inputs=inputs)


def _ensure_tool_versions(manifest: RunManifest) -> None:
    if "external_tools" in manifest.payload:
        return
    manifest.update(external_tools=collect_tool_versions())


def _copy_test_dataset(target_dir: Path) -> None:
    """Copy the bundled test dataset into the target directory."""
    data_root = resources.files("vartracker").joinpath("test_data")
    if hasattr(data_root, "exists") and data_root.exists():
        with resources.as_file(data_root) as data_path:
            shutil.copytree(data_path, target_dir, dirs_exist_ok=True)
    else:
        fallback = Path(__file__).resolve().parent / "test_data"
        if not fallback.exists():
            raise InputValidationError("Bundled test data not found")
        shutil.copytree(fallback, target_dir, dirs_exist_ok=True)


def _make_csv_paths_absolute(csv_path: Path) -> None:
    """Rewrite CSV file paths to be absolute relative to the CSV location."""
    df = pd.read_csv(csv_path)
    base_dir = csv_path.parent.resolve()

    def _resolve(value: str, column: str) -> str:
        if pd.isna(value):
            return value
        text = str(value).strip()
        if not text:
            return ""
        path = Path(text)
        if not path.is_absolute():
            path = (base_dir / path).resolve()
        else:
            path = path.resolve()

        if column == "coverage" and not path.exists():
            candidates = []
            name = path.name
            if ".depth." in name:
                candidates.append(path.with_name(name.replace(".depth", "_depth", 1)))
            if "_depth." in name:
                candidates.append(path.with_name(name.replace("_depth", ".depth", 1)))
            if name.endswith(".coverage.txt"):
                candidates.append(
                    path.with_name(name.replace(".coverage", "_coverage"))
                )
            for candidate in candidates:
                if candidate.exists():
                    path = candidate
                    break
        return str(path)

    for column in FILE_COLUMNS:
        if column in df.columns:
            df[column] = df[column].apply(
                lambda value, col=column: _resolve(value, col)
            )

    df.to_csv(csv_path, index=False)


def _prepare_test_run(args, mode: str):
    """Prepare the temporary dataset used by --test for each mode."""
    temp_dir = tempfile.TemporaryDirectory(prefix="vartracker_demo_data_")
    base = Path(temp_dir.name)
    _copy_test_dataset(base)

    csv_map = {
        "vcf": "test_input_vcf.csv",
        "bam": "test_input_bam.csv",
        "e2e": "test_input_e2e.csv",
    }
    csv_name = csv_map.get(mode, "test_input_vcf.csv")
    csv_path = base / "inputs" / csv_name
    if not csv_path.exists():
        raise InputValidationError(
            f"Test input spreadsheet for '{mode}' mode not found at {csv_path}"
        )
    _make_csv_paths_absolute(csv_path)

    literature_path = base / "mock_literature" / "mock_literature.csv"
    if not literature_path.exists():
        raise InputValidationError("Bundled mock literature database not found")

    args.search_pokay = False
    args.literature_csv = str(literature_path)
    if args.name is None:
        args.name = f"vartracker_{mode}_demo"

    args._test_data_dir = str(base)

    return temp_dir, csv_path


def _configure_vcf_parser(
    parser: argparse.ArgumentParser,
    *,
    include_input_csv: bool,
    input_csv_required: bool = False,
) -> None:
    if include_input_csv:
        if input_csv_required:
            parser.add_argument("input_csv", help="Input CSV file")
        else:
            parser.add_argument("input_csv", nargs="?", help="Input CSV file")

    analysis_group = parser.add_argument_group("Vartracker Analysis Options")
    output_group = parser.add_argument_group("Vartracker Output Options")

    analysis_group.add_argument(
        "-r",
        "--reference",
        action="store",
        required=False,
        default=None,
        help="Reference genome (default: uses packaged SARS-CoV-2 reference)",
    )
    analysis_group.add_argument(
        "-g",
        "--gff3",
        action="store",
        required=False,
        default=None,
        help="GFF3 annotations to use (default: packaged SARS-CoV-2 annotations)",
    )
    analysis_group.add_argument(
        "--allele-frequency-tag",
        action="store",
        required=False,
        default="AF",
        help="INFO tag name for allele frequency (default: AF)",
    )
    analysis_group.add_argument(
        "-m",
        "--min-snv-freq",
        action="store",
        required=False,
        type=float,
        default=0.03,
        help="Minimum allele frequency of SNV variants to keep (default: 0.03)",
    )
    analysis_group.add_argument(
        "-M",
        "--min-indel-freq",
        action="store",
        required=False,
        type=float,
        default=0.1,
        help="Minimum allele frequency of indel variants to keep (default: 0.1)",
    )
    analysis_group.add_argument(
        "-d",
        "--min-depth",
        action="store",
        required=False,
        type=int,
        default=10,
        help="Minimum depth threshold for variant QC (default: 10)",
    )
    analysis_group.add_argument(
        "--sample-cap",
        action="store",
        type=int,
        help="Only analyse samples with sample_number less than or equal to this value",
        default=None,
    )
    analysis_group.add_argument(
        "--literature-csv",
        action="store",
        required=False,
        default=None,
        help="Path to a literature CSV file (see README for file structure)",
    )
    analysis_group.add_argument(
        "--search-pokay",
        action="store_true",
        help='Automatically download and search against the "pokay" SARS-CoV-2 literature database.',
        default=False,
    )
    analysis_group.add_argument(
        "--test",
        action="store_true",
        help="Run vartracker against the bundled demonstration dataset",
        default=False,
    )

    output_group.add_argument(
        "-n",
        "--name",
        action="store",
        required=False,
        default=None,
        help="Optional: add a column to results with the name specified here",
    )
    output_group.add_argument(
        "-o",
        "--outdir",
        action="store",
        required=False,
        default=".",
        help="Output directory for vartracker results (default: current directory)",
    )
    output_group.add_argument(
        "--manifest-level",
        choices=["light", "deep"],
        default="light",
        help="Manifest detail level for run metadata (default: light)",
    )
    output_group.add_argument(
        "-f",
        "--filename",
        action="store",
        required=False,
        default="results.csv",
        help="Output file name (default: results.csv)",
    )

    output_group.add_argument(
        "--debug",
        action="store_true",
        help="Print commands being run for debugging",
        default=False,
    )
    output_group.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep temporary files for debugging",
        default=False,
    )


def _move_action_group_after(
    parser: argparse.ArgumentParser, group_title: str, anchor_title: str
) -> None:
    groups = parser._action_groups
    try:
        group_index = next(
            i for i, grp in enumerate(groups) if grp.title == group_title
        )
        anchor_index = next(
            i for i, grp in enumerate(groups) if grp.title == anchor_title
        )
    except StopIteration:
        return

    group = groups.pop(group_index)
    if group_index < anchor_index:
        anchor_index -= 1
    groups.insert(anchor_index + 1, group)


def _add_vcf_subparser(subparsers):
    description = (
        "Analyse VCF-based longitudinal variant data, generating reports and plots."
    )
    vcf_parser = subparsers.add_parser(
        "vcf",
        help="Analyse VCF inputs",
        description=description,
        formatter_class=HelpFormatter,
        epilog="""
The input CSV file must contain the columns:

sample_name, sample_number, reads1, reads2, bam, vcf, coverage

In VCF mode the `reads1`, `reads2`, and `bam` entries may be blank, but `vcf`
and `coverage` must reference existing files. See the e2e mode if you need to
generate these columns with Snakemake first.
""",
    )

    _configure_vcf_parser(vcf_parser, include_input_csv=True)

    vcf_parser.set_defaults(handler=_run_vcf_command, _subparser=vcf_parser)


def _add_placeholder_subcommand(
    subparsers, name: str, aliases: tuple[str, ...] = ()
) -> None:
    parser = subparsers.add_parser(
        name,
        help=f"Placeholder for upcoming '{name}' functionality",
        aliases=list(aliases),
        description=f"The '{name}' workflow is not yet implemented.",
    )
    parser.set_defaults(handler=_placeholder_handler)


def _placeholder_handler(_args):
    print(
        "Subcommand not yet implemented. Please use 'vartracker vcf' for current functionality."
    )
    return 1


def _add_bam_subparser(subparsers):
    description = (
        "Process BAM inputs through the Snakemake workflow before running vartracker."
    )
    bam_parser = subparsers.add_parser(
        "bam",
        help="Run the BAM preprocessing workflow",
        description=description,
        formatter_class=HelpFormatter,
        epilog="""
Provide a CSV with columns:

sample_name,sample_number,reads1,reads2,bam,vcf,coverage

In BAM mode the `bam` column must point to existing files while `reads1`,
`reads2`, `vcf`, and `coverage` may be blank.
""",
    )

    snk_group = bam_parser.add_argument_group("Snakemake options")
    snk_group.add_argument(
        "--snakemake-outdir",
        default=None,
        help="Output directory for Snakemake artefacts (default: use --outdir)",
    )
    snk_group.add_argument(
        "--cores",
        type=int,
        default=8,
        help="Number of cores for Snakemake execution (default: 8)",
    )
    snk_group.add_argument(
        "--snakemake-dryrun",
        action="store_true",
        help="Perform a Snakemake dry-run and skip vartracker analysis",
        default=False,
    )
    snk_group.add_argument(
        "--rulegraph",
        default=None,
        help="Write the Snakemake rulegraph to a DOT file and exit",
    )
    snk_group.add_argument(
        "--verbose",
        action="store_true",
        help="Print Snakemake shell commands during execution",
        default=False,
    )
    snk_group.add_argument(
        "--redo",
        action="store_true",
        help="Force Snakemake to recompute all targets (forceall)",
        default=False,
    )

    _configure_vcf_parser(bam_parser, include_input_csv=True, input_csv_required=False)
    _move_action_group_after(
        bam_parser, "Snakemake options", "Vartracker Analysis Options"
    )

    bam_parser.set_defaults(handler=_run_bam_command, _subparser=bam_parser)


def _add_spreadsheet_subparser(subparsers):
    description = "Generate template CSV files for vartracker input."
    gen_parser = subparsers.add_parser(
        "spreadsheet",
        help="Generate input spreadsheets from an existing directory of files",
        description=description,
        formatter_class=HelpFormatter,
    )

    gen_parser.add_argument(
        "--mode",
        choices=["vcf", "bam", "e2e"],
        required=True,
        help="What type of pipeline the generated CSV should target",
    )
    gen_parser.add_argument(
        "--dir",
        required=True,
        help="Directory to scan for input files",
    )
    gen_parser.add_argument(
        "--out",
        default="generated_input.csv",
        help="Where to write the generated CSV (default: generated_input.csv)",
    )
    gen_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only report what would be generated without writing the CSV",
    )

    gen_parser.set_defaults(handler=_run_generate_command)


def _add_prepare_reference_subparser(subparsers):
    description = (
        "Build a bcftools-csq-compatible reference bundle from GenBank accessions."
    )
    parser = subparsers.add_parser(
        "reference",
        help="Prepare FASTA/GFF3 reference files from GenBank accessions",
        description=description,
        formatter_class=HelpFormatter,
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--accessions",
        help="Comma-separated GenBank nucleotide accession list",
    )
    group.add_argument(
        "--accession-file",
        help="Path to file with one GenBank accession per line",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for generated reference files",
    )
    parser.add_argument(
        "--prefix",
        default="reference",
        help="Output file prefix (default: reference)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing outputs in --outdir",
        default=False,
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep per-accession GenBank/FASTA/GFF3 intermediate files",
        default=False,
    )
    parser.add_argument(
        "--skip-csq-validation",
        action="store_true",
        help="Skip the bcftools csq smoke test validation",
        default=False,
    )
    parser.set_defaults(handler=_run_prepare_reference_command)


def _run_prepare_command(args):
    if getattr(args, "handler", None) is None or args.command in {"prep", "prepare"}:
        args._subparser.print_help()
        return 1
    return args.handler(args)


def _run_prepare_reference_command(args):
    try:
        accession_list = parse_accessions(
            accessions=args.accessions,
            accession_file=args.accession_file,
        )
        metadata = prepare_reference_bundle(
            accessions=accession_list,
            outdir=args.outdir,
            prefix=args.prefix,
            force=args.force,
            keep_intermediates=args.keep_intermediates,
            skip_csq_validation=args.skip_csq_validation,
            invocation=getattr(args, "_invocation", None),
            argv=getattr(args, "_argv", None),
        )
    except Exception as exc:
        print(f"\nERROR: {exc}\n")
        return 1

    outputs = metadata.get("outputs", {})
    print("Reference bundle prepared successfully.")
    print(f"FASTA: {outputs.get('fasta')}")
    print(f"GFF3: {outputs.get('gff3')}")
    print(f"FASTA index: {outputs.get('fai')}")
    print(f"Metadata: {outputs.get('metadata')}")
    return 0


def _add_prep_subparser(subparsers):
    description = "Prepare inputs and templates for vartracker workflows."
    prep_parser = subparsers.add_parser(
        "prepare",
        help="Prepare inputs for vartracker",
        description=description,
        aliases=["prep"],
        formatter_class=HelpFormatter,
    )
    prep_subparsers = prep_parser.add_subparsers(dest="prep_command")
    _add_spreadsheet_subparser(prep_subparsers)
    _add_prepare_reference_subparser(prep_subparsers)
    prep_parser.set_defaults(handler=_run_prepare_command, _subparser=prep_parser)


def _add_schema_subparser(subparsers):
    description = "Describe schemas for vartracker outputs and literature inputs."
    schema_parser = subparsers.add_parser(
        "schema",
        help="Print schemas for results tables or literature CSV input",
        description=description,
        aliases=["describe-output"],
        formatter_class=HelpFormatter,
    )
    schema_parser.add_argument(
        "schema_target",
        choices=["results", "literature"],
        help="Schema target to describe",
    )

    schema_parser.add_argument(
        "--out",
        help="Optional path to write the schema (CSV or JSON)",
        default=None,
    )
    schema_parser.add_argument(
        "--format",
        choices=["csv", "json"],
        default="csv",
        help="Output format when using --out (default: csv)",
    )

    schema_parser.set_defaults(handler=_run_describe_output_command)


def _run_plot_command(args):
    if getattr(args, "handler", None) is None or args.command == "plot":
        args._subparser.print_help()
        return 1
    return args.handler(args)


def _prepare_standalone_plot_inputs(args):
    results_csv = Path(args.results_csv).expanduser().resolve()
    table = load_results_table(results_csv)
    summary, long_df, sample_names, sample_numbers = prepare_plot_inputs(table)
    shared_kwargs = _collect_shared_plot_kwargs(args)
    summary, long_df = apply_shared_plot_filters(summary, long_df, **shared_kwargs)
    project_name = _project_name_from_results(table, getattr(args, "name", None))
    return (
        results_csv,
        table,
        summary,
        long_df,
        sample_names,
        sample_numbers,
        project_name,
        shared_kwargs,
    )


def _run_plot_trajectory_command(args):
    try:
        (
            results_csv,
            _table,
            summary,
            long_df,
            _sample_names,
            _sample_numbers,
            project_name,
            shared_kwargs,
        ) = _prepare_standalone_plot_inputs(args)
        thresholds = parse_thresholds(args.thresholds)
        if args.crossing_only and not thresholds:
            raise InputValidationError("--crossing-only requires --thresholds")
        if args.crossing_only:
            summary = get_threshold_crossing_variants(
                summary,
                long_df,
                thresholds,
                crossing_rule=args.crossing_rule,
            )
            long_df = long_df[long_df["variant_id"].isin(summary["variant_id"])]
            if summary.empty:
                raise ProcessingError("No variants crossed the requested thresholds")
        selected = _resolve_plot_variants(
            summary,
            explicit_variants=shared_kwargs["variants"],
            top_n=args.top_n,
            prefer_crossing=bool(thresholds),
            thresholds=thresholds,
        )
        output_path = resolve_plot_output_path(
            results_csv,
            out=args.out,
            outdir=args.outdir,
            fmt=args.format,
            filename="variant_trajectory_plot",
        )
        plot_variant_trajectory(
            summary,
            long_df,
            selected_variants=selected,
            output_path=output_path,
            thresholds=thresholds,
            title=args.title
            or (f"{project_name}: variant trajectories" if project_name else None),
            width=args.width,
            height=args.height,
            dpi=args.dpi,
            label_lines=args.label_lines,
            label_threshold_crossers=args.label_threshold_crossers,
            crossing_rule=args.crossing_rule,
            label_mode=args.label_mode,
            sample_axis_mode=args.sample_axis,
        )
        print(f"\nFinished: wrote {output_path}\n")
        return 0
    except (InputValidationError, ProcessingError) as exc:
        print(f"\nERROR: {exc}\n")
        return 1


def _run_plot_turnover_command(args):
    try:
        (
            results_csv,
            _table,
            summary,
            long_df,
            _sample_names,
            _sample_numbers,
            project_name,
            shared_kwargs,
        ) = _prepare_standalone_plot_inputs(args)
        if shared_kwargs["variants"]:
            summary = summary[summary["variant_id"].isin(shared_kwargs["variants"])]
            long_df = long_df[long_df["variant_id"].isin(summary["variant_id"])]
            if summary.empty:
                raise ProcessingError(
                    "None of the requested variants were found in the filtered results"
                )
        output_path = resolve_plot_output_path(
            results_csv,
            out=args.out,
            outdir=args.outdir,
            fmt=args.format,
            filename="variant_turnover_plot",
        )
        plot_variant_turnover(
            summary,
            long_df,
            output_path=output_path,
            title=args.title
            or (f"{project_name}: variant turnover" if project_name else None),
            width=args.width,
            height=args.height,
            dpi=args.dpi,
            count_mode=args.count_mode,
            sample_axis_mode=args.sample_axis,
        )
        print(f"\nFinished: wrote {output_path}\n")
        return 0
    except (InputValidationError, ProcessingError) as exc:
        print(f"\nERROR: {exc}\n")
        return 1


def _run_plot_lifespan_command(args):
    try:
        (
            results_csv,
            _table,
            summary,
            _long_df,
            _sample_names,
            _sample_numbers,
            project_name,
            shared_kwargs,
        ) = _prepare_standalone_plot_inputs(args)
        selected = _resolve_plot_variants(
            summary,
            explicit_variants=shared_kwargs["variants"],
            top_n=args.top_n,
        )
        output_path = resolve_plot_output_path(
            results_csv,
            out=args.out,
            outdir=args.outdir,
            fmt=args.format,
            filename="variant_lifespan_plot",
        )
        plot_variant_lifespan(
            summary,
            selected_variants=selected,
            output_path=output_path,
            title=args.title
            or (f"{project_name}: variant lifespan" if project_name else None),
            width=args.width,
            height=args.height,
            dpi=args.dpi,
            sort_by=args.sort_by,
            annotate_class=args.annotate_class,
            sample_axis_mode=args.sample_axis,
            sample_names=_sample_names,
            sample_numbers=_sample_numbers,
        )
        print(f"\nFinished: wrote {output_path}\n")
        return 0
    except (InputValidationError, ProcessingError) as exc:
        print(f"\nERROR: {exc}\n")
        return 1


def _run_plot_genome_command(args):
    try:
        results_csv = Path(args.results_csv).expanduser().resolve()
        table = load_results_table(results_csv)
        project_name = _project_name_from_results(table, getattr(args, "name", None))
        metadata = load_reference_feature_metadata(results_csv)
        focus_ranges, focus_labels = resolve_focus_regions(
            args.focus_coords, args.focus_region_file
        )
        output_path = resolve_plot_output_path(
            results_csv,
            out=args.out,
            outdir=args.outdir,
            fmt=args.format,
            filename="variant_genome_plot",
        )
        plot_variant_genome(
            table,
            metadata,
            output_path=output_path,
            gene=args.gene,
            aa_scale=args.aa_scale,
            cds_scale=args.cds_scale,
            focus_ranges=focus_ranges,
            focus_labels=focus_labels,
            min_af=args.min_af,
            max_af=args.max_af,
            effects=_parse_csv_option_list(args.effect),
            persistent_only=args.persistent_only,
            new_only=args.new_only,
            include_indels=args.include_indels,
            show_intersections=args.show_intersections,
            title=args.title
            or (
                f"{project_name}: genome variant distribution" if project_name else None
            ),
            width=args.width,
            height=args.height,
            dpi=args.dpi,
        )
        print(f"\nFinished: wrote {output_path}\n")
        return 0
    except (InputValidationError, ProcessingError) as exc:
        print(f"\nERROR: {exc}\n")
        return 1


def _add_plot_heatmap_subparser(subparsers):
    parser = subparsers.add_parser(
        "heatmap",
        aliases=["hm"],
        help="Regenerate the variant heatmap from an existing results CSV",
        description="Read a vartracker results CSV and regenerate the heatmap outputs.",
        formatter_class=HelpFormatter,
    )
    parser.add_argument("results_csv", help="Path to a vartracker results CSV")
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory for regenerated heatmap files (default: results CSV directory)",
    )
    parser.add_argument(
        "--name",
        default=None,
        help="Optional plot title prefix (default: use the `name` column if present)",
    )
    parser.add_argument(
        "--literature-csv",
        default=None,
        help="Optional literature hits CSV to link from the interactive heatmap",
    )
    parser.add_argument(
        "-m",
        "--min-snv-freq",
        action="store",
        required=False,
        type=float,
        default=0.03,
        help="Minimum allele frequency of SNV variants to keep (default: 0.03)",
    )
    parser.add_argument(
        "-M",
        "--min-indel-freq",
        action="store",
        required=False,
        type=float,
        default=0.1,
        help="Minimum allele frequency of indel variants to keep (default: 0.1)",
    )
    heatmap_group = parser.add_argument_group("Heatmap options")
    _add_heatmap_option_arguments(heatmap_group)
    parser.set_defaults(handler=_run_plot_heatmap_command)


def _add_plot_trajectory_subparser(subparsers):
    parser = subparsers.add_parser(
        "trajectory",
        help="Plot selected variant trajectories from an existing results CSV",
        description="Generate allele-frequency trajectory plots from a vartracker results CSV.",
        formatter_class=HelpFormatter,
    )
    parser.add_argument("results_csv", help="Path to a vartracker results CSV")
    filter_group = parser.add_argument_group("Filtering")
    _add_shared_plot_filter_arguments(filter_group)
    parser.add_argument(
        "--top-n",
        type=int,
        default=DEFAULT_TRAJECTORY_TOP_N,
        help=f"Auto-select up to this many variants when --variants is not used (default: {DEFAULT_TRAJECTORY_TOP_N})",
    )
    parser.add_argument(
        "--label-lines",
        action="store_true",
        default=False,
        help="Add end labels to plotted lines",
    )
    parser.add_argument(
        "--label-mode",
        choices=["aa", "nt"],
        default="aa",
        help="Use amino-acid-style labels or nucleotide-level variant labels (default: aa)",
    )
    parser.add_argument(
        "--sample-axis",
        choices=["number", "name"],
        default="number",
        help="Use sample numbers or sample names on the x-axis (default: number)",
    )
    parser.add_argument(
        "--thresholds",
        default="",
        help="Comma-separated AF thresholds to draw, e.g. 0.5,0.9",
    )
    parser.add_argument(
        "--crossing-only",
        action="store_true",
        default=False,
        help="Only keep variants crossing at least one requested threshold",
    )
    parser.add_argument(
        "--label-threshold-crossers",
        action="store_true",
        default=False,
        help="Label variants that cross at least one supplied threshold",
    )
    parser.add_argument(
        "--crossing-rule",
        choices=["at_or_above", "strictly_above"],
        default="at_or_above",
        help="How to evaluate exact threshold matches (default: at_or_above)",
    )
    parser.add_argument("--title", default=None, help="Optional plot title")
    parser.add_argument(
        "--width", type=float, default=8.5, help="Figure width in inches"
    )
    parser.add_argument(
        "--height", type=float, default=5.5, help="Figure height in inches"
    )
    _add_plot_output_arguments(parser)
    parser.set_defaults(handler=_run_plot_trajectory_command)


def _add_plot_turnover_subparser(subparsers):
    parser = subparsers.add_parser(
        "turnover",
        help="Plot longitudinal variant turnover from an existing results CSV",
        description="Generate new-versus-lost turnover summaries from a vartracker results CSV.",
        formatter_class=HelpFormatter,
    )
    parser.add_argument("results_csv", help="Path to a vartracker results CSV")
    filter_group = parser.add_argument_group("Filtering")
    _add_shared_plot_filter_arguments(filter_group)
    parser.add_argument(
        "--count-mode",
        choices=["count", "sum_af"],
        default="count",
        help="Aggregate turnover as counts or summed allele frequencies (default: count)",
    )
    parser.add_argument(
        "--sample-axis",
        choices=["number", "name"],
        default="number",
        help="Use sample numbers or sample names on the x-axis (default: number)",
    )
    parser.add_argument("--title", default=None, help="Optional plot title")
    parser.add_argument(
        "--width", type=float, default=8.5, help="Figure width in inches"
    )
    parser.add_argument(
        "--height", type=float, default=5.0, help="Figure height in inches"
    )
    _add_plot_output_arguments(parser)
    parser.set_defaults(handler=_run_plot_turnover_command)


def _add_plot_lifespan_subparser(subparsers):
    parser = subparsers.add_parser(
        "lifespan",
        help="Plot first-to-last detection spans from an existing results CSV",
        description="Generate a horizontal lifespan plot from a vartracker results CSV.",
        formatter_class=HelpFormatter,
    )
    parser.add_argument("results_csv", help="Path to a vartracker results CSV")
    filter_group = parser.add_argument_group("Filtering")
    _add_shared_plot_filter_arguments(filter_group)
    parser.add_argument(
        "--top-n",
        type=int,
        default=DEFAULT_LIFESPAN_TOP_N,
        help=f"Auto-select up to this many variants when --variants is not used (default: {DEFAULT_LIFESPAN_TOP_N})",
    )
    parser.add_argument(
        "--sort-by",
        choices=["first_seen", "last_seen", "duration", "max_af"],
        default="duration",
        help="Sort plotted variants by this summary metric (default: duration)",
    )
    parser.add_argument(
        "--annotate-class",
        action="store_true",
        default=False,
        help="Append variant/new and persistent/transient classes to labels",
    )
    parser.add_argument(
        "--sample-axis",
        choices=["number", "name"],
        default="number",
        help="Use sample numbers or sample names on the x-axis (default: number)",
    )
    parser.add_argument("--title", default=None, help="Optional plot title")
    parser.add_argument(
        "--width", type=float, default=8.5, help="Figure width in inches"
    )
    parser.add_argument(
        "--height", type=float, default=6.0, help="Figure height in inches"
    )
    _add_plot_output_arguments(parser)
    parser.set_defaults(handler=_run_plot_lifespan_command)


def _add_plot_genome_subparser(subparsers):
    parser = subparsers.add_parser(
        "genome",
        help="Plot collapsed variant allele frequencies along the genome",
        description=(
            "Generate a genome-position summary plot from a vartracker results CSV."
        ),
        formatter_class=HelpFormatter,
    )
    parser.add_argument("results_csv", help="Path to a vartracker results CSV")
    filter_group = parser.add_argument_group("Filtering")
    filter_group.add_argument("--gene", default=None, help="Restrict to a single gene")
    filter_group.add_argument(
        "--effect",
        default="",
        help="Comma-separated effect classes to include (e.g. missense,synonymous)",
    )
    filter_group.add_argument(
        "--min-af",
        type=float,
        default=None,
        help="Minimum collapsed allele frequency to include",
    )
    filter_group.add_argument(
        "--max-af",
        type=float,
        default=None,
        help="Maximum collapsed allele frequency to include",
    )
    filter_group.add_argument(
        "--persistent-only",
        action="store_true",
        default=False,
        help="Only include variants with persistence_status == new_persistent",
    )
    filter_group.add_argument(
        "--new-only",
        action="store_true",
        default=False,
        help="Only include variants with variant_status == new",
    )
    filter_group.add_argument(
        "--include-indels",
        action="store_true",
        default=False,
        help="Include indels as well as SNPs in the genome plot (default: SNPs only)",
    )
    parser.add_argument(
        "--aa-scale",
        action="store_true",
        default=False,
        help="With --gene, plot x-axis in amino-acid coordinates for that gene",
    )
    parser.add_argument(
        "--cds-scale",
        action="store_true",
        default=False,
        help="With --gene, plot x-axis in CDS-relative nucleotide coordinates for that gene",
    )
    parser.add_argument(
        "--focus-coords",
        default="",
        help=(
            "Coordinate ranges to highlight. Use commas for separate ranges or "
            "semicolons to group ranges with the same color, e.g. "
            "150-300,900-1800;50-120 or Name:150-300,900-1800;Other:50-120"
        ),
    )
    parser.add_argument(
        "--focus-region-file",
        default="",
        help=(
            "Optional JSON/CSV/TSV file defining named focus region groups for the "
            "genome plot"
        ),
    )
    parser.add_argument(
        "--show-intersections",
        action="store_true",
        default=False,
        help="Show a compact table of variants intersecting the highlighted focus regions",
    )
    parser.add_argument("--title", default=None, help="Optional plot title")
    parser.add_argument(
        "--width", type=float, default=10.0, help="Figure width in inches"
    )
    parser.add_argument(
        "--height", type=float, default=6.5, help="Figure height in inches"
    )
    _add_plot_output_arguments(parser)
    parser.set_defaults(handler=_run_plot_genome_command)


def _add_plot_subparser(subparsers):
    plot_parser = subparsers.add_parser(
        "plot",
        help="Regenerate plots from existing vartracker outputs",
        description="Regenerate selected plots from an existing vartracker results file.",
        formatter_class=HelpFormatter,
    )
    plot_subparsers = plot_parser.add_subparsers(dest="plot_command")
    _add_plot_heatmap_subparser(plot_subparsers)
    _add_plot_genome_subparser(plot_subparsers)
    _add_plot_trajectory_subparser(plot_subparsers)
    _add_plot_turnover_subparser(plot_subparsers)
    _add_plot_lifespan_subparser(plot_subparsers)
    plot_parser.set_defaults(handler=_run_plot_command, _subparser=plot_parser)


def create_parser():
    """Create and return the top-level argument parser with subcommands."""

    parser = argparse.ArgumentParser(
        description="vartracker: longitudinal variant tracking and summarisation for pathogen sequencing",
        formatter_class=HelpFormatter,
    )
    parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {__version__}"
    )

    subparsers = parser.add_subparsers(dest="command")
    _add_vcf_subparser(subparsers)
    _add_bam_subparser(subparsers)
    _add_e2e_subparser(subparsers)
    _add_plot_subparser(subparsers)
    _add_prep_subparser(subparsers)
    _add_schema_subparser(subparsers)

    return parser


def validate_dependencies(mode: str = "vcf"):
    deps = check_dependencies(mode)
    missing = [tool for tool, available in deps.items() if not available]
    if missing:
        missing_list = ", ".join(missing)
        install_hint = f"conda install -c bioconda {' '.join(missing)}"
        message = (
            f"Dependency error: to run the {mode} mode the following tools must be available: "
            f"{missing_list}"
        )
        tip = f"try installing missing dependencies with (e.g.) `{install_hint}`"
        raise DependencyError(message, mode=mode, missing=missing, tip=tip)


def setup_default_paths(args):
    """Set up default paths for reference files if not provided."""
    if args.reference is None:
        args.reference = get_package_data_path("NC_045512.fasta")

    annotation = getattr(args, "gff3", None)
    if annotation is None:
        annotation = get_package_data_path("NC_045512.gff3")
    args.gff3 = annotation

    return args


def _normalise_rulegraph_path(path: str | None) -> str | None:
    if not path:
        return None
    resolved = str(Path(path).expanduser())
    if not resolved.lower().endswith(".dot"):
        resolved = f"{resolved}.dot"
    return str(Path(resolved).resolve())


def _drop_exact_duplicate_result_rows(table: pd.DataFrame) -> pd.DataFrame:
    deduped = table.drop_duplicates().reset_index(drop=True)
    removed = len(table) - len(deduped)
    if removed:
        print(f"Removed {removed} exact duplicate result rows.")
    return deduped


def main(sysargs=None):
    """Entry point for the vartracker CLI."""
    if sysargs is None:
        sysargs = sys.argv[1:]

    parser = create_parser()

    try:
        args = parser.parse_args(sysargs)
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else 0
        return code

    invocation = "vartracker"
    if sysargs:
        invocation = f"vartracker {' '.join(sysargs)}"
    setattr(args, "_invocation", invocation)
    setattr(args, "_argv", list(sysargs))

    handler = getattr(args, "handler", None)
    if handler is None:
        parser.print_help()
        return 1

    return handler(args)


def _run_vcf_command(args):
    annotation_supplied = getattr(args, "gff3", None) is not None

    with ExitStack() as stack:
        manifest = None
        input_key = "input_csv"
        input_csv_path = args.input_csv

        if args.test:
            temp_dir, csv_path = _prepare_test_run(args, "vcf")
            stack.enter_context(temp_dir)
            input_csv_path = str(csv_path.resolve())
            if args.outdir == ".":
                args.outdir = os.path.join(os.getcwd(), "vartracker_test_results")
        elif input_csv_path is None:
            args._subparser.print_help()
            return 1

        args.input_csv = input_csv_path

        if not getattr(args, "_suppress_logo", False):
            print(get_logo())

        try:
            # Set up default paths
            args = setup_default_paths(args)

            # Create output directory and manifest as early as possible
            os.makedirs(args.outdir, exist_ok=True)
            manifest = _start_manifest(args)
            if manifest is not None:
                _ensure_tool_versions(manifest)
                existing_input = (
                    manifest.payload.get("inputs", {}).get("input_csv", {}).get("path")
                )
                resolved_input = str(Path(args.input_csv).expanduser().resolve())
                if existing_input and existing_input != resolved_input:
                    input_key = "analysis_csv"
                _update_manifest_inputs(
                    manifest,
                    input_csv=args.input_csv,
                    reference=args.reference,
                    gff3=args.gff3,
                    input_key=input_key,
                    manifest_level=args.manifest_level,
                )

            # Validate dependencies
            validate_dependencies()

            try:
                validate_reference_and_annotation(args.reference, args.gff3)
            except ValueError as exc:
                raise InputValidationError(str(exc)) from exc

            gene_lengths_override = None
            if annotation_supplied:
                try:
                    gene_lengths_override = gene_lengths_from_gff3(args.gff3)
                except (
                    Exception
                ) as exc:  # pragma: no cover - parsing errors surfaced to user
                    raise InputValidationError(
                        f"Failed to parse annotation file '{args.gff3}': {exc}"
                    ) from exc

            # Read input CSV
            try:
                input_file = pd.read_csv(args.input_csv)
            except Exception as e:
                raise InputValidationError(
                    f"Could not read input file: {args.input_csv}. {str(e)}"
                )

            if manifest is not None and args.manifest_level == "deep":
                input_files_key = (
                    "analysis_input_files"
                    if input_key == "analysis_csv"
                    else "input_files"
                )
                _update_manifest_inputs(
                    manifest,
                    input_csv=args.input_csv,
                    reference=args.reference,
                    gff3=args.gff3,
                    input_rows=input_file.to_dict(orient="records"),
                    input_key=input_key,
                    input_files_key=input_files_key,
                    manifest_level=args.manifest_level,
                )

            # Resolve potential relative paths for VCF and coverage files
            csv_dir = os.path.dirname(os.path.abspath(args.input_csv))

            def resolve_path(value: str) -> str:
                path_str = str(value)
                if os.path.isabs(path_str):
                    return os.path.abspath(path_str)

                cwd_candidate = os.path.abspath(path_str)
                if os.path.exists(cwd_candidate):
                    return cwd_candidate

                return os.path.abspath(os.path.join(csv_dir, path_str))

            for column in ("vcf", "coverage"):
                if column in input_file.columns:
                    input_file[column] = input_file[column].apply(resolve_path)

            # Validate input file
            validate_input_file(
                input_file,
                optional_empty={"reads1", "reads2", "bam"},
            )

            # Apply sample cap if specified
            if args.sample_cap is not None:
                input_file = input_file[input_file["sample_number"] <= args.sample_cap]

            literature = None
            if args.search_pokay and args.literature_csv is not None:
                raise InputValidationError(
                    "Use either --search-pokay or --literature-csv, not both."
                )

            if args.literature_csv is not None:
                try:
                    literature = pd.read_csv(args.literature_csv)
                except Exception as e:
                    raise InputValidationError(
                        f"Could not read literature CSV: {str(e)}"
                    )
                for column in ("gene", "mutation"):
                    if column not in literature.columns:
                        raise InputValidationError(
                            "Literature CSV must contain 'gene' and 'mutation' columns."
                        )

            # Set up variables
            vcfs = list(input_file["vcf"])
            sample_names = list(input_file["sample_name"])
            covs = list(input_file["coverage"])

            # Warning for single VCF input
            if len(vcfs) == 1:
                print(
                    "\033[93mWarning:\033[0m vartracker was designed for longitudinal comparisons but only one input VCF was provided. Some results may not make sense."
                )

            # Create temporary directory
            if args.keep_temp:
                # Create a persistent temporary directory for debugging
                tempdir = tempfile.mkdtemp(prefix="vartracker_")
                try:
                    # Process everything in the persistent tempdir
                    _process_files(
                        tempdir,
                        args,
                        vcfs,
                        sample_names,
                        covs,
                        input_file,
                        literature,
                        gene_lengths_override,
                    )

                    # Copy tempdir to outdir at the end
                    temp_copy_path = os.path.join(args.outdir, "temp_files")
                    if os.path.exists(temp_copy_path):
                        shutil.rmtree(temp_copy_path)
                    shutil.copytree(tempdir, temp_copy_path)
                    print(f"Temporary files copied to: {temp_copy_path}")

                finally:
                    # Clean up the original tempdir
                    shutil.rmtree(tempdir)
            else:
                # Use context manager for automatic cleanup
                with tempfile.TemporaryDirectory(prefix="vartracker_") as tempdir:
                    _process_files(
                        tempdir,
                        args,
                        vcfs,
                        sample_names,
                        covs,
                        input_file,
                        literature,
                        gene_lengths_override,
                    )

            results_metadata_path = write_results_metadata(
                args.outdir, results_filename=args.filename
            )

            outputs = {
                "results_csv": os.path.join(args.outdir, args.filename),
                "results_metadata": str(results_metadata_path),
                "reference_features": os.path.join(
                    args.outdir, "reference_features.json"
                ),
                "new_mutations_csv": os.path.join(args.outdir, "new_mutations.csv"),
                "persistent_new_mutations_csv": os.path.join(
                    args.outdir, "persistent_new_mutations.csv"
                ),
                "cumulative_mutations_plot": os.path.join(
                    args.outdir, "cumulative_mutations.pdf"
                ),
                "mutations_per_gene_plot": os.path.join(
                    args.outdir, "mutations_per_gene.pdf"
                ),
                "variant_turnover_plot": os.path.join(
                    args.outdir, "variant_turnover_plot.pdf"
                ),
                "variant_genome_plot": os.path.join(
                    args.outdir, "variant_genome_plot.pdf"
                ),
                "variant_allele_frequency_heatmap_html": os.path.join(
                    args.outdir, "variant_allele_frequency_heatmap.html"
                ),
                "variant_allele_frequency_heatmap_pdf": os.path.join(
                    args.outdir, "variant_allele_frequency_heatmap.pdf"
                ),
            }

            if args.search_pokay or args.literature_csv is not None:
                literature_name = args.name if args.name else "sample"
                literature_full = os.path.join(
                    args.outdir, f"{literature_name}.literature_database_hits.full.csv"
                )
                literature_concise = os.path.join(
                    args.outdir,
                    f"{literature_name}.literature_database_hits.concise.csv",
                )
                if os.path.exists(literature_full):
                    outputs["literature_hits_full_csv"] = literature_full
                if os.path.exists(literature_concise):
                    outputs["literature_hits_concise_csv"] = literature_concise
                parsed_literature = os.path.join(args.outdir, "literature_database.csv")
                if os.path.exists(parsed_literature):
                    outputs["literature_database_csv"] = parsed_literature

            if manifest is not None:
                manifest.finish(status="success", outputs=outputs)

            print(f"\nFinished: find results in {args.outdir}\n")
            return 0
        except DependencyError as e:
            _print_dependency_error(e)
            if manifest is not None:
                manifest.finish(status="failed", error=str(e))
            return 1
        except (InputValidationError, ProcessingError) as e:
            print(f"\nERROR: {str(e)}\n")
            if manifest is not None:
                manifest.finish(status="failed", error=str(e))
            return 1
        except Exception as e:
            print(f"\nUnexpected error: {str(e)}\n")
            if args.debug:
                import traceback

                traceback.print_exc()
            if manifest is not None:
                manifest.finish(status="failed", error=str(e))
            return 1


def _add_e2e_subparser(subparsers):
    description = "Run the Snakemake-based end-to-end pipeline and post-process results with vartracker."
    e2e_parser = subparsers.add_parser(
        "end-to-end",
        help="Run the end-to-end workflow (Snakemake + vartracker)",
        description=description,
        formatter_class=HelpFormatter,
        aliases=["e2e"],
    )

    e2e_parser.add_argument(
        "samples_csv",
        nargs="?",
        help="Path to read sample metadata CSV for the Snakemake workflow",
    )

    snk_group = e2e_parser.add_argument_group("Snakemake options")
    snk_group.add_argument(
        "--snakemake-outdir",
        default=None,
        help="Output directory for Snakemake artefacts (default: use --outdir)",
    )
    snk_group.add_argument(
        "--cores",
        type=int,
        default=8,
        help="Number of cores for Snakemake execution (default: 8)",
    )
    snk_group.add_argument(
        "--primer-bed",
        help="Optional primer BED file for amplicon clipping in Snakemake",
    )
    snk_group.add_argument(
        "--snakemake-dryrun",
        action="store_true",
        help="Perform a Snakemake dry-run and skip vartracker analysis",
        default=False,
    )
    snk_group.add_argument(
        "--rulegraph",
        default=None,
        help="Write the Snakemake rulegraph to a DOT file and exit",
    )
    snk_group.add_argument(
        "--verbose",
        action="store_true",
        help="Print Snakemake shell commands during execution",
        default=False,
    )
    snk_group.add_argument(
        "--redo",
        action="store_true",
        help="Force Snakemake to recompute all targets (forceall)",
        default=False,
    )

    _configure_vcf_parser(e2e_parser, include_input_csv=False)
    _move_action_group_after(
        e2e_parser, "Snakemake options", "Vartracker Analysis Options"
    )

    e2e_parser.set_defaults(handler=_run_e2e_command, _subparser=e2e_parser)


def _run_e2e_command(args):
    test_context = None
    manifest = None
    try:
        args = setup_default_paths(args)

        if getattr(args, "test", False):
            test_context, csv_path = _prepare_test_run(args, "e2e")
            args.samples_csv = str(csv_path.resolve())
            if args.outdir == ".":
                args.outdir = os.path.join(os.getcwd(), "vartracker_e2e_test_results")
            if args.snakemake_outdir is None:
                args.snakemake_outdir = args.outdir
        elif args.samples_csv is None:
            print(
                "ERROR: Missing required argument: samples_csv.\n"
                "Usage: vartracker end-to-end <samples_csv> [options]"
            )
            return 1

        os.makedirs(args.outdir, exist_ok=True)
        manifest = _start_manifest(args)
        if manifest is not None:
            _ensure_tool_versions(manifest)
            _update_manifest_inputs(
                manifest,
                input_csv=args.samples_csv,
                reference=args.reference,
                gff3=args.gff3,
                input_key="input_csv",
                manifest_level=args.manifest_level,
            )

        try:
            validate_dependencies("e2e")
        except DependencyError as error:
            _print_dependency_error(error)
            if manifest is not None:
                manifest.finish(status="failed", error=str(error))
            return 1

        snakemake_outdir = args.snakemake_outdir or args.outdir

        try:
            snakemake_input = pd.read_csv(args.samples_csv)
        except Exception as exc:
            raise InputValidationError(
                f"Could not read input file: {args.samples_csv}. {exc}"
            )

        if manifest is not None and args.manifest_level == "deep":
            _update_manifest_inputs(
                manifest,
                input_csv=args.samples_csv,
                reference=args.reference,
                gff3=args.gff3,
                input_rows=snakemake_input.to_dict(orient="records"),
                input_key="input_csv",
                manifest_level=args.manifest_level,
            )

        validate_input_file(
            snakemake_input,
            optional_empty={"bam", "vcf", "coverage", "reads2"},
        )

        print(get_logo())

        rulegraph_path = _normalise_rulegraph_path(args.rulegraph)
        args.rulegraph = rulegraph_path

        if rulegraph_path and args.snakemake_dryrun:
            raise InputValidationError(
                "--rulegraph cannot be combined with --snakemake-dryrun"
            )

        updated_csv = run_e2e_workflow(
            samples_csv=args.samples_csv,
            reference=args.reference,
            outdir=snakemake_outdir,
            cores=args.cores,
            primer_bed=args.primer_bed,
            dryrun=args.snakemake_dryrun,
            force_all=args.redo,
            quiet=not args.verbose,
            mode="reads",
            rulegraph_path=rulegraph_path,
        )

        if rulegraph_path:
            if manifest is not None:
                manifest.finish(
                    status="success",
                    outputs={"snakemake_rulegraph": rulegraph_path},
                )
            return 0

        if args.snakemake_dryrun:
            if manifest is not None:
                manifest.finish(
                    status="success",
                    outputs={"snakemake_outdir": snakemake_outdir},
                )
            return 0

        if not updated_csv or not os.path.exists(updated_csv):
            raise ProcessingError(
                "End-to-end workflow did not produce the expected updated sample sheet"
            )

        if manifest is not None:
            manifest.update(
                outputs={"snakemake_outdir": snakemake_outdir},
            )

        vcf_args = copy.deepcopy(args)
        vcf_args.input_csv = updated_csv
        vcf_args.command = "vcf"
        setattr(vcf_args, "_suppress_logo", True)
        vcf_args.test = False
        setattr(vcf_args, "_manifest", manifest)

        return _run_vcf_command(vcf_args)
    except (InputValidationError, ProcessingError) as exc:
        print(f"\nERROR: {exc}\n")
        if manifest is not None:
            manifest.finish(status="failed", error=str(exc))
        return 1
    except Exception as exc:
        print(f"\nUnexpected error: {exc}\n")
        if getattr(args, "debug", False):
            import traceback

            traceback.print_exc()
        if manifest is not None:
            manifest.finish(status="failed", error=str(exc))
        return 1
    finally:
        if test_context is not None:
            test_context.cleanup()


def _run_bam_command(args):
    test_context = None
    manifest = None
    try:
        args = setup_default_paths(args)

        if getattr(args, "test", False):
            test_context, csv_path = _prepare_test_run(args, "bam")
            args.input_csv = str(csv_path.resolve())
            if args.outdir == ".":
                args.outdir = os.path.join(os.getcwd(), "vartracker_bam_test_results")
            if args.snakemake_outdir is None:
                args.snakemake_outdir = args.outdir
        elif args.input_csv is None:
            print(
                "ERROR: Missing required argument: input_csv.\n"
                "Usage: vartracker bam <input_csv> [options]"
            )
            return 1

        os.makedirs(args.outdir, exist_ok=True)
        manifest = _start_manifest(args)
        if manifest is not None:
            _ensure_tool_versions(manifest)
            _update_manifest_inputs(
                manifest,
                input_csv=args.input_csv,
                reference=args.reference,
                gff3=args.gff3,
                input_key="input_csv",
                manifest_level=args.manifest_level,
            )

        try:
            validate_dependencies("bam")
        except DependencyError as error:
            _print_dependency_error(error)
            if manifest is not None:
                manifest.finish(status="failed", error=str(error))
            return 1

        try:
            bam_input = pd.read_csv(args.input_csv)
        except Exception as exc:
            raise InputValidationError(
                f"Could not read input file: {args.input_csv}. {exc}"
            )

        if manifest is not None and args.manifest_level == "deep":
            _update_manifest_inputs(
                manifest,
                input_csv=args.input_csv,
                reference=args.reference,
                gff3=args.gff3,
                input_rows=bam_input.to_dict(orient="records"),
                input_key="input_csv",
                manifest_level=args.manifest_level,
            )

        validate_input_file(
            bam_input,
            optional_empty={"reads1", "reads2", "vcf", "coverage"},
        )

        snakemake_outdir = args.snakemake_outdir or args.outdir

        print(get_logo())

        rulegraph_path = _normalise_rulegraph_path(args.rulegraph)
        args.rulegraph = rulegraph_path

        if rulegraph_path and args.snakemake_dryrun:
            raise InputValidationError(
                "--rulegraph cannot be combined with --snakemake-dryrun"
            )

        updated_csv = run_e2e_workflow(
            samples_csv=args.input_csv,
            reference=args.reference,
            outdir=snakemake_outdir,
            cores=args.cores,
            primer_bed=None,
            dryrun=args.snakemake_dryrun,
            force_all=args.redo,
            quiet=not args.verbose,
            mode="bam",
            rulegraph_path=rulegraph_path,
        )

        if rulegraph_path:
            if manifest is not None:
                manifest.finish(
                    status="success",
                    outputs={"snakemake_rulegraph": rulegraph_path},
                )
            return 0

        if args.snakemake_dryrun:
            if manifest is not None:
                manifest.finish(
                    status="success",
                    outputs={"snakemake_outdir": snakemake_outdir},
                )
            return 0

        if not updated_csv or not os.path.exists(updated_csv):
            raise ProcessingError(
                "BAM workflow did not produce the expected updated sample sheet"
            )

        if manifest is not None:
            manifest.update(
                outputs={"snakemake_outdir": snakemake_outdir},
            )

        vcf_args = copy.deepcopy(args)
        vcf_args.input_csv = updated_csv
        vcf_args.command = "vcf"
        setattr(vcf_args, "_suppress_logo", True)
        vcf_args.test = False
        setattr(vcf_args, "_manifest", manifest)

        return _run_vcf_command(vcf_args)
    except (InputValidationError, ProcessingError) as exc:
        print(f"\nERROR: {exc}\n")
        if manifest is not None:
            manifest.finish(status="failed", error=str(exc))
        return 1
    except Exception as exc:
        print(f"\nUnexpected error: {exc}\n")
        if getattr(args, "debug", False):
            import traceback

            traceback.print_exc()
        if manifest is not None:
            manifest.finish(status="failed", error=str(exc))
        return 1
    finally:
        if test_context is not None:
            test_context.cleanup()


def _run_plot_heatmap_command(args):
    try:
        results_csv = Path(args.results_csv).expanduser().resolve()
        if not results_csv.exists():
            raise InputValidationError(f"Results CSV not found: {results_csv}")

        outdir = (
            Path(args.outdir).expanduser().resolve()
            if args.outdir
            else results_csv.parent
        )
        outdir.mkdir(parents=True, exist_ok=True)

        table = pd.read_csv(results_csv, keep_default_na=False)
        if table.empty:
            raise InputValidationError("Results CSV is empty")
        if "samples" not in table.columns:
            raise InputValidationError(
                "Results CSV must contain a 'samples' column to regenerate the heatmap"
            )

        sample_names = [
            token.strip()
            for token in str(table.iloc[0]["samples"]).split(" / ")
            if token.strip()
        ]
        if not sample_names:
            raise InputValidationError(
                "Could not determine sample names from the results CSV"
            )

        project_name = args.name
        if project_name is None and "name" in table.columns:
            names = [
                str(value).strip()
                for value in table["name"].unique()
                if str(value).strip()
            ]
            if len(names) == 1:
                project_name = names[0]
        if project_name is None:
            project_name = ""

        literature_df = None
        literature_path = None
        if args.literature_csv:
            literature_path = str(Path(args.literature_csv).expanduser().resolve())
            try:
                literature_df = pd.read_csv(literature_path)
            except Exception as exc:
                raise InputValidationError(
                    f"Could not read literature CSV: {exc}"
                ) from exc

        cli_command = getattr(args, "_invocation", None)
        generate_variant_heatmap(
            table,
            sample_names,
            sample_names,
            str(outdir),
            project_name,
            args.min_snv_freq,
            args.min_indel_freq,
            literature_hits=literature_df,
            literature_table_path=literature_path,
            cli_command=cli_command,
            **_collect_heatmap_kwargs(args),
        )
        print(f"\nFinished: find results in {outdir}\n")
        return 0
    except (InputValidationError, ProcessingError) as exc:
        print(f"\nERROR: {exc}\n")
        return 1
    except Exception as exc:
        print(f"\nUnexpected error: {exc}\n")
        if getattr(args, "debug", False):
            import traceback

            traceback.print_exc()
        return 1


def _run_generate_command(args):
    directory = Path(args.dir).expanduser()
    output_path = Path(args.out).expanduser().resolve()

    try:
        records, warnings = generate_csv(
            args.mode,
            directory,
            output_path,
            dry_run=args.dry_run,
        )
    except Exception as exc:
        print(f"\nERROR: {exc}\n")
        return 1

    if args.dry_run:
        buffer = io.StringIO()
        writer = csv.writer(buffer)
        writer.writerow(
            [
                "sample_name",
                "sample_number",
                "reads1",
                "reads2",
                "bam",
                "vcf",
                "coverage",
            ]
        )
        for record in records:
            writer.writerow(record.as_row())

        table = buffer.getvalue().strip()
        if table:
            print(table)
        note = "(dry run)"
    else:
        note = f"written to {output_path}"

    print(f"Generated input spreadsheet {note}")
    if warnings:
        print("\nWarnings:")
        for entry in warnings:
            print(f" - {entry}")

    return 0


def _run_describe_output_command(args):
    if args.schema_target == "results":
        schema = list(get_results_schema())
        schema_text = render_output_schema_text()
    else:
        schema = list(LITERATURE_SCHEMA)
        header = [
            "Literature schema",
            "Expected columns when using `--literature-csv`:",
            "",
        ]
        lines = [
            f"- {entry['name']}: {entry['description']}" for entry in LITERATURE_SCHEMA
        ]
        schema_text = "\n".join(header + lines)

    if args.out:
        output_path = Path(args.out).expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if args.format == "json":
            output_path.write_text(
                json.dumps(schema, indent=2, sort_keys=True), encoding="utf-8"
            )
        else:
            with output_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["name", "type", "description", "units", "values"],
                )
                writer.writeheader()
                writer.writerows(schema)
        print(f"Wrote schema to {output_path}")
        return 0

    print(schema_text)
    return 0


def _process_files(
    tempdir,
    args,
    vcfs,
    sample_names,
    covs,
    input_file,
    literature,
    gene_lengths,
):
    """Process files in the given temporary directory."""
    print("Pre-processing VCF files for compatibility...")

    cli_command = getattr(args, "_invocation", None)

    literature_df = literature

    if args.search_pokay and literature_df is None:
        download_path = os.path.join(args.outdir, "literature_database.csv")
        exit_code = parse_pokay_module.main([download_path])
        if exit_code != 0 or not os.path.exists(download_path):
            raise ProcessingError("Failed to download literature database")
        literature_df = pd.read_csv(download_path)

    skip_heavy_processing = (
        getattr(args, "test", False)
        and os.environ.get("VARTRACKER_SKIP_BCFTOOLS") == "1"
        and getattr(args, "_test_data_dir", None) is not None
    )

    if skip_heavy_processing:
        precomputed_path = (
            Path(args._test_data_dir) / "precomputed" / "test_results.csv"
        )
        if not precomputed_path.exists():
            raise ProcessingError("Precomputed test results not found")
        table = pd.read_csv(precomputed_path, keep_default_na=False)
    else:
        # Format VCF files
        csq_paths = []
        for vcf, sample in zip(vcfs, sample_names):
            _, csq_path = format_vcf(
                vcf,
                sample,
                tempdir,
                args.min_snv_freq,
                args.min_indel_freq,
                args.reference,
                args.gff3,
                args.debug,
                args.allele_frequency_tag,
            )
            csq_paths.append(csq_path)

        # Prepare file lists for merging
        new_vcfs = csq_paths

        # Write VCF list file
        with open(os.path.join(tempdir, "vcf_list.txt"), "w") as f:
            for vcf in new_vcfs:
                f.write(f"{vcf}\n")

        # Write sample names file
        sample_names_file = os.path.join(tempdir, "sample_names.txt")
        with open(sample_names_file, "w") as f:
            for sample_name in sample_names:
                f.write(f"{sample_name}\n")

        # Merge VCF files
        print("Annotating results...")
        csq_file = os.path.join(tempdir, "vcf_annotated.vcf")
        merge_consequences(tempdir, csq_file, sample_names_file, args.debug)

        # Process VCF and extract variants
        print("Summarising results...")
        table = process_vcf(csq_file, covs, args.min_depth, sample_names)

    # Add sample name column if provided
    if args.name is not None:
        table["name"] = args.name
        pname = args.name
    else:
        pname = ""

    table["sample_number"] = " / ".join(
        [str(value) for value in list(input_file["sample_number"])]
    )

    table = _drop_exact_duplicate_result_rows(table)

    # Write initial results
    outfile = os.path.join(args.outdir, args.filename)
    table.fillna("").to_csv(outfile, index=None)

    # Process joint variants
    if not skip_heavy_processing:
        table = process_joint_variants(outfile)
    else:
        table = pd.read_csv(outfile, keep_default_na=False)

    table = _drop_exact_duplicate_result_rows(table)
    table.fillna("").to_csv(outfile, index=None)

    literature_hits_df = None
    literature_full_csv_path = None
    if (
        args.search_pokay or args.literature_csv is not None
    ) and literature_df is not None:
        new_mutations_subset = table[table["variant_status"] == "new"]
        literature_name = args.name if args.name else "sample"
        literature_hits_df = search_literature(
            new_mutations_subset,
            literature_df,
            args.outdir,
            literature_name,
            args.debug,
        )
        literature_full_csv_path = os.path.join(
            args.outdir, f"{literature_name}.literature_database_hits.full.csv"
        )

    # Generate plots and analysis
    print("Plotting results...")
    generate_cumulative_lineplot(
        table,
        pname,
        input_file["sample_number"],
        os.path.join(args.outdir, "cumulative_mutations.pdf"),
    )

    gene_table = generate_gene_table(table, gene_lengths)
    plot_gene_table(gene_table, pname, args.outdir)
    plot_variant_turnover(
        *apply_shared_plot_filters(
            *prepare_plot_inputs(table)[:2],
            include_synonymous=True,
        ),
        output_path=os.path.join(args.outdir, "variant_turnover_plot.pdf"),
        title=f"{pname}: variant turnover" if pname else None,
    )
    try:
        write_reference_feature_metadata(args.gff3, args.outdir)
        plot_variant_genome(
            table,
            load_reference_feature_metadata(outfile),
            output_path=os.path.join(args.outdir, "variant_genome_plot.pdf"),
            title=f"{pname}: genome variant distribution" if pname else None,
        )
    except (InputValidationError, ProcessingError) as exc:
        print(f"Warning: skipped genome plot generation ({exc})")
    generate_variant_heatmap(
        table,
        sample_names,
        list(input_file["sample_number"]),
        args.outdir,
        pname,
        args.min_snv_freq,
        args.min_indel_freq,
        gene_lengths=gene_lengths,
        literature_hits=literature_hits_df,
        literature_table_path=literature_full_csv_path,
        cli_command=cli_command,
    )

    # Write specialized tables
    new_mutations = table[table.presence_absence.str[0] != "Y"][
        ["gene", "variant", "amino_acid_consequence", "nsp_aa_change"]
    ].reset_index(drop=True)
    new_mutations.to_csv(os.path.join(args.outdir, "new_mutations.csv"), index=None)

    persistent_mutations = table[table.persistence_status == "new_persistent"][
        ["gene", "variant", "amino_acid_consequence", "nsp_aa_change"]
    ].reset_index(drop=True)
    persistent_mutations.to_csv(
        os.path.join(args.outdir, "persistent_new_mutations.csv"), index=None
    )


if __name__ == "__main__":
    sys.exit(main())
