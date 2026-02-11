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

import pandas as pd
from argparse_formatter import FlexiFormatter

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
    search_pokay,
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

_RED = "\033[91m"
_YELLOW = "\033[93m"
_RESET = "\033[0m"


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

    pokay_path = base / "mock_pokay" / "mock_pokay.csv"
    if not pokay_path.exists():
        raise InputValidationError("Bundled mock pokay database not found")

    args.search_pokay = True
    args.download_pokay = False
    args.pokay_csv = str(pokay_path)
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

    vt_group = parser.add_argument_group("Vartracker options")

    vt_group.add_argument(
        "-g",
        "--gff3",
        action="store",
        required=False,
        default=None,
        help="GFF3 annotations to use (default: packaged SARS-CoV-2 annotations)",
    )
    vt_group.add_argument(
        "-m",
        "--min-snv-freq",
        action="store",
        required=False,
        type=float,
        default=0.03,
        help="Minimum allele frequency of SNV variants to keep (default: 0.03)",
    )
    vt_group.add_argument(
        "-M",
        "--min-indel-freq",
        action="store",
        required=False,
        type=float,
        default=0.1,
        help="Minimum allele frequency of indel variants to keep (default: 0.1)",
    )
    vt_group.add_argument(
        "-d",
        "--min-depth",
        action="store",
        required=False,
        type=int,
        default=10,
        help="Minimum depth threshold for variant QC (default: 10)",
    )
    vt_group.add_argument(
        "-n",
        "--name",
        action="store",
        required=False,
        default=None,
        help="Optional: add a column to results with the name specified here",
    )
    vt_group.add_argument(
        "-o",
        "--outdir",
        action="store",
        required=False,
        default=".",
        help="Output directory for vartracker results (default: current directory)",
    )
    vt_group.add_argument(
        "--manifest-level",
        choices=["light", "deep"],
        default="light",
        help="Manifest detail level for run metadata (default: light)",
    )
    vt_group.add_argument(
        "-f",
        "--filename",
        action="store",
        required=False,
        default="results.csv",
        help="Output file name (default: results.csv)",
    )
    vt_group.add_argument(
        "-r",
        "--reference",
        action="store",
        required=False,
        default=None,
        help="Reference genome (default: uses packaged SARS-CoV-2 reference)",
    )
    vt_group.add_argument(
        "--passage-cap",
        action="store",
        type=int,
        help="Cap the number of passages at this number",
        default=None,
    )
    vt_group.add_argument(
        "--debug",
        action="store_true",
        help="Print commands being run for debugging",
        default=False,
    )
    vt_group.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep temporary files for debugging",
        default=False,
    )
    vt_group.add_argument(
        "--pokay-csv",
        action="store",
        required=False,
        default=None,
        help="Path to a pre-parsed pokay CSV file",
    )
    vt_group.add_argument(
        "--search-pokay",
        action="store_true",
        help="Run literature search against the pokay database.",
        default=False,
    )
    vt_group.add_argument(
        "--download-pokay",
        action="store_true",
        help="Automatically download and parse the pokay literature database.",
        default=False,
    )
    vt_group.add_argument(
        "--test",
        action="store_true",
        help="Run vartracker against the bundled demonstration dataset",
        default=False,
    )
    vt_group.add_argument(
        "--allele-frequency-tag",
        action="store",
        required=False,
        default="AF",
        help="INFO tag name for allele frequency (default: AF)",
    )


def _add_vcf_subparser(subparsers):
    description = (
        "Analyse VCF-based longitudinal variant data, generating reports and plots."
    )
    vcf_parser = subparsers.add_parser(
        "vcf",
        help="Analyse VCF inputs",
        description=description,
        formatter_class=FlexiFormatter,
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
        formatter_class=FlexiFormatter,
        epilog="""
Provide a CSV with columns:

sample_name,sample_number,reads1,reads2,bam,vcf,coverage

In BAM mode the `bam` column must point to existing files while `reads1`,
`reads2`, `vcf`, and `coverage` may be blank.
""",
    )

    _configure_vcf_parser(bam_parser, include_input_csv=True, input_csv_required=False)

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

    bam_parser.set_defaults(handler=_run_bam_command, _subparser=bam_parser)


def _add_spreadsheet_subparser(subparsers):
    description = "Generate template CSV files for vartracker input."
    gen_parser = subparsers.add_parser(
        "spreadsheet",
        help="Generate input spreadsheets from an existing directory of files",
        description=description,
        formatter_class=FlexiFormatter,
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


def _run_prep_command(args):
    if getattr(args, "handler", None) is None or args.command == "prep":
        args._subparser.print_help()
        return 1
    return args.handler(args)


def _add_prep_subparser(subparsers):
    description = "Prepare inputs and templates for vartracker workflows."
    prep_parser = subparsers.add_parser(
        "prep",
        help="Prepare inputs for vartracker",
        description=description,
        formatter_class=FlexiFormatter,
    )
    prep_subparsers = prep_parser.add_subparsers(dest="prep_command")
    _add_spreadsheet_subparser(prep_subparsers)
    prep_parser.set_defaults(handler=_run_prep_command, _subparser=prep_parser)


def _add_schema_subparser(subparsers):
    description = "Describe the output schema for vartracker results."
    schema_parser = subparsers.add_parser(
        "schema",
        help="Print the output schema for results tables",
        description=description,
        aliases=["describe-output"],
        formatter_class=FlexiFormatter,
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


def create_parser():
    """Create and return the top-level argument parser with subcommands."""

    parser = argparse.ArgumentParser(
        description="vartracker: longitudinal variant tracking and summarisation for pathogen sequencing",
        formatter_class=FlexiFormatter,
    )
    parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {__version__}"
    )

    subparsers = parser.add_subparsers(dest="command")
    _add_vcf_subparser(subparsers)
    _add_bam_subparser(subparsers)
    _add_e2e_subparser(subparsers)
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

            # Apply passage cap if specified
            if args.passage_cap is not None:
                input_file = input_file[input_file["sample_number"] <= args.passage_cap]

            pokay = None
            if args.search_pokay:
                if args.pokay_csv is not None:
                    try:
                        pokay = pd.read_csv(args.pokay_csv)
                    except Exception as e:
                        raise InputValidationError(
                            f"Could not read pokay database: {str(e)}"
                        )
                elif not args.download_pokay:
                    raise InputValidationError(
                        "--search-pokay requires either --pokay-csv or --download-pokay"
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
                        pokay,
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
                        pokay,
                        gene_lengths_override,
                    )

            results_metadata_path = write_results_metadata(
                args.outdir, results_filename=args.filename
            )

            outputs = {
                "results_csv": os.path.join(args.outdir, args.filename),
                "results_metadata": str(results_metadata_path),
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
                "variant_allele_frequency_heatmap_html": os.path.join(
                    args.outdir, "variant_allele_frequency_heatmap.html"
                ),
                "variant_allele_frequency_heatmap_pdf": os.path.join(
                    args.outdir, "variant_allele_frequency_heatmap.pdf"
                ),
            }

            if args.search_pokay:
                pokay_name = args.name if args.name else "sample"
                pokay_full = os.path.join(
                    args.outdir, f"{pokay_name}.pokay_database_hits.full.csv"
                )
                pokay_concise = os.path.join(
                    args.outdir, f"{pokay_name}.pokay_database_hits.concise.csv"
                )
                if os.path.exists(pokay_full):
                    outputs["pokay_hits_full_csv"] = pokay_full
                if os.path.exists(pokay_concise):
                    outputs["pokay_hits_concise_csv"] = pokay_concise

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
        formatter_class=FlexiFormatter,
        aliases=["e2e"],
    )

    _configure_vcf_parser(e2e_parser, include_input_csv=False)

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
    if args.out:
        output_path = Path(args.out).expanduser().resolve()
        schema = list(get_results_schema())
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

    print(render_output_schema_text())
    return 0


def _process_files(
    tempdir,
    args,
    vcfs,
    sample_names,
    covs,
    input_file,
    pokay,
    gene_lengths,
):
    """Process files in the given temporary directory."""
    print("Pre-processing VCF files for compatibility...")

    cli_command = getattr(args, "_invocation", None)

    pokay_df = pokay

    if args.search_pokay and pokay_df is None and args.download_pokay:
        download_path = os.path.join(tempdir, "pokay_database.csv")
        exit_code = parse_pokay_module.main([download_path])
        if exit_code != 0 or not os.path.exists(download_path):
            raise ProcessingError("Failed to download pokay literature database")
        pokay_df = pd.read_csv(download_path)

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

    # Write initial results
    outfile = os.path.join(args.outdir, args.filename)
    table.fillna("").to_csv(outfile, index=None)

    # Process joint variants
    if not skip_heavy_processing:
        table = process_joint_variants(outfile)
    else:
        table = pd.read_csv(outfile, keep_default_na=False)

    pokay_hits_df = None
    pokay_full_csv_path = None
    if args.search_pokay and pokay_df is not None:
        new_mutations_subset = table[table["variant_status"] == "new"]
        pokay_name = args.name if args.name else "sample"
        pokay_hits_df = search_pokay(
            new_mutations_subset, pokay_df, args.outdir, pokay_name, args.debug
        )
        pokay_full_csv_path = os.path.join(
            args.outdir, f"{pokay_name}.pokay_database_hits.full.csv"
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
    generate_variant_heatmap(
        table,
        sample_names,
        list(input_file["sample_number"]),
        args.outdir,
        pname,
        args.min_snv_freq,
        args.min_indel_freq,
        gene_lengths=gene_lengths,
        pokay_hits=pokay_hits_df,
        pokay_table_path=pokay_full_csv_path,
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
