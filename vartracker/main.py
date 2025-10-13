"""
Main module for vartracker.

Contains the main function and argument parsing for the vartracker command-line interface.
"""

import argparse
import copy
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
from .core import (
    get_logo,
    get_package_data_path,
    validate_input_file,
    check_dependencies,
    DependencyError,
    InputValidationError,
    ProcessingError,
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
from .data import parse_pokay as parse_pokay_module
from .annotation_processing import (
    gene_lengths_from_gff3,
    validate_reference_and_annotation,
)


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

    _configure_vcf_parser(bam_parser, include_input_csv=True, input_csv_required=True)

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


def create_parser():
    """Create and return the top-level argument parser with subcommands."""

    parser = argparse.ArgumentParser(
        description="vartracker: track the persistence (or not) of mutations during long-term passaging",
        formatter_class=FlexiFormatter,
    )
    parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {__version__}"
    )

    subparsers = parser.add_subparsers(dest="command")
    _add_vcf_subparser(subparsers)
    _add_bam_subparser(subparsers)
    _add_e2e_subparser(subparsers)

    return parser


def validate_dependencies():
    """Check and validate required dependencies."""
    deps = check_dependencies()
    missing = [tool for tool, available in deps.items() if not available]

    if missing:
        raise DependencyError(
            f"Missing required tools: {', '.join(missing)}. "
            "Please install bcftools and tabix."
        )


def setup_default_paths(args):
    """Set up default paths for reference files if not provided."""
    if args.reference is None:
        args.reference = get_package_data_path("NC_045512.fasta")

    annotation = getattr(args, "gff3", None)
    if annotation is None:
        annotation = get_package_data_path("NC_045512.gff3")
    args.gff3 = annotation

    return args


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

    handler = getattr(args, "handler", None)
    if handler is None:
        parser.print_help()
        return 1

    return handler(args)


def _run_vcf_command(args):
    annotation_supplied = getattr(args, "gff3", None) is not None

    with ExitStack() as stack:
        input_csv_path = args.input_csv

        if args.test:
            data_root = resources.files("vartracker").joinpath("test_data")
            temp_dir_obj = tempfile.TemporaryDirectory(prefix="vartracker_demo_data_")
            stack.enter_context(temp_dir_obj)

            target_dir = Path(temp_dir_obj.name)

            if hasattr(data_root, "exists") and data_root.exists():
                with resources.as_file(data_root) as data_path:
                    data_path = Path(data_path)
                    shutil.copytree(data_path, target_dir, dirs_exist_ok=True)
            else:
                fallback = Path(__file__).resolve().parent / "test_data"
                if not fallback.exists():
                    raise InputValidationError("Bundled test data not found")
                shutil.copytree(fallback, target_dir, dirs_exist_ok=True)

            base = target_dir
            args._test_data_dir = str(base)
            input_csv_path = str((base / "inputs" / "test_input.csv").resolve())
            if not args.download_pokay:
                args.pokay_csv = str((base / "mock_pokay" / "mock_pokay.csv").resolve())
            args.search_pokay = True
            if args.name is None:
                args.name = "vartracker_demo"
            if args.outdir == ".":
                args.outdir = os.path.join(os.getcwd(), "vartracker_test_results")
        elif input_csv_path is None:
            args._subparser.print_help()
            return 1

        args.input_csv = input_csv_path

        if not getattr(args, "_suppress_logo", False):
            print(get_logo())

        try:
            # Validate dependencies
            validate_dependencies()

            # Set up default paths
            args = setup_default_paths(args)

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

            # Create output directory
            os.makedirs(args.outdir, exist_ok=True)

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

            print(f"\nFinished: find results in {args.outdir}\n")
            return 0
        except (DependencyError, InputValidationError, ProcessingError) as e:
            print(f"\nERROR: {str(e)}\n")
            return 1
        except Exception as e:
            print(f"\nUnexpected error: {str(e)}\n")
            if args.debug:
                import traceback

                traceback.print_exc()
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
    args = setup_default_paths(args)

    snakemake_outdir = args.snakemake_outdir or args.outdir

    try:
        snakemake_input = pd.read_csv(args.samples_csv)
    except Exception as exc:
        raise InputValidationError(
            f"Could not read input file: {args.samples_csv}. {exc}"
        )

    validate_input_file(
        snakemake_input,
        optional_empty={"bam", "vcf", "coverage", "reads2"},
    )

    print(get_logo())

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
    )

    if args.snakemake_dryrun:
        return 0

    if not updated_csv or not os.path.exists(updated_csv):
        raise ProcessingError(
            "End-to-end workflow did not produce the expected updated sample sheet"
        )

    vcf_args = copy.deepcopy(args)
    vcf_args.input_csv = updated_csv
    vcf_args.command = "vcf"
    setattr(vcf_args, "_suppress_logo", True)

    return _run_vcf_command(vcf_args)


def _run_bam_command(args):
    args = setup_default_paths(args)

    if args.input_csv is None:
        print(
            "ERROR: Missing required argument: input_csv.\n"
            "Usage: vartracker bam <input_csv> [options]"
        )
        return 1

    try:
        bam_input = pd.read_csv(args.input_csv)
    except Exception as exc:
        raise InputValidationError(
            f"Could not read input file: {args.input_csv}. {exc}"
        )

    validate_input_file(
        bam_input,
        optional_empty={"reads1", "reads2", "vcf", "coverage"},
    )

    snakemake_outdir = args.snakemake_outdir or args.outdir

    print(get_logo())

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
    )

    if args.snakemake_dryrun:
        return 0

    if not updated_csv or not os.path.exists(updated_csv):
        raise ProcessingError(
            "BAM workflow did not produce the expected updated sample sheet"
        )

    vcf_args = copy.deepcopy(args)
    vcf_args.input_csv = updated_csv
    vcf_args.command = "vcf"
    setattr(vcf_args, "_suppress_logo", True)

    return _run_vcf_command(vcf_args)


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

    # Search pokay database
    if args.search_pokay and pokay_df is not None:
        new_mutations_table = table[table["variant_status"] == "new"]
        pokay_name = args.name if args.name else "sample"
        search_pokay(new_mutations_table, pokay_df, args.outdir, pokay_name, args.debug)


if __name__ == "__main__":
    sys.exit(main())
