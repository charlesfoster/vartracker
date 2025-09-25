"""
Main module for vartracker.

Contains the main function and argument parsing for the vartracker command-line interface.
"""

import argparse
import os
import shutil
import sys
import tempfile

import pandas as pd
from argparse_formatter import FlexiFormatter

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


def create_parser():
    """Create and return the argument parser."""
    parser = argparse.ArgumentParser(
        description="vartracker: track the persistence (or not) of mutations during long-term passaging",
        usage="vartracker [options] <input_csv>",
        formatter_class=FlexiFormatter,
        epilog="""
    The input CSV file should have four columns:
    1. 'vcf': full paths to bgzipped vcf files containing the called variants for each sample of interest
    2. 'coverage': full paths to files containing the genome coverage for each sample. The coverage files should be in the format output by `bedtools genomecov` (e.g., bedtools genomecov -ibam input.bam -d > coverage.txt), whereby the columns are (in order): the reference name, the 1-based position in the genome, and the depth of coverage at that position.
    3. 'sample_name': the name of the sample in the given VCF file
    4. 'sample_number': the sample number. In a longitudinal comparison like long-term passaging, the numbers in the column might go 0, 1, 2, ..., 15.

    Note: the order of values in the input CSV file matters, dictating the order of results in the output CSV file.

    If 'new mutations' are to be searched against the functional 'pokay' database (default) (see vartracker repository README), the path to
    the parsed pokay csv file must be provided.
    """,
    )

    parser.add_argument("input_csv", nargs=1, help="Input CSV file. See below.")
    parser.add_argument(
        "-a",
        "--annotation",
        action="store",
        required=False,
        default=None,  # Will be set to package data if None
        help="Annotations to use in GFF3 format (default: uses packaged SARS-CoV-2 annotations)",
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
    parser.add_argument(
        "-n",
        "--name",
        action="store",
        required=False,
        default=None,
        help="Optional: add a column to results with the name specified here",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        required=False,
        default=".",
        help="Output directory (default: current directory)",
    )
    parser.add_argument(
        "-f",
        "--filename",
        action="store",
        required=False,
        default="results.csv",
        help="Output file name (default: results.csv)",
    )
    parser.add_argument(
        "-r",
        "--reference",
        action="store",
        required=False,
        default=None,  # Will be set to package data if None
        help="Reference genome (default: uses packaged SARS-CoV-2 reference)",
    )
    parser.add_argument(
        "--passage-cap",
        action="store",
        type=int,
        help="Cap the number of passages at this number",
        default=None,
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print commands being run for debugging",
        default=False,
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep temporary files for debugging",
        default=False,
    )
    parser.add_argument(
        "--pokay-csv",
        action="store",
        required=False,
        default=None,
        help="Path to a pre-parsed pokay CSV file",
    )
    parser.add_argument(
        "--search-pokay",
        action="store_true",
        help="Run literature search against the pokay database.",
        default=False,
    )
    parser.add_argument(
        "--download-pokay",
        action="store_true",
        help="Automatically download and parse the pokay literature database.",
        default=False,
    )
    parser.add_argument(
        "--allele-frequency-tag",
        action="store",
        required=False,
        default="AF",
        help="INFO tag name for allele frequency (default: AF)",
    )
    parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {__version__}"
    )

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

    if args.annotation is None:
        args.annotation = get_package_data_path("NC_045512.gff3")

    return args


def main(sysargs=None):
    """Main function for vartracker."""
    if sysargs is None:
        sysargs = sys.argv[1:]

    parser = create_parser()

    if len(sysargs) < 1:
        parser.print_help()
        return 1

    try:
        args = parser.parse_args(sysargs)
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else 0
        return code

    print(get_logo())

    try:
        # Validate dependencies
        validate_dependencies()

        # Set up default paths
        args = setup_default_paths(args)

        # Read input CSV
        try:
            input_file = pd.read_csv(args.input_csv[0])
        except Exception as e:
            raise InputValidationError(
                f"Could not read input file: {args.input_csv[0]}. {str(e)}"
            )

        # Resolve potential relative paths for VCF and coverage files
        csv_dir = os.path.dirname(os.path.abspath(args.input_csv[0]))

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
        validate_input_file(input_file)

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
                    tempdir, args, vcfs, sample_names, covs, input_file, pokay
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
                    tempdir, args, vcfs, sample_names, covs, input_file, pokay
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


def _process_files(tempdir, args, vcfs, sample_names, covs, input_file, pokay):
    """Process files in the given temporary directory."""
    print("Pre-processing VCF files for compatibility...")

    pokay_df = pokay

    if args.search_pokay and pokay_df is None and args.download_pokay:
        download_path = os.path.join(tempdir, "pokay_database.csv")
        exit_code = parse_pokay_module.main([download_path])
        if exit_code != 0 or not os.path.exists(download_path):
            raise ProcessingError("Failed to download pokay literature database")
        pokay_df = pd.read_csv(download_path)

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
            args.annotation,
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
    table = process_vcf(csq_file, covs)

    # Add sample name column if provided
    if args.name is not None:
        table["name"] = args.name
        pname = args.name
    else:
        pname = ""

    # Write initial results
    outfile = os.path.join(args.outdir, args.filename)
    table.to_csv(outfile, index=None)

    # Process joint variants
    table = process_joint_variants(outfile)

    # Generate plots and analysis
    print("Plotting results...")
    generate_cumulative_lineplot(
        table,
        pname,
        input_file["sample_number"],
        os.path.join(args.outdir, "cumulative_mutations.pdf"),
    )

    gene_table = generate_gene_table(table)
    plot_gene_table(gene_table, pname, args.outdir)
    generate_variant_heatmap(
        table,
        sample_names,
        list(input_file["sample_number"]),
        args.outdir,
        pname,
        args.min_snv_freq,
        args.min_indel_freq,
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
