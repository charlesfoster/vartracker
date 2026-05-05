#!/usr/bin/env python3
"""Launch the Snakemake workflow via the Snakemake Python API."""

from __future__ import annotations

import sys
import io
from contextlib import redirect_stdout
from pathlib import Path
from typing import Optional

from snakemake.api import (
    ConfigSettings,
    ResourceSettings,
    SnakemakeApi,
    DAGSettings,
    ExecutionSettings,
    OutputSettings,
)
from snakemake.settings.enums import PrintDag, Quietness


def _normalise_path(value: str | Path) -> str:
    return str(Path(value).expanduser().resolve())


def _first_fasta_record_id(fasta_path: str | Path) -> str | None:
    with Path(fasta_path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                return line[1:].split()[0].strip()
    return None


def _primer_bed_contigs(primer_bed_path: str | Path) -> set[str]:
    contigs = set()
    with Path(primer_bed_path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("track ") or line.startswith("browser "):
                continue
            contigs.add(line.split()[0])
    return contigs


def _validate_primer_bed_reference(
    primer_bed_path: str | Path, reference_path: str | Path
) -> None:
    reference_id = _first_fasta_record_id(reference_path)
    if not reference_id:
        raise ValueError(
            f"Unable to determine reference FASTA ID from: {reference_path}"
        )

    bed_contigs = _primer_bed_contigs(primer_bed_path)
    if not bed_contigs:
        raise ValueError(f"Primer BED file contains no intervals: {primer_bed_path}")

    mismatched = sorted(contig for contig in bed_contigs if contig != reference_id)
    if mismatched:
        joined = ", ".join(mismatched)
        raise ValueError(
            f"Primer BED contig(s) {joined!r} do not match reference FASTA ID "
            f"{reference_id!r}"
        )


def _validate_lofreq_primer_rescue(
    mode: str,
    primer_bed_path: Optional[str],
    min_af: float,
    min_dp: int,
    min_alt_count: int,
    min_qual: float,
    max_ref_count: int,
) -> None:
    if mode not in {"auto", "on", "off"}:
        raise ValueError("lofreq_primer_rescue must be one of: auto, on, off")
    if mode == "on" and not primer_bed_path:
        raise ValueError("--lofreq-primer-rescue on requires --primer-bed")
    if not 0 <= min_af <= 1:
        raise ValueError("lofreq_rescue_min_af must be between 0 and 1")
    if min_dp < 0:
        raise ValueError("lofreq_rescue_min_dp must be >= 0")
    if min_alt_count < 0:
        raise ValueError("lofreq_rescue_min_alt_count must be >= 0")
    if min_qual < 0:
        raise ValueError("lofreq_rescue_min_qual must be >= 0")
    if max_ref_count < 0:
        raise ValueError("lofreq_rescue_max_ref_count must be >= 0")


def run_workflow(
    samples_csv: str | Path,
    reference: str | Path,
    outdir: str | Path = "results",
    cores: int = 8,
    primer_bed: Optional[str | Path] = None,
    ampliconclip_tolerance: int = 1,
    min_depth: int = 10,
    consensus_snp_min_af: float = 0.25,
    consensus_snp_thresh: float = 0.75,
    consensus_indel_thresh: float = 0.75,
    dryrun: bool = False,
    force_all: bool = False,
    quiet: bool = True,
    mode: str = "reads",
    rulegraph_path: Optional[str | Path] = None,
    lofreq_primer_rescue: str = "auto",
    lofreq_rescue_min_af: float = 0.95,
    lofreq_rescue_min_dp: int = 100,
    lofreq_rescue_min_alt_count: int = 95,
    lofreq_rescue_min_qual: float = 100.0,
    lofreq_rescue_max_ref_count: int = 20,
) -> Optional[str]:
    """Run the lofreq variant calling workflow via the Snakemake API.

    Returns the path to the updated sample CSV when the workflow runs, or
    ``None`` if a dry run or rulegraph-only request was made.
    """

    samples_csv = _normalise_path(samples_csv)
    reference = _normalise_path(reference)
    outdir = _normalise_path(outdir)
    primer_bed_path = _normalise_path(primer_bed) if primer_bed else None
    rulegraph_path = _normalise_path(rulegraph_path) if rulegraph_path else None
    if rulegraph_path and not rulegraph_path.lower().endswith(".dot"):
        rulegraph_path = f"{rulegraph_path}.dot"

    # Validate inputs
    if not Path(samples_csv).exists():
        raise FileNotFoundError(f"Sample CSV not found: {samples_csv}")
    if not Path(reference).exists():
        raise FileNotFoundError(f"Reference genome not found: {reference}")
    if primer_bed_path and not Path(primer_bed_path).exists():
        raise FileNotFoundError(f"Primer BED file not found: {primer_bed_path}")
    if primer_bed_path:
        _validate_primer_bed_reference(primer_bed_path, reference)
    _validate_lofreq_primer_rescue(
        lofreq_primer_rescue,
        primer_bed_path,
        lofreq_rescue_min_af,
        lofreq_rescue_min_dp,
        lofreq_rescue_min_alt_count,
        lofreq_rescue_min_qual,
        lofreq_rescue_max_ref_count,
    )

    Path(outdir).mkdir(parents=True, exist_ok=True)

    # Build config for the workflow (paths must be strings)
    config_dict = {
        "samples_csv": samples_csv,
        "reference": reference,
        "outdir": outdir,
        "mode": mode,
        "ampliconclip_tolerance": ampliconclip_tolerance,
        "min_depth": min_depth,
        "consensus_snp_min_af": consensus_snp_min_af,
        "consensus_snp_thresh": consensus_snp_thresh,
        "consensus_indel_thresh": consensus_indel_thresh,
        "lofreq_primer_rescue": lofreq_primer_rescue,
        "lofreq_rescue_min_af": lofreq_rescue_min_af,
        "lofreq_rescue_min_dp": lofreq_rescue_min_dp,
        "lofreq_rescue_min_alt_count": lofreq_rescue_min_alt_count,
        "lofreq_rescue_min_qual": lofreq_rescue_min_qual,
        "lofreq_rescue_max_ref_count": lofreq_rescue_max_ref_count,
    }
    if primer_bed_path:
        config_dict["primer_bed"] = primer_bed_path

    snakefile = Path(__file__).resolve().parent / "Snakefile"
    if not snakefile.exists():
        raise FileNotFoundError(f"Snakefile not found at: {snakefile}")

    workdir = snakefile.parent

    # Snakemake expects an iterable of Quietness enums (or None),
    # not a bare bool. Convert our boolean to the canonical set.
    quiet_setting = frozenset({Quietness.ALL}) if quiet else None

    if rulegraph_path:
        print("Rendering Snakemake rulegraph...")
    else:
        print("Deploying analysis with snakemake...")

    try:
        with SnakemakeApi(
            OutputSettings(
                dryrun=dryrun,
                quiet=quiet_setting,
                printshellcmds=not quiet,
            )
        ) as snakemake_api:
            workflow = snakemake_api.workflow(
                resource_settings=ResourceSettings(cores=cores),
                config_settings=ConfigSettings(config=config_dict),
                snakefile=snakefile,
                workdir=workdir,
            )

            dag = workflow.dag(
                DAGSettings(forceall=force_all, print_dag_as=str(PrintDag.DOT))
            )

            if rulegraph_path:
                output_buffer = io.StringIO()
                with redirect_stdout(output_buffer):
                    dag.printrulegraph()
                Path(rulegraph_path).parent.mkdir(parents=True, exist_ok=True)
                Path(rulegraph_path).write_text(
                    output_buffer.getvalue(), encoding="utf-8"
                )
                print(f"✓ Rulegraph written to: {rulegraph_path}")
                return None

            executor = "dryrun" if dryrun else "local"
            execution_settings = ExecutionSettings()
            dag.execute_workflow(
                executor=executor, execution_settings=execution_settings
            )

    except Exception as exc:  # pragma: no cover - surfaced to CLI
        raise RuntimeError(f"Snakemake execution failed: {exc}") from exc

    if dryrun:
        print("\n✓ Dry run completed. No jobs were executed.")
        return None

    print("✓ Snakemake workflow completed successfully!")
    print(f"Results are in: {outdir}/")
    updated_csv = str(Path(outdir) / "vartracker_execution_spreadsheet.csv")
    print(f"Updated sample sheet: {updated_csv}\n")
    print("Running vartracker summary workflow...")
    return updated_csv


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run lofreq variant calling pipeline")
    parser.add_argument("--samples", required=True, help="Path to samples CSV")
    parser.add_argument("--reference", required=True, help="Path to reference genome")
    parser.add_argument("--outdir", default="results", help="Output directory")
    parser.add_argument("--cores", type=int, default=8, help="Number of cores")
    parser.add_argument(
        "--primer-bed", help="Optional primer BED file for amplicon clipping"
    )
    parser.add_argument(
        "--lofreq-primer-rescue",
        choices=("auto", "on", "off"),
        default="auto",
        help="Primer-overlap rescue mode for LoFreq calls (default: auto)",
    )
    parser.add_argument("--dryrun", action="store_true", help="Perform dry run")

    args = parser.parse_args()

    try:
        run_workflow(
            samples_csv=args.samples,
            reference=args.reference,
            outdir=args.outdir,
            cores=args.cores,
            primer_bed=args.primer_bed,
            lofreq_primer_rescue=args.lofreq_primer_rescue,
            dryrun=args.dryrun,
        )
    except Exception as exc:  # pragma: no cover - CLI surface
        print(f"\n✗ Workflow failed: {exc}", file=sys.stderr)
        sys.exit(1)
