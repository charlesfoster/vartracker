"""
VCF file processing functionality for vartracker.

Contains functions for formatting, merging, and processing VCF files.
"""

import os
import re
import subprocess
import gzip
import shutil
from pathlib import Path
from typing import TypedDict

import pandas as pd
import numpy as np
from cyvcf2 import VCF, Writer
from itertools import dropwhile

from .constants import reformat_csq_notation
from .amino_acids import AminoAcidChange


class _AltEntry(TypedDict):
    index: int
    alt: str
    af: float


def _open_vcf(path: str, mode: str):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def _ensure_format_and_sample(
    vcf_path: str, sample: str, tempdir: str, debug: bool = False
) -> str:
    """Ensure VCF has FORMAT header and a sample column; return path to updated VCF."""

    has_format_header = False
    with _open_vcf(vcf_path, "rt") as src:
        for line in src:
            if line.startswith("##FORMAT="):
                has_format_header = True
                break
            if line.startswith("#CHROM"):
                break

    if has_format_header:
        return vcf_path

    if debug:
        print("Formatting VCF to add GT field and sample column...")

    intermediate_path = os.path.join(
        tempdir, f"{Path(vcf_path).stem}.with_gt_sample.vcf"
    )

    with (
        _open_vcf(vcf_path, "rt") as src,
        open(intermediate_path, "w", encoding="utf-8") as dest,
    ):
        format_added = False
        header_written = False

        for line in src:
            if line.startswith("##"):
                if line.startswith("##FORMAT="):
                    format_added = True
                dest.write(line)
                continue

            if line.startswith("#CHROM"):
                if not format_added:
                    dest.write(
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                    )
                    format_added = True

                header_cols = line.rstrip("\n").split("\t")
                if len(header_cols) <= 8:
                    header_cols.extend(["FORMAT", sample])
                elif header_cols[-1] == "FORMAT":
                    header_cols.append(sample)
                else:
                    header_cols.extend(["FORMAT", sample])
                dest.write("\t".join(header_cols) + "\n")
                header_written = True
                break

            dest.write(line)

        if not header_written:
            raise RuntimeError("VCF header missing '#CHROM' line")

        for data_line in src:
            if not data_line.strip() or data_line.startswith("#"):
                continue
            parts = data_line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            base_fields = parts[:8]
            dest.write("\t".join(base_fields + ["GT", "1"]) + "\n")

    sample_file = os.path.join(tempdir, f"{Path(vcf_path).stem}.sample_names.txt")
    with open(sample_file, "w", encoding="utf-8") as handle:
        handle.write(f"{sample}\n")

    reheadered_path = os.path.join(tempdir, f"{Path(vcf_path).stem}.with_sample.vcf")

    cmd = f"bcftools reheader -s {sample_file} {intermediate_path} > {reheadered_path}"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"Failed to reheader VCF to add sample information: {exc}"
        ) from exc

    return reheadered_path


def _derive_vcf_output_paths(
    vcf_path: str, tempdir: str, sample: str
) -> tuple[str, str, str]:
    """Return paths for cyvcf2 intermediate, csq output, and log files."""

    basename = os.path.basename(vcf_path)
    prefix = re.sub(r"\.vcf(?:\.gz|\.bgz)?$", "", basename, flags=re.IGNORECASE)
    if prefix == basename:
        # Fallback when extension doesn't match expected pattern
        prefix = basename

    out = os.path.join(tempdir, f"{prefix}.cyvcf2.vcf.gz")
    csq_file = os.path.join(tempdir, f"{prefix}.csq.vcf.gz")
    log_file = os.path.join(tempdir, f"{sample}.log")

    return out, csq_file, log_file


def _normalize_header_line(line: str) -> str:
    """Adjust known header cardinalities to satisfy bcftools sanity checks."""

    if line.startswith("##INFO=<ID=AF,") or line.startswith("##FORMAT=<ID=AF,"):
        return re.sub(r"Number=[^,>]+", "Number=A", line, count=1)
    if line.startswith("##INFO=<ID=SB,") or line.startswith("##FORMAT=<ID=SB,"):
        return re.sub(r"Number=[^,>]+", "Number=4", line, count=1)
    return line


def _run_logged_command(cmd: str, log_file: str, context: str) -> None:
    """Run a shell command and surface the tail of the log on failure."""

    try:
        with open(log_file, "a", encoding="utf-8") as err:
            subprocess.run(cmd, shell=True, stderr=err, check=True)
    except subprocess.CalledProcessError as exc:
        log_tail = ""
        if os.path.exists(log_file):
            try:
                with open(
                    log_file, "r", encoding="utf-8", errors="replace"
                ) as err_file:
                    tail_lines = err_file.readlines()[-20:]
                    log_tail = "".join(tail_lines).strip()
            except OSError:
                log_tail = "<unable to read bcftools log>"

        details = log_tail if log_tail else exc.stderr or ""
        message = (
            f"{context} failed with exit status {exc.returncode} for command '{cmd}'."
        )
        if details:
            message = f"{message}\nLast bcftools log lines:\n{details}"
        raise RuntimeError(message) from exc


def _split_multiallelic_records(
    vcf_path: str, tempdir: str, sample: str, log_file: str, debug: bool = False
) -> str:
    """Split multi-ALT records so downstream processing can track each allele."""

    split_path = os.path.join(tempdir, f"{Path(vcf_path).stem}.{sample}.split.vcf")
    cmd = f"bcftools norm -m -any -Ov -o {split_path} {vcf_path}"
    _run_logged_command(cmd, log_file, "bcftools norm")

    if debug:
        print(f"Command: {cmd}")

    return split_path


def rindex(lst, item):
    """Find the last occurrence of an item in a list."""

    def index_ne(x):
        return lst[x] != item

    try:
        return next(dropwhile(index_ne, reversed(range(len(lst)))))
    except StopIteration:
        raise ValueError("rindex(lst, item): item not in list")


def format_vcf(
    vcf,
    sample,
    tempdir,
    min_snv_freq,
    min_indel_freq,
    reference,
    annotation,
    debug,
    allele_frequency_tag="AF",
):
    """
    Format VCF file for compatibility and filter variants before merging.

    Args:
        vcf (str): Path to input VCF file
        sample (str): Sample name
        tempdir (str): Temporary directory path
        min_snv_freq (float): Minimum SNV frequency threshold
        min_indel_freq (float): Minimum indel frequency threshold
        reference (str): Path to reference genome
        annotation (str): Path to annotation file
        debug (bool): Whether to print debug information
        allele_frequency_tag (str): INFO tag name for allele frequency (default: AF)
    """
    out, csq_file, log = _derive_vcf_output_paths(vcf, tempdir, sample)
    raw_out = os.path.join(tempdir, f"{Path(out).stem}.raw.vcf.gz")

    try:
        prepared_vcf = _ensure_format_and_sample(vcf, sample, tempdir, debug)
        split_vcf = _split_multiallelic_records(
            prepared_vcf, tempdir, sample, log, debug
        )

        vcf_mod = VCF(split_vcf, strict_gt=True)
        existing_samples = list(vcf_mod.samples)
        sample_count = len(existing_samples) if existing_samples else 1

        # Validate that the allele frequency tag exists in the VCF
        info_headers = {
            header["ID"]
            for header in vcf_mod.header_iter()
            if header["HeaderType"] == "INFO"
        }

        if allele_frequency_tag not in info_headers:
            raise RuntimeError(
                f"Allele frequency tag '{allele_frequency_tag}' not found in VCF INFO headers. "
                f"Available INFO tags: {', '.join(sorted(info_headers))}"
            )

        # Add INFO headers for DP and AF if they don't exist
        if "DP" not in info_headers:
            vcf_mod.add_info_to_header(
                {
                    "ID": "DP",
                    "Description": "Raw Depth",
                    "Type": "Integer",
                    "Number": "1",
                }
            )

        if "AF" not in info_headers:
            vcf_mod.add_info_to_header(
                {
                    "ID": "AF",
                    "Description": "Allele Frequency",
                    "Type": "Float",
                    "Number": "A",
                }
            )

        if "TYPE" not in info_headers:
            vcf_mod.add_info_to_header(
                {
                    "ID": "TYPE",
                    "Description": "Variant type inferred by vartracker",
                    "Type": "String",
                    "Number": "1",
                }
            )

        # Add FORMAT headers
        vcf_mod.add_format_to_header(
            {"ID": "DP", "Description": "Raw Depth", "Type": "Integer", "Number": "1"}
        )
        vcf_mod.add_format_to_header(
            {
                "ID": "AF",
                "Description": "Allele Frequency",
                "Type": "Float",
                "Number": "A",
            }
        )

        header_lines = []
        for line in vcf_mod.raw_header.strip().splitlines():
            normalized_line = _normalize_header_line(line)
            if line.startswith("#CHROM") and not existing_samples:
                header_lines.append(f"{normalized_line}\t{sample}")
            else:
                header_lines.append(normalized_line)
        header_str = "\n".join(header_lines) + "\n"
        w = Writer.from_string(raw_out, header_str, mode="wz")
        variants = {}

        def _scalar(value):
            if isinstance(value, (list, tuple, np.ndarray)):
                return value[0] if value else 0
            return value

        for v in vcf_mod:
            if not v.ALT:
                continue

            variant_key = (v.CHROM, int(v.POS), v.REF, v.ALT[0])

            # Get allele frequency from the specified tag
            if allele_frequency_tag != "AF":
                # Rename the allele frequency tag to AF for internal processing
                af_value = v.INFO[allele_frequency_tag]
                # Handle floating point precision issues by rounding to 6 decimal places
                if isinstance(af_value, float):
                    af_value = round(af_value, 6)
                elif isinstance(af_value, (list, tuple)) and len(af_value) > 0:
                    af_value = [
                        round(x, 6) if isinstance(x, float) else x for x in af_value
                    ]
                v.INFO["AF"] = af_value
            else:
                # Handle existing AF values that might have precision issues
                af_value = v.INFO["AF"]
                if isinstance(af_value, float):
                    v.INFO["AF"] = round(af_value, 6)
                elif isinstance(af_value, (list, tuple)) and len(af_value) > 0:
                    v.INFO["AF"] = [
                        round(x, 6) if isinstance(x, float) else x for x in af_value
                    ]

            if variant_key in variants:
                if variants[variant_key].INFO["AF"] < v.INFO["AF"]:
                    del variants[variant_key]
                else:
                    continue

            if v.is_indel:
                variant_type = "INDEL"
            elif v.is_snp:
                variant_type = "SNP"
            else:
                variant_type = "OTHER"
            v.INFO["TYPE"] = variant_type

            # Curate variant
            v.ID = v.REF + str(v.POS) + v.ALT[0]
            dp_value = _scalar(getattr(v.INFO, "get", lambda *_: 0)("DP", 0))
            af_value = _scalar(getattr(v.INFO, "get", lambda *_: 0.0)("AF", 0.0))
            dp_array = np.repeat(dp_value, sample_count)
            af_array = np.repeat(af_value, sample_count)
            v.set_format("DP", np.array(dp_array, dtype=int))
            v.set_format("AF", np.array(af_array, dtype=float))
            v.genotypes = [[1, True] for _ in range(sample_count)]
            variants[variant_key] = v

        for variant in variants.values():
            w.write_record(variant)
        w.close()

        # Index the intermediate output file
        cmd = f"bcftools index -f {raw_out}"
        subprocess.run(cmd, shell=True, check=True)

        # Filter variants before the cross-sample merge. Consequence annotation is
        # run only after merging so BCSQ reflects the sample-specific haplotype at
        # each longitudinal timepoint.
        cmd = (
            f'bcftools view -i \'(INFO/AF >= {min_snv_freq} & INFO/TYPE == "SNP") | '
            f'(INFO/AF >= {min_indel_freq} & INFO/TYPE == "INDEL")\' {raw_out} '
            f"-Oz -o {out}"
        )

        _run_logged_command(cmd, log, "bcftools view filter")
        _run_logged_command(f"bcftools index -f {out}", log, "bcftools index")

        if debug:
            print(f"Command: {cmd}")

    except Exception as e:
        raise RuntimeError(f"Error formatting VCF {vcf}: {str(e)}")

    return out, csq_file


def merge_consequences(tempdir, csq_file, sample_names, debug):
    """
    Merge per-sample VCF files before consequence annotation.

    Args:
        tempdir (str): Temporary directory path
        csq_file (str): Output CSQ file path
        sample_names (str): Path to sample names file
        debug (bool): Whether to print debug information
    """
    vcf_list_file = os.path.join(tempdir, "vcf_list.txt")

    # Read the VCF list to determine how many files we have
    with open(vcf_list_file, "r") as f:
        vcf_files = [line.strip() for line in f if line.strip()]

    probe_vcf = VCF(vcf_files[0])
    has_samples = len(probe_vcf.samples) > 0
    probe_vcf.close()

    if len(vcf_files) == 1:
        if has_samples:
            cmd = f"bcftools reheader -s {sample_names} {vcf_files[0]} > {csq_file}"
        else:
            cmd = f"bcftools view -Ov {vcf_files[0]} > {csq_file}"
    else:
        # Handle multiple samples case - merge then reheader
        cmd = (
            f"bcftools merge -0 -Ov -l {vcf_list_file} | "
            f"bcftools reheader -s {sample_names} > {csq_file}"
        )

    if debug:
        print(f"Command: {cmd}")

    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error merging VCF files: {str(e)}")


def _variant_threshold_flag(v) -> str:
    ref = str(v.REF or "")
    alts = [str(alt) for alt in (v.ALT or [])]
    is_indel = any(len(alt) != len(ref) for alt in alts)
    return "--min-indel-freq" if is_indel else "--min-snv-freq"


def _normalise_af_row(sample_values, alt_count: int) -> list[float]:
    values = np.atleast_1d(sample_values).tolist()
    if len(values) < alt_count:
        values.extend([np.nan] * (alt_count - len(values)))
    return values[:alt_count]


def _present_alt_entries(v, sample_values) -> list[_AltEntry]:
    entries: list[_AltEntry] = []
    values = _normalise_af_row(sample_values, len(v.ALT or []))
    for alt_index, raw_value in enumerate(values, start=1):
        try:
            numeric_value = float(raw_value)
        except (TypeError, ValueError):
            continue
        if np.isnan(numeric_value) or numeric_value <= 0:
            continue
        alt = str(v.ALT[alt_index - 1])
        entries.append({"index": alt_index, "alt": alt, "af": numeric_value})
    return entries


def _describe_alt_entries(v, entries: list[_AltEntry]) -> str:
    if not entries:
        return "None"
    ref = str(v.REF or "")
    return ", ".join(f"{ref}>{entry['alt']} AF={entry['af']:.3f}" for entry in entries)


def _build_overflow_events(v, sample_names, af_matrix) -> list[dict[str, object]]:
    events = []
    for sample_name, sample_values in zip(sample_names, af_matrix.tolist()):
        present_entries = _present_alt_entries(v, sample_values)
        if len(present_entries) <= 2:
            continue

        sorted_entries = sorted(
            present_entries, key=lambda entry: (entry["af"], entry["index"])
        )
        drop_count = len(sorted_entries) - 2
        dropped_entries = sorted_entries[:drop_count]
        kept_entries = sorted_entries[drop_count:]
        threshold_suggestion = max(float(entry["af"]) for entry in dropped_entries)

        events.append(
            {
                "sample": sample_name,
                "present": present_entries,
                "dropped": dropped_entries,
                "kept": kept_entries,
                "suggestion": threshold_suggestion,
            }
        )

    return events


def _format_overflow_message(v, overflow_events, *, action: str) -> str:
    threshold_flag = _variant_threshold_flag(v)
    location = f"{v.CHROM}:{v.POS}"

    lines = [
        f"bcftools csq input preparation found more than two ALT alleles present at {location}.",
        f"Action: {action}.",
    ]
    for event in overflow_events:
        lines.append(f"Sample {event['sample']}:")
        lines.append(
            f"Retained ALT alleles after filtering: {_describe_alt_entries(v, event['present'])}."
        )
        lines.append(
            f"To avoid this at the current site, increase {threshold_flag} above {event['suggestion']:.3f}."
        )
        if event["dropped"]:
            lines.append(
                f"Lowest-frequency ALT alleles affected: {_describe_alt_entries(v, event['dropped'])}."
            )
    lines.append("bcftools csq can represent at most two ALT alleles per sample/site.")
    return "\n".join(lines)


def _vcf_record_sort_key(line: str) -> tuple[str, int, str, str]:
    fields = line.rstrip("\n").split("\t")
    chrom = fields[0] if len(fields) > 0 else ""
    try:
        pos = int(fields[1]) if len(fields) > 1 else 0
    except ValueError:
        pos = 0
    ref = fields[3] if len(fields) > 3 else ""
    alt = fields[4] if len(fields) > 4 else ""
    return chrom, pos, ref, alt


def _merge_vcf_outputs(primary_vcf: str, secondary_vcf: str, output_file: str) -> None:
    with open(primary_vcf, "r", encoding="utf-8") as primary_handle:
        primary_lines = primary_handle.readlines()
    with open(secondary_vcf, "r", encoding="utf-8") as secondary_handle:
        secondary_lines = secondary_handle.readlines()

    header_lines = [line for line in primary_lines if line.startswith("#")]
    record_lines = [
        line for line in primary_lines if line.strip() and not line.startswith("#")
    ]
    record_lines.extend(
        line for line in secondary_lines if line.strip() and not line.startswith("#")
    )
    record_lines.sort(key=_vcf_record_sort_key)

    with open(output_file, "w", encoding="utf-8") as out_handle:
        out_handle.writelines(header_lines)
        out_handle.writelines(record_lines)


def annotate_vcf(
    vcf_file,
    output_file,
    reference,
    annotation,
    debug,
    multiallelic_overflow="error",
):
    """Annotate a merged multi-sample VCF with bcftools csq."""
    if multiallelic_overflow not in {"error", "drop-lowest-af", "skip-site"}:
        raise RuntimeError(
            "multiallelic_overflow must be one of: error, drop-lowest-af, skip-site"
        )

    output_dir = os.path.dirname(output_file) or "."
    prefix = Path(output_file).stem
    joined_vcf = os.path.join(output_dir, f"{prefix}.joined_for_csq.vcf")
    prepared_vcf = os.path.join(output_dir, f"{prefix}.prepared_for_csq.vcf")
    skipped_vcf = os.path.join(output_dir, f"{prefix}.skipped_from_csq.vcf")
    annotated_only_vcf = os.path.join(output_dir, f"{prefix}.annotated_only.vcf")

    join_cmd = f"bcftools norm -m +any -Ov -o {joined_vcf} {vcf_file}"
    if debug:
        print(f"Command: {join_cmd}")
    try:
        subprocess.run(join_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Error preparing merged VCF for bcftools csq with multiallelic joins: {str(e)}"
        )

    joined = VCF(joined_vcf, strict_gt=True)
    sample_names = list(joined.samples)
    writer = Writer.from_string(prepared_vcf, joined.raw_header, mode="w")
    skipped_writer = Writer.from_string(skipped_vcf, joined.raw_header, mode="w")
    records_for_csq = 0
    skipped_records = 0

    try:
        for variant in joined:
            try:
                af_matrix = variant.format("AF")
            except (KeyError, AttributeError, ValueError, RuntimeError):
                af_matrix = None

            if af_matrix is None:
                writer.write_record(variant)
                records_for_csq += 1
                continue

            overflow_events = _build_overflow_events(variant, sample_names, af_matrix)
            if overflow_events and multiallelic_overflow == "error":
                raise RuntimeError(
                    _format_overflow_message(
                        variant,
                        overflow_events,
                        action="Stopping analysis before bcftools csq",
                    )
                )

            if overflow_events and multiallelic_overflow == "skip-site":
                print(
                    "Warning: "
                    + _format_overflow_message(
                        variant,
                        overflow_events,
                        action="Skipping consequence calling for this site",
                    )
                )
                skipped_writer.write_record(variant)
                skipped_records += 1
                continue

            updated_af_rows = []
            genotypes = []
            for sample_name, sample_values in zip(sample_names, af_matrix.tolist()):
                row_values = _normalise_af_row(sample_values, len(variant.ALT or []))
                present_entries = _present_alt_entries(variant, row_values)

                if len(present_entries) > 2:
                    sorted_entries = sorted(
                        present_entries,
                        key=lambda entry: (entry["af"], entry["index"]),
                    )
                    dropped_entries = sorted_entries[: len(sorted_entries) - 2]
                    kept_entries = sorted_entries[len(sorted_entries) - 2 :]
                    dropped_indices = {entry["index"] for entry in dropped_entries}

                    for alt_index in dropped_indices:
                        row_values[alt_index - 1] = np.nan

                    print(
                        "Warning: "
                        + _format_overflow_message(
                            variant,
                            [
                                {
                                    "sample": sample_name,
                                    "present": present_entries,
                                    "dropped": dropped_entries,
                                    "kept": kept_entries,
                                    "suggestion": max(
                                        float(entry["af"]) for entry in dropped_entries
                                    ),
                                }
                            ],
                            action="Dropping the lowest-frequency ALT allele(s) for this sample",
                        )
                    )
                    present = sorted(entry["index"] for entry in kept_entries)
                else:
                    present = sorted(entry["index"] for entry in present_entries)

                updated_af_rows.append(row_values)
                if not present:
                    genotypes.append([0, 0, False])
                elif len(present) == 1:
                    genotypes.append([present[0], present[0], False])
                else:
                    genotypes.append([present[0], present[1], False])

            variant.set_format("AF", np.asarray(updated_af_rows, dtype=float))
            variant.genotypes = genotypes
            writer.write_record(variant)
            records_for_csq += 1
    finally:
        writer.close()
        skipped_writer.close()
        joined.close()

    if records_for_csq == 0:
        shutil.copyfile(skipped_vcf, output_file)
        return

    cmd = (
        f"bcftools csq -p R -f {reference} -g {annotation} --force "
        f"-Ov -o {annotated_only_vcf} {prepared_vcf}"
    )

    if debug:
        print(f"Command: {cmd}")

    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error annotating merged VCF with bcftools csq: {str(e)}")

    if skipped_records:
        _merge_vcf_outputs(annotated_only_vcf, skipped_vcf, output_file)
    else:
        shutil.copyfile(annotated_only_vcf, output_file)


def calculate_variant_site_depths(cov_df, v, samples, min_depth: int):
    """
    Calculate depth metrics for variant sites.

    Args:
        cov_df (pd.DataFrame): Coverage data
        v: Variant object from cyvcf2
        samples (list): List of sample names

    Returns:
        dict: Dictionary with depth and QC metrics
    """
    start = v.start + 2  # Account for 0-based coords and wanting first alt base

    if v.is_indel and not v.is_deletion:
        end = v.end + 2
    else:
        end = v.end

    if v.var_type == "snp":
        site_depths = list(cov_df.loc[cov_df["pos"] == start]["depth"])
    else:
        site_depths = (
            cov_df[cov_df["pos"].between(start, end, inclusive="both")]
            .groupby("sample")["depth"]
            .mean()
            .reset_index()
        )
        site_depths["order"] = pd.Categorical(
            site_depths["sample"], categories=samples, ordered=True
        )
        site_depths = [
            int(round(x, 0)) for x in list(site_depths.sort_values("order")["depth"])
        ]

    window_depth = (
        cov_df[cov_df["pos"].between(start - 10, end + 10, inclusive="neither")]
        .groupby("sample")["depth"]
        .mean()
        .reset_index()
    )
    window_depth["order"] = pd.Categorical(
        window_depth["sample"], categories=samples, ordered=True
    )
    window_depth = [
        int(round(x, 0)) for x in list(window_depth.sort_values("order")["depth"])
    ]

    try:
        dp_values = v.format("DP").tolist()
    except (KeyError, AttributeError, ValueError, RuntimeError):
        info_dp = v.INFO.get("DP")
        if isinstance(info_dp, (list, tuple)):
            dp_values = [[int(float(info_dp[0] or 0))]] * len(samples)
        elif info_dp is not None:
            try:
                dp_val = int(float(info_dp))
            except (TypeError, ValueError):
                dp_val = 0
            dp_values = [[dp_val]] * len(samples)
        else:
            dp_values = [[-1]] * len(samples)

    variant_depths = [(x[0] * 0) - 1 if x[0] < 0 else x[0] for x in dp_values]

    variant_qc = []
    for x, y, z in zip(variant_depths, site_depths, window_depth):
        if x >= 0:
            variant_qc.append("P")
        elif x < 0 and y >= min_depth:
            variant_qc.append("P")
        elif v.is_indel and x < 0 and y < min_depth and z >= min_depth:
            variant_qc.append("P")
        elif v.var_type == "snp" and x < 0 and y < min_depth:
            variant_qc.append("F")
        elif x < 0 and y < min_depth and z < min_depth:
            variant_qc.append("F")

    variant_depths = ["M" if x == -1 else x for x in variant_depths]

    result = {
        "site_depth": [str(x) for x in site_depths],
        "window_depth": [str(x) for x in window_depth],
        "variant_depth": [str(x) for x in variant_depths],
        "variant_qc": [str(x) for x in variant_qc],
    }

    return result


def _summarise_sample_trajectory(allele_freqs, samples):
    """Derive trajectory metadata from per-sample allele frequencies."""
    presence_absence = ["N" if x == "." else "Y" for x in allele_freqs]
    variant_status = "new" if allele_freqs[0] == "." else "original"

    if allele_freqs[0] != "." and allele_freqs[-1] == ".":
        persistent_status = "original_lost"
    elif allele_freqs[0] != "." and allele_freqs[-1] != ".":
        persistent_status = "original_retained"
    elif allele_freqs[0] == "." and allele_freqs[-1] != ".":
        persistent_status = "new_persistent"
    elif allele_freqs[0] == "." and allele_freqs[-1] == ".":
        persistent_status = "new_transient"
    else:
        persistent_status = "unknown"

    first_appearance = (
        samples[presence_absence.index("Y")] if "Y" in presence_absence else "None"
    )
    last_appearance = (
        samples[rindex(presence_absence, "Y")] if "Y" in presence_absence else "None"
    )

    return {
        "presence_absence": presence_absence,
        "variant_status": variant_status,
        "persistent_status": persistent_status,
        "first_appearance": first_appearance,
        "last_appearance": last_appearance,
    }


def _extract_scalar_mask(value) -> int:
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return 0
        value = value.flat[0]
    elif isinstance(value, (list, tuple)):
        if not value:
            return 0
        value = value[0]
    try:
        return int(value or 0)
    except (TypeError, ValueError):
        return 0


def _decode_sample_bcsq_annotations(v, samples, annotations):
    """Decode FORMAT/BCSQ bitmasks into sample-specific annotation assignments."""
    if not annotations:
        return {}

    try:
        mask_values = v.format("BCSQ")
    except (KeyError, AttributeError, ValueError, RuntimeError):
        return {}

    if mask_values is None:
        return {}

    decoded = {}
    for sample, sample_mask in zip(samples, mask_values):
        bitmask = _extract_scalar_mask(sample_mask)
        matches = []
        for idx, annotation in enumerate(annotations):
            first_haplotype_bit = 1 << (2 * idx)
            second_haplotype_bit = 1 << (2 * idx + 1)
            if bitmask & first_haplotype_bit or bitmask & second_haplotype_bit:
                matches.append(annotation)
        decoded[sample] = matches

    return decoded


def _mask_allele_frequencies_for_annotation(
    annotation, allele_freqs, samples, sample_bcsq_map
):
    """Keep allele frequencies only for samples where the annotation applies."""
    masked = []
    for sample, allele_freq in zip(samples, allele_freqs):
        if allele_freq == ".":
            masked.append(".")
            continue
        masked.append(
            allele_freq if annotation in sample_bcsq_map.get(sample, []) else "."
        )
    return masked


def _format_frequency_token(value) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "."
    if np.isnan(numeric):
        return "."
    return "{:.3f}".format(numeric)


def _extract_alt_frequency_map(v, samples):
    """Return per-ALT allele-frequency trajectories for a VCF record."""

    alt_map = {str(alt): ["."] * len(samples) for alt in (v.ALT or [])}
    if not alt_map:
        return alt_map

    try:
        allele_freq_arrays = v.format("AF")
    except (KeyError, AttributeError, ValueError, RuntimeError):
        allele_freq_arrays = None

    if allele_freq_arrays is not None:
        for sample_index, sample_values in enumerate(allele_freq_arrays.tolist()):
            values = np.atleast_1d(sample_values).tolist()
            for alt_index, alt in enumerate(v.ALT):
                raw_value = values[alt_index] if alt_index < len(values) else None
                alt_map[str(alt)][sample_index] = _format_frequency_token(raw_value)
        return alt_map

    info_af = v.INFO.get("AF") or v.INFO.get("VAF")
    if isinstance(info_af, (list, tuple)):
        values = list(info_af)
    elif info_af is None:
        values = []
    else:
        values = [info_af]

    for alt_index, alt in enumerate(v.ALT):
        raw_value = values[alt_index] if alt_index < len(values) else None
        token = _format_frequency_token(raw_value)
        alt_map[str(alt)] = [token] * len(samples)

    return alt_map


def _collapse_alt_frequency_map(alt_frequency_map, sample_count: int):
    """Collapse per-ALT frequencies into a site-level trajectory using max AF."""

    collapsed = []
    for sample_index in range(sample_count):
        observed = []
        for allele_freqs in alt_frequency_map.values():
            if sample_index >= len(allele_freqs):
                continue
            token = allele_freqs[sample_index]
            if token == ".":
                continue
            try:
                observed.append(float(token))
            except (TypeError, ValueError):
                continue
        collapsed.append("{:.3f}".format(max(observed)) if observed else ".")
    return collapsed


def _annotation_alt_for_record(v, anno):
    """Infer which ALT allele a BCSQ annotation belongs to for the current record."""

    if not v.ALT:
        return None
    if len(v.ALT) == 1:
        return str(v.ALT[0])
    if len(anno) <= 6:
        return None

    dna_change = str(anno[6] or "").strip()
    if not dna_change:
        return None

    for token in dna_change.split("+"):
        match = re.match(r"^(\d+)([^>]+)>([^>]+)$", token.strip())
        if not match:
            continue
        pos = int(match.group(1))
        alt = match.group(3)
        if pos == int(v.POS) and alt in v.ALT:
            return alt

    return None


def process_vcf(vcf_file, covs, min_depth, sample_names_override=None):
    """
    Process VCF file and extract variant information.

    Args:
        vcf_file (str): Path to VCF file
        covs (list): List of coverage file paths

    Returns:
        pd.DataFrame: DataFrame with variant information
    """
    results = []
    vcf = VCF(vcf_file)
    samples = list(vcf.samples)

    if sample_names_override is not None:
        override = list(sample_names_override)
        if not samples or len(samples) != len(override):
            samples = override

    if not samples:
        raise RuntimeError(
            "No sample columns present in VCF and no sample names could be inferred."
        )

    if len(samples) != len(covs):
        raise RuntimeError(
            "Number of coverage files does not match number of VCF samples."
        )

    # Load coverage data
    cov_list = []
    total_cov_list = []

    for i, cov_file in enumerate(covs):
        try:
            cov = pd.read_csv(cov_file, sep="\t", names=["ref", "pos", "depth"])
            total_cov = "{:.2f}".format(((sum(cov.depth > 10) / cov.shape[0]) * 100))
            total_cov_list.append(total_cov)
            cov["sample"] = samples[i]
            cov_list.append(cov)
        except Exception as e:
            raise RuntimeError(f"Error reading coverage file {cov_file}: {str(e)}")

    cov_df = pd.concat(cov_list).set_index("sample")

    # Process variants
    for v in vcf:
        info = dict(list(v.INFO))
        alt_frequency_map = _extract_alt_frequency_map(v, samples)
        allele_freqs = _collapse_alt_frequency_map(alt_frequency_map, len(samples))
        if not allele_freqs:
            allele_freqs = ["."] * max(1, len(samples))
        trajectory = _summarise_sample_trajectory(allele_freqs, samples)

        depths_qc = calculate_variant_site_depths(cov_df, v, samples, min_depth)
        all_samples_pass_qc = "F" not in depths_qc["variant_qc"]
        proportion_samples_passing_qc = (
            sum(flag == "P" for flag in depths_qc["variant_qc"])
            / len(depths_qc["variant_qc"])
            if depths_qc["variant_qc"]
            else 0.0
        )

        # Process annotations
        if "BCSQ" in info:
            annotations = v.INFO["BCSQ"].split(",")
            sample_bcsq_map = _decode_sample_bcsq_annotations(v, samples, annotations)
            produced_annotation_specific_row = False

            if sample_bcsq_map:
                for annot in annotations:
                    anno = annot.split("|")
                    annotation_alt = _annotation_alt_for_record(v, anno)
                    annotation_allele_freqs = (
                        alt_frequency_map.get(annotation_alt, allele_freqs)
                        if annotation_alt is not None
                        else allele_freqs
                    )
                    masked_allele_freqs = _mask_allele_frequencies_for_annotation(
                        annot, annotation_allele_freqs, samples, sample_bcsq_map
                    )
                    masked_trajectory = _summarise_sample_trajectory(
                        masked_allele_freqs, samples
                    )
                    if "Y" not in masked_trajectory["presence_absence"]:
                        continue

                    result = _process_annotation(
                        v,
                        anno,
                        masked_trajectory["variant_status"],
                        masked_trajectory["persistent_status"],
                        masked_trajectory["presence_absence"],
                        masked_trajectory["first_appearance"],
                        masked_trajectory["last_appearance"],
                        all_samples_pass_qc,
                        proportion_samples_passing_qc,
                        depths_qc,
                        masked_allele_freqs,
                        samples,
                        total_cov_list,
                        annotation_alt,
                    )
                    results.append(result)
                    produced_annotation_specific_row = True

            if not produced_annotation_specific_row:
                for annot in annotations:
                    anno = annot.split("|")
                    annotation_alt = _annotation_alt_for_record(v, anno)
                    annotation_allele_freqs = (
                        alt_frequency_map.get(annotation_alt, allele_freqs)
                        if annotation_alt is not None
                        else allele_freqs
                    )
                    annotation_trajectory = _summarise_sample_trajectory(
                        annotation_allele_freqs, samples
                    )
                    result = _process_annotation(
                        v,
                        anno,
                        annotation_trajectory["variant_status"],
                        annotation_trajectory["persistent_status"],
                        annotation_trajectory["presence_absence"],
                        annotation_trajectory["first_appearance"],
                        annotation_trajectory["last_appearance"],
                        all_samples_pass_qc,
                        proportion_samples_passing_qc,
                        depths_qc,
                        annotation_allele_freqs,
                        samples,
                        total_cov_list,
                        annotation_alt,
                    )
                    results.append(result)
        else:
            # Handle variants without annotations
            if len(alt_frequency_map) > 1:
                for alt, alt_allele_freqs in alt_frequency_map.items():
                    alt_trajectory = _summarise_sample_trajectory(
                        alt_allele_freqs, samples
                    )
                    if "Y" not in alt_trajectory["presence_absence"]:
                        continue
                    result = _create_unannotated_result(
                        v,
                        alt_trajectory["variant_status"],
                        alt_trajectory["persistent_status"],
                        alt_trajectory["presence_absence"],
                        alt_trajectory["first_appearance"],
                        alt_trajectory["last_appearance"],
                        all_samples_pass_qc,
                        proportion_samples_passing_qc,
                        depths_qc,
                        alt_allele_freqs,
                        samples,
                        total_cov_list,
                        alt,
                    )
                    results.append(result)
            else:
                result = _create_unannotated_result(
                    v,
                    trajectory["variant_status"],
                    trajectory["persistent_status"],
                    trajectory["presence_absence"],
                    trajectory["first_appearance"],
                    trajectory["last_appearance"],
                    all_samples_pass_qc,
                    proportion_samples_passing_qc,
                    depths_qc,
                    allele_freqs,
                    samples,
                    total_cov_list,
                )
                results.append(result)

    return pd.DataFrame(results)


def _process_annotation(
    v,
    anno,
    variant_status,
    persistent_status,
    presence_absence,
    first_appearance,
    last_appearance,
    all_samples_pass_qc,
    proportion_samples_passing_qc,
    depths_qc,
    allele_freqs,
    samples,
    total_cov_list,
    alt_allele=None,
):
    """Process a single annotation from bcftools csq."""
    selected_alt = alt_allele or (v.ALT[0] if v.ALT else "")
    if len(anno) == 1 and anno[0].startswith("@"):
        # Joint variant annotation
        return {
            "chrom": v.CHROM,
            "start": v.start + 1,
            "end": v.end,
            "gene": anno[0],
            "ref": v.REF,
            "alt": selected_alt,
            "variant": v.REF + str(v.POS) + selected_alt,
            "amino_acid_consequence": anno[0],
            "nsp_aa_change": anno[0],
            "bcsq_nt_notation": anno[0],
            "bcsq_aa_notation": anno[0],
            "type_of_variant": v.var_type,
            "type_of_change": anno[0],
            "variant_status": variant_status,
            "persistence_status": persistent_status,
            "presence_absence": " / ".join(presence_absence),
            "first_appearance": first_appearance,
            "last_appearance": last_appearance,
            "all_samples_pass_qc": all_samples_pass_qc,
            "proportion_samples_passing_qc": proportion_samples_passing_qc,
            "per_sample_variant_qc": " / ".join(depths_qc["variant_qc"]),
            "aa1_total_properties": anno[0],
            "aa2_total_properties": anno[0],
            "aa1_unique_properties": anno[0],
            "aa2_unique_properties": anno[0],
            "aa1_weight": anno[0],
            "aa2_weight": anno[0],
            "weight_difference": anno[0],
            "alt_freq": " / ".join(allele_freqs),
            "variant_depth": " / ".join(depths_qc["variant_depth"]),
            "variant_site_depth": " / ".join(depths_qc["site_depth"]),
            "variant_window_depth": " / ".join(depths_qc["window_depth"]),
            "samples": " / ".join(samples),
            "total_genome_coverage": " / ".join(total_cov_list),
        }
    else:
        # Regular annotation
        reformatted_aa = (
            reformat_csq_notation(anno[1], anno[5])
            if len(anno) > 5
            else [anno[1] + ":" + anno[0], ""]
        )

        aa_exploration = AminoAcidChange(reformatted_aa[0].replace("*", ""))

        return {
            "chrom": v.CHROM,
            "start": v.start + 1,
            "end": v.end,
            "gene": anno[1],
            "ref": v.REF,
            "alt": selected_alt,
            "variant": v.REF + str(v.POS) + selected_alt,
            "amino_acid_consequence": reformatted_aa[0],
            "nsp_aa_change": reformatted_aa[1],
            "bcsq_nt_notation": anno[6] if len(anno) > 5 else "",
            "bcsq_aa_notation": anno[5] if len(anno) > 5 else "",
            "type_of_variant": v.var_type,
            "type_of_change": anno[0],
            "variant_status": variant_status,
            "persistence_status": persistent_status,
            "presence_absence": " / ".join(presence_absence),
            "first_appearance": first_appearance,
            "last_appearance": last_appearance,
            "all_samples_pass_qc": all_samples_pass_qc,
            "proportion_samples_passing_qc": proportion_samples_passing_qc,
            "per_sample_variant_qc": " / ".join(depths_qc["variant_qc"]),
            "aa1_total_properties": ";".join(aa_exploration.aa1_total_properties),
            "aa2_total_properties": ";".join(aa_exploration.aa2_total_properties),
            "aa1_unique_properties": ";".join(aa_exploration.aa1_unique_properties),
            "aa2_unique_properties": ";".join(aa_exploration.aa2_unique_properties),
            "aa1_weight": aa_exploration.aa1_weight,
            "aa2_weight": aa_exploration.aa2_weight,
            "weight_difference": aa_exploration.weight_difference,
            "alt_freq": " / ".join(allele_freqs),
            "variant_depth": " / ".join(depths_qc["variant_depth"]),
            "variant_site_depth": " / ".join(depths_qc["site_depth"]),
            "variant_window_depth": " / ".join(depths_qc["window_depth"]),
            "samples": " / ".join(samples),
            "total_genome_coverage": " / ".join(total_cov_list),
        }


def _create_unannotated_result(
    v,
    variant_status,
    persistent_status,
    presence_absence,
    first_appearance,
    last_appearance,
    all_samples_pass_qc,
    proportion_samples_passing_qc,
    depths_qc,
    allele_freqs,
    samples,
    total_cov_list,
    alt_allele=None,
):
    """Create result for variants without annotations."""
    selected_alt = alt_allele or (v.ALT[0] if v.ALT else "")
    if v.start + 1 < 266:
        gene = "5' UTR"
    elif v.start + 1 > 29674:
        gene = "3' UTR"
    else:
        gene = "INTERGENIC"

    return {
        "chrom": v.CHROM,
        "start": v.start + 1,
        "end": v.end,
        "gene": gene,
        "ref": v.REF,
        "alt": selected_alt,
        "variant": v.REF + str(v.POS) + selected_alt,
        "amino_acid_consequence": "None",
        "nsp_aa_change": "None",
        "bcsq_nt_notation": "None",
        "bcsq_aa_notation": "None",
        "type_of_variant": v.var_type,
        "type_of_change": "None",
        "variant_status": variant_status,
        "persistence_status": persistent_status,
        "presence_absence": " / ".join(presence_absence),
        "first_appearance": first_appearance,
        "last_appearance": last_appearance,
        "all_samples_pass_qc": all_samples_pass_qc,
        "proportion_samples_passing_qc": proportion_samples_passing_qc,
        "per_sample_variant_qc": " / ".join(depths_qc["variant_qc"]),
        "aa1_total_properties": "None",
        "aa2_total_properties": "None",
        "aa1_unique_properties": "None",
        "aa2_unique_properties": "None",
        "aa1_weight": "None",
        "aa2_weight": "None",
        "weight_difference": "None",
        "alt_freq": " / ".join(allele_freqs),
        "variant_depth": " / ".join(depths_qc["variant_depth"]),
        "variant_site_depth": " / ".join(depths_qc["site_depth"]),
        "variant_window_depth": " / ".join(depths_qc["window_depth"]),
        "samples": " / ".join(samples),
        "total_genome_coverage": " / ".join(total_cov_list),
    }
