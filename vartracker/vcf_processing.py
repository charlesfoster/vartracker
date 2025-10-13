"""
VCF file processing functionality for vartracker.

Contains functions for formatting, merging, and processing VCF files.
"""

import os
import re
import subprocess
import gzip
from pathlib import Path

import pandas as pd
import numpy as np
from cyvcf2 import VCF, Writer
from itertools import dropwhile

from .constants import reformat_csq_notation
from .amino_acids import AminoAcidChange


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
    Format VCF file for compatibility and add consequences.

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

    try:
        prepared_vcf = _ensure_format_and_sample(vcf, sample, tempdir, debug)

        vcf_mod = VCF(prepared_vcf, strict_gt=True)
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
                    "Number": "1",
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
                "Number": "1",
            }
        )

        if existing_samples:
            w = Writer(out, vcf_mod, mode="wz")
        else:
            header_lines = []
            for line in vcf_mod.raw_header.strip().splitlines():
                if line.startswith("#CHROM"):
                    header_lines.append(f"{line}\t{sample}")
                else:
                    header_lines.append(line)
            header_str = "\n".join(header_lines) + "\n"
            w = Writer.from_string(out, header_str, mode="wz")
        variants = {}

        def _scalar(value):
            if isinstance(value, (list, tuple, np.ndarray)):
                return value[0] if value else 0
            return value

        for v in vcf_mod:
            pos_key = str(v.POS)

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

            if pos_key in variants:
                if variants[pos_key].INFO["AF"] < v.INFO["AF"]:
                    del variants[pos_key]
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
            variants[pos_key] = v

        for variant in variants.values():
            w.write_record(variant)
        w.close()

        # Index the output file
        cmd = f"bcftools index -f {out}"
        subprocess.run(cmd, shell=True, check=True)

        # Run csq before merging for better AF filtering
        cmd = (
            f'bcftools view -i \'(INFO/AF >= {min_snv_freq} & INFO/TYPE == "SNP") | '
            f'(INFO/AF >= {min_indel_freq} & INFO/TYPE == "INDEL")\' {out} | '
            f"bcftools csq -f {reference} -g {annotation} --force -Oz -o {csq_file}; "
            f"tabix -f -p vcf {csq_file}"
        )

        try:
            with open(log, "a", encoding="utf-8") as err:
                subprocess.run(cmd, shell=True, stderr=err, check=True)
        except subprocess.CalledProcessError as exc:
            log_tail = ""
            if os.path.exists(log):
                try:
                    with open(log, "r", encoding="utf-8", errors="replace") as err_file:
                        tail_lines = err_file.readlines()[-20:]
                        log_tail = "".join(tail_lines).strip()
                except OSError:
                    log_tail = "<unable to read bcftools log>"

            details = log_tail if log_tail else exc.stderr or ""
            message = f"Command '{cmd}' returned non-zero exit status {exc.returncode}."
            if details:
                message = f"{message}\nLast bcftools log lines:\n{details}"
            raise RuntimeError(message) from exc

        if debug:
            print(f"Command: {cmd}")

    except Exception as e:
        raise RuntimeError(f"Error formatting VCF {vcf}: {str(e)}")

    return out, csq_file


def merge_consequences(tempdir, csq_file, sample_names, debug):
    """
    Merge VCF files with consequences.

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
        try:
            allele_freq_arrays = v.format("AF").tolist()
        except (KeyError, AttributeError, ValueError, RuntimeError):
            allele_freq_arrays = None

        if allele_freq_arrays:
            allele_freqs = [
                "{:.3f}".format(x[0]).replace("nan", ".") for x in allele_freq_arrays
            ]
        else:
            info_af = v.INFO.get("AF") or v.INFO.get("VAF")
            if isinstance(info_af, (list, tuple)):
                values = list(info_af)
            elif info_af is None:
                values = [None] * len(samples)
            else:
                values = [info_af]
            if not values:
                values = [None]
            if len(values) < len(samples):
                values.extend([values[-1]] * (len(samples) - len(values)))
            allele_freqs = []
            for val in values[: len(samples)]:
                if val is None:
                    allele_freqs.append(".")
                    continue
                try:
                    allele_freqs.append("{:.3f}".format(float(val)))
                except (TypeError, ValueError):
                    allele_freqs.append(".")
            if not allele_freqs:
                allele_freqs = ["."] * max(1, len(samples))
        presence_absence = ["N" if x == "." else "Y" for x in allele_freqs]
        variant_status = "new" if allele_freqs[0] == "." else "original"

        # Calculate persistence status
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

        depths_qc = calculate_variant_site_depths(cov_df, v, samples, min_depth)
        overall_variant_qc = "FAIL" if "F" in depths_qc["variant_qc"] else "PASS"
        first_appearance = (
            samples[presence_absence.index("Y")] if "Y" in presence_absence else "None"
        )
        last_appearance = (
            samples[rindex(presence_absence, "Y")]
            if "Y" in presence_absence
            else "None"
        )

        # Process annotations
        if "BCSQ" in info:
            annotations = v.INFO["BCSQ"].split(",")
            for annot in annotations:
                anno = annot.split("|")
                result = _process_annotation(
                    v,
                    anno,
                    variant_status,
                    persistent_status,
                    presence_absence,
                    first_appearance,
                    last_appearance,
                    overall_variant_qc,
                    depths_qc,
                    allele_freqs,
                    samples,
                    total_cov_list,
                )
                results.append(result)
        else:
            # Handle variants without annotations
            result = _create_unannotated_result(
                v,
                variant_status,
                persistent_status,
                presence_absence,
                first_appearance,
                last_appearance,
                overall_variant_qc,
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
    overall_variant_qc,
    depths_qc,
    allele_freqs,
    samples,
    total_cov_list,
):
    """Process a single annotation from bcftools csq."""
    if len(anno) == 1 and anno[0].startswith("@"):
        # Joint variant annotation
        return {
            "chrom": v.CHROM,
            "start": v.start + 1,
            "end": v.end,
            "gene": anno[0],
            "ref": v.REF,
            "alt": v.ALT[0],
            "variant": v.REF + str(v.POS) + v.ALT[0],
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
            "overall_variant_qc": overall_variant_qc,
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
            "alt": v.ALT[0],
            "variant": v.REF + str(v.POS) + v.ALT[0],
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
            "overall_variant_qc": overall_variant_qc,
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
    overall_variant_qc,
    depths_qc,
    allele_freqs,
    samples,
    total_cov_list,
):
    """Create result for variants without annotations."""
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
        "alt": v.ALT[0],
        "variant": v.REF + str(v.POS) + v.ALT[0],
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
        "overall_variant_qc": overall_variant_qc,
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
