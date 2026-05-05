"""Rescue high-confidence LoFreq SNPs filtered at primer-overlap sites."""

from __future__ import annotations

import argparse
import gzip
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO, cast


RESCUE_FILTER = (
    '##FILTER=<ID=RESCUED_PRIMER_OVERLAP,Description="Variant failed default LoFreq '
    'filtering but was retained by primer-overlap near-fixed rescue rule">'
)
PRIMER_INFO = (
    '##INFO=<ID=PRIMER_OVERLAP,Number=0,Type=Flag,Description="Variant overlaps '
    'primer interval">'
)
RESCUED_BY_INFO = (
    '##INFO=<ID=RESCUED_BY,Number=1,Type=String,Description="Reason variant was '
    'rescued">'
)


@dataclass(frozen=True)
class PrimerRescueThresholds:
    """Thresholds applied only to primer-overlap rescue candidates."""

    min_af: float = 0.95
    min_dp: int = 100
    min_alt_count: int = 95
    min_qual: float = 100.0
    max_ref_count: int = 20


@dataclass(frozen=True)
class PrimerRescueResult:
    """Summary of a LoFreq primer rescue run."""

    normal_passed: int
    rescued: int
    discarded: int
    rescued_tsv: str


def _open_text(path: str | Path, mode: str = "rt") -> TextIO:
    text_path = str(path)
    if text_path.endswith(".gz"):
        return cast(TextIO, gzip.open(text_path, mode, encoding="utf-8"))
    return cast(TextIO, open(text_path, mode, encoding="utf-8"))


def load_primers(path: str | Path) -> dict[str, list[tuple[int, int]]]:
    """Load a primer BED file as 0-based half-open intervals."""

    intervals: dict[str, list[tuple[int, int]]] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("track ") or line.startswith("browser "):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            intervals.setdefault(chrom, []).append((start, end))
    return intervals


def overlaps_primer(
    chrom: str, start: int, end: int, intervals: dict[str, list[tuple[int, int]]]
) -> bool:
    """Return whether a half-open interval overlaps a primer interval."""

    for primer_start, primer_end in intervals.get(chrom, []):
        if primer_start < end and start < primer_end:
            return True
    return False


def parse_info(info_text: str) -> dict[str, str | bool]:
    """Parse a VCF INFO string into a simple mapping."""

    info: dict[str, str | bool] = {}
    if info_text in ("", "."):
        return info
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def _first_float(value: str | bool) -> float:
    return float(str(value).split(",", 1)[0])


def _first_int(value: str | bool) -> int:
    return int(str(value).split(",", 1)[0])


def _dp4_counts(value: str | bool) -> list[int] | None:
    parts = str(value).split(",")
    if len(parts) != 4:
        return None
    return [int(part) for part in parts]


def _record_key(fields: list[str]) -> tuple[str, str, str, str]:
    return (fields[0], fields[1], fields[3], fields[4])


def _filter_is_pass(fields: list[str]) -> bool:
    return fields[6] in ("PASS", ".")


def read_default_filter_results(
    path: str | Path,
) -> tuple[set[tuple[str, str, str, str]], dict[tuple[str, str, str, str], str]]:
    """Read records emitted by ``lofreq filter`` and return PASS keys/reasons."""

    pass_keys: set[tuple[str, str, str, str]] = set()
    filter_reasons: dict[tuple[str, str, str, str], str] = {}
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 8:
                key = _record_key(fields)
                filter_reasons[key] = fields[6]
                if _filter_is_pass(fields):
                    pass_keys.add(key)
    return pass_keys, filter_reasons


def rescue_metrics(
    fields: list[str],
    primers: dict[str, list[tuple[int, int]]],
    thresholds: PrimerRescueThresholds,
) -> str | None:
    """Return rescue metrics for a candidate record, or ``None`` if it fails."""

    chrom = fields[0]
    pos = int(fields[1])
    ref = fields[3]
    alt = fields[4]
    qual_text = fields[5]

    if "," in alt:
        return None
    if len(ref) != 1 or len(alt) != 1:
        return None

    start = pos - 1
    end = pos
    if not overlaps_primer(chrom, start, end, primers):
        return None

    if qual_text == ".":
        return None

    info = parse_info(fields[7])
    if "AF" not in info or "DP" not in info or "DP4" not in info:
        return None

    try:
        af = _first_float(info["AF"])
        dp = _first_int(info["DP"])
        qual = float(qual_text)
        counts = _dp4_counts(info["DP4"])
    except ValueError:
        return None

    if counts is None:
        return None

    ref_fwd, ref_rev, alt_fwd, alt_rev = counts
    ref_count = ref_fwd + ref_rev
    alt_count = alt_fwd + alt_rev
    dp4_total = ref_count + alt_count

    if dp4_total == 0:
        return None
    if af < thresholds.min_af:
        return None
    if dp < thresholds.min_dp:
        return None
    if qual < thresholds.min_qual:
        return None
    if alt_count < thresholds.min_alt_count:
        return None
    if alt_count / dp4_total < thresholds.min_af:
        return None
    if ref_count > thresholds.max_ref_count:
        return None
    if not (alt_fwd == 0 or alt_rev == 0):
        return None

    dp4_af = alt_count / dp4_total
    return (
        f"AF={af:g},DP={dp},QUAL={qual:g},"
        f"DP4={ref_fwd}/{ref_rev}/{alt_fwd}/{alt_rev},"
        f"alt_count={alt_count},ref_count={ref_count},dp4_af={dp4_af:g}"
    )


def _add_info_flag(info_text: str, flag: str) -> str:
    if info_text in ("", "."):
        return flag
    items = info_text.split(";")
    if flag not in items:
        items.append(flag)
    return ";".join(items)


def _add_info_value(info_text: str, key: str, value: str) -> str:
    item = f"{key}={value}"
    if info_text in ("", "."):
        return item

    items = []
    replaced = False
    for existing in info_text.split(";"):
        if existing == key or existing.startswith(key + "="):
            items.append(item)
            replaced = True
        else:
            items.append(existing)
    if not replaced:
        items.append(item)
    return ";".join(items)


def _rescued_fields(fields: list[str]) -> list[str]:
    out = fields[:]
    out[6] = "RESCUED_PRIMER_OVERLAP"
    out[7] = _add_info_flag(out[7], "PRIMER_OVERLAP")
    out[7] = _add_info_value(out[7], "RESCUED_BY", "overlap_primer_interval")
    return out


def _pass_fields(fields: list[str]) -> list[str]:
    out = fields[:]
    out[6] = "PASS"
    return out


def _variant_name(fields: list[str]) -> str:
    return f"{fields[3]}{fields[1]}{fields[4]}"


def _write_headers(headers: list[str], output: TextIO) -> None:
    has_rescue_filter = any(
        line.startswith("##FILTER=<ID=RESCUED_PRIMER_OVERLAP,") for line in headers
    )
    has_primer_info = any(
        line.startswith("##INFO=<ID=PRIMER_OVERLAP,") for line in headers
    )
    has_rescued_by_info = any(
        line.startswith("##INFO=<ID=RESCUED_BY,") for line in headers
    )

    for line in headers:
        if line.startswith("#CHROM"):
            if not has_rescue_filter:
                output.write(RESCUE_FILTER + "\n")
            if not has_primer_info:
                output.write(PRIMER_INFO + "\n")
            if not has_rescued_by_info:
                output.write(RESCUED_BY_INFO + "\n")
        output.write(line)


def write_final_vcf(
    raw_vcf: str | Path,
    output_vcf: str | Path,
    rescued_tsv_path: str | Path,
    pass_keys: set[tuple[str, str, str, str]],
    filter_reasons: dict[tuple[str, str, str, str], str],
    primers: dict[str, list[tuple[int, int]]],
    thresholds: PrimerRescueThresholds,
) -> PrimerRescueResult:
    """Write the default-filtered VCF plus any rescued primer-overlap records."""

    normal_passed = 0
    rescued = 0
    discarded = 0
    headers: list[str] = []

    with (
        _open_text(raw_vcf) as raw,
        _open_text(output_vcf, "wt") as output,
        Path(rescued_tsv_path).open("w", encoding="utf-8") as rescued_tsv,
    ):
        rescued_tsv.write("variant\treason_filtered\treason_rescued\tmetrics\n")

        for line in raw:
            if line.startswith("#"):
                headers.append(line)
                continue

            if headers:
                _write_headers(headers, output)
                headers = []

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                discarded += 1
                continue

            key = _record_key(fields)
            if key in pass_keys:
                output.write("\t".join(_pass_fields(fields)) + "\n")
                normal_passed += 1
            else:
                metrics = rescue_metrics(fields, primers, thresholds)
                if metrics is None:
                    discarded += 1
                    continue

                output.write("\t".join(_rescued_fields(fields)) + "\n")
                reason_filtered = filter_reasons.get(
                    key, "not_passed_default_lofreq_filter"
                )
                rescued_tsv.write(
                    f"{_variant_name(fields)}\t{reason_filtered}\t"
                    f"overlap_primer_interval\t{metrics}\n"
                )
                rescued += 1

        if headers:
            _write_headers(headers, output)

    return PrimerRescueResult(
        normal_passed=normal_passed,
        rescued=rescued,
        discarded=discarded,
        rescued_tsv=str(rescued_tsv_path),
    )


def rescue_lofreq_primer_variants(
    raw_vcf: str | Path,
    primers_bed: str | Path,
    output_vcf: str | Path,
    rescued_tsv: str | Path | None = None,
    tmp_dir: str | Path | None = None,
    lofreq: str = "lofreq",
    thresholds: PrimerRescueThresholds | None = None,
) -> PrimerRescueResult:
    """Apply default LoFreq filtering plus the primer-overlap rescue rule."""

    thresholds = thresholds or PrimerRescueThresholds()
    if rescued_tsv is None:
        rescued_tsv = f"{output_vcf}.rescued.tsv"

    primers = load_primers(primers_bed)

    with tempfile.TemporaryDirectory(dir=tmp_dir) as run_tmp_dir:
        filtered_vcf = os.path.join(run_tmp_dir, "lofreq.default_filtered.vcf")
        subprocess.run(
            [lofreq, "filter", "-i", str(raw_vcf), "-o", filtered_vcf],
            check=True,
        )
        pass_keys, filter_reasons = read_default_filter_results(filtered_vcf)

    return write_final_vcf(
        raw_vcf,
        output_vcf,
        rescued_tsv,
        pass_keys,
        filter_reasons,
        primers,
        thresholds,
    )


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Preserve default LoFreq filtering while rescuing near-fixed "
            "primer-overlap SNPs as RESCUED_PRIMER_OVERLAP."
        )
    )
    parser.add_argument(
        "--raw-vcf", required=True, help="Raw LoFreq VCF made with --no-default-filter"
    )
    parser.add_argument("--primers-bed", required=True, help="Primer BED file")
    parser.add_argument("--output-vcf", required=True, help="Final filtered VCF")
    parser.add_argument(
        "--rescued-tsv",
        help="TSV of rescued variants; default: OUTPUT_VCF.rescued.tsv",
    )
    parser.add_argument("--tmp-dir", help="Optional temporary directory")
    parser.add_argument("--lofreq", default="lofreq", help="lofreq executable")
    parser.add_argument(
        "--min-af",
        type=float,
        default=PrimerRescueThresholds.min_af,
        help="Minimum INFO/AF for rescue candidates only (default: 0.95)",
    )
    parser.add_argument(
        "--min-dp",
        type=int,
        default=PrimerRescueThresholds.min_dp,
        help="Minimum INFO/DP for rescue candidates only (default: 100)",
    )
    parser.add_argument(
        "--min-alt-count",
        type=int,
        default=PrimerRescueThresholds.min_alt_count,
        help="Minimum DP4 alt count for rescue candidates only (default: 95)",
    )
    parser.add_argument(
        "--min-qual",
        type=float,
        default=PrimerRescueThresholds.min_qual,
        help="Minimum QUAL for rescue candidates only (default: 100)",
    )
    parser.add_argument(
        "--max-ref-count",
        type=int,
        default=PrimerRescueThresholds.max_ref_count,
        help="Maximum DP4 ref count for rescue candidates only (default: 20)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    thresholds = PrimerRescueThresholds(
        min_af=args.min_af,
        min_dp=args.min_dp,
        min_alt_count=args.min_alt_count,
        min_qual=args.min_qual,
        max_ref_count=args.max_ref_count,
    )
    result = rescue_lofreq_primer_variants(
        raw_vcf=args.raw_vcf,
        primers_bed=args.primers_bed,
        output_vcf=args.output_vcf,
        rescued_tsv=args.rescued_tsv,
        tmp_dir=args.tmp_dir,
        lofreq=args.lofreq,
        thresholds=thresholds,
    )
    print(
        f"normal_passed={result.normal_passed} rescued={result.rescued} "
        f"discarded={result.discarded}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
