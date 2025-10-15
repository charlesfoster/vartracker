"""Utilities for generating vartracker input spreadsheets from directories."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, List, Optional


READ_PATTERN = re.compile(
    r"^(?P<stem>.*?)(?:[_-]R(?P<read>[12]))?\.(?:f(?:ast)?q)(?:\.gz)?$",
    re.IGNORECASE,
)
BAM_PATTERN = re.compile(r"^(?P<stem>.*?)(?:\.bam)$", re.IGNORECASE)
VCF_PATTERN = re.compile(
    r"^(?P<stem>.*?)(?:\.(?:lofreq|filtered|csq))*\.vcf(?:\.gz)?$",
    re.IGNORECASE,
)


@dataclass
class SampleRecord:
    sample_name: str
    sample_number: Optional[int]
    reads1: str = ""
    reads2: str = ""
    bam: str = ""
    vcf: str = ""
    coverage: str = ""
    warnings: List[str] = field(default_factory=list)

    def as_row(self) -> List[str]:
        number = "" if self.sample_number is None else str(self.sample_number)
        return [
            self.sample_name,
            number,
            self.reads1,
            self.reads2,
            self.bam,
            self.vcf,
            self.coverage,
        ]


def _extract_sample_number(stem: str) -> Optional[int]:
    match = re.search(r"(\d+)$", stem)
    return int(match.group(1)) if match else None


def _parse_read_file(path: Path) -> tuple[str, Optional[int], Optional[str]]:
    match = READ_PATTERN.match(path.name)
    if not match:
        return path.stem, None, None
    stem = match.group("stem")
    read = match.group("read")
    sample_number = _extract_sample_number(stem)
    return stem, sample_number, read


def _parse_bam_file(path: Path) -> tuple[str, Optional[int]]:
    match = BAM_PATTERN.match(path.name)
    stem = match.group("stem") if match else path.stem
    return stem, _extract_sample_number(stem)


def _parse_vcf_file(path: Path) -> tuple[str, Optional[int]]:
    match = VCF_PATTERN.match(path.name)
    stem = match.group("stem") if match else path.stem
    return stem, _extract_sample_number(stem)


def _find_files(directory: Path, extensions: Iterable[str]) -> List[Path]:
    return sorted(
        p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in extensions
    )


def generate_from_reads(directory: Path) -> List[SampleRecord]:
    reads = [
        p for p in directory.iterdir() if p.is_file() and READ_PATTERN.match(p.name)
    ]
    if not reads:
        raise FileNotFoundError("No FASTQ files found in the provided directory")

    grouped: dict[str, dict[str, Path]] = {}
    numbers: dict[str, Optional[int]] = {}
    for path in reads:
        sample_name, sample_number, read = _parse_read_file(path)
        grouped.setdefault(sample_name, {})
        numbers.setdefault(sample_name, sample_number)
        if sample_number is not None:
            numbers[sample_name] = sample_number
        if read == "2":
            grouped[sample_name]["R2"] = path
        else:
            grouped[sample_name]["R1"] = path
        if numbers[sample_name] is None:
            numbers[sample_name] = sample_number

    records: List[SampleRecord] = []
    for sample, files in grouped.items():
        sample_number = numbers.get(sample)
        reads1_path = files.get("R1")
        reads2_path = files.get("R2")
        reads1 = str(reads1_path.resolve()) if reads1_path else ""
        reads2 = str(reads2_path.resolve()) if reads2_path else ""
        record = SampleRecord(
            sample_name=sample,
            sample_number=sample_number,
            reads1=reads1,
            reads2=reads2,
        )
        if sample_number is None:
            record.warnings.append(
                f"Could not determine passage number for sample '{sample}'"
            )
        if not reads2:
            record.warnings.append(
                f"No R2 read found for sample '{sample}' (assuming SE)"
            )
        records.append(record)

    return sorted(records, key=lambda r: (r.sample_number or 0, r.sample_name))


def generate_from_bams(directory: Path) -> List[SampleRecord]:
    bam_files = [
        p for p in directory.iterdir() if p.is_file() and BAM_PATTERN.match(p.name)
    ]
    if not bam_files:
        raise FileNotFoundError("No BAM files found in the provided directory")

    records: List[SampleRecord] = []
    for path in bam_files:
        sample, number = _parse_bam_file(path)
        record = SampleRecord(
            sample_name=sample,
            sample_number=number,
            bam=str(path.resolve()),
        )
        if number is None:
            record.warnings.append(
                f"Could not determine passage number for sample '{sample}'"
            )
        records.append(record)

    return sorted(records, key=lambda r: (r.sample_number or 0, r.sample_name))


def generate_from_vcfs(directory: Path) -> List[SampleRecord]:
    vcf_files = [
        p for p in directory.iterdir() if p.is_file() and VCF_PATTERN.match(p.name)
    ]
    if not vcf_files:
        raise FileNotFoundError("No VCF files found in the provided directory")

    records: List[SampleRecord] = []
    for path in vcf_files:
        sample, number = _parse_vcf_file(path)
        record = SampleRecord(
            sample_name=sample,
            sample_number=number,
            vcf=str(path.resolve()),
        )
        if number is None:
            record.warnings.append(
                f"Could not determine passage number for sample '{sample}'"
            )
        records.append(record)

    return sorted(records, key=lambda r: (r.sample_number or 0, r.sample_name))


def write_records(records: List[SampleRecord], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
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


def generate_csv(
    mode: str,
    directory: Path,
    output: Path,
    dry_run: bool = False,
) -> tuple[List[SampleRecord], List[str]]:
    directory = directory.expanduser().resolve()
    if not directory.is_dir():
        raise FileNotFoundError(f"Directory not found: {directory}")

    if mode == "e2e":
        records = generate_from_reads(directory)
    elif mode == "bam":
        records = generate_from_bams(directory)
    else:
        records = generate_from_vcfs(directory)

    warnings: List[str] = []
    for record in records:
        warnings.extend(record.warnings or [])

    if dry_run:
        return records, warnings

    write_records(records, output)
    return records, warnings
