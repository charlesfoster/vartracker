"""Utilities for preparing reference bundles from GenBank accessions."""
# modified based on/inspired by:
# - work by Damien Farrell https://dmnfarrell.github.io/bioinformatics/bcftools-csq-gff-format
# - script bundled with bcftools: https://github.com/samtools/bcftools/blob/develop/misc/gff2gff.py

from __future__ import annotations

import datetime as dt
import hashlib
import json
import re
import shutil
import subprocess
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, SeqFeature

from ._version import __version__
from .provenance import collect_tool_versions, compute_sha256

EUTILS_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_ACCESSION_RE = re.compile(r"^[A-Za-z0-9_.]+$")
_DNA_BASES = ("A", "C", "G", "T")


@dataclass(frozen=True)
class GffFeature:
    """In-memory representation of a single GFF3 feature line."""

    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: tuple[tuple[str, str], ...]

    def to_gff3_line(self) -> str:
        attr_text = ";".join(f"{k}={_gff_escape(v)}" for k, v in self.attributes)
        fields = [
            self.seqid,
            self.source,
            self.feature_type,
            str(self.start),
            str(self.end),
            self.score,
            self.strand,
            self.phase,
            attr_text,
        ]
        return "\t".join(fields)


def _gff_escape(value: str) -> str:
    return urllib.parse.quote(str(value), safe="._:-")


def _md5(path: Path) -> tuple[str | None, str | None]:
    try:
        hasher = hashlib.md5()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                hasher.update(chunk)
        return hasher.hexdigest(), None
    except Exception as exc:  # pragma: no cover - metadata fallback
        return None, str(exc)


def parse_accessions(
    *, accessions: str | None = None, accession_file: str | Path | None = None
) -> list[str]:
    """Parse and validate accession input from CLI values."""
    if bool(accessions) == bool(accession_file):
        raise ValueError(
            "Exactly one of --accessions or --accession-file must be provided."
        )

    raw_values: list[str] = []
    if accessions:
        raw_values = [token.strip() for token in accessions.split(",")]
    else:
        assert accession_file is not None
        file_path = Path(accession_file).expanduser().resolve()
        if not file_path.exists():
            raise FileNotFoundError(f"Accession file not found: {file_path}")
        raw_values = [
            line.strip()
            for line in file_path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.strip().startswith("#")
        ]

    cleaned: list[str] = []
    seen: set[str] = set()
    for item in raw_values:
        if not item:
            continue
        if not _ACCESSION_RE.fullmatch(item):
            raise ValueError(
                f"Invalid accession '{item}'. Accessions may only contain letters, numbers, underscore, and period."
            )
        if item not in seen:
            seen.add(item)
            cleaned.append(item)

    if not cleaned:
        raise ValueError("No accessions were provided after parsing input.")

    return cleaned


def fetch_genbank_record(accession: str, *, timeout_seconds: int = 30) -> str:
    """Download a GenBank flat file for a nucleotide accession from NCBI."""
    query = urllib.parse.urlencode(
        {
            "db": "nuccore",
            "id": accession,
            "rettype": "gbwithparts",
            "retmode": "text",
        }
    )
    url = f"{EUTILS_EFETCH}?{query}"
    request = urllib.request.Request(
        url=url,
        headers={"User-Agent": f"vartracker/{__version__} (reference prepare)"},
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
            payload = response.read().decode("utf-8")
    except Exception as exc:
        raise RuntimeError(
            f"Failed to fetch GenBank record for accession '{accession}': {exc}"
        ) from exc

    if "LOCUS" not in payload:
        raise RuntimeError(
            f"Failed to fetch GenBank record for accession '{accession}': response does not look like GenBank text."
        )
    return payload


def _feature_name(feature: SeqFeature, idx: int) -> str:
    qualifiers = feature.qualifiers
    for key in ("gene", "locus_tag", "product", "protein_id"):
        values = qualifiers.get(key)
        if values and values[0]:
            return str(values[0])
    return f"feature_{idx}"


def _strand_symbol(feature: SeqFeature) -> str:
    strand = feature.location.strand
    if strand == 1:
        return "+"
    if strand == -1:
        return "-"
    return "."


def _feature_parts(feature: SeqFeature) -> list[tuple[int, int]]:
    loc = feature.location
    if isinstance(loc, CompoundLocation):
        parts = loc.parts
    else:
        parts = [loc]
    return [(int(part.start) + 1, int(part.end)) for part in parts]


def _cds_phases(feature: SeqFeature, parts: Sequence[tuple[int, int]]) -> dict[tuple[int, int], str]:
    codon_start_raw = feature.qualifiers.get("codon_start", ["1"])[0]
    try:
        codon_start = int(codon_start_raw)
    except ValueError:
        codon_start = 1
    phase = max(0, min(2, codon_start - 1))

    strand = feature.location.strand
    ordered = sorted(parts, key=lambda x: x[0], reverse=bool(strand == -1))
    out: dict[tuple[int, int], str] = {}
    for start, end in ordered:
        out[(start, end)] = str(phase)
        seg_len = max(0, end - start + 1)
        phase = (3 - ((seg_len - phase) % 3)) % 3
    return out


def convert_genbank_to_features(
    gb_path: str | Path, *, seqid: str
) -> tuple[str, list[GffFeature]]:
    """Convert a GenBank record into sequence + bcftools-csq-friendly GFF features."""
    gb_file = Path(gb_path).expanduser().resolve()
    record = SeqIO.read(str(gb_file), "genbank")
    sequence = str(record.seq).upper()
    if not sequence:
        raise ValueError(f"No sequence found in GenBank record: {gb_file}")

    features: list[GffFeature] = []
    source = "vartracker_prepare"
    cds_index = 1

    for feature in record.features:
        if feature.type != "CDS":
            continue

        name = _feature_name(feature, cds_index)
        strand = _strand_symbol(feature)
        parts = _feature_parts(feature)
        if not parts:
            continue

        gene_id = f"gene:{seqid}_{cds_index}"
        transcript_id = f"transcript:{seqid}_{cds_index}"
        cds_base_id = f"cds:{seqid}_{cds_index}"
        starts = [start for start, _ in parts]
        ends = [end for _, end in parts]
        gene_start = min(starts)
        gene_end = max(ends)

        features.append(
            GffFeature(
                seqid=seqid,
                source=source,
                feature_type="gene",
                start=gene_start,
                end=gene_end,
                score=".",
                strand=strand,
                phase=".",
                attributes=(("ID", gene_id), ("Name", name), ("biotype", "protein_coding")),
            )
        )
        features.append(
            GffFeature(
                seqid=seqid,
                source=source,
                feature_type="mRNA",
                start=gene_start,
                end=gene_end,
                score=".",
                strand=strand,
                phase=".",
                attributes=(
                    ("ID", transcript_id),
                    ("Parent", gene_id),
                    ("Name", name),
                    ("biotype", "protein_coding"),
                ),
            )
        )

        phase_map = _cds_phases(feature, parts)
        for part_idx, (start, end) in enumerate(sorted(parts, key=lambda x: (x[0], x[1])), start=1):
            cds_id = cds_base_id if len(parts) == 1 else f"{cds_base_id}.{part_idx}"
            attrs: list[tuple[str, str]] = [("ID", cds_id), ("Parent", transcript_id), ("gene", name)]
            protein_ids = feature.qualifiers.get("protein_id")
            if protein_ids and protein_ids[0]:
                attrs.append(("protein_id", str(protein_ids[0])))
            features.append(
                GffFeature(
                    seqid=seqid,
                    source=source,
                    feature_type="CDS",
                    start=start,
                    end=end,
                    score=".",
                    strand=strand,
                    phase=phase_map[(start, end)],
                    attributes=tuple(attrs),
                )
            )

        cds_index += 1

    if not any(f.feature_type == "CDS" for f in features):
        raise ValueError(f"No CDS features found in GenBank record: {gb_file}")

    return sequence, features


def _write_fasta_with_index(records: Sequence[tuple[str, str]], fasta_path: Path) -> Path:
    wrap = 60
    index_rows: list[tuple[str, int, int, int, int]] = []

    with fasta_path.open("wb") as handle:
        for seqid, sequence in records:
            header = f">{seqid}\n".encode("ascii")
            handle.write(header)
            offset = handle.tell()

            seq_len = len(sequence)
            line_blen = wrap if seq_len > wrap else max(1, seq_len)
            line_len = line_blen + 1
            for i in range(0, seq_len, wrap):
                handle.write(sequence[i : i + wrap].encode("ascii"))
                handle.write(b"\n")
            index_rows.append((seqid, seq_len, offset, line_blen, line_len))

    fai_path = Path(f"{fasta_path}.fai")
    with fai_path.open("w", encoding="utf-8") as handle:
        for row in index_rows:
            handle.write("\t".join(str(x) for x in row))
            handle.write("\n")

    return fai_path


def _write_merged_gff3(features: Sequence[GffFeature], gff_path: Path) -> None:
    ordered = sorted(
        features,
        key=lambda f: (
            f.seqid,
            f.start,
            f.end,
            f.feature_type,
            tuple(f.attributes),
        ),
    )
    with gff_path.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for feature in ordered:
            handle.write(feature.to_gff3_line())
            handle.write("\n")


def _first_cds(gff_path: Path) -> tuple[str, int]:
    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            if fields[2] != "CDS":
                continue
            return fields[0], int(fields[3])
    raise ValueError("No CDS feature found in merged GFF3 for csq validation.")


def _get_base_at(fasta_path: Path, seqid: str, position_1based: int) -> str:
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        if record.id == seqid:
            if position_1based < 1 or position_1based > len(record.seq):
                raise ValueError(
                    f"Validation position {position_1based} is out of range for {seqid}."
                )
            return str(record.seq[position_1based - 1]).upper()
    raise ValueError(f"Sequence ID '{seqid}' not found in FASTA.")


def _choose_alt(ref_base: str) -> str:
    base = ref_base.upper()
    for candidate in _DNA_BASES:
        if candidate != base:
            return candidate
    return "A"


def _write_dummy_vcf(
    path: Path, *, seqid: str, seq_len: int, pos: int, ref: str, alt: str
) -> None:
    text = (
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID={seqid},length={seq_len}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"{seqid}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n"
    )
    path.write_text(text, encoding="utf-8")


def _fasta_lengths(fasta_path: Path) -> dict[str, int]:
    return {record.id: len(record.seq) for record in SeqIO.parse(str(fasta_path), "fasta")}


def validate_csq_with_dummy_variant(
    *, fasta_path: str | Path, gff_path: str | Path, workdir: str | Path
) -> dict[str, object]:
    """Run a bcftools-csq smoke test against the first CDS feature."""
    fasta = Path(fasta_path).expanduser().resolve()
    gff = Path(gff_path).expanduser().resolve()
    outdir = Path(workdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    if shutil.which("bcftools") is None:
        raise RuntimeError(
            "bcftools is required for csq validation but was not found on PATH. "
            "Install bcftools or rerun with --skip-csq-validation."
        )

    seqid, pos = _first_cds(gff)
    lengths = _fasta_lengths(fasta)
    if seqid not in lengths:
        raise ValueError(
            f"SeqID mismatch: CDS seqid '{seqid}' from GFF3 does not exist in FASTA."
        )
    ref = _get_base_at(fasta, seqid, pos)
    alt = _choose_alt(ref)

    dummy_vcf = outdir / "csq_validation.vcf"
    dummy_out = outdir / "csq_validation.annotated.vcf"
    _write_dummy_vcf(
        dummy_vcf, seqid=seqid, seq_len=lengths[seqid], pos=pos, ref=ref, alt=alt
    )

    cmd = [
        "bcftools",
        "csq",
        "-f",
        str(fasta),
        "-g",
        str(gff),
        str(dummy_vcf),
        "-Ov",
        "-o",
        str(dummy_out),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        stderr_snippet = (result.stderr or "").strip().splitlines()[-10:]
        snippet_text = "\n".join(stderr_snippet) if stderr_snippet else "No stderr output."
        raise RuntimeError(
            "bcftools csq validation failed.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Stderr snippet:\n{snippet_text}"
        )

    found_bcsq = False
    for line in dummy_out.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) >= 8 and "BCSQ=" in fields[7]:
            found_bcsq = True
            break

    if not found_bcsq:
        raise RuntimeError(
            "bcftools csq validation failed: output VCF does not contain BCSQ annotation for the dummy variant."
        )

    return {
        "status": "passed",
        "command": " ".join(cmd),
        "dummy_vcf": str(dummy_vcf),
        "dummy_output_vcf": str(dummy_out),
    }


def prepare_reference_bundle(
    *,
    accessions: Sequence[str],
    outdir: str | Path,
    prefix: str = "reference",
    force: bool = False,
    keep_intermediates: bool = False,
    skip_csq_validation: bool = False,
    invocation: str | None = None,
    argv: Sequence[str] | None = None,
    fetcher: Callable[[str], str] = fetch_genbank_record,
) -> dict[str, object]:
    """Build merged FASTA/GFF3 reference files from GenBank nucleotide accessions."""
    if not accessions:
        raise ValueError("At least one accession is required.")

    output_dir = Path(outdir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = output_dir / f"{prefix}.fa"
    gff_path = output_dir / f"{prefix}.gff3"
    fai_path = output_dir / f"{prefix}.fa.fai"
    metadata_path = output_dir / "prepare_metadata.json"

    for path in (fasta_path, gff_path, fai_path, metadata_path):
        if path.exists() and not force:
            raise FileExistsError(
                f"Output already exists: {path}. Use --force to overwrite existing files."
            )

    intermediates_dir = output_dir / "intermediates"
    intermediates_dir.mkdir(parents=True, exist_ok=True)

    merged_fasta: list[tuple[str, str]] = []
    merged_features: list[GffFeature] = []
    seen_seqids: set[str] = set()
    per_accession_outputs: list[dict[str, str]] = []

    for accession in accessions:
        gb_text = fetcher(accession)
        gb_path = intermediates_dir / f"{accession}.gb"
        gb_path.write_text(gb_text, encoding="utf-8")

        sequence, features = convert_genbank_to_features(gb_path, seqid=accession)
        if accession in seen_seqids:
            raise ValueError(f"Duplicate seqid detected after conversion: {accession}")
        seen_seqids.add(accession)
        merged_fasta.append((accession, sequence))
        merged_features.extend(features)

        per_item: dict[str, str] = {"accession": accession}
        if keep_intermediates:
            fasta_intermediate = intermediates_dir / f"{accession}.fa"
            gff_intermediate = intermediates_dir / f"{accession}.gff3"
            fasta_intermediate.write_text(
                f">{accession}\n{sequence}\n", encoding="utf-8"
            )
            _write_merged_gff3([x for x in features if x.seqid == accession], gff_intermediate)
            per_item.update(
                {
                    "genbank": str(gb_path),
                    "fasta": str(fasta_intermediate),
                    "gff3": str(gff_intermediate),
                }
            )

        per_accession_outputs.append(per_item)
        if not keep_intermediates:
            gb_path.unlink(missing_ok=True)

    _write_fasta_with_index(merged_fasta, fasta_path)
    _write_merged_gff3(merged_features, gff_path)

    if not fai_path.exists():
        raise RuntimeError(f"Failed to create FASTA index at {fai_path}")

    csq_validation: dict[str, object]
    if skip_csq_validation:
        csq_validation = {"status": "skipped", "reason": "Skipped by --skip-csq-validation"}
    else:
        csq_validation = validate_csq_with_dummy_variant(
            fasta_path=fasta_path,
            gff_path=gff_path,
            workdir=intermediates_dir if keep_intermediates else output_dir,
        )

    fa_sha256, fa_sha_err = compute_sha256(fasta_path)
    gff_sha256, gff_sha_err = compute_sha256(gff_path)
    fa_md5, fa_md5_err = _md5(fasta_path)
    gff_md5, gff_md5_err = _md5(gff_path)

    metadata: dict[str, object] = {
        "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
        "vartracker_version": __version__,
        "command": invocation,
        "args": list(argv or []),
        "accessions": list(accessions),
        "outputs": {
            "fasta": str(fasta_path),
            "gff3": str(gff_path),
            "fai": str(fai_path),
            "metadata": str(metadata_path),
        },
        "checksums": {
            "fasta": {"sha256": fa_sha256, "md5": fa_md5},
            "gff3": {"sha256": gff_sha256, "md5": gff_md5},
        },
        "csq_validation": csq_validation,
        "external_tools": collect_tool_versions(("bcftools", "samtools", "tabix")),
        "per_accession_outputs": per_accession_outputs,
        "keep_intermediates": keep_intermediates,
    }
    checksum_errors: dict[str, dict[str, str | None]] = {}
    if fa_sha_err or fa_md5_err:
        checksum_errors["fasta"] = {
            "sha256_error": fa_sha_err,
            "md5_error": fa_md5_err,
        }
    if gff_sha_err or gff_md5_err:
        checksum_errors["gff3"] = {
            "sha256_error": gff_sha_err,
            "md5_error": gff_md5_err,
        }
    if checksum_errors:
        metadata["checksum_errors"] = checksum_errors

    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")
    if not keep_intermediates:
        try:
            intermediates_dir.rmdir()
        except OSError:
            pass
    return metadata


__all__ = [
    "GffFeature",
    "convert_genbank_to_features",
    "fetch_genbank_record",
    "parse_accessions",
    "prepare_reference_bundle",
    "validate_csq_with_dummy_variant",
]
