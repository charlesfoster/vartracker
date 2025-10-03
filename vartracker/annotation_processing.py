from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple
from urllib.parse import unquote


def _first_fasta_record_id(fasta_path: Path) -> str | None:
    with fasta_path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                return line[1:].split()[0].strip()
    return None


def _first_gff_seqid(gff_path: Path) -> str | None:
    with gff_path.open() as fh:
        for raw in fh:
            if not raw or raw.startswith("#"):
                continue
            seqid = raw.split("\t", 1)[0].strip()
            if seqid.startswith(">"):
                # Some GFF files may embed FASTA headers (e.g., ##FASTA sections)
                seqid = seqid.lstrip(">").split()[0]
            return seqid
    return None


def validate_reference_and_annotation(
    reference: str | Path, gff_path: str | Path
) -> None:
    fasta_id = _first_fasta_record_id(Path(reference))
    gff_id = _first_gff_seqid(Path(gff_path))

    if not fasta_id or not gff_id:
        raise ValueError(
            "Unable to determine sequence IDs from reference or annotation"
        )

    if fasta_id != gff_id:
        raise ValueError(
            f"Reference FASTA ID '{fasta_id}' does not match annotation ID '{gff_id}'"
        )


def _parse_attrs(attr: str) -> dict:
    out = {}
    for part in attr.split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = unquote(v)
    return out


def gene_lengths_from_gff3(gff_path: str | Path) -> Dict[str, int]:
    """
    Parse a GFF3 and return {gene_name: CDS_length}, plus 5'/3' UTR and INTERGENIC=1.

    Robust to feature order (e.g., CDS before mRNA). Prefers:
      1) CDS attribute 'gene'
      2) Parent=gene:<id> → gene 'Name'
      3) Parent=transcript:<id> → mRNA 'Name' or its parent gene 'Name'
      4) As last resort, protein_id/tx id (rare).
    """
    gff_path = Path(gff_path)

    # Maps discovered during scan
    gene_id_to_name: Dict[str, str] = {}  # gene-id -> gene symbol/name
    tx_id_to_gene_name: Dict[str, str] = {}  # transcript-id -> gene name

    # Accumulators
    gene_cds_bp: Dict[str, int] = defaultdict(int)
    pending_by_tx: Dict[str, List[int]] = defaultdict(
        list
    )  # CDS lengths waiting for mRNA mapping

    # For UTR computation (whole-genome)
    seq_range: Tuple[int, int] | None = None
    first_cds_start: int | None = None
    last_cds_end: int | None = None

    # --- First pass: parse everything; resolve when possible; queue when not
    with gff_path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue

            if line.startswith("##"):
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            seq_start = int(parts[2])
                            seq_end = int(parts[3])
                            seq_range = (seq_start, seq_end)
                        except ValueError:
                            pass
                continue
            if line.startswith("#"):
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs = line.split(
                "\t"
            )
            start_i, end_i = int(start), int(end)
            a = _parse_attrs(attrs)

            if ftype == "gene":
                gid = (a.get("ID") or "").split(":", 1)[-1]
                gname = a.get("Name") or a.get("gene") or gid
                gene_id_to_name[gid] = gname

            elif ftype in ("mRNA", "transcript"):
                tid = (a.get("ID") or "").split(":", 1)[-1]
                # Prefer mRNA Name, else parent gene's Name, else 'gene' tag if present
                gname = a.get("Name")
                if not gname:
                    parent_gene_id = (a.get("Parent") or "").split(":", 1)[-1]
                    gname = gene_id_to_name.get(parent_gene_id) or a.get("gene") or tid
                tx_id_to_gene_name[tid] = gname

            elif ftype == "CDS":
                # Compute CDS length
                cds_len = end_i - start_i + 1
                # Track for UTRs
                first_cds_start = (
                    start_i
                    if first_cds_start is None
                    else min(first_cds_start, start_i)
                )
                last_cds_end = (
                    end_i if last_cds_end is None else max(last_cds_end, end_i)
                )

                # Try to resolve gene name immediately
                gname = (
                    a.get("gene") or None  # many RSV/Ncbi GFFs include 'gene=' on CDS
                )
                if not gname:
                    parent = a.get("Parent") or ""
                    if parent.startswith("gene:"):
                        gid = parent.split(":", 1)[-1]
                        gname = gene_id_to_name.get(gid)
                    else:
                        tid = parent.split(":", 1)[-1]
                        gname = tx_id_to_gene_name.get(tid)
                        if not gname:
                            # Defer until after we see the mRNA
                            pending_by_tx[tid].append(cds_len)
                            continue  # will add later once tx->gene is known

                if not gname:
                    # Skip unresolvable CDS entries; they will be handled later
                    continue

                gene_cds_bp[str(gname)] += cds_len

    # --- Resolve any pending CDS whose transcripts appeared after
    for tid, lengths in pending_by_tx.items():
        gname = tx_id_to_gene_name.get(tid)
        if not gname:
            # Last resort: use the transcript id as a stand-in (or skip)
            gname = tid
        if not gname:
            continue
        gene_cds_bp[str(gname)] += sum(lengths)

    # --- Add UTRs if we have genome bounds
    if seq_range and first_cds_start is not None and last_cds_end is not None:
        seq_start, seq_end = seq_range
        five_utr = max(0, first_cds_start - seq_start)  # 1..first_cds_start-1
        three_utr = max(0, seq_end - last_cds_end)  # last_cds_end+1..end
        gene_cds_bp["5' UTR"] = five_utr
        gene_cds_bp["3' UTR"] = three_utr

    # Conventional placeholder
    gene_cds_bp["INTERGENIC"] = 1

    return dict(gene_cds_bp)
