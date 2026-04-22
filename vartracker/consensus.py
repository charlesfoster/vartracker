"""Consensus allele/genotype helpers for the Snakemake workflow."""

from __future__ import annotations

IUPAC_CODES: dict[frozenset[str], str] = {
    frozenset({"A", "G"}): "R",
    frozenset({"C", "T"}): "Y",
    frozenset({"G", "C"}): "S",
    frozenset({"A", "T"}): "W",
    frozenset({"G", "T"}): "K",
    frozenset({"A", "C"}): "M",
}


def first_alt_allele(alt: str) -> str:
    """Return the first ALT allele from a possibly multi-allelic VCF field."""
    return str(alt).split(",", 1)[0]


def is_snp_variant(ref: str, alt: str) -> bool:
    """Return True when the first ALT allele represents a SNP."""
    return len(str(ref)) == 1 and len(first_alt_allele(alt)) == 1


def consensus_genotype_for_variant(
    ref: str,
    alt: str,
    *,
    depth: int,
    af: float,
    min_depth: int,
    consensus_snp_min_af: float,
    consensus_snp_thresh: float,
    consensus_indel_thresh: float,
) -> str:
    """Return the genotype bcftools consensus should use for IUPAC output."""
    if depth < min_depth:
        return "./."

    if is_snp_variant(ref, alt):
        if af >= consensus_snp_thresh:
            return "1/1"
        if af >= consensus_snp_min_af:
            return "0/1"
        return "./."

    if af >= consensus_indel_thresh:
        return "1/1"
    return "./."


def simple_consensus_base(
    ref: str,
    alt: str,
    *,
    depth: int,
    af: float,
    min_depth: int,
    consensus_snp_min_af: float,
    consensus_snp_thresh: float,
    consensus_indel_thresh: float,
) -> str:
    """Return REF or ALT for simple consensus tests."""
    if depth < min_depth:
        return ref

    if is_snp_variant(ref, alt):
        if af >= consensus_snp_thresh:
            return first_alt_allele(alt)
        return ref

    if af >= consensus_indel_thresh:
        return first_alt_allele(alt)
    return ref


def iupac_consensus_base(
    ref: str,
    alt: str,
    *,
    depth: int,
    af: float,
    min_depth: int,
    consensus_snp_min_af: float,
    consensus_snp_thresh: float,
    consensus_indel_thresh: float,
) -> str:
    """Return the expected single-site IUPAC consensus symbol for tests."""
    genotype = consensus_genotype_for_variant(
        ref,
        alt,
        depth=depth,
        af=af,
        min_depth=min_depth,
        consensus_snp_min_af=consensus_snp_min_af,
        consensus_snp_thresh=consensus_snp_thresh,
        consensus_indel_thresh=consensus_indel_thresh,
    )
    if genotype == "./.":
        return ref

    first_alt = first_alt_allele(alt)
    if genotype == "1/1":
        return first_alt

    if is_snp_variant(ref, alt):
        return IUPAC_CODES.get(frozenset({ref.upper(), first_alt.upper()}), "N")
    return ref
