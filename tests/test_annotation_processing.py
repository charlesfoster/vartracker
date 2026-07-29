"""Tests for annotation_processing helpers, especially multi-contig handling."""

from __future__ import annotations

from vartracker.annotation_processing import (
    ambiguous_gene_names,
    gene_lengths_from_gff3,
)

SINGLE_CONTIG_GFF3 = """\
##gff-version 3
chrom1\tcustom\tCDS\t1\t300\t.\t+\t0\tID=cds-dnaA;gene=dnaA
chrom1\tcustom\tCDS\t400\t600\t.\t+\t0\tID=cds-repA;gene=repA
"""

MULTI_CONTIG_COLLISION_GFF3 = """\
##gff-version 3
chrom1\tcustom\tCDS\t1\t300\t.\t+\t0\tID=cds-dnaA;gene=dnaA
chrom1\tcustom\tCDS\t400\t600\t.\t+\t0\tID=cds-repA-chrom;gene=repA
plasmid1\tcustom\tCDS\t1\t150\t.\t+\t0\tID=cds-repA-plasmid;gene=repA
plasmid1\tcustom\tCDS\t200\t350\t.\t+\t0\tID=cds-mobA;gene=mobA
"""


def test_gene_lengths_single_contig_unqualified(tmp_path):
    """Single-contig references (e.g. viral genomes) keep plain gene names."""
    gff_path = tmp_path / "single.gff3"
    gff_path.write_text(SINGLE_CONTIG_GFF3)

    lengths = gene_lengths_from_gff3(gff_path)

    assert lengths["dnaA"] == 300
    assert lengths["repA"] == 201
    assert ambiguous_gene_names(gff_path) == set()


def test_gene_lengths_disambiguates_colliding_names_across_contigs(tmp_path):
    """Identically-named genes on different contigs/replicons must not be summed."""
    gff_path = tmp_path / "multi.gff3"
    gff_path.write_text(MULTI_CONTIG_COLLISION_GFF3)

    lengths = gene_lengths_from_gff3(gff_path)

    # Colliding gene name is split per-contig, each keeping its own length.
    assert lengths["repA (chrom1)"] == 201
    assert lengths["repA (plasmid1)"] == 150
    assert "repA" not in lengths

    # Non-colliding gene names on either contig stay unqualified.
    assert lengths["dnaA"] == 300
    assert lengths["mobA"] == 151

    assert ambiguous_gene_names(gff_path) == {"repA"}
