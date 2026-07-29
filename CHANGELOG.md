# Changelog

This file tracks changes made during the NARGAB manuscript revision cycle, in
addition to (not replacing) the project's versioned GitHub releases.

Format loosely follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Fixed

- **Multi-contig gene-name collisions no longer silently merged.** Bacterial
  reference genomes commonly reuse gene names across replicons (e.g. `repA`
  on both a chromosome and a plasmid). Previously, `gene_lengths_from_gff3`
  and `generate_gene_table` aggregated CDS lengths and mutation statistics by
  gene name alone, so two unrelated genes sharing a name on different contigs
  had their lengths and variant counts silently summed into one row.
  - `gene_lengths_from_gff3` (`annotation_processing.py`) now tracks CDS
    lengths per `(contig, gene)` internally and only disambiguates the
    output label (e.g. `repA (chrom1)` / `repA (plasmid1)`) when a gene name
    actually collides across contigs. Gene names confined to a single contig
    are returned unqualified, so single-contig (viral) references and
    existing multi-segment references with uniquely-named segments (e.g.
    8-segment influenza) are unaffected.
  - `generate_gene_table` (`analysis.py`) now splits colliding gene names by
    `chrom` in the same way, using the results table's own `chrom` column.
  - Added `ambiguous_gene_names()` (`annotation_processing.py`), which
    identifies gene names that collide across contigs directly from the
    reference annotation. This closes an edge case where a colliding gene
    only has variants on one of its contigs in a given dataset: without it,
    that gene's real data could be silently dropped when merged against the
    gene-length scaffold, because the scaffold would only contain the
    contig-qualified labels while the variants table only produced the bare,
    unqualified label.
  - Added test coverage: `tests/test_annotation_processing.py` (new file)
    and three new tests in `tests/test_analysis.py` covering the
    single-contig baseline, the colliding-gene-name split, and the
    ambiguous-genes data-loss edge case.
  - No changes to output schema, column names, or behaviour for any
    single-contig or uniquely-named multi-segment reference.
