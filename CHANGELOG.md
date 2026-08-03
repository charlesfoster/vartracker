# Changelog

All notable changes to vartracker are documented in this file, in addition to
(not replacing) the project's versioned GitHub releases.

Format loosely follows [Keep a Changelog](https://keepachangelog.com/).

## [2.3.0] - 2026-07-31

### Added

- `--local-csq` option (`vartracker vcf`/`bam`/`end-to-end`) to switch `bcftools csq` to
  independent, per-variant consequence calling instead of the default joint/compound calling -
  useful for gene-dense, high-variant-density data where unphased, sub-consensus variants can
  cluster in the same gene without genotype evidence they co-occur.
- `--out`/`--outdir` for `vartracker plot heatmap`, to control the output filename and location.
- `--multiallelic-overflow` (`error`/`drop-lowest-af`/`skip-site`) to control how sites with more
  than two surviving ALT alleles in a single sample are handled before `bcftools csq`.
- LoFreq primer-overlap rescue for amplicon data (`--primer-bed`, `--lofreq-primer-rescue`,
  `--lofreq-rescue-*` thresholds), with raw, rescued, and filtered-out LoFreq calls written to
  per-sample audit TSVs.
- Consensus genome generation in the `bam`/`end-to-end` workflows (`--consensus-snp-min-af`,
  `--consensus-snp-thresh`, `--consensus-indel-thresh`), writing both a simple and an
  IUPAC-aware consensus FASTA per sample.
- Gene-name disambiguation across contigs/replicons: `ambiguous_gene_names()` and updates to
  `gene_lengths_from_gff3`/`generate_gene_table` so genes that legitimately reuse a name on
  different contigs (e.g. `repA` on a chromosome and a plasmid) are tracked and reported separately.
- Gene-wise summary plot (`mutations_per_gene.pdf`) is now capped to the top genes by newly
  emerged variants on large (e.g. bacterial) genomes, via `--max-plot-genes`/`--plot-genes`; the
  tabular output is unaffected.
- New standalone plotting subcommands (`plot heatmap`, `plot genome`, `plot trajectory`,
  `plot turnover`, `plot lifespan`) and heatmap filtering options (`--aa-exclude`, `--aa-include`,
  `--only-persistent`, `--only-new`, `--gene-include`/`--gene-exclude`, `--variant-type`, `--qc`,
  `--min-prop-passing-qc`, `--min-persistence`, `--min-max-af`, `--min-sample-af`,
  `--sample-subset`, `--hide-singletons`, `--min-depth`, `--x-labels`, `--title`,
  `--literature-csv`).
- Structured logging for Snakemake pipeline steps.
- Documentation: a guide for building a custom `--literature-csv` for non-SARS-CoV-2 pathogens, a
  QC-column interpretation guide, and a bacterial-genome validation/limitations section.

### Changed

- `vartracker plot heatmap` CLI flags were de-prefixed for standalone use (e.g.
  `--heatmap-aa-exclude` -> `--aa-exclude`); the legacy prefixed forms are kept as aliases.
- The default heatmap now always shows a variant's canonical row, even if that row is
  joint/compound, instead of unconditionally hiding all joint rows; `--include-joint` reveals
  any additional joint/compound annotation-group rows a variant has.
- Per-sample variant QC is now reported directly (`per_sample_variant_qc`, `P`/`F`) rather than
  only as an overall pass/fail, and the default heatmap visually flags failing cells.

### Fixed

- Duplicate `results.csv` rows for the same variant (one per `bcftools csq` annotation group) were
  previously plotted as independent data points, inflating turnover counts; plotting now collapses
  to one row per (variant, sample), taking the union of presence and the maximum allele frequency.
- `joint_variant` could miss rows with a genuinely compound `bcsq_nt_notation` when a variant had
  several sibling rows at the same position; it now flags every such row.
- Heatmap label truncation could make two distinct variants render an identical label, silently
  dropping one of them; colliding labels are now disambiguated with a positional suffix.
- Gene names that collide across contigs/replicons were previously merged into a single row in
  gene-length and gene-wise summary tables; they are now tracked and reported per-contig.
- Percent-encoded gene names read from `bcftools`' `BCSQ` field (e.g. `ercS%27`) were not
  URL-decoded on the variant-calling side, so they never matched the (correctly decoded)
  annotation side and were silently dropped from gene-wise summaries.
- Corrected several low-frequency multiallelic-site handling bugs, so distinct ALT alleles at the
  same position are kept separate through preprocessing and correctly rejoined before
  `bcftools csq`.
- Fixed potential results-row duplication for some samples during heatmap generation.
- Fixed longitudinal variant-rescue edge cases and LoFreq strand-bias rescue thresholds.
