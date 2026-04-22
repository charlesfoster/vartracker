[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![PyPI version](https://img.shields.io/pypi/v/vartracker.svg)](https://pypi.org/project/vartracker/)

```
██    ██  █████  ██████  ████████ ██████   █████   ██████ ██   ██ ███████ ██████
██    ██ ██   ██ ██   ██    ██    ██   ██ ██   ██ ██      ██  ██  ██      ██   ██
██    ██ ███████ ██████     ██    ██████  ███████ ██      █████   █████   ██████
 ██  ██  ██   ██ ██   ██    ██    ██   ██ ██   ██ ██      ██  ██  ██      ██   ██
  ████   ██   ██ ██   ██    ██    ██   ██ ██   ██  ██████ ██   ██ ███████ ██   ██


```
# vartracker

A bioinformatics pipeline to summarise variants called against a reference in a longitudinal study design. Written to investigate longitudinal sequencing data from long-term passaging of SARS-CoV-2. However, with appropriate reference data it can be expanded to other pathogens too.

**Author:** Dr Charles Foster

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Output](#output)
- [What does vartracker do?](#what-does-vartracker-do)
- [Citation](#citation)
- [License](#license)
- [Contributing](#contributing)
- [Support](#support)

## Features

- Track mutation persistence across longitudinal samples
- Comprehensive variant analysis including amino acid consequences
- Built-in SARS-CoV-2 reference data and annotations
- Integration with functional mutation databases (literature)
- Automated plotting and statistical analysis
- Support for both SNPs and indels
- Quality control metrics for variants

## Installation

The simplest, and preferred, installation route is via a conda-compatible package manager (pixi, conda, or mamba).

## Conda/Mamba
`vartracker` and its external bioinformatics dependencies can be installed from the package channels directly:

```bash
mamba create -n vartracker -c conda-forge -c bioconda vartracker
mamba activate vartracker
```

If you prefer `conda`:

```bash
conda create -n vartracker -c conda-forge -c bioconda vartracker
conda activate vartracker
```

### Pixi

You can also install `vartracker` via `pixi` either globally or inside an existing workspace:

```bash
# global installation
pixi global install vartracker

# add into an existing workspace
pixi workspace channel add conda-forge
pixi workspace channel add bioconda
pixi add vartracker
```

### Biocontainers

Every Bioconda package is available as a container image for usage with your preferred container runtime. An example command to pull `vartracker` with `docker`:

```bash
# latest build
docker pull quay.io/biocontainers/vartracker:latest
# specific tag
docker pull quay.io/biocontainers/vartracker:<tag>
```

### Alternative: PyPI (Python-only)

If you want a Python-only install (requires Python 3.11 or newer), you can still install from PyPI. In that case you must provide the required external bioinformatics tools yourself (see below):

```bash
pip install vartracker
```

### External Dependencies

vartracker shells out to a handful of bioinformatics tools. Make sure they are discoverable on `PATH` before running the CLI.
Minimum tested versions are tracked in `docs/DEPENDENCIES.md`.

- **bcftools** and **tabix** – required for all modes
- **samtools**, **lofreq**, **fastp**, **bwa**, and **snakemake** – required for the `bam` and `end-to-end` Snakemake workflows

If you only plan to run `vartracker vcf` against pre-generated VCFs, the first pair is sufficient. The additional tools are needed whenever you ask vartracker to align reads or call variants for you.
Consensus genome generation in the `bam` and `end-to-end` workflows uses `bcftools` and `samtools`; it does not require `bedtools`.

Note: the pinned micromamba environment installs `tabix`/`bgzip` via `htslib`.

#### Installing bcftools and tabix

**On macOS:**
```bash
# Using Homebrew
brew install bcftools htslib samtools fastp bwa
# lofreq is available via bioconda (requires conda/mamba)
conda install -c bioconda lofreq

# Using MacPorts
sudo port install bcftools htslib samtools fastp bwa
```

**On Linux (Ubuntu/Debian):**
```bash
sudo apt-get update
sudo apt-get install bcftools tabix samtools fastp bwa
# lofreq is easiest to install via bioconda on Debian-based systems:
conda install -c bioconda lofreq
```

**On Linux (CentOS/RHEL/Fedora):**
```bash
# CentOS/RHEL with EPEL
sudo yum install epel-release
sudo yum install bcftools htslib samtools fastp bwa

# Fedora
sudo dnf install bcftools htslib samtools fastp bwa
# Install lofreq via bioconda on RPM-based systems:
conda install -c bioconda lofreq
```

**Using conda:**
```bash
conda install -c bioconda bcftools samtools tabix fastp bwa lofreq
```

### Development Installation

For development or to get the latest version (requires Python 3.11+):

```bash
git clone https://github.com/charlesfoster/vartracker.git
cd vartracker
pip install -e .[dev]
pre-commit install
```

### Docker: self-build

Build a container image that bundles Python, vartracker, and all external bioinformatics tools:

```bash
  # released version on Bioconda
  docker build -t vartracker:release .
  # development version
  docker build -t vartracker:dev -f Dockerfile.dev .
```

Docker is a self-contained reproducible option. If you publish the image, record the digest and
set it when running to include it in the run manifest:

```bash
export VARTRACKER_CONTAINER_IMAGE=ghcr.io/your-org/vartracker:2.0.0
export VARTRACKER_CONTAINER_DIGEST=sha256:...
```

Run workflows by mounting your data directory into the container. The command below analyses an input CSV located in the current directory and writes results beside it:

```bash
docker run --rm -v "$(pwd)":/workspace vartracker \
  vcf /workspace/inputs/vcf_inputs.csv \
  --outdir /workspace/results
```

## Quick Start

After installation, `vartracker` will be available as a command-line tool:

```bash
vartracker --help
```

### Typical commands

```bash
# Analyse pre-called VCFs plus coverage files
vartracker vcf path/to/vcf_inputs.csv --outdir results/vcf_run

# Run BAMs through the Snakemake workflow, then summarise variants
vartracker bam path/to/bam_inputs.csv \
  --snakemake-outdir work/bam_pipeline \
  --outdir results/bam_summary

# Start from raw reads (FASTQ) and run the full pipeline
vartracker end-to-end path/to/read_inputs.csv \
  --cores 12 \
  --outdir results/e2e_summary

# Re-plot a heatmap from an existing vartracker results file
vartracker plot heatmap results/results.csv \
  --aa-exclude "*frameshift*" \
  --x-labels sample-number \
  --literature-csv results/sample.literature_database_hits.full.csv \
  --title "Variant allele frequencies"

# Plot whole-dataset turnover from an existing results file
vartracker plot turnover results/results.csv

# Plot collapsed variant frequencies along the genome
vartracker plot genome results/results.csv

# Zoom to a gene region, optionally using amino-acid coordinates
vartracker plot genome results/results.csv --gene F --aa-scale

# Plot selected variant trajectories from an existing results file
vartracker plot trajectory results/results.csv \
  --variants "S:D614G,S:E484K,S:N501Y"

# Plot takeover-style trajectories using AF thresholds
vartracker plot trajectory results/results.csv \
  --thresholds 0.5,0.9 \
  --crossing-only

# Generate a template spreadsheet for a directory of files
vartracker prepare spreadsheet --mode e2e --dir data/passaging --out inputs.csv

# Build a reference FASTA+GFF3 bundle from GenBank accessions
vartracker prepare reference --accessions CY114381,CY114382 --outdir refs/flu --prefix flu_ref

# Exercise the bundled smoke-test dataset
vartracker vcf --test
vartracker bam --test
vartracker end-to-end --test
```

All modes understand `--test`, which copies the example dataset from `vartracker/test_data`
into a temporary directory, resolves relative paths, and runs the appropriate workflow.

Temporary LoFreq note:
- In `bam` and `end-to-end` mode, `vartracker` currently caps `lofreq call-parallel` at 8 threads even if `--cores` is higher.
- This is a temporary workaround for an older Bioconda LoFreq build that can fail during `call-parallel` final filtering when many shards produce an excessively long merged VCF header.
- The cap will be revisited once an updated LoFreq build is available through Bioconda.

### Input Spreadsheets

Every CLI mode reads the same canonical columns:

- `sample_name` (required) – display name for the sample
- `sample_number` (required) – passage/order index used in longitudinal plots
- `reads1`, `reads2` – FASTQ paths (required for `end-to-end`, optional elsewhere). The pipeline runs in single-end mode (leave the `reads2` column empty) but the results are less well tested.
- `bam` – BAM file aligned against the SARS-CoV-2 reference
- `vcf` – bgzipped VCF containing variant calls with depth (`DP`) and allele-frequency tags
- `coverage` – per-base coverage TSV with columns `reference<TAB>position<TAB>depth`

Mode-specific expectations:

- **VCF mode** requires `vcf` and `coverage`, while leaving `reads*`/`bam` empty.
- **BAM mode** requires `bam` and will fill `vcf` + `coverage` during the workflow.
- **End-to-end mode** requires `reads1` (and optionally `reads2`); remaining fields are generated.

The `bam` and `end-to-end` workflows also write two consensus FASTA columns to
the updated Snakemake spreadsheet: `consensus` for a simple consensus and
`iupac_consensus` for an IUPAC-aware consensus. SNPs below
`--consensus-snp-min-af` are ignored, SNPs from `--consensus-snp-min-af` up to
`--consensus-snp-thresh` stay as reference bases in the simple consensus
and become REF+ALT ambiguity codes in the IUPAC consensus, and SNPs at or above
`--consensus-snp-thresh` become ALT bases. Indels are controlled
separately by `--consensus-indel-thresh` in both consensus modes. Low-depth bases
are masked as `N`, except for called deletion intervals so true deletions are not
converted to low-depth masks.

Relative paths are resolved with respect to the CSV location, so you can store the sheet alongside
your sequencing artefacts. The `prepare spreadsheet` subcommand can scaffold a CSV and highlight missing files.

Coverage files can be produced with `samtools depth -aa sample.bam > sample_depth.txt` or
`bedtools genomecov -ibam sample.bam -d`. The file name suffix does not matter; vartracker checks
for both `.depth.txt` and `_depth.txt` patterns when preparing its internal test dataset.

### Mode-specific options

- `vartracker vcf` – accepts core analysis options such as `--min-snv-freq`, `--min-indel-freq`,
  `--allele-frequency-tag`, `--name`, `--outdir`, `--sample-cap`, `--manifest-level`, and literature controls
  (`--search-pokay`, `--literature-csv`). Use `--test` to run the bundled smoke test.
- `vartracker bam` – everything from `vcf`, plus Snakemake options:
  `--snakemake-outdir`, `--cores`, `--snakemake-dryrun`, `--verbose`, `--redo`,
  `--rulegraph`, `--consensus-snp-min-af`, `--consensus-snp-thresh`, and
  `--consensus-indel-thresh`.
- `vartracker end-to-end` – similar to `bam`, with optional amplicon clipping controls:
  `--primer-bed` and `--ampliconclip-tolerance` (default: `1`).
- `vartracker plot heatmap` (`hm`) – regenerate the heatmap from an existing vartracker results CSV, including all heatmap customization filters.
- `vartracker plot genome` – plot SNP positions along the genome or a selected gene region using all observed allele-frequency values for each variant.
- `vartracker plot trajectory` – plot allele-frequency trajectories for a selected or auto-ranked subset of variants, optionally in takeover mode using threshold lines and threshold-based filtering.
- `vartracker plot turnover` – plot new-versus-lost longitudinal turnover from the filtered result set.
- `vartracker plot lifespan` – plot first-to-last detection spans for a selected or auto-ranked subset of variants.

Heatmap filtering:
- `vcf`, `bam`, and `end-to-end` always write the default heatmap. To customize heatmap content after a run, use `vartracker plot heatmap results.csv [options]`.
- By default, all consequence classes are included except joint variants. Use `--include-joint` to show joint variants.
- `--aa-exclude`: comma-separated `type_of_change` patterns to exclude. Wildcards are supported.
- `--aa-include`: comma-separated `type_of_change` patterns to include.
- `--only-persistent`: only include `new_persistent` variants.
- `--only-new`: only include variants with `variant_status == new`.
- `--gene-include` and `--gene-exclude`: comma-separated gene patterns.
- `--variant-type`: comma-separated variant-type patterns such as `snp` or `indel`.
- `--qc`: comma-separated `all_samples_pass_qc` patterns to include. Accepted values include `true`, `false`, `pass`, and `fail`.
- `--min-prop-passing-qc`: minimum fraction of samples that must pass per-sample QC.
- `--min-persistence`: minimum number of included samples in which the variant must be present.
- `--min-max-af`: minimum maximum allele frequency across included samples.
- `--min-sample-af`: minimum allele frequency that must be reached in at least one included sample.
- `--sample-subset`: comma-separated sample-name patterns to plot.
- `--hide-singletons`: hide variants present in only one included sample.
- `--min-depth`: minimum site depth a variant must reach in at least one included sample.
- `--x-labels sample-number`: label heatmap x-axis columns by `sample_number` instead of sample name.
- `--title`: set the heatmap plot title. The default is `Variant allele frequencies`.
- `--literature-csv`: include literature links in the interactive HTML heatmap using a literature hits CSV.
- Example: `--aa-exclude "synonymous,*frameshift*,stop_gained"`

Standalone plot filtering:
- `--gene`, `--effect`, `--min-af`, `--max-af`: restrict the plotted result set before ranking/selection.
- `--variants` or `--variant-file`: explicitly choose variants and preserve that order.
- `--sample-min`, `--sample-max`: restrict the passage/sample-number window.
- `--persistent-only` and `--new-only`: keep only persistent new variants or only variants with `variant_status == new`.
- `trajectory` and `lifespan` auto-select a limited subset by default (`--top-n`) to stay readable.
- `turnover` uses all filtered variants by default and is also written automatically during the main `vcf`/`bam`/`end-to-end` workflows as `variant_turnover_plot.pdf`.
- `genome` uses SNPs only by default, keeps all observed allele-frequency values for each plotted variant, and writes `variant_genome_plot.pdf` during the main workflows.

Standalone plot output:
- `--out`: write to an exact file path.
- `--outdir`: write beside the results CSV or into the chosen directory using deterministic names such as `variant_trajectory_plot.pdf` or `variant_genome_plot.pdf`.
- `--format`: choose `pdf`, `png`, or `svg`.
- `--dpi`: set raster output resolution.

Genome plot options:
- `--gene`: zoom to a single gene region.
- `--aa-scale`: with `--gene`, use amino-acid coordinates on the x-axis.
- `--cds-scale`: with `--gene`, use CDS-relative nucleotide coordinates on the x-axis.
- `--focus-coords`: highlight nucleotide or amino-acid coordinate ranges, depending on the current x-axis mode. Separate color groups with `;`, ranges within a group with `,`, and optionally prefix a group with `Name:`.
- `--focus-region-file`: read named focus region groups from a `.json`, `.csv`, or `.tsv` file for an inset legend.
- `--show-intersections`: add a compact `Region | Variant` table below the genome plot for highlighted-region hits.
- In the genome plot, undetected samples are rendered at the detection threshold rather than zero; by default this floor is `0.03`, or `--min-af` if supplied, and the dashed guide line follows that same threshold.
- `--include-indels`: opt in to plotting indels too. This may be ambiguous or hard to interpret.
- The standalone genome plot auto-discovers `reference_features.json` beside `results.csv`; workflow runs generate this sidecar automatically.

Trajectory threshold mode:
- `--thresholds`: draw horizontal AF threshold lines, e.g. `0.5,0.9`.
- `--crossing-only`: keep only variants crossing at least one supplied threshold.
- `--label-threshold-crossers`: label only threshold-crossing variants to reduce clutter.
- `--crossing-rule`: choose whether threshold equality counts (`at_or_above`) or requires a strict exceedance (`strictly_above`).

Standalone plot examples:
- `vartracker plot genome results.csv`
- `vartracker plot genome results.csv --gene F`
- `vartracker plot genome results.csv --gene F --aa-scale`
- `vartracker plot genome results.csv --gene F --cds-scale --focus-coords "184-210,586-630"`
- `vartracker plot genome results.csv --focus-coords "150-300,900-1800"`
- `vartracker plot genome results.csv --focus-coords "62-69,196-210;31-42,323-332,379-399;254-277"`
- `vartracker plot genome results.csv --gene F --aa-scale --focus-coords "Ø:62-69,196-210;I:31-42,323-332,379-399;II:254-277"`
- `vartracker plot genome results.csv --gene F --aa-scale --focus-coords "Ø:62-69,196-210;I:31-42,323-332,379-399" --show-intersections`
- `vartracker plot genome results.csv --focus-region-file fusion_regions.json`
- `vartracker plot genome results.csv --gene F --aa-scale --focus-coords "50-120,180-220"`
- `vartracker plot turnover results.csv`
- `vartracker plot trajectory results.csv --variants "S:D614G,S:E484K"`
- `vartracker plot trajectory results.csv --thresholds 0.5,0.9`
- `vartracker plot trajectory results.csv --thresholds 0.5,0.9 --crossing-only`
- `vartracker plot trajectory results.csv --thresholds 0.5,0.9 --crossing-only --label-threshold-crossers`
- `vartracker plot lifespan results.csv --top-n 20 --persistent-only`

Note:
- The standalone `plot` commands require `results.csv` files written by current vartracker versions, which now include a slash-separated `sample_number` column for stable passage ordering.

- `vartracker prepare spreadsheet` – specify `--mode` (`vcf`, `bam`, or `e2e`), `--dir` to scan, `--out` for the CSV,
  and `--dry-run` to preview without writing a file.
- `vartracker prepare reference` – build a merged FASTA/GFF3 bundle from GenBank nucleotide accessions.
  Use `--accessions` or `--accession-file`, plus `--outdir`. Optional flags: `--prefix`, `--force`,
  `--keep-intermediates`, `--skip-csq-validation`.

### Using Literature Database

To search mutations against functional databases:

1. **Set up a literature database (optional):**
```bash
parse_pokay pokay_database.csv
```
   This command automatically downloads the required literature files from the
   pokay repository into `pokay_literature/NC_045512` (override with
   `--download-dir`) and writes the processed CSV for downstream analysis.

2. **Run vartracker with literature search:**
```bash
vartracker [mode] input_data.csv --literature-csv pokay_database.csv -o results/
```
   Alternatively, pass `--search-pokay` to automatically download and search
   against the Pokay SARS-CoV-2 literature database.

### Command Line Reference

```
usage: main.py [-h] [-V] {vcf,bam,end-to-end,e2e,prepare,schema} ...

positional arguments:
  {vcf,bam,end-to-end,e2e,prepare,schema}
    vcf                 Analyse VCF inputs
    bam                 Run the BAM preprocessing workflow
    end-to-end (e2e)    Run the end-to-end workflow (Snakemake + vartracker)
    prepare             Prepare inputs and references for vartracker
    schema              Print schemas for results tables or literature CSV input

options:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
```

Use `vartracker <subcommand> --help` to inspect the full list of mode-specific arguments.

### Prepare reference from accessions

Use this workflow to build a `bcftools csq`-ready reference bundle from nucleotide accessions:

```bash
# Comma-separated accessions
# Example with influenza A segments
vartracker prepare reference \
  --accessions CY114381,CY114382,CY114383,CY114384,CY114385,CY114386,CY114387,CY114388 \
  --outdir refs/influenza_a \
  --prefix influenza_a_ref

# One accession per line in a file
vartracker prepare reference \
  --accession-file accessions.txt \
  --outdir refs/
```

Required external tools:

- `bcftools` for csq smoke validation

Outputs:

- `<outdir>/<prefix>.fa`
- `<outdir>/<prefix>.gff3`
- `<outdir>/<prefix>.fa.fai`
- `<outdir>/prepare_metadata.json`

Validation notes:

- Unless `--skip-csq-validation` is supplied, vartracker writes a dummy coding-region VCF variant
  and runs `bcftools csq` against the generated FASTA/GFF3.
- Validation fails fast if `bcftools csq` exits non-zero or if the output VCF does not contain `BCSQ`.

Troubleshooting:

- Accession fetch failures: verify accession spelling and network access to NCBI efetch.
- SeqID mismatch errors: confirm FASTA headers and GFF3 seqids match exactly.
- csq validation failure: inspect the stderr snippet in the error output and confirm `bcftools`
  version and annotation structure.

### Installation Test

After installation you can verify the workflows using the bundled
demonstration dataset:

```bash
vartracker vcf --test --outdir vartracker_vcf_test_results
vartracker bam --test --outdir vartracker_bam_test_results
vartracker end-to-end --test --outdir vartracker_e2e_test_results
```

Each command copies the example dataset, resolves relative paths, checks for
the required external tools, and writes a self-contained set of results.

## Output

vartracker produces several output files:

- **results.csv**: Comprehensive variant analysis with all metrics
- **results_metadata.json**: Output schema version and results metadata
- **new_mutations.csv**: Mutations not present in the first sample
- **persistent_new_mutations.csv**: New mutations that persist to the final sample
- **cumulative_mutations.pdf**: Plot showing mutation accumulation over time
- **mutations_per_gene.pdf**: Gene-wise mutation statistics
- **variant_allele_frequency_heatmap.html**: Interactive heatmap with optional literature annotations
- **variant_allele_frequency_heatmap.pdf**: Heatmap of variant allele frequencies across passages
- **literature_database_hits.*.csv**: Functional annotation results (if literature search used)
- **run_metadata.json**: Provenance manifest capturing inputs, tool versions, and run status

By default the manifest is lightweight. Use `--manifest-level deep` to checksum all referenced
input files (FASTQ/BAM/VCF/coverage) and include file sizes.

### Output schema

The results table schema is documented in `docs/OUTPUT_SCHEMA.md`. You can also print it from the CLI:

```bash
vartracker schema results
```

To write the schema to a file instead, use:

```bash
vartracker schema results --out docs/output_schema.csv
vartracker schema results --out docs/output_schema.json --format json
```

To print the expected literature CSV structure for `--literature-csv`, use:

```bash
vartracker schema literature
```

## What does vartracker do?

The pipeline performs the following analysis:

1. **VCF Standardization**: Normalizes and standardizes input VCF files
2. **Annotation**: Adds amino acid consequences using `bcftools csq`
3. **Variant Merging**: Combines all longitudinal samples
4. **Comprehensive Analysis**: For each variant, determines:
   - Gene location and amino acid consequences
   - Variant type (SNP/indel) and change type (synonymous/missense/etc.)
   - Persistence across samples (new/original, persistent/transient)
   - Quality control metrics
   - Amino acid property changes
   - Allele frequency dynamics

5. **Visualization**: Generates plots for mutation accumulation and gene-wise statistics
6. **Functional Annotation**: (optional) Searches against literature databases for known functional impacts

## Citation

When using vartracker, please cite the software release you used. Citation metadata is provided
in `CITATION.cff`, and GitHub releases are archived on Zenodo (DOI will appear here once minted).

If you use vartracker, please cite the software record on Zenodo:

- Foster, C. (2026). *vartracker* (Version x.y.z). Zenodo. https://doi.org/10.5281/zenodo.XXXXX
- Concept DOI (all versions): https://doi.org/10.5281/zenodo.18452274
Note: a version-specific DOI is minted by Zenodo after each GitHub release.

Also cite relevant methods or data sources, for example:

- Foster CSP, et al. Long-term serial passaging of SARS-CoV-2 reveals signatures of convergent evolution. Journal of Virology. 2025;99: e00363-25. doi:10.1128/jvi.00363-25
- Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10. doi:10.1093/gigascience/giab008
- Danecek P, McCarthy SA. BCFtools/csq: haplotype-aware variant consequences. Bioinformatics. 2017;33: 2037–2039. doi:10.1093/bioinformatics/btx100
- Wilm A, Aw PPK, Bertrand D, Yeo GHT, Ong SH, Wong CH, et al. LoFreq: a sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. Nucleic Acids Res. 2012;40: 11189–11201. doi:10.1093/nar/gks918
- Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34: i884–i890. doi:10.1093/bioinformatics/bty560
- Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv. 2013 [cited 13 Apr 2021]. Available: https://arxiv.org/abs/1303.3997v2
- Mölder F, Jablonski KP, Letcher B, Hall MB, Tomkins-Tinch CH, Sochat V, et al. Sustainable data analysis with Snakemake. F1000Res. 2021;10: 33. doi:10.12688/f1000research.29032.2

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

If you encounter any issues or have questions:

1. Check the [documentation](https://github.com/charlesfoster/vartracker)
2. Search existing [issues](https://github.com/charlesfoster/vartracker/issues)
3. Create a new issue with detailed information about your problem
