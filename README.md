[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![PyPI version](https://badge.fury.io/py/vartracker.svg)](https://badge.fury.io/py/vartracker)

```
██    ██  █████  ██████  ████████ ██████   █████   ██████ ██   ██ ███████ ██████
██    ██ ██   ██ ██   ██    ██    ██   ██ ██   ██ ██      ██  ██  ██      ██   ██
██    ██ ███████ ██████     ██    ██████  ███████ ██      █████   █████   ██████
 ██  ██  ██   ██ ██   ██    ██    ██   ██ ██   ██ ██      ██  ██  ██      ██   ██
  ████   ██   ██ ██   ██    ██    ██   ██ ██   ██  ██████ ██   ██ ███████ ██   ██


```
# vartracker

A bioinformatics pipeline to summarise variants called against a reference in a longitudinal study design. Written to investigate longitudinal sequencing data from long-term passaging of SARS-CoV-2. In theory, it could be expanded for other organisms too.

**Author:** Dr Charles Foster

## Features

- Track mutation persistence across longitudinal samples
- Comprehensive variant analysis including amino acid consequences
- Built-in SARS-CoV-2 reference data and annotations
- Integration with functional mutation databases (pokay)
- Automated plotting and statistical analysis
- Support for both SNPs and indels
- Quality control metrics for variants

## Installation

Requires Python 3.11 or newer.

### From PyPI (Recommended)

```bash
pip install vartracker
```

### External Dependencies

vartracker shells out to a handful of bioinformatics tools. Make sure they are discoverable on `PATH` before running the CLI.

- **bcftools** (>=1.10) and **tabix** – required for all modes
- **samtools**, **lofreq**, **fastp**, and **bwa** – required for the `bam` and `end-to-end` Snakemake workflows

If you only plan to run `vartracker vcf` against pre-generated VCFs, the first pair is sufficient. The additional tools are needed whenever you ask vartracker to align reads or call variants for you.

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

### Docker

Build a container image that bundles Python, vartracker, and all external bioinformatics tools:

```bash
docker build -t vartracker:latest .
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

# Generate a template spreadsheet for a directory of files
vartracker generate --mode e2e --dir data/passaging --out inputs.csv

# Exercise the bundled smoke-test dataset
vartracker vcf --test
vartracker bam --test
vartracker end-to-end --test
```

All modes understand `--test`, which copies the example dataset from `vartracker/test_data`
into a temporary directory, resolves relative paths, and runs the appropriate workflow.

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

Relative paths are resolved with respect to the CSV location, so you can store the sheet alongside
your sequencing artefacts. The `generate` subcommand can scaffold a CSV and highlight missing files.

Coverage files can be produced with `samtools depth -aa sample.bam > sample_depth.txt` or
`bedtools genomecov -ibam sample.bam -d`. The file name suffix does not matter; vartracker checks
for both `.depth.txt` and `_depth.txt` patterns when preparing its internal test dataset.

### Mode-specific options

- `vartracker vcf` – accepts plotting and filtering options such as `--min-snv-freq`, `--min-indel-freq`,
  `--allele-frequency-tag`, `--name`, `--outdir`, `--passage-cap`, and pokay controls (`--search-pokay`,
  `--pokay-csv`, `--download-pokay`). Use `--test` to run the bundled smoke test.
- `vartracker bam` – everything from `vcf`, plus Snakemake options:
  `--snakemake-outdir`, `--cores`, `--snakemake-dryrun`, `--verbose`, `--redo`.
- `vartracker end-to-end` – similar to `bam`, with an optional `--primer-bed` for amplicon clipping.
- `vartracker generate` – specify `--mode` (`vcf`, `bam`, or `e2e`), `--dir` to scan, `--out` for the CSV,
  and `--dry-run` to preview without writing a file.

### Using with pokay Database

To search mutations against functional databases:

1. **Set up pokay database (optional):**
```bash
parse_pokay pokay_database.csv
```
   This command automatically downloads the required literature files from the
   pokay repository into `pokay_literature/NC_045512` (override with
   `--download-dir`) and writes the processed CSV for downstream analysis. You
   can also let `vartracker` download and parse the database on demand with the
   `--download-pokay` flag.

2. **Run vartracker with pokay search:**
```bash
vartracker input_data.csv --search-pokay --pokay-csv pokay_database.csv -o results/
```
   Alternatively, omit `--pokay-csv` and pass `--download-pokay` to fetch the
   database automatically during execution.

### Command Line Reference

```
usage: main.py [-h] [-V] {vcf,bam,end-to-end,e2e,generate} ...

positional arguments:
  {vcf,bam,end-to-end,e2e,generate}
    vcf                 Analyse VCF inputs
    bam                 Run the BAM preprocessing workflow
    end-to-end (e2e)    Run the end-to-end workflow (Snakemake + vartracker)
    generate            Generate input spreadsheets from an existing directory of files

options:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
```

Use `vartracker <subcommand> --help` to inspect the full list of mode-specific arguments.

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
- **new_mutations.csv**: Mutations not present in the first sample
- **persistent_new_mutations.csv**: New mutations that persist to the final sample
- **cumulative_mutations.pdf**: Plot showing mutation accumulation over time
- **mutations_per_gene.pdf**: Gene-wise mutation statistics
- **variant_allele_frequency_heatmap.pdf**: Heatmap of variant allele frequencies across passages
- **pokay_database_hits.*.csv**: Functional annotation results (if pokay used)

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

When using vartracker, please cite:

- The `vartracker` tool:
  - Foster, C. (2025). vartracker: Track the persistence (or not) of mutations during longitudinal sequencing (Version 2.0). GitHub. https://github.com/charlesfoster/vartracker
  - Foster CSP, et al. Long-term serial passaging of SARS-CoV-2 reveals signatures of convergent evolution. Journal of Virology. 2025;99: e00363-25. doi:10.1128/jvi.00363-25
- bcftools: Li H. (2011) A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics, 27, 2987-2993.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

If you encounter any issues or have questions:

1. Check the [documentation](https://github.com/charlesfoster/vartracker)
2. Search existing [issues](https://github.com/charlesfoster/vartracker/issues)
3. Create a new issue with detailed information about your problem
