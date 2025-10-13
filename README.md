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

### From PyPI (Recommended)

```bash
pip install vartracker
```

### External Dependencies

vartracker requires the following external tools to be installed and available in your PATH:

- **bcftools** (>=1.10): For VCF file processing and variant calling consequences
- **tabix**: For indexing compressed files

#### Installing bcftools and tabix

**On macOS:**
```bash
# Using Homebrew
brew install bcftools htslib

# Using MacPorts
sudo port install bcftools htslib
```

**On Linux (Ubuntu/Debian):**
```bash
sudo apt-get update
sudo apt-get install bcftools tabix
```

**On Linux (CentOS/RHEL/Fedora):**
```bash
# CentOS/RHEL with EPEL
sudo yum install epel-release
sudo yum install bcftools htslib

# Fedora
sudo dnf install bcftools htslib
```

**Using conda:**
```bash
conda install -c bioconda bcftools htslib
```

### Development Installation

For development or to get the latest version:

```bash
git clone https://github.com/charlesfoster/vartracker.git
cd vartracker
pip install -e .[dev]
pre-commit install
```

## Quick Start

After installation, `vartracker` will be available as a command-line tool:

```bash
vartracker --help
```

### Basic Usage

```bash
vartracker input_data.csv -o output_directory
```

### Input Format

The input CSV file should have four columns:

1. **vcf**: Full paths to VCF files (bgzipped or uncompressed) containing variants for each sample called against the NC_045512.2 reference genome. The VCF must have depth ("DP") and variant allele frequency tags in the INFO field. The variant allele frequency tag name can be specified via the `vartracker` CLI (see below).
2. **coverage**: Full paths to coverage files. These files are expected to be in TSV format
with three columns (no header): `reference\t1-based position\tdepth`. For example:

```
NC_045512.2	1	0
NC_045512.2	2	0
NC_045512.2	3	0
NC_045512.2	4	0
NC_045512.2	5	0
NC_045512.2	6	0
NC_045512.2	7	0
NC_045512.2	8	0
NC_045512.2	9	0
NC_045512.2	10	0
...
NC_045512.2	250	932
NC_045512.2	251	929
NC_045512.2	252	931
NC_045512.2	253	900
NC_045512.2	254	900
NC_045512.2	255	897
NC_045512.2	256	891
NC_045512.2	257	895
NC_045512.2	258	903
NC_045512.2	259	898
NC_045512.2	260	895
```

Note: this is the format produced by `bedtools genomecov` via the command `bedtools genomecov -ibam input.bam -d > coverage.txt` or `samtools depth` via the command `samtools depth -a input.bam`.

3. **sample_name**: Name of the sample in the VCF file
4. **sample_number**: Sample number (e.g., 0, 1, 2, ..., 15 for passages)
5. **reads1** (optional in VCF mode): Path to read 1 FASTQ (leave blank if not available)
6. **reads2** (optional): Path to read 2 FASTQ (leave blank for single-end data)
7. **bam** (optional in VCF/e2e modes): Pre-processed BAM path if already generated
8. **vcf**: Path to the sample VCF file
9. **coverage**: Path to the per-base coverage file

**Example input CSV:**
```csv
sample_name,sample_number,reads1,reads2,bam,vcf,coverage
Passage_0,0,,,,/path/to/sample0.vcf.gz,/path/to/sample0.coverage
Passage_1,1,,,,/path/to/sample1.vcf.gz,/path/to/sample1.coverage
Passage_2,2,,,,/path/to/sample2.vcf.gz,/path/to/sample2.coverage
```

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

### Command Line Options

```
usage: vartracker [options] <input_csv>

vartracker: track the persistence (or not) of mutations during long-term passaging

positional arguments:
  input_csv             Input CSV file. See below.

options:
  -h, --help            show this help message and exit
  -g, --gff3 GFF3       GFF3 annotations to use (default: packaged SARS-CoV-2 annotations)
  -m, --min-snv-freq MIN_SNV_FREQ
                        Minimum allele frequency of SNV variants to keep (default: 0.03)
  -M, --min-indel-freq MIN_INDEL_FREQ
                        Minimum allele frequency of indel variants to keep (default: 0.1)
  -n, --name NAME       Optional: add a column to results with the name specified here
  -o, --outdir OUTDIR   Output directory (default: current directory)
  --pokay-csv POKAY_CSV
                        Path to a pre-parsed pokay database CSV file
  --search-pokay        Run literature lookups against the pokay database
  --download-pokay      Automatically download and parse the pokay database
  --test                Run vartracker against the bundled demonstration dataset
  -f, --filename FILENAME
                        Output file name (default: results.csv)
  -r, --reference REFERENCE
                        Reference genome (default: uses packaged SARS-CoV-2 reference)
  --passage-cap PASSAGE_CAP
                        Cap the number of passages at this number
  --debug               Print commands being run for debugging
  --keep-temp           Keep temporary files for debugging
  --allele-frequency-tag ALLELE_FREQUENCY_TAG
                        INFO tag name for allele frequency (default: AF)
  -V, --version         show program's version number and exit
```

### Installation Test

After installation you can verify the full pipeline using the bundled
demonstration dataset:

```bash
vartracker --test --outdir vartracker_test_results
```

This command runs vartracker end-to-end with packaged inputs, producing outputs
in the specified directory and confirming that external dependencies (such as
bcftools) are available.

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

- The vartracker tool: `https://github.com/charlesfoster/vartracker`
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
