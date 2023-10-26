[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
```
██    ██  █████  ██████  ████████ ██████   █████   ██████ ██   ██ ███████ ██████  
██    ██ ██   ██ ██   ██    ██    ██   ██ ██   ██ ██      ██  ██  ██      ██   ██ 
██    ██ ███████ ██████     ██    ██████  ███████ ██      █████   █████   ██████  
 ██  ██  ██   ██ ██   ██    ██    ██   ██ ██   ██ ██      ██  ██  ██      ██   ██ 
  ████   ██   ██ ██   ██    ██    ██   ██ ██   ██  ██████ ██   ██ ███████ ██   ██ 
                                                                                  
                                                                                  
```
# vartracker
A simple pipeline to summarise variants called against a reference in a longitudinal study design. Written to investigate longitudinal sequencing data from long-term passaging of SARS-CoV-2.

Author: Dr Charles Foster

# Starting out
To begin with, clone this github repository:

```
git clone https://github.com/charlesfoster/vartracker.git

cd vartracker

```

Next, install dependencies using [mamba](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install):

```
mamba env create -f environment.yml
```

# Before the first run
One useful part of the pipeline is that it searches de novo mutations (i.e., those not in the initial sample) against a compendium of functional SARS-CoV-2 literature.  If this feature is not wanted, the pipeline can be run with the `--skip_search` argument, otherwise before the first run, this compendium must be initialised locally. Steps:

Clone the repo:
```
git clone https://github.com/nodrogluap/pokay.git
```

Get the path to the repo's `data/literature/NC_045512` directory. For example, if saved in `~/Programs`, the path will be `~/Programs/pokay/data/literature/NC_045512`.

Next, run the `parse_pokay.py` script within the `vartracker` repo: `vartracker/util/parse_pokay.py`

```
conda activate vartracker
cd ~/Programs/vartracker
python util/parse_pokay.py  ~/Programs/pokay/data/literature/NC_045512 ~/Programs/vartracker/pokay_database.csv
```

Now, when running `vartracker.py`, the database should be provided to the script (`--pokay ~/Programs/vartracker/pokay_database.csv`) (see below)

# Usage
Activate the environment:

```
conda activate vartracker
```

See usage:

```
python vartracker.py --help
```

Output:

```
██    ██  █████  ██████  ████████ ██████   █████   ██████ ██   ██ ███████ ██████  
██    ██ ██   ██ ██   ██    ██    ██   ██ ██   ██ ██      ██  ██  ██      ██   ██ 
██    ██ ███████ ██████     ██    ██████  ███████ ██      █████   █████   ██████  
 ██  ██  ██   ██ ██   ██    ██    ██   ██ ██   ██ ██      ██  ██  ██      ██   ██ 
  ████   ██   ██ ██   ██    ██    ██   ██ ██   ██  ██████ ██   ██ ███████ ██   ██ 
                                                                                  
                                                                                  

usage: python3 vartracker.py [options] <input_csv> 

vartracker: track the persistence (or not) or mutations during long-term passaging

positional arguments:
  input_csv             Input CSV file. See below.

options:
  -h, --help            show this help message and exit
  -a ANNOTATION, --annotation ANNOTATION
                        Annotations to use in GFF3 format
  -c CHROM_MAP, --chrom_map CHROM_MAP
                        [DEPRECATED] File with map for renaming chrom in VCF files
  -m MIN_SNV_FREQ, --min_snv_freq MIN_SNV_FREQ
                        Minimum allele frequency of SNV variants to keep
  -M MIN_INDEL_FREQ, --min_indel_freq MIN_INDEL_FREQ
                        Minimum allele frequency of indel variants to keep
  -n NAME, --name NAME  Optional: add a column to results with the name specified here
  -o OUTDIR, --outdir OUTDIR
                        Output directory
  -p POKAY, --pokay POKAY
                        Path to csv file of the parsed 'pokay' database
  -f FILENAME, --filename FILENAME
                        Output file name
  -r REFERENCE, --reference REFERENCE
                        Reference genome
  --passage_cap PASSAGE_CAP
                        Cap the number of passages at this number
  --debug               Print commands being run for debugging
  --keep_temp           Keep temporary files for debugging
  --rename_chrom        [DEPRECATED] Rename chromosomes in VCF files. Useful if different reference
                        names were given during VCF calling. Uses 'chrom_map' file. See
                        bcftools/samtools usage guide.
  --skip_search         Skip literature search for 'new mutations'.

The input CSV file should have four columns:
1. 'vcf': full paths to bgzipped vcf files containing the called variants for each sample of interest
2. 'coverage': full paths to files containing the genome coverage for each sample. The coverage files
   should be in the format output by `bedtools genomecov` (e.g., bedtools genomecov -ibam input.bam -d
   > coverage.txt), whereby the columns are (in order): the reference name, the 1-based position in the
   genome, and the depth of coverage at that position.
3. 'sample_name': the name of the sample in the given VCF file
4. 'sample_number': the sample number. In a longitudinal comparison like long-term passaging, the
   numbers in the column might go 0, 1, 2, ..., 15.
 
Note: the order of values in the input CSV file matters, dictating the order of results in the output
CSV file.
 
If 'new mutations' are to be searched against the functional 'pokay' database (default) (see vartracker
repository README), the path to the parsed pokay csv file must be provided.
```

# What does the pipeline do?
All longitudinal VCF files are standardised, amino acid consequences are called, and then are merged. All variants are then interrogated by custom code and written to file. The output file contains a comprehensive summary and analysis of all variants, including:
- which genes they occur in
- different notations for the amino acid changes
- the type of variant (SNP or indel)
- the type of change (synonymous, missense, frameshift etc. [as classified by `bcftools csq`])
- the variant status (new or original)
- whether the variant persists until the end of the longitudinal study (i.e. is present in passage/sample N)
- various QC metrics for the variants and a PASS/FAIL status
- amino acid properties and how they change (charged vs uncharged etc.)
- changes in allele frequency
- changes in total genome coverage
- changes in variant depth
- and more

Additionally, hits against the functional 'pokay' database (see above) are written to file, and plots of mutations per gene and cumulative mutations are generated.

# Credits
* When this pipeline is used, citations should be found for the programs used internally (minimally `bcftools`).
* A link to this repo should also be cited