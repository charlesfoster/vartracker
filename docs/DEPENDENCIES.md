# External dependencies

The table below records minimum tested versions for external tools required by vartracker.
The recommended micromamba environment pins the same versions for reproducibility.

| Tool | Minimum tested version | Pinned in environment.yml | Notes |
| --- | --- | --- | --- |
| bcftools | 1.21 | 1.21 | Required for all modes and `prepare reference` csq validation. |
| tabix | 1.21 | 1.21 (via htslib) | Provided by htslib; required for all modes. |
| samtools | 1.21 | 1.21 | Required for `bam` and `end-to-end` modes. |
| lofreq | 2.1.5 | 2.1.5 | Variant caller used in Snakemake workflow. |
| fastp | 0.23.4 | 0.23.4 | Read trimming for `end-to-end`. |
| bwa | 0.7.17 | 0.7.17 | Read alignment for `end-to-end`. |
| snakemake | 9.0.1 | 9.0.1 | Snakemake API used for workflows. |
| bgzip | 1.21 | 1.21 (via htslib) | Provided by htslib; used for VCF compression. |

If you validate compatibility with older tool versions, update the minimum tested version here and
consider expanding CI coverage.

Python package dependencies are declared in `pyproject.toml` (including `biopython` for
GenBank-to-GFF3/FASTA conversion in `prepare reference`).
