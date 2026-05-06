from pathlib import Path


def test_snakemake_rules_write_logs_under_outdir():
    snakefile = (
        Path(__file__).resolve().parents[1] / "vartracker" / "Snakefile"
    ).read_text(encoding="utf-8")

    assert "/dev/null" not in snakefile
    assert "threads: min(8, max(1, workflow.cores))" in snakefile

    expected_logs = [
        'f"{OUTDIR}/logs/bwa_index.log"',
        'f"{OUTDIR}/{{sample}}/logs/fastp.log"',
        'f"{OUTDIR}/{{sample}}/logs/bwa_mem.log"',
        'f"{OUTDIR}/{{sample}}/logs/ampliconclip.log"',
        'f"{OUTDIR}/{{sample}}/logs/lofreq_indelqual.log"',
        'f"{OUTDIR}/{{sample}}/logs/samtools_depth.log"',
        'f"{OUTDIR}/{{sample}}/logs/lofreq_call.log"',
        'f"{OUTDIR}/{{sample}}/logs/deletion_variants_bed.log"',
        'f"{OUTDIR}/{{sample}}/logs/depth_mask.log"',
        'f"{OUTDIR}/{{sample}}/logs/iupac_genotyped_vcf.log"',
        'f"{OUTDIR}/{{sample}}/logs/bcftools_consensus.log"',
        'f"{OUTDIR}/{{sample}}/logs/bcftools_iupac_consensus.log"',
        'f"{OUTDIR}/logs/update_csv.log"',
    ]

    for expected in expected_logs:
        assert expected in snakefile

    assert "bedtools" not in snakefile
    assert "--both-ends" not in snakefile
    assert "--tolerance {params.tolerance}" in snakefile
    assert "--no-default-filter" in snakefile
    assert "_variants.raw.vcf.gz" in snakefile
    assert "_variants.raw.vcf.gz.tbi" in snakefile
    assert (
        'vcf_raw = temp(f"{OUTDIR}/{{sample}}/{{sample}}_variants.raw.vcf.gz")'
        not in snakefile
    )
    assert (
        'vcf_raw = f"{OUTDIR}/{{sample}}/{{sample}}_variants.raw.vcf.gz"' in snakefile
    )
    assert (
        'vcf_raw_tbi = f"{OUTDIR}/{{sample}}/{{sample}}_variants.raw.vcf.gz.tbi"'
        in snakefile
    )
    assert "_variants.filtered_out.tsv" in snakefile
    assert "LOFREQ_PRIMER_RESCUE_ENABLED" in snakefile
    assert "LOFREQ_RESCUE_MAX_MINOR_ALT_FRACTION" in snakefile
    assert "max_minor_alt_fraction=LOFREQ_RESCUE_MAX_MINOR_ALT_FRACTION" in snakefile
    assert "lofreq_filter_with_audit" in snakefile
    assert "rescue_lofreq_primer_variants" in snakefile
    assert "_variants.rescued.tsv" in snakefile
    assert "_validate_primer_bed_reference(PRIMER_BED, REF)" in snakefile
    assert "_consensus.fasta" in snakefile
    assert "_iupac_consensus.fasta" in snakefile
    assert "df['consensus']" in snakefile
    assert "df['iupac_consensus']" in snakefile
    assert "df['raw_vcf']" in snakefile
    assert "df['lofreq_rescued_tsv']" in snakefile
    assert "df['lofreq_filtered_out_tsv']" in snakefile
