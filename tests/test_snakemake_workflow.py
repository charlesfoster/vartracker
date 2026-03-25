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
        'f"{OUTDIR}/logs/update_csv.log"',
    ]

    for expected in expected_logs:
        assert expected in snakefile
