from pathlib import Path


def test_snakemake_rules_write_logs_under_outdir():
    snakefile = (
        Path(__file__).resolve().parents[1] / "vartracker" / "Snakefile"
    ).read_text(encoding="utf-8")

    assert "/dev/null" not in snakefile
    assert 'return f"{OUTDIR}/logs/{rule_name}.log"' in snakefile
    assert 'return f"{OUTDIR}/{sample}/logs/{rule_name}.log"' in snakefile

    expected_logs = [
        'log:\n            _workflow_log("bwa_index")',
        'lambda w: _sample_log(w.sample, "fastp")',
        'lambda w: _sample_log(w.sample, "bwa_mem")',
        'lambda w: _sample_log(w.sample, "ampliconclip")',
        'lambda w: _sample_log(w.sample, "lofreq_indelqual")',
        'lambda w: _sample_log(w.sample, "samtools_depth")',
        'lambda w: _sample_log(w.sample, "lofreq_call")',
        'log:\n        _workflow_log("update_csv")',
    ]

    for expected in expected_logs:
        assert expected in snakefile
