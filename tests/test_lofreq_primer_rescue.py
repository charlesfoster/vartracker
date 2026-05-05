import subprocess
from pathlib import Path

from vartracker.lofreq_primer_rescue import rescue_lofreq_primer_variants


def test_rescue_lofreq_primer_variants_keeps_pass_and_rescues_overlap(
    monkeypatch, tmp_path
):
    raw_vcf = tmp_path / "raw.vcf"
    raw_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##INFO=<ID=DP4,Number=4,Type=Integer,Description="Strand depths">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tG\t200\t.\tAF=0.98;DP=120;DP4=1,0,118,0\n"
        "chr1\t30\t.\tC\tT\t150\t.\tAF=0.40;DP=120;DP4=30,30,30,30\n",
        encoding="utf-8",
    )
    primers = tmp_path / "primers.bed"
    primers.write_text("chr1\t0\t20\tprimer_1\n", encoding="utf-8")
    default_filtered_vcf = (
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##INFO=<ID=DP4,Number=4,Type=Integer,Description="Strand depths">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tG\t200\tstrandbias\tAF=0.98;DP=120;DP4=1,0,118,0\n"
        "chr1\t30\t.\tC\tT\t150\tPASS\tAF=0.40;DP=120;DP4=30,30,30,30\n"
    )

    def fake_run(cmd, check):
        assert cmd[:4] == ["lofreq", "filter", "-i", str(raw_vcf)]
        assert check is True
        output_path = Path(cmd[cmd.index("-o") + 1])
        output_path.write_text(default_filtered_vcf, encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr("vartracker.lofreq_primer_rescue.subprocess.run", fake_run)

    output_vcf = tmp_path / "final.vcf"
    rescued_tsv = tmp_path / "rescued.tsv"

    result = rescue_lofreq_primer_variants(
        raw_vcf=raw_vcf,
        primers_bed=primers,
        output_vcf=output_vcf,
        rescued_tsv=rescued_tsv,
    )

    assert result.normal_passed == 1
    assert result.rescued == 1
    assert result.discarded == 0

    output_text = output_vcf.read_text(encoding="utf-8")
    assert "##FILTER=<ID=RESCUED_PRIMER_OVERLAP," in output_text
    assert "chr1\t10\t.\tA\tG\t200\tRESCUED_PRIMER_OVERLAP" in output_text
    assert "PRIMER_OVERLAP;RESCUED_BY=overlap_primer_interval" in output_text
    assert "chr1\t30\t.\tC\tT\t150\tPASS" in output_text

    rescued_lines = rescued_tsv.read_text(encoding="utf-8").splitlines()
    assert rescued_lines[0] == "variant\treason_filtered\treason_rescued\tmetrics"
    assert rescued_lines[1].startswith(
        "A10G\tstrandbias\toverlap_primer_interval\tAF=0.98"
    )
