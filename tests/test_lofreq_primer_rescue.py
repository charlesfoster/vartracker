import subprocess
from pathlib import Path

from vartracker.lofreq_primer_rescue import (
    lofreq_filter_with_audit,
    rescue_lofreq_primer_variants,
)


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
        "chr1\t15\t.\tT\tC\t210\t.\tAF=0.98;DP=125;DP4=0,5,114,5\n"
        "chr1\t18\t.\tG\tT\t220\t.\tAF=0.98;DP=125;DP4=0,4,115,5\n"
        "chr1\t30\t.\tC\tT\t150\t.\tAF=0.40;DP=120;DP4=30,30,30,30\n"
        "chr1\t50\t.\tT\tC\t160\t.\tAF=0.96;DP=130;DP4=2,2,60,65\n",
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
        "chr1\t15\t.\tT\tC\t210\tsb_fdr\tAF=0.98;DP=125;DP4=0,5,114,5\n"
        "chr1\t18\t.\tG\tT\t220\tlowqual\tAF=0.98;DP=125;DP4=0,4,115,5\n"
        "chr1\t30\t.\tC\tT\t150\tPASS\tAF=0.40;DP=120;DP4=30,30,30,30\n"
        "chr1\t50\t.\tT\tC\t160\tsb_fdr\tAF=0.96;DP=130;DP4=2,2,60,65\n"
    )

    def fake_run(cmd, check):
        assert cmd[:5] == ["lofreq", "filter", "--print-all", "-i", str(raw_vcf)]
        assert check is True
        output_path = Path(cmd[cmd.index("-o") + 1])
        output_path.write_text(default_filtered_vcf, encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr("vartracker.lofreq_primer_rescue.subprocess.run", fake_run)

    output_vcf = tmp_path / "final.vcf"
    rescued_tsv = tmp_path / "rescued.tsv"
    filtered_out_tsv = tmp_path / "filtered_out.tsv"

    result = rescue_lofreq_primer_variants(
        raw_vcf=raw_vcf,
        primers_bed=primers,
        output_vcf=output_vcf,
        rescued_tsv=rescued_tsv,
        filtered_out_tsv=filtered_out_tsv,
    )

    assert result.normal_passed == 1
    assert result.rescued == 2
    assert result.discarded == 2
    assert result.filtered_out == 2

    output_text = output_vcf.read_text(encoding="utf-8")
    assert "##FILTER=<ID=RESCUED_PRIMER_OVERLAP," in output_text
    assert "chr1\t10\t.\tA\tG\t200\tRESCUED_PRIMER_OVERLAP" in output_text
    assert "chr1\t15\t.\tT\tC\t210\tRESCUED_PRIMER_OVERLAP" in output_text
    assert "chr1\t18\t.\tG\tT" not in output_text
    assert "PRIMER_OVERLAP;RESCUED_BY=overlap_primer_interval" in output_text
    assert "chr1\t30\t.\tC\tT\t150\tPASS" in output_text

    rescued_lines = rescued_tsv.read_text(encoding="utf-8").splitlines()
    assert rescued_lines[0] == "variant\treason_filtered\treason_rescued\tmetrics"
    assert rescued_lines[1].startswith(
        "A10G\tstrandbias\toverlap_primer_interval\tAF=0.98"
    )
    assert rescued_lines[2].startswith("T15C\tsb_fdr\toverlap_primer_interval\tAF=0.98")
    assert "minor_alt_fraction=0.0420168" in rescued_lines[2]

    filtered_out_lines = filtered_out_tsv.read_text(encoding="utf-8").splitlines()
    assert filtered_out_lines[0] == "variant\treason_filtered\tmetrics"
    assert filtered_out_lines[1].startswith("G18T\tlowqual\tAF=0.98")
    assert filtered_out_lines[2].startswith("T50C\tsb_fdr\tAF=0.96")


def test_lofreq_filter_with_audit_writes_pass_vcf_and_filtered_table(
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
        "chr1\t50\t.\tT\tC\t160\t.\tAF=0.96;DP=130;DP4=2,2,60,65\n",
        encoding="utf-8",
    )
    default_filtered_vcf = (
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##INFO=<ID=DP4,Number=4,Type=Integer,Description="Strand depths">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tG\t200\tPASS\tAF=0.98;DP=120;DP4=1,0,118,0\n"
        "chr1\t50\t.\tT\tC\t160\tsb_fdr\tAF=0.96;DP=130;DP4=2,2,60,65\n"
    )

    def fake_run(cmd, check):
        assert cmd[:5] == ["lofreq", "filter", "--print-all", "-i", str(raw_vcf)]
        assert check is True
        output_path = Path(cmd[cmd.index("-o") + 1])
        output_path.write_text(default_filtered_vcf, encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr("vartracker.lofreq_primer_rescue.subprocess.run", fake_run)

    output_vcf = tmp_path / "final.vcf"
    filtered_out_tsv = tmp_path / "filtered_out.tsv"
    result = lofreq_filter_with_audit(
        raw_vcf=raw_vcf,
        output_vcf=output_vcf,
        filtered_out_tsv=filtered_out_tsv,
    )

    assert result.normal_passed == 1
    assert result.filtered_out == 1
    output_text = output_vcf.read_text(encoding="utf-8")
    assert "chr1\t10\t.\tA\tG\t200\tPASS" in output_text
    assert "chr1\t50\t.\tT\tC" not in output_text
    filtered_out_lines = filtered_out_tsv.read_text(encoding="utf-8").splitlines()
    assert filtered_out_lines[0] == "variant\treason_filtered\tmetrics"
    assert filtered_out_lines[1].startswith("T50C\tsb_fdr\tAF=0.96")
