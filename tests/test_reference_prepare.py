"""Tests for reference preparation utilities."""

from __future__ import annotations

import io
import subprocess
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from vartracker import reference_prepare as prep


def _make_genbank_text(seqid: str, gene_name: str = "GENE1") -> str:
    record = SeqRecord(
        Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"),
        id=seqid,
        name=seqid,
        description="synthetic test sequence",
    )
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(
            FeatureLocation(0, 30, strand=1),
            type="CDS",
            qualifiers={
                "gene": [gene_name],
                "codon_start": ["1"],
                "protein_id": [f"{seqid}.P1"],
            },
        )
    ]
    buffer = io.StringIO()
    SeqIO.write(record, buffer, "genbank")
    return buffer.getvalue()


def test_parse_accessions_from_string():
    values = prep.parse_accessions(accessions="CY114381, CY114382, CY114381")
    assert values == ["CY114381", "CY114382"]


def test_parse_accessions_from_file(tmp_path):
    path = tmp_path / "accessions.txt"
    path.write_text("CY114381\n# comment\nCY114382\n", encoding="utf-8")
    values = prep.parse_accessions(accession_file=path)
    assert values == ["CY114381", "CY114382"]


def test_prepare_reference_bundle_merges_outputs(tmp_path):
    data = {
        "ACC1": _make_genbank_text("ACC1", "HA"),
        "ACC2": _make_genbank_text("ACC2", "NA"),
    }

    def fetcher(accession: str) -> str:
        return data[accession]

    meta = prep.prepare_reference_bundle(
        accessions=["ACC1", "ACC2"],
        outdir=tmp_path,
        prefix="flu_ref",
        keep_intermediates=True,
        skip_csq_validation=True,
        fetcher=fetcher,
    )

    fasta_path = Path(meta["outputs"]["fasta"])
    gff_path = Path(meta["outputs"]["gff3"])
    lines = gff_path.read_text(encoding="utf-8").splitlines()
    assert sum(1 for x in lines if x.startswith("##gff-version 3")) == 1
    seqids = {
        line.split("\t")[0] for line in lines if line and not line.startswith("#")
    }
    assert seqids == {"ACC1", "ACC2"}
    fasta_ids = [record.id for record in SeqIO.parse(str(fasta_path), "fasta")]
    assert fasta_ids == ["ACC1", "ACC2"]


def test_validate_csq_with_dummy_variant_success(tmp_path, monkeypatch):
    fasta = tmp_path / "ref.fa"
    gff = tmp_path / "ref.gff3"
    fasta.write_text(">ACC1\nATGATGATGATG\n", encoding="utf-8")
    gff.write_text(
        "##gff-version 3\n"
        "ACC1\tvartracker_prepare\tgene\t1\t9\t.\t+\t.\tID=gene:ACC1_1;Name=GENE1;biotype=protein_coding\n"
        "ACC1\tvartracker_prepare\tmRNA\t1\t9\t.\t+\t.\tID=transcript:ACC1_1;Parent=gene:ACC1_1;Name=GENE1;biotype=protein_coding\n"
        "ACC1\tvartracker_prepare\tCDS\t1\t9\t.\t+\t0\tID=cds:ACC1_1;Parent=transcript:ACC1_1;gene=GENE1\n",
        encoding="utf-8",
    )

    def fake_run(cmd, capture_output, text, check):
        out_idx = cmd.index("-o") + 1
        output = Path(cmd[out_idx])
        output.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "ACC1\t1\t.\tA\tC\t.\tPASS\tBCSQ=missense|GENE1|tx|protein_coding|+|1M>1L|1A>C\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(prep.subprocess, "run", fake_run)
    result = prep.validate_csq_with_dummy_variant(
        fasta_path=fasta, gff_path=gff, workdir=tmp_path
    )
    assert result["status"] == "passed"

    dummy_vcf = Path(result["dummy_vcf"])
    content = dummy_vcf.read_text(encoding="utf-8")
    assert "ACC1\t1\t.\tA\tC" in content


def test_validate_csq_with_dummy_variant_failure(tmp_path, monkeypatch):
    fasta = tmp_path / "ref.fa"
    gff = tmp_path / "ref.gff3"
    fasta.write_text(">ACC1\nATGATGATGATG\n", encoding="utf-8")
    gff.write_text(
        "##gff-version 3\n"
        "ACC1\tvartracker_prepare\tCDS\t1\t9\t.\t+\t0\tID=cds:ACC1_1;Parent=transcript:ACC1_1;gene=GENE1\n",
        encoding="utf-8",
    )

    def fake_run(cmd, capture_output, text, check):
        return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="csq failed")

    monkeypatch.setattr(prep.subprocess, "run", fake_run)
    with pytest.raises(RuntimeError, match="csq validation failed"):
        prep.validate_csq_with_dummy_variant(
            fasta_path=fasta,
            gff_path=gff,
            workdir=tmp_path,
        )


def test_prepare_reference_end_to_end_with_mocked_download(tmp_path, monkeypatch):
    gb_map = {
        "CY114381": _make_genbank_text("CY114381", "PB2"),
        "CY114382": _make_genbank_text("CY114382", "PB1"),
    }

    def fake_fetch(accession: str) -> str:
        return gb_map[accession]

    monkeypatch.setattr(
        prep,
        "validate_csq_with_dummy_variant",
        lambda **kwargs: {"status": "passed", "command": "bcftools csq ..."},
    )

    meta = prep.prepare_reference_bundle(
        accessions=["CY114381", "CY114382"],
        outdir=tmp_path / "refs",
        prefix="flu",
        keep_intermediates=False,
        skip_csq_validation=False,
        fetcher=fake_fetch,
    )

    outputs = meta["outputs"]
    assert Path(outputs["fasta"]).exists()
    assert Path(outputs["gff3"]).exists()
    assert Path(outputs["fai"]).exists()
    assert Path(outputs["metadata"]).exists()
