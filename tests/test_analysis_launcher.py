import pytest

from vartracker.analysis_launcher import (
    _validate_lofreq_primer_rescue,
    _validate_primer_bed_reference,
)


def test_validate_primer_bed_reference_accepts_matching_contig(tmp_path):
    reference = tmp_path / "reference.fasta"
    reference.write_text(">NC_045512.2 reference\nACGT\n", encoding="utf-8")
    primer_bed = tmp_path / "primers.bed"
    primer_bed.write_text(
        "track name=primers\n"
        "NC_045512.2\t0\t20\tprimer_1\n"
        "NC_045512.2\t20\t40\tprimer_2\n",
        encoding="utf-8",
    )

    _validate_primer_bed_reference(primer_bed, reference)


def test_validate_primer_bed_reference_rejects_mismatched_contig(tmp_path):
    reference = tmp_path / "reference.fasta"
    reference.write_text(">NC_045512.2 reference\nACGT\n", encoding="utf-8")
    primer_bed = tmp_path / "primers.bed"
    primer_bed.write_text("MN908947.3\t0\t20\tprimer_1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="do not match reference FASTA ID"):
        _validate_primer_bed_reference(primer_bed, reference)


def test_validate_primer_bed_reference_rejects_empty_bed(tmp_path):
    reference = tmp_path / "reference.fasta"
    reference.write_text(">NC_045512.2 reference\nACGT\n", encoding="utf-8")
    primer_bed = tmp_path / "primers.bed"
    primer_bed.write_text("# no intervals\n", encoding="utf-8")

    with pytest.raises(ValueError, match="contains no intervals"):
        _validate_primer_bed_reference(primer_bed, reference)


def test_validate_lofreq_primer_rescue_on_requires_primer_bed():
    with pytest.raises(ValueError, match="requires --primer-bed"):
        _validate_lofreq_primer_rescue(
            "on",
            None,
            min_af=0.95,
            min_dp=100,
            min_alt_count=95,
            min_qual=100,
            max_ref_count=20,
        )


def test_validate_lofreq_primer_rescue_rejects_invalid_threshold():
    with pytest.raises(ValueError, match="between 0 and 1"):
        _validate_lofreq_primer_rescue(
            "auto",
            None,
            min_af=1.5,
            min_dp=100,
            min_alt_count=95,
            min_qual=100,
            max_ref_count=20,
        )
