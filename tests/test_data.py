"""Tests for vartracker.data module."""

import os

from vartracker.data import get_data_file, get_reference_genome, get_annotation


class TestDataModule:
    """Test data module functions."""

    def test_get_data_file(self):
        """Test get_data_file function."""
        filename = "NC_045512.fasta"
        result = get_data_file(filename)

        assert isinstance(result, str)
        assert filename in result
        assert os.path.exists(result)

    def test_get_reference_genome(self):
        """Test get_reference_genome function."""
        result = get_reference_genome()

        assert isinstance(result, str)
        assert "NC_045512.fasta" in result
        assert os.path.exists(result)

        # Check that it's a valid FASTA file
        with open(result, "r") as f:
            first_line = f.readline().strip()
            assert first_line.startswith(">")

    def test_get_annotation(self):
        """Test get_annotation function."""
        result = get_annotation()

        assert isinstance(result, str)
        assert "NC_045512.gff3" in result
        assert os.path.exists(result)

        # Check that it's a valid GFF3 file
        with open(result, "r") as f:
            content = f.read()
            assert len(content) > 0

    def test_all_data_files_exist(self):
        """Test that all expected data files exist."""
        expected_files = ["NC_045512.fasta", "NC_045512.gff3", "NC_045512.fasta.fai"]

        for filename in expected_files:
            file_path = get_data_file(filename)
            assert os.path.exists(file_path), f"Data file {filename} does not exist"
