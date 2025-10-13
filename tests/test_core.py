"""Tests for vartracker.core module."""

import os
import tempfile
import pytest
import pandas as pd
from unittest.mock import patch
import re

from vartracker.core import (
    get_package_data_path,
    get_logo,
    validate_input_file,
    check_dependencies,
    DependencyError,
    InputValidationError,
    ProcessingError,
)


class TestGetPackageDataPath:
    """Test get_package_data_path function."""

    def test_get_package_data_path_returns_string(self):
        """Test that get_package_data_path returns a string."""
        result = get_package_data_path("NC_045512.fasta")
        assert isinstance(result, str)
        assert "NC_045512.fasta" in result

    def test_get_package_data_path_file_exists(self):
        """Test that returned path points to existing file."""
        result = get_package_data_path("NC_045512.fasta")
        assert os.path.exists(result)


class TestGetLogo:
    """Test get_logo function."""

    @staticmethod
    def _looks_like_block_art(s: str) -> bool:
        # Heuristics: has block/ASCII art chars and is multi-line
        return (re.search(r"[█▓▒░#@*]{3,}", s) is not None) and (s.count("\n") >= 2)

    def test_get_logo_returns_string(self):
        result = get_logo()
        assert isinstance(result, str)
        assert result.strip()  # non-empty
        # Accept either a wordmark or block-art rendering
        assert ("VARTRACKER" in result.upper()) or self._looks_like_block_art(result)


class TestValidateInputFile:
    """Test validate_input_file function."""

    def test_valid_input_file(self):
        """Test validation with a valid input file."""
        with (
            tempfile.NamedTemporaryFile(
                mode="w", suffix=".vcf", delete=False
            ) as vcf_temp,
            tempfile.NamedTemporaryFile(
                mode="w", suffix=".txt", delete=False
            ) as cov_temp,
            tempfile.NamedTemporaryFile(
                mode="w", suffix=".fq", delete=False
            ) as r1_temp,
            tempfile.NamedTemporaryFile(
                mode="w", suffix=".fq", delete=False
            ) as r2_temp,
            tempfile.NamedTemporaryFile(
                mode="w", suffix=".bam", delete=False
            ) as bam_temp,
        ):
            for handle in (vcf_temp, cov_temp, r1_temp, r2_temp, bam_temp):
                handle.write("dummy content")

            df = pd.DataFrame(
                {
                    "sample_name": ["sample1"],
                    "sample_number": [0],
                    "reads1": [r1_temp.name],
                    "reads2": [r2_temp.name],
                    "bam": [bam_temp.name],
                    "vcf": [vcf_temp.name],
                    "coverage": [cov_temp.name],
                }
            )

            # Should not raise any exception
            validate_input_file(df)

        for path in (
            vcf_temp.name,
            cov_temp.name,
            r1_temp.name,
            r2_temp.name,
            bam_temp.name,
        ):
            os.unlink(path)

    def test_missing_columns(self):
        """Test validation with missing required columns."""
        df = pd.DataFrame(
            {
                "sample_name": ["sample1"],
                "sample_number": [0],
            }
        )

        with pytest.raises(ValueError, match="Missing required columns"):
            validate_input_file(df)

    def test_empty_values(self):
        """Test validation with empty values in required columns."""
        df = pd.DataFrame(
            {
                "sample_name": ["sample1"],
                "sample_number": [0],
                "reads1": [""],
                "reads2": ["reads2.fastq"],
                "bam": ["sample1.bam"],
                "vcf": ["sample1.vcf"],
                "coverage": ["sample1_coverage.txt"],
            }
        )

        with pytest.raises(ValueError, match="contains empty values"):
            validate_input_file(df)

    def test_nonexistent_files(self):
        """Test validation with non-existent files."""
        df = pd.DataFrame(
            {
                "sample_name": ["sample1"],
                "sample_number": [0],
                "reads1": ["/nonexistent/reads1.fastq"],
                "reads2": ["/nonexistent/reads2.fastq"],
                "bam": ["/nonexistent/sample1.bam"],
                "vcf": ["/nonexistent/file.vcf"],
                "coverage": ["/nonexistent/file.txt"],
            }
        )

        with pytest.raises(FileNotFoundError):
            validate_input_file(df)


class TestCheckDependencies:
    """Test check_dependencies function."""

    @patch("shutil.which")
    def test_all_dependencies_available(self, mock_which):
        """Test when all dependencies are available."""
        mock_which.return_value = "/usr/bin/tool"  # Simulate tool found

        result = check_dependencies()

        assert isinstance(result, dict)
        assert "bcftools" in result
        assert "tabix" in result
        assert result["bcftools"] is True
        assert result["tabix"] is True

    @patch("shutil.which")
    def test_missing_dependencies(self, mock_which):
        """Test when some dependencies are missing."""

        def mock_which_side_effect(tool):
            if tool == "bcftools":
                return "/usr/bin/bcftools"
            else:  # tabix
                return None

        mock_which.side_effect = mock_which_side_effect

        result = check_dependencies()

        assert result["bcftools"] is True
        assert result["tabix"] is False


class TestExceptions:
    """Test custom exception classes."""

    def test_dependency_error(self):
        """Test DependencyError exception."""
        error = DependencyError("Test dependency error")
        assert str(error) == "Test dependency error"
        assert isinstance(error, Exception)

    def test_input_validation_error(self):
        """Test InputValidationError exception."""
        error = InputValidationError("Test validation error")
        assert str(error) == "Test validation error"
        assert isinstance(error, Exception)

    def test_processing_error(self):
        """Test ProcessingError exception."""
        error = ProcessingError("Test processing error")
        assert str(error) == "Test processing error"
        assert isinstance(error, Exception)
