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
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".vcf", delete=False
        ) as vcf_temp:
            vcf_temp.write("dummy content")
            vcf_path = vcf_temp.name

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as cov_temp:
            cov_temp.write("dummy content")
            cov_path = cov_temp.name

        try:
            df = pd.DataFrame(
                {
                    "vcf": [vcf_path],
                    "coverage": [cov_path],
                    "sample_name": ["sample1"],
                    "sample_number": [0],
                }
            )

            # Should not raise any exception
            validate_input_file(df)
        finally:
            # Cleanup
            os.unlink(vcf_path)
            os.unlink(cov_path)

    def test_missing_columns(self):
        """Test validation with missing required columns."""
        df = pd.DataFrame(
            {
                "vcf": ["test.vcf"],
                "coverage": ["test.txt"],
                # Missing sample_name and sample_number
            }
        )

        with pytest.raises(ValueError, match="Missing required columns"):
            validate_input_file(df)

    def test_empty_values(self):
        """Test validation with empty values in required columns."""
        df = pd.DataFrame(
            {
                "vcf": ["test.vcf"],
                "coverage": [None],  # Empty value
                "sample_name": ["sample1"],
                "sample_number": [0],
            }
        )

        with pytest.raises(ValueError, match="contains empty values"):
            validate_input_file(df)

    def test_nonexistent_files(self):
        """Test validation with non-existent files."""
        df = pd.DataFrame(
            {
                "vcf": ["/nonexistent/file.vcf"],
                "coverage": ["/nonexistent/file.txt"],
                "sample_name": ["sample1"],
                "sample_number": [0],
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
