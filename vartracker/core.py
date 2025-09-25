"""
Core utilities and data handling for vartracker.

Contains utility functions and data processing routines.
"""

import os

from .data import get_data_file


def get_package_data_path(filename):
    """
    Get the path to a data file within the package.

    Args:
        filename (str): Name of the data file

    Returns:
        str: Full path to the data file
    """
    return get_data_file(filename)


def get_logo():
    """Return the vartracker ASCII logo."""
    return """
██    ██  █████  ██████  ████████ ██████   █████   ██████ ██   ██ ███████ ██████
██    ██ ██   ██ ██   ██    ██    ██   ██ ██   ██ ██      ██  ██  ██      ██   ██
██    ██ ███████ ██████     ██    ██████  ███████ ██      █████   █████   ██████
 ██  ██  ██   ██ ██   ██    ██    ██   ██ ██   ██ ██      ██  ██  ██      ██   ██
  ████   ██   ██ ██   ██    ██    ██   ██ ██   ██  ██████ ██   ██ ███████ ██   ██


"""


def validate_input_file(input_file):
    """
    Validate the input CSV file format.

    Args:
        input_file (pd.DataFrame): Input data frame

    Raises:
        ValueError: If required columns are missing
    """
    required_columns = ["vcf", "coverage", "sample_name", "sample_number"]
    missing_columns = [col for col in required_columns if col not in input_file.columns]

    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    # Check for empty values
    for col in required_columns:
        if input_file[col].isna().any():
            raise ValueError(f"Column '{col}' contains empty values")

    # Check file existence for vcf and coverage files
    for _, row in input_file.iterrows():
        if not os.path.exists(row["vcf"]):
            raise FileNotFoundError(f"VCF file not found: {row['vcf']}")
        if not os.path.exists(row["coverage"]):
            raise FileNotFoundError(f"Coverage file not found: {row['coverage']}")


def check_dependencies():
    """
    Check if required external tools are available.

    Returns:
        dict: Dictionary of tool availability
    """
    import shutil

    tools = {
        "bcftools": shutil.which("bcftools") is not None,
        "tabix": shutil.which("tabix") is not None,
    }

    return tools


class VartrackerError(Exception):
    """Base exception class for vartracker."""

    pass


class DependencyError(VartrackerError):
    """Raised when required dependencies are missing."""

    pass


class InputValidationError(VartrackerError):
    """Raised when input validation fails."""

    pass


class ProcessingError(VartrackerError):
    """Raised when processing fails."""

    pass
