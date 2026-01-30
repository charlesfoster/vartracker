"""
Core utilities and data handling for vartracker.

Contains utility functions and data processing routines.
"""

import os

import pandas as pd

from .data import get_data_file
from typing import Iterable, Sequence


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


DEFAULT_INPUT_COLUMNS: Sequence[str] = (
    "sample_name",
    "sample_number",
    "reads1",
    "reads2",
    "bam",
    "vcf",
    "coverage",
)

FILE_COLUMNS = {"reads1", "reads2", "bam", "vcf", "coverage"}


def validate_input_file(
    input_file,
    *,
    required_columns: Iterable[str] | None = None,
    optional_empty: Iterable[str] | None = None,
):
    """
    Validate the input CSV file format.

    Args:
        input_file (pd.DataFrame): Input data frame

    Raises:
        ValueError: If required columns are missing
    """
    required = list(required_columns or DEFAULT_INPUT_COLUMNS)
    optional_empty = set(optional_empty or set())

    missing_columns = [col for col in required if col not in input_file.columns]

    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    # Check for empty values
    for col in required:
        if col in optional_empty:
            continue
        series = input_file[col]
        empty_mask = series.isna() | series.astype(str).str.strip().eq("")
        if empty_mask.any():
            raise ValueError(f"Column '{col}' contains empty values")

    # Check file existence for vcf and coverage files
    for _, row in input_file.iterrows():
        for col in FILE_COLUMNS:
            if col not in input_file.columns:
                continue
            value = row[col]
            if pd.isna(value):
                if col in optional_empty:
                    continue
                raise ValueError(f"Column '{col}' contains empty values")
            path = str(value).strip()
            if not path:
                if col in optional_empty:
                    continue
                raise ValueError(f"Column '{col}' contains empty values")
            if not os.path.exists(path):
                raise FileNotFoundError(f"File not found for column '{col}': {path}")


def check_dependencies(mode: str = "vcf"):
    """Check required external tools based on mode."""

    import shutil

    mode = mode.lower()
    if mode == "vcf":
        requirements = ["bcftools", "tabix"]
    elif mode == "bam":
        requirements = ["samtools", "bcftools", "lofreq", "tabix"]
    elif mode == "e2e":
        requirements = ["samtools", "bcftools", "tabix", "lofreq", "fastp", "bwa"]
    else:
        requirements = ["bcftools", "tabix"]

    return {tool: shutil.which(tool) is not None for tool in requirements}


class VartrackerError(Exception):
    """Base exception class for vartracker."""

    pass


class DependencyError(VartrackerError):
    """Raised when required dependencies are missing."""

    def __init__(
        self,
        message: str,
        *,
        mode: str | None = None,
        missing: Iterable[str] | None = None,
        tip: str | None = None,
    ):
        super().__init__(message)
        self.mode = mode
        self.missing = list(missing or [])
        self.tip = tip


class InputValidationError(VartrackerError):
    """Raised when input validation fails."""

    pass


class ProcessingError(VartrackerError):
    """Raised when processing fails."""

    pass
