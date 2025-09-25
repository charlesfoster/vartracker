"""Data files for vartracker package."""

from pathlib import Path

DATA_DIR = Path(__file__).parent


def get_data_file(filename):
    """Get the full path to a data file."""
    return str(DATA_DIR / filename)


# Convenience functions for getting specific data files
def get_reference_genome():
    """Get path to SARS-CoV-2 reference genome."""
    return get_data_file("NC_045512.fasta")


def get_annotation():
    """Get path to SARS-CoV-2 annotations."""
    return get_data_file("NC_045512.gff3")
