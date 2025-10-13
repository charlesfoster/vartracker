"""
vartracker: Track the persistence (or loss) of mutations during long-term passaging.

A bioinformatics tool for analyzing longitudinal variant data, particularly
designed for SARS-CoV-2 long-term passaging experiments.

Author: Dr Charles Foster
"""

from ._version import __version__

__author__ = "Dr Charles Foster"
__email__ = "github.com/charlesfoster"

from .main import main  # noqa: F401

__all__ = ["__version__", "main"]
