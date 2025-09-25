"""
Amino acid analysis functionality for vartracker.

Contains classes and functions for analyzing amino acid changes and their properties.
"""

import re


class AminoAcidChange:
    """Analyze amino acid changes and their biochemical properties."""

    def __init__(self, aa_change):
        """
        Initialize amino acid change analysis.

        Args:
            aa_change (str): Amino acid change in format like 'D614G'
        """
        self.amino_acids = {
            "I": ["Hydrophobic", "Aliphatic", "Uncharged"],
            "L": ["Hydrophobic", "Aliphatic", "Uncharged"],
            "V": ["Hydrophobic", "Small", "Aliphatic", "Uncharged"],
            "C": ["Hydrophobic", "Small", "Uncharged"],
            "A": ["Hydrophobic", "Small", "Tiny", "Uncharged"],
            "G": ["Hydrophobic", "Small", "Tiny", "Uncharged"],
            "M": ["Hydrophobic", "Uncharged"],
            "F": ["Hydrophobic", "Aromatic", "Uncharged"],
            "Y": ["Hydrophobic", "Polar", "Aromatic", "Uncharged"],
            "W": ["Hydrophobic", "Polar", "Aromatic", "Uncharged"],
            "H": ["Hydrophobic", "Polar", "Aromatic", "Positive", "Charged"],
            "K": ["Hydrophobic", "Polar", "Positive", "Charged"],
            "R": ["Polar", "Positive", "Charged"],
            "E": ["Polar", "Negative", "Charged"],
            "Q": ["Polar", "Uncharged"],
            "D": ["Polar", "Proline", "Negative", "Charged"],
            "N": ["Polar", "Proline", "Uncharged"],
            "S": ["Polar", "Proline", "Tiny", "Uncharged"],
            "T": ["Hydrophobic", "Polar", "Proline", "Uncharged"],
            "P": ["Proline", "Uncharged"],
            "B": ["Polar", "Uncharged"],
            "Z": ["Polar", "Uncharged"],
            "X": [
                "Hydrophobic",
                "Polar",
                "Tiny",
                "Aliphatic",
                "Aromatic",
                "Positive",
                "Negative",
                "Charged",
                "Uncharged",
            ],
            "-": [
                "Hydrophobic",
                "Polar",
                "Tiny",
                "Aliphatic",
                "Aromatic",
                "Positive",
                "Negative",
                "Charged",
                "Uncharged",
            ],
        }

        self.weights = {
            "A": 71.037114,
            "R": 156.101111,
            "N": 114.042927,
            "D": 115.026943,
            "C": 103.009185,
            "E": 129.042593,
            "Q": 128.058578,
            "G": 57.021464,
            "H": 137.058912,
            "I": 113.084064,
            "L": 113.084064,
            "K": 128.094963,
            "M": 131.040485,
            "F": 147.068414,
            "P": 97.052764,
            "S": 87.032028,
            "T": 101.047679,
            "U": 150.95363,
            "W": 186.079313,
            "Y": 163.06332,
            "V": 99.068414,
            "X": 0,
            "B": 113.084064,
            "Z": 0,
        }

        self._parse_amino_acid_change(aa_change)

    def _parse_amino_acid_change(self, aa_change):
        """Parse the amino acid change string and calculate properties."""
        try:
            match = re.match(r"([a-z]+)([0-9]+)([a-z+])", aa_change, re.I)
            if match:
                self.aa1 = match.group(1)
                self.aa2 = match.group(3)
                self.change = aa_change

                self.aa1_total_properties = set(
                    item
                    for sublist in [self.amino_acids.get(x, []) for x in self.aa1]
                    for item in sublist
                )
                self.aa2_total_properties = set(
                    item
                    for sublist in [self.amino_acids.get(x, []) for x in self.aa2]
                    for item in sublist
                )

                self.aa1_unique_properties = list(
                    self.aa1_total_properties - self.aa2_total_properties
                )
                self.aa2_unique_properties = list(
                    self.aa2_total_properties - self.aa1_total_properties
                )

                self.aa1_weight = float(
                    "{:.3f}".format(sum(self.weights.get(x, 0) for x in self.aa1))
                )
                self.aa2_weight = float(
                    "{:.3f}".format(sum(self.weights.get(x, 0) for x in self.aa2))
                )
                self.weight_difference = round(
                    (self.aa1_weight - self.aa2_weight) * -1, 3
                )
            else:
                raise ValueError(f"Invalid amino acid change format: {aa_change}")

        except Exception:
            # Fallback for invalid formats
            self.change = aa_change
            self.aa1 = ""
            self.aa2 = ""
            self.aa1_total_properties = set()
            self.aa2_total_properties = set()
            self.aa1_unique_properties = []
            self.aa2_unique_properties = []
            self.aa1_weight = ""
            self.aa2_weight = ""
            self.weight_difference = 0.0

    def __repr__(self):
        """Return string representation of the amino acid change analysis."""
        return f"AminoAcidChange('{self.change}')"

    def __str__(self):
        """Return detailed string representation."""
        return f"Information for the following change: {self.change}\n{self.__dict__}"
