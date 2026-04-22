import pytest

from vartracker.consensus import iupac_consensus_base, simple_consensus_base


@pytest.mark.parametrize(
    ("af", "simple", "iupac"),
    [
        (0.20, "A", "A"),
        (0.30, "A", "R"),
        (0.50, "A", "R"),
        (0.80, "G", "G"),
    ],
)
def test_snp_consensus_three_tier_thresholds(af, simple, iupac):
    kwargs = {
        "depth": 100,
        "af": af,
        "min_depth": 10,
        "consensus_snp_min_af": 0.25,
        "consensus_snp_thresh": 0.75,
        "consensus_indel_thresh": 0.75,
    }

    assert simple_consensus_base("A", "G", **kwargs) == simple
    assert iupac_consensus_base("A", "G", **kwargs) == iupac


@pytest.mark.parametrize(
    ("af", "simple"),
    [
        (0.74, "AT"),
        (0.75, "A"),
        (0.80, "A"),
    ],
)
def test_indel_af_cutoff_is_separate_from_snp_ambiguity(af, simple):
    kwargs = {
        "depth": 100,
        "af": af,
        "min_depth": 10,
        "consensus_snp_min_af": 0.25,
        "consensus_snp_thresh": 0.75,
        "consensus_indel_thresh": 0.75,
    }

    assert simple_consensus_base("AT", "A", **kwargs) == simple
    assert iupac_consensus_base("AT", "A", **kwargs) == simple
