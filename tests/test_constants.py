"""Tests for reference-coordinate formatting helpers."""

from __future__ import annotations

from vartracker.constants import reformat_csq_notation


def test_reformat_csq_notation_preserves_orf1ab_stop_gained_nsp_change():
    reformatted, nsp_change = reformat_csq_notation("ORF1ab", "889L>889*")

    assert reformatted == "L889*"
    assert nsp_change == "nsp3_PLpro:L71*"
