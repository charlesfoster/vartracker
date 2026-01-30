"""Tests for output schema metadata."""

from pathlib import Path

import importlib.util
import pandas as pd


def _load_schema_module():
    spec = importlib.util.spec_from_file_location(
        "schemas", Path("vartracker") / "schemas.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_results_schema_matches_precomputed_columns():
    module = _load_schema_module()
    schema_columns = [entry["name"] for entry in module.RESULTS_SCHEMA]
    precomputed = pd.read_csv(
        Path("vartracker") / "test_data" / "precomputed" / "test_results.csv"
    )
    assert schema_columns == list(precomputed.columns)
