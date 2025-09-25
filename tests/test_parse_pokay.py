"""Tests for the parse_pokay utility."""

from __future__ import annotations

import json
from io import BytesIO

import pandas as pd

from vartracker.data import parse_pokay


class DummyResponse(BytesIO):
    """Simple context manager around BytesIO for mocking urlopen."""

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


def test_parse_arguments_auto_download(tmp_path):
    download_dir = tmp_path / "remote"
    config = parse_pokay.parse_arguments(
        ["--download-dir", str(download_dir), "pokay.csv"]
    )

    assert config.source_dir == str(download_dir.resolve())
    assert config.output_file.endswith("pokay.csv")
    assert config.should_download is True


def test_parse_arguments_existing_dir(tmp_path):
    data_dir = tmp_path / "pokay"
    data_dir.mkdir()
    (data_dir / "example.txt").write_text("example", encoding="utf-8")

    config = parse_pokay.parse_arguments([str(data_dir), "out.csv"])

    assert config.source_dir == str(data_dir.resolve())
    assert config.should_download is False


def test_parse_arguments_missing_dir_triggers_download(tmp_path):
    data_dir = tmp_path / "missing"

    config = parse_pokay.parse_arguments([str(data_dir), "out.csv"])

    assert config.source_dir == str(data_dir.resolve())
    assert config.should_download is True


def test_download_pokay_files(monkeypatch, tmp_path):
    target_dir = tmp_path / "pokay"

    listing = json.dumps(
        [
            {
                "name": "file1.txt",
                "download_url": "https://example.com/file1.txt",
            },
            {
                "name": "README.md",
                "download_url": "https://example.com/README.md",
            },
        ]
    ).encode("utf-8")

    file_payload = b"mutation data"

    def fake_urlopen(request):
        url = request
        if hasattr(request, "full_url"):
            url = request.full_url
        if url == parse_pokay.GITHUB_CONTENTS_URL:
            return DummyResponse(listing)
        if url == "https://example.com/file1.txt":
            return DummyResponse(file_payload)
        raise AssertionError(f"Unexpected URL requested: {url}")

    monkeypatch.setattr(parse_pokay, "urlopen", fake_urlopen)

    parse_pokay.download_pokay_files(str(target_dir))

    downloaded = target_dir / "file1.txt"
    assert downloaded.exists()
    assert downloaded.read_bytes() == file_payload


def test_build_dataframe(tmp_path):
    file_path = tmp_path / "example.txt"
    file_path.write_text(
        """nsp3_PLpro_category
#Info line
example info http://example.com
Mutation
""",
        encoding="utf-8",
    )

    df = parse_pokay.build_dataframe([str(file_path)])
    assert isinstance(df, pd.DataFrame)
    assert {"gene", "category", "mutation", "information", "reference"} <= set(
        df.columns
    )
