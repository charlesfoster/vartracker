"""Utility for preparing pokay literature data for vartracker."""

from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from typing import Iterable, List, Sequence
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pandas as pd
from urlextract import URLExtract

DEFAULT_DOWNLOAD_SUBDIR = os.path.join("pokay_literature", "NC_045512")
GITHUB_CONTENTS_URL = (
    "https://api.github.com/repos/nodrogluap/pokay/contents/data/literature/NC_045512"
)

# Mapping of product names to nsp names
PRODUCT_TO_NSP = {
    "PLpro": "nsp3",
    "3CLpro": "nsp5",
    "3CL": "nsp5",
    "RdRp": "nsp12",
    "Hel": "nsp13",
    "ExoN": "nsp14",
    "NendoU": "nsp15",
    "2'-O-MT": "nsp16",
}


def parse_mut_file(fname: str) -> pd.DataFrame:
    """Split a list into sub lists, depending on a delimiter."""
    # 1: read the file
    print(f"Working with: {fname}")
    with open(fname, "r") as f:
        thelist = [line.rstrip() for line in f.readlines()]
    # 2: get info from filename, set up URLextractor
    bname = os.path.basename(fname)
    gene = bname.split("_")[0]
    # Update gene name if it matches a product
    gene = PRODUCT_TO_NSP.get(gene, gene)  # Convert to nsp if in mapping
    category = "_".join([x.replace(".txt", "") for x in bname.split("_")[1:]])
    delimiters = ""
    extractor = URLExtract()
    # 3: get info from the file based on: https://stackoverflow.com/questions/45281189/split-list-into-lists-based-on-a-character-occurring-inside-of-an-element
    results = []
    sublist: List[str] = []
    for item in thelist:
        if item in delimiters:
            results.append(sublist)
            sublist = []
        else:
            sublist.append(item)
    if sublist:
        results.append(sublist)
    results = [x for x in results if len(x) > 0]
    # 4 organise the information
    mut_list = []
    for r in results:
        mut = r.pop()
        flat = "".join([x.replace("#", "") for x in r])
        extraction = extractor.find_urls(flat)
        reference = " ; ".join(extraction) if len(extraction) > 0 else "N/A"
        info = {
            "gene": gene,
            "category": category,
            "mutation": mut,
            "information": flat,
            "reference": reference,
        }
        mut_list.append(info)
    # 5: return the organised dictionary
    return pd.DataFrame(mut_list)


def has_pokay_text_files(path: str) -> bool:
    """Return True if the directory contains pokay literature text files."""
    return os.path.isdir(path) and bool(glob.glob(os.path.join(path, "*.txt")))


def _fetch_bytes(url: str, headers: dict[str, str] | None = None) -> bytes:
    """Fetch bytes from the given URL, falling back to curl if SSL is unavailable."""

    request_headers = headers or {}
    request = Request(url, headers=request_headers)

    try:
        with urlopen(request) as response:
            return response.read()
    except HTTPError:
        raise
    except URLError as exc:
        reason = getattr(exc, "reason", "")
        if (
            isinstance(reason, str)
            and "unknown url type" in reason
            and shutil.which("curl")
        ):
            curl_cmd = ["curl", "-fsSL"]
            for key, value in request_headers.items():
                curl_cmd.extend(["-H", f"{key}: {value}"])
            curl_cmd.append(url)
            result = subprocess.run(
                curl_cmd, capture_output=True, check=False
            )  # pragma: no cover - network fallback only
            if result.returncode == 0:
                return result.stdout
            raise RuntimeError(
                f"curl failed to fetch {url}: {result.stderr.decode('utf-8', 'ignore').strip()}"
            ) from exc
        raise


def download_pokay_files(target_dir: str) -> None:
    """Download pokay literature text files into the target directory."""
    os.makedirs(target_dir, exist_ok=True)

    headers = {
        "Accept": "application/vnd.github.v3+json",
        "User-Agent": "vartracker",
    }

    try:
        payload = _fetch_bytes(GITHUB_CONTENTS_URL, headers=headers)
    except HTTPError as exc:  # pragma: no cover - network failure
        raise RuntimeError(
            f"Failed to list pokay literature files: {exc.code} {exc.reason}"
        ) from exc
    except URLError as exc:  # pragma: no cover - network failure
        raise RuntimeError(f"Failed to reach GitHub: {exc.reason}") from exc
    except RuntimeError as exc:  # pragma: no cover - curl failure
        raise RuntimeError(str(exc)) from exc

    try:
        entries = json.loads(payload.decode("utf-8"))
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            "Unexpected response while listing pokay literature files."
        ) from exc

    txt_entries = [entry for entry in entries if entry.get("name", "").endswith(".txt")]
    if not txt_entries:
        raise RuntimeError("No .txt files found in the pokay literature repository.")

    for entry in txt_entries:
        download_url = entry.get("download_url")
        file_name = entry.get("name")
        if not download_url or not file_name:
            continue

        destination = os.path.join(target_dir, file_name)
        try:
            payload = _fetch_bytes(download_url)
        except HTTPError as exc:  # pragma: no cover - network failure
            raise RuntimeError(
                f"Failed to download {file_name}: {exc.code} {exc.reason}"
            ) from exc
        except URLError as exc:  # pragma: no cover - network failure
            raise RuntimeError(f"Failed to download {file_name}: {exc.reason}") from exc
        except RuntimeError as exc:  # pragma: no cover - curl failure
            raise RuntimeError(str(exc)) from exc

        with open(destination, "wb") as handle:
            handle.write(payload)


@dataclass
class ParsePokayConfig:
    """Configuration returned by CLI parsing."""

    source_dir: str
    output_file: str
    should_download: bool


def parse_arguments(argv: Sequence[str]) -> ParsePokayConfig:
    """Parse command-line arguments for parse_pokay."""

    parser = argparse.ArgumentParser(
        description=(
            "Prepare pokay literature annotations. If only an output file is provided "
            "the required text files are downloaded automatically from GitHub."
        )
    )
    parser.add_argument(
        "--download-dir",
        default=DEFAULT_DOWNLOAD_SUBDIR,
        help=(
            "Directory used when downloading the pokay literature files automatically "
            f"(default: {DEFAULT_DOWNLOAD_SUBDIR})."
        ),
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Always re-download the pokay literature files into the target directory.",
    )
    parser.add_argument(
        "positional",
        nargs="*",
        help="Provide OUTPUT_FILE or DATA_DIR OUTPUT_FILE for backwards compatibility.",
    )

    args = parser.parse_args(list(argv))

    positional: List[str] = args.positional
    if len(positional) == 1:
        output_file = positional[0]
        source_dir = args.download_dir
        should_download = True
    elif len(positional) == 2:
        source_dir, output_file = positional
        should_download = args.force_download or not has_pokay_text_files(source_dir)
    else:
        parser.error("Provide OUTPUT_FILE or DATA_DIR OUTPUT_FILE.")

    source_dir = os.path.abspath(source_dir)
    output_file = os.path.abspath(output_file)

    return ParsePokayConfig(
        source_dir=source_dir, output_file=output_file, should_download=should_download
    )


def build_dataframe(file_paths: Iterable[str]) -> pd.DataFrame:
    """Aggregate pokay mutation data from the given text files."""
    df_list = [parse_mut_file(path) for path in file_paths]
    return pd.concat(df_list).drop_duplicates().reset_index(drop=True)


def main(argv: Sequence[str] | None = None) -> int:
    """Main entry point for the parse_pokay command."""

    if argv is None:
        argv = sys.argv[1:]

    config = parse_arguments(argv)

    if config.should_download:
        print(f"Downloading pokay literature files to {config.source_dir}...")
        download_pokay_files(config.source_dir)

    files = sorted(glob.glob(os.path.join(config.source_dir, "*.txt")))
    if not files:
        print(
            f"No .txt files found in {config.source_dir}. Provide a directory containing "
            "pokay literature text files or allow the script to download them.",
            file=sys.stderr,
        )
        return 1

    combined = build_dataframe(files)
    combined.to_csv(config.output_file, index=None)

    print(f"Parsed {len(files)} files, found {len(combined)} unique mutations.")
    print(f"Results written to {config.output_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
