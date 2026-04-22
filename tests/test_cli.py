# tests/test_cli.py
import contextlib
import io
import json

import vartracker
from vartracker.main import main


def _run_main(*args):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        exit_code = main(list(args))
    return exit_code, buf.getvalue()


def test_cli_help_smoke():
    exit_code, out = _run_main("-h")
    assert exit_code == 0
    assert "usage" in out.lower() or "help" in out.lower()


def test_plot_heatmap_help_uses_standalone_option_names():
    exit_code, out = _run_main("plot", "heatmap", "-h")
    assert exit_code == 0
    assert "--aa-exclude" in out
    assert "--include-joint" in out
    assert "--title" in out
    assert "--literature-csv" in out
    assert "--x-labels" in out
    heatmap_options = out.split("Heatmap options:", 1)[1]
    assert "--title" in heatmap_options
    assert "--literature-csv" in heatmap_options
    assert "--x-labels" in heatmap_options
    assert "--heatmap-aa-exclude" not in out
    assert "--heatmap-include-joint" not in out
    assert "--outdir" not in out
    assert "--name" not in out
    assert "--min-snv-freq" not in out
    assert "--min-indel-freq" not in out


def test_cli_version_option():
    exit_code, out = _run_main("--version")
    assert exit_code == 0
    assert vartracker.__version__ in out


def test_cli_describe_output():
    exit_code, out = _run_main("schema", "results")
    assert exit_code == 0
    assert "Output schema" in out


def test_cli_describe_output_writes_json(tmp_path):
    output_path = tmp_path / "schema.json"
    exit_code, out = _run_main(
        "schema", "results", "--out", str(output_path), "--format", "json"
    )
    assert exit_code == 0
    assert output_path.exists()
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert isinstance(payload, list)
    assert payload
    assert "name" in payload[0]


def test_cli_missing_arguments_prints_help():
    exit_code, out = _run_main()
    assert exit_code == 1
    assert "usage" in out.lower()


def test_cli_literature_schema():
    exit_code, out = _run_main("schema", "literature")
    assert exit_code == 0
    assert "Literature schema" in out
    assert "gene:" in out
    assert "reference:" in out


def test_prepare_reference_requires_accessions_source(tmp_path):
    exit_code, out = _run_main(
        "prepare",
        "reference",
        "--outdir",
        str(tmp_path),
    )
    assert exit_code != 0
    assert "one of the arguments --accessions --accession-file is required" in out


def test_prepare_reference_rejects_both_accession_inputs(tmp_path):
    accessions_file = tmp_path / "accessions.txt"
    accessions_file.write_text("CY114381\n", encoding="utf-8")
    exit_code, out = _run_main(
        "prepare",
        "reference",
        "--outdir",
        str(tmp_path),
        "--accessions",
        "CY114381",
        "--accession-file",
        str(accessions_file),
    )
    assert exit_code != 0
    assert "not allowed with argument" in out
