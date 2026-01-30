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


def test_cli_version_option():
    exit_code, out = _run_main("--version")
    assert exit_code == 0
    assert vartracker.__version__ in out


def test_cli_describe_output():
    exit_code, out = _run_main("describe-output")
    assert exit_code == 0
    assert "Output schema" in out


def test_cli_describe_output_writes_json(tmp_path):
    output_path = tmp_path / "schema.json"
    exit_code, out = _run_main(
        "describe-output", "--out", str(output_path), "--format", "json"
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
