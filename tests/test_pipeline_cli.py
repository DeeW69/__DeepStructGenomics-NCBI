"""Exercise sequence commands and NCBI errors without network access."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from unittest.mock import Mock

import pytest
import requests


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_pipeline.py"


@pytest.fixture
def cli(tmp_path, monkeypatch):
    def reject_network(*args, **kwargs):
        pytest.fail("CLI tests must run offline")

    monkeypatch.setattr(requests.Session, "request", reject_network)
    monkeypatch.chdir(tmp_path)
    spec = importlib.util.spec_from_file_location("pipeline_cli_under_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("label", [None, "my_comparison"])
def test_fasta_command_exports_comparison_with_default_or_explicit_label(tmp_path, label):
    reference = tmp_path / "reference sequence.fa"
    mutant = tmp_path / "mutant sequence.fa"
    reference.write_text(">reference DNA example\naT\ngc\n", encoding="utf-8-sig")
    mutant.write_text(">mutant RNA example\nAUGC\n", encoding="utf-8")
    output_dir = tmp_path / "results with spaces"
    command = [
        sys.executable, "-B", str(SCRIPT), "--fasta", str(reference),
        "--mutant-fasta", str(mutant), "--output-dir", str(output_dir),
    ]
    if label:
        command.extend(["--label", label])

    completed = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, timeout=30)

    assert completed.returncode == 0, completed.stderr
    identifier = label or "reference"
    report = json.loads((output_dir / f"{identifier}.json").read_text(encoding="utf-8"))
    assert report["identifier"] == identifier
    assert report["sequence"] == "AUGC"
    assert report["annotations"]["au_ratio"] == 0.5
    assert report["variant_analysis"]["total_differences"] == 0
    assert report["annotations"]["metadata"]["source"] == "fasta"
    assert report["annotations"]["metadata"]["fasta_identifier"] == "reference"
    assert (output_dir / f"{identifier}.md").is_file()
    assert (output_dir / identifier / "visualization" / "visualization_manifest.json").is_file()


@pytest.mark.parametrize("reference_from_fasta", [False, True])
def test_raw_and_fasta_sources_can_be_combined(cli, tmp_path, monkeypatch, reference_from_fasta):
    fasta = tmp_path / "sequence.fa"
    fasta.write_text(">sample\nATGC\n", encoding="utf-8")
    arguments = (
        ["--fasta", str(fasta), "--mutant-sequence", "AUGC"]
        if reference_from_fasta
        else ["--sequence", "AUGC", "--mutant-fasta", str(fasta)]
    )
    output_dir = tmp_path / "reports"
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), *arguments, "--output-dir", str(output_dir)])

    cli.main()

    identifier = "sample" if reference_from_fasta else "custom_sequence"
    report = json.loads((output_dir / f"{identifier}.json").read_text(encoding="utf-8"))
    assert report["variant_analysis"]["total_differences"] == 0


@pytest.mark.parametrize(
    "arguments",
    [
        [],
        ["--sequence", "AUGC", "--fasta", "reference.fa"],
        ["--accession", "NC_TEST", "--fasta", "reference.fa"],
        ["--accession", "NC_TEST", "--sequence", "AUGC"],
        ["--sequence", "AUGC", "--mutant-sequence", "AUGC", "--mutant-fasta", "mutant.fa"],
    ],
)
def test_cli_rejects_ambiguous_or_missing_sources(cli, tmp_path, monkeypatch, capsys, arguments):
    output_dir = tmp_path / "reports"
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), *arguments, "--output-dir", str(output_dir)])

    with pytest.raises(SystemExit) as error:
        cli.main()

    assert error.value.code == 2
    assert "error:" in capsys.readouterr().err
    assert not list(output_dir.rglob("*"))


@pytest.mark.parametrize("endpoint", ["efetch.fcgi", "esummary.fcgi"])
@pytest.mark.parametrize("failure", ["http", "timeout", "connection"])
def test_ncbi_request_failures_stop_before_prediction_and_export(
    cli, tmp_path, monkeypatch, capsys, endpoint, failure
):
    api_key = "fake-key-for-offline-test"
    email = "offline@example.invalid"
    calls = []

    def respond(session, method, url, **kwargs):
        calls.append(url.rsplit("/", 1)[-1])
        prepared = requests.Request(method, url, params=kwargs["params"]).prepare()
        response = requests.Response()
        response.request = prepared
        response.url = prepared.url
        response.status_code = 200
        response._content = b">NC_TEST.1 DNA reference\nATGC\n"
        if url.endswith(endpoint):
            if failure == "http":
                response.status_code = 503
                response._content = b"Service Unavailable"
            elif failure == "timeout":
                raise requests.Timeout(f"Timeout: {prepared.url}", request=prepared)
            else:
                raise requests.ConnectionError(f"Connection failed: {prepared.url}", request=prepared)
        return response

    monkeypatch.setattr(requests.Session, "request", respond)
    predict = Mock(side_effect=AssertionError("Prediction must wait for a complete NCBI record"))
    monkeypatch.setattr("deepstructgenomics.pipeline.predict_secondary_structure", predict)
    output_dir = tmp_path / "reports"
    monkeypatch.setattr(sys, "argv", [
        str(SCRIPT), "--accession", "NC_TEST.1", "--output-dir", str(output_dir),
        "--ncbi-api-key", api_key, "--ncbi-email", email,
    ])

    with pytest.raises(SystemExit) as error:
        cli.main()

    assert "Erreur NCBI" in str(error.value.code)
    captured = capsys.readouterr()
    displayed = str(error.value.code) + captured.out + captured.err
    assert api_key not in displayed
    assert email not in displayed
    assert "Traceback" not in displayed
    assert "Rapports generes" not in displayed
    assert calls == (["efetch.fcgi"] if endpoint == "efetch.fcgi" else ["efetch.fcgi", "esummary.fcgi"])
    predict.assert_not_called()
    assert not list(output_dir.rglob("*"))
    assert not list((tmp_path / "outputs").rglob("*"))


def test_ncbi_empty_fasta_is_a_readable_input_error_without_reports(cli, tmp_path, monkeypatch, capsys):
    response = requests.Response()
    response.status_code = 200
    response._content = b""
    request = Mock(return_value=response)
    monkeypatch.setattr(requests.Session, "request", request)
    output_dir = tmp_path / "reports"
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "--accession", "NC_TEST.1", "--output-dir", str(output_dir)])

    with pytest.raises(SystemExit) as error:
        cli.main()

    assert error.value.code == 2
    captured = capsys.readouterr()
    assert "FASTA vide" in captured.err
    assert "Traceback" not in captured.err
    assert "Rapports generes" not in captured.out
    assert request.call_count == 1
    assert not list(output_dir.rglob("*"))


def test_ncbi_empty_summary_still_exports_the_valid_sequence(cli, tmp_path, monkeypatch, capsys):
    def respond(session, method, url, **kwargs):
        response = requests.Response()
        response.status_code = 200
        if url.endswith("efetch.fcgi"):
            response._content = b">NC_TEST.1 DNA reference\nATGC\n"
        elif url.endswith("esummary.fcgi"):
            response._content = b'{"result": null}'
        else:
            pytest.fail(f"Unexpected NCBI endpoint: {url}")
        return response

    monkeypatch.setattr(requests.Session, "request", respond)
    output_dir = tmp_path / "reports"
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "--accession", "NC_TEST.1", "--output-dir", str(output_dir)])

    cli.main()

    report = json.loads((output_dir / "NC_TEST.1.json").read_text(encoding="utf-8"))
    assert report["sequence"] == "AUGC"
    assert report["description"] == "DNA reference"
    metadata = report["annotations"]["metadata"]
    assert set(metadata) == {"provenance"}
    assert metadata["provenance"]["requested_accession"] == "NC_TEST.1"
    assert (output_dir / "NC_TEST.1.md").is_file()
    assert (output_dir / "NC_TEST.1/visualization/visualization_manifest.json").is_file()
    assert "Rapports generes" in capsys.readouterr().out


@pytest.mark.parametrize(
    "input_kind",
    ["missing", "multiple", "empty", "invalid-encoding", "ambiguous-raw", "empty-mutant", "ambiguous-mutant"],
)
def test_cli_input_errors_are_readable_and_leave_no_reports(cli, tmp_path, monkeypatch, capsys, input_kind):
    fasta = tmp_path / "input.fa"
    arguments = ["--fasta", str(fasta)]
    if input_kind == "multiple":
        fasta.write_text(">one\nAUGC\n>two\nAUGC\n", encoding="utf-8")
    elif input_kind == "empty":
        fasta.write_text(">one\n", encoding="utf-8")
    elif input_kind == "invalid-encoding":
        fasta.write_bytes(b">one\nAT\xffGC\n")
    elif input_kind == "ambiguous-raw":
        arguments = ["--sequence", "AUGNC"]
    elif input_kind == "empty-mutant":
        arguments = ["--sequence", "AUGC", "--mutant-sequence", ""]
    elif input_kind == "ambiguous-mutant":
        arguments = ["--sequence", "AUGC", "--mutant-sequence", "AUGNC"]
    output_dir = tmp_path / "reports"
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), *arguments, "--output-dir", str(output_dir)])

    with pytest.raises(SystemExit) as error:
        cli.main()

    assert error.value.code == 2
    stderr = capsys.readouterr().err
    assert "Erreur :" in stderr
    assert "Traceback" not in stderr
    if input_kind.startswith("ambiguous"):
        assert "N" in stderr and "4" in stderr
    assert not list(output_dir.rglob("*"))


def test_cli_exports_csv_and_applies_thresholds(cli, tmp_path, monkeypatch, capsys):
    output_dir = tmp_path / "reports"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT),
            "--sequence", "GAAAAC",
            "--mutant-sequence", "AAAAAC",
            "--output-dir", str(output_dir),
            "--top-k", "1",
            "--min-abs-delta", "0.65",
            "--base-pair-threshold", "0.65",
        ],
    )

    cli.main()

    stdout = capsys.readouterr().out
    assert "- CSV      :" in stdout
    csv_file = output_dir / "custom_sequence_hotspots.csv"
    assert csv_file.is_file()
    lines = csv_file.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0] == "position,reference,mutant,delta,abs_delta"
    # Only position 1 should be present because top-k=1 and min_abs_delta=0.65
    assert len(lines) == 2
    assert lines[1] == "1,G,A,-0.7,0.7"

    json_report = json.loads((output_dir / "custom_sequence.json").read_text(encoding="utf-8"))
    assert json_report["impact_summary"]["parameters"]["top_k"] == 1
    assert json_report["impact_summary"]["parameters"]["min_abs_delta"] == 0.65
    assert len(json_report["impact_summary"]["hotspots"]) == 1


def test_cli_accepts_delta_threshold_alias(cli, tmp_path, monkeypatch):
    output_dir = tmp_path / "reports"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT),
            "--sequence", "GAAAAC",
            "--mutant-sequence", "AAAAAC",
            "--output-dir", str(output_dir),
            "--delta-threshold", "0.65",
        ],
    )

    cli.main()

    json_report = json.loads((output_dir / "custom_sequence.json").read_text(encoding="utf-8"))
    assert json_report["impact_summary"]["parameters"]["min_abs_delta"] == 0.65


@pytest.mark.parametrize("invalid_arg", [["--top-k", "-1"], ["--min-abs-delta", "-0.1"], ["--base-pair-threshold", "-0.2"]])
def test_cli_rejects_negative_thresholds(cli, tmp_path, monkeypatch, capsys, invalid_arg):
    output_dir = tmp_path / "reports"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT),
            "--sequence", "GAAAAC",
            "--output-dir", str(output_dir),
            *invalid_arg,
        ],
    )

    with pytest.raises(SystemExit) as error:
        cli.main()

    assert error.value.code == 2
    assert "Erreur :" in capsys.readouterr().err
