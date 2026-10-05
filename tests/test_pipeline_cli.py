"""Exercise FASTA commands and actionable CLI errors without network access."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest
import requests


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_pipeline.py"


@pytest.fixture
def cli(tmp_path, monkeypatch):
    def reject_network(*args, **kwargs):
        pytest.fail("FASTA commands must run offline")

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
