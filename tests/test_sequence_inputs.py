"""Protect sequence positions across raw, FASTA and NCBI input paths."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import Mock

import pytest
import requests

from deepstructgenomics.config import PipelineConfig
from deepstructgenomics.data_sources.ncbi_client import (
    load_fasta_record,
    normalize_rna_sequence,
    normalize_user_sequence,
    parse_single_fasta,
)
from deepstructgenomics.pipeline import DeepStructPipeline, PipelineInput


@pytest.fixture
def pipeline(tmp_path, monkeypatch):
    def reject_network(*args, **kwargs):
        pytest.fail("Local sequences must not make HTTP requests")

    monkeypatch.setattr(requests.Session, "request", reject_network)
    return DeepStructPipeline(
        PipelineConfig(cache_dir=tmp_path / "cache", default_output_dir=tmp_path / "reports")
    )


def test_raw_dna_and_rna_keep_all_positions():
    assert normalize_rna_sequence(" aT\tg\r\nc u\vA\f") == "AUGCUA"
    assert normalize_rna_sequence("a\u00a0t\u2028gc") == "AUGC"
    record = normalize_user_sequence(" aT\tg\r\nc u\vA\f", "sample")
    assert record.sequence == "AUGCUA"
    assert record.identifier == "sample"
    assert record.metadata["source"] == "user"


@pytest.mark.parametrize("symbol", ["N", "R", "-", "7", ">", "."])
def test_unsupported_bases_report_the_original_position(symbol):
    with pytest.raises(ValueError) as error:
        normalize_user_sequence(f" a\tU\nG {symbol}C", "invalid")

    message = str(error.value)
    assert symbol in message
    assert "4" in message


@pytest.mark.parametrize("sequence", ["", " \r\n\t", "AU\u200bGC"])
def test_empty_sequences_and_invisible_non_whitespace_are_rejected(sequence):
    with pytest.raises(ValueError):
        normalize_rna_sequence(sequence)


def test_multiline_fasta_preserves_identifier_and_description():
    assert parse_single_fasta("\n>ref.1 example description\r\nat g\r\nc u\n\n") == (
        "ref.1", "example description", "AUGCU"
    )


@pytest.mark.parametrize(
    "text",
    [
        "",
        "ATGC\n",
        ">\nATGC\n",
        "> \t\nATGC\n",
        ">ref\n \t\n",
        ">ref\nATGC\n>other\nATGC\n",
        ">ref\nATGC\n>other\n",
        ">ref\nATNC\n",
    ],
    ids=["empty", "no-header", "no-id", "blank-id", "no-bases", "two-records", "empty-second", "ambiguous"],
)
def test_invalid_fasta_is_never_silently_truncated(text):
    with pytest.raises(ValueError):
        parse_single_fasta(text)


def test_utf8_bom_fasta_with_spaces_in_path_keeps_source_metadata(tmp_path):
    path = tmp_path / "reference with spaces.fa"
    path.write_text(">reference descriptive title\nat\ngc\n", encoding="utf-8-sig")

    record = load_fasta_record(path)
    assert record.identifier == "reference"
    assert record.description == "descriptive title"
    assert record.sequence == "AUGC"
    assert record.metadata["source"] == "fasta"
    assert record.metadata["fasta_identifier"] == "reference"
    assert Path(record.metadata["path"]) == path

    renamed = load_fasta_record(path, label="my_sample")
    assert renamed.identifier == "my_sample"
    assert renamed.metadata["fasta_identifier"] == "reference"


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"accession": "NC_TEST", "sequence": "AUGC"},
        {"accession": "NC_TEST", "fasta_path": "reference.fa"},
        {"sequence": "AUGC", "fasta_path": "reference.fa"},
        {"sequence": "AUGC", "mutant_sequence": "AUGC", "mutant_fasta_path": "mutant.fa"},
        {"sequence": "AUGC", "mutant_sequence": "", "mutant_fasta_path": "mutant.fa"},
    ],
)
def test_python_api_requires_unambiguous_sources(kwargs):
    with pytest.raises(ValueError):
        PipelineInput(**kwargs)


def test_existing_positional_pipeline_input_remains_compatible(pipeline):
    result = pipeline.run(PipelineInput(None, "ATGC", "legacy_sample", "AUGC"))

    assert result.sequence_record.identifier == "legacy_sample"
    assert result.variant_result.total_differences == 0


@pytest.mark.parametrize("mutant", [None, "", " \n", "AUGNC"])
def test_missing_mutant_differs_from_invalid_mutant(pipeline, monkeypatch, mutant):
    if mutant is None:
        result = pipeline.run(PipelineInput(sequence="ATGC"))
        assert result.variant_result is None
        assert result.sequence_record.identifier == "custom_sequence"
        return

    predict = Mock(side_effect=AssertionError("Prediction must wait for valid inputs"))
    monkeypatch.setattr("deepstructgenomics.pipeline.predict_secondary_structure", predict)
    with pytest.raises(ValueError):
        pipeline.run_and_export(PipelineInput(sequence="ATGC", mutant_sequence=mutant))

    predict.assert_not_called()
    assert not list(pipeline.config.default_output_dir.rglob("*"))


@pytest.mark.parametrize("invalid_source", ["raw-reference", "fasta-reference", "fasta-mutant", "missing-mutant"])
def test_invalid_inputs_fail_before_prediction_or_report_creation(pipeline, tmp_path, monkeypatch, invalid_source):
    bad_fasta = tmp_path / "invalid.fa"
    bad_fasta.write_text(">sample\nATGN\n", encoding="utf-8")
    kwargs = {"sequence": "ATGC"}
    if invalid_source == "raw-reference":
        kwargs = {"sequence": "AUGNC"}
    elif invalid_source == "fasta-reference":
        kwargs = {"fasta_path": bad_fasta}
    elif invalid_source == "fasta-mutant":
        kwargs["mutant_fasta_path"] = bad_fasta
    else:
        kwargs["mutant_fasta_path"] = tmp_path / "missing.fa"
    predict = Mock(side_effect=AssertionError("Prediction must wait for valid inputs"))
    monkeypatch.setattr("deepstructgenomics.pipeline.predict_secondary_structure", predict)

    with pytest.raises((ValueError, OSError)):
        pipeline.run_and_export(PipelineInput(**kwargs))

    predict.assert_not_called()
    assert not list(pipeline.config.default_output_dir.rglob("*"))


def test_ncbi_dna_and_equivalent_rna_have_no_false_substitution(pipeline, monkeypatch):
    def respond(session, method, url, **kwargs):
        response = requests.Response()
        response.status_code = 200
        response.url = url
        if url.endswith("efetch.fcgi"):
            response._content = b">NC_TEST.1 DNA reference\nATGC\n"
        elif url.endswith("esummary.fcgi"):
            response._content = json.dumps({"result": {"uids": ["1"], "1": {"title": "DNA reference"}}}).encode()
        else:
            pytest.fail(f"Unexpected NCBI endpoint: {url}")
        return response

    monkeypatch.setattr(requests.Session, "request", respond)
    result = pipeline.run_and_export(PipelineInput(accession="NC_TEST.1", mutant_sequence="AUGC"))

    assert result.sequence_record.sequence == result.mutant_sequence == "AUGC"
    assert result.structure.sequence == "AUGC"
    assert result.annotations["au_ratio"] == 0.5
    assert result.variant_result.total_differences == 0
    payload = json.loads(result.report_paths.json_path.read_text(encoding="utf-8"))
    assert payload["sequence"] == "AUGC"
    assert payload["variant_analysis"]["substitutions"] == []
