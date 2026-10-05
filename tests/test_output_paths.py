"""FASTA identifiers stay metadata, never output directory instructions."""

import json

import pytest

from deepstructgenomics.config import PipelineConfig
from deepstructgenomics.pipeline import DeepStructPipeline, PipelineInput


@pytest.mark.parametrize("identifier", ["../outside", "C:\\outside\\sample", "NUL", "CON.txt", "..."])
def test_fasta_identifiers_preserve_metadata_and_stay_in_output_dir(tmp_path, identifier):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(f">{identifier} Synthetic sequence\nGAAAAC\n", encoding="utf-8")
    output_dir = tmp_path / "reports"
    pipeline = DeepStructPipeline(PipelineConfig(cache_dir=tmp_path / "cache", default_output_dir=output_dir))

    result = pipeline.run_and_export(PipelineInput(fasta_path=fasta))

    payload = json.loads(result.report_paths.json_path.read_text(encoding="utf-8"))
    assert payload["identifier"] == identifier
    assert payload["annotations"]["metadata"]["fasta_identifier"] == identifier
    assert result.report_paths.json_path.parent == output_dir
    assert result.report_paths.markdown_path.parent == output_dir
    assert result.visualization_paths.root_dir.parent.parent == output_dir
    assert result.visualization_paths.root_dir.parent.name == result.report_paths.json_path.stem
    assert result.visualization_paths.wt_structure.is_file()
    assert not (tmp_path / "outside.json").exists()
