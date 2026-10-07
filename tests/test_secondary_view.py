"""Verify the diagram represents exported predictions, including missing mutants."""

import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import pytest
import requests

from deepstructgenomics.config import PipelineConfig
from deepstructgenomics.pipeline import DeepStructPipeline, PipelineInput
from deepstructgenomics.visualization.secondary_view import (
    LOST, build_secondary_figure, load_secondary_data, validate_secondary_data,
)


@pytest.fixture
def pipeline(tmp_path, monkeypatch):
    def reject(*args, **kwargs):
        pytest.fail("Secondary diagrams must be entirely offline")
    monkeypatch.setattr(requests.Session, "request", reject)
    return DeepStructPipeline(PipelineConfig(cache_dir=tmp_path / "cache", default_output_dir=tmp_path / "reports"))


def test_exported_diagram_matches_pairs_scores_and_csv(pipeline):
    result = pipeline.run_and_export(PipelineInput(sequence="GGGGAAAACCCC", mutant_sequence="GGGGAAAACCCA"))
    data = load_secondary_data(result.visualization_paths.manifest_path)
    assert data["reference"]["base_pairs"] == [[0, 11], [1, 10], [2, 9], [3, 8]]
    assert data["mutant"]["base_pairs"] == [[1, 10], [2, 9], [3, 8]]
    folder = result.visualization_paths.root_dir
    pairs = json.loads((folder / "wt_base_pairs.json").read_text())
    assert pairs["base_pairs"] == [[1, 12], [2, 11], [3, 10], [4, 9]]
    with result.report_paths.csv_path.open() as stream:
        for row in csv.DictReader(stream):
            idx = int(row["position"]) - 1
            assert data["mutant"]["scores"][idx] - data["reference"]["scores"][idx] == pytest.approx(float(row["delta"]))
    fig = build_secondary_figure(data)
    try:
        assert len(fig.axes[0].patches) == 4
        assert len(fig.axes[1].patches) == 3
        assert sum(p.get_edgecolor() == to_rgba(LOST) for p in fig.axes[0].patches) == 1
        assert [bar.get_height() for bar in fig.axes[2].patches] == pytest.approx([-0.6] + [0] * 10 + [-0.7])
    finally:
        plt.close(fig)


@pytest.mark.parametrize("mutant", [None, "GGGGAAAACCCC", "GGGGAAAACCCCAAA"])
def test_missing_identical_and_unequal_mutants_remain_distinct(pipeline, mutant):
    result = pipeline.run_and_export(PipelineInput(sequence="GGGGAAAACCCC", mutant_sequence=mutant))
    data = load_secondary_data(result.visualization_paths.manifest_path)
    fig = build_secondary_figure(data)
    try:
        if mutant is None:
            assert data["mutant"] is None
            assert not fig.axes[2].patches
            assert any("indisponible" in text.get_text() for text in fig.axes[2].texts)
        else:
            assert len(fig.axes[2].patches) == 12
            if len(mutant) == 12:
                assert all(bar.get_height() == 0 for bar in fig.axes[2].patches)
            else:
                assert any("sans alignement" in text.get_text() for text in fig.texts)
    finally:
        plt.close(fig)


def test_reference_only_regeneration_does_not_keep_mutant_pairs(pipeline):
    pipeline.run_and_export(PipelineInput(sequence="GGGGAAAACCCC", mutant_sequence="GGGGAAAACCCA"))
    result = pipeline.run_and_export(PipelineInput(sequence="GGGGAAAACCCC"))
    assert load_secondary_data(result.visualization_paths.manifest_path)["mutant"] is None
    assert json.loads((result.visualization_paths.root_dir / "mut_base_pairs.json").read_text())["base_pairs"] == []


def test_inconsistent_pairs_are_rejected_instead_of_drawn(pipeline):
    result = pipeline.run_and_export(PipelineInput(sequence="GGGGAAAACCCC"))
    data = load_secondary_data(result.visualization_paths.manifest_path)
    data["reference"]["base_pairs"] = [[0, 10]]
    with pytest.raises(ValueError, match="incoherents"):
        validate_secondary_data(data)


def test_old_manifest_explains_how_to_regenerate(tmp_path):
    manifest = tmp_path / "old.json"
    manifest.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="regenerer"):
        load_secondary_data(manifest)


def test_secondary_cli_exports_from_relocated_bundle(pipeline, tmp_path):
    result = pipeline.run_and_export(PipelineInput(sequence="GGGGAAAACCCC", mutant_sequence="GGGGAAAACCCA"))
    relocated = tmp_path / "relocated"
    shutil.copytree(result.visualization_paths.root_dir, relocated)
    # Force the basename fallback without removing or editing original reports.
    manifest = relocated / "visualization_manifest.json"
    payload = json.loads(manifest.read_text())
    payload["secondary_structure_file"] = "old/location/secondary_structures.json"
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    script = Path(__file__).resolve().parents[1] / "scripts/view_secondary.py"
    output = tmp_path / "figure.png"
    process = subprocess.run([sys.executable, "-B", str(script), "--manifest", str(manifest),
                              "--export", str(output), "--export-only"],
                             cwd=tmp_path, capture_output=True, text=True, timeout=45)
    assert process.returncode == 0, process.stderr
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
