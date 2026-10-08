"""Offline integration coverage for reports and visualization artifacts."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import Mock

import pytest
import requests
from vtkmodules.vtkRenderingCore import vtkRenderer

from deepstructgenomics.config import PipelineConfig
from deepstructgenomics.pipeline import DeepStructPipeline, PipelineInput
from deepstructgenomics.visualization.io_structures import load_structure
from deepstructgenomics.visualization.score_mapping import (
    compute_delta_scores,
    load_score_table,
    map_scores_to_atoms,
)
from deepstructgenomics.visualization.tk_vtk_overlay import BasePairOptions, OverlayDeltaViewer, OverlayViewer


REFERENCE = "GAAAAC"
MUTANT = "AAAAAC"


@pytest.fixture
def pipeline(tmp_path, monkeypatch):
    def reject_network(*args, **kwargs):
        pytest.fail("User-provided sequences must not make HTTP requests")

    monkeypatch.setattr(requests.Session, "request", reject_network)
    return DeepStructPipeline(
        PipelineConfig(
            cache_dir=tmp_path / "cache",
            default_output_dir=tmp_path / "reports",
        )
    )


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def test_pair_legend_matches_rendered_and_hidden_pairs(pipeline):
    result = pipeline.run_and_export(PipelineInput(sequence=REFERENCE, mutant_sequence=MUTANT))
    paths = result.visualization_paths
    viewer = OverlayDeltaViewer(paths.wt_structure, paths.mutant_structure,
                                paths.wt_score_file, paths.mutant_score_file,
                                base_pair_config=BasePairOptions(enabled=True, wt_pairs=((1, 6),)))
    viewer._build_scene(vtkRenderer())
    assert "BasePairs: on | WT=1 MUT=0" in viewer._status_actor.GetInput()
    for actor in viewer._base_pair_actors:
        actor.SetVisibility(0)
    viewer._refresh_base_pair_status()
    assert "BasePairs: off | WT=1 MUT=0" in viewer._status_actor.GetInput()
    assert viewer._scalar_bar_actor.GetUnconstrainedFontSize()
    lut = viewer._lut_diverging()
    neutral = lut.GetTableValue(lut.GetNumberOfTableValues() // 2)[:3]
    assert min(neutral) > 0.8 and max(neutral) - min(neutral) < 0.1


def test_run_computes_summary_before_export(pipeline):
    result = pipeline.run(
        PipelineInput(sequence=REFERENCE, mutant_sequence=MUTANT, sequence_label="mini")
    )

    assert result.report_paths is None
    assert result.visualization_paths is None
    assert result.impact_summary["identifier"] == "mini"
    assert [entry["position"] for entry in result.impact_summary["hotspots"]] == [1, 6]
    assert not list(pipeline.config.default_output_dir.iterdir())


@pytest.mark.parametrize("mode", ["overlay", "overlay-delta"])
def test_viewer_builds_annotations_and_tooltips_without_a_render_window(pipeline, mode):
    result = pipeline.run_and_export(PipelineInput(sequence=REFERENCE, mutant_sequence=MUTANT))
    artifacts = result.visualization_paths
    if mode == "overlay":
        viewer = OverlayViewer(artifacts.wt_structure, artifacts.mutant_structure, artifacts.mutant_score_file)
    else:
        viewer = OverlayDeltaViewer(
            artifacts.wt_structure, artifacts.mutant_structure,
            artifacts.wt_score_file, artifacts.mutant_score_file,
        )
    renderer = vtkRenderer()
    viewer._build_scene(renderer)
    assert renderer.HasViewProp(viewer._scalar_bar_actor)
    assert renderer.HasViewProp(viewer._status_actor)
    previous_count = renderer.GetViewProps().GetNumberOfItems()
    viewer._configure_tooltip(renderer, Mock(), Mock())
    assert renderer.GetViewProps().GetNumberOfItems() == previous_count + 1


def test_reports_share_summary_and_preserve_existing_artifacts(pipeline):
    result = pipeline.run_and_export(
        PipelineInput(sequence=REFERENCE, mutant_sequence=MUTANT, sequence_label="mini")
    )

    payload = read_json(result.report_paths.json_path)
    artifacts = result.visualization_paths
    summary = payload["impact_summary"]
    assert summary == result.impact_summary == read_json(artifacts.impact_summary_file)
    assert summary["identifier"] == payload["identifier"] == "mini"
    assert summary["generated_at"].endswith("Z")
    assert "note" not in summary

    # A single substitution breaks the terminal pair in this six-base hairpin.
    assert payload["sequence"] == REFERENCE
    assert payload["annotations"]["length"] == 6
    assert payload["structure"]["dot_bracket"] == "(....)"
    assert payload["structure"]["base_pairs"] == [[0, 5]]
    assert payload["variant_analysis"]["total_differences"] == 1
    assert payload["variant_analysis"]["structure_delta_score"] == -1.0
    assert payload["variant_analysis"]["substitutions"] == [
        {"position": 0, "reference": "G", "mutant": "A"}
    ]
    assert result.mutant_structure.dot_bracket == "......"

    hotspots = summary["hotspots"]
    assert [entry["position"] for entry in hotspots] == [1, 6]
    assert [entry["delta"] for entry in hotspots] == pytest.approx([-0.7, -0.6])
    assert [entry["abs_delta"] for entry in hotspots] == pytest.approx([0.7, 0.6])
    assert summary["delta_statistics"]["mean"] == pytest.approx(-1.3 / 6)
    assert summary["delta_statistics"]["negative_ratio"] == pytest.approx(2 / 6)
    assert summary["base_pair_summary"]["sets"]["wt"]["affected_pairs"] == 1

    for path in (
        artifacts.wt_structure,
        artifacts.mutant_structure,
        artifacts.wt_score_file,
        artifacts.mutant_score_file,
    ):
        assert path.is_file()
    manifest = read_json(artifacts.manifest_path)
    assert Path(manifest["impact_summary_file"]) == artifacts.impact_summary_file

    # Re-open the actual PDB/score files just as the viewer does. A shifted PDB
    # chain column would silently map every delta to zero despite correct JSON.
    mutant_atoms = load_structure(artifacts.mutant_structure)
    delta_table = compute_delta_scores(
        load_score_table(artifacts.wt_score_file), load_score_table(artifacts.mutant_score_file)
    )
    rendered_scores = map_scores_to_atoms(mutant_atoms, delta_table, clamp_min=-1.0, clamp_max=1.0)
    assert rendered_scores.tolist() == pytest.approx([-0.7, 0.0, 0.0, 0.0, 0.0, -0.6])

    markdown = result.report_paths.markdown_path.read_text(encoding="utf-8")
    assert "(....)" in markdown
    assert "## Résumé des variations par position" in markdown
    assert "Positions numérotées à partir de 1" in markdown
    assert "| 1 | G | A |" in markdown
    assert "| 1 | -0.700 |" in markdown
    assert "| 6 | -0.600 |" in markdown

    assert result.report_paths.csv_path is not None
    assert result.report_paths.csv_path.is_file()
    csv_text = result.report_paths.csv_path.read_text(encoding="utf-8")
    assert csv_text.splitlines()[0] == "position,reference,mutant,delta,abs_delta"
    assert "1,G,A,-0.7,0.7" in csv_text
    assert "6,C,C,-0.6,0.6" in csv_text
    assert artifacts.hotspots_csv_file is not None
    assert artifacts.hotspots_csv_file.is_file()
    assert Path(manifest["hotspots_csv_file"]) == artifacts.hotspots_csv_file


def test_export_hotspots_csv_empty(tmp_path):
    from deepstructgenomics.reporting.report_generator import export_hotspots_csv
    out_csv = tmp_path / "empty_hotspots.csv"
    export_hotspots_csv([], out_csv)
    assert out_csv.read_text(encoding="utf-8").strip() == "position,reference,mutant,delta,abs_delta"


@pytest.mark.parametrize("mutant_sequence", [None, REFERENCE], ids=["no-mutant", "identical-mutant"])
def test_missing_mutant_is_distinct_from_zero_impact(pipeline, mutant_sequence):
    result = pipeline.run_and_export(
        PipelineInput(sequence=REFERENCE, mutant_sequence=mutant_sequence)
    )

    payload = read_json(result.report_paths.json_path)
    artifacts = result.visualization_paths
    summary = result.impact_summary
    assert summary == payload["impact_summary"] == read_json(artifacts.impact_summary_file)
    assert summary["hotspots"] == []
    assert summary["delta_statistics"]["mean"] == 0.0
    assert summary["delta_statistics"]["nonzero_ratio"] == 0.0
    markdown = result.report_paths.markdown_path.read_text(encoding="utf-8")

    if mutant_sequence is None:
        assert summary["note"] == "mutant not provided"
        assert result.variant_result is None
        assert "variant_analysis" not in payload
        assert artifacts.mutant_structure is None
        assert "Comparaison indisponible : aucune séquence mutante fournie." in markdown
        assert "Aucune position n'atteint le seuil" not in markdown
        # Retain the legacy score-file alias even for reference-only analyses.
        assert artifacts.score_file == artifacts.mutant_score_file
        assert read_json(artifacts.mutant_score_file) == read_json(artifacts.wt_score_file)
    else:
        assert "note" not in summary
        assert payload["variant_analysis"]["total_differences"] == 0
        assert artifacts.mutant_structure.is_file()
        assert "Comparaison indisponible" not in markdown
        assert "Aucune position n'atteint le seuil" in markdown


def test_different_lengths_explain_positional_comparison(pipeline):
    result = pipeline.run_and_export(
        PipelineInput(sequence=REFERENCE, mutant_sequence=REFERENCE + "AAA", align_mutant=False)
    )

    summary = result.impact_summary
    assert summary["note"] == "different sequence lengths; positional comparison without alignment"
    assert summary == read_json(result.report_paths.json_path)["impact_summary"]
    assert summary == read_json(result.visualization_paths.impact_summary_file)
    assert all(1 <= entry["position"] <= len(REFERENCE) for entry in summary["hotspots"])
    markdown = result.report_paths.markdown_path.read_text(encoding="utf-8")
    assert "longueurs différentes" in markdown
    assert "seules les positions communes" in markdown
    assert "sans alignement" in markdown
