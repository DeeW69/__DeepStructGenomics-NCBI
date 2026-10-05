"""Tests for score mapping utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import numpy as np
import pytest

from deepstructgenomics.visualization.io_structures import AtomRecord, MolecularStructure, ResidueKey
from deepstructgenomics.visualization.score_mapping import (
    ScoreTable,
    compute_delta_scores,
    compute_score_statistics,
    detect_hotspots,
    load_score_table,
    map_scores_to_atoms,
    summarize_deltas,
    write_score_table,
)


def _dummy_structure(tmp_path: Path) -> MolecularStructure:
    atoms = [
        AtomRecord(
            serial_number=1,
            name="P",
            element="P",
            residue_name="NTP",
            residue_key=ResidueKey("A", 1, ""),
            coord=np.array([0.0, 0.0, 0.0]),
            b_factor=0.0,
            occupancy=1.0,
            sequence_index=1,
        ),
        AtomRecord(
            serial_number=2,
            name="P",
            element="P",
            residue_name="NTP",
            residue_key=ResidueKey("A", 2, ""),
            coord=np.array([1.0, 0.0, 0.0]),
            b_factor=0.0,
            occupancy=1.0,
            sequence_index=2,
        ),
    ]
    return MolecularStructure(atoms=atoms, source_path=tmp_path / "dummy.pdb", format="pdb", metadata={})


def test_map_scores_prefers_residue_keys(tmp_path):
    structure = _dummy_structure(tmp_path)
    scores = ScoreTable(residue_scores={"A:2:": 0.9}, position_scores={1: 0.4}, metadata={})
    scalars = map_scores_to_atoms(structure, scores)
    assert pytest.approx(scalars[0]) == 0.4
    assert pytest.approx(scalars[1]) == 0.9


def test_score_table_validation(tmp_path):
    payload_path = tmp_path / "scores.json"
    scores = ScoreTable(residue_scores={"A:1:": 0.5}, position_scores={1: 0.5}, metadata={"context": "unit-test"})
    write_score_table(payload_path, scores)
    loaded = load_score_table(payload_path)
    assert loaded.residue_scores["A:1:"] == pytest.approx(0.5)

    invalid_path = tmp_path / "invalid.json"
    invalid_path.write_text('{"position_scores": []}', encoding="utf-8")
    with pytest.raises(ValueError):
        load_score_table(invalid_path)


def test_compute_delta_scores_with_positions():
    wt = ScoreTable(position_scores={1: 0.2, 2: 0.4}, metadata={"reference_identifier": "ref"})
    mut = ScoreTable(position_scores={1: 0.6, 2: -0.5}, metadata={"reference_identifier": "ref"})
    delta = compute_delta_scores(wt, mut)
    assert delta.position_scores[1] == pytest.approx(0.4)
    # clamped to [-1,1]
    assert delta.position_scores[2] == pytest.approx(-0.9)
    assert delta.metadata["position_overlap"] == 2
    assert delta.metadata["strategy"] == "position_first"
    assert delta.metadata["clamp_min"] == pytest.approx(-1.0)
    assert delta.metadata["clamp_max"] == pytest.approx(1.0)


def test_map_scores_supports_negative_range(tmp_path):
    structure = _dummy_structure(tmp_path)
    scores = ScoreTable(position_scores={1: -0.6, 2: 0.8}, metadata={})
    scalars = map_scores_to_atoms(structure, scores, clamp_min=-1.0, clamp_max=1.0)
    assert pytest.approx(scalars[0]) == -0.6
    assert pytest.approx(scalars[1]) == 0.8


def test_compute_score_statistics_handles_signs():
    scores = ScoreTable(position_scores={1: -0.5, 2: 0.0, 3: 0.7}, metadata={})
    stats = compute_score_statistics(scores)
    assert stats["scalar_min"] == pytest.approx(-0.5)
    assert stats["scalar_max"] == pytest.approx(0.7)
    assert stats["nonzero_ratio"] == pytest.approx(2 / 3)
    assert stats["positive_ratio"] == pytest.approx(1 / 3)
    assert stats["negative_ratio"] == pytest.approx(1 / 3)


def _delta_table(values: Dict[int, float]) -> ScoreTable:
    return ScoreTable(position_scores=values, metadata={"context": "delta"})


def test_detect_hotspots_threshold_and_min_run():
    table = _delta_table({1: 0.4, 2: 0.5, 3: 0.1, 4: -0.35, 5: -0.45, 10: 0.9})
    hotspots = detect_hotspots(table, abs_threshold=0.3, min_run=2)
    assert len(hotspots) == 2
    assert hotspots[0]["start"] == 1 and hotspots[0]["end"] == 2
    assert pytest.approx(hotspots[0]["mean_abs_delta"], rel=1e-3) == 0.45
    assert hotspots[1]["start"] == 4 and hotspots[1]["end"] == 5


def test_summarize_deltas_includes_top_positions_and_hotspots():
    table = _delta_table({1: -0.6, 2: 0.05, 3: 0.9, 4: 0.2})
    summary = summarize_deltas(table, abs_threshold=0.2, top_k=2)
    assert summary["max_abs_delta"] == pytest.approx(0.9)
    assert summary["pct_positions_over_threshold"] == pytest.approx(0.75)
    assert len(summary["top_positions"]) == 2
    assert summary["hotspots"] == []
