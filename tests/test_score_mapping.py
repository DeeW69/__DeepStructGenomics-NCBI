"""Tests for score mapping utilities."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from deepstructgenomics.visualization.io_structures import AtomRecord, MolecularStructure, ResidueKey
from deepstructgenomics.visualization.score_mapping import (
    ScoreTable,
    load_score_table,
    map_scores_to_atoms,
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
