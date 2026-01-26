"""Tests for overlay helper utilities."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from deepstructgenomics.visualization.io_structures import AtomRecord, MolecularStructure, ResidueKey
from deepstructgenomics.visualization.tk_vtk_overlay import resolve_point_metadata


def _make_structure() -> MolecularStructure:
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
            sequence_index=3,
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
            sequence_index=None,
        ),
    ]
    return MolecularStructure(atoms=atoms, source_path=Path("dummy.pdb"), format="pdb", metadata={})


def test_resolve_point_metadata_uses_sequence_index():
    structure = _make_structure()
    assert resolve_point_metadata(structure, 0) == (3, "A:1:")
    # fallback to point index + 1 when sequence_index missing
    assert resolve_point_metadata(structure, 1) == (2, "A:2:")
    assert resolve_point_metadata(structure, 99) is None
