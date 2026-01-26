"""Tests for visualization IO helpers."""

from __future__ import annotations

from pathlib import Path

from deepstructgenomics.visualization.io_structures import load_structure


def test_load_structure_from_pdb(tmp_path):
    pdb_content = (
        "ATOM      1  P   NTP A   1      0.000   0.000   0.000  1.00  0.00          P\n"
        "ATOM      2  P   NTP A   2      5.000   0.000   0.000  1.00  0.00          P\n"
        "END\n"
    )
    pdb_path = tmp_path / "mini.pdb"
    pdb_path.write_text(pdb_content, encoding="utf-8")
    structure = load_structure(pdb_path)
    assert structure.format == "pdb"
    assert len(structure.atoms) == 2
    assert structure.atoms[0].residue_key.as_compact() == "A:1:"


def test_load_structure_from_mmcif(tmp_path):
    mmcif_content = """data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
ATOM 1 P P . NTP A A 1 1 ? 0.000 0.000 0.000 1.00 0.00 1
ATOM 2 P P . NTP A A 1 2 ? 1.000 0.000 0.000 1.00 0.00 1
"""
    cif_path = tmp_path / "mini.cif"
    cif_path.write_text(mmcif_content, encoding="utf-8")
    structure = load_structure(cif_path)
    assert structure.format == "mmcif"
    assert len(structure.atoms) == 2
    assert structure.atoms[1].residue_key.as_compact().startswith("A:")
