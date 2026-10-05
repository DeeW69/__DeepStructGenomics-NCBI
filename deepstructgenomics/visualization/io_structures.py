"""Structure I/O helpers for visualization modules."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser


@dataclass(frozen=True)
class ResidueKey:
    """Identifies a residue within a structure."""

    chain_id: str
    resseq: int
    insertion_code: str = ""

    def as_compact(self) -> str:
        """Return a compact string representation chain:resseq:icode."""

        code = self.insertion_code.strip() or ""
        chain = self.chain_id.strip() or "_"
        return f"{chain}:{self.resseq}:{code}"


@dataclass(frozen=True)
class AtomRecord:
    """Minimal atom representation required for visualization."""

    serial_number: int
    name: str
    element: Optional[str]
    residue_name: str
    residue_key: ResidueKey
    coord: np.ndarray
    b_factor: float
    occupancy: float
    sequence_index: Optional[int] = None


@dataclass
class MolecularStructure:
    """In-memory representation of a structure file."""

    atoms: List[AtomRecord]
    source_path: Path
    format: str
    metadata: Dict[str, str]

    def as_numpy(self) -> np.ndarray:
        """Return Nx3 coordinates array."""

        if not self.atoms:
            return np.zeros((0, 3), dtype=float)
        return np.vstack([atom.coord for atom in self.atoms])


def load_structure(path: str | Path, model_index: int = 0) -> MolecularStructure:
    """Load a PDB or mmCIF file into MolecularStructure."""

    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(file_path)
    suffix = file_path.suffix.lower()
    if suffix not in {".pdb", ".cif", ".mmcif"}:
        raise ValueError(f"Format non supporte pour {file_path.name}.")

    if suffix == ".pdb":
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(file_path.stem, str(file_path))
        fmt = "pdb"
    else:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(file_path.stem, str(file_path))
        fmt = "mmcif"

    models = list(structure.get_models())
    if not models:
        raise ValueError(f"Aucun modele trouve dans {file_path}.")
    if model_index >= len(models):
        raise IndexError(f"Modele {model_index} indisponible (max {len(models) - 1}).")
    model = models[model_index]

    atoms: List[AtomRecord] = []
    for chain in model:
        chain_id = (chain.id or "_").strip() or "_"
        for residue in chain:
            hetflag, resseq, icode = residue.id
            if hetflag.strip().upper() == "W":  # skip waters
                continue
            residue_key = ResidueKey(chain_id=chain_id, resseq=int(resseq), insertion_code=icode.strip())
            for atom in residue:
                coord = np.array(atom.coord, dtype=float)
                atoms.append(
                    AtomRecord(
                        serial_number=int(atom.serial_number),
                        name=atom.name.strip(),
                        element=(atom.element.strip() if atom.element else None),
                        residue_name=residue.resname.strip(),
                        residue_key=residue_key,
                        coord=coord,
                        b_factor=float(atom.bfactor),
                        occupancy=float(atom.occupancy),
                        sequence_index=None,
                    )
                )

    metadata = {
        "model_index": str(model_index),
        "atom_count": str(len(atoms)),
    }
    return MolecularStructure(atoms=atoms, source_path=file_path, format=fmt, metadata=metadata)


def generate_coarse_backbone(sequence: str, radius: float = 12.0, rise_per_base: float = 2.8) -> np.ndarray:
    """Generate a coarse 3D backbone as a helix projection."""

    seq = sequence.upper().replace("T", "U")
    n = len(seq)
    if n == 0:
        raise ValueError("Sequence vide pour la generation 3D.")

    coords = np.zeros((n, 3), dtype=float)
    twist = math.pi / 3.2
    for idx in range(n):
        angle = idx * twist
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        z = idx * rise_per_base
        coords[idx] = (x, y, z)
    return coords


def export_sequence_as_pseudo_pdb(
    sequence: str,
    output_path: str | Path,
    chain_id: str = "A",
    residue_name: str = "NTP",
) -> Path:
    """Export a coarse-grained PDB with one pseudo-atom per nucleotide."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    coords = generate_coarse_backbone(sequence)
    lines: List[str] = []
    for idx, coord in enumerate(coords, start=1):
        atom_name = "P"
        resseq = idx
        lines.append(
            "ATOM  {serial:5d} {name:^4} {resname:>3} {chain:1s}{resseq:4d}    "
            "{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfactor:6.2f}          {element:>2s}".format(
                serial=idx,
                name=atom_name,
                resname=residue_name,
                chain=chain_id,
                resseq=resseq,
                x=coord[0],
                y=coord[1],
                z=coord[2],
                occ=1.00,
                bfactor=0.00,
                element="P",
            )
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_visualization_manifest(
    path: str | Path,
    *,
    identifier: str,
    source: str,
    parameters: Dict[str, object],
    wt_structure: Optional[Path],
    mutant_structure: Optional[Path],
    wt_score_file: Optional[Path],
    mutant_score_file: Optional[Path],
    impact_summary_file: Optional[Path] = None,
) -> Path:
    """Write the manifest describing visualization artifacts."""

    manifest = {
        "identifier": identifier,
        "source": source,
        "wt_structure": str(wt_structure) if wt_structure else None,
        "mutant_structure": str(mutant_structure) if mutant_structure else None,
        "wt_score_file": str(wt_score_file) if wt_score_file else None,
        "mutant_score_file": str(mutant_score_file) if mutant_score_file else None,
        "impact_summary_file": str(impact_summary_file) if impact_summary_file else None,
        "parameters": parameters,
    }
    manifest_path = Path(path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest_path


def build_residue_atom_index(structure: MolecularStructure) -> Dict[str, List[int]]:
    """Return mapping residue_key -> indices of atoms."""

    mapping: Dict[str, List[int]] = {}
    for idx, atom in enumerate(structure.atoms):
        key = atom.residue_key.as_compact()
        mapping.setdefault(key, []).append(idx)
    return mapping


def assign_sequence_indices(
    structure: MolecularStructure,
    sequence_indices: Sequence[int],
) -> MolecularStructure:
    """Return a new MolecularStructure with per-atom sequence indices."""

    if len(sequence_indices) != len(structure.atoms):
        raise ValueError("Nombre d'indices sequence != nombre d'atomes.")
    atoms = [
        AtomRecord(
            serial_number=atom.serial_number,
            name=atom.name,
            element=atom.element,
            residue_name=atom.residue_name,
            residue_key=atom.residue_key,
            coord=atom.coord.copy(),
            b_factor=atom.b_factor,
            occupancy=atom.occupancy,
            sequence_index=sequence_indices[idx],
        )
        for idx, atom in enumerate(structure.atoms)
    ]
    return MolecularStructure(atoms=atoms, source_path=structure.source_path, format=structure.format, metadata=structure.metadata)
