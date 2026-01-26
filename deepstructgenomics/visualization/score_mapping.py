"""Score transformation utilities for visualization."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional

import numpy as np

from deepstructgenomics.rna.secondary_structure import SecondaryStructureResult

from .io_structures import MolecularStructure


def _clamp(value: float, minimum: float = 0.0, maximum: float = 1.0) -> float:
    return max(min(value, maximum), minimum)


@dataclass
class ScoreTable:
    """Normalized scores accessible by residue or sequential index."""

    residue_scores: Dict[str, float] = field(default_factory=dict)
    position_scores: Dict[int, float] = field(default_factory=dict)
    metadata: Dict[str, str] = field(default_factory=dict)
    version: int = 1

    def to_dict(self) -> Dict[str, object]:
        return {
            "version": self.version,
            "metadata": self.metadata,
            "residue_scores": self.residue_scores,
            "position_scores": {str(k): v for k, v in self.position_scores.items()},
        }

    @classmethod
    def from_dict(cls, payload: Dict[str, object]) -> "ScoreTable":
        version = int(payload.get("version", 1))
        metadata = {str(k): str(v) for k, v in (payload.get("metadata") or {}).items()}

        raw_positions = payload.get("position_scores")
        if raw_positions is None:
            raw_positions = {}
        elif not isinstance(raw_positions, dict):
            raise ValueError("position_scores doit etre un dictionnaire.")
        position_scores: Dict[int, float] = {}
        for key, value in raw_positions.items():
            idx = int(key)
            position_scores[idx] = _clamp(float(value))

        raw_residue_scores = payload.get("residue_scores")
        if raw_residue_scores is None:
            raw_residue_scores = {}
        elif not isinstance(raw_residue_scores, dict):
            raise ValueError("residue_scores doit etre un dictionnaire.")
        residue_scores = {str(key): _clamp(float(value)) for key, value in raw_residue_scores.items()}

        return cls(residue_scores=residue_scores, position_scores=position_scores, metadata=metadata, version=version)


def write_score_table(path: str | Path, scores: ScoreTable) -> Path:
    """Serialize a ScoreTable to JSON."""

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(scores.to_dict(), indent=2, ensure_ascii=False), encoding="utf-8")
    return target


def load_score_table(path: str | Path) -> ScoreTable:
    """Load and validate a JSON score payload."""

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return ScoreTable.from_dict(payload)


def derive_position_scores_from_structure(structure: SecondaryStructureResult) -> Dict[int, float]:
    """Derive heuristic per-position scores [0,1] from a secondary structure."""

    sequence = structure.sequence
    pair_lookup: Dict[int, int] = {}
    for i, j in structure.base_pairs:
        pair_lookup[i] = j
        pair_lookup[j] = i

    scores: Dict[int, float] = {}
    for idx, base in enumerate(sequence):
        base_class = structure.dot_bracket[idx]
        if base_class == "(" or base_class == ")":
            score = 0.8
        else:
            score = 0.3
        if base in {"G", "C"}:
            score += 0.1
        if idx in pair_lookup:
            distance = abs(pair_lookup[idx] - idx)
            score += min(distance / len(sequence), 0.2)
        scores[idx + 1] = _clamp(score)
    return scores


def map_scores_to_atoms(structure: MolecularStructure, scores: ScoreTable) -> np.ndarray:
    """Map residue/position scores to atom-level scalar array."""

    scalars = np.zeros(len(structure.atoms), dtype=float)
    for idx, atom in enumerate(structure.atoms):
        residue_key = atom.residue_key.as_compact()
        value: Optional[float] = scores.residue_scores.get(residue_key)
        if value is None and atom.sequence_index is not None:
            value = scores.position_scores.get(int(atom.sequence_index))
        if value is None:
            value = 0.0
        scalars[idx] = _clamp(value)
    return scalars


def merge_residue_scores(
    structure: MolecularStructure,
    position_scores: Dict[int, float],
    default_chain: str = "A",
) -> Dict[str, float]:
    """Create residue scores from sequential indices when explicit mappings are absent."""

    residue_scores: Dict[str, float] = {}
    for atom in structure.atoms:
        seq_idx = atom.sequence_index
        if seq_idx is None:
            continue
        if seq_idx not in position_scores:
            continue
        key = atom.residue_key.as_compact()
        residue_scores.setdefault(key, position_scores[seq_idx])
    if not residue_scores and position_scores:
        # fallback for coarse models with implicit residues
        for idx, value in position_scores.items():
            key = f"{default_chain}:{idx}:"
            residue_scores[key] = value
    return residue_scores
