"""Score transformation utilities for visualization."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional

import numpy as np

from deepstructgenomics.rna.secondary_structure import SecondaryStructureResult

from .io_structures import MolecularStructure


DEBUG = os.getenv("DSG_VIZ_DEBUG") == "1"


def _clamp(value: float) -> float:
    return max(min(value, 1.0), 0.0)


def _clamp_range(value: float, minimum: float, maximum: float) -> float:
    if minimum > maximum:
        minimum, maximum = maximum, minimum
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
            position_scores[idx] = float(value)

        raw_residue_scores = payload.get("residue_scores")
        if raw_residue_scores is None:
            raw_residue_scores = {}
        elif not isinstance(raw_residue_scores, dict):
            raise ValueError("residue_scores doit etre un dictionnaire.")
        residue_scores = {str(key): float(value) for key, value in raw_residue_scores.items()}

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


def map_scores_to_atoms(
    structure: MolecularStructure,
    scores: ScoreTable,
    *,
    clamp_min: float = 0.0,
    clamp_max: float = 1.0,
) -> np.ndarray:
    """Map residue/position scores to atom-level scalar array."""

    scalars = np.zeros(len(structure.atoms), dtype=float)
    for idx, atom in enumerate(structure.atoms):
        residue_key = atom.residue_key.as_compact()
        value: Optional[float] = scores.residue_scores.get(residue_key)
        if value is None and atom.sequence_index is not None:
            value = scores.position_scores.get(int(atom.sequence_index))
        if value is None:
            value = 0.0
        scalars[idx] = _clamp_range(value, clamp_min, clamp_max)
    if DEBUG and len(scalars):
        nonzero = int((np.abs(scalars) > 1e-9).sum())
        pct = (nonzero / len(scalars)) * 100.0
        print(f"[map_scores_to_atoms] atoms={len(scalars)} nonzero={nonzero} ({pct:.1f}%)")
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


def compute_score_statistics(scores: ScoreTable) -> Dict[str, float]:
    """Compute descriptive stats for manifest logging."""

    values = list(scores.position_scores.values()) or list(scores.residue_scores.values())
    if not values:
        return {
            "value_count": 0,
            "scalar_min": 0.0,
            "scalar_max": 0.0,
            "nonzero_ratio": 0.0,
            "positive_ratio": 0.0,
            "negative_ratio": 0.0,
        }
    n_values = len(values)
    scalar_min = float(min(values))
    scalar_max = float(max(values))
    nonzero = sum(1 for val in values if abs(val) > 1e-9)
    positive = sum(1 for val in values if val > 0)
    negative = sum(1 for val in values if val < 0)
    return {
        "value_count": n_values,
        "scalar_min": scalar_min,
        "scalar_max": scalar_max,
        "nonzero_ratio": nonzero / n_values,
        "positive_ratio": positive / n_values,
        "negative_ratio": negative / n_values,
    }


def compute_delta_scores(
    wt_scores: ScoreTable,
    mutant_scores: ScoreTable,
    *,
    clamp_min: float = -1.0,
    clamp_max: float = 1.0,
) -> ScoreTable:
    """Return delta (mutant - wt) score table."""

    delta_positions: Dict[int, float] = {}
    delta_residues: Dict[str, float] = {}

    pos_keys = set(wt_scores.position_scores.keys()) & set(mutant_scores.position_scores.keys())
    if pos_keys:
        for idx in sorted(pos_keys):
            delta = mutant_scores.position_scores[idx] - wt_scores.position_scores[idx]
            delta_positions[idx] = _clamp_range(delta, clamp_min, clamp_max)

    res_keys = set(wt_scores.residue_scores.keys()) & set(mutant_scores.residue_scores.keys())
    if res_keys:
        for key in sorted(res_keys):
            delta = mutant_scores.residue_scores[key] - wt_scores.residue_scores[key]
            delta_residues[key] = _clamp_range(delta, clamp_min, clamp_max)

    if not delta_positions and not delta_residues:
        raise ValueError("Impossible de calculer le delta : aucun recouvrement entre les scores.")

    metadata = {
        "context": "delta",
        "reference_identifier": mutant_scores.metadata.get("reference_identifier")
        or wt_scores.metadata.get("reference_identifier", ""),
        "strategy": "position_first" if delta_positions else "residue_only",
        "position_overlap": len(delta_positions),
        "residue_overlap": len(delta_residues),
        "clamp_min": clamp_min,
        "clamp_max": clamp_max,
    }
    return ScoreTable(residue_scores=delta_residues, position_scores=delta_positions, metadata=metadata, version=wt_scores.version)
