"""Pure helpers for overlay link construction (no VTK dependency)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np


def build_position_mapping(
    wt_sequence_indices: Sequence[Optional[int]],
    mut_sequence_indices: Sequence[Optional[int]],
) -> Dict[int, Tuple[int, int]]:
    """Return mapping index -> (wt_point_idx, mut_point_idx).

    Sequence indices have priority; fallback uses 1-based point order.
    """

    def _canonical_lookup(indices: Sequence[Optional[int]]) -> Dict[int, int]:
        lookup: Dict[int, int] = {}
        for idx, seq_idx in enumerate(indices):
            key = int(seq_idx) if seq_idx is not None else (idx + 1)
            if key not in lookup:
                lookup[key] = idx
        return lookup

    wt_lookup = _canonical_lookup(wt_sequence_indices)
    mut_lookup = _canonical_lookup(mut_sequence_indices)

    mapping: Dict[int, Tuple[int, int]] = {}
    for key in sorted(set(wt_lookup.keys()) & set(mut_lookup.keys())):
        mapping[key] = (wt_lookup[key], mut_lookup[key])
    return mapping


def compute_link_scalars(
    mode: str,
    *,
    wt_scores: Optional[Dict[int, float]] = None,
    mut_scores: Optional[Dict[int, float]] = None,
    clamp_min: float,
    clamp_max: float,
) -> Dict[int, float]:
    """Return per-position scalar values to color WT->MUT links."""

    if mut_scores is None:
        mut_scores = {}
    lo, hi = (clamp_min, clamp_max) if clamp_min <= clamp_max else (clamp_max, clamp_min)

    def _clamp(value: float) -> float:
        return max(min(value, hi), lo)

    if mode == "overlay":
        base = mut_scores
    elif mode == "overlay-delta":
        base = {}
        if wt_scores is not None:
            common = set(mut_scores.keys()) & set(wt_scores.keys())
            for key in common:
                base[key] = mut_scores[key] - wt_scores[key]
        else:
            base = mut_scores
    else:
        raise ValueError(f"Unsupported mode: {mode}")

    return {key: _clamp(value) for key, value in base.items()}


def filter_links_by_threshold(values: Dict[int, float], threshold: Optional[float]) -> Set[int]:
    """Return set of positions that satisfy |value| >= threshold (if provided)."""

    if threshold is None or threshold <= 0:
        return set(values.keys())
    minimum = abs(threshold)
    return {key for key, value in values.items() if abs(value) >= minimum}


def compute_displacements(
    wt_coords: np.ndarray,
    mut_coords: np.ndarray,
    mapping: Dict[int, Tuple[int, int]],
) -> Dict[int, float]:
    """Return per-position Euclidean distances between WT and mutant."""

    distances: Dict[int, float] = {}
    if wt_coords.ndim != 2 or mut_coords.ndim != 2:
        return distances
    n_wt = wt_coords.shape[0]
    n_mut = mut_coords.shape[0]
    for position, (wt_idx, mut_idx) in mapping.items():
        if wt_idx < 0 or wt_idx >= n_wt or mut_idx < 0 or mut_idx >= n_mut:
            continue
        wt_point = wt_coords[wt_idx]
        mut_point = mut_coords[mut_idx]
        if np.isnan(wt_point).any() or np.isnan(mut_point).any():
            continue
        distances[position] = float(np.linalg.norm(mut_point - wt_point))
    return distances


def summarize_displacements(distances: Dict[int, float]) -> Dict[str, float]:
    """Return descriptive statistics for displacement values."""

    if not distances:
        return {"count": 0.0, "mean": 0.0, "median": 0.0, "p95": 0.0, "max": 0.0, "min": 0.0}
    values = np.array(list(distances.values()), dtype=float)
    return {
        "count": float(len(values)),
        "mean": float(values.mean()),
        "median": float(np.median(values)),
        "p95": float(np.percentile(values, 95)),
        "max": float(values.max()),
        "min": float(values.min()),
    }


def filter_links_by_distance(
    mapping: Dict[int, Tuple[int, int]],
    distances: Dict[int, float],
    min_distance: float = 0.0,
    max_distance: Optional[float] = None,
) -> Dict[int, Tuple[int, int]]:
    """Return mapping filtered based on min/max distance thresholds."""

    if not mapping:
        return {}
    min_distance = max(0.0, float(min_distance))
    max_distance = float(max_distance) if max_distance is not None else None
    filtered: Dict[int, Tuple[int, int]] = {}
    for position, indices in mapping.items():
        value = distances.get(position)
        if value is None:
            continue
        if value < min_distance:
            continue
        if max_distance is not None and value > max_distance:
            continue
        filtered[position] = indices
    return filtered


def normalize_base_pairs(pairs_iterable: Iterable[Iterable[int]]) -> List[Tuple[int, int]]:
    """Normalize iterable of base-pair indices to sorted tuples with i < j."""

    normalized: Set[Tuple[int, int]] = set()
    for pair in pairs_iterable:
        if pair is None:
            continue
        if len(pair) != 2:
            continue
        i, j = pair
        try:
            a = int(i)
            b = int(j)
        except (TypeError, ValueError):
            continue
        if a == b:
            continue
        if a > b:
            a, b = b, a
        if a <= 0 or b <= 0:
            continue
        normalized.add((a, b))
    return sorted(normalized)


def load_base_pairs_json(path: str | Path) -> List[Tuple[int, int]]:
    """Load base-pair pairs stored as {"base_pairs": [[i,j], ...]}."""

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    raw_pairs = payload.get("base_pairs")
    if raw_pairs is None:
        raise ValueError("Le fichier JSON doit contenir la cle 'base_pairs'.")
    if not isinstance(raw_pairs, list):
        raise ValueError("'base_pairs' doit etre une liste de paires.")
    return normalize_base_pairs(raw_pairs)


def compute_base_pair_delta(
    wt_pairs: Set[Tuple[int, int]],
    mut_pairs: Set[Tuple[int, int]],
) -> Dict[str, Set[Tuple[int, int]]]:
    """Return dictionary with gained, lost, and shared pair sets."""

    shared = wt_pairs & mut_pairs
    gained = mut_pairs - wt_pairs
    lost = wt_pairs - mut_pairs
    return {"shared": shared, "gained": gained, "lost": lost}


def base_pair_delta_value(
    scores_wt: Optional[Dict[int, float]],
    scores_mut: Optional[Dict[int, float]],
    pair: Tuple[int, int],
) -> Optional[float]:
    """Return max(|mut-wt|) on the pair indices, or None if unavailable."""

    if scores_wt is None or scores_mut is None:
        return None
    i, j = pair
    deltas: List[float] = []
    for position in (int(i), int(j)):
        wt_value = scores_wt.get(position)
        mut_value = scores_mut.get(position)
        if wt_value is None or mut_value is None:
            return None
        deltas.append(abs(mut_value - wt_value))
    return max(deltas) if deltas else None


def filter_base_pairs_by_delta(
    scores_wt: Optional[Dict[int, float]],
    scores_mut: Optional[Dict[int, float]],
    pairs: Sequence[Tuple[int, int]],
    threshold: float,
) -> List[Tuple[int, int]]:
    """Filter base pairs based on max(|mut - wt|) >= threshold.

    If threshold <= 0 or score dictionaries are missing, the original list is returned.
    Missing per-position entries fall back to a conservative keep (no filtering).
    """

    if threshold is None or threshold <= 0:
        return list(pairs)
    if not scores_wt or not scores_mut:
        return list(pairs)
    filtered: List[Tuple[int, int]] = []
    for pair in pairs:
        delta = base_pair_delta_value(scores_wt, scores_mut, pair)
        if delta is None:
            filtered.append(pair)
            continue
        if delta >= threshold:
            filtered.append(pair)
    return filtered


def build_base_pair_segments(
    coords: np.ndarray,
    pairs: Sequence[Tuple[int, int]],
) -> np.ndarray:
    """Return array of point-to-point segments for given base pairs.

    Output shape: (M, 2, 3) with 1-based indices referencing coords[i - 1].
    """

    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        return np.zeros((0, 2, 3), dtype=float)
    n_points = coords.shape[0]
    if n_points == 0 or not pairs:
        return np.zeros((0, 2, 3), dtype=float)
    segments: List[np.ndarray] = []
    for i, j in pairs:
        start_idx = int(i) - 1
        end_idx = int(j) - 1
        if start_idx < 0 or end_idx < 0 or start_idx >= n_points or end_idx >= n_points:
            continue
        p1 = coords[start_idx]
        p2 = coords[end_idx]
        if np.isnan(p1).any() or np.isnan(p2).any():
            continue
        segments.append(np.stack((p1, p2), axis=0))
    if not segments:
        return np.zeros((0, 2, 3), dtype=float)
    return np.stack(segments, axis=0).astype(float, copy=False)
