"""Pure helpers to compute impact summaries (no heavy dependencies)."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Dict, List, Mapping, Optional, Sequence, Tuple


def compute_hotspots(
    delta_scores: Mapping[int, float],
    top_k: int = 10,
    min_abs_delta: float = 0.1,
) -> List[Dict[str, float]]:
    """Return top-k positions sorted by |delta| descending."""

    if not delta_scores or top_k <= 0:
        return []
    threshold = max(0.0, float(min_abs_delta))
    filtered: List[Tuple[int, float]] = []
    for position, value in delta_scores.items():
        value_f = float(value)
        if abs(value_f) >= threshold:
            filtered.append((int(position), value_f))
    filtered.sort(key=lambda item: abs(item[1]), reverse=True)
    limit = min(len(filtered), int(top_k))
    return [
        {"position": pos, "delta": value, "abs_delta": abs(value)}
        for pos, value in filtered[:limit]
    ]


def summarize_delta(delta_scores: Mapping[int, float]) -> Dict[str, float]:
    """Return descriptive statistics for a delta-score dictionary."""

    if not delta_scores:
        return {
            "min": 0.0,
            "max": 0.0,
            "mean": 0.0,
            "median": 0.0,
            "nonzero_ratio": 0.0,
            "positive_ratio": 0.0,
            "negative_ratio": 0.0,
        }
    values = [float(v) for v in delta_scores.values()]
    sorted_values = sorted(values)
    count = len(values)
    mid = count // 2
    if count % 2 == 0:
        median = (sorted_values[mid - 1] + sorted_values[mid]) / 2.0
    else:
        median = sorted_values[mid]
    eps = 1e-9
    nonzero_ratio = sum(1 for value in values if abs(value) > eps) / count
    positive_ratio = sum(1 for value in values if value > eps) / count
    negative_ratio = sum(1 for value in values if value < -eps) / count
    return {
        "min": float(min(values)),
        "max": float(max(values)),
        "mean": float(sum(values) / count),
        "median": float(median),
        "nonzero_ratio": float(nonzero_ratio),
        "positive_ratio": float(positive_ratio),
        "negative_ratio": float(negative_ratio),
    }


def base_pair_hotspots(
    pairs: Sequence[Tuple[int, int]],
    delta_scores: Mapping[int, float],
    threshold: float,
) -> Dict[str, float]:
    """Return counts of base pairs impacted by |delta| >= threshold."""

    total = len(pairs)
    if total == 0:
        return {"total_pairs": 0, "affected_pairs": 0, "ratio": 0.0, "threshold": float(threshold)}
    cutoff = max(0.0, float(threshold))
    affected = 0
    for raw_i, raw_j in pairs:
        i = int(raw_i)
        j = int(raw_j)
        delta_i = abs(float(delta_scores.get(i, 0.0)))
        delta_j = abs(float(delta_scores.get(j, 0.0)))
        if max(delta_i, delta_j) >= cutoff:
            affected += 1
    ratio = affected / total if total else 0.0
    return {
        "total_pairs": total,
        "affected_pairs": affected,
        "ratio": float(ratio),
        "threshold": float(cutoff),
    }


def build_impact_summary(
    identifier: str,
    generated_at: datetime,
    delta_scores: Mapping[int, float],
    *,
    top_k: int = 10,
    min_abs_delta: float = 0.1,
    note: Optional[str] = None,
    base_pairs: Optional[Mapping[str, Sequence[Tuple[int, int]]]] = None,
    base_pair_threshold: float = 0.2,
) -> Dict[str, object]:
    """Aggregate statistics, hotspots and optional base-pair insights."""

    stats = summarize_delta(delta_scores)
    hotspots = compute_hotspots(delta_scores, top_k=top_k, min_abs_delta=min_abs_delta)
    timestamp = _format_timestamp(generated_at)
    summary: Dict[str, object] = {
        "identifier": identifier,
        "generated_at": timestamp,
        "delta_statistics": stats,
        "hotspots": hotspots,
        "parameters": {
            "top_k": int(top_k),
            "min_abs_delta": float(min_abs_delta),
        },
    }
    if note:
        summary["note"] = note
    if base_pairs:
        base_summary: Dict[str, Dict[str, float]] = {}
        for label, pair_list in base_pairs.items():
            if not pair_list:
                continue
            base_summary[label] = base_pair_hotspots(pair_list, delta_scores, base_pair_threshold)
        if base_summary:
            summary["base_pair_summary"] = {
                "threshold": float(base_pair_threshold),
                "sets": base_summary,
            }
    return summary


def _format_timestamp(moment: datetime) -> str:
    """Return an ISO8601 timestamp suffixed with Z (UTC)."""

    if moment.tzinfo is None:
        utc_dt = moment
    else:
        utc_dt = moment.astimezone(timezone.utc).replace(tzinfo=None)
    return utc_dt.isoformat(timespec="seconds") + "Z"
