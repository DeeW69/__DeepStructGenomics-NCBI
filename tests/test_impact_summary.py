"""Tests for impact summary helpers."""

from __future__ import annotations

from datetime import datetime

from deepstructgenomics.analysis.impact_summary import (
    base_pair_hotspots,
    build_impact_summary,
    compute_hotspots,
    summarize_delta,
)


def test_compute_hotspots_sorting_and_threshold():
    scores = {1: 0.05, 2: -0.5, 3: 0.21, 4: -0.18}
    hotspots = compute_hotspots(scores, top_k=3, min_abs_delta=0.15)
    assert [entry["position"] for entry in hotspots] == [2, 3, 4]
    assert hotspots[0]["abs_delta"] == 0.5
    assert hotspots[-1]["delta"] == -0.18


def test_summarize_delta_sign_ratios():
    scores = {1: -0.2, 2: 0.0, 3: 0.3, 4: 0.45}
    summary = summarize_delta(scores)
    assert summary["min"] == -0.2
    assert summary["max"] == 0.45
    assert summary["median"] == 0.15
    assert summary["positive_ratio"] == 0.5
    assert summary["negative_ratio"] == 0.25
    assert summary["nonzero_ratio"] == 0.75


def test_no_mutant_summary_shape():
    generated_at = datetime(2024, 1, 1, 12, 0, 0)
    summary = build_impact_summary(
        "run-1",
        generated_at,
        {},
        note="mutant not provided",
        base_pairs={"wt": [(1, 2), (3, 4)]},
        base_pair_threshold=0.3,
    )
    assert summary["identifier"] == "run-1"
    assert summary["hotspots"] == []
    assert summary["delta_statistics"]["mean"] == 0.0
    assert summary["note"] == "mutant not provided"
    base_info = summary["base_pair_summary"]
    assert base_info["threshold"] == 0.3
    assert base_info["sets"]["wt"]["total_pairs"] == 2
    assert base_info["sets"]["wt"]["affected_pairs"] == 0
