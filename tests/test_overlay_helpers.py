"""Tests for overlay helper utilities."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from deepstructgenomics.visualization.io_structures import AtomRecord, MolecularStructure, ResidueKey
from deepstructgenomics.visualization.overlay_helpers import (
    build_base_pair_segments,
    build_position_mapping,
    compute_base_pair_delta,
    compute_displacements,
    compute_link_scalars,
    filter_base_pairs_by_delta,
    filter_links_by_distance,
    filter_links_by_threshold,
    load_base_pairs_json,
    normalize_base_pairs,
    summarize_displacements,
)
from deepstructgenomics.visualization.tk_vtk_overlay import build_status_text, resolve_point_metadata


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


def test_build_status_text_contains_expected_lines():
    text = build_status_text(
        mode_label="overlay",
        score_files=["mut_scores.json"],
        scalar_range=(0.1, 0.9),
        clamp_range=(0.0, 1.0),
        legend_line="WT=gray transparent | MUT=colored score",
        extra_lines=("Extra line",),
    )
    assert "Mode: overlay | Scores: mut_scores.json" in text
    assert "Scalar range: 0.100..0.900 | Clamp: [0.000,1.000]" in text
    assert "Legend: WT=gray transparent | MUT=colored score" in text
    assert "Extra line" in text
    last_line = text.splitlines()[-1]
    assert "Keys: R reset" in last_line
    assert last_line.endswith("P base pairs")


def test_build_position_mapping_prefers_sequence_index():
    wt_indices = [1, 2, None, 4]
    mut_indices = [1, None, 3, 4]
    mapping = build_position_mapping(wt_indices, mut_indices)
    # seq indices present in both -> mapping uses those
    assert mapping[1] == (0, 0)
    assert mapping[4] == (3, 3)
    # fallback on 1-based index when seq index missing on one side
    assert mapping[3] == (2, 2)


def test_compute_link_scalars_overlay_and_delta():
    mut_scores = {1: 0.2, 2: 1.2}
    overlay_values = compute_link_scalars(
        "overlay",
        mut_scores=mut_scores,
        clamp_min=0.0,
        clamp_max=1.0,
    )
    assert overlay_values[1] == 0.2
    # clamped to 1.0
    assert overlay_values[2] == 1.0

    wt_scores = {1: 0.1, 2: -0.2}
    mut_scores_delta = {1: 0.5, 2: 0.2}
    delta_values = compute_link_scalars(
        "overlay-delta",
        wt_scores=wt_scores,
        mut_scores=mut_scores_delta,
        clamp_min=-1.0,
        clamp_max=1.0,
    )
    # (0.5 - 0.1) = 0.4
    assert delta_values[1] == 0.4
    # (0.2 - -0.2) = 0.4
    assert delta_values[2] == 0.4


def test_filter_links_by_threshold_applies_absolute_value():
    values = {1: 0.1, 2: -0.3, 3: 0.25}
    assert filter_links_by_threshold(values, None) == {1, 2, 3}
    assert filter_links_by_threshold(values, 0.2) == {2, 3}


def test_load_base_pairs_normalization_dedup(tmp_path: Path):
    raw_pairs = [[5, 2], [3, 3], ["7", "9"], [-1, 4], [9, 7], [2, 5]]
    normalized = normalize_base_pairs(raw_pairs)
    assert normalized == [(2, 5), (7, 9)]

    data = {"base_pairs": raw_pairs}
    path = tmp_path / "pairs.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    loaded = load_base_pairs_json(path)
    assert loaded == [(2, 5), (7, 9)]

    delta = compute_base_pair_delta({(1, 2), (2, 5)}, {(2, 5), (4, 6)})
    assert delta["shared"] == {(2, 5)}
    assert delta["lost"] == {(1, 2)}
    assert delta["gained"] == {(4, 6)}


def test_build_base_pair_segments_bounds():
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    pairs = [(1, 2), (2, 4)]  # second pair out of range
    segments = build_base_pair_segments(coords, pairs)
    assert segments.shape == (1, 2, 3)
    assert tuple(segments[0, 0]) == (0.0, 0.0, 0.0)
    assert tuple(segments[0, 1]) == (1.0, 0.0, 0.0)


def test_filter_base_pairs_by_delta_threshold():
    pairs = [(1, 2), (2, 3)]
    scores_wt = {1: 0.1, 2: 0.2, 3: 0.3}
    scores_mut = {1: 0.18, 2: 0.21, 3: 0.65}

    # threshold <= 0 keeps everything
    assert filter_base_pairs_by_delta(scores_wt, scores_mut, pairs, 0.0) == pairs

    filtered = filter_base_pairs_by_delta(scores_wt, scores_mut, pairs, 0.2)
    # Only (2,3) keeps max(|mut-wt|) >= 0.2
    assert filtered == [(2, 3)]

    # Missing score maps -> no filtering
    assert filter_base_pairs_by_delta(None, scores_mut, pairs, 0.3) == pairs
    assert filter_base_pairs_by_delta(scores_wt, None, pairs, 0.3) == pairs

    # Missing per-position entry => conservatively keep
    partial_wt = {1: 0.1}
    partial_mut = {1: 0.3}
    assert filter_base_pairs_by_delta(partial_wt, partial_mut, pairs, 0.3) == pairs


def test_compute_displacements_and_summary_and_distance_filter():
    wt_coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ]
    )
    mut_coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 3.0, 0.0],
        ]
    )
    mapping = {1: (0, 0), 2: (1, 1), 3: (2, 2)}
    distances = compute_displacements(wt_coords, mut_coords, mapping)
    assert distances[1] == 0.0
    assert distances[2] == 1.0
    assert distances[3] == 3.0

    summary = summarize_displacements(distances)
    assert summary["count"] == 3
    assert summary["max"] == 3.0
    assert summary["median"] == 1.0
    assert summary["p95"] >= 2.8

    filtered = filter_links_by_distance(mapping, distances, min_distance=0.5, max_distance=2.5)
    assert list(filtered.keys()) == [2]
