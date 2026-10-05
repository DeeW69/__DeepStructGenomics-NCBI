"""Tests for export helper utilities."""

from __future__ import annotations

import pytest

from deepstructgenomics.visualization.export_helpers import compute_camera_angles, parse_export_view


def test_parse_export_view_normalizes_and_validates():
    assert parse_export_view("ISO") == "iso"
    assert parse_export_view(None) == "auto"
    with pytest.raises(ValueError):
        parse_export_view("diagonal")


def test_compute_camera_angles_returns_expected_pairs():
    assert compute_camera_angles("auto") == (0.0, 0.0)
    assert compute_camera_angles("top") == (0.0, -90.0)
    assert compute_camera_angles("side") == (90.0, 0.0)
