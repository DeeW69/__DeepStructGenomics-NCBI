"""Utilities for export presets (camera view, validation)."""

from __future__ import annotations

from typing import Tuple

VALID_EXPORT_VIEWS = {"auto", "iso", "top", "side", "front"}
CAMERA_PRESETS = {
    "iso": (45.0, 35.0),
    "top": (0.0, -90.0),
    "side": (90.0, 0.0),
    "front": (0.0, 0.0),
}


def parse_export_view(value: str | None) -> str:
    """Normalize export view preset."""

    if not value:
        return "auto"
    preset = value.strip().lower()
    if preset not in VALID_EXPORT_VIEWS:
        raise ValueError(f"Export view '{value}' invalide (attendu: {sorted(VALID_EXPORT_VIEWS)}).")
    return preset


def compute_camera_angles(preset: str) -> Tuple[float, float]:
    """Return (azimuth, elevation) delta for a preset."""

    preset = parse_export_view(preset)
    if preset == "auto":
        return (0.0, 0.0)
    return CAMERA_PRESETS[preset]
