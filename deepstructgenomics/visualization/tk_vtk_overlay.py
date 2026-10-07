"""Overlay viewers highlighting WT vs mutant differences."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import vtkmodules.vtkRenderingOpenGL2  # noqa: F401
from vtkmodules.vtkCommonCore import vtkFloatArray, vtkLookupTable, vtkPoints, vtkUnsignedCharArray
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkLine, vtkPolyData, vtkPolyLine
from vtkmodules.vtkFiltersCore import vtkGlyph3D, vtkTubeFilter
from vtkmodules.vtkFiltersSources import vtkSphereSource
from vtkmodules.vtkRenderingAnnotation import vtkScalarBarActor
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkPointPicker,
    vtkPolyDataMapper,
    vtkRenderer,
    vtkTextActor,
)

from .io_structures import MolecularStructure, load_structure
from .overlay_helpers import (
    base_pair_delta_value,
    build_base_pair_segments,
    build_position_mapping,
    compute_displacements,
    compute_link_scalars,
    filter_base_pairs_by_delta,
    filter_links_by_distance,
    filter_links_by_threshold,
    summarize_displacements,
)
from .score_mapping import (
    ScoreTable,
    compute_delta_scores,
    compute_score_statistics,
    load_score_table,
    map_scores_to_atoms,
)
from .vtk_window import VTKStandaloneWindow

DEBUG = os.getenv("DSG_VIZ_DEBUG") == "1"

BASE_PAIR_COLOR = (0.82, 0.82, 0.88)


@dataclass
class BasePairOptions:
    enabled: bool = False
    mode: str = "lines"
    threshold: float = 0.0
    wt_pairs: Sequence[Tuple[int, int]] = field(default_factory=tuple)
    mut_pairs: Sequence[Tuple[int, int]] = field(default_factory=tuple)


def build_status_text(
    *,
    mode_label: str,
    score_files: Sequence[str],
    scalar_range: Tuple[float, float],
    clamp_range: Tuple[float, float],
    legend_line: str,
    extra_lines: Iterable[str] = (),
) -> str:
    """Return a newline-separated status text for overlay viewers."""

    scalar_min, scalar_max = scalar_range
    clamp_min, clamp_max = clamp_range
    files_text = " + ".join(score_files) if score_files else "n/a"
    lines = [
        f"Mode: {mode_label} | Scores: {files_text}",
        f"Scalar range: {scalar_min:.3f}..{scalar_max:.3f} | Clamp: [{clamp_min:.3f},{clamp_max:.3f}]",
        f"Legend: {legend_line}",
    ]
    for extra in extra_lines:
        if extra:
            lines.append(extra)
    lines.append("Keys: R reset | W WT | M MUT | B backbone | L links | P base pairs")
    return "\n".join(lines)


def resolve_point_metadata(structure: MolecularStructure, point_id: int) -> Optional[Tuple[int, str]]:
    """Return (position_index, residue_key) for a glyph point."""

    if point_id < 0 or point_id >= len(structure.atoms):
        return None
    atom = structure.atoms[point_id]
    position = atom.sequence_index if atom.sequence_index is not None else (point_id + 1)
    return position, atom.residue_key.as_compact()


class OverlayBase:
    """Common helpers shared by overlay viewers."""

    def __init__(
        self,
        *,
        title: str,
        width: int,
        height: int,
        show_backbone: bool = True,
        enable_tooltip: Optional[bool] = None,
        tooltip_verbose: Optional[bool] = None,
        enable_links: bool = False,
        link_threshold: Optional[float] = None,
        base_pair_config: Optional[BasePairOptions] = None,
        link_metric: str = "score",
        link_min_distance: Optional[float] = None,
        link_max_distance: Optional[float] = None,
    ) -> None:
        self.window = VTKStandaloneWindow(title=title, width=width, height=height)
        self.show_backbone = show_backbone
        tooltip_enabled = DEBUG if enable_tooltip is None else enable_tooltip
        verbose_flag = DEBUG if tooltip_verbose is None else tooltip_verbose
        self.enable_tooltip = tooltip_enabled
        if tooltip_enabled:
            self._tooltip_mode = "full" if verbose_flag else "minimal"
        else:
            self._tooltip_mode = "off"
        self.enable_links = bool(enable_links)
        self.link_threshold = link_threshold
        self.base_pair_config = base_pair_config
        self.link_metric = (link_metric or "score").lower()
        self.link_min_distance = max(0.0, float(link_min_distance or 0.0))
        if link_max_distance is not None and link_max_distance >= 0:
            self.link_max_distance: Optional[float] = max(float(link_max_distance), self.link_min_distance)
        else:
            self.link_max_distance = None
        self._tooltip_payload: Optional[Dict[str, object]] = None
        self._interaction_actor: Optional[vtkActor] = None
        self._wt_actor: Optional[vtkActor] = None
        self._mut_actor: Optional[vtkActor] = None
        self._backbone_actor: Optional[vtkActor] = None
        self._links_actor: Optional[vtkActor] = None
        self._base_pair_actors: List[Optional[vtkActor]] = []
        self._base_pair_summary: Optional[Dict[str, object]] = None
        self._link_distance_lookup: Dict[int, float] = {}
        self._link_distance_stats: Optional[Dict[str, float]] = None
        self._scalar_bar_actor: Optional[vtkScalarBarActor] = None
        self._status_actor: Optional[vtkTextActor] = None

    def start(
        self,
        *,
        export_path: Optional[str | Path] = None,
        export_only: bool = False,
        export_scale: int = 1,
        export_view: str = "auto",
        export_hide_ui: bool = False,
        export_background: str = "dark",
    ) -> None:
        self.window.run(
            self._build_scene,
            post_setup=self._post_setup,
            export_path=export_path,
            export_only=export_only,
            export_scale=export_scale,
            export_view=export_view,
            export_hide_ui=export_hide_ui,
            export_background=export_background,
            ui_elements_provider=self._ui_actors,
        )

    def _build_scene(self, renderer: vtkRenderer) -> None:  # pragma: no cover - abstract
        raise NotImplementedError

    def _post_setup(self, renderer: vtkRenderer, render_window, interactor) -> None:
        self._configure_key_bindings(renderer, render_window, interactor)
        if self.enable_tooltip:
            self._configure_tooltip(renderer, render_window, interactor)

    def _configure_tooltip(self, renderer: vtkRenderer, render_window, interactor) -> None:
        payload = self._tooltip_payload
        target_actor = self._interaction_actor
        if not payload or target_actor is None:
            return
        picker = vtkPointPicker()
        picker.PickFromListOn()
        picker.SetTolerance(0.01)
        picker.AddPickList(target_actor)
        interactor.SetPicker(picker)

        text = vtkTextActor()
        text.VisibilityOff()
        textprop = text.GetTextProperty()
        textprop.SetFontSize(16)
        textprop.SetColor(0.95, 0.95, 0.95)
        renderer.AddViewProp(text)

        last_message: Optional[str] = None

        def on_mouse_move(obj, event) -> None:
            nonlocal last_message
            x, y = interactor.GetEventPosition()
            picker.Pick(x, y, 0, renderer)
            point_id = picker.GetPointId()
            if point_id < 0:
                if text.GetVisibility():
                    text.VisibilityOff()
                    render_window.Render()
                last_message = None
                return
            message = self._format_tooltip_text(point_id)
            if not message:
                if text.GetVisibility():
                    text.VisibilityOff()
                    render_window.Render()
                last_message = None
                return
            text.SetInput(message)
            text.SetPosition(x + 15, y + 15)
            if not text.GetVisibility():
                text.VisibilityOn()
            if last_message != message:
                render_window.Render()
                last_message = message

        interactor.AddObserver("MouseMoveEvent", on_mouse_move)

    def _configure_key_bindings(self, renderer: vtkRenderer, render_window, interactor) -> None:
        def on_key_press(obj, event) -> None:
            key = interactor.GetKeySym()
            if not key:
                return
            key = key.lower()
            if key == "r":
                renderer.ResetCamera()
                render_window.Render()
            elif key == "w":
                self._toggle_actor(self._wt_actor, render_window)
            elif key == "m":
                self._toggle_actor(self._mut_actor, render_window)
            elif key == "b":
                self._toggle_actor(self._backbone_actor, render_window)
            elif key == "l":
                self._toggle_actor(self._links_actor, render_window)
            elif key == "p":
                self._toggle_actor_list(self._base_pair_actors, render_window)
                self._refresh_base_pair_status()
                render_window.Render()

        interactor.AddObserver("KeyPressEvent", on_key_press)

    @staticmethod
    def _toggle_actor(actor: Optional[vtkActor], render_window) -> None:
        if actor is None:
            return
        actor.SetVisibility(0 if actor.GetVisibility() else 1)
        render_window.Render()

    @staticmethod
    def _toggle_actor_list(actors: Sequence[Optional[vtkActor]], render_window) -> None:
        toggled = False
        for actor in actors:
            if actor is None:
                continue
            actor.SetVisibility(0 if actor.GetVisibility() else 1)
            toggled = True
        if toggled:
            render_window.Render()

    @staticmethod
    def _sequence_indices(structure: MolecularStructure) -> Sequence[Optional[int]]:
        return [atom.sequence_index for atom in structure.atoms]

    @staticmethod
    def _scalars_by_position(structure: MolecularStructure, scalars: np.ndarray) -> Dict[int, float]:
        lookup: Dict[int, float] = {}
        count = min(len(structure.atoms), len(scalars))
        for idx in range(count):
            meta = resolve_point_metadata(structure, idx)
            if not meta:
                continue
            position, _ = meta
            lookup[position] = float(scalars[idx])
        return lookup

    def _ui_actors(self) -> Sequence[object]:
        actors: List[object] = []
        if self._scalar_bar_actor:
            actors.append(self._scalar_bar_actor)
        if self._status_actor:
            actors.append(self._status_actor)
        return actors

    def _base_pair_status_line(self) -> str:
        config = self.base_pair_config
        if not config or not getattr(config, "enabled", False):
            return ""
        summary = self._base_pair_summary or {}
        state = "on" if any(actor.GetVisibility() for actor in self._base_pair_actors if actor) else "off"
        wt_count = int(summary.get("wt", 0))
        mut_count = int(summary.get("mut", 0))
        mode_label = str(summary.get("mode", config.mode or "lines"))
        threshold = float(summary.get("threshold", config.threshold or 0.0))
        return (
            f"BasePairs: {state} | WT={wt_count} MUT={mut_count} | "
            f"mode={mode_label} | threshold={threshold:.2f} | toggle=P"
        )

    def _refresh_base_pair_status(self) -> None:
        if self._status_actor is None:
            return
        lines = self._status_actor.GetInput().splitlines()
        line = self._base_pair_status_line()
        lines = [line if text.startswith("BasePairs:") else text for text in lines]
        if line and not any(text.startswith("BasePairs:") for text in lines):
            lines.append(line)
        self._status_actor.SetInput("\n".join(lines))

    def _build_links_actor(
        self,
        wt_structure: MolecularStructure,
        mut_structure: MolecularStructure,
        mapping: Dict[int, Tuple[int, int]],
        positions: Iterable[int],
        scalar_values: Dict[int, float],
        lut: vtkLookupTable,
        scalar_range: Tuple[float, float],
    ) -> Optional[vtkActor]:
        point_pairs: list[Tuple[int, int, float]] = []
        for position in positions:
            if position not in mapping or position not in scalar_values:
                continue
            wt_idx, mut_idx = mapping[position]
            if wt_idx >= len(wt_structure.atoms) or mut_idx >= len(mut_structure.atoms):
                continue
            value = scalar_values[position]
            point_pairs.append((wt_idx, mut_idx, value))
        if not point_pairs:
            return None

        points = vtkPoints()
        cells = vtkCellArray()
        scalars_arr = vtkFloatArray()
        scalars_arr.SetName("link_scalars")

        for wt_idx, mut_idx, value in point_pairs:
            wt_coord = wt_structure.atoms[wt_idx].coord
            mut_coord = mut_structure.atoms[mut_idx].coord
            start_id = points.InsertNextPoint(float(wt_coord[0]), float(wt_coord[1]), float(wt_coord[2]))
            end_id = points.InsertNextPoint(float(mut_coord[0]), float(mut_coord[1]), float(mut_coord[2]))
            line = vtkLine()
            line.GetPointIds().SetId(0, start_id)
            line.GetPointIds().SetId(1, end_id)
            cells.InsertNextCell(line)
            scalars_arr.InsertNextValue(float(value))
            scalars_arr.InsertNextValue(float(value))

        poly = vtkPolyData()
        poly.SetPoints(points)
        poly.SetLines(cells)
        poly.GetPointData().SetScalars(scalars_arr)

        mapper = vtkPolyDataMapper()
        mapper.SetInputData(poly)
        mapper.SetScalarRange(*scalar_range)
        mapper.SetLookupTable(lut)
        mapper.ScalarVisibilityOn()

        actor = vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetLineWidth(2.0)
        return actor

    def _try_create_links(
        self,
        renderer: vtkRenderer,
        wt_structure: MolecularStructure,
        mut_structure: MolecularStructure,
        scalar_values: Dict[int, float],
        scalar_range: Tuple[float, float],
        lut: vtkLookupTable,
        mode: str,
        *,
        wt_scalar_values: Optional[Dict[int, float]] = None,
    ) -> Optional[vtkActor]:
        self._link_distance_lookup = {}
        self._link_distance_stats = None
        if not self.enable_links:
            return None

        mapping = build_position_mapping(self._sequence_indices(wt_structure), self._sequence_indices(mut_structure))
        if not mapping:
            if DEBUG:
                print("[links] skipped (no shared positions between WT and mutant)")
            return None

        need_distance = self.link_metric == "distance" or self.link_min_distance > 0.0 or self.link_max_distance is not None
        if need_distance:
            wt_coords = self._coords_array_from_structure(wt_structure)
            mut_coords = self._coords_array_from_structure(mut_structure)
            distance_lookup = compute_displacements(wt_coords, mut_coords, mapping)
            self._link_distance_lookup = distance_lookup
        else:
            distance_lookup = {}

        if self.link_min_distance > 0.0 or self.link_max_distance is not None:
            mapping = filter_links_by_distance(mapping, distance_lookup, self.link_min_distance, self.link_max_distance)
            if not mapping:
                if DEBUG:
                    print("[links] skipped (distance filtering removed all positions)")
                return None

        active_lut = lut
        active_range = scalar_range
        if self.link_metric == "distance":
            active_values = {pos: distance_lookup[pos] for pos in mapping if pos in distance_lookup}
            if not active_values:
                if DEBUG:
                    print("[links] skipped (distance metric lacks values)")
                return None
            min_val = min(active_values.values())
            max_val = max(active_values.values())
            if abs(max_val - min_val) < 1e-6:
                max_val = min_val + 1e-3
            active_range = (min_val, max_val)
            active_lut = self._lut_distance(max_val)
            positions = sorted(active_values.keys())
            self._link_distance_stats = summarize_displacements(active_values)
        else:
            values = compute_link_scalars(
                mode,
                wt_scores=wt_scalar_values,
                mut_scores=scalar_values,
                clamp_min=scalar_range[0],
                clamp_max=scalar_range[1],
            )
            allowed_positions = filter_links_by_threshold(values, self.link_threshold)
            active_values = {pos: values[pos] for pos in mapping if pos in allowed_positions}
            if not active_values:
                if DEBUG:
                    print("[links] skipped (threshold or data removed all positions)")
                return None
            positions = sorted(active_values.keys())

        actor = self._build_links_actor(
            wt_structure,
            mut_structure,
            mapping,
            positions,
            active_values,
            active_lut,
            active_range,
        )
        if actor:
            renderer.AddActor(actor)
            if self.link_metric == "distance":
                self._append_distance_status_line()
        else:
            if DEBUG:
                print("[links] no actor created after filtering shared positions")
        return actor

    def _append_distance_status_line(self) -> None:
        if not self.enable_links or self.link_metric != "distance":
            return
        stats = self._link_distance_stats
        if not stats or not self._status_actor:
            return
        line = "links metric: distance (min={min_val:.2f}, max={max_val:.2f})".format(
            min_val=stats.get("min", 0.0),
            max_val=stats.get("max", 0.0),
        )
        current = self._status_actor.GetInput() or ""
        if line in current:
            return
        sep = "\n" if current else ""
        self._status_actor.SetInput(current + sep + line)

    @staticmethod
    def _coords_array_from_structure(structure: MolecularStructure) -> np.ndarray:
        positions: List[int] = []
        for idx in range(len(structure.atoms)):
            meta = resolve_point_metadata(structure, idx)
            if not meta:
                continue
            positions.append(meta[0])
        if not positions:
            return np.zeros((0, 3), dtype=float)
        max_pos = max(positions)
        coords = np.full((max_pos, 3), np.nan, dtype=float)
        for idx in range(len(structure.atoms)):
            meta = resolve_point_metadata(structure, idx)
            if not meta:
                continue
            pos, _ = meta
            atom = structure.atoms[idx]
            coords[pos - 1] = np.array(atom.coord, dtype=float)
        return coords

    def _setup_base_pairs(
        self,
        renderer: vtkRenderer,
        *,
        wt_structure: Optional[MolecularStructure],
        mut_structure: Optional[MolecularStructure],
        wt_scores: Optional[Dict[int, float]] = None,
        mut_scores: Optional[Dict[int, float]] = None,
        lut: Optional[vtkLookupTable] = None,
        scalar_range: Optional[Tuple[float, float]] = None,
    ) -> None:
        self._base_pair_actors = []
        self._base_pair_summary = None
        config = self.base_pair_config
        if not config or not config.enabled:
            return

        mode_label = (config.mode or "lines").lower()
        style = "arc" if mode_label == "arcs" else "line"
        threshold = max(0.0, float(config.threshold or 0.0))
        summary = {"wt": 0, "mut": 0, "mode": "arcs" if style == "arc" else "lines", "threshold": threshold}
        should_color = (
            threshold > 0
            and wt_scores
            and mut_scores
            and lut is not None
            and scalar_range is not None
        )

        def render_dataset(
            structure: Optional[MolecularStructure],
            pairs: Sequence[Tuple[int, int]],
            summary_key: str,
        ) -> None:
            if not structure or not pairs:
                return
            coords = self._coords_array_from_structure(structure)
            if not len(coords):
                return
            filtered_pairs = filter_base_pairs_by_delta(wt_scores, mut_scores, pairs, threshold)
            segments = build_base_pair_segments(coords, filtered_pairs)
            if segments.shape[0] == 0:
                return
            scalar_values: Optional[List[float]] = None
            if should_color:
                scalar_values = self._pair_delta_values(filtered_pairs, wt_scores, mut_scores)
            actor = self._create_base_pair_actor(
                segments,
                style=style,
                color=BASE_PAIR_COLOR,
                lut=lut if scalar_values is not None else None,
                scalar_range=scalar_range if scalar_values is not None else None,
                scalar_values=scalar_values,
            )
            if actor:
                renderer.AddActor(actor)
                self._base_pair_actors.append(actor)
                summary[summary_key] = int(segments.shape[0])

        render_dataset(wt_structure, list(config.wt_pairs or ()), "wt")
        render_dataset(mut_structure, list(config.mut_pairs or ()), "mut")
        self._base_pair_summary = summary
        self._refresh_base_pair_status()

    @staticmethod
    def _paths_from_segments(
        segments: np.ndarray,
        style: str,
    ) -> List[List[Tuple[float, float, float]]]:
        if segments.size == 0:
            return []
        if style == "arc":
            return [OverlayBase._arc_path(tuple(seg[0]), tuple(seg[1])) for seg in segments]
        return [[tuple(seg[0]), tuple(seg[1])] for seg in segments]

    @staticmethod
    def _pair_delta_values(
        pairs: Sequence[Tuple[int, int]],
        scores_wt: Optional[Dict[int, float]],
        scores_mut: Optional[Dict[int, float]],
    ) -> Optional[List[float]]:
        if not pairs or not scores_wt or not scores_mut:
            return None
        values: List[float] = []
        for pair in pairs:
            delta = base_pair_delta_value(scores_wt, scores_mut, pair)
            if delta is None:
                return None
            values.append(delta)
        return values

    @staticmethod
    def _arc_path(
        start: Tuple[float, float, float],
        end: Tuple[float, float, float],
    ) -> List[Tuple[float, float, float]]:
        start_v = np.array(start, dtype=float)
        end_v = np.array(end, dtype=float)
        chord = end_v - start_v
        chord_norm = np.linalg.norm(chord)
        if chord_norm < 1e-6:
            return [tuple(start), tuple(end)]
        midpoint = (start_v + end_v) / 2.0
        up = np.array([0.0, 0.0, 1.0])
        offset = np.cross(chord, up)
        offset_norm = np.linalg.norm(offset)
        if offset_norm < 1e-6:
            up = np.array([0.0, 1.0, 0.0])
            offset = np.cross(chord, up)
            offset_norm = np.linalg.norm(offset)
        if offset_norm < 1e-6:
            return [tuple(start), tuple(end)]
        offset_dir = offset / offset_norm
        height = min(20.0, chord_norm * 0.25)
        control = midpoint + offset_dir * height
        return [tuple(start_v), tuple(control), tuple(end_v)]

    def _create_base_pair_actor(
        self,
        segments: np.ndarray,
        *,
        style: str,
        color: Tuple[float, float, float],
        lut: Optional[vtkLookupTable] = None,
        scalar_range: Optional[Tuple[float, float]] = None,
        scalar_values: Optional[Sequence[float]] = None,
    ) -> Optional[vtkActor]:
        if segments.size == 0:
            return None
        paths = self._paths_from_segments(segments, style)
        if not paths:
            return None
        points = vtkPoints()
        cells = vtkCellArray()
        for path in paths:
            if len(path) < 2:
                continue
            polyline = vtkPolyLine()
            polyline.GetPointIds().SetNumberOfIds(len(path))
            for idx, coord in enumerate(path):
                pid = points.InsertNextPoint(float(coord[0]), float(coord[1]), float(coord[2]))
                polyline.GetPointIds().SetId(idx, pid)
            cells.InsertNextCell(polyline)
        if cells.GetNumberOfCells() == 0:
            return None
        poly = vtkPolyData()
        poly.SetPoints(points)
        poly.SetLines(cells)
        mapper = vtkPolyDataMapper()
        mapper.SetInputData(poly)
        actor = vtkActor()
        actor.SetMapper(mapper)

        use_scalars = (
            scalar_values is not None
            and lut is not None
            and scalar_range is not None
            and len(scalar_values) == cells.GetNumberOfCells()
        )
        if use_scalars:
            scalars = vtkFloatArray()
            scalars.SetName("base_pair_delta")
            for value in scalar_values:
                scalars.InsertNextValue(float(value))
            poly.GetCellData().SetScalars(scalars)
            mapper.SetScalarRange(*scalar_range)
            mapper.SetLookupTable(lut)
            mapper.SetScalarModeToUseCellData()
            mapper.ScalarVisibilityOn()
            actor.GetProperty().SetOpacity(0.45)
        else:
            mapper.ScalarVisibilityOff()
            actor.GetProperty().SetColor(*color)
            actor.GetProperty().SetOpacity(0.35)
        actor.GetProperty().SetLineWidth(1.1)
        return actor

    def _format_tooltip_text(self, point_id: int) -> Optional[str]:
        payload = self._tooltip_payload
        if not payload:
            return None
        meta = resolve_point_metadata(payload["structure"], point_id)
        if not meta:
            return None
        pos, residue = meta
        scalars = payload["scalars"]
        if point_id >= len(scalars):
            return None
        label = str(payload.get("label", "score"))
        value = scalars[point_id]
        lines = [
            f"pos={pos}",
            f"residue={residue}",
            f"{label}={value:.4f}",
        ]
        if getattr(self, "_tooltip_mode", "minimal") == "full":
            wt_val = OverlayBase._score_at_position(payload.get("wt_scores"), pos)
            mut_val = OverlayBase._score_at_position(payload.get("mut_scores"), pos)
            label_key = label.lower()
            if label_key == "delta":
                wt_text = f"{wt_val:.4f}" if wt_val is not None else "NA"
                mut_text = f"{mut_val:.4f}" if mut_val is not None else "NA"
                lines.append(f"WT={wt_text} MUT={mut_text}")
            else:
                if mut_val is not None:
                    lines.append(f"mut_raw={mut_val:.4f}")
                if wt_val is not None:
                    lines.append(f"wt_raw={wt_val:.4f}")
        clamp = payload.get("clamp")
        if clamp:
            lines.append(f"clamp=[{clamp[0]:.3f},{clamp[1]:.3f}]")
        scalar_range = payload.get("scalar_range")
        if scalar_range:
            lines.append(f"scalar=[{scalar_range[0]:.3f},{scalar_range[1]:.3f}]")
        if self.link_metric == "distance":
            distance_lookup = payload.get("distance_lookup")
            if isinstance(distance_lookup, dict):
                value = distance_lookup.get(pos)
                if value is not None:
                    lines.append(f"distance={value:.3f}")
        return "\n".join(lines)

    @staticmethod
    def _score_at_position(table: Optional[ScoreTable], pos: int) -> Optional[float]:
        if not table:
            return None
        if pos in table.position_scores:
            return table.position_scores[pos]
        key = f"A:{pos}:"
        return table.residue_scores.get(key)

    @staticmethod
    def _points_poly(structure: MolecularStructure) -> vtkPolyData:
        coords = structure.as_numpy()
        points = vtkPoints()
        for x, y, z in coords:
            points.InsertNextPoint(float(x), float(y), float(z))
        poly = vtkPolyData()
        poly.SetPoints(points)
        return poly

    @staticmethod
    def _solid_color(poly: vtkPolyData, color: Tuple[int, int, int]) -> None:
        colors = vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)
        colors.SetName("solid_rgb")
        for _ in range(poly.GetNumberOfPoints()):
            colors.InsertNextTypedTuple(color)
        poly.GetPointData().SetScalars(colors)

    @staticmethod
    def _lut_blue_to_red() -> vtkLookupTable:
        lut = vtkLookupTable()
        lut.SetNumberOfTableValues(256)
        lut.SetRange(0.0, 1.0)
        lut.Build()
        for i in range(256):
            t = i / 255.0
            r = t
            g = 0.15 * (1.0 - t)
            b = 1.0 - t
            lut.SetTableValue(i, r, g, b, 1.0)
        return lut

    @staticmethod
    def _lut_diverging() -> vtkLookupTable:
        lut = vtkLookupTable()
        lut.SetNumberOfTableValues(257)
        lut.SetRange(-1.0, 1.0)
        lut.Build()
        for i in range(257):
            t = (i / 256.0) * 2.0 - 1.0
            neutral = np.array([0.88, 0.90, 0.93])
            endpoint = np.array([0.90, 0.27, 0.13] if t >= 0 else [0.16, 0.55, 0.94])
            r, g, b = neutral * (1.0 - abs(t)) + endpoint * abs(t)
            lut.SetTableValue(i, r, g, b, 1.0)
        return lut

    @staticmethod
    def _lut_distance(max_value: float) -> vtkLookupTable:
        lut = vtkLookupTable()
        lut.SetNumberOfTableValues(256)
        lut.SetRange(0.0, max(1e-3, float(max_value)))
        lut.Build()
        for i in range(256):
            t = i / 255.0
            r = min(1.0, 0.2 + 0.8 * t)
            g = 0.9 - 0.7 * t
            b = 0.4 + 0.5 * (1.0 - t)
            lut.SetTableValue(i, r, max(0.0, g), min(1.0, b), 1.0)
        return lut

    @staticmethod
    def _make_glyph_actor(
        poly: vtkPolyData,
        sphere_radius: float,
        opacity: float,
        scalars: Optional[np.ndarray] = None,
        solid_rgb: Optional[Tuple[int, int, int]] = None,
        scalar_range: Optional[Tuple[float, float]] = None,
        lut: Optional[vtkLookupTable] = None,
        lut_builder: Optional[Callable[[], vtkLookupTable]] = None,
    ) -> vtkActor:
        sphere = vtkSphereSource()
        sphere.SetRadius(float(sphere_radius))
        sphere.SetThetaResolution(16)
        sphere.SetPhiResolution(16)

        if scalars is not None:
            arr = vtkFloatArray()
            arr.SetName("scores")
            for value in scalars:
                arr.InsertNextValue(float(value))
            poly.GetPointData().SetScalars(arr)
        elif solid_rgb is not None:
            OverlayBase._solid_color(poly, solid_rgb)

        glyph = vtkGlyph3D()
        glyph.SetInputData(poly)
        glyph.SetSourceConnection(sphere.GetOutputPort())
        glyph.ScalingOff()
        glyph.Update()

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(glyph.GetOutputPort())

        if scalars is not None:
            if scalar_range is None:
                scalar_range = (float(np.min(scalars)), float(np.max(scalars)))
            if abs(scalar_range[1] - scalar_range[0]) < 1e-6:
                scalar_range = (scalar_range[0], scalar_range[0] + 1e-3)
            mapper.SetScalarRange(*scalar_range)
            lut_obj = lut or (lut_builder() if lut_builder else OverlayBase._lut_blue_to_red())
            if scalar_range and lut_obj is not None:
                lut_obj.SetRange(*scalar_range)
            mapper.SetLookupTable(lut_obj)
            mapper.ScalarVisibilityOn()
        else:
            mapper.SetColorModeToDirectScalars()
            mapper.ScalarVisibilityOn()

        actor = vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetOpacity(float(opacity))
        return actor

    @staticmethod
    def _create_backbone_actor(structure: MolecularStructure, radius: float = 0.25) -> Optional[vtkActor]:
        coords = structure.as_numpy()
        if len(coords) < 2:
            return None
        points = vtkPoints()
        for coord in coords:
            points.InsertNextPoint(float(coord[0]), float(coord[1]), float(coord[2]))
        polyline = vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(len(coords))
        for idx in range(len(coords)):
            polyline.GetPointIds().SetId(idx, idx)
        cells = vtkCellArray()
        cells.InsertNextCell(polyline)
        poly = vtkPolyData()
        poly.SetPoints(points)
        poly.SetLines(cells)
        tube = vtkTubeFilter()
        tube.SetInputData(poly)
        tube.SetNumberOfSides(8)
        tube.SetRadius(radius)
        tube.CappingOn()
        tube.Update()
        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(tube.GetOutputPort())
        actor = vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0.9, 0.9, 0.9)
        actor.GetProperty().SetOpacity(0.25)
        return actor

    def _apply_overlay_annotations(
        self,
        renderer: vtkRenderer,
        *,
        lut: vtkLookupTable,
        scalar_title: str,
        status_text: str,
    ) -> None:
        bar = vtkScalarBarActor()
        bar.SetLookupTable(lut)
        bar.SetTitle(scalar_title)
        bar.SetNumberOfLabels(5)
        bar.SetWidth(0.10)
        bar.SetHeight(0.32)
        bar.SetPosition(0.88, 0.18)
        bar.SetMaximumWidthInPixels(110)
        bar.SetUnconstrainedFontSize(True)
        bar.SetLabelFormat("%.2f")
        title_prop = bar.GetTitleTextProperty()
        title_prop.SetFontSize(16)
        title_prop.SetColor(0.95, 0.95, 0.95)
        title_prop.ItalicOff()
        title_prop.ShadowOff()
        label_prop = bar.GetLabelTextProperty()
        label_prop.SetFontSize(14)
        label_prop.SetColor(0.9, 0.9, 0.9)
        label_prop.ItalicOff()
        label_prop.BoldOff()
        label_prop.ShadowOff()
        renderer.AddViewProp(bar)
        self._scalar_bar_actor = bar

        text_actor = vtkTextActor()
        text_actor.SetInput(status_text)
        text_actor.SetPosition(20, 20)
        text_prop = text_actor.GetTextProperty()
        text_prop.SetFontSize(16)
        text_prop.SetColor(0.95, 0.95, 0.95)
        text_prop.SetLineSpacing(1.2)
        renderer.AddViewProp(text_actor)
        self._status_actor = text_actor


class OverlayViewer(OverlayBase):
    """Superpose WT et mutant coloré par scores absolus [0..1]."""

    def __init__(
        self,
        wt_structure: str | Path,
        mutant_structure: str | Path,
        score_file: str | Path,
        *,
        show_backbone: bool = True,
        enable_tooltip: Optional[bool] = None,
        tooltip_verbose: Optional[bool] = None,
        enable_links: bool = False,
        link_threshold: Optional[float] = None,
        base_pair_config: Optional[BasePairOptions] = None,
        link_metric: Optional[str] = None,
        link_min_distance: Optional[float] = None,
        link_max_distance: Optional[float] = None,
    ) -> None:
        super().__init__(
            title="DeepStructGenomics - Overlay WT/Mutant",
            width=1200,
            height=850,
            show_backbone=show_backbone,
            enable_tooltip=enable_tooltip,
            tooltip_verbose=tooltip_verbose,
            enable_links=enable_links,
            link_threshold=link_threshold,
            base_pair_config=base_pair_config,
            link_metric=link_metric or "score",
            link_min_distance=link_min_distance,
            link_max_distance=link_max_distance,
        )
        self.wt_path = Path(wt_structure)
        self.mut_path = Path(mutant_structure)
        self.score_path = Path(score_file)

    def _build_scene(self, renderer: vtkRenderer) -> None:
        wt = load_structure(self.wt_path)
        mut = load_structure(self.mut_path)
        scores = load_score_table(self.score_path)

        wt_poly = self._points_poly(wt)
        mut_poly = self._points_poly(mut)

        wt_actor = self._make_glyph_actor(
            wt_poly,
            sphere_radius=1.05,
            opacity=0.22,
            scalars=None,
            solid_rgb=(170, 170, 210),
        )

        mut_scalars = map_scores_to_atoms(mut, scores)
        scalar_range: Tuple[float, float]
        if len(mut_scalars):
            smin = float(mut_scalars.min())
            smax = float(mut_scalars.max())
            if DEBUG:
                print(f"[overlay] mutant scalars range: {smin:.3f} -> {smax:.3f}")
            if abs(smax - smin) < 1e-6:
                smax = smin + 1e-3
            scalar_range = (smin, smax)
        else:
            scalar_range = (0.0, 1.0)
        lut = self._lut_blue_to_red()
        mut_actor = self._make_glyph_actor(
            mut_poly,
            sphere_radius=1.15,
            opacity=0.95,
            scalars=mut_scalars,
            solid_rgb=None,
            scalar_range=scalar_range,
            lut=lut,
        )

        renderer.AddActor(wt_actor)
        renderer.AddActor(mut_actor)
        self._wt_actor = wt_actor
        self._mut_actor = mut_actor
        self._backbone_actor = None
        if self.show_backbone:
            backbone = self._create_backbone_actor(mut)
            if backbone:
                renderer.AddActor(backbone)
                self._backbone_actor = backbone
        renderer.SetBackground(0.05, 0.05, 0.1)
        clamp_range = (0.0, 1.0)
        extra_lines: List[str] = []
        if self.enable_links:
            extra_lines.append("Links: WT -> MUT (same sequence position), color follows active LUT")
        base_pair_fn = getattr(self, "_base_pair_status_line", None)
        base_pair_line = base_pair_fn() if callable(base_pair_fn) else ""
        if base_pair_line:
            extra_lines.append(base_pair_line)
        status_text = build_status_text(
            mode_label="overlay",
            score_files=[self.score_path.name],
            scalar_range=scalar_range,
            clamp_range=clamp_range,
            legend_line="WT=gray transparent | MUT=colored score",
            extra_lines=tuple(extra_lines),
        )
        self._apply_overlay_annotations(renderer, lut=lut, scalar_title="score", status_text=status_text)

        self._interaction_actor = mut_actor
        self._tooltip_payload = {
            "structure": mut,
            "scalars": mut_scalars,
            "label": "score",
            "mut_scores": scores,
            "clamp": clamp_range,
            "scalar_range": scalar_range,
        }
        if self.link_metric == "distance" and self._link_distance_lookup:
            self._tooltip_payload["distance_lookup"] = self._link_distance_lookup

        if self.enable_links and len(mut_scalars):
            self._links_actor = self._try_create_links(
                renderer=renderer,
                wt_structure=wt,
                mut_structure=mut,
                scalar_values=self._scalars_by_position(mut, mut_scalars),
                scalar_range=scalar_range,
                lut=lut,
                mode="overlay",
            )
        else:
            self._links_actor = None

        self._setup_base_pairs(
            renderer,
            wt_structure=wt,
            mut_structure=mut,
            wt_scores=None,
            mut_scores=scores.position_scores,
        )


class OverlayDeltaViewer(OverlayBase):
    """Superpose WT (tube gris) et mutant coloré par delta (mutant - WT)."""

    def __init__(
        self,
        wt_structure: str | Path,
        mutant_structure: str | Path,
        wt_scores: str | Path,
        mutant_scores: str | Path,
        *,
        show_backbone: bool = True,
        enable_tooltip: Optional[bool] = None,
        tooltip_verbose: Optional[bool] = None,
        enable_links: bool = False,
        link_threshold: Optional[float] = None,
        base_pair_config: Optional[BasePairOptions] = None,
        link_metric: Optional[str] = None,
        link_min_distance: Optional[float] = None,
        link_max_distance: Optional[float] = None,
    ) -> None:
        super().__init__(
            title="DeepStructGenomics - Overlay Delta",
            width=1250,
            height=900,
            show_backbone=show_backbone,
            enable_tooltip=enable_tooltip,
            tooltip_verbose=tooltip_verbose,
            enable_links=enable_links,
            link_threshold=link_threshold,
            base_pair_config=base_pair_config,
            link_metric=link_metric or "delta",
            link_min_distance=link_min_distance,
            link_max_distance=link_max_distance,
        )
        self.wt_path = Path(wt_structure)
        self.mut_path = Path(mutant_structure)
        self.wt_scores_path = Path(wt_scores)
        self.mut_scores_path = Path(mutant_scores)

    def _build_scene(self, renderer: vtkRenderer) -> None:
        wt_structure = load_structure(self.wt_path)
        mut_structure = load_structure(self.mut_path)
        wt_scores = load_score_table(self.wt_scores_path)
        mut_scores = load_score_table(self.mut_scores_path)
        delta_scores = compute_delta_scores(wt_scores, mut_scores)
        stats = compute_score_statistics(delta_scores)
        mut_score_scalars = map_scores_to_atoms(mut_structure, mut_scores, clamp_min=0.0, clamp_max=1.0)
        mut_score_lookup = self._scalars_by_position(mut_structure, mut_score_scalars)
        if len(mut_score_scalars):
            score_min = float(mut_score_scalars.min())
            score_max = float(mut_score_scalars.max())
            if abs(score_max - score_min) < 1e-6:
                score_max = score_min + 1e-3
            score_range = (score_min, score_max)
        else:
            score_range = (0.0, 1.0)

        wt_poly = self._points_poly(wt_structure)
        mut_poly = self._points_poly(mut_structure)

        wt_actor = self._make_glyph_actor(
            wt_poly,
            sphere_radius=1.05,
            opacity=0.22,
            scalars=None,
            solid_rgb=(170, 170, 210),
        )

        delta_scalars = map_scores_to_atoms(mut_structure, delta_scores, clamp_min=-1.0, clamp_max=1.0)
        if len(delta_scalars):
            raw_min = float(delta_scalars.min())
            raw_max = float(delta_scalars.max())
            abs_max = max(abs(raw_min), abs(raw_max), 1e-3)
            scalar_range = (-abs_max, abs_max)
        else:
            scalar_range = (-1.0, 1.0)

        if DEBUG:
            pos_pct = stats["positive_ratio"] * 100.0
            neg_pct = stats["negative_ratio"] * 100.0
            print(
                "[overlay-delta] range {:.3f}->{:.3f} | +{:.1f}% / -{:.1f}%".format(
                    scalar_range[0], scalar_range[1], pos_pct, neg_pct
                )
            )

        lut = self._lut_diverging()
        mut_actor = self._make_glyph_actor(
            mut_poly,
            sphere_radius=1.2,
            opacity=0.98,
            scalars=delta_scalars,
            solid_rgb=None,
            scalar_range=scalar_range,
            lut=lut,
        )

        renderer.AddActor(wt_actor)
        renderer.AddActor(mut_actor)
        self._wt_actor = wt_actor
        self._mut_actor = mut_actor
        self._backbone_actor = None
        if self.show_backbone:
            backbone = self._create_backbone_actor(mut_structure)
            if backbone:
                renderer.AddActor(backbone)
                self._backbone_actor = backbone
        renderer.SetBackground(0.05, 0.05, 0.1)
        status_text = build_status_text(
            mode_label="overlay-delta",
            score_files=[self.wt_scores_path.name, self.mut_scores_path.name],
            scalar_range=scalar_range,
            clamp_range=(-1.0, 1.0),
            legend_line="WT=gray transparent | MUT=colored delta (mut-wt)",
            extra_lines=("Schema geometrique illustratif - pas une conformation moleculaire", *self._delta_status_lines()),
        )
        self._apply_overlay_annotations(renderer, lut=lut, scalar_title="delta", status_text=status_text)

        self._interaction_actor = mut_actor
        scalar_lookup = self._scalars_by_position(mut_structure, delta_scalars)
        self._tooltip_payload = {
            "structure": mut_structure,
            "scalars": delta_scalars,
            "label": "delta",
            "wt_scores": wt_scores,
            "mut_scores": mut_scores,
            "clamp": (-1.0, 1.0),
            "scalar_range": scalar_range,
        }
        if self.link_metric == "distance" and self._link_distance_lookup:
            self._tooltip_payload["distance_lookup"] = self._link_distance_lookup

        if self.enable_links and len(delta_scalars):
            link_values = scalar_lookup
            link_range = scalar_range
            link_mode = "overlay-delta"
            link_lut = lut
            if self.link_metric == "score":
                link_values = mut_score_lookup
                link_range = score_range
                link_mode = "overlay"
                link_lut = self._lut_blue_to_red()
            self._links_actor = self._try_create_links(
                renderer=renderer,
                wt_structure=wt_structure,
                mut_structure=mut_structure,
                scalar_values=link_values,
                scalar_range=link_range,
                lut=link_lut,
                mode=link_mode,
            )
        else:
            self._links_actor = None

        self._setup_base_pairs(
            renderer,
            wt_structure=wt_structure,
            mut_structure=mut_structure,
            wt_scores=wt_scores.position_scores,
            mut_scores=mut_scores.position_scores,
            lut=lut,
            scalar_range=scalar_range,
        )

    def _delta_status_lines(self) -> Tuple[str, ...]:
        lines: List[str] = ["Delta: 0=neutral | blue=decrease red=increase"]
        if self.enable_links:
            lines.append("Links: WT -> MUT (same sequence position), color follows active LUT")
        base_pair_fn = getattr(self, "_base_pair_status_line", None)
        base_pair_line = base_pair_fn() if callable(base_pair_fn) else ""
        if base_pair_line:
            lines.append(base_pair_line)
        return tuple(lines)
