"""Overlay viewers highlighting WT vs mutant differences."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Callable, Dict, Optional, Tuple

import numpy as np
import vtkmodules.vtkRenderingOpenGL2  # noqa: F401
from vtkmodules.vtkCommonCore import vtkFloatArray, vtkLookupTable, vtkPoints, vtkUnsignedCharArray
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData, vtkPolyLine
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
from .score_mapping import (
    ScoreTable,
    compute_delta_scores,
    compute_score_statistics,
    load_score_table,
    map_scores_to_atoms,
)
from .vtk_window import VTKStandaloneWindow

DEBUG = os.getenv("DSG_VIZ_DEBUG") == "1"


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
    ) -> None:
        self.window = VTKStandaloneWindow(title=title, width=width, height=height)
        self.show_backbone = show_backbone
        self.enable_tooltip = DEBUG if enable_tooltip is None else enable_tooltip
        self._tooltip_payload: Optional[Dict[str, object]] = None
        self._interaction_actor: Optional[vtkActor] = None

    def start(
        self,
        *,
        export_path: Optional[str | Path] = None,
        export_only: bool = False,
        export_scale: int = 1,
    ) -> None:
        self.window.run(
            self._build_scene,
            post_setup=self._post_setup if self.enable_tooltip else None,
            export_path=export_path,
            export_only=export_only,
            export_scale=export_scale,
        )

    def _build_scene(self, renderer: vtkRenderer) -> None:  # pragma: no cover - abstract
        raise NotImplementedError

    def _post_setup(self, renderer: vtkRenderer, render_window, interactor) -> None:
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
        renderer.AddActor2D(text)

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
        label = payload.get("label", "Score")
        value = scalars[point_id]
        lines = [
            f"pos={pos}",
            f"residue={residue}",
            f"{label}={value:.4f}",
        ]
        wt_val = OverlayBase._score_at_position(payload.get("wt_scores"), pos)
        mut_val = OverlayBase._score_at_position(payload.get("mut_scores"), pos)
        if label == "Δ":
            wt_text = f"{wt_val:.4f}" if wt_val is not None else "NA"
            mut_text = f"{mut_val:.4f}" if mut_val is not None else "NA"
            lines.append(f"WT={wt_text} MUT={mut_text}")
        elif mut_val is not None:
            lines.append(f"WT={'NA'} MUT={mut_val:.4f}")
        clamp = payload.get("clamp")
        if clamp:
            lines.append(f"range=[{clamp[0]:.3f},{clamp[1]:.3f}]")
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
        lut.SetNumberOfTableValues(256)
        lut.SetRange(-1.0, 1.0)
        lut.Build()
        for i in range(256):
            t = (i / 255.0) * 2.0 - 1.0
            if t >= 0:
                r = min(1.0, 0.35 + 0.65 * t)
                g = 0.2 * (1.0 - t)
                b = 0.2 * (1.0 - t)
            else:
                r = 0.2 * (1.0 + t)
                g = 0.35 * (1.0 + t)
                b = min(1.0, 0.5 - t * 0.5)
            lut.SetTableValue(i, r, g, b, 1.0)
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

    @staticmethod
    def _add_legend(
        renderer: vtkRenderer,
        lut: vtkLookupTable,
        title: str,
        subtitle: str = "",
    ) -> None:
        bar = vtkScalarBarActor()
        bar.SetLookupTable(lut)
        bar.SetTitle(title)
        bar.SetNumberOfLabels(5)
        bar.SetWidth(0.12)
        bar.SetHeight(0.45)
        bar.SetPosition(0.86, 0.10)
        renderer.AddActor2D(bar)

        text = vtkTextActor()
        lines = "WT: gris transparent\nMutant: couleur"
        if subtitle:
            lines += f"\n{subtitle}"
        text.SetInput(lines)
        text.SetPosition(20, 20)
        textprop = text.GetTextProperty()
        textprop.SetFontSize(18)
        textprop.SetColor(0.95, 0.95, 0.95)
        renderer.AddActor2D(text)


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
    ) -> None:
        super().__init__(
            title="DeepStructGenomics - Overlay WT/Mutant",
            width=1200,
            height=850,
            show_backbone=show_backbone,
            enable_tooltip=enable_tooltip,
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
        scalar_range = None
        lut = None
        subtitle = ""
        if len(mut_scalars):
            smin = float(mut_scalars.min())
            smax = float(mut_scalars.max())
            if DEBUG:
                print(f"[overlay] mutant scalars range: {smin:.3f} -> {smax:.3f}")
            if abs(smax - smin) < 1e-6:
                smax = smin + 1e-3
            scalar_range = (smin, smax)
            lut = self._lut_blue_to_red()
            subtitle = f"WT = gris transparent\nMutant = couleur = score mutant\nPlage affichee: {smin:.3f} -> {smax:.3f}"
        else:
            subtitle = "WT = gris transparent\nMutant = couleur = score mutant\nPlage affichee: 0.000 -> 1.000"
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
        if self.show_backbone:
            backbone = self._create_backbone_actor(mut)
            if backbone:
                renderer.AddActor(backbone)
        renderer.SetBackground(0.05, 0.05, 0.1)
        if lut is not None:
            info = subtitle + "\nLUT etiree sur la plage ci-dessus."
            self._add_legend(renderer, lut, "Score mutant (0..1)", info)

        self._interaction_actor = mut_actor
        clamp_range = scalar_range if scalar_range else (0.0, 1.0)
        self._tooltip_payload = {
            "structure": mut,
            "scalars": mut_scalars,
            "label": "score",
            "mut_scores": scores,
            "clamp": clamp_range,
        }


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
    ) -> None:
        super().__init__(
            title="DeepStructGenomics - Overlay Delta",
            width=1250,
            height=900,
            show_backbone=show_backbone,
            enable_tooltip=enable_tooltip,
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
        if self.show_backbone:
            backbone = self._create_backbone_actor(mut_structure)
            if backbone:
                renderer.AddActor(backbone)
        renderer.SetBackground(0.05, 0.05, 0.1)
        subtitle = (
            "WT = gris transparent\n"
            "Mutant = couleur = delta (mutant - WT)\n"
            "Bleu: delta < 0 | Rouge: delta > 0 | Gris ~0\n"
            f"Plage symetrique: {scalar_range[0]:.3f} -> {scalar_range[1]:.3f}"
        )
        self._add_legend(renderer, lut, "Delta score (mutant - WT)", subtitle)

        self._interaction_actor = mut_actor
        self._tooltip_payload = {
            "structure": mut_structure,
            "scalars": delta_scalars,
            "label": "delta",
            "wt_scores": wt_scores,
            "mut_scores": mut_scores,
            "clamp": scalar_range,
        }
