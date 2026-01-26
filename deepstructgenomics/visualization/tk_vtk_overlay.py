"""Viewer overlaying WT vs mutant structures with scalar coloring."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import numpy as np

try:  # pragma: no cover - optional dependency
    from vtkmodules.vtkCommonCore import vtkFloatArray, vtkPoints, vtkUnsignedCharArray
    from vtkmodules.vtkCommonDataModel import vtkPolyData
    from vtkmodules.vtkFiltersCore import vtkGlyph3D
    from vtkmodules.vtkFiltersSources import vtkSphereSource
    from vtkmodules.vtkRenderingCore import vtkActor, vtkLookupTable, vtkPolyDataMapper
except ImportError as exc:  # pragma: no cover
    raise ImportError("VTK est requis pour la visualisation 3D (pip install vtk).") from exc

from .io_structures import MolecularStructure, load_structure
from .score_mapping import ScoreTable, load_score_table, map_scores_to_atoms
from .tk_vtk_minimal import BaseTkVtkApp


class OverlayViewer(BaseTkVtkApp):
    """Superpose WT et mutant avec coloration du mutant via un score [0,1]."""

    def __init__(self, wt_structure: str | Path, mutant_structure: str | Path, score_file: str | Path) -> None:
        super().__init__(title="DeepStructGenomics - Overlay WT/Mutant")
        self.wt = load_structure(wt_structure)
        self.mutant = load_structure(mutant_structure)
        self.scores = load_score_table(score_file)
        self._add_structures()

    def _build_points(self, structure: MolecularStructure) -> vtkPolyData:
        coords = structure.as_numpy()
        points = vtkPoints()
        for x, y, z in coords:
            points.InsertNextPoint(float(x), float(y), float(z))
        poly = vtkPolyData()
        poly.SetPoints(points)
        return poly

    def _add_structures(self) -> None:
        wt_poly = self._build_points(self.wt)
        mutant_poly = self._build_points(self.mutant)

        sphere = vtkSphereSource()
        sphere.SetRadius(1.1)
        sphere.SetThetaResolution(12)
        sphere.SetPhiResolution(12)

        wt_actor = self._build_actor(
            poly=wt_poly,
            sphere=sphere,
            color=(0.6, 0.6, 0.9),
            opacity=0.35,
            scalars=None,
        )
        mutant_scalars = map_scores_to_atoms(self.mutant, self.scores)
        mutant_actor = self._build_actor(
            poly=mutant_poly,
            sphere=sphere,
            color=(0.9, 0.4, 0.4),
            opacity=0.85,
            scalars=mutant_scalars,
        )

        self.renderer.AddActor(wt_actor)
        self.renderer.AddActor(mutant_actor)
        self.renderer.SetBackground(0.07, 0.07, 0.09)
        self.renderer.ResetCamera()

    def _build_actor(
        self,
        poly: vtkPolyData,
        sphere: vtkSphereSource,
        color: Tuple[float, float, float],
        opacity: float,
        scalars: Optional[np.ndarray],
    ) -> vtkActor:
        glyph = vtkGlyph3D()
        glyph.SetInputData(poly)
        glyph.SetSourceConnection(sphere.GetOutputPort())
        glyph.ScalingOff()
        if scalars is not None and len(scalars):
            scalar_array = vtkFloatArray()
            scalar_array.SetName("scores")
            for value in scalars:
                scalar_array.InsertNextValue(float(value))
            poly.GetPointData().SetScalars(scalar_array)
        glyph.Update()

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(glyph.GetOutputPort())

        if scalars is not None and len(scalars):
            mapper.SetScalarRange(0.0, 1.0)
            mapper.SetLookupTable(self._build_lookup_table())
        else:
            mapper.SetColorModeToDirectScalars()
            solid_colors = vtkUnsignedCharArray()
            solid_colors.SetNumberOfComponents(3)
            solid_colors.SetName("solid")
            rgb = tuple(int(channel * 255) for channel in color)
            for _ in range(poly.GetNumberOfPoints()):
                solid_colors.InsertNextTypedTuple(rgb)
            poly.GetPointData().SetScalars(solid_colors)

        actor = vtkActor()
        actor.SetMapper(mapper)
        if scalars is None or not len(scalars):
            actor.GetProperty().SetColor(*color)
        actor.GetProperty().SetOpacity(opacity)
        return actor

    def _build_lookup_table(self) -> vtkLookupTable:
        lut = vtkLookupTable()
        lut.SetNumberOfTableValues(256)
        lut.Build()
        for i in range(256):
            t = i / 255.0
            lut.SetTableValue(i, 0.1 + 0.9 * t, 0.2 * (1 - t), 0.9 * (1 - t), 1.0)
        return lut
