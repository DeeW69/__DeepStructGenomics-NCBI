from __future__ import annotations

from pathlib import Path
from typing import Optional

from vtkmodules.vtkFiltersSources import vtkCubeSource
from vtkmodules.vtkRenderingAnnotation import vtkAxesActor
from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkRenderer

from .vtk_window import VTKStandaloneWindow


class MinimalVTKViewer:
    def __init__(self) -> None:
        self.window = VTKStandaloneWindow(title="DeepStructGenomics - Minimal", width=900, height=700)

    @staticmethod
    def _build_scene(renderer: vtkRenderer) -> None:
        # 1) Cube (géométrie sûre)
        cube = vtkCubeSource()
        cube.SetXLength(50.0)
        cube.SetYLength(50.0)
        cube.SetZLength(50.0)
        cube.Update()

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(cube.GetOutputPort())

        actor = vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0.2, 0.8, 0.8)  # bien visible
        actor.GetProperty().SetSpecular(0.3)
        actor.GetProperty().SetSpecularPower(20)

        renderer.AddActor(actor)

        # 2) Axes (debug visuel : si tu vois ça, la scène est OK)
        axes = vtkAxesActor()
        axes.SetTotalLength(80, 80, 80)
        renderer.AddActor(axes)

        renderer.SetBackground(0.05, 0.05, 0.1)

    def start(self, *, export_path: Optional[str | Path] = None, export_only: bool = False) -> None:
        self.window.run(self._build_scene, export_path=export_path, export_only=export_only)
