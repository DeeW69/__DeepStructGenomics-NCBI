"""Tkinter + VTK scaffold for desktop visualization."""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk
from typing import Tuple

try:  # pragma: no cover - optional heavy dependency
    import vtkmodules.vtkInteractionStyle  # noqa: F401
    from vtkmodules.tk.vtkTkRenderWindowInteractor import vtkTkRenderWindowInteractor
    from vtkmodules.vtkFiltersSources import vtkSphereSource
    from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkRenderWindow, vtkRenderer
except ImportError as exc:  # pragma: no cover - import-time guard
    raise ImportError("VTK est requis pour la visualisation 3D (pip install vtk).") from exc


class BaseTkVtkApp:
    """Base class that embeds a VTK renderer in a Tkinter window."""

    def __init__(self, title: str = "DeepStructGenomics Viewer", size: Tuple[int, int] = (900, 700)) -> None:
        self.title = title
        self.size = size
        self.root = tk.Tk()
        self.root.title(self.title)
        self.root.geometry(f"{size[0]}x{size[1]}")
        self.frame = ttk.Frame(self.root)
        self.frame.pack(fill=tk.BOTH, expand=True)
        self.renderer = vtkRenderer()
        self.render_window = vtkRenderWindow()
        self.render_window.AddRenderer(self.renderer)
        self.vtk_widget = vtkTkRenderWindowInteractor(self.frame, rw=self.render_window, width=size[0], height=size[1])
        self.vtk_widget.pack(fill=tk.BOTH, expand=True)

    def start(self) -> None:
        """Start Tk main loop."""

        self.vtk_widget.Initialize()
        self.renderer.ResetCamera()
        self.render_window.Render()
        self.root.mainloop()


class MinimalVTKViewer(BaseTkVtkApp):
    """Minimal viewer showcasing camera / interaction controls."""

    def __init__(self) -> None:
        super().__init__(title="DeepStructGenomics - Viewer minimal")
        self._build_scene()

    def _build_scene(self) -> None:
        """Create a demo actor."""

        sphere_source = vtkSphereSource()
        sphere_source.SetRadius(10.0)
        sphere_source.SetThetaResolution(32)
        sphere_source.SetPhiResolution(32)

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(sphere_source.GetOutputPort())

        actor = vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0.2, 0.8, 0.8)

        self.renderer.AddActor(actor)
        self.renderer.SetBackground(0.05, 0.05, 0.1)
        self.renderer.ResetCamera()
