from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

import vtkmodules.vtkInteractionStyle  # noqa: F401
import vtkmodules.vtkRenderingOpenGL2  # noqa: F401

from vtkmodules.vtkIOImage import vtkPNGWriter
from vtkmodules.vtkRenderingCore import (
    vtkRenderWindow,
    vtkRenderer,
    vtkRenderWindowInteractor,
    vtkWindowToImageFilter,
)


@dataclass
class VTKStandaloneWindow:
    title: str = "DeepStructGenomics - VTK"
    width: int = 1000
    height: int = 700

    def run(
        self,
        build_scene: Callable[[vtkRenderer], None],
        *,
        post_setup: Optional[Callable[[vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor], None]] = None,
        export_path: Optional[str | Path] = None,
        export_only: bool = False,
        export_scale: int = 1,
    ) -> None:
        renderer = vtkRenderer()
        renderer.SetBackground(0.05, 0.05, 0.1)

        render_window = vtkRenderWindow()
        render_window.SetWindowName(self.title)
        render_window.SetSize(self.width, self.height)
        render_window.AddRenderer(renderer)

        interactor = vtkRenderWindowInteractor()
        interactor.SetRenderWindow(render_window)

        build_scene(renderer)
        if post_setup:
            post_setup(renderer, render_window, interactor)

        renderer.ResetCamera()
        render_window.Render()

        if export_path:
            self.export_png(render_window, Path(export_path), export_scale=export_scale)
            if export_only:
                return

        interactor.Initialize()
        interactor.Start()

    @staticmethod
    def export_png(render_window: vtkRenderWindow, output_path: Path, *, export_scale: int = 1) -> None:
        """Capture the current render window into a PNG image."""

        output_path.parent.mkdir(parents=True, exist_ok=True)
        w2i = vtkWindowToImageFilter()
        w2i.SetInput(render_window)
        scale = max(1, int(export_scale))
        w2i.SetScale(scale, scale, 1)
        w2i.Update()

        writer = vtkPNGWriter()
        writer.SetFileName(str(output_path))
        writer.SetInputConnection(w2i.GetOutputPort())
        writer.Write()
