from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Optional, Sequence, Tuple

import vtkmodules.vtkInteractionStyle  # noqa: F401
import vtkmodules.vtkRenderingOpenGL2  # noqa: F401

from vtkmodules.vtkRenderingCore import vtkCamera
from vtkmodules.vtkIOImage import vtkPNGWriter
from vtkmodules.vtkRenderingCore import (
    vtkRenderWindow,
    vtkRenderer,
    vtkRenderWindowInteractor,
    vtkWindowToImageFilter,
)

from .export_helpers import compute_camera_angles

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
        export_view: str = "auto",
        export_hide_ui: bool = False,
        export_background: str = "dark",
        ui_elements_provider: Optional[Callable[[], Sequence[object]]] = None,
    ) -> None:
        renderer = vtkRenderer()
        renderer.SetBackground(0.05, 0.05, 0.1)

        render_window = vtkRenderWindow()
        render_window.SetWindowName(self.title)
        render_window.SetSize(self.width, self.height)
        render_window.AddRenderer(renderer)
        if export_only:
            render_window.SetOffScreenRendering(1)

        interactor = vtkRenderWindowInteractor()
        interactor.SetRenderWindow(render_window)

        build_scene(renderer)
        if post_setup:
            post_setup(renderer, render_window, interactor)

        renderer.ResetCamera()
        render_window.Render()

        if export_path:
            self._export_scene(
                renderer=renderer,
                render_window=render_window,
                export_path=Path(export_path),
                export_scale=export_scale,
                export_view=export_view,
                export_hide_ui=export_hide_ui,
                export_background=export_background,
                ui_elements_provider=ui_elements_provider,
                export_only=export_only,
            )
            if export_only:
                render_window.Finalize()
                return

        interactor.Initialize()
        interactor.Start()

    @staticmethod
    def export_png(
        render_window: vtkRenderWindow,
        output_path: Path,
        *,
        export_scale: int = 1,
        transparent: bool = False,
    ) -> None:
        """Capture the current render window into a PNG image."""

        output_path.parent.mkdir(parents=True, exist_ok=True)
        w2i = vtkWindowToImageFilter()
        w2i.SetInput(render_window)
        if transparent:
            w2i.SetInputBufferTypeToRGBA()
        else:
            w2i.SetInputBufferTypeToRGB()
        w2i.ReadFrontBufferOff()
        scale = max(1, int(export_scale))
        try:
            w2i.SetScale(scale, scale)
        except TypeError:
            # Windows VTK wheels sometimes expose SetScale(scale) only.
            w2i.SetScale(scale)
        w2i.Update()

        writer = vtkPNGWriter()
        writer.SetFileName(str(output_path))
        writer.SetInputConnection(w2i.GetOutputPort())
        writer.Write()

    def _export_scene(
        self,
        renderer: vtkRenderer,
        render_window: vtkRenderWindow,
        *,
        export_path: Path,
        export_scale: int,
        export_view: str,
        export_hide_ui: bool,
        export_background: str,
        ui_elements_provider: Optional[Callable[[], Sequence[object]]],
        export_only: bool,
    ) -> None:
        camera = renderer.GetActiveCamera()
        camera_state = self._snapshot_camera(camera)

        applied_view = export_view != "auto"
        if applied_view:
            azimuth, elevation = compute_camera_angles(export_view)
            camera.Azimuth(azimuth)
            camera.Elevation(elevation)
            renderer.ResetCameraClippingRange()

        hidden_elements: list[Tuple[object, int]] = []
        if export_hide_ui and ui_elements_provider:
            for actor in ui_elements_provider() or []:
                if actor is None:
                    continue
                if hasattr(actor, "GetVisibility") and hasattr(actor, "SetVisibility"):
                    original = actor.GetVisibility()
                    hidden_elements.append((actor, original))
                    actor.SetVisibility(0)

        original_bg = renderer.GetBackground()
        original_alpha = renderer.GetBackgroundAlpha()
        transparent = export_background == "transparent"
        if export_background == "white":
            renderer.SetBackground(1.0, 1.0, 1.0)
            renderer.SetBackgroundAlpha(1.0)
        elif transparent:
            render_window.SetAlphaBitPlanes(1)
            render_window.SetMultiSamples(0)
            renderer.SetBackground(0.0, 0.0, 0.0)
            renderer.SetBackgroundAlpha(0.0)
        elif export_background == "dark":
            renderer.SetBackground(0.05, 0.05, 0.1)
            renderer.SetBackgroundAlpha(1.0)

        render_window.Render()
        self.export_png(
            render_window,
            export_path,
            export_scale=export_scale,
            transparent=transparent,
        )

        if hidden_elements:
            for actor, visibility in hidden_elements:
                actor.SetVisibility(visibility)

        if export_background in {"white", "transparent"}:
            renderer.SetBackground(*original_bg)
            renderer.SetBackgroundAlpha(original_alpha)

        if applied_view and not export_only:
            self._restore_camera(camera, camera_state)
            renderer.ResetCameraClippingRange()

        render_window.Render()

    @staticmethod
    def _snapshot_camera(camera: vtkCamera) -> Dict[str, Tuple[float, float, float] | float]:
        return {
            "position": tuple(camera.GetPosition()),
            "focal_point": tuple(camera.GetFocalPoint()),
            "view_up": tuple(camera.GetViewUp()),
            "clipping_range": tuple(camera.GetClippingRange()),
            "view_angle": camera.GetViewAngle(),
        }

    @staticmethod
    def _restore_camera(camera: vtkCamera, state: Dict[str, Tuple[float, float, float] | float]) -> None:
        camera.SetPosition(*state["position"])
        camera.SetFocalPoint(*state["focal_point"])
        camera.SetViewUp(*state["view_up"])
        camera.SetClippingRange(*state["clipping_range"])
        camera.SetViewAngle(float(state["view_angle"]))
