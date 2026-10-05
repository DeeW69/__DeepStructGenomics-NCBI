from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Tuple

import vtkmodules.vtkRenderingOpenGL2  # noqa: F401
from vtkmodules.vtkCommonCore import vtkPoints, vtkUnsignedCharArray
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData, vtkPolyLine
from vtkmodules.vtkFiltersCore import vtkGlyph3D, vtkTubeFilter
from vtkmodules.vtkFiltersSources import vtkSphereSource
from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkRenderer

from .io_structures import MolecularStructure, load_structure
from .vtk_window import VTKStandaloneWindow

NUCLEOTIDE_COLORS: Dict[str, Tuple[int, int, int]] = {
    "A": (230, 140, 70),
    "C": (60, 160, 230),
    "G": (80, 200, 90),
    "U": (230, 90, 140),
    "T": (230, 90, 140),
}


class MoleculeViewer:
    """Affiche une structure (PDB/mmCIF) en glyphes sphériques (standalone VTK)."""

    def __init__(self, structure_path: str | Path, sphere_radius: float = 1.2, show_backbone: bool = False) -> None:
        self.structure_path = Path(structure_path)
        self.sphere_radius = float(sphere_radius)
        self.show_backbone = show_backbone
        self.window = VTKStandaloneWindow(title="DeepStructGenomics - Molecule", width=1100, height=800)

    def _build_scene(self, renderer: vtkRenderer) -> None:
        structure = load_structure(self.structure_path)
        coords = structure.as_numpy()
        if coords.size == 0:
            raise ValueError("Structure vide: aucun atome trouvé.")

        points = vtkPoints()
        colors = vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)
        colors.SetName("atom_colors")

        for atom in structure.atoms:
            x, y, z = atom.coord
            points.InsertNextPoint(float(x), float(y), float(z))

            resname = atom.residue_name.strip()
            # pour les pseudo-PDB on a "NTP", sinon on essaie de dériver une lettre
            nucleotide = resname[-1].upper() if resname else "N"
            rgb = NUCLEOTIDE_COLORS.get(nucleotide, (200, 200, 200))
            colors.InsertNextTypedTuple(rgb)

        poly = vtkPolyData()
        poly.SetPoints(points)
        poly.GetPointData().SetScalars(colors)

        sphere = vtkSphereSource()
        sphere.SetRadius(self.sphere_radius)
        sphere.SetThetaResolution(16)
        sphere.SetPhiResolution(16)

        glyph = vtkGlyph3D()
        glyph.SetInputData(poly)
        glyph.SetSourceConnection(sphere.GetOutputPort())
        glyph.ScalingOff()
        glyph.Update()

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(glyph.GetOutputPort())
        mapper.SetColorModeToDirectScalars()
        mapper.ScalarVisibilityOn()

        actor = vtkActor()
        actor.SetMapper(mapper)

        renderer.AddActor(actor)
        if self.show_backbone:
            backbone = self._create_backbone_actor(structure)
            if backbone:
                renderer.AddActor(backbone)
        renderer.SetBackground(0.05, 0.05, 0.1)

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
            export_path=export_path,
            export_only=export_only,
            export_scale=export_scale,
            export_view=export_view,
            export_hide_ui=export_hide_ui,
            export_background=export_background,
        )

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
