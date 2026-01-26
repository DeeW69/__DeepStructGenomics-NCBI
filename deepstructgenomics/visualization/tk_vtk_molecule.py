"""Tkinter window embedding a VTK molecule viewer."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

try:  # pragma: no cover - optional heavy dependency
    from vtkmodules.vtkCommonCore import vtkPoints, vtkUnsignedCharArray
    from vtkmodules.vtkCommonDataModel import vtkPolyData
    from vtkmodules.vtkFiltersCore import vtkGlyph3D
    from vtkmodules.vtkFiltersSources import vtkSphereSource
    from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper
except ImportError as exc:  # pragma: no cover
    raise ImportError("VTK est requis pour la visualisation 3D (pip install vtk).") from exc

from .io_structures import load_structure
from .tk_vtk_minimal import BaseTkVtkApp


NUCLEOTIDE_COLORS: Dict[str, Tuple[float, float, float]] = {
    "A": (0.9, 0.5, 0.2),
    "C": (0.2, 0.6, 0.9),
    "G": (0.4, 0.8, 0.4),
    "U": (0.9, 0.3, 0.5),
    "T": (0.9, 0.3, 0.5),
}


class MoleculeViewer(BaseTkVtkApp):
    """Display atoms from a PDB/mmCIF file using glyph spheres."""

    def __init__(self, structure_path: str | Path, sphere_radius: float = 1.2) -> None:
        super().__init__(title="DeepStructGenomics - Molecule 3D")
        self.structure = load_structure(structure_path)
        self.sphere_radius = sphere_radius
        self._add_atoms()

    def _add_atoms(self) -> None:
        molecule = self.structure
        coords = molecule.as_numpy()
        if coords.size == 0:
            raise ValueError("Structure vide: aucun atome trouve.")

        points = vtkPoints()
        for x, y, z in coords:
            points.InsertNextPoint(float(x), float(y), float(z))

        colors = vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)
        colors.SetName("atom_colors")
        for atom in molecule.atoms:
            resname = atom.residue_name.strip()
            nucleotide = resname[-1] if resname else "N"
            color = NUCLEOTIDE_COLORS.get(nucleotide.upper(), (0.8, 0.8, 0.8))
            rgb = tuple(int(channel * 255) for channel in color)
            colors.InsertNextTypedTuple(rgb)

        poly_data = vtkPolyData()
        poly_data.SetPoints(points)
        poly_data.GetPointData().SetScalars(colors)

        sphere = vtkSphereSource()
        sphere.SetRadius(self.sphere_radius)
        sphere.SetThetaResolution(12)
        sphere.SetPhiResolution(12)

        glyph = vtkGlyph3D()
        glyph.SetInputData(poly_data)
        glyph.SetSourceConnection(sphere.GetOutputPort())
        glyph.ScalingOff()
        glyph.Update()

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(glyph.GetOutputPort())
        mapper.SetColorModeToDirectScalars()

        actor = vtkActor()
        actor.SetMapper(mapper)

        self.renderer.AddActor(actor)
        self.renderer.SetBackground(0.08, 0.08, 0.12)
        self.renderer.ResetCamera()
