"""Qt host for the existing VTK overlay scene, with explicit user controls."""
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QApplication, QCheckBox, QHBoxLayout, QPlainTextEdit, QVBoxLayout, QWidget

from .widgets import button, label


class VTKView(QWidget):
    selected = Signal(int)

    def __init__(self, manifest, data):
        super().__init__()
        if QApplication.platformName() in ("offscreen", "minimal"):
            raise ValueError("Un affichage natif avec OpenGL est requis pour VTK ; utilisez la vue 2D hors écran.")
        # Lazy import: starting the app or viewing 2D does not create an OpenGL context.
        from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
        from vtkmodules.vtkRenderingCore import vtkCellPicker, vtkRenderer
        from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera
        from deepstructgenomics.visualization.io_structures import resolve_manifest_bundle
        from deepstructgenomics.visualization.tk_vtk_overlay import BasePairOptions, OverlayDeltaViewer
        bundle = resolve_manifest_bundle(manifest)
        layout = QVBoxLayout(self)
        layout.addWidget(label("Visualisation schématique · géométrie fondée sur la longueur, pas une conformation moléculaire.", "badge"))
        controls = QHBoxLayout()
        self.viewer = OverlayDeltaViewer(
            bundle["wt_structure"], bundle["mutant_structure"],
            bundle["wt_score_file"], bundle["mutant_score_file"],
            base_pair_config=BasePairOptions(enabled=True,
                wt_pairs=[(i + 1, j + 1) for i, j in data["reference"]["base_pairs"]],
                mut_pairs=[(i + 1, j + 1) for i, j in data["mutant"]["base_pairs"]]),
        )
        for title, key in (("Référence", "wt"), ("Mutant / Δ", "mutant"),
                           ("Squelette", "backbone"), ("Paires", "pairs")):
            check = QCheckBox(title)
            check.setChecked(True)
            check.toggled.connect(lambda value, layer=key: self.toggle(layer, value))
            controls.addWidget(check)
        controls.addStretch()
        controls.addWidget(button("Recentrer", self.reset))
        layout.addLayout(controls)
        self.widget = QVTKRenderWindowInteractor(self)
        layout.addWidget(self.widget, 1)
        self.renderer = vtkRenderer()
        self.widget.GetRenderWindow().AddRenderer(self.renderer)
        details = self.viewer.attach_scene(self.renderer)
        self.renderer.SetBackground(0.08, 0.14, 0.20)
        self.widget.SetInteractorStyle(vtkInteractorStyleTrackballCamera())
        self.picker = vtkCellPicker()
        self.picker.PickFromListOn()
        self.picker.AddPickList(self.viewer.pick_actor)
        self.widget.AddObserver("LeftButtonPressEvent", self.pick)
        self.widget.Initialize()
        self.reset()
        layout.addWidget(label("Bleu : diminution du score   ·   Blanc : aucun changement   ·   Orange : augmentation"))
        technical = QCheckBox("Informations techniques")
        layout.addWidget(technical)
        self.details = QPlainTextEdit(details)
        self.details.setReadOnly(True)
        self.details.setMaximumHeight(100)
        self.details.hide()
        technical.toggled.connect(self.details.setVisible)
        layout.addWidget(self.details)

    def toggle(self, layer, visible):
        self.viewer.set_layer_visible(layer, visible)
        self.widget.GetRenderWindow().Render()

    def reset(self):
        self.renderer.ResetCamera()
        self.widget.GetRenderWindow().Render()

    def pick(self, *_):
        x, y = self.widget.GetEventPosition()
        if self.picker.Pick(x, y, 0, self.renderer):
            position = self.viewer.position_near(self.picker.GetPickPosition())
            if position is not None:
                self.selected.emit(position)

    def shutdown(self):
        self.widget.Finalize()
