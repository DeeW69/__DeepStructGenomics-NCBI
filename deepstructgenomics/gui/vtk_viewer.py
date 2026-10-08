"""Qt host for the existing VTK overlay scene, with explicit user controls."""
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QApplication, QCheckBox, QHBoxLayout, QPlainTextEdit, QVBoxLayout, QWidget

from .widgets import button, label
from deepstructgenomics.alignment.mapping import ComparisonMapping


class VTKView(QWidget):
    selected = Signal(int)

    def __init__(self, manifest, data):
        super().__init__()
        if QApplication.platformName() in ("offscreen", "minimal"):
            raise ValueError("Un affichage natif avec OpenGL est requis pour VTK ; utilisez la vue 2D hors écran.")
        # Lazy import: starting the app or viewing 2D does not create an OpenGL context.
        from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
        from vtkmodules.vtkRenderingCore import vtkBillboardTextActor3D, vtkCellPicker, vtkRenderer
        from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera
        from deepstructgenomics.visualization.io_structures import load_structure, resolve_manifest_bundle
        from deepstructgenomics.visualization.tk_vtk_overlay import BasePairOptions, OverlayDeltaViewer
        bundle = resolve_manifest_bundle(manifest)
        self.mapping = ComparisonMapping.from_data(data)
        layout = QVBoxLayout(self)
        layout.addWidget(label("Visualisation schématique · géométrie fondée sur la longueur, pas une conformation moléculaire.", "badge"))
        controls = QHBoxLayout()
        self.viewer = OverlayDeltaViewer(
            bundle["wt_structure"], bundle["mutant_structure"],
            bundle["wt_score_file"], bundle["mutant_score_file"],
            glyph_scale=.58,
            alignment=self.mapping.alignment,
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
        self.viewer.pick_actor.GetProperty().SetSpecular(0.)
        self.viewer.pick_actor.GetProperty().SetInterpolationToPhong()
        self.viewer.pick_actor.GetProperty().SetAmbient(.35)
        self.viewer.pick_actor.GetProperty().SetDiffuse(.65)
        self.base_labels = []
        structure = load_structure(bundle["mutant_structure"])
        if len(structure.atoms) <= 80:
            for index, coord in enumerate(structure.as_numpy()):
                text = vtkBillboardTextActor3D()
                base = data["mutant"]["sequence"][index]
                column = self.mapping.columns[self.mapping.to_column["mutant"][index + 1] - 1]
                wt_base = column.reference_base
                text.SetInput(f"{index + 1} {wt_base}>{base}" if wt_base and wt_base != base else f"{index + 1} {base}")
                text.SetPosition(float(coord[0]), float(coord[1]) + 1.1, float(coord[2]) + .6)
                text.GetTextProperty().SetFontSize(13)
                text.GetTextProperty().SetColor(.86, .92, .97)
                self.renderer.AddActor(text)
                self.base_labels.append(text)
        letters = QCheckBox("Bases")
        letters.setChecked(bool(self.base_labels))
        letters.setEnabled(bool(self.base_labels))
        letters.setToolTip("Lettres affichées pour les séquences de 80 bases maximum.")
        letters.toggled.connect(self.toggle_labels)
        controls.insertWidget(4, letters)
        self.renderer.SetBackground(0.08, 0.14, 0.20)
        self.widget.SetInteractorStyle(vtkInteractorStyleTrackballCamera())
        self.picker = vtkCellPicker()
        self.picker.PickFromListOn()
        self.picker.AddPickList(self.viewer.pick_actor)
        self.widget.AddObserver("LeftButtonPressEvent", self.pick)
        self.widget.Initialize()
        self.reset()
        layout.addWidget(label("Bleu : diminution · Blanc : aucun changement · Rouge : augmentation · Gris : pas d'homologue"))
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

    def toggle_labels(self, visible):
        for actor in self.base_labels:
            actor.SetVisibility(visible)
        self.widget.GetRenderWindow().Render()

    def pick(self, *_):
        x, y = self.widget.GetEventPosition()
        if self.picker.Pick(x, y, 0, self.renderer):
            position = self.viewer.position_near(self.picker.GetPickPosition())
            if position is not None:
                self.selected.emit(self.mapping.to_column["mutant"][position])

    def shutdown(self):
        self.widget.Finalize()
