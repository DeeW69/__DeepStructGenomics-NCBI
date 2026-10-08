"""Comparison page consuming exported structures, scores and reports."""
from pathlib import Path

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QFileDialog, QHBoxLayout, QHeaderView, QSpinBox, QSplitter,
    QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget,
)
from .services import position_details, read_rna
from .widgets import FigureView, Metrics, button, label


class ComparisonPage(QWidget):
    failed = Signal(str)

    def __init__(self):
        super().__init__()
        self.manifest = None
        self.data = None
        self.vtk_view = None
        layout = QVBoxLayout(self)
        layout.setContentsMargins(26, 22, 26, 20)
        layout.addWidget(label("Comparaison structurelle", "heading"))
        self.subtitle = label("Ouvrez un résultat ou lancez une analyse ARN.", "subheading")
        layout.addWidget(self.subtitle)
        self.metrics = Metrics()
        layout.addWidget(self.metrics)
        self.tabs = QTabWidget()
        self.figure_view = FigureView()
        self.tabs.addTab(self.figure_view, "Structure secondaire")
        self.table = QTableWidget(0, 8)
        self.table.setHorizontalHeaderLabels(["Position", "WT", "Mutant", "Score WT", "Score MUT", "Δ", "Contexte WT", "Contexte MUT"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.cellClicked.connect(lambda row, _: self.position.setValue(row + 1))
        self.tabs.addTab(self.table, "Séquence")
        delta_page = QWidget()
        delta_layout = QVBoxLayout(delta_page)
        delta_layout.setContentsMargins(0, 0, 0, 0)
        self.delta_view = FigureView()
        delta_layout.addWidget(self.delta_view, 1)
        self.pair_changes = label("")
        delta_layout.addWidget(self.pair_changes)
        delta_layout.addWidget(label("Énergie MFE indisponible : le moteur Nussinov pondéré ne la calcule pas."))
        self.tabs.addTab(delta_page, "Delta structural")
        self.vtk_host = QWidget()
        self.vtk_layout = QVBoxLayout(self.vtk_host)
        self.vtk_button = button("Charger la vue 3D schématique", self.load_vtk)
        self.vtk_layout.addWidget(self.vtk_button)
        self.tabs.addTab(self.vtk_host, "3D schématique")
        self.splitter = QSplitter()
        self.splitter.addWidget(self.tabs)
        panel = QWidget()
        panel.setMinimumWidth(200)
        panel.setMaximumWidth(270)
        detail_layout = QVBoxLayout(panel)
        detail_layout.addWidget(label("Inspecter une position", "subheading"))
        self.position = QSpinBox()
        self.position.setMinimum(1)
        self.position.valueChanged.connect(self.show_position)
        detail_layout.addWidget(self.position)
        self.detail = label("Aucune séquence chargée.")
        detail_layout.addWidget(self.detail)
        detail_layout.addStretch()
        detail_layout.addWidget(label("Scores heuristiques, sans unité. Les paires représentent une prédiction secondaire.", "badge"))
        self.splitter.addWidget(panel)
        self.splitter.setStretchFactor(0, 1)
        self.splitter.setStretchFactor(1, 0)
        self.splitter.setSizes([850, 230])
        layout.addWidget(self.splitter, 1)
        actions = QHBoxLayout()
        self.export_button = button("Exporter la structure 2D PNG…", self.export_png)
        self.tabs.currentChanged.connect(lambda index: self.export_button.setText(
            "Exporter le delta PNG…" if index == 2 else "Exporter la structure 2D PNG…"))
        self.export_button.setEnabled(False)
        actions.addWidget(self.export_button)
        actions.addWidget(label("Rapports JSON, Markdown et CSV disponibles via « Dossier des résultats »."))
        layout.addLayout(actions)

    def load(self, manifest):
        from deepstructgenomics.visualization.secondary_view import build_secondary_figure
        from deepstructgenomics.visualization.secondary_diagram import build_structure_diagram
        data, metrics = read_rna(manifest)
        figure = build_structure_diagram(data, compact=True)
        delta_figure = build_secondary_figure(data, compact=True)
        self.shutdown()
        self.manifest, self.data = Path(manifest), data
        self.summary = (f"{len(data['reference']['sequence'])} bases · {metrics['wt_pairs']} paires WT"
                        + (f" / {metrics['mut_pairs']} MUT · {metrics['lost']} perdue(s)"
                           f" · |Δ| max {metrics['max_abs_delta']:.3f}" if data.get("mutant") else " · référence seule"))
        self.vtk_button.show()
        self.vtk_button.setEnabled(data.get("mutant") is not None)
        self.vtk_button.setText("Charger la vue 3D schématique" if data.get("mutant") else "Vue superposée : mutant requis")
        self.subtitle.setText(str(data.get("identifier", "ARN")) + " · Nussinov pondéré · positions depuis 1")
        if data.get("mutant") and len(data["mutant"]["sequence"]) != len(data["reference"]["sequence"]):
            self.subtitle.setText(self.subtitle.text() + " · longueurs différentes, sans alignement")
        self.metrics.set_values([
            ("Paires WT", metrics["wt_pairs"]), ("Paires MUT", metrics["mut_pairs"]),
            ("Paires perdues", metrics["lost"]), ("Paires gagnées", metrics["gained"]),
            ("Δ absolu maximal", f"{metrics['max_abs_delta']:.3f}" if metrics["max_abs_delta"] is not None else None),
        ])
        self.figure_view.set_figure(figure)
        self.delta_view.set_figure(delta_figure)
        self.figure_view.canvas.mpl_connect("button_press_event", self.pick_position)
        self.delta_view.canvas.mpl_connect("button_press_event", self.pick_position)
        self.contexts = list(figure.rna_contexts.values())
        self.subtitle.setText(self.subtitle.text() + " · " + figure.rna_layout_note)
        if data.get("mutant"):
            wt_pairs = set(map(tuple, data["reference"]["base_pairs"]))
            mut_pairs = set(map(tuple, data["mutant"]["base_pairs"]))
            def pair_text(pairs):
                ordered = sorted(pairs)
                return (", ".join(f"{i + 1}–{j + 1}" for i, j in ordered[:12]) or "aucune") + ("…" if len(ordered) > 12 else "")
            self.pair_changes.setText(f"Paires perdues : {pair_text(wt_pairs - mut_pairs)}\nNouvelles paires : {pair_text(mut_pairs - wt_pairs)}")
        else:
            self.pair_changes.setText("Comparaison indisponible sans mutant.")
        n = max(len(data["reference"]["sequence"]), len(data["mutant"]["sequence"]) if data.get("mutant") else 0)
        self.position.setMaximum(n)
        self.table.setRowCount(n)
        for index in range(n):
            detail = position_details(data, index + 1)
            wt, mut = detail["reference"], detail["mutant"]
            values = [index + 1, wt["base"] if wt else "—", mut["base"] if mut else "—",
                      f"{wt['score']:.3f}" if wt else "—", f"{mut['score']:.3f}" if mut else "—",
                      f"{detail['delta']:+.3f}" if detail["delta"] is not None else "—",
                      self.contexts[0][index] if index < len(self.contexts[0]) else "—",
                      self.contexts[1][index] if len(self.contexts) > 1 and index < len(self.contexts[1]) else "—"]
            for column, value in enumerate(values):
                self.table.setItem(index, column, QTableWidgetItem(str(value)))
        self.position.setValue(1)
        self.show_position(1)
        self.export_button.setEnabled(True)
        self.tabs.setCurrentIndex(0)

    def pick_position(self, event):
        if getattr(event.canvas, "toolbar", None) and event.canvas.toolbar.mode:
            return
        if event.canvas.figure is self.figure_view.figure:
            from deepstructgenomics.visualization.secondary_diagram import picked_position
            position = picked_position(self.figure_view.figure, event)
            if position is not None:
                self.position.setValue(position)
        elif event.inaxes and event.xdata is not None:
            self.position.setValue(max(1, min(self.position.maximum(), round(event.xdata))))

    def show_position(self, position):
        if not self.data:
            return
        detail = position_details(self.data, position)
        lines = [f"Position {position}"]
        for column, (key, title) in enumerate((("reference", "Référence"), ("mutant", "Mutant"))):
            entry = detail[key]
            if entry:
                pair = f"{position} ↔ {entry['partner']}" if entry["partner"] else "aucune"
                lines += ["", f"{title} : {entry['base']}", f"Score : {entry['score']:.3f}", f"Paire : {pair}"]
                lines.append(self.contexts[column][position - 1])
            else:
                lines += ["", f"{title} : indisponible"]
        if detail["delta"] is not None:
            lines += ["", f"Variation : {detail['delta']:+.3f}"]
            wt, mut = detail["reference"], detail["mutant"]
            if wt["partner"] and not mut["partner"]:
                lines += ["", "La paire prédite à cette position disparaît dans le mutant."]
            elif mut["partner"] and not wt["partner"]:
                lines += ["", "Une paire prédite apparaît dans le mutant."]
        self.detail.setText("\n".join(lines))
        from deepstructgenomics.visualization.secondary_diagram import select_position
        select_position(self.figure_view.figure, position)

    def load_vtk(self):
        if self.vtk_view or not self.data or not self.data.get("mutant"):
            return
        try:
            from .vtk_viewer import VTKView
            self.vtk_view = VTKView(self.manifest, self.data)
            self.vtk_view.selected.connect(self.position.setValue)
            self.vtk_layout.addWidget(self.vtk_view)
            self.vtk_button.hide()
        except Exception as exc:
            self.failed.emit(f"Vue VTK indisponible : {exc}")

    def export_png(self):
        path, _ = QFileDialog.getSaveFileName(self, "Exporter la comparaison", "comparaison.png", "PNG (*.png)")
        if path:
            from deepstructgenomics.visualization.secondary_view import build_secondary_figure
            from deepstructgenomics.visualization.secondary_diagram import build_structure_diagram
            from matplotlib import pyplot as plt
            figure = build_secondary_figure(self.data) if self.tabs.currentIndex() == 2 else build_structure_diagram(self.data)
            try:
                figure.savefig(path, dpi=180)
            except OSError as exc:
                self.failed.emit(str(exc))
            finally:
                plt.close(figure)

    def shutdown(self):
        if self.vtk_view:
            self.vtk_view.shutdown()
            self.vtk_layout.removeWidget(self.vtk_view)
            self.vtk_view.deleteLater()
            self.vtk_view = None
