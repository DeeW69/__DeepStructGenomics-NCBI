"""Hi-C region form and views of the existing advanced analysis report."""
import json
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QFormLayout, QHBoxLayout,
    QLineEdit, QScrollArea, QSpinBox, QSplitter, QVBoxLayout, QWidget,
)
from .widgets import FigureView, Metrics, StateCard, button, label


def hic_figure(report, mode="Contacts équilibrés"):
    figure = Figure(figsize=(9, 6), constrained_layout=True, facecolor="white")
    region = report["region"]
    if mode == "Reconstruction 3D":
        ax = figure.add_subplot(111, projection="3d")
        reconstruction = report.get("reconstruction") or {}
        coords = reconstruction.get("coordinates", [])
        if coords:
            xyz = np.array([[row[key] for key in ("x", "y", "z")] for row in coords])
            ax.scatter(*xyz.T, c=[row["bin"] for row in coords], cmap="viridis", s=32)
            for i in range(len(coords) - 1):
                if coords[i + 1]["bin"] == coords[i]["bin"] + 1:
                    ax.plot(*xyz[i:i + 2].T, color="#96a8b8", linewidth=1)
        else:
            ax.text2D(.1, .5, "Reconstruction indisponible", transform=ax.transAxes)
        ax.set(title=f"MDS inférée · {reconstruction.get('status', 'indisponible')}",
               xlabel="x (u.a.)", ylabel="y (u.a.)", zlabel="z (u.a.)")
    else:
        ax = figure.add_subplot(111)
        matrix = np.array(report["balanced_matrix"], dtype=float)
        masked = report["normalization"]["masked_bins"]
        matrix[masked, :] = np.nan
        matrix[:, masked] = np.nan
        start, end, step = region["start"], region["end"], region["bin_size"]
        image = ax.imshow(matrix, origin="lower", cmap="magma", extent=(start, end, start, end))
        figure.colorbar(image, ax=ax, label="Poids équilibré")
        domains = report.get("domains") or {}
        for boundary in domains.get("selected_boundaries", []):
            ax.axvline(start + boundary * step, color="cyan", linewidth=.8)
            ax.axhline(start + boundary * step, color="cyan", linewidth=.8)
        if mode == "Boucles candidates":
            loops = [row for row in (report.get("loops") or {}).get("tests", []) if row["significant"]]
            ax.scatter([start + (row["bin1"] + .5) * step for row in loops],
                       [start + (row["bin2"] + .5) * step for row in loops],
                       facecolors="none", edgecolors="cyan", s=80, linewidths=1.5)
        ax.set(title=f"{region['chromosome']} · {mode}", xlabel="Position (bp, base 0)", ylabel="Position (bp, base 0)")
    return figure


class HicPage(QWidget):
    requested = Signal(dict)
    failed = Signal(str)

    def __init__(self):
        super().__init__()
        self.report = None
        self.result_path = None
        layout = QVBoxLayout(self)
        layout.setContentsMargins(26, 22, 26, 20)
        layout.addWidget(label("Hi-C Explorer", "heading"))
        layout.addWidget(label("Contacts cis, boucles et domaines candidats · reconstruction inférée", "subheading"))
        split = QSplitter()
        form_widget = QWidget()
        self.form_widget = form_widget
        form_widget.setMaximumWidth(335)
        form_layout = QVBoxLayout(form_widget)
        self.path = QLineEdit()
        self.path.setPlaceholderText("Fichier TSV de contacts bruts")
        form_layout.addWidget(self.path)
        form_layout.addWidget(button("Choisir les contacts…", self.choose_file))
        form = QFormLayout()
        self.assembly = QLineEdit()
        self.chromosome = QLineEdit()
        self.start = QSpinBox()
        self.end = QSpinBox()
        self.bin_size = QSpinBox()
        for field in (self.start, self.end, self.bin_size):
            field.setRange(0, 2_000_000_000)
            field.setSingleStep(10000)
        self.bin_size.setMinimum(1)
        self.bin_size.setValue(10000)
        self.end.setValue(320000)
        for title, field in (("Assemblage", self.assembly), ("Chromosome", self.chromosome),
                             ("Début (bp)", self.start), ("Fin exclue (bp)", self.end), ("Bin (bp)", self.bin_size)):
            form.addRow(title, field)
        form_layout.addLayout(form)
        self.missing = QCheckBox("Les contacts absents sont\ndes zéros observés")
        form_layout.addWidget(self.missing)
        form_layout.addWidget(label("Région alignée de 400 bins maximum. Fournir les comptes bruts complets, pas une sélection de contacts."))
        form_layout.addWidget(button("Analyser la région", self.submit, True))
        self.demo = button("Remplir l'exemple synthétique", self.fill_demo)
        form_layout.addWidget(self.demo)
        form_layout.addWidget(label("SciPy requis (extra research). FDR 0,05 · 999 permutations · graine 46 · fenêtre 3 bins."))
        form_layout.addStretch()
        form_scroll = QScrollArea()
        form_scroll.setWidgetResizable(True)
        form_scroll.setMaximumWidth(350)
        form_scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        form_scroll.setWidget(form_widget)
        split.addWidget(form_scroll)
        results = QWidget()
        results_layout = QVBoxLayout(results)
        self.metrics = Metrics()
        results_layout.addWidget(self.metrics)
        self.mode = QComboBox()
        self.mode.addItems(["Contacts équilibrés", "Boucles candidates", "Reconstruction 3D"])
        self.mode.currentTextChanged.connect(self.redraw)
        results_layout.addWidget(self.mode)
        self.figure_view = FigureView()
        self.figure_view.empty.set_state("empty", "Aucun contact analysé", "Importez des contacts TSV et définissez une région, ou remplissez l'exemple synthétique.",
                                         [("Remplir l'exemple", self.fill_demo)])
        results_layout.addWidget(self.figure_view, 1)
        self.result_state = StateCard()
        self.result_state.hide()
        self.status = self.result_state.body
        results_layout.addWidget(self.result_state)
        self.mode.setEnabled(False)
        self.export = button("Exporter la figure PNG…", self.export_png)
        self.export.setEnabled(False)
        results_layout.addWidget(self.export)
        split.addWidget(results)
        split.setStretchFactor(1, 1)
        layout.addWidget(split, 1)
        layout.addWidget(label("Méthodes exploratoires : q-values dépendantes des modèles nuls, domaines candidats et coordonnées 3D en unités arbitraires."))

    def choose_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Contacts Hi-C", "", "TSV (*.tsv);;Tous (*)")
        if path:
            self.path.setText(path)

    def fill_demo(self):
        path = Path(__file__).resolve().parents[2] / "data" / "examples" / "hic_region.tsv"
        if not path.is_file():
            self.failed.emit("Exemple disponible dans le dépôt source : data/examples/hic_region.tsv.")
            return
        self.path.setText(str(path))
        self.assembly.setText("synthetic")
        self.chromosome.setText("chrSynthetic")
        self.start.setValue(0)
        self.end.setValue(320000)
        self.bin_size.setValue(10000)
        self.missing.setChecked(True)

    def parameters(self):
        return {"path": self.path.text(), "assembly": self.assembly.text().strip(),
                "chromosome": self.chromosome.text().strip(), "start": self.start.value(),
                "end": self.end.value(), "bin_size": self.bin_size.value(),
                "missing_as_zero": self.missing.isChecked()}

    def submit(self):
        self.requested.emit(self.parameters())

    def load(self, path):
        report = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(report, dict) or report.get("schema") != 1 or not isinstance(report.get("region"), dict):
            raise ValueError("Ce fichier n'est pas un rapport Hi-C avancé compatible.")
        matrix = np.asarray(report.get("balanced_matrix"), dtype=float)
        if (matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]
                or not 2 <= len(matrix) <= 400 or not np.isfinite(matrix).all()
                or report.get("status") not in ("ok", "not_converged")
                or report["region"].get("bins") != len(matrix)):
            raise ValueError("Matrice ou statut Hi-C invalide.")
        # Build before replacing the currently displayed report.
        figure = hic_figure(report, self.mode.currentText())
        self.report, self.result_path = report, Path(path)
        loops, domains = report.get("loops"), report.get("domains")
        self.metrics.set_values([
            ("Bins", report["region"]["bins"]),
            ("Boucles candidates", sum(row["significant"] for row in loops["tests"]) if loops else None),
            ("Frontières", len(domains["selected_boundaries"]) if domains else None),
        ])
        self.result_state.set_state("success" if report["status"] == "ok" else "warning",
                                    "Contacts analysés" if report["status"] == "ok" else "Équilibrage non convergé",
                                    "Résultats exploratoires disponibles dans les vues ci-dessus." if report["status"] == "ok"
                                    else "Aucun test ni reconstruction exécuté. Vérifiez la région et les contacts fournis.")
        self.mode.setEnabled(True)
        self.figure_view.set_figure(figure)
        self.export.setEnabled(True)

    def redraw(self, *_):
        if self.report:
            self.figure_view.set_figure(hic_figure(self.report, self.mode.currentText()))

    def export_png(self):
        path, _ = QFileDialog.getSaveFileName(self, "Exporter Hi-C", "hic.png", "PNG (*.png)")
        if path:
            try:
                self.figure_view.figure.savefig(path, dpi=180)
            except OSError as exc:
                self.failed.emit(str(exc))
