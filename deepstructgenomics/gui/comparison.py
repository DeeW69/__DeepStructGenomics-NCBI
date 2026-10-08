"""Comparison page consuming exported structures, scores and reports."""
from pathlib import Path

from PySide6.QtCore import Signal
from deepstructgenomics.alignment.mapping import ComparisonMapping
from PySide6.QtWidgets import (
    QFileDialog, QHBoxLayout, QHeaderView, QSpinBox, QSplitter,
    QTableView, QTabWidget, QVBoxLayout, QWidget,
)
from .services import read_rna
from .widgets import FigureView, Metrics, StateCard, button, label
from .inspector import RnaInspector
from .alignment_table import AlignmentTableModel


class ComparisonPage(QWidget):
    failed = Signal(str)
    new_analysis = Signal()
    search_requested = Signal()

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
        self.alignment_note = label("", "subheading")
        layout.addWidget(self.alignment_note)
        self.empty = StateCard("Aucune séquence sélectionnée", "Recherchez NCBI, importez un FASTA ou saisissez une séquence.",
                               actions=[("Rechercher NCBI", self.search_requested.emit), ("Nouvelle analyse", self.new_analysis.emit)])
        layout.addWidget(self.empty)
        self.metrics = Metrics()
        layout.addWidget(self.metrics)
        self.tabs = QTabWidget()
        self.figure_view = FigureView()
        self.tabs.addTab(self.figure_view, "Structure secondaire")
        self.table = QTableView()
        self.table.verticalHeader().hide()
        self.table.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QTableView.SelectionMode.SingleSelection)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setResizeContentsPrecision(50)
        self.table.setEditTriggers(QTableView.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.clicked.connect(lambda index: self.position.setValue(index.row() + 1))
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
        self.technical = label("")
        self.tabs.addTab(self.technical, "Paramètres utilisés")
        self.splitter = QSplitter()
        self.splitter.addWidget(self.tabs)
        self.inspector = RnaInspector()
        self.position = QSpinBox()
        self.position.setMinimum(1)
        self.position.valueChanged.connect(self.show_position)
        self.position.setAccessibleName("Position ARN à inspecter, numérotation depuis 1")
        self.inspector.layout.insertWidget(1, self.position)
        self.splitter.addWidget(self.inspector)
        self.splitter.setStretchFactor(0, 1)
        self.splitter.setStretchFactor(1, 0)
        self.splitter.setSizes([820, 290])
        layout.addWidget(self.splitter, 1)
        self.splitter.hide()
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
        self.mapping = ComparisonMapping.from_data(data)
        run = data.get("alignment_run")
        if run:
            matrix, config = run["matrix"], run["configuration"]
            self.technical.setText(f"Configuration enregistrée pour cette analyse · {run['status']}\n\n"
                f"WT : {matrix['wt_length']:,} nt · MUT : {matrix['mut_length']:,} nt\n"
                f"Matrice théorique WT × MUT : {matrix['cells']:,} cellules\n"
                f"Limite : {config['max_cells']:,} cellules · utilisation : {100 * matrix['cells'] / config['max_cells']:.6g} %\n"
                f"Longueur maximale : {config['max_length']:,} nt\n"
                f"Match {config['match_score']:+g} · mismatch {config['mismatch_score']:+g} · ouverture {config['gap_open_score']:+g} · extension {config['gap_extend_score']:+g}\n"
                f"Politique de dépassement : {config['overflow_policy']}\n\nLa matrice théorique ne mesure pas la mémoire allouée.\n" + run.get("reason", ""))
        else:
            self.technical.setText("Ancien résultat : configuration effective non enregistrée. Aucun paramètre actuel n'est substitué aux paramètres historiques.")
        self.empty.hide()
        self.splitter.show()
        self.summary = (f"{len(data['reference']['sequence'])} bases · {metrics['wt_pairs']} paires WT"
                        + (f" / {metrics['mut_pairs']} MUT · {metrics['lost']} perdue(s)"
                           f" · |Δ| max {metrics['max_abs_delta']:.3f}" if data.get("mutant") else " · référence seule"))
        self.vtk_button.show()
        self.vtk_button.setEnabled(data.get("mutant") is not None)
        self.vtk_button.setText("Charger la vue 3D schématique" if data.get("mutant") else "Vue superposée : mutant requis")
        self.subtitle.setText(str(data.get("identifier", "ARN")) + (" · WT vs mutant" if data.get("mutant") else " · WT seul")
                              + " · Nussinov pondéré · positions depuis 1")
        alignment = self.mapping.alignment
        if alignment:
            c = alignment.counts
            self.alignment_note.setText(f"Alignement WT/MUT : actif · WT {len(data['reference']['sequence'])} nt / MUT {len(data['mutant']['sequence'])} nt · "
                f"Identité {alignment.identity:.1%} (gaps inclus) · {c['substitution']} substitution(s), {c['insertion']} insertion(s), {c['deletion']} délétion(s)"
                + ("\n" + " ".join(alignment.warnings) if alignment.warnings else ""))
        else:
            self.alignment_note.setText("Alignement WT/MUT : non calculé — comparaison brute par position. Les indels peuvent décaler les correspondances." if data.get("mutant") else "Référence seule : comparaison indisponible")
        self.metrics.set_values([
            ("Bases WT", len(data["reference"]["sequence"])),
            ("Paires WT", metrics["wt_pairs"]), ("Paires MUT", metrics["mut_pairs"]),
            ("Perdues / supprimées" if metrics["deleted_pairs"] else "Paires perdues", metrics["lost"]),
            ("Nouvelles (+ indels)" if metrics["inserted_pairs"] else "Paires gagnées", metrics["gained"]),
            ("Δ absolu maximal", f"{metrics['max_abs_delta']:.3f}" if metrics["max_abs_delta"] is not None else None),
        ])
        self.figure_view.set_figure(figure)
        self.delta_view.set_figure(delta_figure)
        self.figure_view.canvas.mpl_connect("button_press_event", self.pick_position)
        self.delta_view.canvas.mpl_connect("button_press_event", self.pick_position)
        self.contexts = list(figure.rna_contexts.values())
        self.subtitle.setText(self.subtitle.text() + " · " + figure.rna_layout_note)
        if data.get("mutant"):
            changes = self.mapping.classify_pairs(data["reference"]["base_pairs"], data["mutant"]["base_pairs"])
            def pair_text(pairs):
                ordered = sorted(pairs)
                return (", ".join(f"{i}–{j}" for i, j in ordered[:12]) or "aucune") + ("…" if len(ordered) > 12 else "")
            self.pair_changes.setText(("Colonnes d'alignement" if alignment else "Positions brutes") +
                f" · Perdues : {pair_text(changes['lost'])} · Nouvelles : {pair_text(changes['gained'])}"
                + (f"\nSupprimées avec une base : {pair_text(changes['deleted'])} · Nouvelles avec insertion : {pair_text(changes['inserted'])}" if alignment else ""))
        else:
            self.pair_changes.setText("Comparaison indisponible sans mutant.")
        n = len(self.mapping.columns)
        self.position.setMaximum(n)
        old_model = self.table.model()
        self.table.setModel(AlignmentTableModel(data, self.mapping, self.contexts, self.table))
        if old_model:
            old_model.deleteLater()
        self.position.setValue(1)
        self.show_position(1)
        self.export_button.setEnabled(True)
        self.tabs.setCurrentIndex(0)
        self.tabs.setTabText(1, "Alignement WT/MUT" if alignment else "Séquence")

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
        self.inspector.display(self.data, position, self.contexts)
        self.table.selectRow(position - 1)
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
