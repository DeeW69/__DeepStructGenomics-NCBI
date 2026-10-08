"""Search, inspect and explicitly select a versioned NCBI reference."""
from PySide6.QtCore import Signal, Qt
from PySide6.QtWidgets import (QCheckBox, QFormLayout, QHBoxLayout, QHeaderView,
    QLineEdit, QPlainTextEdit, QSplitter, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)
from .widgets import StateCard, button, label


class NcbiPage(QWidget):
    requested = Signal(str, dict)
    chosen = Signal(dict)

    def __init__(self):
        super().__init__()
        self.search_parameters = None
        self.result = None
        self.preview = None
        layout = QVBoxLayout(self)
        layout.setContentsMargins(26, 22, 26, 20)
        layout.setSpacing(12)
        layout.addWidget(label("Rechercher une séquence NCBI", "heading"))
        layout.addWidget(label("Gène, accession ou mots-clés → aperçu → référence WT", "subheading"))
        self.query = QLineEdit()
        self.query.setPlaceholderText("BRCA1, TP53, HOTAIR ou une accession…")
        self.organism = QLineEdit()
        self.organism.setPlaceholderText("Organisme facultatif, ex. Homo sapiens")
        form = QFormLayout()
        form.addRow("Recherche", self.query)
        form.addRow("Organisme", self.organism)
        layout.addLayout(form)
        filters = QHBoxLayout()
        self.rna_only = QCheckBox("ARN uniquement")
        self.rna_only.setChecked(True)
        filters.addWidget(self.rna_only)
        filters.addStretch()
        filters.addWidget(button("Rechercher", self.search, True))
        layout.addLayout(filters)
        self.query.returnPressed.connect(self.search)
        self.state_card = StateCard("Trouvez votre référence WT", "Saisissez un gène, une accession ou des mots-clés, puis cliquez sur Rechercher.")
        self.note = self.state_card.body
        layout.addWidget(self.state_card)
        split = QSplitter(Qt.Orientation.Vertical)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["Accession", "Organisme", "Type", "Longueur", "Description"])
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeMode.Stretch)
        self.table.itemSelectionChanged.connect(self.selection_changed)
        self.table.cellDoubleClicked.connect(lambda *_: self.fetch_preview())
        split.addWidget(self.table)
        panel = QWidget()
        preview_layout = QVBoxLayout(panel)
        self.details = label("Sélectionnez une notice, puis chargez son aperçu.", "badge")
        preview_layout.addWidget(self.details)
        self.sequence = QPlainTextEdit()
        self.sequence.setReadOnly(True)
        self.sequence.setPlaceholderText("La séquence choisie apparaîtra ici, sans lancer de prédiction.")
        preview_layout.addWidget(self.sequence)
        split.addWidget(panel)
        split.setSizes([280, 200])
        layout.addWidget(split, 1)
        actions = QHBoxLayout()
        self.previous = button("Page précédente", lambda: self.search(-1))
        self.next = button("Page suivante", lambda: self.search(1))
        self.preview_button = button("Charger l'aperçu", self.fetch_preview)
        self.use = button("Utiliser cette séquence comme WT", self.use_sequence, True)
        for control in (self.previous, self.next, self.preview_button, self.use):
            control.setEnabled(False)
            actions.addWidget(control)
        layout.addLayout(actions)
        self.query.textEdited.connect(self.invalidate)
        self.organism.textEdited.connect(self.invalidate)
        self.rna_only.toggled.connect(self.invalidate)

    def invalidate(self, *_):
        self.result = None
        self.search_parameters = None
        self.table.setRowCount(0)
        self.selection_changed()
        self.previous.setEnabled(False)
        self.next.setEnabled(False)
        self.state_card.set_state("empty", "Recherche à lancer", "Critères modifiés : relancer la recherche.")

    def search(self, direction=0):
        if direction and self.result and self.search_parameters:
            parameters = {**self.search_parameters, "start": max(0, self.result["start"] + direction * self.result["page_size"])}
        else:
            if not self.query.text().strip():
                self.state_card.set_state("empty", "Aucun terme saisi", "Saisir un gène, une accession ou des mots-clés.")
                return
            parameters = {"term": self.query.text().strip(), "organism": self.organism.text().strip(),
                          "rna_only": self.rna_only.isChecked(), "start": 0, "page_size": 20}
        self.search_parameters = parameters
        self.table.setRowCount(0)
        self.selection_changed()
        self.previous.setEnabled(False)
        self.next.setEnabled(False)
        self.state_card.set_state("empty", "Recherche demandée", "Les résultats remplaceront la liste précédente.")
        self.requested.emit("ncbi_search", parameters)

    def show_results(self, result):
        self.result = result
        self.table.setRowCount(len(result["rows"]))
        for index, row in enumerate(result["rows"]):
            for column, key in enumerate(("accession", "organism", "molecule", "length", "description")):
                item = QTableWidgetItem(str(row[key]))
                item.setToolTip(str(row[key]))
                self.table.setItem(index, column, item)
        count, start = len(result["rows"]), result["start"]
        self.state_card.set_state("success" if count else "empty", "Sélectionnez une notice" if count else "Aucun résultat",
                                  f"{start + 1}–{start + count} sur {result['total']} notices · {result['term']}" if count
                                  else "Essayez un autre terme ou élargissez les filtres.")
        self.note.setToolTip(result.get("translated_query", result["query"]))
        if result.get("unmatched_terms"):
            self.note.setText(self.note.text() + " · Termes non reconnus : " + ", ".join(result["unmatched_terms"]))
        self.previous.setEnabled(start > 0)
        self.next.setEnabled(start + result["page_size"] < result["total"])

    def selected(self):
        index = self.table.currentRow()
        return self.result["rows"][index] if self.result and 0 <= index < len(self.result["rows"]) else None

    def selection_changed(self):
        self.preview = None
        self.use.setEnabled(False)
        self.sequence.clear()
        row = self.selected()
        self.preview_button.setEnabled(row is not None)
        self.details.setText(f"{row['accession']} · {row['organism']} · {row['length']} nt\n{row['description']}" if row
                             else "Sélectionnez une notice, puis chargez son aperçu.")
        if row:
            self.state_card.set_state("empty", "Notice sélectionnée", "Chargez l'aperçu pour valider cette référence avant de l'utiliser.")

    def fetch_preview(self):
        row = self.selected()
        if row:
            self.preview = None
            self.use.setEnabled(False)
            self.state_card.set_state("empty", "Aperçu demandé", "Chargement de la séquence et de sa provenance…")
            self.requested.emit("ncbi_preview", {"accession": row["accession"]})

    def show_preview(self, result):
        row = self.selected()
        if not row or row["accession"] != result["accession"]:
            return
        self.preview = result
        sequence = result["sequence"]
        # Bound widget size for large records, while preserving the full cached sequence.
        excerpt = sequence[:5000]
        self.sequence.setPlainText("\n".join(excerpt[i:i + 80] for i in range(0, len(excerpt), 80)))
        self.details.setText(f"{result['accession']} · {row['organism']} · {len(sequence)} nt\n{result['description']}")
        self.state_card.set_state("success", "Aperçu prêt · sélectionnez cette référence WT", "T converti en U · la sélection ne lance pas l'analyse."
                          + (" Affichage limité aux 5 000 premières bases." if len(sequence) > 5000 else "")
                          + (" ARN long : le calcul Nussinov peut être très coûteux et reste annulable." if len(sequence) > 500 else ""))
        self.use.setEnabled(True)

    def use_sequence(self):
        if self.preview:
            self.chosen.emit(self.preview)
