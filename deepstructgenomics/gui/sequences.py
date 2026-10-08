"""Source selection form feeding PipelineInput without duplicating validation."""
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout,
    QLineEdit, QPlainTextEdit, QStackedWidget, QVBoxLayout, QWidget,
)
from .widgets import StateCard, button, label
from .alignment_settings import AlignmentSettings
from deepstructgenomics.alignment.config import load_alignment_config, load_sequence_inputs


class SequencesPage(QWidget):
    requested = Signal(dict)
    search_requested = Signal()
    failed = Signal(str)

    def __init__(self, alignment_config=None, sequence_inputs=None):
        super().__init__()
        self.ncbi_preview = None
        layout = QVBoxLayout(self)
        layout.setContentsMargins(26, 22, 26, 20)
        layout.setSpacing(14)
        layout.addWidget(label("Nouvelle analyse ARN", "heading"))
        layout.addWidget(label("Choisissez la référence, puis ajoutez un mutant pour comparer les appariements.", "subheading"))
        self.empty = StateCard("Choisissez votre référence WT", "Saisissez une séquence ci-dessous, importez un FASTA ou recherchez NCBI.",
                               actions=[("Rechercher NCBI", self.search_requested.emit)])
        layout.addWidget(self.empty)
        self.source = QComboBox()
        self.source.addItems(["Saisie directe", "Fichier FASTA", "Accession NCBI"])
        layout.addWidget(self.source)
        self.sources = QStackedWidget()
        self.sequence = QPlainTextEdit()
        self.sequence.setPlaceholderText("Séquence WT · A, C, G, U ou T")
        self.sequence.setMaximumHeight(125)
        self.sources.addWidget(self.sequence)
        fasta_widget = QWidget()
        fasta_layout = QHBoxLayout(fasta_widget)
        self.fasta = QLineEdit()
        self.fasta.setPlaceholderText("FASTA contenant une seule séquence")
        fasta_layout.addWidget(self.fasta)
        fasta_layout.addWidget(button("Parcourir…", lambda: self.choose_file(self.fasta)))
        self.sources.addWidget(fasta_widget)
        ncbi_widget = QWidget()
        form = QFormLayout(ncbi_widget)
        self.accession = QLineEdit()
        self.accession.setPlaceholderText("NR_000027.1 · accession nuccore exacte")
        self.email = QLineEdit()
        self.email.setPlaceholderText("Contact E-utilities (facultatif)")
        self.refresh = QCheckBox("Rafraîchir le cache NCBI")
        form.addRow("Accession", self.accession)
        form.addRow("E-mail", self.email)
        form.addRow(self.refresh)
        form.addRow(button("Rechercher un gène ou une séquence…", self.search_requested.emit))
        self.sources.addWidget(ncbi_widget)
        self.sources.setMaximumHeight(190)
        layout.addWidget(self.sources)
        self.source.currentIndexChanged.connect(self.sources.setCurrentIndex)
        form = QFormLayout()
        self.analysis_label = QLineEdit()
        self.analysis_label.setPlaceholderText("Facultatif · sinon identifiant de la source")
        form.addRow("Nom de l'analyse", self.analysis_label)
        layout.addLayout(form)
        self.reference_note = label("")
        self.reference_note.hide()
        self.source.currentIndexChanged.connect(
            lambda index: self.reference_note.setVisible(index == 2 and self.ncbi_preview is not None))
        layout.addWidget(self.reference_note)
        self.accession.textChanged.connect(self.invalidate_preview)
        mutant_box = QGroupBox("Mutant facultatif")
        mutant_layout = QVBoxLayout(mutant_box)
        self.mutant_source = QComboBox()
        self.mutant_source.addItems(["Sans mutant", "Séquence mutante", "FASTA mutant"])
        mutant_layout.addWidget(self.mutant_source)
        self.mutant = QPlainTextEdit()
        self.mutant.setPlaceholderText("Séquence mutante complète")
        self.mutant.setMaximumHeight(95)
        mutant_layout.addWidget(self.mutant)
        mutant_row = QHBoxLayout()
        self.mutant_fasta = QLineEdit()
        self.mutant_browse = button("Parcourir…", lambda: self.choose_file(self.mutant_fasta))
        mutant_row.addWidget(self.mutant_fasta)
        mutant_row.addWidget(self.mutant_browse)
        mutant_layout.addLayout(mutant_row)
        self.mutant_source.currentIndexChanged.connect(self.update_mutant)
        self.update_mutant(0)
        layout.addWidget(mutant_box)
        self.advanced = AlignmentSettings(alignment_config or load_alignment_config())
        self.align_mutant = self.advanced.enabled
        layout.addWidget(self.advanced)
        layout.addWidget(label("Méthode : Nussinov pondéré · scores heuristiques · aucune énergie MFE calculée.", "badge"))
        actions = QHBoxLayout()
        actions.addWidget(button("Remplir la démo C12A", self.fill_demo))
        actions.addStretch()
        actions.addWidget(button("Lancer l'analyse", self.submit, True))
        layout.addLayout(actions)
        layout.addStretch()
        self.source.currentIndexChanged.connect(self.update_empty)
        self.sequence.textChanged.connect(self.update_empty)
        self.fasta.textChanged.connect(self.update_empty)
        self.accession.textChanged.connect(self.update_empty)
        self.update_empty()
        inputs = load_sequence_inputs() if sequence_inputs is None else sequence_inputs
        if inputs.get("sequence"):
            self.sequence.setPlainText(inputs["sequence"])
        elif inputs.get("fasta_path"):
            self.fasta.setText(str(inputs["fasta_path"]))
            self.source.setCurrentIndex(1)
        if inputs.get("mutant_sequence"):
            self.mutant.setPlainText(inputs["mutant_sequence"])
            self.mutant_source.setCurrentIndex(1)
        elif inputs.get("mutant_fasta_path"):
            self.mutant_fasta.setText(str(inputs["mutant_fasta_path"]))
            self.mutant_source.setCurrentIndex(2)
        for signal in (self.sequence.textChanged, self.mutant.textChanged, self.source.currentIndexChanged,
                       self.mutant_source.currentIndexChanged):
            signal.connect(self.update_matrix)
        self.update_matrix()

    def update_matrix(self, *_):
        wt = len("".join(self.sequence.toPlainText().split())) if self.source.currentIndex() == 0 else None
        if self.source.currentIndex() == 2 and self.ncbi_preview:
            wt = len(self.ncbi_preview["sequence"])
        mut = len("".join(self.mutant.toPlainText().split())) if self.mutant_source.currentIndex() == 1 else None
        self.advanced.set_lengths(wt, mut)

    def update_empty(self, *_):
        value = (self.sequence.toPlainText(), self.fasta.text(), self.accession.text())[self.source.currentIndex()]
        self.empty.setVisible(not value.strip())

    def update_mutant(self, index):
        self.mutant.setVisible(index == 1)
        self.mutant_fasta.setVisible(index == 2)
        self.mutant_browse.setVisible(index == 2)

    def choose_file(self, field):
        path, _ = QFileDialog.getOpenFileName(self, "Séquence FASTA", "", "FASTA (*.fa *.fasta *.fna);;Tous (*)")
        if path:
            field.setText(path)

    def fill_demo(self):
        self.source.setCurrentIndex(0)
        self.sequence.setPlainText("GGGGAAAACCCC")
        self.mutant_source.setCurrentIndex(1)
        self.mutant.setPlainText("GGGGAAAACCCA")
        self.analysis_label.setText("demo_C12A")

    def invalidate_preview(self, *_):
        self.ncbi_preview = None
        self.reference_note.hide()

    def use_ncbi(self, preview):
        self.source.setCurrentIndex(2)
        self.accession.setText(preview["accession"])
        self.refresh.setChecked(False)
        self.ncbi_preview = preview
        length = len(preview["sequence"])
        self.reference_note.setText(f"Référence sélectionnée : {preview['accession']} · {length} nt\n{preview['description']}"
                                   + ("\nARN long : calcul Nussinov potentiellement coûteux, annulable." if length > 500 else ""))
        self.reference_note.show()
        self.update_matrix()

    def parameters(self):
        index = self.source.currentIndex()
        source = ({"sequence": self.sequence.toPlainText()}, {"fasta_path": self.fasta.text()},
                  {"accession": self.accession.text(), "ncbi_email": self.email.text() or None,
                   "refresh_cache": self.refresh.isChecked()})[index]
        if index != 2 and self.analysis_label.text().strip():
            source["sequence_label"] = self.analysis_label.text().strip()
        if index == 2 and self.ncbi_preview:
            source["expected_ncbi_sha256"] = self.ncbi_preview["sha256"]
        if self.mutant_source.currentIndex() == 1:
            source["mutant_sequence"] = self.mutant.toPlainText()
        elif self.mutant_source.currentIndex() == 2:
            source["mutant_fasta_path"] = self.mutant_fasta.text()
        source["alignment_config"] = self.advanced.configuration().to_dict()
        return source

    def submit(self):
        try:
            parameters = self.parameters()
        except ValueError as exc:
            self.failed.emit(str(exc))
            self.advanced.expand()
            return
        self.requested.emit(parameters)
