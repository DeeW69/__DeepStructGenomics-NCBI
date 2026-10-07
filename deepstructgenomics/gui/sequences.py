"""Source selection form feeding PipelineInput without duplicating validation."""
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout,
    QLineEdit, QPlainTextEdit, QStackedWidget, QVBoxLayout, QWidget,
)
from .widgets import button, label


class SequencesPage(QWidget):
    requested = Signal(dict)

    def __init__(self):
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(26, 22, 26, 20)
        layout.setSpacing(14)
        layout.addWidget(label("Nouvelle analyse ARN", "heading"))
        layout.addWidget(label("Choisissez la référence, puis ajoutez un mutant pour comparer les appariements.", "subheading"))
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
        self.sources.addWidget(ncbi_widget)
        self.sources.setMaximumHeight(150)
        layout.addWidget(self.sources)
        self.source.currentIndexChanged.connect(self.sources.setCurrentIndex)
        form = QFormLayout()
        self.analysis_label = QLineEdit()
        self.analysis_label.setPlaceholderText("Facultatif · sinon identifiant de la source")
        form.addRow("Nom de l'analyse", self.analysis_label)
        layout.addLayout(form)
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
        layout.addWidget(label("Méthode : Nussinov pondéré · scores heuristiques · aucune énergie MFE calculée.", "badge"))
        actions = QHBoxLayout()
        actions.addWidget(button("Remplir la démo C12A", self.fill_demo))
        actions.addStretch()
        actions.addWidget(button("Lancer l'analyse", self.submit, True))
        layout.addLayout(actions)
        layout.addStretch()

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

    def parameters(self):
        index = self.source.currentIndex()
        source = ({"sequence": self.sequence.toPlainText()}, {"fasta_path": self.fasta.text()},
                  {"accession": self.accession.text(), "ncbi_email": self.email.text() or None,
                   "refresh_cache": self.refresh.isChecked()})[index]
        if index != 2 and self.analysis_label.text().strip():
            source["sequence_label"] = self.analysis_label.text().strip()
        if self.mutant_source.currentIndex() == 1:
            source["mutant_sequence"] = self.mutant.toPlainText()
        elif self.mutant_source.currentIndex() == 2:
            source["mutant_fasta_path"] = self.mutant_fasta.text()
        return source

    def submit(self):
        self.requested.emit(self.parameters())
