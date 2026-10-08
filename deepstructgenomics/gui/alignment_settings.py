"""Collapsible per-analysis settings backed by the core configuration."""
from math import isqrt

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QCheckBox, QComboBox, QFormLayout, QLineEdit, QToolButton, QVBoxLayout, QWidget

from deepstructgenomics.alignment.config import AlignmentConfig, alignment_cell_count
from .widgets import button, label


class NumericField(QLineEdit):
    """Keep exact loaded values; Qt spin boxes otherwise round or clamp them."""
    def __init__(self, integer=False):
        super().__init__()
        self.integer = integer

    def setValue(self, value):
        self.setText(str(value))

    def value(self):
        try:
            return int(self.text()) if self.integer else float(self.text())
        except ValueError:
            raise ValueError("Paramètres avancés : renseigner une limite entière positive et des scores numériques finis.") from None


class AlignmentSettings(QWidget):
    def __init__(self, config):
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.toggle = QToolButton()
        self.toggle.setText("Paramètres avancés · alignement WT/MUT")
        self.toggle.setCheckable(True)
        self.toggle.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.toggle.setArrowType(Qt.ArrowType.RightArrow)
        layout.addWidget(self.toggle)
        self.panel = QWidget()
        form = QFormLayout(self.panel)
        self.enabled = QCheckBox("Activer l'alignement WT/MUT")
        form.addRow(self.enabled)
        self.profile = QComboBox()
        for name, cells in (("Standard — 4 M", AlignmentConfig().max_cells), ("Étendu — 25 M", 25_000_000),
                            ("Très étendu — 50 M", 50_000_000), ("Personnalisé", None)):
            self.profile.addItem(name, cells)
        form.addRow("Profil", self.profile)
        self.fields = {}
        for key, title in (("max_cells", "Limite maximale (cellules)"), ("max_length", "Longueur maximale (nt)"),
                           ("match_score", "Match"), ("mismatch_score", "Mismatch"),
                           ("gap_open_score", "Ouverture de gap"), ("gap_extend_score", "Extension de gap")):
            widget = NumericField(integer=key in ("max_cells", "max_length"))
            form.addRow(title, widget)
            self.fields[key] = widget
            widget.textChanged.connect(self.update_usage)
        self.overflow = QComboBox()
        self.overflow.addItem("Rejeter (défaut)", "reject")
        self.overflow.addItem("Comparer par position brute si dépassement", "positional")
        form.addRow("Dépassement", self.overflow)
        self.usage = label("")
        form.addRow(self.usage)
        form.addRow(label("La limite ne mesure pas la RAM. Le repli positionnel peut décaler les indels.\nCes réglages concernent uniquement l'alignement, pas le coût du repliement Nussinov."))
        form.addRow(button("Restaurer les valeurs par défaut", lambda: self.set_config(AlignmentConfig())))
        layout.addWidget(self.panel)
        self.panel.hide()
        self.toggle.toggled.connect(self.expand)
        self.profile.currentIndexChanged.connect(self.choose_profile)
        self.lengths = None
        self.set_config(config)

    def expand(self, checked=True):
        self.toggle.setChecked(checked)
        self.panel.setVisible(checked)
        self.toggle.setArrowType(Qt.ArrowType.DownArrow if checked else Qt.ArrowType.RightArrow)

    def set_config(self, config):
        self.enabled.setChecked(config.enabled)
        self.overflow.setCurrentIndex(self.overflow.findData(config.overflow_policy))
        self.profile.blockSignals(True)
        profile = self.profile.findData(config.max_cells)
        defaults = AlignmentConfig()
        if any(getattr(config, key) != getattr(defaults, key) for key in self.fields if key != "max_cells"):
            profile = -1
        self.profile.setCurrentIndex(profile if profile >= 0 else 3)
        for key, widget in self.fields.items():
            widget.setValue(getattr(config, key))
        self.profile.blockSignals(False)
        self.set_editability()
        self.update_usage()

    def choose_profile(self):
        cells = self.profile.currentData()
        if cells is not None:
            defaults = AlignmentConfig()
            for key, widget in self.fields.items():
                widget.setValue(cells if key == "max_cells" else getattr(defaults, key))
        self.set_editability()
        self.update_usage()

    def set_editability(self):
        for widget in self.fields.values():
            widget.setEnabled(self.profile.currentData() is None)

    def configuration(self):
        return AlignmentConfig(enabled=self.enabled.isChecked(), overflow_policy=self.overflow.currentData(),
                               **{key: widget.value() for key, widget in self.fields.items()})

    def set_lengths(self, wt=None, mut=None):
        self.lengths = (wt, mut) if wt is not None and mut is not None else None
        self.update_usage()

    def update_usage(self, *_):
        if "max_cells" not in self.fields or not hasattr(self, "usage"):
            return
        try:
            limit = self.fields["max_cells"].value()
            if limit <= 0:
                raise ValueError
        except ValueError:
            self.usage.setText("Limite invalide : saisir un entier strictement positif.")
            return
        text = f"Limite : {limit:,} cellules · ≈ {isqrt(limit):,} × {isqrt(limit):,} nt à tailles égales."
        if getattr(self, "lengths", None):
            wt, mut = self.lengths
            cells = alignment_cell_count(wt, mut)
            text += f"\nAnalyse actuelle : WT {wt:,} × MUT {mut:,} = {cells:,} cellules ({100 * cells / limit:.6g} % de la limite)."
        else:
            text += "\nMatrice actuelle : disponible après lecture des deux sources."
        self.usage.setText(text)
