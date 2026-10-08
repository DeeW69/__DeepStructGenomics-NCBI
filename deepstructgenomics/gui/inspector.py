"""Readable per-base inspection, independent of the selected visual tab."""
from html import escape

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QFormLayout, QFrame, QScrollArea, QVBoxLayout, QWidget

from .services import inspection_details
from .widgets import StateCard, label


class RnaInspector(QScrollArea):
    def __init__(self):
        super().__init__()
        self.setWidgetResizable(True)
        self.setFrameShape(QFrame.Shape.NoFrame)
        self.setMinimumWidth(245)
        self.setMaximumWidth(330)
        content = QWidget()
        self.layout = QVBoxLayout(content)
        self.layout.setContentsMargins(8, 0, 4, 0)
        self.heading = label("Inspecter une position", "subheading")
        self.layout.addWidget(self.heading)
        self.fields = {}
        for key, title in (("reference", "WT"), ("mutant", "MUT")):
            card = QFrame()
            card.setObjectName("card")
            rows = QFormLayout(card)
            rows.addRow(label(title, "stateTitle"))
            self.fields[key] = {}
            for field, caption in (("base", "Base"), ("score", "Score"), ("context", "Contexte"), ("partner", "Partenaire")):
                value = label("—")
                rows.addRow(caption, value)
                self.fields[key][field] = value
            self.layout.addWidget(card)
        self.delta = label("Δ indisponible", "inspectorDelta")
        self.layout.addWidget(self.delta)
        self.badges = [StateCard() for _ in range(3)]
        for badge in self.badges:
            badge.hide()
            self.layout.addWidget(badge)
        self.local_title = label("Séquence locale")
        self.layout.addWidget(self.local_title)
        self.local = label("")
        self.local.setObjectName("localSequence")
        self.local.setTextFormat(Qt.TextFormat.RichText)
        self.local.setWordWrap(False)
        self.layout.addWidget(self.local)
        self.note = label("Scores heuristiques, sans unité. Paires prédites ; positions depuis 1.", "subheading")
        self.layout.addWidget(self.note)
        self.layout.addStretch()
        self.setWidget(content)

    def display(self, data, position, contexts):
        detail = inspection_details(data, position)
        self.heading.setText(f"Position {position}")
        for column, key in enumerate(("reference", "mutant")):
            entry = detail[key]
            values = {"base": entry["base"] if entry else "Indisponible",
                      "score": f"{entry['score']:.3f}" if entry else "—",
                      "context": contexts[column][position - 1] if entry else "—",
                      "partner": entry["partner_label"] if entry else "—"}
            for field, value in values.items():
                self.fields[key][field].setText(value)
        self.delta.setText(f"Δ  {detail['delta']:+.3f}" if detail["delta"] is not None else "Δ indisponible")
        for index, badge in enumerate(self.badges):
            if index < len(detail["changes"]):
                state, title, text = detail["changes"][index]
                badge.set_state(state, title, "" if state == "substitution" else text)
            else:
                badge.hide()
        start, end = detail["local_start"], detail["local_end"]
        self.local_title.setText(f"Séquence locale · {start}–{end}")
        rows = []
        for key, title in (("reference", "WT "), ("mutant", "MUT")):
            bases = [f'<span style="background-color:#127c82;color:white"><b>{escape(base)}</b></span>'
                     if start + i == position else escape(base) for i, base in enumerate(detail["local"][key])]
            rows.append(title + "  " + " ".join(bases))
        self.local.setText('<pre style="font-family:Consolas,monospace">' + "\n".join(rows) + "</pre>")
        if not data.get("mutant"):
            note = "Référence seule : comparaison indisponible."
        elif not detail["reference"] or not detail["mutant"]:
            note = "Position absente d'une séquence : comparaison indisponible."
        elif len(data["reference"]["sequence"]) != len(data["mutant"]["sequence"]):
            note = "Longueurs différentes : comparaison par position, sans alignement."
        elif not detail["changes"]:
            note = "Base non appariée dans les deux prédictions."
        else:
            note = "Paires prédites ; positions depuis 1."
        self.note.setText(note + " Scores heuristiques, sans unité.")
