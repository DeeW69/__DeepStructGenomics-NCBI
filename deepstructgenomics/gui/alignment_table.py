"""Virtual table: only visible cells are materialized by Qt."""
from functools import lru_cache
from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt
from PySide6.QtGui import QColor
from .services import position_details


class AlignmentTableModel(QAbstractTableModel):
    headers = ["Colonne", "Pos. WT", "WT", "Pos. MUT", "MUT", "Variation", "Score WT", "Score MUT", "Δ", "Contexte WT", "Contexte MUT"]

    def __init__(self, data, mapping, contexts, parent=None):
        super().__init__(parent)
        self.source, self.mapping, self.contexts = data, mapping, contexts
        self.headers = list(self.headers)
        if not mapping.alignment:
            self.headers[0] = "Position"
        self.row_values = lru_cache(maxsize=256)(self._row_values)

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.mapping.columns)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.headers)

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Horizontal:
            return self.headers[section]

    def _row_values(self, row):
        detail = position_details(self.source, row + 1)
        wt, mut = detail["reference"], detail["mutant"]
        operation = {"match": "Identique", "substitution": "Substitution", "insertion": "Insertion",
                     "deletion": "Délétion", "positional": "Sans alignement"}[detail["operation"]]
        return [str(row + 1), str(wt["position"]) if wt else "—", wt["base"] if wt else "—",
                str(mut["position"]) if mut else "—", mut["base"] if mut else "—", operation,
                f"{wt['score']:.3f}" if wt else "—", f"{mut['score']:.3f}" if mut else "—",
                f"{detail['delta']:+.3f}" if detail["delta"] is not None else "—",
                self.contexts[0][wt["position"] - 1] if wt else "—",
                self.contexts[1][mut["position"] - 1] if mut else "—"]

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        if role == Qt.ItemDataRole.DisplayRole:
            return self.row_values(index.row())[index.column()]
        if role == Qt.ItemDataRole.BackgroundRole:
            color = {"substitution": "#fff2dc", "insertion": "#f2eafa", "deletion": "#e7f1fb"}.get(self.mapping.columns[index.row()].operation)
            return QColor(color) if color else None
