"""Shared presentation widgets; figures are supplied by the analysis views."""
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame, QHBoxLayout, QLabel, QProgressBar, QPushButton, QVBoxLayout, QWidget,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT


def label(text, name=None):
    widget = QLabel(text)
    widget.setTextFormat(Qt.TextFormat.PlainText)
    widget.setWordWrap(True)
    if name:
        widget.setObjectName(name)
    return widget


def button(text, callback, primary=False):
    widget = QPushButton(text)
    if primary:
        widget.setObjectName("primary")
    widget.clicked.connect(callback)
    return widget


def page(title, subtitle):
    widget = QWidget()
    layout = QVBoxLayout(widget)
    layout.setContentsMargins(26, 22, 26, 20)
    layout.setSpacing(14)
    layout.addWidget(label(title, "heading"))
    layout.addWidget(label(subtitle, "subheading"))
    return widget, layout


class Metrics(QWidget):
    def __init__(self):
        super().__init__()
        self.row = QHBoxLayout(self)
        self.row.setContentsMargins(0, 0, 0, 0)

    def set_values(self, values):
        while self.row.count():
            old = self.row.takeAt(0).widget()
            old.hide()
            old.deleteLater()
        for title, value in values:
            card = QFrame()
            card.setObjectName("card")
            layout = QVBoxLayout(card)
            layout.addWidget(label("—" if value is None else str(value), "metric"))
            layout.addWidget(label(title))
            self.row.addWidget(card)


class StateCard(QFrame):
    """One vocabulary for empty, busy, success, warning and error states."""
    def __init__(self, title="", text="", state="empty", actions=()):
        super().__init__()
        self.setObjectName("stateCard")
        self._initialized = False
        layout = QVBoxLayout(self)
        layout.setContentsMargins(14, 10, 14, 10)
        layout.setSpacing(6)
        self.title = label("", "stateTitle")
        self.body = label("")
        layout.addWidget(self.title)
        layout.addWidget(self.body)
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setAccessibleName("Opération en cours, progression non mesurée")
        layout.addWidget(self.progress)
        self.action_row = QWidget()
        self.actions = QHBoxLayout(self.action_row)
        self.actions.setContentsMargins(0, 0, 0, 0)
        self.controls = []
        for i in range(3):
            control = button("", lambda checked=False, index=i: self.trigger(index), i == 0)
            self.controls.append(control)
            self.actions.addWidget(control)
        self.actions.addStretch()
        layout.addWidget(self.action_row)
        self.set_state(state, title, text, actions)
        self._initialized = True

    def trigger(self, index):
        if index < len(self.callbacks):
            self.callbacks[index]()

    def set_state(self, state, title, text, actions=()):
        self.state = state
        self.setProperty("state", state)
        self.style().unpolish(self)
        self.style().polish(self)
        self.title.setText(title)
        self.body.setText(text)
        self.body.setVisible(bool(text))
        self.progress.setVisible(state == "busy")
        self.callbacks = [callback for _, callback in actions]
        self.action_row.setVisible(bool(actions))
        for i, control in enumerate(self.controls):
            control.setVisible(i < len(actions))
            if i < len(actions):
                control.setText(actions[i][0])
        if self._initialized:
            self.show()


class FigureView(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.canvas = None
        self.figure = None
        self.empty = StateCard("Aucune analyse", "Lancez une analyse ou ouvrez un résultat enregistré.")
        self.layout.addWidget(self.empty)

    def set_figure(self, figure):
        old_figure = self.figure
        while self.layout.count():
            old = self.layout.takeAt(0).widget()
            old.hide()
            old.deleteLater()
        if old_figure:
            from matplotlib import pyplot as plt
            plt.close(old_figure)
        self.figure = figure
        self.canvas = FigureCanvasQTAgg(figure)
        self.layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        self.layout.addWidget(self.canvas, 1)
        self.canvas.draw_idle()
