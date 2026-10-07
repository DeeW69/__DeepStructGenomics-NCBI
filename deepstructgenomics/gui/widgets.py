"""Shared presentation widgets; figures are supplied by the analysis views."""
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget,
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
            self.row.takeAt(0).widget().deleteLater()
        for title, value in values:
            card = QFrame()
            card.setObjectName("card")
            layout = QVBoxLayout(card)
            layout.addWidget(label("—" if value is None else str(value), "metric"))
            layout.addWidget(label(title))
            self.row.addWidget(card)


class FigureView(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.canvas = None
        self.figure = None
        self.layout.addWidget(label("Les résultats de votre analyse apparaîtront ici."))

    def set_figure(self, figure):
        old_figure = self.figure
        while self.layout.count():
            self.layout.takeAt(0).widget().deleteLater()
        if old_figure:
            from matplotlib import pyplot as plt
            plt.close(old_figure)
        self.figure = figure
        self.canvas = FigureCanvasQTAgg(figure)
        self.layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        self.layout.addWidget(self.canvas, 1)
        self.canvas.draw_idle()
