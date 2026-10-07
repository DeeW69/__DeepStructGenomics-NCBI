"""Entry point for the optional Qt desktop application."""
import argparse
from pathlib import Path
import sys


def apply_theme(app):
    from PySide6.QtGui import QColor, QPalette
    app.setStyle("Fusion")
    palette = QPalette()
    for role, color in ((QPalette.ColorRole.Window, "#f3f6fa"),
                        (QPalette.ColorRole.Base, "#ffffff"),
                        (QPalette.ColorRole.WindowText, "#20374c"),
                        (QPalette.ColorRole.Text, "#20374c"),
                        (QPalette.ColorRole.Button, "#ffffff"),
                        (QPalette.ColorRole.ButtonText, "#20374c"),
                        (QPalette.ColorRole.Highlight, "#127c82"),
                        (QPalette.ColorRole.HighlightedText, "#ffffff")):
        palette.setColor(role, QColor(color))
    app.setPalette(palette)
    theme = Path(__file__).parent / "styles" / "theme.qss"
    app.setStyleSheet(theme.read_text(encoding="utf-8"))


def main(argv=None):
    parser = argparse.ArgumentParser(description="DeepStructGenomics — espace d'analyse desktop")
    parser.add_argument("--workspace", type=Path, default=Path("outputs_gui"))
    parser.add_argument("--manifest", type=Path, help="Ouvrir un résultat ARN existant.")
    args = parser.parse_args(argv)
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        parser.exit(2, "Interface facultative absente : installer requirements-gui-lock.txt.\n")
    from .main_window import MainWindow
    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName("DeepStructGenomics")
    app.setOrganizationName("DeeW69")
    apply_theme(app)
    window = MainWindow(args.workspace)
    window.show()
    if args.manifest:
        window.open_rna(args.manifest)
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
