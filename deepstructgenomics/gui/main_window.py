"""Four-page desktop shell; scientific work executes in a cancellable process."""
import json
from pathlib import Path
import sys

from PySide6.QtCore import QProcess, QProcessEnvironment, Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QButtonGroup, QFileDialog, QFrame, QHBoxLayout, QListWidget, QListWidgetItem,
    QMainWindow, QProgressBar, QPushButton, QStackedWidget, QVBoxLayout, QWidget,
)
from .comparison import ComparisonPage
from .hic import HicPage
from .sequences import SequencesPage
from .ncbi import NcbiPage
from .services import RecentStore, read_rna
from .widgets import Metrics, button, label, page


class MainWindow(QMainWindow):
    def __init__(self, workspace):
        super().__init__()
        self.workspace = Path(workspace).expanduser().resolve()
        self.history = RecentStore(self.workspace / "recent.json")
        self.current_directory = None
        self.result_directories = {}
        self.process = None
        self.cancelled = False
        self.stdout = bytearray()
        self.stderr = bytearray()
        self.setWindowTitle("DeepStructGenomics · RNA & genomic interactions")
        self.resize(1420, 920)
        self.setMinimumSize(1100, 760)
        central = QWidget()
        row = QHBoxLayout(central)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(0)
        sidebar = QFrame()
        sidebar.setObjectName("sidebar")
        sidebar.setFixedWidth(225)
        sidebar_layout = QVBoxLayout(sidebar)
        sidebar_layout.setContentsMargins(14, 27, 14, 20)
        sidebar_layout.addWidget(label("DeepStruct\nGenomics", "brand"))
        sidebar_layout.addWidget(label("RNA structure & genomic\ninteraction analysis workspace"))
        sidebar_layout.addSpacing(35)
        self.navigation = QButtonGroup(self)
        for index, title in enumerate(("01   Accueil", "02   Séquences ARN", "03   Comparaison", "04   Hi-C Explorer")):
            nav = QPushButton(title)
            nav.setCheckable(True)
            self.navigation.addButton(nav, index)
            sidebar_layout.addWidget(nav)
        sidebar_layout.addStretch()
        sidebar_layout.addWidget(label("Moteur scientifique v0.4.6\nInterface en développement"))
        sidebar_layout.addWidget(label("DeeW69"))
        row.addWidget(sidebar)
        body = QVBoxLayout()
        body.setContentsMargins(0, 0, 0, 0)
        toolbar = QHBoxLayout()
        toolbar.setContentsMargins(26, 14, 26, 0)
        toolbar.addWidget(label("ESPACE DE RECHERCHE"))
        toolbar.addStretch()
        toolbar.addWidget(button("Ouvrir un résultat…", self.choose_result))
        self.folder_button = button("Dossier des résultats", self.open_directory)
        self.folder_button.setEnabled(False)
        toolbar.addWidget(self.folder_button)
        body.addLayout(toolbar)
        self.pages = QStackedWidget()
        self.dashboard, home = page("Votre espace d'analyse", "Des séquences aux appariements, des contacts aux interactions génomiques.")
        hero = QFrame()
        hero.setObjectName("card")
        hero_layout = QVBoxLayout(hero)
        hero_layout.setContentsMargins(24, 24, 24, 24)
        hero_layout.addWidget(label("Explorer une nouvelle hypothèse", "heading"))
        hero_layout.addWidget(label("Importez votre ARN ou vos contacts Hi-C. Retrouvez les résultats, leurs paramètres et les exports dans un même espace."))
        actions = QHBoxLayout()
        actions.addWidget(button("Rechercher NCBI", self.open_ncbi, True))
        actions.addWidget(button("Importer FASTA", self.import_fasta))
        actions.addWidget(button("Saisie directe", lambda: self.open_sequences(0)))
        hero_layout.addLayout(actions)
        secondary_actions = QHBoxLayout()
        secondary_actions.addWidget(button("Essayer la démo C12A", self.start_demo))
        secondary_actions.addWidget(button("Explorer Hi-C", lambda: self.navigate(3)))
        secondary_actions.addStretch()
        hero_layout.addLayout(secondary_actions)
        home.addWidget(hero)
        home.addWidget(label("Analyses récentes", "subheading"))
        self.recents = QListWidget()
        self.recents.itemDoubleClicked.connect(self.open_recent)
        home.addWidget(self.recents, 1)
        home.addWidget(button("Ouvrir l'analyse sélectionnée", self.open_selected_recent))
        self.latest_rna = None
        self.latest_title = label("Dernière analyse ARN consultée", "subheading")
        home.addWidget(self.latest_title)
        self.latest_metrics = Metrics()
        home.addWidget(self.latest_metrics)
        self.latest_button = button("Ouvrir cette analyse ARN", self.open_latest_rna)
        self.latest_button.setEnabled(False)
        home.addWidget(self.latest_button)
        home.addWidget(label(f"Résultats enregistrés dans : {self.workspace}"))
        self.sequences = SequencesPage()
        self.comparison = ComparisonPage()
        self.hic = HicPage()
        self.ncbi = NcbiPage()
        for widget in (self.dashboard, self.sequences, self.comparison, self.hic, self.ncbi):
            self.pages.addWidget(widget)
        self.sequences.requested.connect(lambda parameters: self.start_job("rna", parameters))
        self.hic.requested.connect(lambda parameters: self.start_job("hic", parameters))
        self.sequences.search_requested.connect(self.open_ncbi)
        self.ncbi.requested.connect(self.start_ncbi_job)
        self.ncbi.chosen.connect(self.choose_ncbi_reference)
        self.comparison.failed.connect(self.show_error)
        self.hic.failed.connect(self.show_error)
        body.addWidget(self.pages, 1)
        self.message = label("")
        self.message.setContentsMargins(26, 4, 26, 4)
        self.message.hide()
        body.addWidget(self.message)
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.hide()
        body.addWidget(self.progress)
        self.cancel_button = button("Annuler le calcul", self.cancel_job)
        self.cancel_button.hide()
        body.addWidget(self.cancel_button)
        row.addLayout(body, 1)
        self.setCentralWidget(central)
        self.navigation.idClicked.connect(self.navigate)
        self.navigate(0)
        self.refresh_history()
        self.statusBar().showMessage("Prêt · calculs locaux, réseau uniquement pour NCBI")

    def navigate(self, index):
        self.pages.setCurrentIndex(index)
        self.navigation.button(1 if index == 4 else index).setChecked(True)
        if index in (2, 3):
            self.current_directory = self.result_directories.get("rna" if index == 2 else "hic")
            self.folder_button.setEnabled(self.current_directory is not None)

    def refresh_history(self):
        self.recents.clear()
        entries = self.history.load()
        for entry in entries:
            kind = "ARN" if entry["kind"] == "rna" else "Hi-C"
            item = QListWidgetItem(f"{entry['label']}  ·  {kind}\n{entry.get('summary', 'Résultat enregistré · ouvrir pour consulter')}")
            item.setToolTip(entry["path"])
            item.setData(Qt.ItemDataRole.UserRole, entry)
            self.recents.addItem(item)
        if not self.recents.count():
            item = QListWidgetItem("Aucune analyse récente. Lancez la démo ou ouvrez un résultat existant.")
            item.setFlags(Qt.ItemFlag.NoItemFlags)
            self.recents.addItem(item)
        if self.history.warning:
            self.show_error(self.history.warning)
        self.latest_rna = None
        self.latest_button.setEnabled(False)
        self.latest_metrics.set_values([])
        self.latest_title.setText("Aucune analyse ARN consultable pour le moment")
        for entry in entries:
            if entry["kind"] != "rna":
                continue
            try:
                data, metrics = read_rna(entry["path"])
            except (ValueError, OSError, KeyError, TypeError):
                continue
            self.latest_rna = entry
            self.latest_title.setText(f"Dernière analyse ARN consultée · {entry['label']}")
            self.latest_metrics.set_values([
                ("Bases WT", len(data["reference"]["sequence"])), ("Paires WT", metrics["wt_pairs"]),
                ("Paires MUT", metrics["mut_pairs"]), ("Perdues", metrics["lost"]),
                ("Nouvelles", metrics["gained"]), ("|Δ| maximal", f"{metrics['max_abs_delta']:.3f}" if metrics["max_abs_delta"] is not None else None)])
            self.latest_button.setEnabled(True)
            break

    def open_sequences(self, source):
        self.sequences.source.setCurrentIndex(source)
        self.navigate(1)

    def import_fasta(self):
        self.open_sequences(1)
        self.sequences.choose_file(self.sequences.fasta)

    def open_ncbi(self):
        self.navigate(4)
        self.ncbi.query.setFocus()

    def start_ncbi_job(self, kind, parameters):
        self.start_job(kind, {**parameters, "ncbi_email": self.sequences.email.text() or None})

    def choose_ncbi_reference(self, preview):
        self.sequences.use_ncbi(preview)
        self.navigate(1)
        self.statusBar().showMessage("Référence WT sélectionnée · ajoutez un mutant ou lancez l'analyse")

    def open_latest_rna(self):
        if self.latest_rna:
            self.open_rna(self.latest_rna["path"])

    def start_demo(self):
        self.sequences.fill_demo()
        self.navigate(1)
        self.start_job("rna", self.sequences.parameters())

    def start_job(self, kind, parameters):
        if self.process is not None:
            self.show_error("Un calcul est déjà en cours. Attendez sa fin ou annulez-le.")
            return
        self.message.hide()
        self.stdout.clear()
        self.stderr.clear()
        self.cancelled = False
        self.active_kind = kind
        process = QProcess(self)
        self.process = process
        # Works both from the checkout and from an installed wheel, independent of cwd.
        environment = QProcessEnvironment.systemEnvironment()
        environment.insert("PYTHONPATH", str(Path(__file__).resolve().parents[2]))
        environment.insert("PYTHONUTF8", "1")
        process.setProcessEnvironment(environment)
        process.setProgram(sys.executable)
        process.setArguments(["-m", "deepstructgenomics.gui.worker"])
        process.readyReadStandardOutput.connect(self.read_stdout)
        process.readyReadStandardError.connect(self.read_stderr)
        payload = json.dumps({"kind": kind, "workspace": str(self.workspace), "parameters": parameters}).encode("utf-8")
        process.started.connect(lambda: (process.write(payload), process.closeWriteChannel()))
        process.finished.connect(self.job_finished)
        process.errorOccurred.connect(self.process_error)
        self.progress.show()
        self.cancel_button.show()
        self.sequences.setEnabled(False)
        self.ncbi.setEnabled(False)
        title = {"ncbi_search": "Recherche NCBI", "ncbi_preview": "Aperçu NCBI"}.get(kind, f"Analyse {kind.upper()}")
        self.statusBar().showMessage(f"{title} en cours… Vous pouvez parcourir les résultats précédents.")
        process.start()

    def read_stdout(self):
        if self.process:
            self.stdout.extend(bytes(self.process.readAllStandardOutput()))

    def read_stderr(self):
        if self.process:
            self.stderr.extend(bytes(self.process.readAllStandardError()))

    def process_error(self, error):
        if error == QProcess.ProcessError.FailedToStart:
            self.show_error("Impossible de démarrer le processus d'analyse.")
            if self.active_kind.startswith("ncbi_"):
                self.ncbi.note.setText("Démarrage impossible. Réessayez.")
            self.finish_process()

    def job_finished(self, code, _status):
        self.read_stdout()
        self.read_stderr()
        was_cancelled = self.cancelled
        kind = self.active_kind
        self.finish_process()
        if was_cancelled:
            if kind.startswith("ncbi_"):
                self.ncbi.note.setText("Opération annulée. Vous pouvez relancer la recherche ou l'aperçu.")
            self.statusBar().showMessage("Calcul annulé. Aucun résultat ajouté à l'historique.")
            return
        try:
            result = json.loads(self.stdout.decode("utf-8"))
            if code or "error" in result:
                raise ValueError(result.get("error", "Le calcul n'a pas abouti."))
            if kind == "ncbi_search":
                self.ncbi.show_results(result)
            elif kind == "ncbi_preview":
                self.ncbi.show_preview(result)
            else:
                self.open_entry(result)
            self.statusBar().showMessage("Opération terminée" if kind.startswith("ncbi_") else "Analyse terminée · résultats enregistrés")
        except (ValueError, OSError, KeyError, TypeError) as exc:
            self.show_error(str(exc) if self.stdout else "Le processus d'analyse s'est arrêté sans résultat.")
            if kind.startswith("ncbi_"):
                self.ncbi.note.setText("NCBI indisponible ou réponse invalide. Vérifiez les critères puis réessayez.")

    def finish_process(self):
        if self.process:
            self.process.deleteLater()
            self.process = None
        self.progress.hide()
        self.cancel_button.hide()
        self.sequences.setEnabled(True)
        self.ncbi.setEnabled(True)

    def cancel_job(self):
        if self.process:
            self.cancelled = True
            self.process.kill()

    def choose_result(self):
        path, _ = QFileDialog.getOpenFileName(self, "Manifeste ARN ou rapport Hi-C", str(self.workspace), "JSON (*.json)")
        if path:
            try:
                with Path(path).open(encoding="utf-8") as stream:
                    data = json.load(stream)
                if isinstance(data, dict) and "balanced_matrix" in data:
                    self.open_entry({"kind": "hic", "path": path, "directory": str(Path(path).parent),
                                     "label": data.get("region", {}).get("chromosome", "Hi-C")})
                else:
                    self.open_rna(path)
            except (ValueError, OSError, KeyError, TypeError) as exc:
                self.show_error(f"Résultat illisible : {exc}")

    def open_rna(self, path):
        path = Path(path).resolve()
        report_root = path.parent.parent.parent
        directory = report_root if (report_root / f"{path.parent.parent.name}.json").is_file() else path.parent
        try:
            self.open_entry({"kind": "rna", "path": str(path), "label": path.parent.parent.name,
                             "directory": str(directory)})
        except (ValueError, OSError, KeyError, TypeError) as exc:
            self.show_error(f"Résultat ARN illisible : {exc}")

    def open_entry(self, entry):
        if entry["kind"] == "rna":
            self.comparison.load(entry["path"])
            entry = {**entry, "label": str(self.comparison.data.get("identifier", "ARN")),
                     "summary": self.comparison.summary}
            self.navigate(2)
        else:
            self.hic.load(entry["path"])
            report = self.hic.report
            region = report["region"]
            entry = {**entry, "summary": f"{region['start']}–{region['end']} bp · {region['bins']} bins · "
                     + ("équilibrage convergé" if report["status"] == "ok" else "équilibrage non convergé")}
            self.navigate(3)
        self.current_directory = Path(entry["directory"])
        self.result_directories[entry["kind"]] = self.current_directory
        self.folder_button.setEnabled(True)
        self.message.hide()
        try:
            self.history.add(entry)
        except OSError:
            self.show_error("Résultat ouvert, mais l'historique ne peut pas être enregistré.")
        self.refresh_history()

    def open_selected_recent(self):
        item = self.recents.currentItem()
        if item:
            self.open_recent(item)

    def open_recent(self, item):
        entry = item.data(Qt.ItemDataRole.UserRole)
        if entry:
            try:
                self.open_entry(entry)
            except (ValueError, OSError, KeyError, TypeError) as exc:
                self.show_error(f"Résultat indisponible : {exc}")

    def open_directory(self):
        if self.current_directory:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.current_directory)))

    def show_error(self, text):
        self.message.setText(text)
        self.message.setStyleSheet("color: #9c3428; background: #fff0e9; padding: 10px;")
        self.message.show()
        self.statusBar().showMessage("Action non terminée · voir le message")

    def closeEvent(self, event):
        if self.process:
            self.show_error("Annulez le calcul en cours ou attendez sa fin avant de fermer.")
            event.ignore()
            return
        self.comparison.shutdown()
        event.accept()
