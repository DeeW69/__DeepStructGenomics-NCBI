"""Four-page desktop shell; scientific work executes in a cancellable process."""
import json
from pathlib import Path
import sys

from PySide6.QtCore import QProcess, QProcessEnvironment, Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QButtonGroup, QFileDialog, QFrame, QHBoxLayout, QListWidget, QListWidgetItem,
    QMainWindow, QPushButton, QScrollArea, QStackedWidget, QVBoxLayout, QWidget,
)
from .comparison import ComparisonPage
from .hic import HicPage
from .sequences import SequencesPage
from .ncbi import NcbiPage
from .services import RecentStore, read_rna
from .widgets import Metrics, StateCard, button, label, page


class MainWindow(QMainWindow):
    def __init__(self, workspace, alignment_config=None, sequence_inputs=None):
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
        sidebar_layout.addWidget(label("v0.5.0-rc1 · en préparation\nVersion candidate desktop"))
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
        self.recent_button = button("Ouvrir l'analyse sélectionnée", self.open_selected_recent)
        home.addWidget(self.recent_button)
        self.empty_history = StateCard("Aucune analyse récente", "Commencez par une séquence ou essayez la démo hors ligne.",
                                       actions=[("Essayer la démo C12A", self.start_demo)])
        home.addWidget(self.empty_history)
        self.latest_rna = None
        self.latest_title = label("Dernière analyse ARN consultée", "subheading")
        home.addWidget(self.latest_title)
        self.latest_metrics = Metrics()
        home.addWidget(self.latest_metrics)
        self.latest_button = button("Ouvrir cette analyse ARN", self.open_latest_rna)
        self.latest_button.setEnabled(False)
        home.addWidget(self.latest_button)
        home.addWidget(label(f"Résultats enregistrés dans : {self.workspace}"))
        self.sequences = SequencesPage(alignment_config, sequence_inputs)
        self.sequences.failed.connect(self.show_error)
        self.comparison = ComparisonPage()
        self.hic = HicPage()
        self.ncbi = NcbiPage()
        for index, widget in enumerate((self.dashboard, self.sequences, self.comparison, self.hic, self.ncbi)):
            if index in (0, 1):
                scroll = QScrollArea()
                scroll.setWidgetResizable(True)
                scroll.setFrameShape(QFrame.Shape.NoFrame)
                scroll.setWidget(widget)
                self.pages.addWidget(scroll)
            else:
                self.pages.addWidget(widget)
        self.sequences.requested.connect(lambda parameters: self.start_job("rna", parameters))
        self.hic.requested.connect(lambda parameters: self.start_job("hic", parameters))
        self.sequences.search_requested.connect(self.open_ncbi)
        self.comparison.search_requested.connect(self.open_ncbi)
        self.comparison.new_analysis.connect(lambda: self.open_sequences(0))
        self.ncbi.requested.connect(self.start_ncbi_job)
        self.ncbi.chosen.connect(self.choose_ncbi_reference)
        self.comparison.failed.connect(self.show_error)
        self.hic.failed.connect(self.show_error)
        body.addWidget(self.pages, 1)
        self.notice = StateCard()
        self.notice.hide()
        self.message = self.notice.body
        body.addWidget(self.notice)
        self.operation = StateCard()
        self.operation.hide()
        self.progress = self.operation.progress
        self.cancel_button = self.operation.controls[0]
        body.addWidget(self.operation)
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
        self.empty_history.setVisible(not entries)
        self.recents.setVisible(bool(entries))
        self.recent_button.setVisible(bool(entries))
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
                ("Paires MUT", metrics["mut_pairs"]), ("Absentes MUT" if metrics["deleted_pairs"] else "Perdues", metrics["lost"]),
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
        self.notice.set_state("success", "Référence WT sélectionnée", f"{preview['accession']} · {len(preview['sequence'])} nt. Ajoutez un mutant ou lancez l'analyse.",
                              [("Fermer", self.notice.hide)])

    def open_latest_rna(self):
        if self.latest_rna:
            self.open_rna(self.latest_rna["path"])

    def start_demo(self):
        self.sequences.fill_demo()
        self.navigate(1)
        self.sequences.submit()

    def start_job(self, kind, parameters):
        if self.process is not None:
            self.show_error("Un calcul est déjà en cours. Attendez sa fin ou annulez-le.")
            return
        self.notice.hide()
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
        self.sequences.setEnabled(False)
        self.ncbi.setEnabled(False)
        self.hic.form_widget.setEnabled(False)
        title = {"ncbi_search": "Recherche NCBI", "ncbi_preview": "Aperçu NCBI"}.get(kind, f"Analyse {kind.upper()}")
        if kind == "rna":
            source = parameters.get("accession") or parameters.get("sequence_label") or parameters.get("fasta_path") or "Séquence saisie"
            length = len("".join(parameters.get("sequence", "").split()))
            context = f"{source}" + (f" · {length} nt" if length else "")
            detail = "Validation, prédiction et export des résultats."
        elif kind == "hic":
            context = f"{parameters.get('chromosome', '')} · {parameters.get('start', 0)}–{parameters.get('end', 0)} bp"
            detail = "Analyse des contacts et export des résultats."
        else:
            context = parameters.get("term") or parameters.get("accession", "NCBI")
            detail = "Interrogation du service NCBI." if kind == "ncbi_search" else "Téléchargement et validation de la séquence."
        self.operation.set_state("busy", title + " en cours", f"{context}\n{detail} Durée restante non estimée.",
                                 [("Annuler", self.cancel_job)])
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
            self.show_job_error("Impossible de démarrer le processus d'analyse.", self.active_kind)
            if self.active_kind.startswith("ncbi_"):
                self.ncbi.state_card.set_state("error", "Démarrage impossible", "Réessayez la recherche ou l'aperçu.")
            self.finish_process()

    def job_finished(self, code, _status):
        self.read_stdout()
        self.read_stderr()
        was_cancelled = self.cancelled
        kind = self.active_kind
        self.finish_process()
        if was_cancelled:
            if kind.startswith("ncbi_"):
                self.ncbi.state_card.set_state("warning", "Opération annulée", "Vous pouvez relancer la recherche ou l'aperçu.")
            self.statusBar().showMessage("Calcul annulé. Aucun résultat ajouté à l'historique.")
            self.notice.set_state("warning", "Opération annulée", "Aucun résultat ajouté à l'historique. Les résultats précédents restent consultables.",
                                  [("Revenir aux paramètres", lambda: self.navigate(self.source_page(kind))), ("Fermer", self.notice.hide)])
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
                if self.notice.state != "error" or self.notice.isHidden():
                    self.notice.set_state("success", "Analyse terminée", "Résultats enregistrés. Explorez les vues ou exportez une figure.",
                                          [("Fermer", self.notice.hide)])
            self.statusBar().showMessage("Opération terminée" if kind.startswith("ncbi_") else "Analyse terminée · résultats enregistrés")
        except (ValueError, OSError, KeyError, TypeError) as exc:
            self.show_job_error(str(exc) if self.stdout else "Le processus d'analyse s'est arrêté sans résultat.", kind)
            if kind.startswith("ncbi_"):
                self.ncbi.state_card.set_state("error", "Aucun résultat utilisable", "Consultez l'erreur ci-dessous, puis réessayez.")

    def finish_process(self):
        if self.process:
            self.process.deleteLater()
            self.process = None
        self.operation.hide()
        self.sequences.setEnabled(True)
        self.ncbi.setEnabled(True)
        self.hic.form_widget.setEnabled(True)

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
        self.notice.hide()
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
        self.notice.set_state("error", "Action non terminée", text, [("Fermer", self.notice.hide)])
        self.statusBar().showMessage("Action non terminée · voir le message")

    @staticmethod
    def source_page(kind):
        return 4 if kind.startswith("ncbi_") else 3 if kind == "hic" else 1

    def show_job_error(self, text, kind):
        if kind == "rna" and "Alignement trop coûteux" in text:
            def edit_limit():
                self.navigate(1)
                self.sequences.advanced.expand()
                self.notice.hide()
            def positional():
                self.sequences.align_mutant.setChecked(False)
                edit_limit()
            self.notice.set_state("error", "Limite d'alignement dépassée", text,
                [("Utiliser comparaison positionnelle", positional), ("Modifier la limite", edit_limit), ("Annuler", self.notice.hide)])
            return
        title = "Recherche interrompue" if kind.startswith("ncbi_") else "Analyse interrompue"
        action = "Modifier la recherche" if kind == "ncbi_search" else "Choisir une autre notice" if kind == "ncbi_preview" else "Corriger la séquence" if kind == "rna" else "Corriger les paramètres"
        self.notice.set_state("error", title, text,
                              [(action, lambda: self.navigate(self.source_page(kind))), ("Fermer", self.notice.hide)])
        self.statusBar().showMessage(title + " · corrigez les paramètres puis relancez")

    def closeEvent(self, event):
        if self.process:
            self.show_error("Annulez le calcul en cours ou attendez sa fin avant de fermer.")
            event.ignore()
            return
        self.comparison.shutdown()
        event.accept()
