"""Offscreen UI integration: real pipeline, errors, selection and cancellation."""
import os
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from pathlib import Path
import time

import pytest
pytest.importorskip("PySide6")
from PySide6.QtCore import QProcess
from PySide6.QtWidgets import QApplication
from deepstructgenomics.gui.main_window import MainWindow
from deepstructgenomics.gui.services import run_rna


@pytest.fixture(scope="module")
def app():
    application = QApplication.instance() or QApplication([])
    from deepstructgenomics.gui.app import apply_theme
    apply_theme(application)
    return application


@pytest.fixture
def window(app, tmp_path):
    widget = MainWindow(tmp_path)
    yield widget
    if widget.process:
        widget.cancel_job()
        wait_for(app, lambda: widget.process is None)
    widget.close()
    widget.deleteLater()
    app.processEvents()


def wait_for(app, condition, timeout=30):
    deadline = time.monotonic() + timeout
    while not condition() and time.monotonic() < deadline:
        app.processEvents()
        time.sleep(.01)
    app.processEvents()
    assert condition(), "Desktop operation timed out"


def test_demo_process_to_comparison_and_recent_reopen(app, window):
    window.start_demo()
    assert window.process is not None
    wait_for(app, lambda: window.process is None)
    assert window.pages.currentIndex() == 2, window.message.text()
    assert window.comparison.table.rowCount() == 12
    assert window.comparison.table.item(11, 5).text() == "-0.700"
    assert window.comparison.tabs.currentIndex() == 0
    assert window.comparison.tabs.tabText(0) == "Structure secondaire"
    assert window.comparison.table.item(4, 6).text() == "Boucle terminale"
    window.comparison.position.setValue(12)
    inspector = window.comparison.inspector
    assert inspector.fields["reference"]["base"].text() == "C"
    assert inspector.fields["reference"]["partner"].text() == "G1"
    assert inspector.fields["mutant"]["partner"].text() == "—"
    assert inspector.badges[0].title.text() == "Paire perdue"
    assert inspector.badges[1].title.text() == "Substitution C → A"
    assert inspector.delta.text() == "Δ  -0.700"
    figure = window.comparison.figure_view.figure
    for axis, coordinates in figure.rna_coordinates.items():
        assert figure.rna_selections[axis].get_offsets()[0].tolist() == coordinates[11].tolist()
    window.navigate(0)
    window.open_recent(window.recents.item(0))
    assert window.pages.currentIndex() == 2
    assert len(window.history.load()) == 1
    assert window.latest_rna is not None and window.latest_button.isEnabled()
    window.open_latest_rna()
    assert window.pages.currentIndex() == 2


def test_ncbi_search_preview_selection_and_invalidation(app, window, monkeypatch):
    calls = []
    monkeypatch.setattr(window, "start_job", lambda kind, params: calls.append((kind, params)))
    window.open_ncbi()
    page = window.ncbi
    page.query.setText("HOTAIR")
    page.search()
    assert calls[-1][0] == "ncbi_search" and calls[-1][1]["start"] == 0
    result = {"term": "HOTAIR", "query": "HOTAIR", "start": 0, "page_size": 20, "total": 21,
              "rows": [{"accession": "NR_TEST.1", "organism": "Homo sapiens", "molecule": "ncRNA",
                        "length": 12, "description": "Reference"}]}
    page.show_results(result)
    assert page.next.isEnabled() and not page.previous.isEnabled()
    page.table.selectRow(0)
    assert page.preview_button.isEnabled() and not page.use.isEnabled()
    page.fetch_preview()
    assert calls[-1][0] == "ncbi_preview"
    preview = {"accession": "NR_TEST.1", "description": "Reference", "sequence": "GGGGAAAACCCC", "sha256": "digest"}
    page.show_preview({**preview, "accession": "wrong"})
    assert not page.use.isEnabled()
    page.show_preview(preview)
    assert page.state_card.state == "success"
    page.use_sequence()
    assert window.pages.currentIndex() == 1 and window.sequences.source.currentIndex() == 2
    assert window.sequences.parameters()["expected_ncbi_sha256"] == "digest"
    assert window.history.load() == []  # Browsing is not an analysis.
    window.sequences.source.setCurrentIndex(0)
    assert window.sequences.reference_note.isHidden()
    window.sequences.accession.setText("OTHER.1")
    assert window.sequences.ncbi_preview is None
    page.search(1)
    assert calls[-1][1]["start"] == 20 and calls[-1][1]["term"] == "HOTAIR"
    page.invalidate()
    assert page.table.rowCount() == 0 and page.preview is None and not page.next.isEnabled()
    assert page.state_card.state == "empty"


def test_empty_states_guide_to_inputs_and_inspector_resets(app, window, tmp_path):
    assert not window.empty_history.isHidden()
    assert window.comparison.splitter.isHidden()
    assert not window.comparison.export_button.isEnabled()
    assert not window.hic.mode.isEnabled()
    window.comparison.empty.controls[0].click()
    assert window.pages.currentIndex() == 4
    window.comparison.empty.controls[1].click()
    assert window.pages.currentIndex() == 1
    window.sequences.sequence.setPlainText("ACGU")
    assert window.sequences.empty.isHidden()
    entry = run_rna(tmp_path, {"sequence": "GGGGAAAACCCC", "mutant_sequence": "GGGGAAAACCCA"})
    window.open_entry(entry)
    window.comparison.position.setValue(12)
    inspector = window.comparison.inspector
    assert not inspector.badges[0].isHidden()
    entry = run_rna(tmp_path, {"sequence": "GGGGAAAACCCC"})
    window.open_entry(entry)
    assert all(badge.isHidden() for badge in inspector.badges)
    assert inspector.delta.text() == "Δ indisponible"
    assert inspector.fields["mutant"]["base"].text() == "Indisponible"
    assert window.comparison.empty.isHidden()
    assert window.empty_history.isHidden()


def test_ncbi_errors_use_recoverable_states(app, window):
    # Exercise worker completion routing without making a network call.
    import json
    window.active_kind = "ncbi_search"
    window.stdout.extend(json.dumps({"error": "NCBI temporairement indisponible"}).encode())
    window.job_finished(1, QProcess.ExitStatus.NormalExit)
    assert window.notice.state == "error" and window.ncbi.state_card.state == "error"
    window.notice.controls[0].click()
    assert window.pages.currentIndex() == 4 and window.ncbi.isEnabled()
    assert window.history.load() == []


def test_invalid_input_visible_and_retry_possible(app, window):
    window.start_job("rna", {"sequence": "ACG!"})
    assert window.operation.state == "busy" and not window.operation.isHidden()
    assert not window.hic.form_widget.isEnabled()
    wait_for(app, lambda: window.process is None)
    assert window.message.text()
    assert window.sequences.isEnabled()
    assert window.history.load() == []
    assert window.notice.state == "error" and window.operation.isHidden()
    assert window.notice.controls[0].text() == "Corriger la séquence"
    window.notice.controls[0].click()
    assert window.pages.currentIndex() == 1
    window.notice.controls[1].click()
    assert window.notice.isHidden()
    window.start_demo()
    wait_for(app, lambda: window.process is None)
    assert window.pages.currentIndex() == 2
    assert window.notice.state == "success"


def test_cancel_preserves_previous_result_and_close_waits(app, window):
    window.start_demo()
    wait_for(app, lambda: window.process is None)
    previous = window.comparison.manifest
    window.start_job("rna", {"sequence": "G" * 500 + "C" * 500})
    wait_for(app, lambda: window.process.state() == QProcess.ProcessState.Running)
    assert not window.close()
    window.cancel_job()
    wait_for(app, lambda: window.process is None)
    assert window.comparison.manifest == previous
    assert len(window.history.load()) == 1
    assert window.sequences.isEnabled()
    assert window.notice.state == "warning" and window.operation.isHidden()
    assert window.hic.form_widget.isEnabled()


def test_missing_result_and_no_mutant(app, window, tmp_path):
    window.open_rna(tmp_path / "missing.json")
    assert "illisible" in window.message.text()
    entry = run_rna(tmp_path, {"sequence": "GGGGAAAACCCC"})
    window.open_entry(entry)
    assert not window.comparison.vtk_button.isEnabled()
    assert window.comparison.table.item(0, 5).text() == "—"
    window.open_rna(entry["path"])
    assert window.current_directory == Path(entry["directory"])


def test_hic_process_to_heatmap_and_3d(app, window):
    pytest.importorskip("scipy")
    window.hic.fill_demo()
    assert Path(window.hic.path.text()).is_file()
    parameters = window.hic.parameters()
    parameters["permutations"] = 9
    window.start_job("hic", parameters)
    wait_for(app, lambda: window.process is None)
    assert window.pages.currentIndex() == 3, window.message.text()
    assert window.hic.report["status"] == "ok"
    assert window.hic.report["region"]["bins"] == 32
    window.hic.mode.setCurrentText("Reconstruction 3D")
    assert window.hic.figure_view.figure.axes[0].name == "3d"
    hic_directory = window.current_directory
    window.start_demo()
    wait_for(app, lambda: window.process is None)
    rna_directory = window.current_directory
    assert hic_directory != rna_directory
    window.navigate(3)
    assert window.current_directory == hic_directory
    window.navigate(2)
    assert window.current_directory == rna_directory


def test_export_figure_and_offscreen_vtk_fallback(app, window, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QFileDialog
    entry = run_rna(tmp_path, {"sequence": "GGGGAAAACCCC", "mutant_sequence": "GGGGAAAACCCA"})
    window.open_entry(entry)
    target = tmp_path / "comparison.png"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(target), "PNG"))
    window.comparison.export_png()
    assert target.read_bytes().startswith(b"\x89PNG")
    window.comparison.load_vtk()
    assert "OpenGL" in window.message.text()
    assert window.comparison.figure_view.figure is not None


def test_embedded_vtk_scene_uses_existing_geometry(tmp_path):
    from vtkmodules.vtkRenderingCore import vtkRenderer
    from deepstructgenomics.visualization.io_structures import resolve_manifest_bundle
    from deepstructgenomics.visualization.tk_vtk_overlay import OverlayDeltaViewer
    entry = run_rna(tmp_path, {"sequence": "GGGGAAAACCCC", "mutant_sequence": "GGGGAAAACCCA"})
    bundle = resolve_manifest_bundle(entry["path"])
    viewer = OverlayDeltaViewer(bundle["wt_structure"], bundle["mutant_structure"],
                                bundle["wt_score_file"], bundle["mutant_score_file"])
    details = viewer.attach_scene(vtkRenderer())
    assert "Scalar range" in details
    assert not viewer._status_actor.GetVisibility()
    viewer.set_layer_visible("mutant", False)
    assert not viewer.pick_actor.GetVisibility()
    coords = viewer._coords_array_from_structure(viewer._tooltip_payload["structure"])
    assert viewer.position_near(coords[11]) == 12
