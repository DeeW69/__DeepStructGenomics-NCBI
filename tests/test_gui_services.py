"""Desktop operations must preserve pipeline semantics without requiring Qt."""
import json
from pathlib import Path
import subprocess
import sys

import pytest

from deepstructgenomics.gui.services import RecentStore, position_details, read_rna, run_rna


def test_desktop_rna_results_and_position_conventions(tmp_path):
    entry = run_rna(tmp_path, {"sequence": "GGGGAAAACCCC", "mutant_sequence": "GGGGAAAACCCA"})
    data, metrics = read_rna(entry["path"])
    assert metrics["wt_pairs"] == 4 and metrics["mut_pairs"] == 3
    assert metrics["lost"] == 1 and metrics["gained"] == 0
    assert metrics["changed"] == 2
    assert metrics["mean_delta"] == pytest.approx(-1.3 / 12)
    assert metrics["max_abs_delta"] == pytest.approx(.7)
    detail = position_details(data, 12)
    assert detail["reference"]["base"] == "C"
    assert detail["reference"]["partner"] == 1
    assert detail["mutant"]["base"] == "A"
    assert detail["mutant"]["partner"] is None
    assert detail["delta"] == pytest.approx(-.7)
    second = run_rna(tmp_path, {"sequence": "GGGGAAAACCCC"})
    assert second["directory"] != entry["directory"]
    _, metrics = read_rna(second["path"])
    assert metrics["max_abs_delta"] is None and metrics["lost"] is None


def test_longer_mutant_does_not_invent_reference_base(tmp_path):
    entry = run_rna(tmp_path, {"sequence": "GGGGAAAACCCC", "mutant_sequence": "GGGGAAAACCCCA"})
    data, _ = read_rna(entry["path"])
    detail = position_details(data, 13)
    assert detail["reference"] is None and detail["delta"] is None
    assert detail["mutant"]["base"] == "A"


def test_recent_store_recovers_and_deduplicates(tmp_path):
    store = RecentStore(tmp_path / "recent.json")
    store.path.write_text("{broken", encoding="utf-8")
    assert store.load() == [] and store.warning
    entry = {"kind": "rna", "path": "unavailable.json", "directory": "results", "label": "ARN"}
    store.add(entry)
    store.add(entry)
    assert store.load() == [entry]


def test_worker_success_and_validation_error(tmp_path):
    request = {"kind": "rna", "workspace": str(tmp_path), "parameters": {"sequence": "GGGGAAAACCCC"}}
    command = [sys.executable, "-m", "deepstructgenomics.gui.worker"]
    process = subprocess.run(command, input=json.dumps(request), text=True, capture_output=True, timeout=30)
    assert process.returncode == 0, process.stderr
    assert Path(json.loads(process.stdout)["path"]).is_file()
    request["parameters"]["sequence"] = "ACG!"
    process = subprocess.run(command, input=json.dumps(request), text=True, capture_output=True, timeout=30)
    assert process.returncode == 1
    assert "error" in json.loads(process.stdout)


def test_fasta_and_ncbi_sources_delegate_to_pipeline(tmp_path, monkeypatch):
    from deepstructgenomics.data_sources.ncbi_client import NCBIClient, normalize_user_sequence
    fasta = tmp_path / "reference.fasta"
    fasta.write_text(">local_record\nGGGGAAAACCCC\n", encoding="utf-8")
    entry = run_rna(tmp_path, {"fasta_path": str(fasta)})
    assert entry["label"] == "local_record"
    calls = []

    def fetch(self, accession, *, refresh=False):
        calls.append((accession, refresh))
        return normalize_user_sequence("GGGGAAAACCCC", accession)

    monkeypatch.setattr(NCBIClient, "fetch_sequence", fetch)
    entry = run_rna(tmp_path, {"accession": "TEST.1", "refresh_cache": True, "ncbi_email": "example@example.org"})
    assert entry["label"] == "TEST.1"
    assert calls == [("TEST.1", True)]
