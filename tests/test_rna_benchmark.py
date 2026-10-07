import importlib.util
import json
from pathlib import Path
import zipfile

import pytest

from deepstructgenomics.rna.evaluation import dot_bracket_pairs, evaluate_reference_set
from deepstructgenomics.rna.vienna import ViennaPredictor

ROOT = Path(__file__).resolve().parents[1]


def test_bundled_experimental_references_are_traceable():
    dataset = json.loads((ROOT / "data/benchmarks/bprna_pdb.json").read_text())
    assert len(dataset["records"]) == 30
    assert len({r["sequence"] for r in dataset["records"]}) == 30
    assert dataset["selection"]["pdb_records"] == 669
    for record in dataset["records"]:
        assert record["id"].startswith("bpRNA_PDB_")
        assert len(record["provenance"]["sha256"]) == 64
        dot_bracket_pairs(record["reference"], len(record["sequence"]))


def test_preparation_selects_before_prediction_and_reproduces(tmp_path):
    spec = importlib.util.spec_from_file_location("prepare", ROOT / "scripts/prepare_rna_benchmark.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    archive = tmp_path / "records.zip"
    with zipfile.ZipFile(archive, "w") as z:
        z.writestr("stFiles/bpRNA_PDB_1.st", "#Name: test\n" + "A" * 20 + "\n" + "." * 20)
        z.writestr("stFiles/bpRNA_PDB_2.st", "#Name: test\n" + "A" * 20 + "\n" + "." * 20)
    a = module.prepare(archive, tmp_path / "a.json", size=1)
    b = module.prepare(archive, tmp_path / "b.json", size=1)
    assert a == b
    assert a["selection"]["excluded"]["duplicate_sequence"] == 1


def test_missing_vienna_has_actionable_error(monkeypatch):
    def missing(_):
        raise ImportError("missing")
    monkeypatch.setattr("deepstructgenomics.rna.vienna.import_module", missing)
    with pytest.raises(ValueError, match="requirements-research"):
        ViennaPredictor()


def test_live_vienna_agrees_with_archived_benchmark():
    pytest.importorskip("RNA")
    report = evaluate_reference_set(ROOT / "data/benchmarks/bprna_pdb.json", viennarna=True)
    archived = json.loads((ROOT / "docs/benchmarks/v0.4.5_rna.json").read_text())
    assert report["method_provenance"]["viennarna_mfe"]["executed"]
    assert report["micro_average"] == archived["micro_average"]
    assert [r["predictions"] for r in report["records"]] == [r["predictions"] for r in archived["records"]]
