import json

import pytest

from deepstructgenomics.rna.evaluation import dot_bracket_pairs, evaluate_reference_set, pair_metrics


@pytest.mark.parametrize("structure", ["(()", "())", "[..]", ")..("])
def test_invalid_structure(structure):
    with pytest.raises(ValueError):
        dot_bracket_pairs(structure, len(structure))


def test_exact_metrics_and_undefined():
    metrics = pair_metrics({(0, 8), (1, 7)}, {(0, 8), (2, 6)})
    assert metrics == {"tp": 1, "fp": 1, "fn": 1, "precision": .5, "recall": .5, "f1": .5}
    assert pair_metrics(set(), set())["f1"] is None


def test_external_method_and_micro_counts(tmp_path):
    path = tmp_path / "references.json"
    record = {"id": "hairpin", "sequence": "GGGGAAAACCCC", "reference": "((((....))))",
              "predictions": {"external-v1": "............"}}
    path.write_text(json.dumps({"source": "synthetic", "records": [record]}))
    result = evaluate_reference_set(path)
    assert result["micro_average"]["nussinov"]["f1"] == 1
    assert result["micro_average"]["external-v1"]["fn"] == 4
    assert len(result["dataset_sha256"]) == 64
    path.write_text(json.dumps({"source": "synthetic", "records": [record, record]}))
    with pytest.raises(ValueError, match="duplique"):
        evaluate_reference_set(path)


def test_reject_partial_method_coverage(tmp_path):
    path = tmp_path / "references.json"
    path.write_text(json.dumps({"source": "synthetic", "records": [
        {"id": "a", "sequence": "AAAA", "reference": "....", "predictions": {"method": "...."}},
        {"id": "b", "sequence": "AAAA", "reference": "...."}]}))
    with pytest.raises(ValueError, match="couvrir"):
        evaluate_reference_set(path)
