"""Alignment coordinates, conservative structural interpretation and persistence."""
from itertools import product

import pytest

from deepstructgenomics.alignment import AlignmentResult, align_sequences
from deepstructgenomics.alignment.mapping import ComparisonMapping


@pytest.mark.parametrize("wt,mut,counts", [
    ("ACGU", "ACGU", (0, 0, 0)), ("ACGU", "ACGA", (1, 0, 0)), ("ACGU", "UCGC", (2, 0, 0)),
    ("ACGU", "GACGU", (0, 1, 0)), ("ACGU", "ACGGU", (0, 1, 0)), ("ACGU", "ACGUA", (0, 1, 0)),
    ("GACGU", "ACGU", (0, 0, 1)), ("ACGGU", "ACGU", (0, 0, 1)), ("ACGUA", "ACGU", (0, 0, 1)),
    ("ACGU", "ACGGA", (1, 1, 0)), ("ACGGU", "ACGA", (1, 0, 1)),
    ("ACGUACGUACGUACGU", "CGUACGUACGUACGUA", (0, 1, 1)),
])
def test_variant_cases(wt, mut, counts):
    result = align_sequences(wt, mut)
    assert tuple(result.counts[k] for k in ("substitution", "insertion", "deletion")) == counts
    assert result.aligned_reference.replace("-", "") == wt
    assert result.aligned_mutant.replace("-", "") == mut
    assert sum(result.counts.values()) == len(result.columns)
    assert AlignmentResult.from_dict(result.to_dict(), wt, mut) == result


def test_repeated_base_is_deterministic_and_ambiguity_is_visible():
    result = align_sequences("ACGU", "ACGGU")
    assert result.aligned_reference == "ACG-U"
    assert result.columns[3].mutant_position == 4 and result.columns[3].reference_position is None
    assert result.warnings
    assert align_sequences("ACGGU", "ACGU").aligned_mutant == "ACG-U"


def test_mapping_properties_for_short_sequences():
    sequences = ["".join(s) for n in (1, 2, 3) for s in product("AC", repeat=n)]
    for wt, mut in product(sequences, repeat=2):
        result = align_sequences(wt, mut)
        assert list(result.ref_to_alignment) == list(range(1, len(wt) + 1))
        assert list(result.mut_to_alignment) == list(range(1, len(mut) + 1))
        for pos, other in result.ref_to_mut.items():
            if other is not None:
                assert result.mut_to_ref[other] == pos
                assert result.ref_to_alignment[pos] == result.mut_to_alignment[other]
        for c in result.columns:
            assert c.reference_position is not None or c.mutant_position is not None
            assert c.reference_base == (wt[c.reference_position - 1] if c.reference_position else None)
            assert c.mutant_base == (mut[c.mutant_position - 1] if c.mutant_position else None)


def test_pair_comparison_uses_homologous_positions():
    alignment = align_sequences("ACGU", "GACGU")
    mapping = ComparisonMapping("ACGU", "GACGU", alignment)
    pairs = mapping.classify_pairs([(0, 3)], [(1, 4)])
    assert pairs["conserved"] == {(2, 5)}
    assert not any(pairs[k] for k in ("lost", "gained", "deleted", "inserted"))
    changed = mapping.classify_pairs([(0, 3)], [(1, 3)])
    assert changed["lost"] == {(2, 5)} and changed["gained"] == {(2, 4)}
    inserted = mapping.classify_pairs([], [(0, 4)])
    assert inserted["inserted"] == {(1, 5)} and not inserted["gained"]
    mapping = ComparisonMapping("GACGU", "ACGU", align_sequences("GACGU", "ACGU"))
    removed = mapping.classify_pairs([(0, 4)], [])
    assert removed["deleted"] == {(1, 5)} and not removed["lost"]
    assert mapping.deltas([.5] * 5, [.2] * 4) == {i: pytest.approx(-.3) for i in range(2, 6)}


def test_limits_before_aligner_and_corrupted_artifact(monkeypatch):
    monkeypatch.setattr("deepstructgenomics.alignment.pairwise.PairwiseAligner", lambda **kwargs: pytest.fail("limit checked first"))
    with pytest.raises(ValueError, match="coûteux"):
        align_sequences("A" * 2001, "A" * 2001)
    with pytest.raises(ValueError, match="coûteux"):
        align_sequences("A" * 10001, "A")


def test_corrupted_artifact_is_not_silently_recomputed():
    result = align_sequences("ACGU", "ACGGU")
    with pytest.raises(ValueError, match="incohérent"):
        AlignmentResult.from_dict(result.to_dict(), "AAAA", "ACGGU")
    payload = result.to_dict()
    payload["counts"]["match"] = 100
    with pytest.raises(ValueError, match="incohérent"):
        AlignmentResult.from_dict(payload, "ACGU", "ACGGU")


def test_export_reopen_and_legacy_fallback(tmp_path, monkeypatch):
    import json
    from pathlib import Path
    from deepstructgenomics.gui.services import run_rna, read_rna, inspection_details
    entry = run_rna(tmp_path, {"sequence": "ACGU", "mutant_sequence": "ACGGU"})
    import sys
    from scripts import view_3d
    monkeypatch.setattr(sys, "argv", ["view_3d.py", "--manifest", entry["path"], "--mode", "overlay", "--links"])
    with pytest.raises(SystemExit, match="Connecteurs WT/MUT indisponibles"):
        view_3d.main()
    data, metrics = read_rna(entry["path"])
    assert data["alignment"]["counts"]["insertion"] == 1
    detail = inspection_details(data, 4)
    assert detail["reference"] is None and detail["mutant"]["position"] == 4
    assert detail["delta"] is None and detail["changes"][0][1] == "Insertion"
    assert detail["local"]["reference"] == list("ACG-U")
    assert inspection_details(data, 5)["reference"]["position"] == 4
    assert metrics["changed"] <= 4
    manifest = json.loads(Path(entry["path"]).read_text())
    secondary = Path(manifest["secondary_structure_file"])
    payload = json.loads(secondary.read_text())
    del payload["alignment"]
    del payload["alignment_run"]
    secondary.write_text(json.dumps(payload), encoding="utf-8")
    before = secondary.read_bytes()
    old, _ = read_rna(entry["path"])
    assert ComparisonMapping.from_data(old).alignment is None
    assert secondary.read_bytes() == before
    assert inspection_details(old, 4)["reference"]["base"] == "U"
    raw = run_rna(tmp_path, {"sequence": "ACGU", "mutant_sequence": "ACGGU", "align_mutant": False})
    assert read_rna(raw["path"])[0]["alignment"] is None


def test_figures_pick_both_directions_and_gap_has_one_highlight():
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from types import SimpleNamespace
    from deepstructgenomics.visualization.secondary_diagram import build_structure_diagram, picked_position, select_position
    from deepstructgenomics.visualization.secondary_view import build_secondary_figure
    def structure(seq, pairs, dot):
        return {"sequence": seq, "base_pairs": pairs, "dot_bracket": dot, "scores": [.5] * len(seq)}
    data = {"version": 1, "reference": structure("ACGU", [[0, 3]], "(..)"),
            "mutant": structure("GACGU", [[1, 4]], ".(..)"),
            "alignment": align_sequences("ACGU", "GACGU").to_dict()}
    fig = build_structure_diagram(data, compact=True)
    try:
        fig.canvas.draw()
        for axis, native in ((fig.axes[0], 0), (fig.axes[1], 1)):
            x, y = axis.transData.transform(fig.rna_coordinates[axis][native])
            assert picked_position(fig, SimpleNamespace(inaxes=axis, x=x, y=y)) == 2
        select_position(fig, 2)
        assert fig.rna_selections[fig.axes[0]].get_offsets()[0].tolist() == fig.rna_coordinates[fig.axes[0]][0].tolist()
        assert fig.rna_selections[fig.axes[1]].get_offsets()[0].tolist() == fig.rna_coordinates[fig.axes[1]][1].tolist()
        select_position(fig, 1)
        assert len(fig.rna_selections[fig.axes[0]].get_offsets()) == 0
        assert len(fig.rna_selections[fig.axes[1]].get_offsets()) == 1
    finally:
        plt.close(fig)
    arcs = build_secondary_figure(data)
    try:
        assert len(arcs.axes[2].patches) == 4
        assert [bar.get_x() + bar.get_width() / 2 for bar in arcs.axes[2].patches] == [2, 3, 4, 5]
    finally:
        plt.close(arcs)


def test_long_alignment_within_bound_and_no_mutant_guard(tmp_path, monkeypatch):
    assert align_sequences("ACGU" * 100, "ACGU" * 100).identity == 1
    from deepstructgenomics.gui.services import run_rna
    monkeypatch.setattr("deepstructgenomics.pipeline.align_sequences", lambda *args: pytest.fail("no mutant"))
    run_rna(tmp_path, {"sequence": "ACGU"})
