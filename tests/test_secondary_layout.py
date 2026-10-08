"""Topology, honest comparison and picking for the stem/loop RNA renderer."""
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
import pytest

from deepstructgenomics.visualization.secondary_layout import comparison_layouts, layout_secondary
from deepstructgenomics.visualization.secondary_diagram import build_structure_diagram, picked_position, select_position


def structure(dot):
    stack, pairs = [], []
    for index, char in enumerate(dot):
        if char == "(":
            stack.append(index)
        elif char == ")":
            pairs.append([stack.pop(), index])
    assert not stack
    return {"sequence": ("ACGU" * len(dot))[:len(dot)], "dot_bracket": dot,
            "base_pairs": pairs, "scores": [.5] * len(dot)}


def test_hairpin_ladder_loop_and_every_backbone_edge():
    payload = structure("((((....))))")
    layout = layout_secondary(12, payload["base_pairs"])
    coords = layout.coordinates
    assert np.allclose(coords[:4, 0], coords[0, 0])
    assert np.allclose(coords[8:, 0], coords[8, 0])
    assert coords[5, 1] > coords[4, 1] > coords[3, 1]
    assert np.linalg.norm(np.diff(coords, axis=0), axis=1) == pytest.approx([1.4] * 11)
    assert layout.contexts[4:8] == ["Boucle terminale"] * 4
    assert layout.contexts[0] == "Tige appariée"


@pytest.mark.parametrize("dot, expected", [
    ("....", "Région non appariée"),
    ("((..((...))..))", "Boucle interne"),
    ("((..((...))))", "Renflement"),
    ("((..((...))..((...))..))", "Jonction multibranche"),
    ("...((...))..((...))...", "Boucle terminale"),
    (".", "Région non appariée"),
])
def test_loop_contexts_and_finite_distinct_positions(dot, expected):
    payload = structure(dot)
    layout = layout_secondary(len(dot), payload["base_pairs"])
    assert layout.coordinates.shape == (len(dot), 2)
    assert np.isfinite(layout.coordinates).all()
    assert len(np.unique(np.round(layout.coordinates, 7), axis=0)) == len(dot)
    assert expected in layout.contexts


def test_terminal_pair_loss_keeps_reference_frame_but_not_false_context():
    wt, mut = structure("((((....))))"), structure(".(((....))).")
    left, right, note = comparison_layouts(wt, mut)
    assert np.array_equal(left.coordinates, right.coordinates)
    assert "fixées" in note
    assert left.contexts[0] == "Tige appariée"
    assert right.contexts[0] == "Région non appariée"
    assert wt["base_pairs"] != mut["base_pairs"]


def test_incompatible_pairs_get_independent_layouts_and_no_union_edges():
    wt, mut = structure("((....))...."), structure("....((....))")
    left, right, note = comparison_layouts(wt, mut)
    assert "propres" in note
    assert not np.array_equal(left.coordinates, right.coordinates)
    data = {"version": 1, "reference": wt, "mutant": mut}
    fig = build_structure_diagram(data)
    try:
        for axis, payload in zip(fig.axes, (wt, mut)):
            assert len(axis.collections[1].get_segments()) == len(payload["base_pairs"])
    finally:
        plt.close(fig)


def test_deep_stem_no_recursion_and_crossings_rejected():
    pairs = [(i, 2402 - i) for i in range(1200)]
    layout = layout_secondary(2403, pairs)
    assert np.isfinite(layout.coordinates).all()
    with pytest.raises(ValueError, match="pseudonœuds"):
        layout_secondary(8, [(0, 5), (2, 7)])
    with pytest.raises(ValueError, match="incompatibles"):
        layout_secondary(8, [(0, 7), (0, 6)])


def test_picking_uses_both_screen_coordinates_and_selection_is_synchronized():
    fig = build_structure_diagram({"version": 1, "reference": structure("((((....))))"),
                                  "mutant": structure(".(((....))).")}, compact=True)
    try:
        fig.canvas.draw()
        axis = fig.axes[0]
        coords = fig.rna_coordinates[axis]
        # These bases share x, so the old x-rounding picker cannot work.
        for index in (0, 3, 11):
            x, y = axis.transData.transform(coords[index])
            assert picked_position(fig, SimpleNamespace(inaxes=axis, x=x, y=y)) == index + 1
        select_position(fig, 12)
        for ax in fig.rna_coordinates:
            assert np.allclose(fig.rna_selections[ax].get_offsets()[0], fig.rna_coordinates[ax][11])
        assert picked_position(fig, SimpleNamespace(inaxes=axis, x=-100, y=-100)) is None
    finally:
        plt.close(fig)


def test_longer_mutant_missing_reference_position_and_reference_only():
    wt, mut = structure("((....))"), structure("((....))...")
    left, right, note = comparison_layouts(wt, mut)
    assert len(left.coordinates) == 8 and len(right.coordinates) == 11
    assert "sans alignement" in note
    fig = build_structure_diagram({"version": 1, "reference": wt, "mutant": mut})
    try:
        select_position(fig, 11)
        assert fig.rna_selections[fig.axes[0]].get_offsets().shape == (0, 2)
    finally:
        plt.close(fig)
    fig = build_structure_diagram({"version": 1, "reference": wt, "mutant": None})
    try:
        assert len(fig.rna_coordinates) == 1
        assert any("sans mutant" in text.get_text() for text in fig.texts)
    finally:
        plt.close(fig)


def test_dense_labels_are_sized_before_first_draw_and_reappear_on_zoom():
    payload = structure("(" * 60 + "...." + ")" * 60)
    fig = build_structure_diagram({"version": 1, "reference": payload, "mutant": None})
    try:
        fig.canvas.draw()
        axis = fig.axes[0]
        bases = [text for text in axis.texts if text.get_text() in ("A", "C", "G", "U")]
        assert not any(text.get_visible() for text in bases)
        x, y = fig.rna_coordinates[axis][0]
        axis.set_xlim(x - 2, x + 3)
        axis.set_ylim(y - 1, y + 5)
        fig.canvas.draw()
        assert all(text.get_visible() and text.get_clip_on() for text in bases)
    finally:
        plt.close(fig)
