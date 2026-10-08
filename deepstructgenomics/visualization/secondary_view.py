"""Draw exported secondary structures without inventing molecular coordinates."""

from __future__ import annotations

import json
import math
from pathlib import Path

from .io_structures import resolve_manifest_bundle

INK = "#19334b"
MUTED = "#64748b"
COMMON = "#28866c"
LOST = "#287ec0"
GAINED = "#c84c48"
DECREASE = "#287ec0"
INCREASE = "#c84c48"


def load_secondary_data(manifest: str | Path) -> dict:
    bundle = resolve_manifest_bundle(manifest)
    path = bundle.get("secondary_structure_file")
    if not path or not path.is_file():
        raise ValueError("Structures secondaires absentes : regenerer les resultats avec le pipeline actuel.")
    data = json.loads(path.read_text(encoding="utf-8"))
    validate_secondary_data(data)
    return data


def validate_secondary_data(data: dict) -> None:
    if not isinstance(data, dict) or data.get("version") != 1:
        raise ValueError("Format de structures secondaires non pris en charge.")
    for label in ("reference", "mutant"):
        structure = data.get(label)
        if label == "mutant" and structure is None:
            continue
        if not isinstance(structure, dict):
            raise ValueError(f"Structure {label} manquante.")
        sequence = structure.get("sequence")
        if not isinstance(sequence, str) or not sequence or any(base not in "ACGU" for base in sequence):
            raise ValueError(f"Sequence {label} invalide.")
        n = len(sequence)
        dot = structure.get("dot_bracket")
        scores = structure.get("scores")
        pairs = structure.get("base_pairs")
        if not isinstance(dot, str) or len(dot) != n or any(ch not in ".()" for ch in dot):
            raise ValueError(f"Dot-bracket {label} invalide.")
        if not isinstance(scores, list) or len(scores) != n or any(
            not isinstance(v, (int, float)) or not math.isfinite(v) or not 0 <= v <= 1 for v in scores
        ):
            raise ValueError(f"Scores {label} invalides.")
        if not isinstance(pairs, list) or any(
            not isinstance(pair, (list, tuple)) or len(pair) != 2
            or any(type(pos) is not int for pos in pair) or not 0 <= pair[0] < pair[1] < n
            for pair in pairs
        ):
            raise ValueError(f"Appariements {label} invalides.")
        stack, parsed = [], set()
        for pos, char in enumerate(dot):
            if char == "(":
                stack.append(pos)
            elif char == ")":
                if not stack:
                    raise ValueError(f"Dot-bracket {label} desequilibre.")
                parsed.add((stack.pop(), pos))
        if stack or parsed != {tuple(pair) for pair in pairs} or len(parsed) != len(pairs):
            raise ValueError(f"Appariements et dot-bracket {label} incoherents.")


def build_secondary_figure(data: dict, *, compact: bool = False):
    """Return a Matplotlib figure showing exactly the exported pairs and scores."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Arc, Patch

    validate_secondary_data(data)
    wt, mut = data["reference"], data.get("mutant")
    wt_pairs = {tuple(pair) for pair in wt["base_pairs"]}
    mut_pairs = {tuple(pair) for pair in mut["base_pairs"]} if mut else set()
    lost = wt_pairs - mut_pairs if mut else set()
    gained = mut_pairs - wt_pairs if mut else set()
    length = max(len(wt["sequence"]), len(mut["sequence"]) if mut else 0)
    shared = min(len(wt["sequence"]), len(mut["sequence"])) if mut else 0
    changed = {idx for idx in range(length) if mut and
               wt["sequence"][idx:idx + 1] != mut["sequence"][idx:idx + 1]}

    fig = plt.figure(figsize=(13, 8), facecolor="white")
    grid = fig.add_gridspec(2, 2, left=0.075, right=0.97, bottom=0.15,
                           top=0.76, hspace=0.62, wspace=0.18, height_ratios=[1.2, 1])
    fig.text(0.06, 0.94, "ARN : comparer les appariements prédits", fontsize=21, weight="bold", color=INK)
    identifier = str(data.get("identifier", ""))[:70]
    detail = f"{identifier}  ·  WT : {len(wt['sequence'])} bases"
    detail += f"  ·  mutant : {len(mut['sequence'])} bases" if mut else "  ·  référence seule"
    if mut and len(wt["sequence"]) != len(mut["sequence"]):
        detail += "  ·  sans alignement"
    fig.text(0.06, 0.895, detail, fontsize=11, color=MUTED)
    fig.legend(handles=[Patch(color=COMMON, label="Paire conservée" if mut else "Paire prédite"),
                        Patch(color=LOST, label=f"Perdue ({len(lost)})"),
                        Patch(color=GAINED, label=f"Gagnée ({len(gained)})"),
                        Patch(facecolor="#fff2d6", edgecolor=LOST, label="Base différente")],
               loc="upper left", bbox_to_anchor=(0.055, 0.86), ncol=4, frameon=False, fontsize=10)

    for column, (structure, title, special, color) in enumerate((
        (wt, "Référence (WT)", lost, LOST), (mut, "Mutant", gained, GAINED)
    )):
        ax = fig.add_subplot(grid[0, column])
        ax.axis("off")
        if structure is None:
            ax.text(0.5, 0.5, "Aucun mutant fourni", ha="center", transform=ax.transAxes, color=MUTED)
            continue
        sequence = structure["sequence"]
        for i, j in structure["base_pairs"]:
            marked = (i, j) in special
            ax.add_patch(Arc(((i + j) / 2 + 1, 0), width=j - i, height=(j - i) * 0.43,
                             theta1=0, theta2=180, color=color if marked else COMMON,
                             linewidth=2.8 if marked else 1.8,
                             linestyle="--" if marked and column == 0 else "-"))
        positions = list(range(1, len(sequence) + 1))
        ax.plot(positions, [0] * len(sequence), color="#cbd5df", linewidth=1, zorder=1)
        ax.scatter(positions, [0] * len(sequence), s=max(10, min(360, 4400 / length)),
                   c=["#fff2d6" if i in changed else "#e7f1f2" for i in range(len(sequence))],
                   edgecolors=[LOST if i in changed else "#e7f1f2" for i in range(len(sequence))], zorder=3)
        step = max(1, math.ceil(length / 24))
        for i, base in enumerate(sequence):
            if length <= 40:
                ax.text(i + 1, 0, base, ha="center", va="center", fontsize=10, weight="bold", color=INK)
            if i % step == 0 or i == len(sequence) - 1:
                ax.annotate(str(i + 1), (i + 1, 0), xytext=(0, -19), textcoords="offset points",
                            ha="center", color=MUTED, fontsize=9)
        ax.set_xlim(0.3, length + 0.7)
        ax.set_ylim(-max(0.6, length * 0.03), max(2.6, length * 0.24))
        ax.set_title(f"{title}  ·  {len(structure['base_pairs'])} paires", loc="left", color=INK, weight="bold", fontsize=13)

    ax = fig.add_subplot(grid[1, :])
    if not mut:
        ax.axis("off")
        ax.text(0.5, 0.5, "Comparaison indisponible : fournir une séquence mutante.",
                ha="center", transform=ax.transAxes, color=MUTED)
    else:
        positions = list(range(1, shared + 1))
        delta = [mut["scores"][i] - wt["scores"][i] for i in range(shared)]
        ax.bar(positions, delta, width=0.58, color=[DECREASE if v < 0 else INCREASE for v in delta], zorder=3)
        zeros = [p for p, value in zip(positions, delta) if abs(value) < 1e-9]
        ax.scatter(zeros, [0] * len(zeros), s=15, color="#afbfcc", zorder=4)
        if shared <= 40:
            for position, value in zip(positions, delta):
                if abs(value) > 1e-9:
                    ax.text(position, value + (-0.045 if value < 0 else 0.025), f"{value:+.3f}",
                            ha="center", va="top" if value < 0 else "bottom", fontsize=10,
                            weight="bold", color=DECREASE if value < 0 else INCREASE)
        count = sum(abs(value) > 1e-9 for value in delta)
        ax.set_title(f"Δ score = mutant − WT  ·  {count}/{shared} positions communes modifiées",
                     loc="left", color=INK, fontsize=13, pad=12)
        ax.set_xlim(0.4, length + 0.6)
        ax.set_ylim(min(-0.2, min(delta) - 0.18), max(0.2, max(delta) + 0.18))
        ax.set_xticks(list(range(1, length + 1, max(1, math.ceil(length / 24)))))
        ax.set_xlabel("Position dans la séquence (depuis 1)", color=INK)
        ax.set_ylabel("Δ score", color=INK)
        ax.axhline(0, color=MUTED, linewidth=0.8)
        ax.grid(axis="y", color="#e8edf2", zorder=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.text(0.06, 0.055, "Schéma de structure secondaire · scores heuristiques · aucune conformation 3D prédite",
             color=MUTED, fontsize=10)
    fig.text(0.97, 0.055, "DeeW69", ha="right", color=MUTED, fontsize=10)
    if compact:
        # The desktop shell supplies the title, identifier and scientific note.
        for text in list(fig.texts):
            text.remove()
        fig.legends[0].remove()
        fig.legend(handles=[Patch(color=COMMON, label="Conservée" if mut else "Prédite"),
                            Patch(color=LOST, label=f"Perdue ({len(lost)})"),
                            Patch(color=GAINED, label=f"Gagnée ({len(gained)})"),
                            Patch(facecolor="#fff2d6", edgecolor=LOST, label="Base différente")],
                   loc="upper center", bbox_to_anchor=(0.5, .99), ncol=2, frameon=False, fontsize=9)
        grid.update(top=.82, bottom=.12, left=.10, right=.97, hspace=.72)
        for axis in fig.axes:
            axis.set_title(axis.get_title(), fontsize=10)
            axis.set_title(axis.get_title(loc="left"), loc="left", fontsize=10)
        if mut:
            ax.set_title(f"Δ score (mutant − WT) · {count}/{shared} positions modifiées",
                         loc="left", fontsize=10)
    return fig
