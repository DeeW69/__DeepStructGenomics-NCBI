#!/usr/bin/env python3
"""Rebuild the README figure from the bundled example and pipeline results."""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from deepstructgenomics.config import PipelineConfig
from deepstructgenomics.pipeline import DeepStructPipeline, PipelineInput
from deepstructgenomics.rna.secondary_structure import SecondaryStructureResult
from deepstructgenomics.visualization.score_mapping import (
    ScoreTable,
    compute_delta_scores,
    derive_position_scores_from_structure,
)

INK = "#19334b"
TEAL = "#187f86"
ORANGE = "#c95636"
MUTED = "#64748b"


def draw_structure(ax, structure: SecondaryStructureResult, title: str, lost_pairs: set) -> None:
    """Draw only base pairs actually predicted for this structure."""
    sequence = structure.sequence
    for i, j in structure.base_pairs:
        lost = (i, j) in lost_pairs
        ax.add_patch(
            Arc(
                ((i + j) / 2 + 1, 0),
                width=j - i,
                height=(j - i) * 0.43,
                theta1=0,
                theta2=180,
                color=ORANGE if lost else TEAL,
                linewidth=2.5 if lost else 1.8,
            )
        )
    positions = list(range(1, len(sequence) + 1))
    ax.scatter(positions, [0] * len(sequence), s=375, color="#e7f1f2", zorder=3)
    for position, base in enumerate(sequence, start=1):
        ax.text(position, 0, base, ha="center", va="center", color=INK, fontsize=11, weight="bold")
        ax.text(position, -0.32, str(position), ha="center", va="center", color=MUTED, fontsize=9)
    ax.set_title(f"{title}  ·  {len(structure.base_pairs)} paires", loc="left", fontsize=13, weight="bold", color=INK, pad=14)
    ax.text(0.5, -0.03, structure.dot_bracket, transform=ax.transAxes, ha="center", va="top", family="monospace", fontsize=12, color=MUTED)
    ax.set_xlim(0.4, len(sequence) + 0.6)
    ax.set_ylim(-0.6, max(2.6, len(sequence) * 0.23))
    ax.axis("off")


def main() -> None:
    parser = argparse.ArgumentParser(description="Regenerer la figure de la demo, sans reseau ni fenetre.")
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "docs" / "assets" / "demo_delta.png")
    args = parser.parse_args()
    example = json.loads((PROJECT_ROOT / "data" / "examples" / "rna_variant.json").read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="deepstruct-demo-") as temporary_dir:
        work_dir = Path(temporary_dir)
        pipeline = DeepStructPipeline(PipelineConfig(default_output_dir=work_dir, cache_dir=work_dir / ".cache"))
        result = pipeline.run(
            PipelineInput(
                sequence=example["reference_sequence"],
                sequence_label=example["identifier"],
                mutant_sequence=example["mutant_sequence"],
            )
        )

    wt = result.structure
    mutant = result.mutant_structure
    delta = compute_delta_scores(
        ScoreTable(position_scores=derive_position_scores_from_structure(wt)),
        ScoreTable(position_scores=derive_position_scores_from_structure(mutant)),
    ).position_scores
    lost_pairs = set(wt.base_pairs) - set(mutant.base_pairs)
    substitutions = ", ".join(
        f"{sub['reference']}{sub['position'] + 1}{sub['mutant']}"
        for sub in result.variant_result.substitutions
    )

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    fig = plt.figure(figsize=(11.5, 7.0), facecolor="white")
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.1], left=0.085, right=0.965, bottom=0.15, top=0.79, hspace=0.66, wspace=0.18)
    fig.text(0.06, 0.94, "Une substitution, deux positions affectées", fontsize=21, weight="bold", color=INK)
    fig.text(0.06, 0.89, f"ARN synthétique · {len(wt.sequence)} nucléotides · {substitutions} · comparaison référence / mutant", fontsize=12, color=MUTED)
    draw_structure(fig.add_subplot(grid[0, 0]), wt, "Référence (WT)", lost_pairs)
    draw_structure(fig.add_subplot(grid[0, 1]), mutant, "Mutant", set())
    fig.text(0.06, 0.85, f"En orange : paires WT perdues dans le mutant ({len(lost_pairs)})", fontsize=10, color=ORANGE)

    ax = fig.add_subplot(grid[1, :])
    positions = list(delta)
    values = list(delta.values())
    colors = [ORANGE if value < 0 else TEAL for value in values]
    ax.bar(positions, values, width=0.58, color=colors, zorder=3)
    zeros = [position for position, value in delta.items() if abs(value) < 1e-9]
    ax.scatter(zeros, [0] * len(zeros), color="#b5c4cf", s=24, zorder=4)
    for position, value in delta.items():
        if abs(value) > 1e-9:
            ax.text(position, value - 0.035 if value < 0 else value + 0.035, f"{value:+.3f}".replace(".", ","), ha="center", va="top" if value < 0 else "bottom", weight="bold", color=ORANGE if value < 0 else TEAL, fontsize=11)
    changed = sum(abs(value) > 1e-9 for value in values)
    ax.set_title(f"Δ score = mutant − WT  ·  {changed}/{len(values)} positions modifiées", loc="left", color=INK, fontsize=13, weight="bold", pad=12)
    ax.axhline(0, color=MUTED, linewidth=0.8)
    ax.set_xlim(0.4, len(wt.sequence) + 0.6)
    ax.set_ylim(min(-0.2, min(values) - 0.16), max(0.15, max(values) + 0.15))
    ax.set_xticks(positions)
    ax.set_xlabel("Position dans la séquence (depuis 1)", color=INK, labelpad=8)
    ax.set_ylabel("Δ score", color=INK)
    ax.tick_params(colors=MUTED, length=0, pad=7)
    ax.grid(axis="y", color="#e8edf2", linewidth=0.8, zorder=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.text(0.06, 0.055, "Scores heuristiques d'appariement · aucune interprétation clinique", color=MUTED, fontsize=10)
    fig.text(0.965, 0.055, "DeeW69", ha="right", color=MUTED, fontsize=10)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=150, facecolor=fig.get_facecolor(), metadata={"Author": "DeeW69"})
    plt.close(fig)
    print(f"Figure : {args.output}")
    print(f"Paires WT/mutant : {len(wt.base_pairs)}/{len(mutant.base_pairs)} ; deltas : {delta}")


if __name__ == "__main__":
    main()
