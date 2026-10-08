"""Topology-aware RNA diagrams and selection shared by desktop and CLI."""
import math

import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.figure import Figure

from .secondary_layout import comparison_layouts, layout_secondary
from deepstructgenomics.alignment.mapping import ComparisonMapping
from .secondary_view import COMMON, DECREASE, GAINED, INCREASE, INK, LOST, MUTED, validate_secondary_data

BASE_COLORS = {"A": "#e4bc52", "U": "#7eb8dc", "G": "#92c5a3", "C": "#bd9fd4"}


class StructureFigure(Figure):
    """Size bases before painting, including the very first PNG export."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.rna_style_updates = []

    def draw(self, renderer):
        if hasattr(self, "rna_compact_grid"):
            # Reserve physical space for legends/titles even in a short desktop pane.
            height = self.get_figheight() * 72
            self.rna_compact_grid.update(top=max(.45, 1 - 85 / height), bottom=min(.35, 62 / height))
            self.legends[0].set_bbox_to_anchor((.5, 1 - 3 / height))
            self.legends[1].set_bbox_to_anchor((.5, 1 - 25 / height))
            self.texts[-2].set_y(40 / height)
            self.texts[-1].set_y(9 / height)
        for axis, update in self.rna_style_updates:
            axis.apply_aspect()
            update(axis)
        super().draw(renderer)


def build_structure_diagram(data, *, compact=False):
    """Draw only the exported pairs; never refold or infer a molecular shape."""
    from matplotlib import pyplot as plt
    validate_secondary_data(data)
    wt, mut = data["reference"], data.get("mutant")
    mapping = ComparisonMapping.from_data(data)
    layouts = comparison_layouts(wt, mut)
    if mapping.alignment and (mapping.alignment.counts["insertion"] or mapping.alignment.counts["deletion"]):
        layouts = (layout_secondary(len(wt["sequence"]), wt["base_pairs"]),
                   layout_secondary(len(mut["sequence"]), mut["base_pairs"]), "Sélection par homologue · dispositions propres à chaque structure")
    wt_pairs = set(map(tuple, wt["base_pairs"]))
    mut_pairs = set(map(tuple, mut["base_pairs"])) if mut else set()
    changes = mapping.classify_pairs(wt_pairs, mut_pairs)
    common = changes["conserved"] if mut else mapping.pairs(wt_pairs, "reference")
    lost, gained = (changes["lost"] | changes["deleted"], changes["gained"] | changes["inserted"]) if mut else (set(), set())
    delta_by_column = mapping.deltas(wt["scores"], mut["scores"]) if mut else {}
    fig = plt.figure(figsize=(12, 8), facecolor="white", FigureClass=StructureFigure)
    fig.rna_coordinates = {}
    fig.rna_selections = {}
    fig.rna_column_positions = {}
    fig.rna_contexts = {}
    fig.rna_layout_note = layouts[2]
    grid = fig.add_gridspec(1, 2, left=.055, right=.97, bottom=.15 if compact else .18,
                           top=.80 if compact else .71, wspace=.17)
    if compact:
        fig.rna_compact_grid = grid
    if not compact:
        fig.text(.055, .945, "Structures secondaires ARN", fontsize=21, color=INK, weight="bold")
        fig.text(.055, .902, str(data.get("identifier", "ARN"))[:65] + " · " + layouts[2], fontsize=10, color=MUTED)
        fig.text(.055, .035, "Schéma des appariements prédits · distances sans signification physique · DeeW69", fontsize=9, color=MUTED)
    base_legend = [Line2D([], [], marker="o", linestyle="", markersize=9,
                         markerfacecolor=color, markeredgecolor="white", label=base)
                   for base, color in BASE_COLORS.items()]
    fig.legend(handles=base_legend, loc="upper center", bbox_to_anchor=(.5, .99 if compact else .875),
               ncol=4, frameon=False, fontsize=10, handletextpad=.3, columnspacing=1.7)
    pair_legend = [Line2D([], [], color=COMMON, linewidth=2, label="Conservée" if mut else "Prédite")]
    if mut:
        pair_legend += [Line2D([], [], color=LOST, linewidth=2, linestyle="--", label=f"Perdue ({len(lost)})"),
                        Line2D([], [], color=GAINED, linewidth=2, label=f"Nouvelle ({len(gained)})")]
    if changes["deleted"] or changes["inserted"]:
        pair_legend[1].set_label(f"Perdue/supprimée ({len(lost)})")
    fig.legend(handles=pair_legend, loc="upper center", bbox_to_anchor=(.5, .93 if compact else .825),
               ncol=3, frameon=False, fontsize=9)
    all_coords = np.vstack([entry.coordinates for entry in layouts[:2] if entry is not None])
    bounds = all_coords.min(axis=0) - 1.4, all_coords.max(axis=0) + 1.4
    for column, (structure, layout, title) in enumerate(zip((wt, mut), layouts[:2], ("Référence · WT", "Mutant"))):
        ax = fig.add_subplot(grid[0, column])
        ax.set_axis_off()
        ax.set_aspect("equal")
        if structure is None:
            ax.text(.5, .5, "Aucun mutant fourni", ha="center", transform=ax.transAxes, color=MUTED)
            continue
        coords = layout.coordinates
        n = len(coords)
        ax.set_xlim(bounds[0][0], bounds[1][0])
        ax.set_ylim(bounds[0][1], bounds[1][1])
        ax.set_title(f"{title} · {len(structure['base_pairs'])} paires", color=INK, fontsize=11, pad=18)
        if n > 1:
            ax.add_collection(LineCollection(np.stack([coords[:-1], coords[1:]], axis=1),
                                             colors="#aebcc7", linewidths=1.3, zorder=1))
        pairs = [tuple(pair) for pair in structure["base_pairs"]]
        side = "reference" if column == 0 else "mutant"
        columns = mapping.to_column[side]
        pair_columns = [(columns[i + 1], columns[j + 1]) for i, j in pairs]
        ax.add_collection(LineCollection([coords[[i, j]] for i, j in pairs],
            colors=[COMMON if pair in common else LOST if column == 0 else GAINED for pair in pair_columns],
            linestyles=["--" if column == 0 and pair in lost else "-" for pair in pair_columns], linewidths=2.2, zorder=2))
        # Labels shrink with dense layouts; the toolbar zoom reveals individual bases.
        span = max(float(np.ptp(coords[:, 0])), float(np.ptp(coords[:, 1])), 5.)
        size = max(12., min(360., 12000. / span ** 2))
        delta = np.array([delta_by_column.get(columns[i + 1], 0.) for i in range(n)])
        changed = [i for i in range(n) if mut and (lambda c: c.reference_base != c.mutant_base)(mapping.columns[columns[i + 1] - 1])]
        affected = np.flatnonzero(np.abs(delta) > 1e-9)
        halo = None
        if len(affected):
            halo = ax.scatter(*coords[affected].T, s=size * 2.3, c=[DECREASE if delta[i] < 0 else INCREASE for i in affected],
                              alpha=.18, edgecolors="none", zorder=2)
        bases = ax.scatter(*coords.T, s=size, c=[BASE_COLORS[base] for base in structure["sequence"]],
                           edgecolors=[INK if i in changed else "white" for i in range(n)], linewidths=1.2, zorder=3)
        bases.set_gid(f"bases-{'wt' if column == 0 else 'mutant'}")
        texts = [ax.text(x, y, base, ha="center", va="center", color=INK, weight="bold", fontsize=9, zorder=4, clip_on=True)
                 for (x, y), base in zip(coords, structure["sequence"])]
        step = max(1, math.ceil(n / 12))
        for i in range(n):
            if i not in (0, n - 1) and (i + 1) % step:
                continue
            previous, following = coords[max(0, i - 1)], coords[min(n - 1, i + 1)]
            direction = following - previous
            normal = np.array([-direction[1], direction[0]])
            normal /= max(np.linalg.norm(normal), 1e-9)
            text = f"{i + 1}" + (" · 5′/3′" if n == 1 else " · 5′" if i == 0 else " · 3′" if i == n - 1 else "")
            offset = (-8 if i == 0 else 8, -21) if i in (0, n - 1) else normal * 14
            alignment = "right" if i == 0 and n > 1 else "left" if i == n - 1 and n > 1 else "center"
            ax.annotate(text, coords[i], xytext=offset, textcoords="offset points",
                        ha=alignment, va="center", fontsize=7, color=MUTED)
        selection = ax.scatter([], [], s=size * 2.6, facecolors="none", edgecolors=INK, linewidths=1.8, zorder=5)
        fig.rna_coordinates[ax] = coords
        fig.rna_column_positions[ax] = [columns[i + 1] for i in range(n)]
        fig.rna_selections[ax] = selection
        fig.rna_contexts[ax] = layout.contexts

        def resize_labels(_axis, axis=ax, labels=texts, dots=bases, ring=selection, aura=halo):
            pixels = np.linalg.norm(axis.transData.transform((1., 0.)) - axis.transData.transform((0., 0.)))
            diameter = min(21., pixels * .62 * 72 / fig.dpi)
            dots.set_sizes([max(2., diameter ** 2)])
            ring.set_sizes([max(7., diameter ** 2 * 2.0)])
            if aura is not None:
                aura.set_sizes([max(5., diameter ** 2 * 2.3)])
            for text in labels:
                text.set_visible(diameter >= 9)
                text.set_fontsize(max(5, min(11, diameter * .48)))
        ax.callbacks.connect("xlim_changed", resize_labels)
        ax.callbacks.connect("ylim_changed", resize_labels)
        fig.rna_style_updates.append((ax, resize_labels))
    if mut:
        caption = f"{len(common)} paire(s) conservée(s) · {len(lost)} perdue(s) · {len(gained)} nouvelle(s)"
    else:
        caption = "Comparaison indisponible sans mutant"
    fig.text(.5, .09 if compact else .12, caption, ha="center", fontsize=10, color=INK)
    fig.text(.5, .035 if compact else .075,
             ("Contour sombre : base différente\nHalo bleu/rouge : diminution/augmentation du score" if compact
              else "Contour sombre : base différente · halo bleu/rouge : diminution/augmentation du score"),
             ha="center", fontsize=8, color=MUTED)
    return fig


def select_position(figure, position):
    for axis, coords in figure.rna_coordinates.items():
        columns = figure.rna_column_positions[axis]
        native = columns.index(position) if position in columns else None
        offsets = coords[native:native + 1] if native is not None else np.empty((0, 2))
        figure.rna_selections[axis].set_offsets(offsets)
    figure.canvas.draw_idle()


def picked_position(figure, event, tolerance=16):
    """Pick in display pixels, not by rounding x as in the old arc diagram."""
    if event.inaxes not in figure.rna_coordinates or event.x is None or event.y is None:
        return None
    coords = event.inaxes.transData.transform(figure.rna_coordinates[event.inaxes])
    distance = np.linalg.norm(coords - [event.x, event.y], axis=1)
    index = int(np.argmin(distance))
    return figure.rna_column_positions[event.inaxes][index] if distance[index] <= tolerance else None
