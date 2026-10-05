#!/usr/bin/env python3
"""Tkinter + VTK viewer launcher."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from deepstructgenomics.visualization.export_helpers import parse_export_view
from deepstructgenomics.visualization.overlay_helpers import load_base_pairs_json
from deepstructgenomics.visualization.tk_vtk_minimal import MinimalVTKViewer
from deepstructgenomics.visualization.tk_vtk_molecule import MoleculeViewer
from deepstructgenomics.visualization.tk_vtk_overlay import BasePairOptions, OverlayDeltaViewer, OverlayViewer


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Visualisation 3D Tkinter + VTK (DeepStructGenomics).")
    parser.add_argument(
        "--mode",
        choices=["minimal", "molecule", "overlay", "overlay-delta"],
        required=True,
        help="Type de viewer a afficher.",
    )
    parser.add_argument("--structure", help="Chemin vers une structure (PDB/mmCIF) pour le mode molecule.")
    parser.add_argument("--wt", help="Structure de reference pour le mode overlay.")
    parser.add_argument("--mut", help="Structure mutante pour le mode overlay.")
    parser.add_argument("--scores", help="Fichier JSON de scores (mode overlay).")
    parser.add_argument("--wt-scores", help="Fichier JSON des scores WT (mode overlay-delta).")
    parser.add_argument("--mut-scores", help="Fichier JSON des scores mutant (mode overlay-delta).")
    parser.add_argument(
        "--backbone",
        action="store_true",
        help="Affiche le backbone (molecule).",
    )
    parser.add_argument(
        "--no-backbone",
        action="store_true",
        help="Desactive le backbone pour overlay/overlay-delta.",
    )
    parser.add_argument("--export", help="Chemin PNG a generer avant ouverture de la fenetre.")
    parser.add_argument(
        "--export-only",
        action="store_true",
        help="Capture le PNG puis quitte sans lancer l'interacteur.",
    )
    parser.add_argument(
        "--tooltips",
        action="store_true",
        help="Active les tooltips (mode compact) meme si DSG_VIZ_DEBUG!=1.",
    )
    parser.add_argument(
        "--tooltips-verbose",
        action="store_true",
        help="Active les tooltips detailles (sans logs) meme si DSG_VIZ_DEBUG!=1.",
    )
    parser.add_argument(
        "--export-scale",
        type=int,
        default=1,
        help="Facteur de sur-echantillonnage pour --export (1 = taille fenetre).",
    )
    parser.add_argument(
        "--export-view",
        default="auto",
        help="Preset camera pour les exports (auto, iso, top, side, front).",
    )
    parser.add_argument(
        "--export-hide-ui",
        action="store_true",
        help="Cache les elements UI (legende/barre) uniquement pendant l'export.",
    )
    parser.add_argument(
        "--export-background",
        choices=["dark", "white", "transparent"],
        default="dark",
        help="Fond utilise pour le PNG exporte (dark = valeur par defaut interactive).",
    )
    parser.add_argument(
        "--links",
        action="store_true",
        help="Active les liens WT->MUT (segments positionnels).",
    )
    parser.add_argument(
        "--no-links",
        action="store_true",
        help="Force la desactivation des liens meme si --links est present.",
    )
    parser.add_argument(
        "--link-threshold",
        type=float,
        default=None,
        help="Filtre les liens selon |valeur| >= seuil (score ou delta).",
    )
    parser.add_argument(
        "--link-metric",
        choices=["score", "delta", "distance"],
        help="Choisit la metrique des liens (score, delta ou distance).",
    )
    parser.add_argument(
        "--link-min-distance",
        type=float,
        default=0.0,
        help="Distance minimale pour afficher un lien (en unite VTK).",
    )
    parser.add_argument(
        "--link-max-distance",
        type=float,
        default=None,
        help="Distance maximale pour afficher un lien (en unite VTK).",
    )
    parser.add_argument("--base-pairs-wt", help="JSON des paires de bases WT (1-based).")
    parser.add_argument("--base-pairs-mut", help="JSON des paires de bases mutant (1-based).")
    parser.add_argument(
        "--base-pairs-mode",
        choices=["lines", "arcs"],
        default="lines",
        help="Aspect des connecteurs de paires (segments ou arcs).",
    )
    parser.add_argument(
        "--base-pairs-threshold",
        type=float,
        default=0.0,
        help="Seuil delta pour filtrer les paires (overlay-delta uniquement).",
    )
    parser.add_argument(
        "--no-base-pairs",
        action="store_true",
        help="Desactive les connecteurs de paires meme si des JSON sont fournis.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.export_only and not args.export:
        raise SystemExit("--export-only requiert --export.")
    if args.export_scale < 1:
        raise SystemExit("--export-scale doit etre >= 1.")
    export_path = Path(args.export) if args.export else None
    try:
        export_view = parse_export_view(args.export_view)
    except ValueError as exc:  # pragma: no cover - validation
        raise SystemExit(str(exc))
    export_hide_ui = args.export_hide_ui
    export_background = args.export_background
    if args.tooltips_verbose:
        tooltip_enabled: Optional[bool] = True
        tooltip_verbose: Optional[bool] = True
    elif args.tooltips:
        tooltip_enabled = True
        tooltip_verbose = False
    else:
        tooltip_enabled = None
        tooltip_verbose = None
    links_enabled = args.links and not args.no_links
    link_threshold = args.link_threshold
    link_metric = args.link_metric
    if not link_metric:
        link_metric = "delta" if args.mode == "overlay-delta" else "score"
    if link_metric == "delta" and args.mode != "overlay-delta":
        raise SystemExit("--link-metric delta requiert --mode overlay-delta.")
    link_min_distance = max(0.0, float(args.link_min_distance or 0.0))
    if args.link_max_distance is not None:
        if args.link_max_distance < link_min_distance:
            raise SystemExit("--link-max-distance doit etre >= --link-min-distance.")
        link_max_distance: Optional[float] = float(args.link_max_distance)
    else:
        link_max_distance = None
    base_pair_config: Optional[BasePairOptions]
    if args.no_base_pairs:
        base_pair_config = None
    else:
        wt_pairs: List[Tuple[int, int]] = []
        mut_pairs: List[Tuple[int, int]] = []
        if args.base_pairs_wt:
            wt_pairs = load_base_pairs_json(args.base_pairs_wt)
        if args.base_pairs_mut:
            mut_pairs = load_base_pairs_json(args.base_pairs_mut)
        if wt_pairs or mut_pairs:
            base_pair_config = BasePairOptions(
                enabled=True,
                mode=args.base_pairs_mode,
                threshold=max(0.0, float(args.base_pairs_threshold or 0.0)),
                wt_pairs=tuple(wt_pairs),
                mut_pairs=tuple(mut_pairs),
            )
        else:
            base_pair_config = None

    if args.mode == "minimal":
        viewer = MinimalVTKViewer()
    elif args.mode == "molecule":
        if not args.structure:
            raise SystemExit("--structure est requis pour le mode molecule.")
        viewer = MoleculeViewer(args.structure, show_backbone=args.backbone)
    elif args.mode == "overlay":
        if not (args.wt and args.mut and args.scores):
            raise SystemExit("--wt, --mut et --scores sont requis pour le mode overlay.")
        viewer = OverlayViewer(
            args.wt,
            args.mut,
            args.scores,
            show_backbone=not args.no_backbone,
            enable_tooltip=tooltip_enabled,
            tooltip_verbose=tooltip_verbose,
            enable_links=links_enabled,
            link_threshold=link_threshold,
            link_metric=link_metric,
            link_min_distance=link_min_distance,
            link_max_distance=link_max_distance,
            base_pair_config=base_pair_config,
        )
    else:  # overlay-delta
        if not (args.wt and args.mut and args.wt_scores and args.mut_scores):
            raise SystemExit("--wt, --mut, --wt-scores et --mut-scores sont requis pour overlay-delta.")
        viewer = OverlayDeltaViewer(
            args.wt,
            args.mut,
            args.wt_scores,
            args.mut_scores,
            show_backbone=not args.no_backbone,
            enable_tooltip=tooltip_enabled,
            tooltip_verbose=tooltip_verbose,
            enable_links=links_enabled,
            link_threshold=link_threshold,
            link_metric=link_metric,
            link_min_distance=link_min_distance,
            link_max_distance=link_max_distance,
            base_pair_config=base_pair_config,
        )
    viewer.start(
        export_path=export_path,
        export_only=args.export_only,
        export_scale=args.export_scale,
        export_view=export_view,
        export_hide_ui=export_hide_ui,
        export_background=export_background,
    )


if __name__ == "__main__":
    main()
