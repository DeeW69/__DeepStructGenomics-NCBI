#!/usr/bin/env python3
"""Tkinter + VTK viewer launcher."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from deepstructgenomics.visualization.tk_vtk_minimal import MinimalVTKViewer
from deepstructgenomics.visualization.tk_vtk_molecule import MoleculeViewer
from deepstructgenomics.visualization.tk_vtk_overlay import OverlayDeltaViewer, OverlayViewer


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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.export_only and not args.export:
        raise SystemExit("--export-only requiert --export.")
    export_path = Path(args.export) if args.export else None

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
        )
    viewer.start(export_path=export_path, export_only=args.export_only)


if __name__ == "__main__":
    main()
