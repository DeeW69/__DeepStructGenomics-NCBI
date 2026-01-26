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
from deepstructgenomics.visualization.tk_vtk_overlay import OverlayViewer


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Visualisation 3D Tkinter + VTK (DeepStructGenomics).")
    parser.add_argument(
        "--mode",
        choices=["minimal", "molecule", "overlay"],
        required=True,
        help="Type de viewer a afficher.",
    )
    parser.add_argument("--structure", help="Chemin vers une structure (PDB/mmCIF) pour le mode molecule.")
    parser.add_argument("--wt", help="Structure de reference pour le mode overlay.")
    parser.add_argument("--mut", help="Structure mutante pour le mode overlay.")
    parser.add_argument("--scores", help="Fichier JSON de scores (mode overlay).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.mode == "minimal":
        viewer = MinimalVTKViewer()
    elif args.mode == "molecule":
        if not args.structure:
            raise SystemExit("--structure est requis pour le mode molecule.")
        viewer = MoleculeViewer(args.structure)
    else:
        if not (args.wt and args.mut and args.scores):
            raise SystemExit("--wt, --mut et --scores sont requis pour le mode overlay.")
        viewer = OverlayViewer(args.wt, args.mut, args.scores)
    viewer.start()


if __name__ == "__main__":
    main()
