#!/usr/bin/env python3
"""Show actual predicted base pairs and score differences from a result bundle."""

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def main():
    parser = argparse.ArgumentParser(description="Comparer les appariements ARN predits, sans geometrie 3D artificielle.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--export", type=Path, help="Enregistrer la figure (par exemple comparaison.png).")
    parser.add_argument("--export-only", action="store_true", help="Exporter sans ouvrir de fenetre.")
    parser.add_argument("--view", choices=["structure", "arcs"], default="structure",
                        help="Tiges/boucles (defaut), ou diagramme d'arcs et deltas.")
    args = parser.parse_args()
    if args.export_only and not args.export:
        parser.error("--export-only requiert --export")

    import matplotlib
    if args.export_only:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from deepstructgenomics.visualization.secondary_view import build_secondary_figure, load_secondary_data
    from deepstructgenomics.visualization.secondary_diagram import build_structure_diagram

    try:
        draw = build_structure_diagram if args.view == "structure" else build_secondary_figure
        figure = draw(load_secondary_data(args.manifest))
        if args.export:
            args.export.parent.mkdir(parents=True, exist_ok=True)
            figure.savefig(args.export, dpi=150, facecolor="white")
            print(f"Figure : {args.export}")
        if not args.export_only:
            plt.show()
        plt.close(figure)
    except (ValueError, OSError) as exc:
        parser.exit(2, f"Erreur : {exc}\n")


if __name__ == "__main__":
    main()
