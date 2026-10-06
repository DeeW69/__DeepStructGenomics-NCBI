#!/usr/bin/env python3
"""CLI wrapper for the DeepStructGenomics pipeline."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import requests

# Ensure the repository root is importable when running from source.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from deepstructgenomics.config import NCBIConfig, PipelineConfig
from deepstructgenomics.pipeline import DeepStructPipeline, PipelineInput


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pipeline sequence -> structure -> impact base sur les donnees NCBI."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--accession", help="Identifiant NCBI (nuccore).")
    source.add_argument("--sequence", help="Sequence ADN/ARN brute (ACGTU).")
    source.add_argument("--fasta", type=Path, help="Fichier FASTA local contenant une seule sequence ADN/ARN.")
    parser.add_argument("--label", help="Nom associe a --sequence ou --fasta (defaut : identifiant FASTA ou custom_sequence).")
    mutant = parser.add_mutually_exclusive_group()
    mutant.add_argument("--mutant-sequence", help="Sequence mutante pour comparaison.")
    mutant.add_argument("--mutant-fasta", type=Path, help="Fichier FASTA local contenant une seule sequence mutante.")
    parser.add_argument("--output-dir", default="outputs", help="Dossier de sortie des rapports.")
    parser.add_argument("--top-k", type=int, default=10, help="Nombre maximal de positions les plus modifiees a inclure dans les rapports (defaut : 10).")
    parser.add_argument(
        "--min-abs-delta",
        "--delta-threshold",
        dest="min_abs_delta",
        type=float,
        default=0.1,
        help="Seuil minimal |delta| pour retenir une position modifiee (defaut : 0.1).",
    )
    parser.add_argument(
        "--base-pair-threshold",
        type=float,
        default=0.2,
        help="Seuil |delta| pour comptabiliser une paire de bases affectee (defaut : 0.2).",
    )
    parser.add_argument("--ncbi-email", help="Email requis par les E-utilities.")
    parser.add_argument("--ncbi-api-key", help="Cle API optionnelle pour augmenter les quotas.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        request = PipelineInput(
            accession=args.accession,
            sequence=args.sequence,
            sequence_label=args.label,
            mutant_sequence=args.mutant_sequence,
            fasta_path=args.fasta,
            mutant_fasta_path=args.mutant_fasta,
            top_k=args.top_k,
            min_abs_delta=args.min_abs_delta,
            base_pair_threshold=args.base_pair_threshold,
        )
        ncbi_config = NCBIConfig(email=args.ncbi_email, api_key=args.ncbi_api_key)
        pipeline = DeepStructPipeline(PipelineConfig(ncbi=ncbi_config))
        result = pipeline.run_and_export(request, output_dir=Path(args.output_dir))
    except requests.RequestException:
        raise SystemExit("Erreur NCBI : requete impossible. Verifier l'accession et la connexion reseau.") from None
    except (ValueError, OSError) as exc:
        print(f"Erreur : {exc}", file=sys.stderr)
        raise SystemExit(2) from None
    print("Rapports generes :")
    print(f"- JSON     : {result.report_paths.json_path}")
    print(f"- Markdown : {result.report_paths.markdown_path}")
    if result.report_paths.csv_path:
        print(f"- CSV      : {result.report_paths.csv_path}")
    if result.visualization_paths and result.visualization_paths.impact_summary_file:
        print(f"- Impact   : {result.visualization_paths.impact_summary_file}")


if __name__ == "__main__":
    main()
