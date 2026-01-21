#!/usr/bin/env python3
"""CLI wrapper for the DeepStructGenomics pipeline."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

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
    source.add_argument("--sequence", help="Sequence ARN brute (AUGC).")
    parser.add_argument("--label", default="custom_sequence", help="Nom associe a la sequence fournie.")
    parser.add_argument("--mutant-sequence", help="Sequence mutante pour comparaison.")
    parser.add_argument("--output-dir", default="outputs", help="Dossier de sortie des rapports.")
    parser.add_argument("--ncbi-email", help="Email requis par les E-utilities.")
    parser.add_argument("--ncbi-api-key", help="Cle API optionnelle pour augmenter les quotas.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ncbi_config = NCBIConfig(email=args.ncbi_email, api_key=args.ncbi_api_key)
    pipeline = DeepStructPipeline(PipelineConfig(ncbi=ncbi_config))
    request = PipelineInput(
        accession=args.accession,
        sequence=args.sequence,
        sequence_label=args.label,
        mutant_sequence=args.mutant_sequence,
    )
    result = pipeline.run_and_export(request, output_dir=Path(args.output_dir))
    print("Rapports generes :")
    print(f"- JSON     : {result.report_paths.json_path}")
    print(f"- Markdown : {result.report_paths.markdown_path}")


if __name__ == "__main__":
    main()
