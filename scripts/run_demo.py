#!/usr/bin/env python3
"""Run the bundled synthetic RNA comparison without network or GUI calls."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from deepstructgenomics.config import PipelineConfig
from deepstructgenomics.pipeline import DeepStructPipeline, PipelineInput


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Demo hors ligne : comparaison de deux petits ARN synthetiques."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs_demo"),
        help="Dossier de sortie (defaut : outputs_demo).",
    )
    args = parser.parse_args()
    example_path = PROJECT_ROOT / "data" / "examples" / "rna_variant.json"
    example = json.loads(example_path.read_text(encoding="utf-8"))

    pipeline = DeepStructPipeline(
        PipelineConfig(
            default_output_dir=args.output_dir,
            cache_dir=args.output_dir / ".cache",
        )
    )
    result = pipeline.run_and_export(
        PipelineInput(
            sequence=example["reference_sequence"],
            sequence_label=example["identifier"],
            mutant_sequence=example["mutant_sequence"],
        ),
        output_dir=args.output_dir,
    )

    print("Demo synthetique hors ligne : substitution C12A (positions depuis 1).")
    print(f"WT     : {result.sequence_record.sequence}")
    print(f"         {result.structure.dot_bracket}")
    print(f"Mutant : {result.mutant_sequence}")
    print(f"         {result.mutant_structure.dot_bracket}")
    print("Positions les plus modifiees (delta = score mutant - score WT) :")
    for hotspot in result.impact_summary["hotspots"]:
        print(f"- position {hotspot['position']}: {hotspot['delta']:+.3f}")
    print("Scores heuristiques, sans interpretation clinique.")
    print(f"Rapport Markdown : {result.report_paths.markdown_path}")
    print(f"Rapport JSON     : {result.report_paths.json_path}")
    print(f"Artefacts        : {result.visualization_paths.root_dir}")


if __name__ == "__main__":
    main()
