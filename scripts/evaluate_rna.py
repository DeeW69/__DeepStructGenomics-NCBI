#!/usr/bin/env python3
"""Compare RNA predictions against a JSON reference dataset."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from deepstructgenomics.rna.evaluation import evaluate_reference_set


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = evaluate_reference_set(args.reference)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    except (ValueError, OSError) as exc:
        parser.exit(2, f"Erreur : {exc}\n")
    print(args.output)


if __name__ == "__main__":
    main()
