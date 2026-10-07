#!/usr/bin/env python3
"""Export descriptive Hi-C contact summaries from binned TSV data."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from deepstructgenomics.analysis.hic import summarize_contacts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contacts", type=Path, required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--bin-size", type=int, required=True)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = summarize_contacts(args.contacts, assembly=args.assembly,
                                    bin_size=args.bin_size, top_k=args.top_k)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "hic_summary.json").write_text(
            json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
        lines = ["# Exploration Hi-C", "", f"Assemblage : {result['assembly']}",
                 f"Taille des bins : {args.bin_size} bp", f"Contacts : {result['contact_rows']}",
                 f"Poids total : {result['total_weight']}", f"Poids cis : {result['cis_weight']}",
                 f"Poids trans : {result['trans_weight']}", "", result["limitations"], "",
                 "| Chromosome 1 | Début 1 | Chromosome 2 | Début 2 | Poids |",
                 "| --- | --- | --- | --- | --- |"]
        lines += [f"| {r['chrom1']} | {r['start1']} | {r['chrom2']} | {r['start2']} | {r['count']} |"
                  for r in result["top_contacts"]]
        (args.output_dir / "hic_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    except (ValueError, OSError) as exc:
        parser.exit(2, f"Erreur : {exc}\n")
    print(args.output_dir / "hic_summary.json")


if __name__ == "__main__":
    main()
