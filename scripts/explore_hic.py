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
    parser.add_argument("--advanced", action="store_true", help="Equilibrage, tests et reconstruction d'une region cis.")
    parser.add_argument("--chromosome")
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    parser.add_argument("--missing-as-zero", action="store_true", help="Declarer que les contacts absents sont des zeros observes.")
    parser.add_argument("--window-bins", type=int, default=3)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=46)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--min-separation", type=int, default=2)
    parser.add_argument("--distance-exponent", type=float, default=1 / 3)
    parser.add_argument("--max-iterations", type=int, default=2000)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    try:
        if args.advanced:
            if args.start is None or args.end is None or not args.chromosome:
                raise ValueError("--advanced requiert --chromosome, --start et --end.")
            from deepstructgenomics.analysis.hic_advanced import analyze_region
            from deepstructgenomics.reporting.hic_report import export_hic_analysis
            result = analyze_region(args.contacts, assembly=args.assembly, chromosome=args.chromosome,
                                    start=args.start, end=args.end, bin_size=args.bin_size,
                                    missing_as_zero=args.missing_as_zero, window=args.window_bins,
                                    permutations=args.permutations, seed=args.seed, fdr=args.fdr,
                                    min_separation=args.min_separation, exponent=args.distance_exponent,
                                    tolerance=args.tolerance, max_iterations=args.max_iterations)
            # Avoid mixing stale successful tables with a later unsuccessful run.
            args.output_dir.mkdir(parents=True, exist_ok=False)
            export_hic_analysis(result, args.output_dir, plot=not args.no_plot)
            print(args.output_dir / "hic_analysis.json")
            if result["status"] != "ok":
                raise SystemExit(2)
            return
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
