"""Shared command entry points for alignment-only and synthetic benchmarks."""
import argparse
import json
from pathlib import Path

from .benchmark import generate_test_pair, measure_alignment
from .config import add_alignment_arguments, config_from_args, load_sequence_inputs
from deepstructgenomics.data_sources.ncbi_client import load_fasta_record, normalize_rna_sequence


def benchmark_main():
    return main(benchmark=True)


def main(argv=None, *, benchmark=False):
    parser = argparse.ArgumentParser(description="Alignement WT/MUT seul ; génération synthétique sans calcul de structure ARN.")
    add_alignment_arguments(parser)
    parser.add_argument("--show-config", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--length" if benchmark else "--generate-length", dest="length", type=int)
    parser.add_argument("--substitutions", "--mutations", type=int, default=1)
    parser.add_argument("--insertions", type=int, default=0)
    parser.add_argument("--deletions", type=int, default=0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--pattern", choices=("repetitive", "deterministic-random"), default="deterministic-random")
    for flag in ("wt", "mut", "wt-fasta", "mut-fasta"):
        parser.add_argument("--" + flag)
    args = parser.parse_args(argv)
    try:
        config = config_from_args(args)
        if args.show_config:
            print(json.dumps(config.to_dict(), indent=2))
            if not args.length and not any((args.wt, args.mut, args.wt_fasta, args.mut_fasta)):
                return 0
        if args.length is not None:
            if any((args.wt, args.mut, args.wt_fasta, args.mut_fasta)):
                parser.error("Ne pas mélanger génération synthétique et données utilisateur.")
            wt, mut = generate_test_pair(args.length, args.substitutions, args.insertions, args.deletions, args.seed, args.pattern)
        else:
            if benchmark:
                parser.error("Le benchmark requiert --length (données synthétiques uniquement).")
            inputs = load_sequence_inputs(args.env_file, cli={"sequence": args.wt, "mutant_sequence": args.mut,
                                                              "fasta_path": args.wt_fasta, "mutant_fasta_path": args.mut_fasta})
            def resolve(seq, fasta):
                if inputs.get(fasta):
                    return load_fasta_record(inputs[fasta]).sequence
                if not inputs.get(seq):
                    raise ValueError("Fournir WT et MUT : séquences, FASTA ou configuration externe.")
                return normalize_rna_sequence(inputs[seq])
            wt, mut = resolve("sequence", "fasta_path"), resolve("mutant_sequence", "mutant_fasta_path")
        summary, result = measure_alignment(wt, mut, config)
        summary["synthetic"] = args.length is not None
        if summary["synthetic"]:
            summary["generation"] = {key: getattr(args, key) for key in ("length", "substitutions", "insertions", "deletions", "seed", "pattern")}
        print("DeepStructGenomics — " + ("BENCHMARK SYNTHÉTIQUE (aucune structure ARN)" if summary["synthetic"] else "Alignement WT/MUT"))
        for key, value in summary.items():
            if key not in ("configuration", "generation"):
                print(f"{key:18} : {value}")
        print("Mémoire : pic du processus entier en MiB, imports inclus ; cells = WT × MUT, pas des octets.")
        if not result:
            print("Alignement non calculé. Une comparaison brute par position peut décaler les indels.")
        if args.verbose and result:
            print(result.aligned_reference, result.aligned_mutant, sep="\n")
        if args.json_output:
            args.json_output.parent.mkdir(parents=True, exist_ok=True)
            args.json_output.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        return 2 if summary["status"] == "rejected" else 0
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
