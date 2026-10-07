#!/usr/bin/env python3
"""Build a reproducible bpRNA PDB subset from the original stFiles.zip archive."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random
import re
import sys
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from deepstructgenomics.rna.evaluation import dot_bracket_pairs

SOURCE_URL = "https://bprna.cgrb.oregonstate.edu/bpRNA_1m/stFiles.zip"


def prepare(archive, output, *, size=30, seed=45):
    if size < 1:
        raise ValueError("Taille du corpus positive requise.")
    archive = Path(archive)
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    candidates, seen, exclusions = [], set(), Counter()
    with zipfile.ZipFile(archive) as z:
        names = sorted((n for n in z.namelist() if re.fullmatch(r"(?:.*/)?bpRNA_PDB_\d+\.st", n)),
                       key=lambda n: int(re.search(r"_(\d+)\.st$", n)[1]))
        for name in names:
            raw = z.read(name)
            lines = [line.strip() for line in raw.decode("utf-8").splitlines() if line and not line.startswith("#")]
            if len(lines) < 2:
                raise ValueError(f"Fichier incomplet : {name}")
            sequence, structure = lines[:2]
            reason = None
            if not 20 <= len(sequence) <= 120:
                reason = "length_outside_20_120"
            elif set(sequence) - set("ACGU"):
                reason = "non_ACGU_sequence"
            elif set(structure) - set(".()"):
                reason = "pseudoknot_or_extended_notation"
            elif sequence in seen:
                reason = "duplicate_sequence"
            if reason:
                exclusions[reason] += 1
                continue
            dot_bracket_pairs(structure, len(sequence))
            seen.add(sequence)
            candidates.append({"id": Path(name).stem, "sequence": sequence, "reference": structure,
                               "provenance": {"archive_member": name, "sha256": hashlib.sha256(raw).hexdigest(),
                                              "origin": "PDB-derived secondary annotation in bpRNA-1m v1.0"}})
    if size > len(candidates):
        raise ValueError(f"Seulement {len(candidates)} references eligibles.")
    records = sorted(random.Random(seed).sample(candidates, size), key=lambda r: int(r["id"].rsplit("_", 1)[1]))
    dataset = {"source": {"database": "bpRNA-1m v1.0, PDB subset", "url": SOURCE_URL,
                          "publication": "https://doi.org/10.1093/nar/gky285", "archive_sha256": digest},
               "selection": {"seed": seed, "sample_size": size, "pdb_records": len(names),
                             "eligible_unique_sequences": len(candidates), "excluded": dict(exclusions),
                             "policy": "Numeric ID order; 20-120 nt; ACGU; .() only; exact sequence deduplication; seeded random sample before prediction.",
                             "limitations": "Small experimental-structure-derived sample, not family-balanced or homology-filtered; no claim of independence from ViennaRNA parameter fitting."},
               "records": records}
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(dataset, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return dataset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--size", type=int, default=30)
    parser.add_argument("--seed", type=int, default=45)
    args = parser.parse_args()
    try:
        dataset = prepare(args.archive, args.output, size=args.size, seed=args.seed)
    except (ValueError, OSError, zipfile.BadZipFile) as exc:
        parser.exit(2, f"Erreur : {exc}\n")
    print(json.dumps(dataset["selection"], indent=2))


if __name__ == "__main__":
    main()
