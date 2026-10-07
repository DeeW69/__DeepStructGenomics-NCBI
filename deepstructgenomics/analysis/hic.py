"""Exploratory sparse contact summaries; no 3D reconstruction or bias correction."""

import csv
from datetime import datetime, timezone
import hashlib
import io
import math
from pathlib import Path


def summarize_contacts(path: str | Path, *, assembly: str, bin_size: int, top_k: int = 10):
    """Read binned TSV: chrom1,start1,chrom2,start2,count (tab separated).

    Starts are zero-based, aligned to bin_size; bins are half-open intervals.
    Duplicate or mirrored bin pairs are rejected to prevent accidental doubling.
    """
    if not isinstance(assembly, str) or not assembly.strip():
        raise ValueError("Assemblage genomique requis.")
    if bin_size <= 0 or top_k < 0:
        raise ValueError("bin_size doit etre positif et top_k non negatif.")
    raw = Path(path).read_bytes()
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
    columns = ["chrom1", "start1", "chrom2", "start2", "count"]
    if reader.fieldnames != columns:
        raise ValueError("Colonnes TSV attendues : " + ", ".join(columns))
    contacts, seen, coverage, distances = [], set(), {}, {}
    total = cis = trans = diagonal = 0.0
    for line, row in enumerate(reader, 2):
        try:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Nombre de colonnes incorrect")
            a = (row["chrom1"].strip(), int(row["start1"]))
            b = (row["chrom2"].strip(), int(row["start2"]))
            count = float(row["count"])
            if (not a[0] or not b[0] or a[1] < 0 or b[1] < 0
                    or a[1] % bin_size or b[1] % bin_size
                    or not math.isfinite(count) or count < 0):
                raise ValueError("Coordonnees ou compte invalides")
            a, b = sorted((a, b))
            if (a, b) in seen:
                raise ValueError("Contact duplique ou symetrique")
            seen.add((a, b))
        except (ValueError, TypeError) as exc:
            raise ValueError(f"Ligne {line}: {exc}") from exc
        total += count
        if a[0] == b[0]:
            cis += count
            distance = abs(a[1] - b[1])
            distances[distance] = distances.get(distance, 0.0) + count
            if a == b:
                diagonal += count
        else:
            trans += count
        for anchor in set((a, b)):
            coverage[anchor] = coverage.get(anchor, 0.0) + count
        contacts.append({"chrom1": a[0], "start1": a[1], "chrom2": b[0], "start2": b[1], "count": count})
    if not contacts:
        raise ValueError("Aucun contact dans le fichier.")
    if not math.isfinite(total):
        raise ValueError("Somme des contacts hors limites numeriques.")
    top = sorted(contacts, key=lambda row: (-row["count"], row["chrom1"], row["start1"], row["chrom2"], row["start2"]))[:top_k]
    return {
        "source": str(path), "sha256": hashlib.sha256(raw).hexdigest(),
        "assembly": assembly.strip(), "bin_size": bin_size,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "coordinates": "0-based, half-open bins [start, start + bin_size)",
        "normalization": "none; input weights retained; fractions are descriptive only",
        "contact_rows": len(contacts), "total_weight": total,
        "cis_weight": cis, "trans_weight": trans, "diagonal_weight": diagonal,
        "cis_fraction": cis / total if total else None,
        "top_contacts": top,
        "bin_coverage": [{"chrom": key[0], "start": key[1], "weight": value}
                         for key, value in sorted(coverage.items())],
        "cis_distance_weights": [{"distance_bp": key, "weight": value} for key, value in sorted(distances.items())],
        "limitations": "Exploration only: no ICE/KR balancing, statistical loop calling, TAD detection or 3D reconstruction. Diagonal counted once in coverage.",
    }
