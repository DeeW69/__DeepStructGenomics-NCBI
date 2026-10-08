"""Variant comparison utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple


@dataclass
class VariantImpactResult:
    """Summarizes the delta between a reference and mutant sequence/structure."""

    total_differences: int
    substitutions: List[Dict[str, int | str]]
    structure_delta_score: float
    commentary: str
    coordinate_system: str = "sequence_position_0based"


def compare_sequences(
    reference_seq: str,
    mutant_seq: str,
    reference_pairs: Sequence[Tuple[int, int]],
    mutant_pairs: Sequence[Tuple[int, int]],
    *, alignment=None,
) -> VariantImpactResult:
    """Simple heuristic impact estimation baseline."""

    ref_len = len(reference_seq)
    mut_len = len(mutant_seq)
    max_len = max(ref_len, mut_len)
    substitutions: List[Dict[str, int | str]] = []

    if alignment:
        from deepstructgenomics.alignment.mapping import ComparisonMapping
        changes = [{"position": c.index - 1, "reference": c.reference_base or "-", "mutant": c.mutant_base or "-"}
                   for c in alignment.columns if c.operation != "match"]
        pairs = ComparisonMapping(reference_seq, mutant_seq, alignment).classify_pairs(reference_pairs, mutant_pairs)
        delta = len(pairs["gained"]) + len(pairs["inserted"]) - len(pairs["lost"]) - len(pairs["deleted"])
        return VariantImpactResult(len(changes), changes, float(delta),
            "Comparaison des appariements via l'alignement ; ce delta n'est pas une énergie ni un effet biologique mesuré.",
            "alignment_column_0based")

    for idx in range(max_len):
        ref_base = reference_seq[idx] if idx < ref_len else "-"
        mut_base = mutant_seq[idx] if idx < mut_len else "-"
        if ref_base != mut_base:
            substitutions.append({"position": idx, "reference": ref_base, "mutant": mut_base})

    ref_pair_set = {tuple(sorted(pair)) for pair in reference_pairs}
    mut_pair_set = {tuple(sorted(pair)) for pair in mutant_pairs}
    lost_pairs = len(ref_pair_set - mut_pair_set)
    gained_pairs = len(mut_pair_set - ref_pair_set)
    structure_delta_score = gained_pairs - lost_pairs

    commentary = (
        "Variant potentiellement stabilisant."
        if structure_delta_score > 0
        else "Variant potentiellement destabilisant."
        if structure_delta_score < 0
        else "Impact structural neutre au premier ordre."
    )

    return VariantImpactResult(
        total_differences=len(substitutions),
        substitutions=substitutions,
        structure_delta_score=float(structure_delta_score),
        commentary=commentary,
    )
