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


def compare_sequences(
    reference_seq: str,
    mutant_seq: str,
    reference_pairs: Sequence[Tuple[int, int]],
    mutant_pairs: Sequence[Tuple[int, int]],
) -> VariantImpactResult:
    """Simple heuristic impact estimation baseline."""

    ref_len = len(reference_seq)
    mut_len = len(mutant_seq)
    max_len = max(ref_len, mut_len)
    substitutions: List[Dict[str, int | str]] = []

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
