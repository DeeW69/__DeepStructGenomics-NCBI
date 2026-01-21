"""Simple RNA secondary structure module (Nussinov-like baseline)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from deepstructgenomics.config import RNAConfig

PAIRINGS = {
    ("A", "U"),
    ("U", "A"),
    ("G", "C"),
    ("C", "G"),
    ("G", "U"),
    ("U", "G"),
}


@dataclass
class SecondaryStructureResult:
    """Container for RNA structure predictions."""

    sequence: str
    dot_bracket: str
    base_pairs: Sequence[Tuple[int, int]]
    stability_score: float
    metrics: Dict[str, float]


def predict_secondary_structure(sequence: str, config: RNAConfig | None = None) -> SecondaryStructureResult:
    """Run a simplified Nussinov algorithm to produce a dot-bracket structure."""

    cfg = config or RNAConfig()
    sequence = sequence.upper().replace("T", "U")
    n = len(sequence)
    if n == 0:
        raise ValueError("La sequence fournie est vide.")

    dp = [[0.0 for _ in range(n)] for _ in range(n)]
    decision: List[List[Tuple[str, int] | None]] = [[None for _ in range(n)] for _ in range(n)]

    for length in range(1, n):
        for i in range(n - length):
            j = i + length
            best = dp[i + 1][j]
            decision[i][j] = ("skip_i", i + 1)

            if dp[i][j - 1] > best:
                best = dp[i][j - 1]
                decision[i][j] = ("skip_j", j - 1)

            if (j - i) > cfg.min_loop_length and can_pair(sequence[i], sequence[j]):
                pair_score = dp[i + 1][j - 1] + pairing_score(sequence[i], sequence[j], cfg)
                if pair_score > best:
                    best = pair_score
                    decision[i][j] = ("pair", j)

            for k in range(i + 1, j):
                combined = dp[i][k] + dp[k + 1][j]
                if combined > best:
                    best = combined
                    decision[i][j] = ("split", k)

            dp[i][j] = best

    base_pairs = _traceback(decision, sequence, cfg)
    dot_bracket = _pairs_to_dot_bracket(n, base_pairs)
    metrics = _compute_structure_metrics(sequence, base_pairs)
    stability_score = dp[0][n - 1] if n > 1 else 0.0

    return SecondaryStructureResult(
        sequence=sequence,
        dot_bracket=dot_bracket,
        base_pairs=base_pairs,
        stability_score=stability_score,
        metrics=metrics,
    )


def _traceback(
    decisions: List[List[Tuple[str, int] | None]],
    sequence: str,
    config: RNAConfig,
) -> List[Tuple[int, int]]:
    """Recover the base pairs from the decision matrix."""

    pairs: List[Tuple[int, int]] = []

    def walk(i: int, j: int):
        if i >= j or decisions[i][j] is None:
            return
        action, value = decisions[i][j]
        if action == "skip_i":
            walk(i + 1, j)
        elif action == "skip_j":
            walk(i, j - 1)
        elif action == "pair":
            if (j - i) > config.min_loop_length and can_pair(sequence[i], sequence[value]):
                pairs.append((i, value))
                walk(i + 1, value - 1)
            walk(value + 1, j)
        elif action == "split":
            walk(i, value)
            walk(value + 1, j)

    walk(0, len(sequence) - 1)
    return sorted(pairs)


def _pairs_to_dot_bracket(length: int, pairs: Sequence[Tuple[int, int]]) -> str:
    """Convert list of base pairs to dot-bracket notation."""

    structure = ["." for _ in range(length)]
    for i, j in pairs:
        structure[i] = "("
        structure[j] = ")"
    return "".join(structure)


def _compute_structure_metrics(sequence: str, base_pairs: Sequence[Tuple[int, int]]) -> Dict[str, float]:
    """Derive simple descriptive metrics."""

    n = len(sequence)
    gc_count = sum(1 for base in sequence if base in {"G", "C"})
    metrics = {
        "length": float(n),
        "gc_content": gc_count / n if n else 0.0,
        "paired_fraction": (2 * len(base_pairs)) / n if n else 0.0,
    }

    # Estimate loop count by counting '.' segments in dot-bracket
    dot_bracket = _pairs_to_dot_bracket(n, base_pairs)
    loop_count = sum(1 for segment in dot_bracket.split("(") for part in segment.split(")") if part and all(ch == "." for ch in part))
    metrics["loop_count"] = float(loop_count)
    return metrics


def can_pair(base_a: str, base_b: str) -> bool:
    """Return True if the bases can pair under canonical/au/gu rules."""

    return (base_a, base_b) in PAIRINGS


def pairing_score(base_a: str, base_b: str, config: RNAConfig) -> float:
    """Heuristic energy contribution for a pairing."""

    pair = (base_a, base_b)
    if pair in {("G", "C"), ("C", "G")}:
        return config.gc_weight
    if pair in {("A", "U"), ("U", "A")}:
        return config.au_weight
    return config.gu_weight
