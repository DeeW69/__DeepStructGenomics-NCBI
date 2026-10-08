"""Bounded global affine alignment; no Qt, no alternate-traceback enumeration."""

import Bio
from Bio.Align import PairwiseAligner

from deepstructgenomics.data_sources.ncbi_client import normalize_rna_sequence
from .models import AlignmentResult
from .config import AlignmentConfig, alignment_metadata, overflow_reason


def align_sequences(reference: str, mutant: str, config: AlignmentConfig | None = None) -> AlignmentResult | None:
    config = config or AlignmentConfig()
    reference, mutant = normalize_rna_sequence(reference), normalize_rna_sequence(mutant)
    if not config.enabled:
        return None
    reason = overflow_reason(config, len(reference), len(mutant))
    if reason:
        if config.overflow_policy == "positional":
            return None
        raise ValueError(reason)
    aligner = PairwiseAligner(mode="global", match_score=config.match_score, mismatch_score=config.mismatch_score,
                             open_gap_score=config.gap_open_score, extend_gap_score=config.gap_extend_score)
    aligner.epsilon = 1e-6
    candidates = iter(aligner.align(reference[::-1], mutant[::-1]))
    best = next(candidates)
    ambiguous = next(candidates, None) is not None
    warnings = ("Plusieurs alignements optimaux : la position des indels est conventionnelle, notamment dans les répétitions.",) if ambiguous else ()
    return AlignmentResult(str(best[0])[::-1], str(best[1])[::-1], best.score, config.scoring, Bio.__version__, warnings,
                           alignment_metadata(config, len(reference), len(mutant), "aligned"), ambiguous)
