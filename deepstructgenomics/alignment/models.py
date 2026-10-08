"""All exposed positions and columns are 1-based; a gap maps to None."""
from dataclasses import asdict, dataclass
from functools import cached_property
from itertools import groupby
import math


@dataclass(frozen=True)
class AlignmentColumn:
    index: int
    reference_position: int | None
    mutant_position: int | None
    reference_base: str | None
    mutant_base: str | None
    operation: str


@dataclass(frozen=True)
class AlignmentResult:
    aligned_reference: str
    aligned_mutant: str
    score: float
    parameters: dict
    implementation_version: str
    warnings: tuple[str, ...] = ()
    run_metadata: dict | None = None
    ambiguous: bool = False

    def __post_init__(self):
        if (not self.aligned_reference or len(self.aligned_reference) != len(self.aligned_mutant)
                or not math.isfinite(self.score)
                or any(base not in "ACGU-" for base in self.aligned_reference + self.aligned_mutant)
                or any(a == b == "-" for a, b in zip(self.aligned_reference, self.aligned_mutant))
                or not self.aligned_reference.replace("-", "") or not self.aligned_mutant.replace("-", "")):
            raise ValueError("Alignement invalide : bases, colonnes ou score incohérents.")

    @cached_property
    def columns(self) -> tuple[AlignmentColumn, ...]:
        ref, mut = 0, 0
        result = []
        for index, (a, b) in enumerate(zip(self.aligned_reference, self.aligned_mutant), 1):
            ref += a != "-"
            mut += b != "-"
            operation = "insertion" if a == "-" else "deletion" if b == "-" else "match" if a == b else "substitution"
            result.append(AlignmentColumn(index, ref if a != "-" else None, mut if b != "-" else None,
                                          a if a != "-" else None, b if b != "-" else None, operation))
        return tuple(result)

    @cached_property
    def ref_to_mut(self):
        return {c.reference_position: c.mutant_position for c in self.columns if c.reference_position is not None}

    @cached_property
    def mut_to_ref(self):
        return {c.mutant_position: c.reference_position for c in self.columns if c.mutant_position is not None}

    @cached_property
    def ref_to_alignment(self):
        return {c.reference_position: c.index for c in self.columns if c.reference_position is not None}

    @cached_property
    def mut_to_alignment(self):
        return {c.mutant_position: c.index for c in self.columns if c.mutant_position is not None}

    @cached_property
    def counts(self):
        return {kind: sum(c.operation == kind for c in self.columns) for kind in ("match", "substitution", "insertion", "deletion")}

    @property
    def identity(self):
        """Matches / all alignment columns, including gaps."""
        return self.counts["match"] / len(self.columns)

    @property
    def cigar(self):
        symbols = {"match": "=", "substitution": "X", "insertion": "I", "deletion": "D"}
        return "".join(f"{sum(1 for _ in group)}{symbols[kind]}" for kind, group in groupby(c.operation for c in self.columns))

    def to_dict(self):
        return {"schema": 1, "coordinates": "1-based", "algorithm": "global-affine",
                "implementation": "Bio.Align.PairwiseAligner", "implementation_version": self.implementation_version,
                "tie_break": "first optimal traceback on reversed sequences; reverse columns back",
                "parameters": self.parameters, "score": self.score, "identity": self.identity,
                "counts": self.counts, "cigar": self.cigar, "warnings": list(self.warnings),
                "run_metadata": self.run_metadata, "ambiguous": self.ambiguous,
                "columns": [asdict(column) for column in self.columns],
                "aligned_reference": self.aligned_reference, "aligned_mutant": self.aligned_mutant}

    @classmethod
    def from_dict(cls, payload, reference, mutant):
        try:
            if payload["schema"] != 1 or payload["coordinates"] != "1-based":
                raise ValueError
            result = cls(payload["aligned_reference"], payload["aligned_mutant"], float(payload["score"]),
                         payload["parameters"], payload["implementation_version"], tuple(payload["warnings"]),
                         payload.get("run_metadata"), payload.get("ambiguous", bool(payload["warnings"])))
            if result.run_metadata:
                from .config import validate_alignment_metadata, overflow_reason
                config = validate_alignment_metadata(result.run_metadata, len(reference), len(mutant))
                if (config.scoring != result.parameters or not config.enabled
                        or result.run_metadata["status"] != "aligned" or overflow_reason(config, len(reference), len(mutant))):
                    raise ValueError
            if (result.aligned_reference.replace("-", "") != reference or result.aligned_mutant.replace("-", "") != mutant
                    or result.counts != payload["counts"] or result.cigar != payload["cigar"]
                    or not math.isclose(result.identity, payload["identity"])
                    or payload["columns"] != [asdict(column) for column in result.columns]):
                raise ValueError
            return result
        except (KeyError, TypeError, ValueError, AttributeError) as exc:
            raise ValueError("Alignement enregistré incohérent avec les séquences. Rouvrez un résultat valide ou relancez l'analyse.") from exc
