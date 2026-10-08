"""Single mapping adapter used by reports and views. Pairs supplied in base 0."""
from .models import AlignmentColumn, AlignmentResult


def alignment_for(data):
    payload = data.get("alignment")
    if payload is None:
        return None
    return AlignmentResult.from_dict(payload, data["reference"]["sequence"], data["mutant"]["sequence"])


class ComparisonMapping:
    def __init__(self, reference, mutant=None, alignment=None):
        self.alignment = alignment
        if alignment:
            self.columns = alignment.columns
        else:
            self.columns = tuple(AlignmentColumn(i + 1, i + 1 if i < len(reference) else None,
                i + 1 if mutant is not None and i < len(mutant) else None,
                reference[i] if i < len(reference) else None, mutant[i] if mutant is not None and i < len(mutant) else None,
                "positional") for i in range(max(len(reference), len(mutant or ""))))
        self.to_column = {
            "reference": {c.reference_position: c.index for c in self.columns if c.reference_position is not None},
            "mutant": {c.mutant_position: c.index for c in self.columns if c.mutant_position is not None}}

    @classmethod
    def from_data(cls, data):
        if "_comparison_mapping" in data:
            return data["_comparison_mapping"]
        return cls(data["reference"]["sequence"], data["mutant"]["sequence"] if data.get("mutant") else None, alignment_for(data))

    def pairs(self, pairs, side):
        mapping = self.to_column[side]
        return {(mapping[i + 1], mapping[j + 1]) for i, j in pairs}

    def classify_pairs(self, reference_pairs, mutant_pairs):
        wt, mut = self.pairs(reference_pairs, "reference"), self.pairs(mutant_pairs, "mutant")
        removed = {p for p in wt - mut if any(self.columns[c - 1].mutant_position is None for c in p)} if self.alignment else set()
        inserted = {p for p in mut - wt if any(self.columns[c - 1].reference_position is None for c in p)} if self.alignment else set()
        return {"conserved": wt & mut, "lost": wt - mut - removed, "gained": mut - wt - inserted,
                "deleted": removed, "inserted": inserted}

    def deltas(self, wt_scores, mut_scores):
        return {c.index: mut_scores[c.mutant_position - 1] - wt_scores[c.reference_position - 1]
                for c in self.columns if c.reference_position is not None and c.mutant_position is not None}
