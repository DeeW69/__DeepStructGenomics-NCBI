"""Optional ViennaRNA MFE predictor, isolated from the default pipeline."""

from importlib import import_module


class ViennaPredictor:
    def __init__(self):
        try:
            self.rna = import_module("RNA")
        except ImportError as exc:
            raise ValueError("ViennaRNA absent : installer requirements-research-lock.txt.") from exc
        self.model = self.rna.md()
        self.model.temperature = 37.0
        self.model.dangles = 2
        self.model.noLP = 0
        self.provenance = {
            "software": "ViennaRNA", "version": self.rna.__version__,
            "algorithm": "fold_compound.mfe", "temperature_celsius": 37.0,
            "dangles": 2, "noLP": False, "parameters": "library default Turner 2004",
            "constraints": "none", "executed": True,
        }

    def predict(self, sequence):
        structure, energy = self.rna.fold_compound(sequence, self.model).mfe()
        return structure, float(energy)
