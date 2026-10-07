"""Reproducible base-pair evaluation against explicit reference structures."""

from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

from deepstructgenomics.config import RNAConfig
from deepstructgenomics.data_sources.ncbi_client import normalize_rna_sequence
from deepstructgenomics.rna.secondary_structure import can_pair, predict_secondary_structure


def dot_bracket_pairs(structure: str, length: int):
    if not isinstance(structure, str) or len(structure) != length:
        raise ValueError("Longueur de structure incompatible avec la sequence.")
    stack, pairs = [], set()
    for index, symbol in enumerate(structure):
        if symbol == "(":
            stack.append(index)
        elif symbol == ")" and stack:
            pairs.add((stack.pop(), index))
        elif symbol != ".":
            raise ValueError("Dot-bracket invalide ; seuls .() sans pseudonoeuds sont acceptes.")
    if stack:
        raise ValueError("Parentheses non equilibrees.")
    return pairs


def pair_metrics(reference, prediction):
    tp = len(reference & prediction)
    fp, fn = len(prediction - reference), len(reference - prediction)
    return {"tp": tp, "fp": fp, "fn": fn,
            "precision": tp / (tp + fp) if tp + fp else None,
            "recall": tp / (tp + fn) if tp + fn else None,
            "f1": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else None}


def greedy_pairs(sequence, min_loop_length=3):
    """Leftmost base paired with the rightmost compatible base, recursively nested."""
    pairs, intervals = set(), [(0, len(sequence) - 1)]
    while intervals:
        left, right = intervals.pop()
        if left >= right:
            continue
        partner = next((j for j in range(right, left + min_loop_length, -1)
                        if can_pair(sequence[left], sequence[j])), None)
        if partner is None:
            intervals.append((left + 1, right))
        else:
            pairs.add((left, partner))
            intervals.extend([(left + 1, partner - 1), (partner + 1, right)])
    return pairs


def evaluate_reference_set(path: str | Path, *, config=None, viennarna=False):
    """Evaluate Nussinov, greedy and all-unpaired controls plus supplied predictions.

    Every record must contain all external methods, preventing unequal comparisons.
    Undefined metrics are JSON null, including pair-free reference/prediction pairs.
    """
    raw = Path(path).read_bytes()
    dataset = json.loads(raw)
    if not isinstance(dataset, dict) or not isinstance(dataset.get("records"), list) or not dataset["records"]:
        raise ValueError("Jeu de reference vide ou invalide.")
    if not dataset.get("source"):
        raise ValueError("La provenance source du jeu est requise.")
    config = config or RNAConfig()
    predictor = None
    if viennarna:
        from deepstructgenomics.rna.vienna import ViennaPredictor
        predictor = ViennaPredictor()
    method_provenance = dict(dataset.get("method_provenance", {}))
    if predictor:
        method_provenance["viennarna_mfe"] = predictor.provenance
    records, seen, external_methods = [], set(), None
    totals = {}
    for record in dataset["records"]:
        if not isinstance(record, dict) or not isinstance(record.get("id"), str) or not record["id"]:
            raise ValueError("Identifiant de reference requis.")
        identifier = record["id"]
        if identifier in seen:
            raise ValueError(f"Identifiant duplique : {identifier}")
        seen.add(identifier)
        if not isinstance(record.get("sequence"), str):
            raise ValueError(f"Sequence requise : {identifier}")
        sequence = normalize_rna_sequence(record["sequence"])
        reference = dot_bracket_pairs(record.get("reference"), len(sequence))
        external = record.get("predictions", {})
        if not isinstance(external, dict):
            raise ValueError("predictions doit etre un objet methode:dot-bracket.")
        if external_methods is None:
            external_methods = set(external)
        if set(external) != external_methods:
            raise ValueError("Chaque methode externe doit couvrir toutes les references.")
        predictions = {
            "nussinov": set(predict_secondary_structure(sequence, config).base_pairs),
            "greedy": greedy_pairs(sequence, config.min_loop_length),
            "unpaired": set(),
        }
        if set(external) & predictions.keys():
            raise ValueError("Nom de methode externe reserve.")
        energy = None
        if predictor:
            if "viennarna_mfe" in external:
                raise ValueError("viennarna_mfe est reserve a la methode executee.")
            structure, energy = predictor.predict(sequence)
            predictions["viennarna_mfe"] = dot_bracket_pairs(structure, len(sequence))
        predictions.update({name: dot_bracket_pairs(value, len(sequence)) for name, value in external.items()})
        metrics = {}
        for name, pairs in predictions.items():
            metrics[name] = pair_metrics(reference, pairs)
            counts = totals.setdefault(name, {"tp": 0, "fp": 0, "fn": 0})
            for key in counts:
                counts[key] += metrics[name][key]
        structures = {}
        for name, pairs in predictions.items():
            symbols = ["."] * len(sequence)
            for left, right in pairs:
                symbols[left], symbols[right] = "(", ")"
            structures[name] = "".join(symbols)
        records.append({"id": identifier, "length": len(sequence), "methods": metrics,
                        "sequence": sequence, "reference": record["reference"],
                        "provenance": record.get("provenance", {}),
                        "predictions": structures, "viennarna_mfe_kcal_mol": energy})
    micro = {}
    for name, counts in totals.items():
        tp, fp, fn = (counts[key] for key in ("tp", "fp", "fn"))
        micro[name] = {**counts,
                       "precision": tp / (tp + fp) if tp + fp else None,
                       "recall": tp / (tp + fn) if tp + fn else None,
                       "f1": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else None}
    return {"source": dataset["source"], "dataset_sha256": hashlib.sha256(raw).hexdigest(),
            "generated_at": datetime.now(timezone.utc).isoformat(), "parameters": asdict(config),
            "method_provenance": method_provenance,
            "records": records, "micro_average": micro,
            "metric_policy": "Exact base pairs; undefined ratios are null; no pseudoknots."}
