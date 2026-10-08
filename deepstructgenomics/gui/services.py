"""Desktop operations delegated to the existing scientific pipeline (no Qt)."""
from datetime import datetime, timezone
import json
from pathlib import Path
from uuid import uuid4
import hashlib

from deepstructgenomics.config import NCBIConfig, PipelineConfig
from deepstructgenomics.pipeline import DeepStructPipeline, PipelineInput
from deepstructgenomics.data_sources.ncbi_client import NCBIClient
from deepstructgenomics.visualization.secondary_view import load_secondary_data
from deepstructgenomics.alignment.mapping import ComparisonMapping


def new_run_directory(workspace, kind):
    root = Path(workspace).expanduser().resolve()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
    directory = root / f"{kind}_{stamp}_{uuid4().hex[:8]}"
    directory.mkdir(parents=True, exist_ok=False)
    return directory


def run_rna(workspace, parameters):
    parameters = dict(parameters)
    email = parameters.pop("ncbi_email", None)
    expected = parameters.pop("expected_ncbi_sha256", None)
    directory = new_run_directory(workspace, "rna")
    pipeline = DeepStructPipeline(PipelineConfig(
        ncbi=NCBIConfig(email=email), cache_dir=Path(workspace).resolve() / "cache",
        default_output_dir=directory,
    ))
    if expected:
        record = pipeline.ncbi_client.fetch_sequence(parameters["accession"], refresh=parameters.get("refresh_cache", False))
        if hashlib.sha256(record.sequence.encode()).hexdigest() != expected:
            raise ValueError("La séquence NCBI a changé depuis l'aperçu. Sélectionner de nouveau la référence.")
        parameters["refresh_cache"] = False
    result = pipeline.run_and_export(PipelineInput(**parameters), directory)
    return {"kind": "rna", "path": str(result.visualization_paths.manifest_path),
            "directory": str(directory), "label": result.sequence_record.identifier}


def search_ncbi(workspace, parameters):
    parameters = dict(parameters)
    client = NCBIClient(NCBIConfig(email=parameters.pop("ncbi_email", None)))
    return client.search_sequences(**parameters)


def preview_ncbi(workspace, parameters):
    client = NCBIClient(NCBIConfig(email=parameters.get("ncbi_email")), cache_dir=Path(workspace).resolve() / "cache")
    accession = parameters["accession"]
    record = client.fetch_sequence(accession)
    return {"accession": accession, "identifier": record.identifier, "description": record.description,
            "sequence": record.sequence, "metadata": record.metadata,
            "sha256": hashlib.sha256(record.sequence.encode()).hexdigest()}


def run_hic(workspace, parameters):
    from deepstructgenomics.analysis.hic_advanced import analyze_region
    from deepstructgenomics.reporting.hic_report import export_hic_analysis
    report = analyze_region(**parameters)
    directory = new_run_directory(workspace, "hic")
    # Qt/Matplotlib rendering stays on the main thread.
    export_hic_analysis(report, directory, plot=False)
    return {"kind": "hic", "path": str(directory / "hic_analysis.json"),
            "directory": str(directory), "label": report["region"]["chromosome"]}


def read_rna(path):
    """Consume the same exported values as the standalone secondary viewer."""
    data = load_secondary_data(path)
    wt, mut = data["reference"], data.get("mutant")
    wt_pairs = set(map(tuple, wt["base_pairs"]))
    mut_pairs = set(map(tuple, mut["base_pairs"])) if mut else set()
    mapping = ComparisonMapping.from_data(data)
    pairs = mapping.classify_pairs(wt["base_pairs"], mut["base_pairs"] if mut else [])
    delta = list(mapping.deltas(wt["scores"], mut["scores"]).values()) if mut else []
    metrics = {
        "wt_pairs": len(wt_pairs), "mut_pairs": len(mut_pairs) if mut else None,
        "lost": len(pairs["lost"]) + len(pairs["deleted"]) if mut else None,
        "gained": len(pairs["gained"]) + len(pairs["inserted"]) if mut else None,
        "deleted_pairs": len(pairs["deleted"]) if mut else None,
        "inserted_pairs": len(pairs["inserted"]) if mut else None,
        "changed": sum(abs(value) > 1e-9 for value in delta) if mut else None,
        "mean_delta": sum(delta) / len(delta) if delta else None,
        "max_abs_delta": max(map(abs, delta)) if delta else None,
    }
    return data, metrics


def position_details(data, position):
    """Inspect one 1-based alignment column, or explicit positional fallback."""
    result = {"position": position}
    mapping = ComparisonMapping.from_data(data)
    column = mapping.columns[position - 1]
    for key in ("reference", "mutant"):
        structure = data.get(key)
        native = column.reference_position if key == "reference" else column.mutant_position
        if not structure or native is None:
            result[key] = None
            continue
        partner = next((j + 1 if i == native - 1 else i + 1
                        for i, j in structure["base_pairs"] if native - 1 in (i, j)), None)
        result[key] = {"base": structure["sequence"][native - 1], "position": native,
                       "score": structure["scores"][native - 1], "partner": partner}
    result["operation"] = column.operation
    wt, mut = result["reference"], result["mutant"]
    result["delta"] = mut["score"] - wt["score"] if wt and mut else None
    return result


def inspection_details(data, position, radius=4):
    """Describe positional changes; absence never implies a lost pair or an indel."""
    detail = position_details(data, position)
    mapping = ComparisonMapping.from_data(data)
    changes = []
    wt, mut = detail["reference"], detail["mutant"]
    for key in ("reference", "mutant"):
        entry = detail[key]
        if entry:
            sequence = data[key]["sequence"]
            partner = entry["partner"]
            entry["partner_label"] = f"{sequence[partner - 1]}{partner}" if partner else "—"
    if wt and mut:
        def pair_text(entry):
            return f"{entry['base']}{entry['position']} ↔ {entry['partner_label']}"
        wp = mapping.to_column["reference"].get(wt["partner"])
        mp = mapping.to_column["mutant"].get(mut["partner"])
        if wp is not None and wp == mp:
            changes.append(("conserved", "Conservée", pair_text(wt) + " / " + pair_text(mut)))
        else:
            if wt["partner"] is not None:
                deleted = mapping.alignment and mapping.columns[wp - 1].mutant_position is None
                changes.append(("lost", "Paire supprimée avec la base" if deleted else "Paire perdue", pair_text(wt) + " disparaît."))
            if mut["partner"] is not None:
                inserted = mapping.alignment and mapping.columns[mp - 1].reference_position is None
                changes.append(("gained", "Nouvelle paire avec insertion" if inserted else "Nouvelle paire", pair_text(mut) + " apparaît."))
        if wt["base"] != mut["base"]:
            same_length = mapping.alignment or len(data["reference"]["sequence"]) == len(data["mutant"]["sequence"])
            changes.append(("substitution", f"Substitution {wt['base']} → {mut['base']}" if same_length
                            else f"Différence {wt['base']} → {mut['base']}",
                            "Bases mises en correspondance par alignement." if mapping.alignment else "Comparaison à la même position, sans alignement."))
    elif mapping.alignment:
        entry = wt or mut
        changes.append(("lost" if wt else "gained", "Délétion" if wt else "Insertion",
                        f"{entry['base']}{entry['position']} · pas d'homologue."))
        if entry["partner"] is not None:
            changes.append(("lost" if wt else "gained", "Paire supprimée avec la base" if wt else "Nouvelle paire avec insertion",
                            f"{entry['base']}{entry['position']} ↔ {entry['partner_label']}"))
    detail["changes"] = changes
    start = max(1, position - radius)
    end = min(len(mapping.columns), position + radius)
    detail["local_start"], detail["local_end"] = start, end
    gap = "-" if mapping.alignment else "—"
    detail["local"] = {"reference": [c.reference_base or gap for c in mapping.columns[start - 1:end]],
                       "mutant": [c.mutant_base or gap for c in mapping.columns[start - 1:end]]}
    return detail


class RecentStore:
    """Only remember result paths; never persist input sequences or credentials."""
    def __init__(self, path):
        self.path = Path(path)
        self.warning = ""

    def load(self):
        self.warning = ""
        try:
            entries = json.loads(self.path.read_text(encoding="utf-8"))
            if not isinstance(entries, list):
                raise ValueError("Historique invalide")
            return [row for row in entries if isinstance(row, dict)
                    and row.get("kind") in ("rna", "hic")
                    and all(isinstance(row.get(k), str) for k in ("path", "label", "directory"))][:20]
        except FileNotFoundError:
            return []
        except (OSError, ValueError):
            self.warning = "Historique illisible ; les fichiers de résultats restent disponibles."
            return []

    def add(self, entry):
        entries = [entry] + [row for row in self.load() if row["path"] != entry["path"]]
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(".tmp")
        temporary.write_text(json.dumps(entries[:20], ensure_ascii=False, indent=2), encoding="utf-8")
        temporary.replace(self.path)
