"""Desktop operations delegated to the existing scientific pipeline (no Qt)."""
from datetime import datetime, timezone
import json
from pathlib import Path
from uuid import uuid4

from deepstructgenomics.config import NCBIConfig, PipelineConfig
from deepstructgenomics.pipeline import DeepStructPipeline, PipelineInput
from deepstructgenomics.visualization.secondary_view import load_secondary_data


def new_run_directory(workspace, kind):
    root = Path(workspace).expanduser().resolve()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
    directory = root / f"{kind}_{stamp}_{uuid4().hex[:8]}"
    directory.mkdir(parents=True, exist_ok=False)
    return directory


def run_rna(workspace, parameters):
    parameters = dict(parameters)
    email = parameters.pop("ncbi_email", None)
    directory = new_run_directory(workspace, "rna")
    pipeline = DeepStructPipeline(PipelineConfig(
        ncbi=NCBIConfig(email=email), cache_dir=Path(workspace).resolve() / "cache",
        default_output_dir=directory,
    ))
    result = pipeline.run_and_export(PipelineInput(**parameters), directory)
    return {"kind": "rna", "path": str(result.visualization_paths.manifest_path),
            "directory": str(directory), "label": result.sequence_record.identifier}


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
    delta = [b - a for a, b in zip(wt["scores"], mut["scores"])] if mut else []
    metrics = {
        "wt_pairs": len(wt_pairs), "mut_pairs": len(mut_pairs) if mut else None,
        "lost": len(wt_pairs - mut_pairs) if mut else None,
        "gained": len(mut_pairs - wt_pairs) if mut else None,
        "changed": sum(abs(value) > 1e-9 for value in delta) if mut else None,
        "mean_delta": sum(delta) / len(delta) if delta else None,
        "max_abs_delta": max(map(abs, delta)) if delta else None,
    }
    return data, metrics


def position_details(data, position):
    """Use biological positions at the UI boundary, including unequal lengths."""
    result = {"position": position}
    for key in ("reference", "mutant"):
        structure = data.get(key)
        if not structure or not 1 <= position <= len(structure["sequence"]):
            result[key] = None
            continue
        partner = next((j + 1 if i == position - 1 else i + 1
                        for i, j in structure["base_pairs"] if position - 1 in (i, j)), None)
        result[key] = {"base": structure["sequence"][position - 1],
                       "score": structure["scores"][position - 1], "partner": partner}
    wt, mut = result["reference"], result["mutant"]
    result["delta"] = mut["score"] - wt["score"] if wt and mut else None
    return result


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
