"""Sequential FASTA batches with isolated outputs and explicit per-entry failures."""

from pathlib import Path
import json

from deepstructgenomics.data_sources.ncbi_client import parse_single_fasta
from deepstructgenomics.pipeline import DeepStructPipeline, PipelineInput
from deepstructgenomics.reporting.report_generator import export_report


def run_fasta_batch(path: str | Path, output_dir: str | Path, *, pipeline=None,
                    top_k=10, min_abs_delta=0.1, base_pair_threshold=0.2):
    """Write one report bundle per record and a summary; continue invalid records.

    Input before the first header and empty files are fatal framing errors.
    Record numbering, rather than user identifiers, determines output folders.
    """
    path, output_dir = Path(path), Path(output_dir)
    pipeline = pipeline or DeepStructPipeline()
    entries = []
    current = []
    with path.open(encoding="utf-8-sig") as stream:
        for line in stream:
            if not line.strip():
                continue
            if line.lstrip().startswith(">"):
                if current:
                    entries.append("".join(current))
                current = [line]
            elif not current:
                raise ValueError("FASTA : donnees avant le premier en-tete.")
            else:
                current.append(line)
    if current:
        entries.append("".join(current))
    if not entries:
        raise ValueError("Fichier FASTA vide.")
    # Refuse reuse: otherwise failures could leave apparently valid old reports.
    output_dir.mkdir(parents=True, exist_ok=False)
    summary = {"source": str(path), "entries": [], "succeeded": 0, "failed": 0}
    for number, text in enumerate(entries, 1):
        folder = output_dir / f"entry_{number:06d}"
        folder.mkdir()
        item = {"index": number, "header": text.splitlines()[0][1:].strip(),
                "directory": folder.name}
        try:
            identifier, description, sequence = parse_single_fasta(text)
            request = PipelineInput(sequence=sequence, sequence_label=identifier,
                                    top_k=top_k, min_abs_delta=min_abs_delta,
                                    base_pair_threshold=base_pair_threshold)
            result = pipeline.run(request)
            result.sequence_record.description = description
            result.sequence_record.metadata = {
                "source": "fasta", "path": str(path), "fasta_identifier": identifier,
                "entry_index": number,
            }
            result.annotations["metadata"] = result.sequence_record.metadata
            result.report_paths = export_report(result, folder)
            result.visualization_paths = pipeline._export_visualization_bundle(result, folder)
            item.update(status="ok", identifier=identifier,
                        report=str(result.report_paths.json_path.relative_to(output_dir)))
            summary["succeeded"] += 1
        except ValueError as exc:
            item.update(status="error", error=str(exc))
            (folder / "error.json").write_text(json.dumps(item, indent=2), encoding="utf-8")
            summary["failed"] += 1
        summary["entries"].append(item)
    (output_dir / "batch_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary
