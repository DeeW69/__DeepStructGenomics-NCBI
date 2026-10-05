"""High-level orchestration for the DeepStructGenomics pipeline."""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

from deepstructgenomics.analysis.impact_summary import build_impact_summary
from deepstructgenomics.config import PipelineConfig
from deepstructgenomics.data_sources.ncbi_client import (
    NCBIClient,
    SequenceRecord,
    load_fasta_record,
    normalize_user_sequence,
)
from deepstructgenomics.reporting.report_generator import ReportPaths, export_report, safe_output_stem
from deepstructgenomics.rna.secondary_structure import (
    SecondaryStructureResult,
    predict_secondary_structure,
)
from deepstructgenomics.variants.variant_analysis import VariantImpactResult, compare_sequences
from deepstructgenomics.visualization.io_structures import (
    export_sequence_as_pseudo_pdb,
    generate_coarse_backbone,
    write_visualization_manifest,
)
from deepstructgenomics.visualization.overlay_helpers import compute_displacements, summarize_displacements
from deepstructgenomics.visualization.score_mapping import (
    ScoreTable,
    compute_delta_scores,
    compute_score_statistics,
    derive_position_scores_from_structure,
    write_score_table,
)


@dataclass
class PipelineInput:
    """User-facing description of what should be processed."""

    accession: Optional[str] = None
    sequence: Optional[str] = None
    sequence_label: Optional[str] = None
    mutant_sequence: Optional[str] = None
    fasta_path: Optional[str | Path] = None
    mutant_fasta_path: Optional[str | Path] = None

    def __post_init__(self) -> None:
        sources = (self.accession, self.sequence, self.fasta_path)
        if sum(value is not None for value in sources) != 1:
            raise ValueError("Fournir exactement une source : accession NCBI, sequence ou fichier FASTA.")
        if self.accession is not None:
            self.accession = self.accession.strip()
            if not self.accession:
                raise ValueError("L'identifiant NCBI ne peut pas etre vide.")
        if self.mutant_sequence is not None and self.mutant_fasta_path is not None:
            raise ValueError("Fournir soit une sequence mutante, soit un fichier FASTA mutant, pas les deux.")
        if self.sequence_label is not None:
            self.sequence_label = self.sequence_label.strip()
            if not self.sequence_label:
                raise ValueError("Le label ne peut pas etre vide.")


@dataclass
class PipelineResult:
    """Structured outputs for downstream consumers."""

    sequence_record: SequenceRecord
    annotations: Dict[str, Any]
    structure: SecondaryStructureResult
    variant_result: Optional[VariantImpactResult]
    mutant_sequence: Optional[str]
    mutant_structure: Optional[SecondaryStructureResult]
    generated_at: datetime
    report_paths: Optional[ReportPaths] = None
    visualization_paths: Optional["VisualizationArtifactPaths"] = None
    impact_summary: Optional[Dict[str, Any]] = None


@dataclass
class VisualizationArtifactPaths:
    """Holds locations of generated visualization files."""

    root_dir: Path
    wt_structure: Optional[Path]
    mutant_structure: Optional[Path]
    wt_score_file: Optional[Path]
    mutant_score_file: Optional[Path]
    manifest_path: Optional[Path]
    impact_summary_file: Optional[Path] = None

    @property
    def score_file(self) -> Optional[Path]:  # backward-compatible alias
        return self.mutant_score_file


class DeepStructPipeline:
    """Coordinates data retrieval, structure prediction and reporting."""

    def __init__(self, config: Optional[PipelineConfig] = None) -> None:
        self.config = config or PipelineConfig()
        self.ncbi_client = NCBIClient(self.config.ncbi)
        self.config.cache_dir.mkdir(parents=True, exist_ok=True)
        self.config.default_output_dir.mkdir(parents=True, exist_ok=True)

    def run(self, request: PipelineInput) -> PipelineResult:
        """Processes a sequence and returns the core artifacts."""

        sequence_record = self._resolve_sequence(request)
        mutant_record: Optional[SequenceRecord] = None
        if request.mutant_fasta_path is not None:
            mutant_record = load_fasta_record(request.mutant_fasta_path)
        elif request.mutant_sequence is not None:
            mutant_record = normalize_user_sequence(request.mutant_sequence, f"{sequence_record.identifier}_mut")

        # Resolve and validate both inputs before starting the folding algorithm.
        annotations = self._annotate_sequence(sequence_record)
        structure = predict_secondary_structure(sequence_record.sequence, self.config.rna)

        variant_result: Optional[VariantImpactResult] = None
        mutant_sequence: Optional[str] = None
        mutant_structure: Optional[SecondaryStructureResult] = None
        if mutant_record is not None:
            mutant_structure = predict_secondary_structure(mutant_record.sequence, self.config.rna)
            variant_result = compare_sequences(
                sequence_record.sequence,
                mutant_record.sequence,
                structure.base_pairs,
                mutant_structure.base_pairs,
            )
            mutant_sequence = mutant_record.sequence

        result = PipelineResult(
            sequence_record=sequence_record,
            annotations=annotations,
            structure=structure,
            variant_result=variant_result,
            mutant_sequence=mutant_sequence,
            mutant_structure=mutant_structure,
            generated_at=datetime.now(timezone.utc),
        )
        result.impact_summary = self._summarize_impact(result)
        return result

    def run_and_export(self, request: PipelineInput, output_dir: Optional[Path] = None) -> PipelineResult:
        """Runs the pipeline then writes JSON/Markdown reports."""

        result = self.run(request)
        target_dir = Path(output_dir) if output_dir else self.config.default_output_dir
        target_dir.mkdir(parents=True, exist_ok=True)
        report_paths = export_report(result, target_dir)
        result.report_paths = report_paths
        visualization_paths = self._export_visualization_bundle(result, target_dir)
        result.visualization_paths = visualization_paths
        return result

    def _resolve_sequence(self, request: PipelineInput) -> SequenceRecord:
        if request.accession is not None:
            return self.ncbi_client.fetch_sequence(request.accession)
        if request.fasta_path is not None:
            return load_fasta_record(request.fasta_path, label=request.sequence_label)
        assert request.sequence is not None
        return normalize_user_sequence(request.sequence, request.sequence_label or "custom_sequence")

    @staticmethod
    def _summarize_impact(result: PipelineResult) -> Dict[str, Any]:
        """Compute one summary shared by the reports and visualization bundle."""

        delta_scores: Dict[int, float] = {}
        base_pair_sets: Dict[str, Sequence[Tuple[int, int]]] = {
            "wt": [(i + 1, j + 1) for i, j in result.structure.base_pairs],
        }
        note = "mutant not provided"
        if result.mutant_structure is not None:
            wt_scores = ScoreTable(position_scores=derive_position_scores_from_structure(result.structure))
            mutant_scores = ScoreTable(position_scores=derive_position_scores_from_structure(result.mutant_structure))
            delta_scores = compute_delta_scores(wt_scores, mutant_scores).position_scores
            base_pair_sets["mutant"] = [(i + 1, j + 1) for i, j in result.mutant_structure.base_pairs]
            note = None if delta_scores else "delta scores unavailable"
            if len(result.structure.sequence) != len(result.mutant_structure.sequence):
                note = "different sequence lengths; positional comparison without alignment"

        return build_impact_summary(
            result.sequence_record.identifier,
            result.generated_at,
            delta_scores,
            top_k=10,
            min_abs_delta=0.1,
            note=note,
            base_pairs=base_pair_sets,
            base_pair_threshold=0.2,
        )

    @staticmethod
    def _annotate_sequence(record: SequenceRecord) -> Dict[str, Any]:
        """Compute sequence-level descriptors (length, GC%, etc.)."""

        seq = record.sequence.upper()
        length = len(seq)
        gc_count = sum(1 for base in seq if base in {"G", "C"})
        au_count = sum(1 for base in seq if base in {"A", "U"})
        annotations = {
            "length": length,
            "gc_content": gc_count / length if length else 0.0,
            "au_ratio": au_count / length if length else 0.0,
            "metadata": record.metadata,
        }
        return annotations

    def _export_visualization_bundle(
        self,
        result: PipelineResult,
        output_dir: Path,
    ) -> VisualizationArtifactPaths:
        """Generate pseudo-3D structures + score payload for the viewer."""

        run_id = safe_output_stem(result.sequence_record.identifier)
        viz_dir = output_dir / run_id / "visualization"
        viz_dir.mkdir(parents=True, exist_ok=True)

        wt_structure_path = export_sequence_as_pseudo_pdb(
            result.structure.sequence,
            viz_dir / "wt_structure.pdb",
        )

        wt_position_scores = derive_position_scores_from_structure(result.structure)
        wt_residue_scores = {f"A:{idx}:": value for idx, value in wt_position_scores.items()}
        wt_table = ScoreTable(
            position_scores=wt_position_scores,
            residue_scores=wt_residue_scores,
            metadata={"context": "reference", "reference_identifier": result.sequence_record.identifier},
        )
        wt_score_path = write_score_table(viz_dir / "wt_scores.json", wt_table)
        wt_stats = compute_score_statistics(wt_table)

        mutant_structure_path: Optional[Path] = None
        mutant_score_path: Optional[Path] = None
        mutant_stats: Optional[Dict[str, float]] = None
        delta_stats: Optional[Dict[str, float]] = None
        delta_range_mode: Optional[str] = None
        displacement_stats: Optional[Dict[str, float]] = None
        if result.mutant_structure and result.mutant_sequence:
            mutant_structure_path = export_sequence_as_pseudo_pdb(
                result.mutant_structure.sequence,
                viz_dir / "mutant_structure.pdb",
            )
            mut_position_scores = derive_position_scores_from_structure(result.mutant_structure)
            mut_residue_scores = {f"A:{idx}:": value for idx, value in mut_position_scores.items()}
            mutant_table = ScoreTable(
                position_scores=mut_position_scores,
                residue_scores=mut_residue_scores,
                metadata={"context": "mutant", "reference_identifier": result.sequence_record.identifier},
            )
            mutant_score_path = write_score_table(viz_dir / "mutant_scores.json", mutant_table)
            mutant_stats = compute_score_statistics(mutant_table)

            delta_table = compute_delta_scores(wt_table, mutant_table)
            delta_stats = compute_score_statistics(delta_table)
            delta_range_mode = "dynamic_symmetric"

            wt_coords = generate_coarse_backbone(result.structure.sequence)
            mut_coords = generate_coarse_backbone(result.mutant_structure.sequence)
            max_len = min(len(wt_coords), len(mut_coords))
            if max_len:
                mapping = {idx + 1: (idx, idx) for idx in range(max_len)}
                distances = compute_displacements(wt_coords, mut_coords, mapping)
                displacement_stats = summarize_displacements(distances)
        else:
            # maintain backward compatibility by keeping a mutant_scores.json identical to WT
            mutant_score_path = write_score_table(viz_dir / "mutant_scores.json", wt_table)
        impact_summary_path = viz_dir / "impact_summary.json"
        impact_summary_path.write_text(
            json.dumps(result.impact_summary, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

        manifest_path = write_visualization_manifest(
            viz_dir / "visualization_manifest.json",
            identifier=result.sequence_record.identifier,
            source=result.sequence_record.metadata.get("source", "NCBI"),
            parameters={
                "generated_at": result.generated_at.replace(tzinfo=None).isoformat() + "Z",
                "score_context": "reference",
                "note": "Projection 3D contrainte derivee de la structure secondaire.",
                "score_statistics_wt": wt_stats,
                "score_statistics_mutant": mutant_stats,
                "score_statistics_delta": delta_stats,
                "delta_range_mode": delta_range_mode,
                "displacement_statistics": displacement_stats,
            },
            wt_structure=wt_structure_path,
            mutant_structure=mutant_structure_path,
            wt_score_file=wt_score_path,
            mutant_score_file=mutant_score_path,
            impact_summary_file=impact_summary_path,
        )

        return VisualizationArtifactPaths(
            root_dir=viz_dir,
            wt_structure=wt_structure_path,
            mutant_structure=mutant_structure_path,
            wt_score_file=wt_score_path,
            mutant_score_file=mutant_score_path,
            manifest_path=manifest_path,
            impact_summary_file=impact_summary_path,
        )
