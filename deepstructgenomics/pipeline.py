"""High-level orchestration for the DeepStructGenomics pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

from deepstructgenomics.config import PipelineConfig
from deepstructgenomics.data_sources.ncbi_client import (
    NCBIClient,
    SequenceRecord,
    normalize_user_sequence,
)
from deepstructgenomics.reporting.report_generator import ReportPaths, export_report
from deepstructgenomics.rna.secondary_structure import (
    SecondaryStructureResult,
    predict_secondary_structure,
)
from deepstructgenomics.variants.variant_analysis import VariantImpactResult, compare_sequences


@dataclass
class PipelineInput:
    """User-facing description of what should be processed."""

    accession: Optional[str] = None
    sequence: Optional[str] = None
    sequence_label: str = "custom_sequence"
    mutant_sequence: Optional[str] = None

    def __post_init__(self) -> None:
        if not self.accession and not self.sequence:
            raise ValueError("Un identifiant NCBI ou une sequence doivent etre fournis.")


@dataclass
class PipelineResult:
    """Structured outputs for downstream consumers."""

    sequence_record: SequenceRecord
    annotations: Dict[str, Any]
    structure: SecondaryStructureResult
    variant_result: Optional[VariantImpactResult]
    generated_at: datetime
    report_paths: Optional[ReportPaths] = None


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
        annotations = self._annotate_sequence(sequence_record)
        structure = predict_secondary_structure(sequence_record.sequence, self.config.rna)

        variant_result: Optional[VariantImpactResult] = None
        if request.mutant_sequence:
            mutant = normalize_user_sequence(request.mutant_sequence, f"{sequence_record.identifier}_mut")
            mutant_structure = predict_secondary_structure(mutant.sequence, self.config.rna)
            variant_result = compare_sequences(
                sequence_record.sequence,
                mutant.sequence,
                structure.base_pairs,
                mutant_structure.base_pairs,
            )

        return PipelineResult(
            sequence_record=sequence_record,
            annotations=annotations,
            structure=structure,
            variant_result=variant_result,
            generated_at=datetime.utcnow(),
        )

    def run_and_export(self, request: PipelineInput, output_dir: Optional[Path] = None) -> PipelineResult:
        """Runs the pipeline then writes JSON/Markdown reports."""

        result = self.run(request)
        target_dir = Path(output_dir) if output_dir else self.config.default_output_dir
        target_dir.mkdir(parents=True, exist_ok=True)
        report_paths = export_report(result, target_dir)
        result.report_paths = report_paths
        return result

    def _resolve_sequence(self, request: PipelineInput) -> SequenceRecord:
        if request.accession:
            return self.ncbi_client.fetch_sequence(request.accession)
        assert request.sequence is not None
        return normalize_user_sequence(request.sequence, request.sequence_label)

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
