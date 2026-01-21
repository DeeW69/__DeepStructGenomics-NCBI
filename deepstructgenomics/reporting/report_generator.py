"""Serialization helpers for DeepStructGenomics."""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Optional

if TYPE_CHECKING:  # pragma: no cover - only used for typing
    from deepstructgenomics.pipeline import PipelineResult


@dataclass
class ReportPaths:
    """Holds the locations of the exported artifacts."""

    json_path: Path
    markdown_path: Path


def export_report(result: "PipelineResult", output_dir: Path) -> ReportPaths:
    """Write JSON and Markdown representations to disk."""

    payload = _build_payload(result)
    base_name = result.sequence_record.identifier.replace("|", "_")
    json_path = output_dir / f"{base_name}.json"
    markdown_path = output_dir / f"{base_name}.md"

    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    markdown_path.write_text(_render_markdown(payload), encoding="utf-8")
    return ReportPaths(json_path=json_path, markdown_path=markdown_path)


def _build_payload(result: "PipelineResult") -> Dict[str, Any]:
    """Convert PipelineResult into a serializable dictionary."""

    structure_block = {
        "dot_bracket": result.structure.dot_bracket,
        "base_pairs": list(result.structure.base_pairs),
        "stability_score": result.structure.stability_score,
        "metrics": result.structure.metrics,
    }
    payload: Dict[str, Any] = {
        "identifier": result.sequence_record.identifier,
        "description": result.sequence_record.description,
        "sequence": result.sequence_record.sequence,
        "annotations": result.annotations,
        "structure": structure_block,
        "generated_at": result.generated_at.replace(tzinfo=None).isoformat() + "Z",
    }
    if result.variant_result:
        payload["variant_analysis"] = {
            "total_differences": result.variant_result.total_differences,
            "substitutions": result.variant_result.substitutions,
            "structure_delta_score": result.variant_result.structure_delta_score,
            "commentary": result.variant_result.commentary,
        }
    return payload


def _render_markdown(payload: Dict[str, Any]) -> str:
    """Create a human-readable scientific summary."""

    lines = []
    lines.append(f"# Rapport structurel - {payload['identifier']}")
    lines.append("")
    lines.append(f"**Description** : {payload.get('description', 'N/A')}")
    annotations = payload["annotations"]
    lines.append(
        f"- Longueur : {annotations['length']} nt  \n"
        f"- GC% : {annotations['gc_content']*100:.2f}  \n"
        f"- AU ratio : {annotations['au_ratio']*100:.2f} %"
    )
    lines.append("")
    struct = payload["structure"]
    lines.append("## Structure secondaire (dot-bracket)")
    lines.append("")
    lines.append("```text")
    lines.append(payload["sequence"])
    lines.append(struct["dot_bracket"])
    lines.append("```")
    lines.append("")
    lines.append(
        f"Stabilite heuristique : {struct['stability_score']:.3f} - "
        f"Fraction appariee : {struct['metrics']['paired_fraction']:.2f}"
    )
    lines.append("")

    if "variant_analysis" in payload:
        var = payload["variant_analysis"]
        lines.append("## Impact des variants")
        lines.append("")
        lines.append(
            f"- Nombre de substitutions : {var['total_differences']}\n"
            f"- Score delta structure : {var['structure_delta_score']:.2f}\n"
            f"- Interpretation : {var['commentary']}"
        )
        if var["substitutions"]:
            lines.append("")
            lines.append("| Position | Reference | Mutant |")
            lines.append("|----------|-----------|--------|")
            for sub in var["substitutions"][:50]:
                lines.append(f"| {sub['position']} | {sub['reference']} | {sub['mutant']} |")
            if len(var["substitutions"]) > 50:
                lines.append("")
                lines.append(f"_({len(var['substitutions']) - 50} substitutions supplementaires non listees)_")
        lines.append("")

    lines.append("## Metadonnees")
    lines.append("")
    for key, value in annotations.get("metadata", {}).items():
        lines.append(f"- **{key}** : {value}")
    lines.append("")
    lines.append(f"_rapport genere le {payload['generated_at']}_")
    return "\n".join(lines)
