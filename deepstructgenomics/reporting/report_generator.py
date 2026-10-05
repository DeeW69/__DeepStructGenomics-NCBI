"""Serialization helpers for DeepStructGenomics."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict

if TYPE_CHECKING:  # pragma: no cover - only used for typing
    from deepstructgenomics.pipeline import PipelineResult


@dataclass
class ReportPaths:
    """Holds the locations of the exported artifacts."""

    json_path: Path
    markdown_path: Path


def safe_output_stem(identifier: str) -> str:
    """Keep a sequence identifier within one portable filename component."""

    stem = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", identifier).rstrip(" .") or "sequence"
    reserved = {"CON", "PRN", "AUX", "NUL", "CONIN$", "CONOUT$"}
    reserved.update(f"{prefix}{number}" for prefix in ("COM", "LPT") for number in range(1, 10))
    if stem.split(".", 1)[0].upper() in reserved:
        stem = "_" + stem
    return stem


def export_report(result: "PipelineResult", output_dir: Path) -> ReportPaths:
    """Write JSON and Markdown representations to disk."""

    payload = _build_payload(result)
    base_name = safe_output_stem(result.sequence_record.identifier)
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
    if result.impact_summary is not None:
        payload["impact_summary"] = result.impact_summary
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
            lines.append("| Position (1-based) | Reference | Mutant |")
            lines.append("|----------|-----------|--------|")
            for sub in var["substitutions"][:50]:
                lines.append(f"| {sub['position'] + 1} | {sub['reference']} | {sub['mutant']} |")
            if len(var["substitutions"]) > 50:
                lines.append("")
                lines.append(f"_({len(var['substitutions']) - 50} substitutions supplementaires non listees)_")
        lines.append("")

    if "impact_summary" in payload:
        lines.extend(_render_impact_summary(payload["impact_summary"]))

    lines.append("## Metadonnees")
    lines.append("")
    for key, value in annotations.get("metadata", {}).items():
        lines.append(f"- **{key}** : {value}")
    lines.append("")
    lines.append(f"_rapport genere le {payload['generated_at']}_")
    return "\n".join(lines)


def _render_impact_summary(summary: Dict[str, Any]) -> list[str]:
    """Explain positional score changes without treating them as measured effects."""

    lines = ["## Résumé des variations par position", ""]
    note = summary.get("note")
    if note == "mutant not provided":
        return lines + ["Comparaison indisponible : aucune séquence mutante fournie.", ""]
    if note == "delta scores unavailable":
        return lines + ["Comparaison indisponible : aucun score delta disponible.", ""]
    if note:
        lines.extend([
            "Les séquences ont des longueurs différentes : seules les positions communes "
            "sont comparées, sans alignement. Les insertions et délétions ne sont pas réalignées.",
            "",
        ])

    lines.extend([
        "Positions numérotées à partir de 1. Delta = score mutant - score référence.",
        "Ces scores heuristiques décrivent la composition et les appariements prédits, "
        "pas un effet biologique mesuré.",
        "",
    ])
    stats = summary["delta_statistics"]
    lines.extend([
        f"- Delta minimum / maximum : {stats['min']:+.3f} / {stats['max']:+.3f}",
        f"- Delta moyen / médian : {stats['mean']:+.3f} / {stats['median']:+.3f}",
        f"- Positions avec un delta non nul : {stats['nonzero_ratio']:.1%}",
        "",
    ])
    parameters = summary["parameters"]
    threshold = parameters["min_abs_delta"]
    lines.extend([
        f"### Positions les plus modifiées (jusqu'à {parameters['top_k']}, |delta| ≥ {threshold:g})",
        "",
    ])
    if summary["hotspots"]:
        lines.extend(["| Position | Delta | Amplitude |", "|----------|-------|-----------|"])
        for hotspot in summary["hotspots"]:
            lines.append(f"| {hotspot['position']} | {hotspot['delta']:+.3f} | {hotspot['abs_delta']:.3f} |")
    else:
        lines.append(f"Aucune position n'atteint le seuil |delta| ≥ {threshold:g}.")
    lines.append("")

    pairs = summary.get("base_pair_summary")
    if pairs:
        lines.extend([
            f"Paires dont au moins une position atteint |delta| ≥ {pairs['threshold']:g} "
            "(ce comptage ne décrit pas les paires perdues ou gagnées) :",
            "",
        ])
        for label, values in pairs["sets"].items():
            display_name = {"wt": "Référence", "mutant": "Mutant"}.get(label, label)
            lines.append(
                f"- {display_name} : {values['affected_pairs']}/{values['total_pairs']} "
                f"({values['ratio']:.1%})"
            )
        lines.append("")
    return lines
