"""NCBI client and validated RNA sequence inputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

import requests

from deepstructgenomics.config import NCBIConfig


@dataclass
class SequenceRecord:
    """Represents a nucleotide sequence fetched from NCBI or supplied by the user."""

    identifier: str
    description: str
    sequence: str
    metadata: Dict[str, Any]


class NCBIClient:
    """Minimal wrapper around the NCBI REST endpoints."""

    def __init__(self, config: Optional[NCBIConfig] = None) -> None:
        self.config = config or NCBIConfig()
        self._session = requests.Session()

    def fetch_sequence(self, accession: str) -> SequenceRecord:
        """Fetches a nucleotide sequence (FASTA + metadata) from NCBI."""

        fasta_text = self._call_endpoint(
            "efetch.fcgi",
            {
                "db": self.config.database,
                "id": accession,
                "rettype": "fasta",
                "retmode": "text",
            },
        )
        seq_id, description, sequence = self._parse_fasta(fasta_text)
        metadata = self._fetch_summary(accession)
        return SequenceRecord(
            identifier=seq_id or accession,
            description=description or metadata.get("title", ""),
            sequence=sequence,
            metadata=metadata,
        )

    def _fetch_summary(self, accession: str) -> Dict[str, Any]:
        """Retrieve metadata, tolerating empty or malformed successful responses.

        Transport and HTTP failures still propagate to the caller.
        """

        raw_json = self._call_endpoint(
            "esummary.fcgi",
            {
                "db": self.config.database,
                "id": accession,
                "retmode": "json",
            },
        )
        try:
            payload = raw_json.json()
        except ValueError:
            return {}
        if not isinstance(payload, dict):
            return {}
        result = payload.get("result", {})
        if not isinstance(result, dict):
            return {}
        uid = next((key for key in result if key != "uids"), accession)
        summary = result.get(uid, {})
        return summary if isinstance(summary, dict) else {}

    def _call_endpoint(self, endpoint: str, params: Dict[str, Any]):
        """Generic helper that attaches shared parameters to every request."""

        query = {
            **params,
            "tool": self.config.tool_name,
        }
        if self.config.email:
            query["email"] = self.config.email
        if self.config.api_key:
            query["api_key"] = self.config.api_key

        url = f"{self.config.base_url}/{endpoint}"
        response = self._session.get(url, params=query, timeout=self.config.timeout)
        response.raise_for_status()
        return response

    @staticmethod
    def _parse_fasta(fasta_text: requests.Response | str) -> tuple[str, str, str]:
        """Parse one NCBI record with the same validation as local inputs."""

        text = fasta_text.text if isinstance(fasta_text, requests.Response) else fasta_text
        return parse_single_fasta(text)


def normalize_rna_sequence(sequence: str) -> str:
    """Normalize DNA/RNA without silently removing unsupported bases."""

    compact = "".join(base for base in sequence if not base.isspace())
    if not compact:
        raise ValueError("Sequence vide : fournir au moins une base A, C, G, T ou U.")
    for position, base in enumerate(compact, start=1):
        if base not in "ACGTUacgtu":
            raise ValueError(
                f"Caractere non pris en charge {base!r} a la position {position} "
                "(numerotation depuis 1, sans espaces). Bases acceptees : A, C, G, T, U."
            )
    return compact.upper().replace("T", "U")


def parse_single_fasta(text: str) -> tuple[str, str, str]:
    """Read exactly one nonempty FASTA record and return normalized RNA."""

    lines = [line.strip() for line in text.lstrip("\ufeff").splitlines() if line.strip()]
    if not lines:
        raise ValueError("Fichier FASTA vide.")
    if not lines[0].startswith(">"):
        raise ValueError("En-tete FASTA manquant : la premiere ligne doit commencer par >.")
    header = lines[0][1:].strip()
    if not header:
        raise ValueError("En-tete FASTA vide : un identifiant est requis apres >.")
    if any(line.startswith(">") for line in lines[1:]):
        raise ValueError("Plusieurs sequences FASTA detectees : fournir une seule sequence par fichier.")
    parts = header.split(None, 1)
    identifier = parts[0]
    description = parts[1] if len(parts) > 1 else ""
    return identifier, description, normalize_rna_sequence("".join(lines[1:]))


def load_fasta_record(path: str | Path, label: Optional[str] = None) -> SequenceRecord:
    """Load one UTF-8 FASTA file while preserving its identifier and provenance."""

    file_path = Path(path)
    text = file_path.read_text(encoding="utf-8-sig")
    try:
        identifier, description, sequence = parse_single_fasta(text)
    except ValueError as exc:
        raise ValueError(f"{file_path}: {exc}") from exc
    return SequenceRecord(
        identifier=label if label is not None else identifier,
        description=description,
        sequence=sequence,
        metadata={"source": "fasta", "path": str(file_path), "fasta_identifier": identifier},
    )


def normalize_user_sequence(sequence: str, label: str) -> SequenceRecord:
    """Create a SequenceRecord from an arbitrary RNA sequence supplied by the user."""

    return SequenceRecord(
        identifier=label,
        description=f"User provided sequence ({label})",
        sequence=normalize_rna_sequence(sequence),
        metadata={"source": "user"},
    )
