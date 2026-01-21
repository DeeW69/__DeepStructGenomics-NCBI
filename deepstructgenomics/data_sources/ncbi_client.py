"""Lightweight client to interact with the NCBI E-utilities."""

from __future__ import annotations

from dataclasses import dataclass
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
        metadata = self._fetch_summary(accession)
        seq_id, description, sequence = self._parse_fasta(fasta_text)
        return SequenceRecord(
            identifier=seq_id or accession,
            description=description or metadata.get("title", ""),
            sequence=sequence,
            metadata=metadata,
        )

    def _fetch_summary(self, accession: str) -> Dict[str, Any]:
        """Retrieve structured metadata for an accession."""

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
            result = payload.get("result", {})
            uid = next((key for key in result.keys() if key != "uids"), accession)
            summary = result.get(uid, {})
        except ValueError:
            summary = {}
        return summary

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
    def _parse_fasta(fasta_text: requests.Response | str):
        """Parse the FASTA payload and return (identifier, description, sequence)."""

        text = fasta_text.text if isinstance(fasta_text, requests.Response) else fasta_text
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        if not lines:
            raise ValueError("FASTA payload is empty.")
        header = lines[0]
        seq_lines = lines[1:]
        if not header.startswith(">"):
            raise ValueError("Invalid FASTA header")
        header = header[1:]
        parts = header.split(None, 1)
        seq_id = parts[0]
        desc = parts[1] if len(parts) > 1 else ""
        sequence = "".join(seq_lines).replace("U", "T").upper()
        return seq_id, desc, sequence


def normalize_user_sequence(sequence: str, label: str) -> SequenceRecord:
    """Create a SequenceRecord from an arbitrary RNA sequence supplied by the user."""

    cleaned = "".join(base for base in sequence.upper() if base in {"A", "C", "G", "T", "U"})
    if not cleaned:
        raise ValueError("Aucune base valide trouvee dans la sequence fournie.")
    cleaned = cleaned.replace("T", "U")
    return SequenceRecord(
        identifier=label,
        description=f"User provided sequence ({label})",
        sequence=cleaned,
        metadata={"source": "user"},
    )
