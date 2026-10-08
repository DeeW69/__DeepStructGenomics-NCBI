"""NCBI client and validated RNA sequence inputs."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
import tempfile
import time
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

    def __init__(self, config: Optional[NCBIConfig] = None, cache_dir: Optional[Path] = None) -> None:
        self.config = config or NCBIConfig()
        self.cache_dir = Path(cache_dir) if cache_dir is not None else None
        self._session = requests.Session()
        self._last_call = 0.0

    def search_sequences(self, term: str, *, organism: str = "", rna_only: bool = True,
                         start: int = 0, page_size: int = 20) -> Dict[str, Any]:
        """Search nuccore and retrieve one bounded page of document summaries."""
        term, organism = term.strip(), organism.strip()
        if not term or len(term) > 1000:
            raise ValueError("Saisir une recherche NCBI de 1 à 1000 caractères.")
        if start < 0 or not 1 <= page_size <= 50:
            raise ValueError("Pagination NCBI invalide.")
        query = f"({term})"
        if organism:
            if any(char in organism for char in '"[]'):
                raise ValueError("Saisir un nom d'organisme, sans guillemets ni crochets.")
            query += f' AND "{organism}"[Organism]'
        if rna_only:
            properties = ("mrna", "ncrna", "rrna", "trna", "snrna", "snorna", "scrna", "transcribed_rna")
            query += " AND (" + " OR ".join(f"biomol_{value}[PROP]" for value in properties) + ")"

        def read_json(endpoint, parameters):
            try:
                payload = self._call_endpoint(endpoint, parameters).json()
            except requests.exceptions.JSONDecodeError as exc:
                raise ValueError("Réponse NCBI illisible. Réessayer la recherche.") from exc
            if not isinstance(payload, dict) or payload.get("error"):
                raise ValueError("NCBI n'a pas pu traiter la recherche.")
            return payload

        payload = read_json("esearch.fcgi", {"db": self.config.database, "term": query,
                            "retmode": "json", "retstart": start, "retmax": page_size})
        search = payload.get("esearchresult")
        if not isinstance(search, dict) or search.get("ERROR"):
            raise ValueError("Requête NCBI invalide ou non reconnue.")
        errors = search.get("errorlist") or {}
        if not isinstance(errors, dict) or errors.get("fieldsnotfound"):
            raise ValueError("Champ de recherche NCBI non reconnu.")
        try:
            total = int(search["count"])
            ids = search["idlist"]
            if total < 0 or not isinstance(ids, list) or any(not isinstance(uid, str) or not uid.isdigit() for uid in ids):
                raise ValueError
        except (ValueError, KeyError, TypeError) as exc:
            raise ValueError("Pagination reçue de NCBI invalide.") from exc
        rows = []
        if ids:
            summaries = read_json("esummary.fcgi", {"db": self.config.database, "id": ",".join(ids[:page_size]), "retmode": "json"}).get("result")
            if not isinstance(summaries, dict):
                raise ValueError("Métadonnées NCBI indisponibles.")
            for uid in ids[:page_size]:
                summary = summaries.get(uid)
                if not isinstance(summary, dict) or summary.get("error") or not summary.get("accessionversion"):
                    raise ValueError("Une notice NCBI est incomplète ; relancer la recherche.")
                try:
                    length = int(summary["slen"])
                    if length < 1:
                        raise ValueError
                except (ValueError, KeyError, TypeError) as exc:
                    raise ValueError("Longueur de séquence NCBI invalide.") from exc
                rows.append({"accession": str(summary["accessionversion"]),
                             "organism": str(summary.get("organism", "Non renseigné")),
                             "molecule": str(summary.get("biomol") or summary.get("moltype") or "Non renseigné"),
                             "length": length, "description": str(summary.get("title", ""))})
        return {"term": term, "query": query, "translated_query": search.get("querytranslation", query),
                "total": total, "start": start, "page_size": page_size, "rows": rows,
                "unmatched_terms": errors.get("phrasesnotfound", [])}

    def fetch_sequence(self, accession: str, *, refresh: bool = False) -> SequenceRecord:
        """Fetches a nucleotide sequence (FASTA + metadata) from NCBI."""

        key = [self.config.base_url, self.config.database, accession]
        cache_path = None
        if self.cache_dir is not None:
            digest = hashlib.sha256(json.dumps(key).encode()).hexdigest()
            cache_path = self.cache_dir / f"{digest}.json"
            if cache_path.exists() and not refresh:
                try:
                    payload = json.loads(cache_path.read_text(encoding="utf-8"))
                    record = SequenceRecord(**payload["record"])
                    checksum = hashlib.sha256(record.sequence.encode()).hexdigest()
                    if (payload["schema"] != 1 or payload["key"] != key
                            or checksum != record.metadata["provenance"]["sha256"]
                            or normalize_rna_sequence(record.sequence) != record.sequence):
                        raise ValueError("Cache incoherent")
                    record.metadata["provenance"]["cache_hit"] = True
                    return record
                except (ValueError, KeyError, TypeError, AttributeError):
                    pass  # A corrupt entry is replaced only after a successful fetch.

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
        record = SequenceRecord(
            identifier=seq_id or accession,
            description=description or metadata.get("title", ""),
            sequence=sequence,
            metadata=metadata,
        )
        if cache_path is not None:
            record.metadata["provenance"] = {
                "source": "NCBI", "base_url": self.config.base_url,
                "database": self.config.database, "requested_accession": accession,
                "resolved_identifier": record.identifier,
                "retrieved_at": datetime.now(timezone.utc).isoformat(),
                "sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "cache_hit": False,
            }
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            fd, temporary = tempfile.mkstemp(dir=cache_path.parent, suffix=".tmp")
            try:
                with os.fdopen(fd, "w", encoding="utf-8") as stream:
                    json.dump({"schema": 1, "key": key, "record": asdict(record)}, stream)
                os.replace(temporary, cache_path)
            finally:
                if os.path.exists(temporary):
                    os.unlink(temporary)
        return record

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
        # E-utilities permit 3 requests/s without a key (10 with a key).
        delay = (0.11 if self.config.api_key else 0.35) - (time.monotonic() - self._last_call)
        if delay > 0:
            time.sleep(delay)
        self._last_call = time.monotonic()
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
