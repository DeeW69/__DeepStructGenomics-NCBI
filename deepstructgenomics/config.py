"""Configuration data classes for the DeepStructGenomics pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class NCBIConfig:
    """Holds the default parameters required to query the NCBI E-utilities."""

    base_url: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    database: str = "nuccore"
    tool_name: str = "DeepStructGenomics"
    email: Optional[str] = None
    api_key: Optional[str] = None
    timeout: int = 15  # seconds


@dataclass
class RNAConfig:
    """Parameters used by the RNA secondary structure module."""

    min_loop_length: int = 3
    gc_weight: float = 0.5
    au_weight: float = 0.4
    gu_weight: float = 0.2


@dataclass
class PipelineConfig:
    """Global pipeline settings."""

    ncbi: NCBIConfig = field(default_factory=NCBIConfig)
    rna: RNAConfig = field(default_factory=RNAConfig)
    cache_dir: Path = field(default_factory=lambda: Path("data/cache"))
    default_output_dir: Path = field(default_factory=lambda: Path("outputs"))

