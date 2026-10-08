"""Immutable alignment settings; environment loading belongs to entry points."""
from dataclasses import asdict, dataclass
import math
import os
from pathlib import Path


@dataclass(frozen=True)
class AlignmentConfig:
    enabled: bool = True
    max_length: int = 10_000
    max_cells: int = 4_000_000
    match_score: float = 2
    mismatch_score: float = -1
    gap_open_score: float = -5
    gap_extend_score: float = -1
    overflow_policy: str = "reject"

    def __post_init__(self):
        if type(self.enabled) is not bool:
            raise ValueError("enabled doit être un booléen.")
        for name in ("max_length", "max_cells"):
            if type(getattr(self, name)) is not int or getattr(self, name) <= 0:
                raise ValueError(f"{name} doit être un entier strictement positif.")
        for name in ("match_score", "mismatch_score", "gap_open_score", "gap_extend_score"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"{name} doit être un nombre fini.")
        if self.overflow_policy not in ("reject", "positional"):
            raise ValueError("overflow_policy : choisir reject ou positional.")

    @property
    def scoring(self):
        return {key: getattr(self, key + "_score") for key in ("match", "mismatch", "gap_open", "gap_extend")}

    def to_dict(self):
        return asdict(self)


def alignment_cell_count(wt_length, mut_length):
    """Theoretical work indicator WT × MUT, not allocated bytes or DP borders."""
    if any(type(n) is not int or n < 0 for n in (wt_length, mut_length)):
        raise ValueError("Les longueurs doivent être des entiers positifs ou nuls.")
    return wt_length * mut_length


def alignment_metadata(config, wt_length, mut_length, status):
    return {"configuration": config.to_dict(), "status": status,
            "matrix": {"wt_length": wt_length, "mut_length": mut_length,
                       "cells": alignment_cell_count(wt_length, mut_length)}}


def overflow_reason(config, wt_length, mut_length):
    cells = alignment_cell_count(wt_length, mut_length)
    if max(wt_length, mut_length) > config.max_length or cells > config.max_cells:
        return (f"Alignement trop coûteux : WT {wt_length:,} nt × MUT {mut_length:,} nt = {cells:,} cellules. "
                f"Limite : {config.max_cells:,} cellules ; longueur maximale : {config.max_length:,} nt. "
                "Modifiez la limite ou choisissez explicitement la comparaison positionnelle.")
    return None


def validate_alignment_metadata(metadata, wt_length, mut_length):
    try:
        config = AlignmentConfig(**metadata["configuration"])
        expected = alignment_metadata(config, wt_length, mut_length, metadata["status"])
        if metadata["matrix"] != expected["matrix"] or metadata["status"] not in ("aligned", "disabled", "reference_only", "overflow_positional"):
            raise ValueError
        return config
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Configuration d'alignement enregistrée incohérente.") from exc


ENV_KEYS = {"enabled": "DSG_ALIGNMENT_ENABLED", "max_length": "DSG_MAX_ALIGNMENT_LENGTH",
            "max_cells": "DSG_MAX_ALIGNMENT_CELLS", "match_score": "DSG_ALIGNMENT_MATCH",
            "mismatch_score": "DSG_ALIGNMENT_MISMATCH", "gap_open_score": "DSG_ALIGNMENT_GAP_OPEN",
            "gap_extend_score": "DSG_ALIGNMENT_GAP_EXTEND", "overflow_policy": "DSG_ALIGNMENT_OVERFLOW_POLICY"}
INPUT_KEYS = {"sequence": "DSG_WT_SEQUENCE", "mutant_sequence": "DSG_MUT_SEQUENCE",
              "fasta_path": "DSG_WT_FASTA", "mutant_fasta_path": "DSG_MUT_FASTA"}


def read_env_file(path=None):
    """Small literal KEY=VALUE format; no shell expansion or execution."""
    target = Path(path) if path is not None else Path.cwd() / ".env"
    if path is None and not target.exists():
        return {}
    values = {}
    for number, line in enumerate(target.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        key, sep, value = line.partition("=")
        key, value = key.strip(), value.strip()
        if not key.startswith("DSG_"):
            continue
        if not sep:
            raise ValueError(f"{target}:{number} : attendu DSG_NOM=valeur.")
        if value.startswith(("'", '"')):
            if len(value) < 2 or value[-1] != value[0]:
                raise ValueError(f"{target}:{number} : guillemets non fermés.")
            value = value[1:-1]
        values[key] = value
    return values


def _convert(name, value):
    try:
        if name == "enabled":
            if value.lower() not in ("true", "false", "1", "0"):
                raise ValueError
            return value.lower() in ("true", "1")
        if name in ("max_length", "max_cells"):
            return int(value)
        if name.endswith("_score"):
            return float(value)
        return value
    except (ValueError, AttributeError):
        raise ValueError(f"{ENV_KEYS[name]} : valeur invalide {value!r}.") from None


def load_alignment_config(env_file=None, *, environ=None, cli=None, gui=None):
    values = {}
    for layer in (read_env_file(env_file), os.environ if environ is None else environ):
        values.update({name: _convert(name, layer[key]) for name, key in ENV_KEYS.items() if key in layer})
    for layer in (cli, gui):
        values.update({key: value for key, value in (layer or {}).items() if value is not None})
    return AlignmentConfig(**values)


def load_sequence_inputs(env_file=None, *, environ=None, cli=None, gui=None):
    values = {}
    for layer in (read_env_file(env_file), os.environ if environ is None else environ):
        _merge_inputs(values, {name: layer[key] for name, key in INPUT_KEYS.items() if layer.get(key)})
    for layer in (cli, gui):
        _merge_inputs(values, {key: value for key, value in (layer or {}).items() if value is not None})
    return values


def _merge_inputs(values, layer):
    for group in (("sequence", "fasta_path", "accession", "batch_fasta"), ("mutant_sequence", "mutant_fasta_path")):
        selected = [key for key in group if key in layer]
        if len(selected) > 1:
            raise ValueError(f"Une seule source par séquence : {', '.join(selected)}.")
        if selected:
            for key in group:
                values.pop(key, None)
            values[selected[0]] = layer[selected[0]]


def add_alignment_arguments(parser):
    parser.add_argument("--env-file", type=Path, help="Fichier .env externe (défaut : répertoire courant).")
    switch = parser.add_mutually_exclusive_group()
    switch.add_argument("--no-alignment", dest="enabled", action="store_false", default=None)
    switch.add_argument("--alignment", dest="enabled", action="store_true")
    for name in ("max_length", "max_cells"):
        parser.add_argument("--" + name.replace("_", "-"), type=int)
    for name in ("match", "mismatch", "gap_open", "gap_extend"):
        parser.add_argument("--" + name.replace("_", "-"), dest=name + "_score", type=float)
    parser.add_argument("--overflow-policy", choices=("reject", "positional"))


def config_from_args(args):
    return load_alignment_config(args.env_file, cli={key: getattr(args, key, None) for key in ENV_KEYS})
