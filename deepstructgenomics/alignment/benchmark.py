"""Synthetic fixtures and manual diagnostics; never invokes RNA folding."""
import os
import random
import sys
from time import perf_counter

from .config import alignment_cell_count, overflow_reason
from .pairwise import align_sequences


def generate_test_pair(length, substitutions=1, insertions=0, deletions=0, seed=42, pattern="deterministic-random"):
    counts = (length, substitutions, insertions, deletions)
    if any(type(n) is not int or n < 0 for n in counts) or length == 0 or deletions >= length or substitutions + deletions > length:
        raise ValueError("Tailles invalides : WT non vide, substitutions + délétions ≤ longueur, au moins une base conservée.")
    rng = random.Random(seed)
    if pattern not in ("deterministic-random", "repetitive"):
        raise ValueError("Motif synthétique inconnu.")
    wt = "".join(rng.choice("ACGU") for _ in range(length)) if pattern == "deterministic-random" else ("ACGU" * ((length + 3) // 4))[:length]
    positions = rng.sample(range(length), substitutions + deletions)
    changed, removed = set(positions[:substitutions]), set(positions[substitutions:])
    mut = [rng.choice("".join(b for b in "ACGU" if b != base)) if i in changed else base
           for i, base in enumerate(wt) if i not in removed]
    for _ in range(insertions):
        mut.insert(rng.randrange(len(mut) + 1), rng.choice("ACGU"))
    return wt, "".join(mut)


def peak_process_memory_mb():
    """Peak process RSS/working set, INCLUDING interpreter/imports (not a delta)."""
    try:
        if os.name == "nt":
            import ctypes
            from ctypes import wintypes
            class Counters(ctypes.Structure):
                _fields_ = [("cb", wintypes.DWORD), ("faults", wintypes.DWORD)] + [
                    (name, ctypes.c_size_t) for name in ("peak", "working", "peak_paged", "paged", "peak_nonpaged", "nonpaged", "pagefile", "peak_pagefile")]
            counters = Counters()
            counters.cb = ctypes.sizeof(counters)
            kernel = ctypes.WinDLL("kernel32", use_last_error=True)
            kernel.GetCurrentProcess.restype = wintypes.HANDLE
            psapi = ctypes.WinDLL("psapi", use_last_error=True)
            psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
            if not psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
                return None
            return counters.peak / 1024**2
        import resource
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return peak / (1024**2 if sys.platform == "darwin" else 1024)
    except (ImportError, OSError, AttributeError):
        return None


def measure_alignment(wt, mut, config):
    summary = {"wt_length": len(wt), "mut_length": len(mut), "cells": alignment_cell_count(len(wt), len(mut)),
               "max_cells": config.max_cells, "configuration": config.to_dict()}
    start = perf_counter()
    try:
        result = align_sequences(wt, mut, config)
    except ValueError as exc:
        if not overflow_reason(config, len(wt), len(mut)):
            raise
        result = None
        summary.update(status="rejected", reason=str(exc))
    else:
        summary["status"] = "aligned" if result else "disabled" if not config.enabled else "overflow_positional"
    summary.update(time_seconds=perf_counter() - start, peak_process_mb=peak_process_memory_mb())
    if result:
        summary.update(score=result.score, identity=result.identity, ambiguous=result.ambiguous, **result.counts)
    return summary, result
