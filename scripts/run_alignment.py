#!/usr/bin/env python3
"""Align user sequences or a synthetic test pair, without RNA folding."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from deepstructgenomics.alignment.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
