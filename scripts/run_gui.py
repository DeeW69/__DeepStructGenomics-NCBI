"""Launch the optional desktop workspace from a source checkout."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from deepstructgenomics.gui.app import main

if __name__ == "__main__":
    raise SystemExit(main())
