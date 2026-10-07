import json
import subprocess
import sys
from pathlib import Path

import pytest

from deepstructgenomics.batch import run_fasta_batch


def test_batch_duplicates_invalid_and_provenance(tmp_path):
    fasta = tmp_path / "batch.fa"
    fasta.write_text(">../same description\nACGT\n>bad\nACNT\n>../same\nGGGGAAAACCCC\n")
    output = tmp_path / "results"
    summary = run_fasta_batch(fasta, output)
    assert (summary["succeeded"], summary["failed"]) == (2, 1)
    for item in (summary["entries"][0], summary["entries"][2]):
        report = json.loads((output / item["report"]).read_text())
        assert report["annotations"]["metadata"]["entry_index"] == item["index"]
    assert (output / "entry_000002/error.json").exists()
    with pytest.raises(FileExistsError):
        run_fasta_batch(fasta, output)


@pytest.mark.parametrize("text", ["", "ACGU\n>x\nACGU"])
def test_bad_framing(tmp_path, text):
    fasta = tmp_path / "batch.fa"
    fasta.write_text(text)
    with pytest.raises(ValueError):
        run_fasta_batch(fasta, tmp_path / "results")


def test_batch_cli_partial_failure(tmp_path):
    fasta = tmp_path / "batch.fa"
    fasta.write_text(">good\nACGU\n>bad\nN\n")
    script = Path(__file__).resolve().parents[1] / "scripts/run_pipeline.py"
    result = subprocess.run([sys.executable, str(script), "--batch-fasta", str(fasta),
                             "--output-dir", str(tmp_path / "results")], capture_output=True)
    assert result.returncode == 2
    assert (tmp_path / "results/batch_summary.json").exists()
