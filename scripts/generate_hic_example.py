#!/usr/bin/env python3
"""Generate a seeded synthetic Hi-C region with two domains and one added loop."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def generate(output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(46)
    n, bin_size = 32, 10000
    bias = np.exp(rng.normal(0, 0.25, n))
    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            same_domain = (i < 16) == (j < 16)
            mean = (100 if same_domain else 3) / np.sqrt(j - i) * bias[i] * bias[j]
            if (i, j) == (5, 10):
                mean *= 15
            rows.append(["chrSynthetic", i * bin_size, "chrSynthetic", j * bin_size, int(rng.poisson(mean))])
    with (output_dir / "hic_region.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom1", "start1", "chrom2", "start2", "count"])
        writer.writerows(rows)
    truth = {"source": "Synthetic Poisson generator; no experimental Hi-C data",
             "seed": 46, "bins": n, "bin_size": bin_size, "boundary_bin": 16,
             "added_loop": [5, 10], "within_domain_amplitude": 100, "cross_domain_amplitude": 3,
             "distance_decay_exponent": 0.5, "bias_log_sd": 0.25, "loop_multiplier": 15,
             "warning": "Domain structure violates a homogeneous distance-only loop null; extra loop calls are not biological discoveries."}
    (output_dir / "hic_region_truth.json").write_text(json.dumps(truth, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("data/examples"))
    generate(parser.parse_args().output_dir)
