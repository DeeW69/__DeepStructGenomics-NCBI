import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from deepstructgenomics.analysis.hic_advanced import (
    analyze_region, balance_contacts, benjamini_hochberg, call_domains, call_loops, reconstruct_3d,
)

ROOT = Path(__file__).resolve().parents[1]


def test_balancing_recovers_known_bias_and_masks_empty_bins():
    bias = np.array([1., 2., 3., 4.])
    raw = (np.ones((4, 4)) - np.eye(4)) * bias[:, None] * bias[None, :]
    raw = np.pad(raw, ((0, 1), (0, 1)))
    balanced, info = balance_contacts(raw)
    assert info["converged"]
    assert info["masked_bins"] == [4]
    np.testing.assert_allclose(balanced[:4, :4].sum(axis=1), 1, rtol=1e-6)
    np.testing.assert_allclose(balanced[:4, :4], (np.ones((4, 4)) - np.eye(4)) / 3, atol=1e-6)
    assert np.all(balanced[4] == 0)


def test_nonconvergent_support_is_reported():
    raw = np.array([[0., 1, 1], [1, 0, 0], [1, 0, 0]])
    result, info = balance_contacts(raw, max_iterations=20)
    assert not info["converged"]
    assert np.isfinite(result).all()


@pytest.mark.parametrize("raw", [np.eye(3), [[0, -1], [-1, 0]], [[0, 2], [1, 0]], [[0, np.nan], [np.nan, 0]]])
def test_bad_matrix_rejected(raw):
    with pytest.raises(ValueError):
        balance_contacts(raw)


def test_bh_known_probabilities_and_empty():
    assert benjamini_hochberg([]) == []
    np.testing.assert_allclose(benjamini_hochberg([.01, .04, .03, .2]), [.04, .0533333333333333, .0533333333333333, .2])


def test_loop_enrichment_and_fractional_counts_rejected():
    pytest.importorskip("scipy")
    raw = np.full((15, 15), 10.)
    np.fill_diagonal(raw, 0)
    raw[3, 6] = raw[6, 3] = 1000
    result = call_loops(raw, np.ones(15))
    significant = [(r["bin1"], r["bin2"]) for r in result["tests"] if r["significant"]]
    assert significant == [(3, 6)]
    with pytest.raises(ValueError, match="entiers"):
        call_loops(raw + .5, np.ones(15))


def test_seeded_domain_test_finds_inserted_boundary():
    raw = np.ones((24, 24))
    raw[:12, :12] = raw[12:, 12:] = 50
    np.fill_diagonal(raw, 0)
    balanced, _ = balance_contacts(raw)
    a = call_domains(balanced, permutations=499, seed=46, window=3)
    b = call_domains(balanced, permutations=499, seed=46, window=3)
    assert a == b
    assert 12 in a["selected_boundaries"]
    assert a["candidate_domains"][0]["end_bin"] == 12


def test_mds_recovers_euclidean_geometry_up_to_scale_and_rotation():
    pytest.importorskip("scipy")
    points = np.array([[0., 0, 0], [1, 0, 0], [0, 2, 0], [0, 0, 3]])
    distance = np.linalg.norm(points[:, None] - points[None, :], axis=2)
    raw = np.zeros_like(distance)
    np.divide(1, distance ** 3, out=raw, where=distance > 0)
    result = reconstruct_3d(raw)
    assert result["status"] == "ok"
    coords = np.array([[p["x"], p["y"], p["z"]] for p in result["coordinates"]])
    fitted = np.linalg.norm(coords[:, None] - coords[None, :], axis=2)
    expected = distance / np.median(distance[np.triu_indices(4, 1)])
    np.testing.assert_allclose(fitted, expected, atol=1e-8)
    assert result["stress"] < 1e-8


def test_disconnected_graph_not_fabricated():
    pytest.importorskip("scipy")
    matrix = np.zeros((4, 4))
    matrix[0, 1] = matrix[1, 0] = matrix[2, 3] = matrix[3, 2] = 1
    assert reconstruct_3d(matrix)["status"] == "disconnected"


def test_explicit_region_and_missingness_required():
    with pytest.raises(ValueError, match="missing-as-zero"):
        analyze_region("unused", assembly="test", chromosome="chr1", start=0, end=100, bin_size=10)
    with pytest.raises(ValueError, match="400"):
        analyze_region("unused", assembly="test", chromosome="chr1", start=0, end=5000, bin_size=10, missing_as_zero=True)


def test_advanced_cli_exports_and_refuses_reuse(tmp_path):
    pytest.importorskip("scipy")
    command = [sys.executable, str(ROOT / "scripts/explore_hic.py"), "--contacts", str(ROOT / "data/examples/hic_region.tsv"),
               "--assembly", "synthetic", "--bin-size", "10000", "--advanced", "--chromosome", "chrSynthetic",
               "--start", "0", "--end", "320000", "--missing-as-zero", "--permutations", "99", "--no-plot",
               "--output-dir", str(tmp_path / "analysis")]
    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    report = json.loads((tmp_path / "analysis/hic_analysis.json").read_text())
    assert report["status"] == "ok"
    assert report["reconstruction"]["status"] == "ok"
    assert (tmp_path / "analysis/coordinates.tsv").exists()
    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert result.returncode == 2


def test_overdispersion_avoids_calling_domain_contacts_as_loops():
    pytest.importorskip("scipy")
    report = analyze_region(ROOT / "data/examples/hic_region.tsv", assembly="synthetic",
                            chromosome="chrSynthetic", start=0, end=320000, bin_size=10000,
                            missing_as_zero=True, permutations=99)
    selected = [(r["bin1"], r["bin2"]) for r in report["loops"]["tests"] if r["significant"]]
    assert selected == [(5, 10)]
    assert any(r["dispersion"] > 0 for r in report["loops"]["tests"])


def test_nonconvergence_prevents_statistics(tmp_path):
    pytest.importorskip("scipy")
    path = tmp_path / "star.tsv"
    path.write_text("chrom1\tstart1\tchrom2\tstart2\tcount\nchr1\t0\tchr1\t10\t1\nchr1\t0\tchr1\t20\t1\n")
    result = analyze_region(path, assembly="synthetic", chromosome="chr1", start=0, end=30,
                            bin_size=10, missing_as_zero=True, window=1, max_iterations=10)
    assert result["status"] == "not_converged"
    assert result["loops"] is None and result["reconstruction"] is None
