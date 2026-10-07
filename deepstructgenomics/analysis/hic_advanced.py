"""Small cis-region analysis: balancing, exploratory tests and distance-based MDS.

This is a transparent research baseline, not an implementation of a published
loop caller. Statistical assumptions and numerical diagnostics travel with results.
"""

import csv
from datetime import datetime, timezone
import hashlib
import io
from pathlib import Path

import numpy as np

from deepstructgenomics.analysis.hic import summarize_contacts


def _matrix(value):
    matrix = np.asarray(value, dtype=float)
    if (matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or len(matrix) < 2
            or not np.isfinite(matrix).all() or np.any(matrix < 0)
            or not np.allclose(matrix, matrix.T, rtol=1e-12, atol=1e-12)):
        raise ValueError("Matrice symetrique carree, finie, non negative requise.")
    return matrix


def balance_contacts(matrix, *, tolerance=1e-6, max_iterations=2000):
    """Symmetric multiplicative equal-visibility balancing with damped updates.

    Diagonal contacts are excluded; zero-coverage bins are masked. Scaling vectors
    are multiplicative: balanced[i,j] = input[i,j] * weight[i] * weight[j].
    Active row sums target one. Non-convergence is a result, never hidden.
    """
    matrix = _matrix(matrix).copy()
    if not np.isfinite(tolerance) or tolerance <= 0 or max_iterations < 1:
        raise ValueError("Tolerance positive et au moins une iteration requises.")
    np.fill_diagonal(matrix, 0)
    active = matrix.sum(axis=1) > 0
    if np.count_nonzero(active) < 2:
        raise ValueError("Au moins deux bins avec contacts hors diagonale requis.")
    work = matrix[np.ix_(active, active)]
    weights = np.ones(len(work))
    converged, error = False, float("inf")
    for iteration in range(max_iterations + 1):
        balanced = work * weights[:, None] * weights[None, :]
        rows = balanced.sum(axis=1)
        if not np.isfinite(rows).all() or np.any(rows <= 0):
            raise ValueError("Equilibrage numeriquement instable.")
        mean = float(rows.mean())
        error = float(np.max(np.abs(rows / mean - 1)))
        if error <= tolerance:
            converged = True
            break
        if iteration < max_iterations:
            log_weights = np.log(weights) - 0.5 * np.log(rows / mean)
            log_weights -= log_weights.mean()
            if np.max(np.abs(log_weights)) > 200:
                break  # Infeasible support can drive weights to infinity.
            weights = np.exp(log_weights)
    weights /= np.sqrt(mean)
    full_weights = np.zeros(len(matrix))
    full_weights[active] = weights
    result = matrix * full_weights[:, None] * full_weights[None, :]
    return result, {"method": "symmetric damped equal-visibility balancing",
                    "converged": converged, "iterations": iteration,
                    "max_relative_row_error": error, "tolerance": tolerance,
                    "diagonal": "excluded", "target_active_row_sum": 1.0,
                    "masked_bins": np.flatnonzero(~active).tolist(),
                    "weights": full_weights.tolist()}


def benjamini_hochberg(pvalues):
    values = np.asarray(pvalues, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all() or np.any((values < 0) | (values > 1)):
        raise ValueError("Probabilites entre 0 et 1 requises.")
    if not len(values):
        return []
    order = np.argsort(values, kind="stable")
    adjusted = values[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum(1, np.minimum.accumulate(adjusted[::-1])[::-1])
    result = np.empty(len(values))
    result[order] = adjusted
    return result.tolist()


def call_loops(raw, weights, *, fdr=0.05, min_separation=2, min_background_pixels=5,
               min_background_count=20):
    """One-sided negative-binomial enrichment conditional on distance/bias.

    For each pixel, fit its distance-stratum rate on all OTHER usable pixels.
    Zero pixels are included. Excess variance of exposure-corrected background
    rates estimates gamma-Poisson dispersion by moments (Poisson if zero).
    Rate/dispersion/bias estimation uncertainty is not modeled.
    """
    try:
        from scipy.stats import nbinom, poisson
    except ImportError as exc:
        raise ValueError("SciPy absent : installer requirements-research-lock.txt.") from exc
    raw = _matrix(raw)
    weights = np.asarray(weights, dtype=float)
    if (weights.shape != (len(raw),) or not np.isfinite(weights).all() or np.any(weights < 0)
            or not np.isfinite(fdr) or not 0 < fdr < 1 or min_separation < 1
            or min_background_pixels < 2 or min_background_count < 1):
        raise ValueError("Parametres de tests de boucles invalides.")
    if np.any(raw != np.floor(raw)):
        raise ValueError("Les tests de comptage necessitent des comptes bruts entiers.")
    active = weights > 0
    candidates, skipped = [], 0
    for distance in range(min_separation, len(raw)):
        left = np.arange(len(raw) - distance)
        right = left + distance
        keep = active[left] & active[right]
        left, right = left[keep], right[keep]
        if len(left) - 1 < min_background_pixels:
            skipped += len(left)
            continue
        observed = raw[left, right]
        exposure = 1 / (weights[left] * weights[right])
        for offset, (i, j, count) in enumerate(zip(left, right, observed)):
            other_count = float(observed.sum() - count)
            other_exposure = float(exposure.sum() - exposure[offset])
            if other_count < min_background_count or other_exposure <= 0:
                skipped += 1
                continue
            expected = other_count / other_exposure * float(exposure[offset])
            others = np.arange(len(observed)) != offset
            rate = other_count / other_exposure
            rates = observed[others] / exposure[others]
            sampling_noise = rate * float(np.mean(1 / exposure[others]))
            dispersion = max(0.0, (float(np.var(rates, ddof=1)) - sampling_noise) / rate ** 2)
            if dispersion > 1e-8:
                shape = 1 / dispersion
                pvalue = float(nbinom.sf(count - 1, shape, shape / (shape + expected)))
            else:
                pvalue = float(poisson.sf(count - 1, expected))
            if not np.isfinite(expected) or not np.isfinite(pvalue):
                raise ValueError("Modele de comptage numeriquement instable.")
            candidates.append({"bin1": int(i), "bin2": int(j), "count": int(count),
                               "expected": expected, "dispersion": dispersion,
                               "enrichment": float(count / expected), "p_value": pvalue})
    for row, qvalue in zip(candidates, benjamini_hochberg([r["p_value"] for r in candidates])):
        row.update(q_value=qvalue, significant=qvalue <= fdr and row["enrichment"] > 1)
    return {"model": "Negative binomial (Poisson at zero dispersion); raw counts; distance and balancing bias; leave-one-pixel-out moments",
            "multiple_testing": "Benjamini-Hochberg over all tested pixels in this region",
            "fdr_threshold": fdr, "tested_pixels": len(candidates), "untested_pixels": skipped,
            "min_separation_bins": min_separation, "min_background_pixels": min_background_pixels,
            "min_background_count": min_background_count, "tests": candidates,
            "limitations": "Exploratory q-values: estimated mean/dispersion/bias, few background pixels, dependent tests and no biological replicates; no calibrated FDR guarantee."}


def call_domains(balanced, *, window=3, permutations=999, seed=46, fdr=0.05):
    """Insulation minima with a distance-stratified permutation null and BH.

    Each full cross-boundary square is compared to randomized contact matrices.
    Values, including zeros, shuffle only within the same genomic separation.
    This preserves distance distributions, but not margins or spatial dependence.
    """
    balanced = _matrix(balanced)
    if window < 1 or permutations < 1 or not np.isfinite(fdr) or not 0 < fdr < 1:
        raise ValueError("Fenetre/permutations positives et FDR entre 0 et 1 requis.")
    active = balanced.sum(axis=1) > 0
    cuts = [i for i in range(window, len(balanced) - window + 1)
            if active[i - window:i + window].all()]
    def scores(matrix):
        return np.array([matrix[i - window:i, i:i + window].mean() for i in cuts])
    observed = scores(balanced)
    exceedances = np.zeros(len(cuts), dtype=int)
    rng = np.random.default_rng(seed)
    strata = []
    for distance in range(1, len(balanced)):
        left = np.arange(len(balanced) - distance)
        right = left + distance
        keep = active[left] & active[right]
        strata.append((left[keep], right[keep], balanced[left[keep], right[keep]]))
    for _ in range(permutations if cuts else 0):
        shuffled = np.zeros_like(balanced)
        for left, right, values in strata:
            values = rng.permutation(values)
            shuffled[left, right] = values
            shuffled[right, left] = values
        exceedances += scores(shuffled) <= observed
    pvalues = (exceedances + 1) / (permutations + 1)
    qvalues = benjamini_hochberg(pvalues)
    positive = observed[observed > 0]
    median = float(np.median(positive)) if len(positive) else 0
    tests, selected = [], []
    for index, cut in enumerate(cuts):
        value = float(observed[index])
        neighbors = [observed[j] for j in range(max(0, index - 1), min(len(cuts), index + 2))
                     if j != index and abs(cuts[j] - cut) == 1]
        minimum = bool(neighbors) and all(value <= neighbor for neighbor in neighbors)
        row = {"boundary_bin": cut, "insulation": value,
               "log2_insulation_over_median": float(np.log2(value / median)) if value > 0 and median > 0 else None,
               "p_value": float(pvalues[index]), "q_value": qvalues[index], "local_minimum": minimum,
               "selected": False}
        tests.append(row)
    # Keep the strongest minimum in a window; deterministic ties use genomic order.
    for row in sorted(tests, key=lambda r: (r["q_value"], r["insulation"], r["boundary_bin"])):
        if row["local_minimum"] and row["q_value"] <= fdr and all(abs(row["boundary_bin"] - b) >= window for b in selected):
            row["selected"] = True
            selected.append(row["boundary_bin"])
    edges = [0, *sorted(selected), len(balanced)]
    return {"method": "square insulation; within-distance permutation lower-tail test",
            "window_bins": window, "permutations": permutations, "seed": seed, "fdr_threshold": fdr,
            "multiple_testing": "Benjamini-Hochberg over all full-window boundaries in this region",
            "tests": tests, "selected_boundaries": sorted(selected),
            "candidate_domains": [{"start_bin": a, "end_bin": b, "active_fraction": float(active[a:b].mean())}
                                  for a, b in zip(edges, edges[1:])],
            "limitations": "Candidate intervals, not validated TADs; exchangeability within distance ignores spatial dependence and does not preserve balanced margins."}


def reconstruct_3d(balanced, *, exponent=1 / 3):
    """Contact-derived shortest-path distances followed by classical MDS in 3D."""
    try:
        from scipy.sparse.csgraph import shortest_path
    except ImportError as exc:
        raise ValueError("SciPy absent : installer requirements-research-lock.txt.") from exc
    balanced = _matrix(balanced)
    if not np.isfinite(exponent) or exponent <= 0:
        raise ValueError("Exposant de distance positif requis.")
    active = np.flatnonzero(balanced.sum(axis=1) > 0)
    matrix = balanced[np.ix_(active, active)]
    np.fill_diagonal(matrix, 0)
    if len(active) < 2:
        return {"status": "insufficient_contacts", "coordinates": []}
    graph = np.full_like(matrix, np.inf)
    positive = matrix > 0
    scale = float(np.median(matrix[positive]))
    with np.errstate(over="ignore", invalid="ignore"):
        graph[positive] = (matrix[positive] / scale) ** (-exponent)
    if not np.isfinite(graph[positive]).all():
        raise ValueError("Exposant incompatible avec les poids : distances non finies.")
    np.fill_diagonal(graph, 0)
    distances = shortest_path(graph, directed=False)
    if not np.isfinite(distances).all():
        return {"status": "disconnected", "coordinates": [], "active_bins": active.tolist(),
                "reason": "No artificial edges inserted between disconnected components."}
    # Normalize distance scale before squaring to keep eigendecomposition stable.
    distances /= np.median(distances[np.triu_indices(len(distances), 1)])
    centering = np.eye(len(active)) - np.ones((len(active), len(active))) / len(active)
    gram = -0.5 * centering @ (distances ** 2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    indices = np.argsort(eigenvalues)[::-1][:3]
    coords = np.zeros((len(active), 3))
    for axis, index in enumerate(indices):
        if eigenvalues[index] > 0:
            coords[:, axis] = eigenvectors[:, index] * np.sqrt(eigenvalues[index])
            if coords[np.argmax(np.abs(coords[:, axis])), axis] < 0:
                coords[:, axis] *= -1
    fitted = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
    stress = float(np.sqrt(np.sum((fitted - distances) ** 2) / np.sum(distances ** 2)))
    absolute = float(np.abs(eigenvalues).sum())
    return {"status": "ok", "method": "inverse-contact shortest paths + classical MDS",
            "distance_exponent": exponent, "units": "arbitrary; median target distance = 1",
            "stress": stress, "negative_eigenvalue_fraction": float(-eigenvalues[eigenvalues < 0].sum() / absolute) if absolute else 0,
            "coordinates": [{"bin": int(index), "x": float(x), "y": float(y), "z": float(z)}
                            for index, (x, y, z) in zip(active, coords)],
            "limitations": "Inferred population consensus under a contact-distance assumption; no physical distance calibration or experimental 3D validation. Orientation is arbitrary."}


def analyze_region(path, *, assembly, chromosome, start, end, bin_size,
                   missing_as_zero=False, window=3, permutations=999, seed=46, fdr=0.05,
                   min_separation=2, exponent=1 / 3, tolerance=1e-6, max_iterations=2000):
    """Analyze an explicitly bounded cis region of at most 400 bins."""
    if (not missing_as_zero or bin_size <= 0 or start < 0 or end <= start
            or start % bin_size or end % bin_size or not chromosome):
        raise ValueError("Region alignee et --missing-as-zero explicite requis ; les contacts absents representent des zeros observes.")
    n = (end - start) // bin_size
    if not 2 <= n <= 400:
        raise ValueError("L'analyse dense accepte de 2 a 400 bins ; reduire la region ou la resolution.")
    if window < 1 or 2 * window > n or not 1 <= permutations <= 10000:
        raise ValueError("Fenetre complete requise et permutations entre 1 et 10000.")
    try:
        import scipy
    except ImportError as exc:
        raise ValueError("SciPy absent : installer requirements-research-lock.txt.") from exc
    summary = summarize_contacts(path, assembly=assembly, bin_size=bin_size, top_k=0)
    raw_bytes = Path(path).read_bytes()
    if hashlib.sha256(raw_bytes).hexdigest() != summary["sha256"]:
        raise ValueError("Fichier modifie pendant la lecture.")
    matrix = np.zeros((n, n), dtype=float)
    included = 0
    for row in csv.DictReader(io.StringIO(raw_bytes.decode("utf-8-sig")), delimiter="\t"):
        a, b = int(row["start1"]), int(row["start2"])
        if row["chrom1"].strip() == row["chrom2"].strip() == chromosome and start <= a < end and start <= b < end:
            count = float(row["count"])
            if count != int(count) or count > 1e12:
                raise ValueError("Comptes bruts entiers <= 10^12 requis pour les statistiques.")
            i, j = (a - start) // bin_size, (b - start) // bin_size
            matrix[i, j] = matrix[j, i] = count
            included += 1
    balanced, normalization = balance_contacts(matrix, tolerance=tolerance, max_iterations=max_iterations)
    report = {"schema": 1, "generated_at": datetime.now(timezone.utc).isoformat(),
              "source": str(path), "sha256": summary["sha256"], "assembly": assembly,
              "region": {"chromosome": chromosome, "start": start, "end": end, "bin_size": bin_size, "bins": n},
              "included_rows": included, "excluded_rows": summary["contact_rows"] - included,
              "missing_contacts": "zero by explicit user declaration; zero-coverage bins masked",
              "normalization": normalization, "balanced_matrix": balanced.tolist()}
    if not normalization["converged"]:
        report.update(status="not_converged", loops=None, domains=None, reconstruction=None)
        return report
    report["software"] = {"numpy": np.__version__, "scipy": scipy.__version__}
    report.update(status="ok",
                  loops=call_loops(matrix, normalization["weights"], fdr=fdr, min_separation=min_separation),
                  domains=call_domains(balanced, window=window, permutations=permutations, seed=seed, fdr=fdr),
                  reconstruction=reconstruct_3d(balanced, exponent=exponent))
    return report
