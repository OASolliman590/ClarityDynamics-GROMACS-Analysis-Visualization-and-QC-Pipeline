"""Pairwise ligand comparisons on replicate-level values.

Replicates are the independent samples (Constitution II), so every test uses one value
per replicate: the replicate mean of an MD metric, or the replicate mean of MM/PBSA
ΔTOTAL. For each metric and each pair of ligands this reports:

- Welch's t-test (unequal variances);
- a permutation test on the difference of means, exact when the number of
  relabellings is small (otherwise Monte-Carlo with a fixed seed). With 3 vs 3
  replicates the smallest attainable two-sided p is 0.1 - noted in the output;
- Hedges' g (bias-corrected standardised difference);
- Holm and Benjamini-Hochberg adjusted p-values over all pairs of that metric.
"""

from __future__ import annotations

import itertools
import math
from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

EXACT_LIMIT = 20000
COLUMNS = [
    "Metric", "LigandA", "LigandB", "n_A", "n_B", "Mean_A", "Mean_B", "Diff", "Welch_t", "Welch_p",
    "Perm_p", "Hedges_g", "p_holm", "p_bh", "Direction", "Note",
]


def welch(a: Sequence[float], b: Sequence[float]) -> Dict[str, float]:
    from scipy import stats  # noqa: WPS433 - scipy is a core dependency

    a, b = np.asarray(a, float), np.asarray(b, float)
    if np.var(a, ddof=1) == 0 and np.var(b, ddof=1) == 0:
        same = math.isclose(a.mean(), b.mean())
        return {"t": 0.0 if same else math.copysign(math.inf, a.mean() - b.mean()), "p": 1.0 if same else 0.0}
    result = stats.ttest_ind(a, b, equal_var=False)
    return {"t": float(result.statistic), "p": float(result.pvalue)}


def permutation_p(a: Sequence[float], b: Sequence[float], n_permutations: int = 10000, seed: int = 0) -> float:
    """Two-sided p for |mean(a) - mean(b)| under random relabelling (observed counted)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    pooled = np.concatenate([a, b])
    n_a, n = len(a), len(pooled)
    observed = abs(a.mean() - b.mean())
    tol = 1e-12 * max(1.0, abs(observed))
    total = math.comb(n, n_a)
    if total <= EXACT_LIMIT:
        hits = 0
        for idx in itertools.combinations(range(n), n_a):
            mask = np.zeros(n, bool)
            mask[list(idx)] = True
            if abs(pooled[mask].mean() - pooled[~mask].mean()) >= observed - tol:
                hits += 1
        return hits / total
    rng = np.random.default_rng(seed)
    perms = np.argsort(rng.random((int(n_permutations), n)), axis=1)
    shuffled = pooled[perms]
    diffs = np.abs(shuffled[:, :n_a].mean(axis=1) - shuffled[:, n_a:].mean(axis=1))
    return float((np.sum(diffs >= observed - tol) + 1) / (len(diffs) + 1))


def hedges_g(a: Sequence[float], b: Sequence[float]) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    n_a, n_b = len(a), len(b)
    pooled_var = ((n_a - 1) * np.var(a, ddof=1) + (n_b - 1) * np.var(b, ddof=1)) / (n_a + n_b - 2)
    diff = a.mean() - b.mean()
    if pooled_var <= 0:
        return 0.0 if diff == 0 else math.copysign(math.inf, diff)
    correction = 1.0 - 3.0 / (4.0 * (n_a + n_b) - 9.0)
    return float(diff / math.sqrt(pooled_var) * correction)


def holm(pvalues: Sequence[float]) -> List[float]:
    p = np.asarray(pvalues, float)
    out = np.full(len(p), np.nan)
    valid = np.where(np.isfinite(p))[0]
    m = len(valid)
    if m == 0:
        return out.tolist()
    order = valid[np.argsort(p[valid])]
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, min(1.0, (m - rank) * p[idx]))
        out[idx] = running
    return out.tolist()


def benjamini_hochberg(pvalues: Sequence[float]) -> List[float]:
    p = np.asarray(pvalues, float)
    out = np.full(len(p), np.nan)
    valid = np.where(np.isfinite(p))[0]
    m = len(valid)
    if m == 0:
        return out.tolist()
    order = valid[np.argsort(p[valid])]
    running = 1.0
    for rank in range(m - 1, -1, -1):
        idx = order[rank]
        running = min(running, p[idx] * m / (rank + 1))
        out[idx] = min(1.0, running)
    return out.tolist()


def pairwise_tests(
    values: Mapping[str, Mapping[str, Sequence[float]]],
    directions: Optional[Mapping[str, str]] = None,
    n_permutations: int = 10000,
    seed: int = 0,
) -> pd.DataFrame:
    """``values[metric][ligand] -> replicate values``; one row per metric and ligand pair."""
    directions = directions or {}
    rows: List[Dict] = []
    for metric in sorted(values):
        groups = {lig: [float(v) for v in vals if v is not None and np.isfinite(v)]
                  for lig, vals in values[metric].items()}
        metric_rows = []
        for lig_a, lig_b in itertools.combinations(sorted(groups), 2):
            a, b = groups[lig_a], groups[lig_b]
            row = {
                "Metric": metric, "LigandA": lig_a, "LigandB": lig_b, "n_A": len(a), "n_B": len(b),
                "Mean_A": float(np.mean(a)) if a else np.nan, "Mean_B": float(np.mean(b)) if b else np.nan,
                "Direction": directions.get(metric, "two-sided"),
            }
            row["Diff"] = row["Mean_A"] - row["Mean_B"]
            if len(a) < 2 or len(b) < 2:
                row.update({"Welch_t": np.nan, "Welch_p": np.nan, "Perm_p": np.nan, "Hedges_g": np.nan,
                            "Note": "insufficient replicates (need >= 2 per ligand)"})
            else:
                w = welch(a, b)
                perm = permutation_p(a, b, n_permutations=n_permutations, seed=seed)
                min_p = 2.0 / math.comb(len(a) + len(b), len(a)) if math.comb(len(a) + len(b), len(a)) <= EXACT_LIMIT else None
                note = f"smallest attainable permutation p = {min_p:.3g}" if min_p and min_p > 0.01 else ""
                row.update({"Welch_t": w["t"], "Welch_p": w["p"], "Perm_p": perm,
                            "Hedges_g": hedges_g(a, b), "Note": note})
            metric_rows.append(row)
        adjusted_holm = holm([r["Welch_p"] for r in metric_rows])
        adjusted_bh = benjamini_hochberg([r["Welch_p"] for r in metric_rows])
        for row, ph, pb in zip(metric_rows, adjusted_holm, adjusted_bh):
            row["p_holm"], row["p_bh"] = ph, pb
        rows.extend(metric_rows)
    return pd.DataFrame(rows, columns=COLUMNS)
