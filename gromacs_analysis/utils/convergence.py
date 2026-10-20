"""Sampling-sufficiency checks for time series and principal-component projections.

- Block-averaged standard error (Flyvbjerg & Petersen, J. Chem. Phys. 1989, 91, 461):
  the series is halved repeatedly by averaging neighbours; the plateau (largest) SEM
  over levels with at least ``min_blocks`` blocks is reported, and the effective number
  of independent samples is ``N_eff = var / SEM^2``.
- Half agreement: the production window's first and second halves should agree within
  their uncertainties; ``z = |m1 - m2| / sqrt(SEM1^2 + SEM2^2)``. A large z means drift.
- Cosine content (Hess, Phys. Rev. E 2002, 65, 031910): for a projection ``p_i`` onto
  PC i, ``c_i = (2/T) (sum_t cos(i*pi*(t+1/2)/T) p_i(t))^2 / sum_t p_i(t)^2``. Values
  near 1 mean the PC looks like random diffusion (the trajectory has not converged).
- Half-split RMSIP: overlap of the essential subspaces of the first and second half of
  the frames; low values mean the dominant motions are not yet sampled consistently.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np

MIN_POINTS = 8

DEFAULT_SERIES_THRESHOLDS = {"neff_min": 20.0, "half_z_max": 3.0}
DEFAULT_PCA_THRESHOLDS = {"cosine_max": 0.5, "half_rmsip_min": 0.5}


def block_average(values, min_blocks: int = 4) -> Dict[str, float]:
    """Plateau block-averaged SEM. Empty dict when the series is too short."""
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    n = len(data)
    if n < 2 * min_blocks:
        return {}
    variance = float(np.var(data, ddof=1))
    naive_sem = np.sqrt(variance / n) if variance > 0 else 0.0
    best = naive_sem
    blocked = data
    while len(blocked) >= 2 * min_blocks:
        if len(blocked) % 2:
            blocked = blocked[:-1]
        blocked = 0.5 * (blocked[0::2] + blocked[1::2])
        best = max(best, float(np.sqrt(np.var(blocked, ddof=1) / len(blocked))))
    if best <= 0 or variance <= 0:
        return {"sem": 0.0, "n_effective": float(n), "statistical_inefficiency": 1.0}
    n_eff = variance / best ** 2
    return {
        "sem": float(best),
        "n_effective": float(min(n, n_eff)),
        "statistical_inefficiency": float(max(1.0, n / n_eff)),
    }


def _sem(values: np.ndarray) -> float:
    block = block_average(values)
    if block:
        return block["sem"]
    if len(values) < 2:
        return float("nan")
    return float(np.std(values, ddof=1) / np.sqrt(len(values)))


def half_agreement(values) -> Dict[str, float]:
    """Difference of the two halves' means and its z-score using each half's SEM."""
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    half = len(data) // 2
    if half < 2:
        return {"half_diff": float("nan"), "half_z": float("nan")}
    first, second = data[:half], data[half:2 * half]
    diff = float(second.mean() - first.mean())
    scale = float(np.hypot(_sem(first), _sem(second)))
    if not np.isfinite(scale):
        z = float("nan")
    elif scale == 0:
        z = 0.0 if diff == 0 else float("inf")
    else:
        z = abs(diff) / scale
    return {"half_diff": diff, "half_z": z}


def series_convergence(values, thresholds: Optional[Dict[str, float]] = None) -> Dict:
    """N, mean, block SEM, N_eff, half agreement and flags for one replicate's series."""
    limits = {**DEFAULT_SERIES_THRESHOLDS, **(thresholds or {})}
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    row: Dict = {"N": int(len(data)), "Mean": float(data.mean()) if len(data) else float("nan")}
    if len(data) < MIN_POINTS:
        row.update({"SEM_block": float("nan"), "N_eff": float("nan"), "Half_Diff": float("nan"),
                    "Half_Z": float("nan"), "Flags": "too_short"})
        return row
    block = block_average(data)
    halves = half_agreement(data)
    row.update({
        "SEM_block": block.get("sem", float("nan")),
        "N_eff": block.get("n_effective", float(len(data))),
        "Half_Diff": halves["half_diff"],
        "Half_Z": halves["half_z"],
    })
    flags: List[str] = []
    if np.isfinite(row["N_eff"]) and row["N_eff"] < limits["neff_min"]:
        flags.append("low_neff")
    if np.isfinite(row["Half_Z"]) and row["Half_Z"] > limits["half_z_max"]:
        flags.append("drift")
    row["Flags"] = ";".join(flags)
    return row


def cosine_content(projection, pc_index: int) -> float:
    """Hess cosine content of a 1-D projection for PC number ``pc_index`` (1-based)."""
    p = np.asarray(projection, dtype=float)
    p = p[np.isfinite(p)]
    t_total = len(p)
    if t_total < 2:
        return float("nan")
    p = p - p.mean()
    denom = float(np.sum(p ** 2))
    if denom == 0:
        return 0.0
    t = np.arange(t_total)
    num = float(np.sum(np.cos(pc_index * np.pi * (t + 0.5) / t_total) * p)) ** 2
    return float(np.clip(2.0 / t_total * num / denom, 0.0, 1.0))


def _top_vectors(flat: np.ndarray, k: int) -> np.ndarray:
    centred = flat - flat.mean(axis=0)
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    return vt[:k].T


def half_split_rmsip(flat, k: int = 10) -> float:
    """RMSIP between the top-k PCs of the first and second half of ``flat`` (frames, dims)."""
    data = np.asarray(flat, dtype=float)
    half = len(data) // 2
    k = int(min(k, data.shape[1], half - 1)) if half > 1 else 0
    if k < 1:
        return float("nan")
    a = _top_vectors(data[:half], k)
    b = _top_vectors(data[half:2 * half], k)
    overlap = a.T @ b
    return float(np.sqrt(np.sum(overlap ** 2) / k))


def projection_convergence(projection, n_pcs: int = 3) -> Dict[str, float]:
    """Cosine content of the first ``n_pcs`` columns of a (frames, k) projection."""
    proj = np.asarray(projection, dtype=float)
    if proj.ndim == 1:
        proj = proj[:, None]
    return {f"Cos_PC{i + 1}": (cosine_content(proj[:, i], i + 1) if i < proj.shape[1] else float("nan"))
            for i in range(n_pcs)}
