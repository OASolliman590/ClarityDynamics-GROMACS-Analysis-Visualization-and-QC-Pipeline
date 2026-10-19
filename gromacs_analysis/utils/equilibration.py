"""Automated equilibration detection (Chodera, J. Chem. Theory Comput. 2016, 12, 1799).

The production window starts at the index t0 that maximises the effective number of
uncorrelated samples ``N_eff(t0) = (N - t0) / g(t0)``, where ``g`` is the statistical
inefficiency of the series from t0 onward. Discarding too little keeps biased initial
data (inflating variance); discarding too much throws away good samples, so N_eff has
a maximum. Unlike threshold rules this needs no unit-dependent parameters.

This is a dependency-free implementation of the same estimator used by
``pymbar.timeseries``; when ``use_pymbar=True`` and pymbar is installed, pymbar is used.
"""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np

MIN_POINTS = 10


def _autocorrelation(dx: np.ndarray) -> np.ndarray:
    """Normalised autocorrelation of a centred series via FFT (C[0] = 1)."""
    n = len(dx)
    size = 1 << (2 * n - 1).bit_length()
    spectrum = np.fft.rfft(dx, size)
    acov = np.fft.irfft(spectrum * np.conj(spectrum), size)[:n]
    acov /= np.arange(n, 0, -1)          # unbiased: divide by (N - t)
    return acov / acov[0]


def statistical_inefficiency(x, mintime: int = 3) -> float:
    """g = 1 + 2 Σ_t (1 - t/N) C(t), summed until C(t) first drops to <= 0 (after ``mintime``)."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 2:
        return 1.0
    dx = x - x.mean()
    if np.allclose(dx, 0):
        return 1.0
    acf = _autocorrelation(dx)
    g = 1.0
    for t in range(1, n - 1):
        c = acf[t]
        if c <= 0 and t > mintime:
            break
        g += 2.0 * c * (1.0 - t / n)
    return max(1.0, float(g))


def detect_equilibration(
    x,
    n_candidates: int = 100,
    max_fraction: float = 0.5,
    use_pymbar: bool = False,
) -> Tuple[int, float, float]:
    """Return ``(t0, g, n_eff)`` for the production start maximising N_eff.

    Candidate starts are spread evenly over the first ``max_fraction`` of the series.
    Series shorter than 10 points return ``(0, 1.0, nan)`` ("insufficient data").
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < MIN_POINTS or not np.any(np.isfinite(x)):
        return 0, 1.0, float("nan")

    if use_pymbar:
        try:
            from pymbar import timeseries  # noqa: WPS433
            t0, g, n_eff = timeseries.detect_equilibration(x)
            return int(t0), float(g), float(n_eff)
        except ImportError:
            pass

    last = max(0, min(n - MIN_POINTS, int(n * max_fraction)))
    candidates = np.unique(np.linspace(0, last, max(2, int(n_candidates))).astype(int))
    best = (0, 1.0, -math.inf)
    for t0 in candidates:
        g = statistical_inefficiency(x[t0:])
        n_eff = (n - t0) / g
        if n_eff > best[2]:
            best = (int(t0), g, n_eff)
    return best[0], float(best[1]), float(best[2])
