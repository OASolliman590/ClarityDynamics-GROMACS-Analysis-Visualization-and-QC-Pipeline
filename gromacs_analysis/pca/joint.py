"""Joint principal component analysis in one shared essential subspace.

Per-system PCA gives each system its own eigenvectors, so "PC1" means a different
motion in every system and projections cannot be compared quantitatively. Here all
systems' frames are fitted to one common reference, one covariance matrix is built
from all of them, and every system is projected on the same eigenvectors
(Amadei et al. 1993; van Aalten et al. 1997). Similarity of the individual systems'
dominant motions is quantified by the root mean square inner product (RMSIP,
Amadei et al. 1999).

The numeric core works on plain arrays shaped (frames, atoms, 3); MDAnalysis is only
needed to read trajectories (``load_coordinates``).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class JointPCAResult:
    """Shared basis, per-system projections and subspace overlaps."""

    system_names: List[str]
    n_atoms: int
    eigenvalues: np.ndarray            # (k,) descending
    eigenvectors: np.ndarray           # (3N, k), sign-normalised
    explained_variance: np.ndarray     # (k,) fraction of total variance
    total_variance: float
    projections: Dict[str, np.ndarray]  # system -> (frames, k)
    rmsip_pairwise: np.ndarray         # (S, S)
    rmsip_to_joint: np.ndarray         # (S,)
    mean_structure: np.ndarray         # (N, 3)
    parameters: Dict = field(default_factory=dict)


def kabsch_fit(mobile: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Optimal rigid superposition of ``mobile`` (N, 3) onto ``reference`` (N, 3)."""
    mob_c = mobile - mobile.mean(axis=0)
    ref_c = reference - reference.mean(axis=0)
    h = mob_c.T @ ref_c
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    correction = np.diag([1.0, 1.0, d])
    rotation = vt.T @ correction @ u.T
    return mob_c @ rotation.T + reference.mean(axis=0)


def _fit_all(frames: np.ndarray, reference: np.ndarray) -> np.ndarray:
    return np.stack([kabsch_fit(frame, reference) for frame in frames])


def iterative_fit(
    coords_by_system: Mapping[str, np.ndarray],
    iterations: int = 2,
) -> Dict[str, np.ndarray]:
    """Fit every frame of every system to a common reference.

    Pass 1 fits to frame 0 of the first system; later passes fit to the mean
    structure of all fitted frames, removing the bias of an arbitrary reference.
    """
    names = list(coords_by_system)
    reference = np.asarray(coords_by_system[names[0]][0], dtype=float)
    fitted: Dict[str, np.ndarray] = {}
    for _ in range(max(1, int(iterations))):
        fitted = {name: _fit_all(np.asarray(coords_by_system[name], dtype=float), reference) for name in names}
        reference = np.concatenate([fitted[name] for name in names]).mean(axis=0)
    return fitted


def _normalise_signs(vectors: np.ndarray) -> np.ndarray:
    """Make the largest-magnitude component of each eigenvector positive (deterministic output)."""
    idx = np.argmax(np.abs(vectors), axis=0)
    signs = np.sign(vectors[idx, np.arange(vectors.shape[1])])
    signs[signs == 0] = 1.0
    return vectors * signs


def _pca(flat: np.ndarray, k: int):
    centred = flat - flat.mean(axis=0)
    cov = centred.T @ centred / max(len(flat) - 1, 1)
    values, vectors = np.linalg.eigh(cov)
    order = np.argsort(values)[::-1]
    values = np.clip(values[order], 0.0, None)
    vectors = _normalise_signs(vectors[:, order][:, :k])
    return values, vectors


def rmsip(a: np.ndarray, b: np.ndarray, k: int) -> float:
    """Root mean square inner product of the first ``k`` columns of two orthonormal bases."""
    k = min(k, a.shape[1], b.shape[1])
    overlap = a[:, :k].T @ b[:, :k]
    return float(np.sqrt(np.sum(overlap ** 2) / k))


def joint_pca(
    coords_by_system: Mapping[str, np.ndarray],
    n_components: int = 10,
    n_overlap: int = 10,
    fit_iterations: int = 2,
) -> JointPCAResult:
    """Joint PCA of several systems on one shared basis.

    Args:
        coords_by_system: system name -> coordinates (frames, atoms, 3). All systems must
            use the same atom selection (same atoms, same order).
        n_components: number of shared PCs to project on and report.
        n_overlap: number of PCs used for RMSIP.
        fit_iterations: fitting passes (>= 1).
    """
    if not coords_by_system:
        raise ValueError("joint PCA needs at least one system")
    names = list(coords_by_system)
    shapes = {name: np.asarray(coords_by_system[name]).shape for name in names}
    for name, shape in shapes.items():
        if len(shape) != 3 or shape[2] != 3:
            raise ValueError(f"System {name}: coordinates must be (frames, atoms, 3), got {shape}")
        if shape[0] < 2:
            raise ValueError(f"System {name}: at least 2 frames are required, got {shape[0]}")
    atom_counts = {name: shape[1] for name, shape in shapes.items()}
    if len(set(atom_counts.values())) > 1:
        detail = ", ".join(f"{name}: {count} atoms" for name, count in atom_counts.items())
        raise ValueError(f"All systems must select the same atoms for joint PCA ({detail})")
    if len(names) == 1:
        logger.warning("Joint PCA with a single system: comparison needs two or more systems.")

    n_atoms = next(iter(atom_counts.values()))
    dim = 3 * n_atoms
    total_frames = sum(shape[0] for shape in shapes.values())
    k = int(max(1, min(n_components, dim, total_frames - 1)))

    fitted = iterative_fit(coords_by_system, iterations=fit_iterations)
    flat = {name: fitted[name].reshape(len(fitted[name]), dim) for name in names}
    stacked = np.concatenate([flat[name] for name in names])
    mean = stacked.mean(axis=0)

    values, vectors = _pca(stacked, k)
    total = float(values.sum())
    eigenvalues = values[:k]
    projections = {name: (flat[name] - mean) @ vectors for name in names}

    k_overlap = int(max(1, min(n_overlap, dim, min(shape[0] for shape in shapes.values()) - 1)))
    own_bases = {name: _pca(flat[name], k_overlap)[1] for name in names}
    pairwise = np.array([[rmsip(own_bases[a], own_bases[b], k_overlap) for b in names] for a in names])
    to_joint = np.array([rmsip(own_bases[name], vectors, k_overlap) for name in names])

    return JointPCAResult(
        system_names=names,
        n_atoms=n_atoms,
        eigenvalues=eigenvalues,
        eigenvectors=vectors,
        explained_variance=eigenvalues / total if total > 0 else np.zeros_like(eigenvalues),
        total_variance=total,
        projections=projections,
        rmsip_pairwise=pairwise,
        rmsip_to_joint=to_joint,
        mean_structure=mean.reshape(n_atoms, 3),
        parameters={
            "n_components": k,
            "n_overlap": k_overlap,
            "fit_iterations": int(fit_iterations),
            "frames": {name: int(shape[0]) for name, shape in shapes.items()},
            "n_atoms": n_atoms,
        },
    )


def load_coordinates(
    topology: Path,
    trajectories: Sequence[Path],
    selection: str = "protein and name CA",
    stride: int = 1,
    max_frames: Optional[int] = None,
) -> np.ndarray:
    """Read selected coordinates (Å) from one or more trajectories with MDAnalysis."""
    try:
        import MDAnalysis as mda  # noqa: WPS433
    except ImportError as exc:  # pragma: no cover - exercised only without MDAnalysis
        raise ImportError("MDAnalysis is required to read trajectories for joint PCA") from exc

    blocks: List[np.ndarray] = []
    for trajectory in trajectories:
        universe = mda.Universe(str(topology), str(trajectory))
        atoms = universe.select_atoms(selection)
        if atoms.n_atoms == 0:
            raise ValueError(f"Selection '{selection}' matched no atoms in {topology}")
        blocks.append(np.array([atoms.positions.copy() for _ in universe.trajectory[:: max(1, int(stride))]]))
    coords = np.concatenate(blocks)
    if max_frames is not None and len(coords) > max_frames:
        # Evenly spaced subsample keeps the whole time span represented.
        idx = np.linspace(0, len(coords) - 1, int(max_frames)).round().astype(int)
        coords = coords[idx]
    return coords
