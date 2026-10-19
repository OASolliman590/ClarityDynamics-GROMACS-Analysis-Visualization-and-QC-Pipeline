"""Write small, physically plausible protein-ligand systems to disk for integration tests.

Two helical-ish chains (A, B) of glycine-like residues with backbone atoms and a
three-atom ligand ``LIG``. Residue numbers restart in each chain, so chain-aware
code paths are exercised. The ligand can drift away linearly so distances have a
known trend. Coordinates are in Å (MDAnalysis convention).
"""

from __future__ import annotations

from pathlib import Path
from typing import Tuple

import numpy as np


def build_system(
    out_dir: Path,
    n_res_per_chain: int = 8,
    n_frames: int = 20,
    drift_per_frame: float = 0.5,
    noise: float = 0.1,
    seed: int = 0,
    name: str = "sys",
) -> Tuple[Path, Path]:
    """Return (topology.pdb, trajectory.xtc) paths."""
    import MDAnalysis as mda

    rng = np.random.default_rng(seed)
    atom_names = ["N", "CA", "C", "O"]
    chains = ["A", "B"]
    n_protein_res = n_res_per_chain * len(chains)
    n_atoms = n_protein_res * len(atom_names) + 3

    resindex = np.concatenate([np.repeat(np.arange(n_protein_res), len(atom_names)), np.full(3, n_protein_res)])
    segindex = np.array([i // n_res_per_chain for i in range(n_protein_res)] + [2])
    u = mda.Universe.empty(n_atoms, n_residues=n_protein_res + 1, n_segments=3,
                           atom_resindex=resindex, residue_segindex=segindex, trajectory=True)
    u.add_TopologyAttr("names", atom_names * n_protein_res + ["C1", "C2", "O1"])
    u.add_TopologyAttr("types", [n[0] for n in atom_names] * n_protein_res + ["C", "C", "O"])
    u.add_TopologyAttr("elements", [n[0] for n in atom_names] * n_protein_res + ["C", "C", "O"])
    u.add_TopologyAttr("resnames", ["GLY"] * n_protein_res + ["LIG"])
    u.add_TopologyAttr("resids", list(range(1, n_res_per_chain + 1)) * len(chains) + [1])
    u.add_TopologyAttr("segids", chains + ["L"])
    u.add_TopologyAttr("chainIDs", [chains[i // (n_res_per_chain * 4)] for i in range(n_protein_res * 4)] + ["L"] * 3)
    u.add_TopologyAttr("masses", [14.0, 12.0, 12.0, 16.0] * n_protein_res + [12.0, 12.0, 16.0])
    u.add_TopologyAttr("occupancies", np.ones(n_atoms))
    u.add_TopologyAttr("tempfactors", np.zeros(n_atoms))
    u.add_TopologyAttr("altLocs", [""] * n_atoms)
    u.add_TopologyAttr("icodes", [""] * n_protein_res + [""])

    base = np.zeros((n_atoms, 3))
    for c, offset in enumerate((0.0, 10.0)):
        for r in range(n_res_per_chain):
            idx = (c * n_res_per_chain + r) * 4
            center = np.array([offset, 3.8 * r, 0.0])
            base[idx:idx + 4] = center + np.array([[-1.2, 0, 0], [0, 0, 0], [1.2, 0.3, 0], [1.8, 1.2, 0]])
    ligand_start = np.array([5.0, 3.8 * (n_res_per_chain / 2), 3.0])
    base[-3:] = ligand_start + np.array([[0, 0, 0], [1.4, 0, 0], [2.0, 1.1, 0]])

    out_dir.mkdir(parents=True, exist_ok=True)
    topology = out_dir / f"{name}.pdb"
    trajectory = out_dir / f"{name}.xtc"
    u.dimensions = np.array([60.0, 60.0, 60.0, 90.0, 90.0, 90.0])
    u.atoms.positions = base
    u.atoms.write(str(topology))
    with mda.Writer(str(trajectory), n_atoms) as writer:
        for frame in range(n_frames):
            coords = base + rng.normal(0, noise, base.shape)
            coords[-3:] += np.array([0.0, 0.0, drift_per_frame * frame])
            u.atoms.positions = coords
            u.trajectory.ts.time = frame * 10.0  # ps
            u.dimensions = np.array([60.0, 60.0, 60.0, 90.0, 90.0, 90.0])
            writer.write(u.atoms)
    return topology, trajectory
