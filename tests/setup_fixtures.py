"""Build small CHARMM-GUI-like project trees for setup/discovery tests."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence


def psf(resnames: Sequence[str]) -> str:
    lines = ["PSF", "", "       1 !NTITLE", " REMARKS test", "", f"{len(resnames):8d} !NATOM"]
    for i, res in enumerate(resnames, start=1):
        lines.append(f"{i:8d} PROA {i:<8d} {res:<6} CA  CT1  0.000 12.011 0")
    lines += ["", "       0 !NBOND: bonds", ""]
    return "\n".join(lines) + "\n"


def pdb(resnames: Sequence[str]) -> str:
    lines = []
    for i, res in enumerate(resnames, start=1):
        record = "ATOM  " if res in {"ALA", "GLY", "HSD", "LYS"} else "HETATM"
        lines.append(f"{record}{i:5d}  CA  {res:<4}A{i:4d}    {i:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00           C")
    return "\n".join(lines) + "\nEND\n"


def gro(resnames: Sequence[str]) -> str:
    lines = ["test system", f"{len(resnames):5d}"]
    for i, res in enumerate(resnames, start=1):
        lines.append(f"{i:5d}{res:<5}{'CA':>5}{i:5d}{i * 0.1:8.3f}{0.0:8.3f}{0.0:8.3f}")
    lines.append("   5.00000   5.00000   5.00000")
    return "\n".join(lines) + "\n"


PROTEIN = ["ALA", "GLY", "HSD", "LYS"]
SOLVENT = ["TIP3", "SOD", "CLA"]


def write(path: Path, text: str = "", size: int = 0) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if size:
        path.write_bytes(b"x" * size)
    else:
        path.write_text(text)


def build_project(root: Path) -> Path:
    """Standard multi-ligand project.

    - HOLO_LigA_1..3: psf with INH + NAD (two ligands), raw and noPBC xtc, metric files, MMPBSA in rep 1
    - HOLO_LigB_rep1..2: gro with LGB, single xtc
    - APO_1..3: pdb with protein + HEM cofactor shared with LigB (not a ligand)
    - LigC/rep1, LigC/rep2: nested layout, pdb with LGC
    - notes/ (no trajectory) and topology_only/ (no trajectory): skipped
    """
    for rep in (1, 2, 3):
        d = root / f"HOLO_LigA_{rep}"
        write(d / "step3_input.psf", psf(PROTEIN + ["INH", "NAD"] + SOLVENT))
        write(d / "step3_input.pdb", pdb(PROTEIN + ["INH", "NAD"]))
        write(d / "step5_production.xtc", size=500)
        write(d / "step5_production_noPBC.xtc", size=400)
        write(d / "RMSD_apo_.dat", "1\n2\n")
        write(d / "RMSF_.dat", "1 1.0\n")
        if rep == 1:
            write(d / "FINAL_RESULTS_MMPBSA.dat", "Delta (Complex - Receptor - Ligand):\n")
    for rep in (1, 2):
        d = root / f"HOLO_LigB_rep{rep}"
        write(d / "md.gro", gro(PROTEIN + ["HEM", "LGB"] + SOLVENT))
        write(d / "md.xtc", size=300)
    for rep in (1, 2, 3):
        d = root / f"APO_{rep}"
        write(d / "step3_input.pdb", pdb(PROTEIN + ["HEM", "SOD"]))
        write(d / "prod.xtc", size=300)
        write(d / "analysis" / "RMSD_apo_.dat", "1\n")
    for rep in (1, 2):
        d = root / "LigC" / f"rep{rep}"
        write(d / "system.pdb", pdb(PROTEIN + ["LGC"]))
        write(d / "traj.dcd", size=100)
    write(root / "notes" / "readme.txt", "x")
    write(root / "topology_only" / "model.pdb", pdb(PROTEIN))
    return root


def names(systems: Iterable) -> list:
    return [s.key for s in systems]
