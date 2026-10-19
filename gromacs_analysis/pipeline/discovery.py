"""Automatic discovery of per-ligand systems, replicates and ligand residues.

Used by the ``prolif_batch`` and ``ttclust_batch`` stages so that a project laid
out as ``<root>/HOLO_<ligand>_<rep>/`` is analysed per ligand without listing
every system by hand, and so that complexes with several ligands (e.g. an
inhibitor plus a cofactor) are fingerprinted once per ligand.
"""

from __future__ import annotations

import logging
import re
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

logger = logging.getLogger(__name__)

PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # Protonation / force-field variants
    "HSD", "HSE", "HSP", "HID", "HIE", "HIP", "HISD", "HISE", "HISH", "CYX", "CYM", "ASH",
    "GLH", "LYN", "ASPP", "GLUP", "LSN", "SEC", "PYL", "MSE",
    # Caps
    "ACE", "NME", "NMA", "NH2",
}
SOLVENT_ION_RESNAMES = {
    "SOL", "WAT", "HOH", "TIP3", "TIP3P", "TIP4", "TIP4P", "TIP5", "SPC", "SPCE", "T3P", "T4P",
    "NA", "NA+", "SOD", "K", "K+", "POT", "CL", "CL-", "CLA", "MG", "MG2", "CA", "CAL", "ZN",
    "ZN2", "CS", "LI", "RB", "BR", "IOD", "F", "NA1", "CL1",
}
LIPID_RESNAMES = {"POPC", "POPE", "POPG", "POPS", "DPPC", "DMPC", "DOPC", "DOPE", "CHL1", "CHOL", "PA", "PC", "OL"}


def _pattern_to_regex(pattern: str) -> re.Pattern:
    """``HOLO_{ligand}_{rep}`` -> regex with named groups (other text matched literally)."""
    out = ""
    pos = 0
    for match in re.finditer(r"\{(\w+)\}", pattern):
        out += re.escape(pattern[pos:match.start()])
        name = match.group(1)
        out += rf"(?P<{name}>\d+)" if name == "rep" else rf"(?P<{name}>[^/]+?)"
        pos = match.end()
    out += re.escape(pattern[pos:])
    return re.compile(rf"^{out}$")


def _pattern_to_glob(pattern: str) -> str:
    return re.sub(r"\{\w+\}", "*", pattern)


def discover_systems(
    root: Path,
    pattern: str,
    topology: str,
    trajectory: str,
    include: Optional[Sequence[str]] = None,
    exclude: Optional[Sequence[str]] = None,
) -> List[Dict]:
    """Group replicate directories by ligand.

    Args:
        root: Directory searched for systems.
        pattern: Relative directory pattern with ``{ligand}`` and optionally
            ``{rep}``, e.g. ``"HOLO_{ligand}_{rep}"`` or ``"{ligand}/rep{rep}"``.
        topology: Topology file name inside each replicate directory (the first
            replicate's topology is used for the system).
        trajectory: Trajectory file name inside each replicate directory.

    Returns:
        One entry per ligand: ``{"name", "ligand_name", "topology", "trajectories", "directories"}``
        with replicates sorted numerically. Directories missing the trajectory are skipped
        with a warning.
    """
    root = Path(root)
    if "{ligand}" not in pattern:
        raise ValueError("auto_discover.pattern must contain a {ligand} placeholder")
    regex = _pattern_to_regex(pattern)
    include_set = {str(v) for v in include} if include else None
    exclude_set = {str(v) for v in exclude or []}

    grouped: "OrderedDict[str, List]" = OrderedDict()
    for candidate in sorted(root.glob(_pattern_to_glob(pattern))):
        if not candidate.is_dir():
            continue
        rel = candidate.relative_to(root).as_posix()
        match = regex.match(rel)
        if not match:
            continue
        ligand = match.group("ligand")
        if include_set is not None and ligand not in include_set:
            continue
        if ligand in exclude_set:
            continue
        rep = int(match.groupdict().get("rep") or 1)
        traj = candidate / trajectory
        if not traj.exists():
            logger.warning(f"Skipping {candidate}: trajectory {trajectory} not found")
            continue
        grouped.setdefault(ligand, []).append((rep, candidate))

    systems: List[Dict] = []
    for ligand, reps in grouped.items():
        reps.sort(key=lambda item: item[0])
        dirs = [path for _, path in reps]
        topo = next((d / topology for d in dirs if (d / topology).exists()), None)
        if topo is None:
            logger.warning(f"Skipping ligand {ligand}: topology {topology} not found in {dirs[0]}")
            continue
        systems.append({
            "name": ligand,
            "ligand_name": ligand,
            "topology": str(topo),
            "trajectories": [str(d / trajectory) for d in dirs],
            "directories": [str(d) for d in dirs],
        })
    if not systems:
        raise ValueError(
            f"auto_discover found no systems under {root} matching '{pattern}' "
            f"with {trajectory} and {topology}"
        )
    logger.info("Auto-discovered %d system(s): %s", len(systems), ", ".join(s["name"] for s in systems))
    return systems


def _resnames_from_text_topology(path: Path) -> List[str]:
    suffix = path.suffix.lower()
    names: "OrderedDict[str, None]" = OrderedDict()
    lines = path.read_text(errors="replace").splitlines()
    if suffix == ".pdb":
        for line in lines:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 20:
                names[line[17:21].strip()] = None
    elif suffix == ".gro":
        for line in lines[2:-1]:
            if len(line) >= 10:
                names[line[5:10].strip()] = None
    elif suffix == ".psf":
        in_atoms = False
        for line in lines:
            if "!NATOM" in line:
                in_atoms = True
                continue
            if in_atoms:
                if not line.strip() or "!" in line:
                    break
                parts = line.split()
                if len(parts) >= 4:
                    names[parts[3]] = None
    else:
        raise ValueError(f"Cannot read residue names from {path.suffix} without MDAnalysis")
    return list(names)


def detect_ligands(
    topology: Path,
    exclude: Iterable[str] = (),
    include_lipids: bool = False,
) -> List[Dict[str, str]]:
    """Residue names that are neither protein, solvent, ions (nor lipids) -> one ligand each.

    Returns ``[{"name": RESNAME, "selection": "resname RESNAME"}, ...]`` in topology order.
    """
    topology = Path(topology)
    resnames: List[str]
    try:
        import MDAnalysis as mda  # noqa: WPS433

        universe = mda.Universe(str(topology))
        resnames = list(OrderedDict.fromkeys(str(r) for r in universe.residues.resnames))
    except Exception:
        resnames = _resnames_from_text_topology(topology)

    skip = {name.upper() for name in exclude} | PROTEIN_RESNAMES | SOLVENT_ION_RESNAMES
    if not include_lipids:
        skip |= LIPID_RESNAMES
    ligands = [
        {"name": name, "selection": f"resname {name}"}
        for name in resnames
        if name and name.upper() not in skip
    ]
    return ligands


def normalise_ligands(spec, topology: Optional[Path], exclude: Iterable[str] = ()) -> List[Dict[str, str]]:
    """Resolve a ``ligands`` setting: ``"auto"``, a list of names, or a list of {name, selection}."""
    if spec in (None, "", []):
        return []
    if isinstance(spec, str):
        if spec.lower() != "auto":
            return [{"name": spec, "selection": f"resname {spec}"}]
        if topology is None:
            raise ValueError("ligands: auto requires a topology")
        found = detect_ligands(topology, exclude=exclude)
        if not found:
            raise ValueError(f"ligands: auto found no ligand residues in {topology}")
        return found
    out = []
    for item in spec:
        if isinstance(item, dict):
            name = str(item.get("name") or item.get("resname") or "").strip()
            if not name:
                raise ValueError(f"Ligand entry needs a name: {item}")
            out.append({"name": name, "selection": item.get("selection") or f"resname {name}"})
        else:
            out.append({"name": str(item), "selection": f"resname {item}"})
    return out
