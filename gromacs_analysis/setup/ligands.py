"""Ligand residue detection from topology files without transferring them.

Each command prints the unique residue names of one topology, so even a
multi-megabyte solvated system costs one short round trip on a remote host.
"""

from __future__ import annotations

import shlex
from typing import Iterable, List, Optional

from ..pipeline.discovery import LIPID_RESNAMES, PROTEIN_RESNAMES, SOLVENT_ION_RESNAMES

# Formats whose residue names can be read as text, in order of preference.
TEXT_TOPOLOGIES = (".psf", ".pdb", ".gro")


def resname_command(path: str) -> Optional[str]:
    """Shell command printing one residue name per line for a text topology."""
    quoted = shlex.quote(path)
    lower = path.lower()
    if lower.endswith(".pdb"):
        return f"awk '/^(ATOM|HETATM)/ {{print substr($0,18,4)}}' {quoted} | sort -u"
    if lower.endswith(".gro"):
        # Atom lines: 5-char resid, 5-char resname; the title, count and box lines are skipped.
        return (
            f"awk 'NR>2 && length($0)>=44 && substr($0,1,5) ~ /^ *[0-9]+$/ "
            f"{{print substr($0,6,5)}}' {quoted} | sort -u"
        )
    if lower.endswith(".psf"):
        return (
            f"awk '/!NATOM/ {{on=1; next}} on && /!/ {{exit}} on && NF>=4 {{print $4}}' {quoted} | sort -u"
        )
    return None


def parse_resnames(output: str) -> List[str]:
    names = []
    for line in output.splitlines():
        name = line.strip()
        if name and name not in names:
            names.append(name)
    return names


def is_ligand_resname(name: str, include_lipids: bool = False) -> bool:
    upper = name.strip().upper()
    if not upper or upper in PROTEIN_RESNAMES or upper in SOLVENT_ION_RESNAMES:
        return False
    if not include_lipids and upper in LIPID_RESNAMES:
        return False
    return True


def filter_ligands(resnames: Iterable[str], shared: Iterable[str] = ()) -> List[str]:
    """Ligand candidates: not protein/water/ion/lipid and not present in the APO reference."""
    shared_upper = {s.strip().upper() for s in shared}
    return sorted(
        {n.strip() for n in resnames if is_ligand_resname(n) and n.strip().upper() not in shared_upper}
    )


def pick_text_topology(paths: Iterable[str]) -> Optional[str]:
    by_ext = {}
    for path in paths:
        for ext in TEXT_TOPOLOGIES:
            if path.lower().endswith(ext) and ext not in by_ext:
                by_ext[ext] = path
    for ext in TEXT_TOPOLOGIES:
        if ext in by_ext:
            return by_ext[ext]
    return None
