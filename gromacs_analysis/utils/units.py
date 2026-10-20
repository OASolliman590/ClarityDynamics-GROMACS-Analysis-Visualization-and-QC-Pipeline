"""Unit detection for GROMACS ``.xvg`` metrics and conversion to display units.

GROMACS tools write lengths in nm and areas in nm²; the VMD helper writes Å and Å².
The unit a file declares in its ``@ yaxis label`` header is converted to the unit
shown on the plot axis (taken from the metric's ``ylabel``), and every conversion
is recorded so it can be reported.
"""

from __future__ import annotations

import re
from typing import Dict, Optional, Tuple

LENGTH_NM = "nm"
AREA_NM2 = "nm^2"
LENGTH_A = "Å"
AREA_A2 = "Å^2"

_FACTORS: Dict[Tuple[str, str], float] = {
    (LENGTH_NM, LENGTH_A): 10.0,
    (LENGTH_A, LENGTH_NM): 0.1,
    (AREA_NM2, AREA_A2): 100.0,
    (AREA_A2, AREA_NM2): 0.01,
}


def normalise_unit(text: Optional[str]) -> Optional[str]:
    """Map unit spellings (incl. GROMACS ``nm\\S2\\N`` formatting) to canonical names.

    Returns None when ``text`` names no unit, and the stripped text itself when it
    names a unit that is not a length or area.
    """
    if text is None:
        return None
    raw = str(text).strip()
    if not raw:
        return None
    t = raw.replace("\\S2\\N", "^2").replace("\\S2", "^2").replace("²", "^2").replace("**2", "^2")
    t = t.replace(" ", "")
    lower = t.lower()
    if lower in {"nm", "nanometer", "nanometers", "nanometre", "nanometres"}:
        return LENGTH_NM
    if lower in {"nm^2", "nm2", "nm*nm"}:
        return AREA_NM2
    if t in {"Å", "A", "Ang", "Angstrom", "angstrom", "Angstroms", "angstroms", "\\cE\\C", "\\AA"}:
        return LENGTH_A
    if t in {"Å^2", "A^2", "A2", "Angstrom^2", "angstrom^2", "\\cE\\C^2"}:
        return AREA_A2
    return raw


def _unit_from_label(label: str) -> Optional[str]:
    """Unit in a label such as ``RMSD (nm)``, ``Area (nm\\S2\\N)`` or a bare ``Å``."""
    if label is None:
        return None
    match = re.search(r"\(([^()]*)\)\s*$", label.strip())
    if match:
        return normalise_unit(match.group(1))
    canonical = normalise_unit(label)
    if canonical in {LENGTH_NM, AREA_NM2, LENGTH_A, AREA_A2}:
        return canonical
    return None


def parse_xvg_units(content: str) -> Tuple[Optional[str], Optional[str]]:
    """Return (x_unit, y_unit) declared in an ``.xvg`` header (None when absent)."""
    x_unit = y_unit = None
    x_match = re.search(r'^@\s*xaxis\s+label\s+"([^"]*)"', content, re.MULTILINE)
    y_match = re.search(r'^@\s*yaxis\s+label\s+"([^"]*)"', content, re.MULTILINE)
    if x_match:
        x_unit = _unit_from_label(x_match.group(1))
    if y_match:
        y_unit = _unit_from_label(y_match.group(1))
    return x_unit, y_unit


def display_unit_from_label(label: Optional[str]) -> Optional[str]:
    """Display unit implied by a plot axis label: 'Å', 'RMSD (Å)', 'Å²', 'Area (nm²)', ..."""
    if not label:
        return None
    unit = _unit_from_label(label)
    if unit in {LENGTH_NM, AREA_NM2, LENGTH_A, AREA_A2}:
        return unit
    # Labels such as "Å" with trailing text: look for a unit token anywhere.
    for token, canonical in (("Å²", AREA_A2), ("Å^2", AREA_A2), ("nm²", AREA_NM2), ("nm^2", AREA_NM2),
                             ("Å", LENGTH_A), ("nm", LENGTH_NM)):
        if re.search(rf"(?<![A-Za-z]){re.escape(token)}(?![A-Za-z^²])", label):
            return canonical
    return None


def conversion_factor(declared: Optional[str], display: Optional[str]) -> Tuple[float, str]:
    """Factor converting values in ``declared`` to ``display`` units, and a status.

    Status is one of ``converted``, ``matched``, ``assumed`` (nothing declared, values
    left unchanged) or ``unrecognised`` (declared unit cannot be converted to the display unit).
    """
    if declared is None:
        return 1.0, "assumed"
    if display is None:
        return 1.0, "unrecognised" if declared not in {LENGTH_NM, AREA_NM2, LENGTH_A, AREA_A2} else "matched"
    if declared == display:
        return 1.0, "matched"
    factor = _FACTORS.get((declared, display))
    if factor is None:
        return 1.0, "unrecognised"
    return factor, "converted"
