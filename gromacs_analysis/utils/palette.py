"""Shared, colour-blind-validated plotting palette.

Categorical hues are assigned to systems in a fixed order (never cycled) so
that a ligand keeps its colour across every figure of a run. The APO/reference
system is drawn in a neutral grey so ligand-bound systems stand out.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

RGB = Tuple[int, int, int]

# Validated categorical order (adjacent-pair CVD ΔE >= 8 on a light surface).
CATEGORICAL_HEX: List[str] = [
    "#2a78d6",  # blue
    "#eb6834",  # orange
    "#1baf7a",  # aqua
    "#eda100",  # yellow
    "#e87ba4",  # magenta
    "#008300",  # green
    "#6250d6",  # violet
    "#e34948",  # red
]
APO_HEX = "#6f6e69"        # neutral reference system
TEXT_HEX = "#0b0b0b"
MUTED_HEX = "#8a8984"      # annotations, excluded regions

# Line dashes used as a secondary encoding once all eight hues are taken.
DASHES = ["solid", "dash", "dot", "dashdot"]

# Perceptually uniform colour scale for continuous surfaces (FEL, heatmaps).
SEQUENTIAL_SCALE = "Viridis"
OCCUPANCY_SCALE = [
    [0.0, "#f4f6fa"],
    [0.25, "#b7d3f6"],
    [0.5, "#6da7ec"],
    [0.75, "#2a78d6"],
    [1.0, "#104281"],
]


def hex_to_rgb(value: str) -> RGB:
    value = value.lstrip("#")
    return tuple(int(value[i:i + 2], 16) for i in (0, 2, 4))  # type: ignore[return-value]


def rgb_to_hex(rgb: Sequence[int]) -> str:
    return "#{:02x}{:02x}{:02x}".format(*[int(max(0, min(255, c))) for c in rgb[:3]])


def rgba(color, alpha: float) -> str:
    rgb = hex_to_rgb(color) if isinstance(color, str) else tuple(color)
    return f"rgba({rgb[0]},{rgb[1]},{rgb[2]},{alpha})"


def categorical(index: int) -> str:
    """Hue for the index-th series; beyond eight hues the caller adds a dash."""
    return CATEGORICAL_HEX[index % len(CATEGORICAL_HEX)]


def dash_for(index: int) -> str:
    return DASHES[(index // len(CATEGORICAL_HEX)) % len(DASHES)]


def assign_system_colors(
    names: Iterable[str],
    is_apo: Optional[Dict[str, bool]] = None,
    explicit: Optional[Dict[str, Sequence[int]]] = None,
) -> Dict[str, RGB]:
    """Map system names to RGB in fixed order: APO -> grey, ligands -> categorical slots."""
    is_apo = is_apo or {}
    explicit = explicit or {}
    colors: Dict[str, RGB] = {}
    slot = 0
    for name in names:
        if name in explicit and explicit[name] is not None:
            colors[name] = tuple(int(c) for c in explicit[name])[:3]  # type: ignore[assignment]
            continue
        if is_apo.get(name, "APO" in str(name).upper()):
            colors[name] = hex_to_rgb(APO_HEX)
            continue
        colors[name] = hex_to_rgb(categorical(slot))
        slot += 1
    return colors


def mix(color, other, amount: float) -> str:
    """Blend ``color`` toward ``other`` (0 = color, 1 = other)."""
    a = hex_to_rgb(color) if isinstance(color, str) else tuple(color)
    b = hex_to_rgb(other) if isinstance(other, str) else tuple(other)
    return rgb_to_hex([round(x + (y - x) * amount) for x, y in zip(a, b)])


def single_hue_scale(color) -> List[List]:
    """Light-to-dark ramp of one hue for ordered values such as simulation time.

    The light end stays clearly visible on a white surface (no near-white start).
    """
    return [
        [0.0, mix(color, "#ffffff", 0.55)],
        [0.5, mix(color, "#ffffff", 0.0)],
        [1.0, mix(color, "#000000", 0.45)],
    ]
