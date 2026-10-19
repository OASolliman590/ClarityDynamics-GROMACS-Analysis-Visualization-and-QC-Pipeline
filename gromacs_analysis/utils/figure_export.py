"""Consistent static export for Plotly figures.

Plotly lays text out at the figure's pixel size, so enlarging the canvas
(``width * scale``) shrinks fonts relative to the plot. Instead we keep the
layout size and pass ``scale`` as a pure resolution multiplier, and figures that
size themselves to their content (e.g. heatmaps with many residues) keep their
own width/height.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple


def figure_size(fig, default_width: int, default_height: int) -> Tuple[int, int]:
    layout = getattr(fig, "layout", None)
    width = getattr(layout, "width", None) if layout is not None else None
    height = getattr(layout, "height", None) if layout is not None else None
    return int(width or default_width), int(height or default_height)


def write_static(
    fig,
    path: Path,
    fmt: Optional[str] = None,
    width: int = 1200,
    height: int = 600,
    scale: float = 2,
) -> None:
    """Write SVG/PNG/PDF at the figure's own size (or the defaults) and ``scale``."""
    fmt = fmt or Path(path).suffix.lstrip(".")
    w, h = figure_size(fig, width, height)
    fig.write_image(str(path), format=fmt, width=w, height=h, scale=scale)


def scaled_height(n_rows: int, per_row: int = 22, base: int = 220, minimum: int = 420, maximum: int = 4000) -> int:
    """Height that gives every row (residue, cluster, ...) a readable band."""
    return int(max(minimum, min(maximum, base + per_row * max(0, n_rows))))


def scaled_width(n_cols: int, per_col: int = 28, base: int = 320, minimum: int = 900, maximum: int = 4000) -> int:
    return int(max(minimum, min(maximum, base + per_col * max(0, n_cols))))
