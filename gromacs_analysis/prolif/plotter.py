"""Plotting helpers for ProLIF outputs."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Optional

import plotly.graph_objects as go

from ..utils.figure_export import write_static

logger = logging.getLogger(__name__)


class ProlifPlotter:
    def __init__(self, plot_config, protein_name: str):
        self.plot_config = plot_config
        self.protein_name = protein_name

    def create_heatmap(
        self,
        z,
        x_labels,
        y_labels,
        title: str,
        colorscale: str = "Viridis",
        zmin: Optional[float] = None,
        zmax: Optional[float] = None,
    ) -> go.Figure:
        fig = go.Figure(
            data=[
                go.Heatmap(
                    z=z,
                    x=x_labels,
                    y=y_labels,
                    colorscale=colorscale,
                    zmin=zmin,
                    zmax=zmax,
                )
            ]
        )
        fig.update_layout(
            title=title,
            font_family=self.plot_config.font_family,
            font_size=self.plot_config.font_size,
            template=self.plot_config.template,
            xaxis_title="Residue" if x_labels else None,
            yaxis_title="Ligand" if y_labels else None,
        )
        return fig

    def save_figure(self, fig, output_dir: Path, name: str, formats: Optional[Iterable[str]] = None) -> None:
        output_dir.mkdir(parents=True, exist_ok=True)
        formats = list(formats or self.plot_config.save_formats)

        normalized_fig = self._normalize_figure(fig)
        for fmt in formats:
            file_path = output_dir / f"{name}.{fmt}"
            try:
                if fmt == "html":
                    if hasattr(normalized_fig, "write_html"):
                        normalized_fig.write_html(str(file_path))
                    elif hasattr(normalized_fig, "to_html"):
                        file_path.write_text(normalized_fig.to_html())
                    elif hasattr(normalized_fig, "save"):
                        normalized_fig.save(str(file_path))
                    else:
                        logger.warning(f"Cannot export HTML for {name}")
                elif fmt in ("svg", "png"):
                    if hasattr(normalized_fig, "write_image"):
                        write_static(
                            normalized_fig,
                            file_path,
                            fmt,
                            self.plot_config.width,
                            self.plot_config.height,
                            self.plot_config.scale,
                        )
                    elif hasattr(normalized_fig, "savefig"):
                        normalized_fig.savefig(str(file_path), format=fmt, bbox_inches="tight")
                    elif hasattr(normalized_fig, "save"):
                        normalized_fig.save(str(file_path))
                    else:
                        logger.warning(f"Cannot export {fmt} for {name}")
                logger.info(f"Saved {fmt} plot: {file_path}")
            except Exception as exc:
                if fmt == "svg" and hasattr(normalized_fig, "savefig"):
                    try:
                        normalized_fig.savefig(str(file_path), format=fmt, bbox_inches="tight")
                        logger.info(f"Saved SVG via matplotlib fallback: {file_path}")
                        continue
                    except Exception:
                        pass
                logger.error(f"Error saving {fmt} plot {file_path}: {exc}")

    def save_matplotlib(self, fig, output_dir: Path, name: str, fmt: str = "svg") -> None:
        output_dir.mkdir(parents=True, exist_ok=True)
        file_path = output_dir / f"{name}.{fmt}"
        fig.savefig(str(file_path), format=fmt, bbox_inches="tight")

    def save_tanimoto_heatmap(self, matrix, labels, output_dir: Path, name: str) -> None:
        fig = self.create_heatmap(
            matrix,
            labels,
            labels,
            title=f"{self.protein_name} · ProLIF Tanimoto Similarity",
            colorscale="RdYlBu",
            zmin=0,
            zmax=1,
        )
        self.save_figure(fig, output_dir, name)

    @staticmethod
    def _normalize_figure(fig):
        # ProLIF barcode returns matplotlib Axes; use underlying Figure for saving.
        if hasattr(fig, "figure") and not hasattr(fig, "savefig"):
            return fig.figure
        return fig
