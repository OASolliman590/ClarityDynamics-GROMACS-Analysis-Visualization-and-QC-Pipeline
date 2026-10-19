"""
Plotter for ligand-protein distance analysis.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

from .config import DistanceConfig
from ..utils.mpl_fallback import init_matplotlib
from ..utils.figure_export import StaticExportUnavailable, write_static

logger = logging.getLogger(__name__)


class DistancePlotter:
    def __init__(self, config: DistanceConfig):
        self.config = config
        pio.templates.default = config.plot_config.template

    def create_per_residue_plot(self, time: np.ndarray, residue: Dict) -> go.Figure:
        fig = go.Figure()
        mean = residue["mean"]
        std = residue["std"]
        label = residue["label"]

        fig.add_trace(go.Scatter(
            x=time,
            y=mean,
            mode="lines",
            name=label,
            line=dict(width=2),
        ))

        fig.add_trace(go.Scatter(
            x=np.concatenate([time, time[::-1]]),
            y=np.concatenate([mean + std, (mean - std)[::-1]]),
            fill='toself',
            opacity=0.3,
            line=dict(color='rgba(0,0,0,0)'),
            showlegend=False,
            hoverinfo='skip'
        ))

        fig.update_layout(
            title=f"Distance between Residue {label} and Ligand",
            xaxis_title="Time (ns)",
            yaxis_title="Distance (A)",
            font_family=self.config.plot_config.font_family,
            font_size=self.config.plot_config.font_size,
            template=self.config.plot_config.template,
        )

        return fig

    def create_combined_plot(self, time: np.ndarray, residues: List[Dict]) -> go.Figure:
        fig = go.Figure()

        for residue in residues:
            mean = residue["mean"]
            std = residue["std"]
            label = residue["label"]

            fig.add_trace(go.Scatter(
                x=time,
                y=mean,
                mode="lines",
                name=label,
                line=dict(width=2),
            ))

            if self.config.combined_show_std:
                fig.add_trace(go.Scatter(
                    x=np.concatenate([time, time[::-1]]),
                    y=np.concatenate([mean + std, (mean - std)[::-1]]),
                    fill='toself',
                    opacity=0.15,
                    line=dict(color='rgba(0,0,0,0)'),
                    showlegend=False,
                    hoverinfo='skip'
                ))

        fig.update_layout(
            title=f"Residue-Ligand Distance · {self.config.protein_name}-{self.config.ligand_name}",
            xaxis_title="Time (ns)",
            yaxis_title="Distance (A)",
            font_family=self.config.plot_config.font_family,
            font_size=self.config.plot_config.font_size,
            template=self.config.plot_config.template,
            legend=dict(
                orientation='h',
                y=-0.4,
                x=0.5,
                xanchor='center'
            )
        )

        return fig

    def save_plot(self, fig: go.Figure, output_path: Path, name: str,
                  fallback: Optional[Dict] = None) -> None:
        output_path.mkdir(parents=True, exist_ok=True)

        for fmt in self.config.plot_config.save_formats:
            file_path = output_path / f"{name}.{fmt}"
            try:
                if fmt == "html":
                    fig.write_html(str(file_path))
                elif fmt in ("svg", "png"):
                    write_static(fig, file_path, fmt, self.config.plot_config.width, self.config.plot_config.height, self.config.plot_config.scale)
                logger.info(f"Saved {fmt} plot: {file_path}")
            except Exception as exc:
                if fmt == "svg" and fallback:
                    if self._save_svg_fallback(file_path, fallback):
                        logger.info(f"Saved SVG via matplotlib fallback: {file_path}")
                        continue
                if isinstance(exc, StaticExportUnavailable):
                    logger.debug("Static export unavailable; skipped %s", file_path)
                    continue
                logger.error(f"Error saving {fmt} plot {file_path}: {exc}")

    def _save_svg_fallback(self, file_path: Path, fallback: Dict) -> bool:
        kind = fallback.get("kind")
        try:
            if kind == "single":
                return self._save_single_svg(file_path, fallback.get("time"), fallback.get("residue"))
            if kind == "combined":
                return self._save_combined_svg(file_path, fallback.get("time"), fallback.get("residues"))
        except Exception as exc:
            logger.error(f"Matplotlib fallback failed for {file_path}: {exc}")
        return False

    def _save_single_svg(self, file_path: Path, time: np.ndarray, residue: Dict) -> bool:
        if time is None or residue is None:
            return False

        _, plt = init_matplotlib()
        fig, ax = plt.subplots(
            figsize=(self.config.plot_config.width / 100, self.config.plot_config.height / 100),
            dpi=100 * self.config.plot_config.scale,
        )
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size

        mean = residue["mean"]
        std = residue["std"]
        label = residue["label"]

        ax.plot(time, mean, linewidth=2)
        ax.fill_between(time, mean - std, mean + std, alpha=0.3)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Distance (A)")
        ax.set_title(f"Distance between Residue {label} and Ligand")
        ax.grid(True, alpha=0.2)
        fig.tight_layout()
        fig.savefig(file_path, format="svg", bbox_inches="tight")
        plt.close(fig)
        return True

    def _save_combined_svg(self, file_path: Path, time: np.ndarray, residues: List[Dict]) -> bool:
        if time is None or not residues:
            return False

        _, plt = init_matplotlib()
        fig, ax = plt.subplots(
            figsize=(self.config.plot_config.width / 100, self.config.plot_config.height / 100),
            dpi=100 * self.config.plot_config.scale,
        )
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size

        for residue in residues:
            mean = residue["mean"]
            std = residue["std"]
            label = residue["label"]
            ax.plot(time, mean, linewidth=2, label=label)
            if self.config.combined_show_std:
                ax.fill_between(time, mean - std, mean + std, alpha=0.15)

        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Distance (A)")
        ax.set_title(f"Residue-Ligand Distance · {self.config.protein_name}-{self.config.ligand_name}")
        ax.grid(True, alpha=0.2)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, frameon=False)
        fig.subplots_adjust(bottom=0.25)
        fig.savefig(file_path, format="svg", bbox_inches="tight")
        plt.close(fig)
        return True


class DistanceComparisonPlotter:
    def __init__(self, plot_config, protein_name: str):
        self.plot_config = plot_config
        self.protein_name = protein_name
        pio.templates.default = plot_config.template

    def create_residue_comparison_plot(
        self,
        time: np.ndarray,
        residue_label: str,
        series: Dict[str, Dict[str, np.ndarray]],
        color_map: Dict[str, Tuple[int, int, int]],
        show_std: bool = False,
    ) -> go.Figure:
        fig = go.Figure()

        for ligand_name, stats in series.items():
            mean = stats["mean"]
            std = stats["std"]
            color = color_map.get(ligand_name)
            line_kwargs = dict(width=2)
            if color:
                line_kwargs["color"] = f"rgb({color[0]},{color[1]},{color[2]})"

            fig.add_trace(go.Scatter(
                x=time,
                y=mean,
                mode="lines",
                name=ligand_name,
                line=line_kwargs,
            ))

            if show_std:
                fig.add_trace(go.Scatter(
                    x=np.concatenate([time, time[::-1]]),
                    y=np.concatenate([mean + std, (mean - std)[::-1]]),
                    fill='toself',
                    opacity=0.2,
                    line=dict(color='rgba(0,0,0,0)'),
                    showlegend=False,
                    hoverinfo='skip'
                ))

        fig.update_layout(
            title=f"{self.protein_name} · {residue_label} Distance Comparison",
            xaxis_title="Time (ns)",
            yaxis_title="Distance (A)",
            font_family=self.plot_config.font_family,
            font_size=self.plot_config.font_size,
            template=self.plot_config.template,
            legend=dict(
                orientation='h',
                y=-0.4,
                x=0.5,
                xanchor='center'
            )
        )

        return fig

    def create_overall_comparison_plot(
        self,
        time: np.ndarray,
        series: Dict[str, Dict[str, np.ndarray]],
        color_map: Dict[str, Tuple[int, int, int]],
        show_std: bool = False,
    ) -> go.Figure:
        fig = go.Figure()

        for ligand_name, stats in series.items():
            mean = stats["mean"]
            std = stats["std"]
            color = color_map.get(ligand_name)
            line_kwargs = dict(width=2)
            if color:
                line_kwargs["color"] = f"rgb({color[0]},{color[1]},{color[2]})"

            fig.add_trace(go.Scatter(
                x=time,
                y=mean,
                mode="lines",
                name=ligand_name,
                line=line_kwargs,
            ))

            if show_std:
                fig.add_trace(go.Scatter(
                    x=np.concatenate([time, time[::-1]]),
                    y=np.concatenate([mean + std, (mean - std)[::-1]]),
                    fill='toself',
                    opacity=0.2,
                    line=dict(color='rgba(0,0,0,0)'),
                    showlegend=False,
                    hoverinfo='skip'
                ))

        fig.update_layout(
            title=f"{self.protein_name} · Overall Distance Comparison",
            xaxis_title="Time (ns)",
            yaxis_title="Distance (A)",
            font_family=self.plot_config.font_family,
            font_size=self.plot_config.font_size,
            template=self.plot_config.template,
            legend=dict(
                orientation='h',
                y=-0.4,
                x=0.5,
                xanchor='center'
            )
        )

        return fig

    def create_distance_heatmap(
        self,
        time: np.ndarray,
        residue_labels: List[str],
        matrix: np.ndarray,
        title: str,
        colorscale: str = "Viridis",
    ) -> go.Figure:
        fig = go.Figure(
            data=[
                go.Heatmap(
                    z=matrix,
                    x=residue_labels,
                    y=time,
                    colorscale=colorscale,
                    coloraxis="coloraxis",
                    zsmooth="best",
                )
            ]
        )
        fig.update_layout(
            title=title,
            xaxis_title="Residue",
            yaxis_title="Time (ns)",
            font_family=self.plot_config.font_family,
            font_size=self.plot_config.font_size,
            template=self.plot_config.template,
        )
        return fig

    def create_delta_heatmap(
        self,
        residue_labels: List[str],
        delta_values: np.ndarray,
        title: str,
        colorscale: str = "RdBu",
    ) -> go.Figure:
        fig = go.Figure(
            data=[
                go.Heatmap(
                    z=[delta_values],
                    x=residue_labels,
                    y=["Δ"],
                    colorscale=colorscale,
                    coloraxis="coloraxis",
                    zsmooth="best",
                )
            ]
        )
        fig.update_layout(
            title=title,
            xaxis_title="Residue",
            yaxis_title="Delta",
            font_family=self.plot_config.font_family,
            font_size=self.plot_config.font_size,
            template=self.plot_config.template,
        )
        return fig

    def save_plot(self, fig: go.Figure, output_path: Path, name: str,
                  fallback: Optional[Dict] = None) -> None:
        output_path.mkdir(parents=True, exist_ok=True)

        for fmt in self.plot_config.save_formats:
            file_path = output_path / f"{name}.{fmt}"
            try:
                if fmt == "html":
                    fig.write_html(str(file_path))
                elif fmt in ("svg", "png"):
                    write_static(
                        fig, file_path, fmt, self.plot_config.width, self.plot_config.height, self.plot_config.scale
                    )
                logger.info(f"Saved {fmt} plot: {file_path}")
            except Exception as exc:
                if fmt == "svg" and fallback:
                    if self._save_svg_fallback(file_path, fallback):
                        logger.info(f"Saved SVG via matplotlib fallback: {file_path}")
                        continue
                if isinstance(exc, StaticExportUnavailable):
                    logger.debug("Static export unavailable; skipped %s", file_path)
                    continue
                logger.error(f"Error saving {fmt} plot {file_path}: {exc}")

    def _save_svg_fallback(self, file_path: Path, fallback: Dict) -> bool:
        kind = fallback.get("kind")
        try:
            if kind == "compare":
                return self._save_compare_svg(
                    file_path,
                    fallback.get("time"),
                    fallback.get("series"),
                    fallback.get("color_map"),
                    fallback.get("title"),
                    fallback.get("show_std", False),
                )
        except Exception as exc:
            logger.error(f"Matplotlib fallback failed for {file_path}: {exc}")
        return False

    def _save_compare_svg(
        self,
        file_path: Path,
        time: np.ndarray,
        series: Dict[str, Dict[str, np.ndarray]],
        color_map: Dict[str, Tuple[int, int, int]],
        title: str,
        show_std: bool,
    ) -> bool:
        if time is None or not series:
            return False

        _, plt = init_matplotlib()
        fig, ax = plt.subplots(
            figsize=(self.plot_config.width / 100, self.plot_config.height / 100),
            dpi=100 * self.plot_config.scale,
        )
        plt.rcParams['font.family'] = self.plot_config.font_family
        plt.rcParams['font.size'] = self.plot_config.font_size

        for ligand_name, stats in series.items():
            mean = stats["mean"]
            std = stats["std"]
            color = color_map.get(ligand_name)
            mpl_color = None
            if color:
                mpl_color = (color[0] / 255, color[1] / 255, color[2] / 255)
            ax.plot(time, mean, linewidth=2, label=ligand_name, color=mpl_color)
            if show_std:
                ax.fill_between(time, mean - std, mean + std, alpha=0.2, color=mpl_color)

        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Distance (A)")
        ax.set_title(title)
        ax.grid(True, alpha=0.2)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, frameon=False)
        fig.subplots_adjust(bottom=0.25)
        fig.savefig(file_path, format="svg", bbox_inches="tight")
        plt.close(fig)
        return True
