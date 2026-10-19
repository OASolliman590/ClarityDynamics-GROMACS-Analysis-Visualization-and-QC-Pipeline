"""
Cross-ligand PCA comparison analyzer.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio

from .config import PCACompareConfig, PCAConfig
from .processor import PCAProcessor
from .parser import PCAParser
from ..utils.mpl_fallback import init_matplotlib

logger = logging.getLogger(__name__)

_PALETTE = [
    (31, 119, 180),
    (255, 127, 14),
    (44, 160, 44),
    (214, 39, 40),
    (148, 103, 189),
    (140, 86, 75),
    (227, 119, 194),
    (127, 127, 127),
    (188, 189, 34),
    (23, 190, 207),
]


class PCAComparisonAnalyzer:
    def __init__(self, config: PCACompareConfig):
        self.config = config
        self.parser = PCAParser()
        pio.templates.default = config.plot_config.template

    def run_analysis(self) -> Dict:
        logger.info("=" * 60)
        logger.info("Starting PCA Comparison")
        logger.info("=" * 60)

        self._create_output_dirs()

        systems_data = self._load_systems()

        plots_dir = self.config.output_dir / "plots"
        data_dir = self.config.output_dir / "data"
        plots_created: List[str] = []

        # 2D projection comparisons.
        # Each system's projections are onto its *own* eigenvectors, so PC axes are
        # neither the same collective motions nor sign-consistent across systems.
        # A quantitative comparison needs a common basis: run one gmx covar on the
        # concatenated (identically fitted) trajectories and project every system
        # onto those eigenvectors, or compare subspaces via RMSIP (gmx anaeig -over).
        if self.config.compare_pairs and len(systems_data) > 1:
            logger.warning(
                "PCA projection overlays use system-specific eigenvectors; axes are not a common "
                "basis. Use a joint covariance (concatenated trajectories) for quantitative comparison."
            )
        for pair in self.config.compare_pairs:
            proj_series = self._collect_projection_series(systems_data, pair)
            if not proj_series:
                continue
            fig = self._create_projection_plot(pair, proj_series)
            name = f"pca_compare_proj_{pair}"
            self._save_plot(fig, plots_dir, name, fallback={
                "kind": "proj",
                "pair": pair,
                "series": proj_series,
            })
            plots_created.append(name)

        # Scree / cumulative variance comparison
        if self.config.compare_scree:
            scree_series = self._collect_scree_series(systems_data)
            if scree_series:
                fig = self._create_scree_plot(scree_series)
                name = "pca_compare_scree"
                self._save_plot(fig, plots_dir, name, fallback={
                    "kind": "scree",
                    "series": scree_series,
                })
                plots_created.append(name)
                self._save_scree_csv(data_dir, scree_series)

        # Cosine content comparison
        if self.config.compare_cosine:
            cosine_series = self._collect_cosine_series(systems_data)
            if cosine_series:
                fig = self._create_cosine_plot(cosine_series)
                name = "pca_compare_cosine"
                self._save_plot(fig, plots_dir, name, fallback={
                    "kind": "cosine",
                    "series": cosine_series,
                })
                plots_created.append(name)

        logger.info("=" * 60)
        logger.info("PCA Comparison Complete")
        logger.info("=" * 60)
        logger.info(f"Output directory: {self.config.output_dir}")
        logger.info(f"Plots created: {len(plots_created)}")

        return {
            "success": True,
            "output_dir": str(self.config.output_dir),
            "plots_created": plots_created,
        }

    def _create_output_dirs(self) -> None:
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        (self.config.output_dir / "plots").mkdir(parents=True, exist_ok=True)
        (self.config.output_dir / "data").mkdir(parents=True, exist_ok=True)

    def _load_systems(self) -> Dict[str, Dict]:
        systems_data: Dict[str, Dict] = {}
        for cfg in self.config.systems:
            processor = PCAProcessor(cfg)
            files = processor.auto_detect_files()
            ligand_name = cfg.ligand_name or files.get("ligand_name") or Path(cfg.base_dir).name
            systems_data[ligand_name] = {"config": cfg, "files": files}
        return systems_data

    def _collect_projection_series(self, systems_data: Dict[str, Dict], pair: str) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        series: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
        for ligand_name, data in systems_data.items():
            files = data["files"]
            proj_file = files.get("proj_2d_holo", {}).get(pair)
            if not proj_file:
                continue
            _, x, y, _, _ = self.parser.parse_xvg_single(proj_file)
            if not x or not y:
                continue
            x_arr = np.array(x, dtype=float)
            y_arr = np.array(y, dtype=float)
            x_arr, y_arr = self._downsample(x_arr, y_arr)
            series[ligand_name] = (x_arr, y_arr)
        return series

    def _collect_scree_series(self, systems_data: Dict[str, Dict]) -> Dict[str, Dict[str, np.ndarray]]:
        series: Dict[str, Dict[str, np.ndarray]] = {}
        for ligand_name, data in systems_data.items():
            eigen_file = data["files"].get("eigenvals")
            if not eigen_file:
                continue
            _, pcs, eigenvals, _, _ = self.parser.parse_xvg_single(eigen_file)
            if not pcs or not eigenvals:
                continue
            vals = np.array(eigenvals, dtype=float)
            if np.sum(vals) == 0:
                continue
            frac = vals / np.sum(vals)
            cum = np.cumsum(frac)
            series[ligand_name] = {
                "pcs": np.array(pcs, dtype=float),
                "eigenvals": vals,
                "fraction": frac,
                "cumulative": cum,
            }
        return series

    def _collect_cosine_series(self, systems_data: Dict[str, Dict]) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        series: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
        for ligand_name, data in systems_data.items():
            cc_file = data["files"].get("cc_holo")
            if not cc_file:
                continue
            _, pcs, values, _, _ = self.parser.parse_xvg_single(cc_file)
            if not pcs or not values:
                continue
            series[ligand_name] = (np.array(pcs, dtype=float), np.array(values, dtype=float))
        return series

    def _downsample(self, x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        if self.config.downsample and self.config.downsample > 1:
            return x[:: self.config.downsample], y[:: self.config.downsample]
        if self.config.max_points and len(x) > self.config.max_points:
            step = max(1, len(x) // self.config.max_points)
            return x[::step], y[::step]
        return x, y

    def _color_for_index(self, idx: int) -> Tuple[int, int, int]:
        return _PALETTE[idx % len(_PALETTE)]

    def _create_projection_plot(self, pair: str, series: Dict[str, Tuple[np.ndarray, np.ndarray]]) -> go.Figure:
        fig = go.Figure()
        for idx, (ligand_name, (x, y)) in enumerate(series.items()):
            color = self._color_for_index(idx)
            fig.add_trace(go.Scatter(
                x=x,
                y=y,
                mode="markers",
                name=ligand_name,
                marker=dict(
                    size=3,
                    opacity=0.55,
                    color=f"rgb({color[0]},{color[1]},{color[2]})"
                )
            ))

        fig.add_annotation(
            text=(
                "Each system projected on its own eigenvectors; axes are not a common basis "
                "(qualitative comparison only)."
            ),
            xref="paper",
            yref="paper",
            x=0.5,
            y=1.02,
            showarrow=False,
            font=dict(size=10),
        )
        fig.update_layout(
            title=f"PCA Projection PC{pair[0]} vs PC{pair[1]} · {self.config.protein_name}",
            xaxis_title=f"PC{pair[0]}",
            yaxis_title=f"PC{pair[1]}",
            font_family=self.config.plot_config.font_family,
            font_size=self.config.plot_config.font_size,
            template=self.config.plot_config.template,
            legend=dict(
                orientation='h',
                y=-0.35,
                x=0.5,
                xanchor='center'
            )
        )
        return fig

    def _create_scree_plot(self, series: Dict[str, Dict[str, np.ndarray]]) -> go.Figure:
        fig = go.Figure()
        for idx, (ligand_name, data) in enumerate(series.items()):
            color = self._color_for_index(idx)
            fig.add_trace(go.Scatter(
                x=data["pcs"],
                y=data["cumulative"],
                mode="lines+markers",
                name=ligand_name,
                line=dict(color=f"rgb({color[0]},{color[1]},{color[2]})", width=2),
            ))

        fig.update_layout(
            title=f"Cumulative Variance Explained · {self.config.protein_name}",
            xaxis_title="PC Index",
            yaxis_title="Cumulative Variance",
            font_family=self.config.plot_config.font_family,
            font_size=self.config.plot_config.font_size,
            template=self.config.plot_config.template,
            legend=dict(
                orientation='h',
                y=-0.35,
                x=0.5,
                xanchor='center'
            )
        )
        return fig

    def _create_cosine_plot(self, series: Dict[str, Tuple[np.ndarray, np.ndarray]]) -> go.Figure:
        fig = go.Figure()
        for idx, (ligand_name, (pcs, values)) in enumerate(series.items()):
            color = self._color_for_index(idx)
            fig.add_trace(go.Scatter(
                x=pcs,
                y=values,
                mode="lines+markers",
                name=ligand_name,
                line=dict(color=f"rgb({color[0]},{color[1]},{color[2]})", width=2),
            ))

        fig.update_layout(
            title=f"Cosine Content · {self.config.protein_name}",
            xaxis_title="PC Index",
            yaxis_title="Cosine Content",
            font_family=self.config.plot_config.font_family,
            font_size=self.config.plot_config.font_size,
            template=self.config.plot_config.template,
            legend=dict(
                orientation='h',
                y=-0.35,
                x=0.5,
                xanchor='center'
            )
        )
        return fig

    def _save_scree_csv(self, data_dir: Path, series: Dict[str, Dict[str, np.ndarray]]) -> None:
        for ligand_name, data in series.items():
            df = pd.DataFrame({
                "pc_index": data["pcs"],
                "eigenvalue": data["eigenvals"],
                "fraction": data["fraction"],
                "cumulative": data["cumulative"],
            })
            safe_name = ligand_name.replace(" ", "_")
            df.to_csv(data_dir / f"scree_{safe_name}.csv", index=False)

    def _save_plot(self, fig: go.Figure, output_path: Path, name: str,
                   fallback: Optional[Dict] = None) -> None:
        output_path.mkdir(parents=True, exist_ok=True)

        for fmt in self.config.plot_config.save_formats:
            file_path = output_path / f"{name}.{fmt}"
            try:
                if fmt == "html":
                    fig.write_html(str(file_path))
                elif fmt in ("svg", "png"):
                    fig.write_image(
                        str(file_path),
                        format=fmt,
                        width=self.config.plot_config.width * self.config.plot_config.scale,
                        height=self.config.plot_config.height * self.config.plot_config.scale,
                    )
                logger.info(f"Saved {fmt} plot: {file_path}")
            except Exception as exc:
                if fmt == "svg" and fallback:
                    if self._save_svg_fallback(file_path, fallback):
                        logger.info(f"Saved SVG via matplotlib fallback: {file_path}")
                        continue
                logger.error(f"Error saving {fmt} plot {file_path}: {exc}")

    def _save_svg_fallback(self, file_path: Path, fallback: Dict) -> bool:
        kind = fallback.get("kind")
        try:
            if kind == "proj":
                return self._save_projection_svg(
                    file_path,
                    fallback.get("pair"),
                    fallback.get("series", {}),
                )
            if kind == "scree":
                return self._save_scree_svg(
                    file_path,
                    fallback.get("series", {}),
                )
            if kind == "cosine":
                return self._save_cosine_svg(
                    file_path,
                    fallback.get("series", {}),
                )
        except Exception as exc:
            logger.error(f"Matplotlib fallback failed for {file_path}: {exc}")
        return False

    def _save_projection_svg(self, file_path: Path, pair: str,
                             series: Dict[str, Tuple[np.ndarray, np.ndarray]]) -> bool:
        if not series:
            return False
        _, plt = init_matplotlib()
        fig, ax = plt.subplots(
            figsize=(self.config.plot_config.width / 100, self.config.plot_config.height / 100),
            dpi=100 * self.config.plot_config.scale,
        )
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size

        for idx, (ligand_name, (x, y)) in enumerate(series.items()):
            color = self._color_for_index(idx)
            ax.scatter(x, y, s=5, alpha=0.55, label=ligand_name, color=(color[0]/255, color[1]/255, color[2]/255))

        ax.set_xlabel(f"PC{pair[0]}")
        ax.set_ylabel(f"PC{pair[1]}")
        ax.set_title(f"PCA Projection PC{pair[0]} vs PC{pair[1]} · {self.config.protein_name}")
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, frameon=False)
        fig.subplots_adjust(bottom=0.25)
        fig.savefig(file_path, format="svg", bbox_inches="tight")
        plt.close(fig)
        return True

    def _save_scree_svg(self, file_path: Path, series: Dict[str, Dict[str, np.ndarray]]) -> bool:
        if not series:
            return False
        _, plt = init_matplotlib()
        fig, ax = plt.subplots(
            figsize=(self.config.plot_config.width / 100, self.config.plot_config.height / 100),
            dpi=100 * self.config.plot_config.scale,
        )
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size

        for idx, (ligand_name, data) in enumerate(series.items()):
            color = self._color_for_index(idx)
            ax.plot(data["pcs"], data["cumulative"], marker='o', linewidth=2, label=ligand_name,
                    color=(color[0]/255, color[1]/255, color[2]/255))

        ax.set_xlabel("PC Index")
        ax.set_ylabel("Cumulative Variance")
        ax.set_title(f"Cumulative Variance Explained · {self.config.protein_name}")
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, frameon=False)
        fig.subplots_adjust(bottom=0.25)
        fig.savefig(file_path, format="svg", bbox_inches="tight")
        plt.close(fig)
        return True

    def _save_cosine_svg(self, file_path: Path, series: Dict[str, Tuple[np.ndarray, np.ndarray]]) -> bool:
        if not series:
            return False
        _, plt = init_matplotlib()
        fig, ax = plt.subplots(
            figsize=(self.config.plot_config.width / 100, self.config.plot_config.height / 100),
            dpi=100 * self.config.plot_config.scale,
        )
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size

        for idx, (ligand_name, (pcs, values)) in enumerate(series.items()):
            color = self._color_for_index(idx)
            ax.plot(pcs, values, marker='o', linewidth=2, label=ligand_name,
                    color=(color[0]/255, color[1]/255, color[2]/255))

        ax.set_xlabel("PC Index")
        ax.set_ylabel("Cosine Content")
        ax.set_title(f"Cosine Content · {self.config.protein_name}")
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, frameon=False)
        fig.subplots_adjust(bottom=0.25)
        fig.savefig(file_path, format="svg", bbox_inches="tight")
        plt.close(fig)
        return True
