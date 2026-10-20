"""Joint PCA stage: shared-basis projections, overlap and figures."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .config import PCAJointConfig
from .joint import JointPCAResult, joint_pca, load_coordinates
from ..utils.figure_export import StaticExportUnavailable, write_static
from ..utils.palette import APO_HEX, OCCUPANCY_SCALE, categorical, mix, rgba

logger = logging.getLogger(__name__)


class PCAJointAnalyzer:
    """Run joint PCA for all configured systems and write data + figures."""

    def __init__(self, config: PCAJointConfig):
        self.config = config
        self.colors = self._assign_colors()
        self.segments: Dict[str, List[int]] = {}

    def _assign_colors(self) -> Dict[str, str]:
        colors: Dict[str, str] = {}
        slot = 0
        for system in self.config.systems:
            if system.color:
                colors[system.name] = system.color
            elif system.is_apo:
                colors[system.name] = APO_HEX
            else:
                colors[system.name] = categorical(slot)
                slot += 1
        return colors

    def _resolve(self, path: Path) -> Path:
        if path.is_absolute() or self.config.base_dir is None:
            return path
        return (self.config.base_dir / path).resolve()

    def load(self) -> Dict[str, np.ndarray]:
        """Read every replicate separately so per-replicate checks know the segment boundaries."""
        coords = {}
        for system in self.config.systems:
            per_traj_max = None
            if self.config.max_frames_per_system:
                per_traj_max = max(2, int(self.config.max_frames_per_system) // len(system.trajectories))
            blocks = [
                load_coordinates(
                    self._resolve(system.topology),
                    [self._resolve(t)],
                    selection=self.config.selection,
                    stride=self.config.stride,
                    max_frames=per_traj_max,
                )
                for t in system.trajectories
            ]
            self.segments[system.name] = [len(b) for b in blocks]
            coords[system.name] = np.concatenate(blocks)
            logger.info("Joint PCA: %s -> %d frames x %d atoms", system.name, *coords[system.name].shape[:2])
        return coords

    def run_analysis(self, coords_by_system: Optional[Dict[str, np.ndarray]] = None) -> Dict:
        out = self.config.output_dir
        data_dir, plots_dir = out / "data", out / "plots"
        data_dir.mkdir(parents=True, exist_ok=True)
        plots_dir.mkdir(parents=True, exist_ok=True)

        coords = coords_by_system if coords_by_system is not None else self.load()
        result = joint_pca(
            coords,
            n_components=self.config.n_components,
            n_overlap=self.config.n_overlap,
            fit_iterations=self.config.fit_iterations,
        )
        written = self._write_data(result, data_dir)
        written.append(self._write_convergence(result, data_dir))
        figures = []
        figures += self._save(self._projection_figure(result), plots_dir, "joint_projection_12")
        figures += self._save(self._density_figure(result), plots_dir, "joint_density_12")
        figures += self._save(self._rmsip_figure(result), plots_dir, "rmsip_heatmap")
        figures += self._save(self._scree_figure(result), plots_dir, "joint_scree")
        return {
            "success": True,
            "output_dir": str(out),
            "data_files": [str(p) for p in written],
            "figures": [str(p) for p in figures],
            "explained_variance_pc1_pc2": [float(v) for v in result.explained_variance[:2]],
        }

    # ------------------------------------------------------------------ data

    def _write_data(self, result: JointPCAResult, data_dir: Path) -> List[Path]:
        k = len(result.eigenvalues)
        pcs = [f"PC{i + 1}" for i in range(k)]
        paths = []
        eig = pd.DataFrame({
            "PC": np.arange(1, k + 1),
            "Eigenvalue": result.eigenvalues,
            "ExplainedVariance": result.explained_variance,
            "Cumulative": np.cumsum(result.explained_variance),
        })
        paths.append(data_dir / "eigenvalues.csv")
        eig.to_csv(paths[-1], index=False)
        for name, proj in result.projections.items():
            frame = pd.DataFrame(proj, columns=pcs)
            frame.insert(0, "Frame", np.arange(len(proj)))
            paths.append(data_dir / f"projections_{name}.csv")
            frame.to_csv(paths[-1], index=False)
        paths.append(data_dir / "rmsip_pairwise.csv")
        pd.DataFrame(result.rmsip_pairwise, index=result.system_names, columns=result.system_names).to_csv(paths[-1])
        paths.append(data_dir / "rmsip_to_joint.csv")
        pd.DataFrame({"System": result.system_names, "RMSIP": result.rmsip_to_joint}).to_csv(paths[-1], index=False)
        paths.append(data_dir / "eigenvectors.npy")
        np.save(paths[-1], result.eigenvectors)
        params = {
            **result.parameters,
            "selection": self.config.selection,
            "stride": self.config.stride,
            "max_frames_per_system": self.config.max_frames_per_system,
            "coordinate_unit": "Å (MDAnalysis)",
            "systems": result.system_names,
        }
        paths.append(data_dir / "parameters.json")
        paths[-1].write_text(json.dumps(params, indent=2) + "\n")
        return paths

    def _write_convergence(self, result: JointPCAResult, data_dir: Path) -> Path:
        """Cosine content of PC1-3 per replicate segment and half-split RMSIP per system."""
        from ..utils.convergence import DEFAULT_PCA_THRESHOLDS, projection_convergence  # noqa: WPS433

        limits = {**DEFAULT_PCA_THRESHOLDS, **(self.config.convergence or {})}
        rows = []
        for name, proj in result.projections.items():
            lengths = self.segments.get(name) or [len(proj)]
            if sum(lengths) != len(proj):
                lengths = [len(proj)]
            half = result.half_rmsip.get(name, float("nan"))
            start = 0
            for index, length in enumerate(lengths, start=1):
                segment = proj[start:start + length]
                start += length
                row = {"System": name, "Replicate": index, "Frames": int(length)}
                row.update(projection_convergence(segment[:, :3]))
                row["Half_RMSIP"] = half
                flags = []
                if np.isfinite(row["Cos_PC1"]) and row["Cos_PC1"] > limits["cosine_max"]:
                    flags.append("pc1_diffusive")
                if np.isfinite(half) and half < limits["half_rmsip_min"]:
                    flags.append("unstable_subspace")
                if length < 8:
                    flags.append("too_short")
                row["Flags"] = ";".join(flags)
                rows.append(row)
        path = data_dir / "convergence.csv"
        frame = pd.DataFrame(rows, columns=["System", "Replicate", "Frames", "Cos_PC1", "Cos_PC2", "Cos_PC3",
                                            "Half_RMSIP", "Flags"])
        frame.to_csv(path, index=False, float_format="%.4g")
        flagged = frame[frame["Flags"].astype(str).str.len() > 0]
        if not flagged.empty:
            logger.warning("Joint PCA sampling checks flagged %d of %d replicate(s) (see %s)",
                           len(flagged), len(frame), path)
        return path

    # ------------------------------------------------------------------ figures

    def _layout(self, fig: go.Figure, title: str, **kwargs) -> go.Figure:
        pc = self.config.plot_config
        fig.update_layout(
            title=title,
            template=pc.template,
            font_family=pc.font_family,
            font_size=pc.font_size,
            legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.18, yanchor="top"),
            **kwargs,
        )
        return fig

    def _pc_title(self, result: JointPCAResult, i: int) -> str:
        return f"PC{i + 1} ({100 * result.explained_variance[i]:.1f}% of variance, Å)"

    def _projection_figure(self, result: JointPCAResult) -> go.Figure:
        fig = go.Figure()
        for name in result.system_names:
            proj = result.projections[name]
            color = self.colors[name]
            fig.add_trace(go.Scattergl(
                x=proj[:, 0], y=proj[:, 1] if proj.shape[1] > 1 else np.zeros(len(proj)),
                mode="markers", name=name, legendgroup=name, showlegend=False,
                marker=dict(size=4, color=rgba(color, 0.35)),
                hovertemplate=f"{name}<br>frame %{{customdata}}<extra></extra>",
                customdata=np.arange(len(proj)),
            ))
        for name in result.system_names:
            proj = result.projections[name]
            fig.add_trace(go.Scatter(
                x=[proj[:, 0].mean()], y=[proj[:, 1].mean() if proj.shape[1] > 1 else 0.0],
                mode="markers", name=name, legendgroup=name, showlegend=True,
                marker=dict(size=15, symbol="diamond", color=mix(self.colors[name], "#000000", 0.25),
                            line=dict(color="white", width=2)),
            ))
        fig.update_xaxes(title=self._pc_title(result, 0))
        if len(result.eigenvalues) > 1:
            fig.update_yaxes(title=self._pc_title(result, 1))
        return self._layout(
            fig,
            f"{self.config.protein_name} · all systems on shared PC1–PC2 (joint PCA, diamonds = means)",
        )

    def _density_figure(self, result: JointPCAResult) -> go.Figure:
        fig = go.Figure()
        for name in result.system_names:
            proj = result.projections[name]
            if proj.shape[1] < 2:
                continue
            color = self.colors[name]
            fig.add_trace(go.Histogram2dContour(
                x=proj[:, 0], y=proj[:, 1], name=name, showlegend=True, showscale=False,
                ncontours=6, histnorm="probability density",
                contours=dict(coloring="none", showlabels=False),
                line=dict(color=color, width=2),
            ))
        fig.update_xaxes(title=self._pc_title(result, 0))
        if len(result.eigenvalues) > 1:
            fig.update_yaxes(title=self._pc_title(result, 1))
        return self._layout(fig, f"{self.config.protein_name} · sampling density on shared PC1–PC2")

    def _rmsip_figure(self, result: JointPCAResult) -> go.Figure:
        names = result.system_names
        k = result.parameters.get("n_overlap")
        fig = go.Figure(go.Heatmap(
            z=result.rmsip_pairwise, x=names, y=names, zmin=0, zmax=1,
            colorscale=OCCUPANCY_SCALE, xgap=2, ygap=2,
            text=np.round(result.rmsip_pairwise, 2), texttemplate="%{text:.2f}",
            colorbar=dict(title="RMSIP"),
            hovertemplate="%{y} vs %{x}: %{z:.3f}<extra></extra>",
        ))
        fig.update_yaxes(autorange="reversed")
        size = max(420, 160 + 90 * len(names))
        return self._layout(
            fig,
            f"Essential-subspace overlap<br><sup>RMSIP over the first {k} PCs · 1 = identical dominant motions</sup>",
            width=size + 160, height=size + 40,
        )

    def _scree_figure(self, result: JointPCAResult) -> go.Figure:
        k = len(result.eigenvalues)
        xs = np.arange(1, k + 1)
        color = categorical(0)
        fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08)
        fig.add_trace(go.Scatter(x=xs, y=result.eigenvalues, mode="lines+markers", name="Eigenvalue",
                                 line=dict(color=color, width=2)), row=1, col=1)
        fig.add_trace(go.Scatter(x=xs, y=100 * np.cumsum(result.explained_variance), mode="lines+markers",
                                 name="Cumulative variance", line=dict(color=color, width=2, dash="dot")),
                      row=2, col=1)
        fig.update_yaxes(title="Eigenvalue (Å²)", rangemode="tozero", row=1, col=1)
        fig.update_yaxes(title="Cumulative (%)", range=[0, 100], row=2, col=1)
        fig.update_xaxes(title="Joint eigenvector", dtick=1, row=2, col=1)
        return self._layout(fig, "Joint PCA eigenvalue spectrum", showlegend=False,
                            height=max(self.config.plot_config.height, 700))

    def _save(self, fig: go.Figure, plots_dir: Path, name: str) -> List[Path]:
        pc = self.config.plot_config
        written: List[Path] = []
        for fmt in pc.save_formats:
            path = plots_dir / f"{name}.{fmt}"
            try:
                if fmt == "html":
                    fig.write_html(str(path))
                else:
                    write_static(fig, path, fmt, pc.width, pc.height, pc.scale)
                written.append(path)
            except StaticExportUnavailable:
                logger.debug("Static export unavailable; skipped %s", path)
            except Exception as exc:  # noqa: BLE001
                logger.error("Could not save %s: %s", path, exc)
        return written


__all__ = ["PCAJointAnalyzer"]
