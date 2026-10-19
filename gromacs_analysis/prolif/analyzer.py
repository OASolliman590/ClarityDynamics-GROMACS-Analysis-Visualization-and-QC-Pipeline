"""ProLIF analysis orchestration."""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from .config import ProlifConfig, ProlifCompareConfig
from .processor import ProlifProcessor
from .plotter import ProlifPlotter
from ..utils.figure_export import scaled_height, scaled_width
from ..utils.palette import OCCUPANCY_SCALE, categorical

logger = logging.getLogger(__name__)


class ProlifAnalyzer:
    def __init__(self, config: ProlifConfig):
        self.config = config
        self.processor = ProlifProcessor(config)
        self.plotter = ProlifPlotter(config.plot_config, config.protein_name)

    def run_analysis(self) -> Dict:
        logger.info("=" * 60)
        logger.info("Starting ProLIF Interaction Fingerprinting")
        logger.info("=" * 60)

        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        data_dir = self.config.output_dir / "data"
        plots_dir = self.config.output_dir / "plots"
        data_dir.mkdir(parents=True, exist_ok=True)
        plots_dir.mkdir(parents=True, exist_ok=True)

        results = self.processor.run()
        replicate_outputs = results.get("replicates", [])
        if not replicate_outputs:
            raise ValueError("No ProLIF replicates were processed.")

        occupancy_by_replicate: List[pd.DataFrame] = []
        for rep in replicate_outputs:
            rep_idx = rep["replicate"]
            fp = rep["fingerprint"]
            rep_df = fp.to_dataframe()
            rep_df.to_csv(data_dir / f"fingerprint_rep{rep_idx}.csv")

            if self.config.save_pickle:
                fp.to_pickle(str(data_dir / f"fingerprint_rep{rep_idx}.pkl"))
                if rep_idx == 1:
                    fp.to_pickle(str(data_dir / "fingerprint.pkl"))

            rep_occupancy = self._compute_occupancy(rep_df, self.config.interaction_types)
            rep_occupancy["Replicate"] = rep_idx
            rep_occupancy["Trajectory"] = rep["trajectory"]
            occupancy_by_replicate.append(rep_occupancy)

            if "barcode" in self.config.outputs:
                self._save_barcode_outputs(fp, rep_df, rep_idx, plots_dir, data_dir)

        replicate_occ = self._combine_replicate_occupancies(occupancy_by_replicate)
        replicate_occ.to_csv(data_dir / "occupancy_replicates.csv", index=False)

        summary_occ = self._summarize_replicate_occupancy(replicate_occ)
        summary_occ.to_csv(data_dir / "occupancy_summary.csv", index=False)
        if "occupancy" in self.config.outputs:
            self._save_occupancy_figures(summary_occ, plots_dir)

        # Backward-compatible occupancy output, now with mean/SD/replicate columns.
        summary_occ.rename(
            columns={
                "OccupancyMean": "Occupancy",
                "OccupancySD": "Occupancy_SD",
                "ReplicateCount": "N_replicates",
            }
        ).to_csv(data_dir / "occupancy.csv", index=False)

        first_rep = replicate_outputs[0]
        first_fp = first_rep["fingerprint"]
        first_df = first_fp.to_dataframe()
        first_df.to_csv(data_dir / "fingerprint.csv")

        if "lignetwork" in self.config.outputs:
            try:
                fp_for_lignetwork = first_rep.get("fingerprint_count") or first_fp
                if self.config.lignetwork_kind == "frame":
                    fig = fp_for_lignetwork.plot_lignetwork(
                        first_rep["ligand_mol"],
                        kind="frame",
                        frame=0,
                        display_all=self.config.lignetwork_display_all,
                    )
                else:
                    # Interactions present in >= threshold of frames, not just one snapshot.
                    fig = fp_for_lignetwork.plot_lignetwork(
                        first_rep["ligand_mol"],
                        kind="aggregate",
                        threshold=self.config.lignetwork_threshold,
                        display_all=self.config.lignetwork_display_all,
                    )
                self._save_plot(fig, plots_dir, "lignetwork")
            except Exception as exc:
                logger.warning(f"Failed to generate lignetwork plot: {exc}")

        if "tanimoto" in self.config.outputs:
            self._save_tanimoto(first_fp, first_df, plots_dir)

        self._write_chemistry_requirements(
            data_dir / "chemistry_requirements.md",
            results.get("chemistry_requirements", {}),
        )
        chemistry_summary = results.get("chemistry_summary", {})
        self._write_chemistry_summary(
            data_dir / "chemistry_summary.csv",
            [rep.get("chemistry", {}) for rep in replicate_outputs],
        )
        if chemistry_summary:
            pd.DataFrame([chemistry_summary]).to_csv(data_dir / "chemistry_overview.csv", index=False)

        return {
            "success": True,
            "output_dir": str(self.config.output_dir),
            "occupancy_file": str(data_dir / "occupancy.csv"),
            "occupancy_replicates_file": str(data_dir / "occupancy_replicates.csv"),
            "occupancy_summary_file": str(data_dir / "occupancy_summary.csv"),
        }

    def _save_plot(self, fig, plots_dir: Path, name: str) -> None:
        self.plotter.save_figure(fig, plots_dir, name)

    def _save_tanimoto(self, fp, df: pd.DataFrame, plots_dir: Path) -> None:
        try:
            from rdkit import DataStructs  # noqa: WPS433
        except ImportError:
            logger.warning("RDKit not installed; skipping tanimoto heatmap")
            return

        bitvectors = fp.to_bitvectors()
        if not bitvectors:
            return

        labels = list(df.index) if df.index is not None else list(range(len(bitvectors)))
        matrix = []
        for bv in bitvectors:
            matrix.append(DataStructs.BulkTanimotoSimilarity(bv, bitvectors))
        matrix = np.array(matrix, dtype=float)

        self.plotter.save_tanimoto_heatmap(matrix, labels, plots_dir, "tanimoto_similarity")

    def _compute_occupancy(self, df: pd.DataFrame, interaction_types: List[str]) -> pd.DataFrame:
        if df.empty:
            return pd.DataFrame(columns=["Interaction", "Residue", "Occupancy"])

        interactions_lower = {name.lower(): name for name in interaction_types}

        rows = []
        for col in df.columns:
            interaction, residue = self._parse_column(col, interactions_lower)
            if not residue:
                continue
            occ = float(df[col].mean()) * 100.0
            rows.append({
                "Interaction": interaction or "Unknown",
                "Residue": residue,
                "Occupancy": occ,
            })

        if not rows:
            return pd.DataFrame(columns=["Interaction", "Residue", "Occupancy"])

        out = pd.DataFrame(rows)
        return out.sort_values(["Interaction", "Residue"]).reset_index(drop=True)

    @staticmethod
    def _combine_replicate_occupancies(occupancies: List[pd.DataFrame]) -> pd.DataFrame:
        if not occupancies:
            return pd.DataFrame(columns=["Interaction", "Residue", "Occupancy", "Replicate", "Trajectory"])
        return pd.concat(occupancies, ignore_index=True)

    @staticmethod
    def _summarize_replicate_occupancy(occupancy_df: pd.DataFrame) -> pd.DataFrame:
        if occupancy_df.empty:
            return pd.DataFrame(
                columns=["Interaction", "Residue", "OccupancyMean", "OccupancySD", "ReplicateCount"]
            )
        out = (
            occupancy_df
            .groupby(["Interaction", "Residue"], as_index=False)["Occupancy"]
            .agg(["mean", "std", "count"])
            .reset_index()
            .rename(
                columns={
                    "mean": "OccupancyMean",
                    "std": "OccupancySD",
                    "count": "ReplicateCount",
                }
            )
        )
        out["OccupancySD"] = out["OccupancySD"].fillna(0.0)
        return out.sort_values(["Interaction", "Residue"]).reset_index(drop=True)

    @staticmethod
    def _write_chemistry_requirements(path: Path, requirements: Dict) -> None:
        lines = ["# ProLIF Ligand Chemistry Requirements", ""]
        if not requirements:
            lines.append("- No requirements metadata was provided.")
        else:
            for key, value in requirements.items():
                lines.append(f"- `{key}`: {value}")
        path.write_text("\n".join(lines) + "\n")

    @staticmethod
    def _write_chemistry_summary(path: Path, chemistry_rows: List[Dict]) -> None:
        if not chemistry_rows:
            pd.DataFrame(columns=["replicate"]).to_csv(path, index=False)
            return
        pd.DataFrame(chemistry_rows).to_csv(path, index=False)

    def _parse_column(self, col, interactions_lower: Dict[str, str]) -> Tuple[str, str]:
        # ProLIF columns are (ligand residue, protein residue, interaction).
        if isinstance(col, tuple) and len(col) >= 3:
            interaction = str(col[-1])
            return interactions_lower.get(interaction.lower(), interaction), str(col[1])
        items = list(col) if isinstance(col, tuple) else [col]

        interaction = None
        residue = None
        for item in items:
            label = str(item)
            key = label.lower()
            if key in interactions_lower and interaction is None:
                interaction = interactions_lower[key]
                continue
            if residue is None and self._looks_like_residue(label):
                residue = label

        if residue is None:
            for item in items:
                label = str(item)
                if self._looks_like_residue(label):
                    residue = label
                    break

        return interaction or "", residue or ""

    @staticmethod
    def _looks_like_residue(label: str) -> bool:
        if "." in label:
            label = label.split(".")[0]
        return bool(re.match(r"^[A-Z0-9]{1,4}?-?\d+", label))

    def _barcode_plot_kwargs(self, n_rows: int = 0) -> Dict:
        width, height = self.config.barcode_figsize
        # Give every interaction row ~0.22 in so residue labels never overlap.
        height = max(height, 1.5 + 0.22 * n_rows)
        return {
            "n_frame_ticks": self.config.barcode_n_frame_ticks,
            "residues_tick_location": self.config.barcode_residues_tick_location,
            "xlabel": self.config.barcode_xlabel,
            "figsize": (width, min(height, 60.0)),
            "dpi": self.config.barcode_dpi,
        }

    def _save_occupancy_figures(self, summary: pd.DataFrame, plots_dir: Path) -> None:
        """Residue x interaction occupancy heatmap and top-residue bar chart (mean ± SD over replicates)."""
        if summary.empty:
            return
        name = f"{self.config.protein_name} · {self.config.ligand_name}"
        font = self.config.plot_config.font_size

        table = summary.pivot_table(index="Residue", columns="Interaction", values="OccupancyMean", aggfunc="max")
        sd = summary.pivot_table(index="Residue", columns="Interaction", values="OccupancySD", aggfunc="max")
        order = table.max(axis=1).sort_values(ascending=False).index
        table, sd = table.loc[order].fillna(0.0), sd.reindex(index=order, columns=table.columns).fillna(0.0)

        heat = go.Figure(go.Heatmap(
            z=table.values,
            x=list(table.columns),
            y=list(table.index),
            customdata=sd.values,
            zmin=0,
            zmax=100,
            colorscale=OCCUPANCY_SCALE,
            xgap=2,
            ygap=2,
            colorbar=dict(title="Occupancy (%)"),
            hovertemplate="%{y} · %{x}<br>%{z:.1f} ± %{customdata:.1f} %<extra></extra>",
        ))
        heat.update_yaxes(autorange="reversed")
        heat.update_layout(
            title=f"{name} · interaction occupancy",
            xaxis_title="Interaction",
            yaxis_title="Residue",
            template=self.config.plot_config.template,
            font_family=self.config.plot_config.font_family,
            font_size=font,
            width=scaled_width(len(table.columns), per_col=110, base=360, minimum=700, maximum=2200),
            height=scaled_height(len(table.index), per_row=max(22, font + 6)),
            margin=dict(l=140, r=40, t=90, b=90),
        )
        self.plotter.save_figure(heat, plots_dir, "occupancy_heatmap")

        per_residue = (
            summary.sort_values("OccupancyMean", ascending=False)
            .drop_duplicates("Residue")
            .head(self.config.occupancy_top_n)
            .iloc[::-1]
        )
        bars = go.Figure(go.Bar(
            x=per_residue["OccupancyMean"],
            y=per_residue["Residue"],
            orientation="h",
            error_x=dict(type="data", array=per_residue["OccupancySD"], visible=True, thickness=1.5,
                         color="rgba(0,0,0,0.45)"),
            marker=dict(color=categorical(0), line=dict(color="white", width=2)),
            customdata=per_residue["Interaction"],
            hovertemplate="%{y}: %{x:.1f}% (strongest: %{customdata})<extra></extra>",
        ))
        bars.update_layout(
            title=f"{name} · top {len(per_residue)} interacting residues",
            xaxis=dict(title="Occupancy (% of frames, mean ± SD over replicates)", range=[0, 110]),
            yaxis_title="Residue",
            template=self.config.plot_config.template,
            font_family=self.config.plot_config.font_family,
            font_size=font,
            width=max(900, self.config.plot_config.width),
            height=scaled_height(len(per_residue), per_row=max(24, font + 8)),
            margin=dict(l=140, r=60, t=90, b=90),
            showlegend=False,
        )
        self.plotter.save_figure(bars, plots_dir, "occupancy_top_residues")

    def _save_barcode_outputs(
        self,
        fp,
        rep_df: pd.DataFrame,
        rep_idx: int,
        plots_dir: Path,
        data_dir: Path,
    ) -> None:
        # 1) Keep the native ProLIF barcode as static figure for compatibility.
        barcode_fig = fp.plot_barcode(**self._barcode_plot_kwargs(len(rep_df.columns)))
        self.plotter.save_figure(barcode_fig, plots_dir, f"barcode_rep{rep_idx}", formats=["svg"])

        # 2) Generate readable interactive barcode HTML with residue-level filtering.
        residue_binary, qc = self._build_residue_frame_barcode(rep_df)
        qc.to_csv(data_dir / f"barcode_qc_rep{rep_idx}.csv", index=False)
        html_fig = self._make_interactive_barcode_figure(residue_binary, rep_idx)
        # Save canonical HTML name and keep a legacy alias for compatibility with older references.
        self.plotter.save_figure(html_fig, plots_dir, f"barcode_rep{rep_idx}", formats=["html"])
        self.plotter.save_figure(html_fig, plots_dir, f"barcode_rep{rep_idx}_interactive", formats=["html"])

    def _build_residue_frame_barcode(self, rep_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        if rep_df.empty:
            return (
                pd.DataFrame(index=rep_df.index),
                pd.DataFrame(columns=["Residue", "InteractingFrames", "OccupancyPercent"]),
            )

        binary_df = rep_df.astype(np.uint8)
        residue_map: Dict = {}
        for col in binary_df.columns:
            residue = self._extract_protein_residue_from_column(col)
            if not residue:
                continue
            residue_map[col] = residue
        if not residue_map:
            return (
                pd.DataFrame(index=rep_df.index),
                pd.DataFrame(columns=["Residue", "InteractingFrames", "OccupancyPercent"]),
            )

        grouped = binary_df[list(residue_map.keys())].T
        grouped["Residue"] = [residue_map[c] for c in grouped.index]
        residue_binary = grouped.groupby("Residue", sort=True).max().T

        interacting_frames = residue_binary.sum(axis=0).astype(int)
        occupancy_percent = (interacting_frames / max(1, len(residue_binary.index)) * 100.0).round(3)
        qc = pd.DataFrame(
            {
                "Residue": interacting_frames.index,
                "InteractingFrames": interacting_frames.values,
                "OccupancyPercent": occupancy_percent.values,
            }
        ).sort_values(["InteractingFrames", "Residue"], ascending=[False, True]).reset_index(drop=True)

        if self.config.barcode_only_interacting_residues:
            keep = set(qc[qc["InteractingFrames"] >= self.config.barcode_min_interaction_frames]["Residue"])
            residue_binary = residue_binary[[c for c in residue_binary.columns if c in keep]]
            qc = qc[qc["Residue"].isin(keep)].reset_index(drop=True)

        if self.config.barcode_max_residues is not None and len(qc) > self.config.barcode_max_residues:
            keep = set(qc.head(self.config.barcode_max_residues)["Residue"])
            residue_binary = residue_binary[[c for c in residue_binary.columns if c in keep]]
            qc = qc.head(self.config.barcode_max_residues).reset_index(drop=True)

        return residue_binary, qc

    def _make_interactive_barcode_figure(self, residue_binary: pd.DataFrame, rep_idx: int) -> go.Figure:
        if residue_binary.empty:
            fig = go.Figure()
            fig.add_annotation(
                text="No interacting residues after barcode QC filters.",
                x=0.5,
                y=0.5,
                xref="paper",
                yref="paper",
                showarrow=False,
            )
            fig.update_layout(
                title=f"{self.config.protein_name} · {self.config.ligand_name} · Barcode Rep{rep_idx}",
                template=self.config.plot_config.template,
                font_family=self.config.plot_config.font_family,
                font_size=self.config.plot_config.font_size,
                height=420,
                width=1200,
            )
            return fig

        residues = list(residue_binary.columns)
        x_labels = [str(i) for i in residue_binary.index.tolist()]
        z = residue_binary.T.values
        height = min(2200, max(460, 24 * len(residues) + 220))
        width = min(2600, max(1200, 6 * len(x_labels) + 300))

        fig = go.Figure(
            data=[
                go.Heatmap(
                    z=z,
                    x=x_labels,
                    y=residues,
                    zmin=0,
                    zmax=1,
                    colorscale=[[0.0, "#F1F3F5"], [1.0, "#2563EB"]],
                    colorbar={"title": "Interaction"},
                    hovertemplate="Frame %{x}<br>Residue %{y}<br>Interaction %{z}<extra></extra>",
                )
            ]
        )
        fig.update_layout(
            title=f"{self.config.protein_name} · {self.config.ligand_name} · Barcode Rep{rep_idx} (Interacting Residues)",
            template=self.config.plot_config.template,
            font_family=self.config.plot_config.font_family,
            font_size=self.config.plot_config.font_size,
            xaxis_title=self.config.barcode_xlabel,
            yaxis_title="Residue",
            height=height,
            width=width,
            margin={"l": 120, "r": 40, "t": 90, "b": 80},
            xaxis={
                "tickmode": "auto",
                "nticks": max(10, min(30, self.config.barcode_n_frame_ticks)),
                "showgrid": False,
            },
            yaxis={
                "tickfont": {"size": max(10, self.config.plot_config.font_size - 3)},
            },
        )
        fig.update_yaxes(autorange="reversed")
        return fig

    def _extract_protein_residue_from_column(self, col) -> str:
        if isinstance(col, tuple) and len(col) >= 3:
            return str(col[1]).split(".")[0]
        if isinstance(col, tuple):
            for item in col:
                label = str(item)
                if self._looks_like_residue(label):
                    return label.split(".")[0]
            if len(col) >= 2 and str(col[1]):
                return str(col[1]).split(".")[0]
        label = str(col)
        if self._looks_like_residue(label):
            return label.split(".")[0]
        return ""


class ProlifComparisonAnalyzer:
    def __init__(self, config: ProlifCompareConfig):
        self.config = config
        self.plotter = ProlifPlotter(config.plot_config, config.protein_name)

    def run_analysis(self) -> Dict:
        logger.info("=" * 60)
        logger.info("Starting ProLIF Comparison")
        logger.info("=" * 60)

        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        data_dir = self.config.output_dir / "data"
        plots_dir = self.config.output_dir / "plots"
        data_dir.mkdir(parents=True, exist_ok=True)
        plots_dir.mkdir(parents=True, exist_ok=True)

        systems_data: Dict[str, pd.DataFrame] = {}
        for cfg in self.config.systems:
            occ_summary = cfg.output_dir / "data" / "occupancy_summary.csv"
            occ_path = occ_summary if occ_summary.exists() else (cfg.output_dir / "data" / "occupancy.csv")
            if not occ_path.exists():
                logger.warning(f"Missing occupancy file for {cfg.ligand_name}: {occ_path}")
                continue
            systems_data[cfg.ligand_name] = pd.read_csv(occ_path)

        if not systems_data:
            raise ValueError("No ProLIF occupancy files found for comparison")

        summary_df = self._build_compare_table(systems_data)
        summary_path = data_dir / "prolif_compare.csv"
        summary_df.to_csv(summary_path, index=False)
        summary_stats_path = data_dir / "prolif_compare_mean_sd.csv"
        summary_df.to_csv(summary_stats_path, index=False)

        # Heatmaps per interaction (ligands x residues)
        if summary_df.empty:
            return {
                "success": True,
                "output_dir": str(self.config.output_dir),
                "summary_file": str(summary_path),
                "summary_stats_file": str(summary_stats_path),
            }

        interactions = sorted(summary_df["Interaction"].unique())
        ligands = list(systems_data.keys())
        for interaction in interactions:
            subset = summary_df[summary_df["Interaction"] == interaction]
            residues = subset["Residue"].unique().tolist()
            if not residues:
                continue
            matrix = np.zeros((len(ligands), len(residues)))
            for i, lig in enumerate(ligands):
                lig_rows = subset[subset["Ligand"] == lig]
                occ_map = dict(zip(lig_rows["Residue"], lig_rows["Occupancy"]))
                for j, residue in enumerate(residues):
                    matrix[i, j] = occ_map.get(residue, 0.0)

            fig = self.plotter.create_heatmap(
                matrix,
                residues,
                ligands,
                title=f"{self.config.protein_name} · ProLIF {interaction} occupancy",
                colorscale=OCCUPANCY_SCALE,
                zmin=0,
                zmax=100,
            )
            fig.update_traces(xgap=2, ygap=2, colorbar=dict(title="Occupancy (%)"),
                              hovertemplate="%{y} · %{x}<br>%{z:.1f} %<extra></extra>")
            fig.update_layout(
                width=scaled_width(len(residues), per_col=34, base=320),
                height=scaled_height(len(ligands), per_row=48, base=260, minimum=380),
                xaxis=dict(tickangle=-60),
                margin=dict(l=140, r=40, t=90, b=140),
            )
            self.plotter.save_figure(fig, plots_dir, f"prolif_{interaction}_compare")

        return {
            "success": True,
            "output_dir": str(self.config.output_dir),
            "summary_file": str(summary_path),
            "summary_stats_file": str(summary_stats_path),
        }

    def _build_compare_table(self, systems_data: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        rows = []
        interaction_sets = {}
        for lig, df in systems_data.items():
            for _, row in df.iterrows():
                interaction = row.get("Interaction", "Unknown")
                residue = row.get("Residue", "")
                if not residue:
                    continue
                interaction_sets.setdefault(interaction, {}).setdefault(lig, set()).add(residue)

        for interaction, lig_map in interaction_sets.items():
            residue_sets = list(lig_map.values())
            if not residue_sets:
                continue
            if self.config.compare_mode == "intersection":
                residues = sorted(set.intersection(*residue_sets))
            else:
                residues = sorted(set.union(*residue_sets))

            for residue in residues:
                for lig, df in systems_data.items():
                    subset = df[(df["Interaction"] == interaction) & (df["Residue"] == residue)]
                    if subset.empty:
                        occupancy = 0.0
                        occupancy_sd = 0.0
                        n_replicates = 0
                    else:
                        if "OccupancyMean" in subset.columns:
                            occupancy = float(subset["OccupancyMean"].iloc[0])
                            sd_val = subset["OccupancySD"].iloc[0] if "OccupancySD" in subset.columns else 0.0
                            rep_val = subset["ReplicateCount"].iloc[0] if "ReplicateCount" in subset.columns else 1
                            occupancy_sd = float(0.0 if pd.isna(sd_val) else sd_val)
                            n_replicates = int(1 if pd.isna(rep_val) else rep_val)
                        else:
                            occupancy = float(subset["Occupancy"].iloc[0])
                            sd_val = subset["Occupancy_SD"].iloc[0] if "Occupancy_SD" in subset.columns else 0.0
                            rep_val = subset["N_replicates"].iloc[0] if "N_replicates" in subset.columns else 1
                            occupancy_sd = float(0.0 if pd.isna(sd_val) else sd_val)
                            n_replicates = int(1 if pd.isna(rep_val) else rep_val)
                    rows.append({
                        "Interaction": interaction,
                        "Residue": residue,
                        "Ligand": lig,
                        "Occupancy": occupancy,
                        "OccupancySD": occupancy_sd,
                        "ReplicateCount": n_replicates,
                    })

        return pd.DataFrame(
            rows,
            columns=[
                "Interaction",
                "Residue",
                "Ligand",
                "Occupancy",
                "OccupancySD",
                "ReplicateCount",
            ],
        )
