"""
Cross-ligand distance comparison analyzer.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

from .config import DistanceCompareConfig
from .processor import DistanceProcessor
from .plotter import DistanceComparisonPlotter

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


def _residue_sort_key(resid):
    """Sort plain numbers and chain-qualified 'A:45' identifiers together."""
    text = str(resid)
    chain, _, number = text.rpartition(":")
    try:
        return (chain, int(number))
    except ValueError:
        return (chain, 0)


class DistanceComparisonAnalyzer:
    """Compare distance profiles across multiple ligands."""

    def __init__(self, config: DistanceCompareConfig):
        self.config = config
        self.plotter = DistanceComparisonPlotter(config.plot_config, config.protein_name)

    def run_analysis(self) -> Dict:
        logger.info("=" * 60)
        logger.info("Starting Distance Comparison")
        logger.info("=" * 60)

        self._create_output_dirs()

        systems_data: Dict[str, Dict] = {}
        for dist_cfg in self.config.systems:
            if self.config.include_overall and not dist_cfg.include_overall:
                dist_cfg.include_overall = True
            data = DistanceProcessor(dist_cfg).load_data()
            systems_data[dist_cfg.ligand_name] = data

        if not systems_data:
            raise ValueError("No distance data available for comparison")

        time_ref, min_len = self._align_time(systems_data)

        residue_ids = self._resolve_residue_ids(systems_data)
        residue_labels = self._build_residue_labels(systems_data)

        color_map = self._build_color_map(list(systems_data.keys()))

        plots_dir = self.config.output_dir / "plots"
        data_dir = self.config.output_dir / "data"

        plots_created: List[str] = []
        data_files: List[Path] = []

        if self.config.per_residue_plots:
            for resid in residue_ids:
                label = residue_labels.get(resid, f"RES{resid}")
                series = self._collect_residue_series(systems_data, resid, min_len)
                if not series:
                    continue

                fig = self.plotter.create_residue_comparison_plot(
                    time_ref,
                    label,
                    series,
                    color_map,
                    show_std=self.config.show_std
                )

                name = f"distance_compare_{str(label).replace(':', '_')}"
                self.plotter.save_plot(
                    fig,
                    plots_dir,
                    name,
                    fallback={
                        "kind": "compare",
                        "time": time_ref,
                        "series": series,
                        "color_map": color_map,
                        "title": f"{self.config.protein_name} · {label} Distance Comparison",
                        "show_std": self.config.show_std,
                    },
                )
                plots_created.append(name)

                csv_path = self._save_residue_csv(data_dir, time_ref, label, series)
                data_files.append(csv_path)

        collective_series = self._collect_collective_series(systems_data, residue_ids, min_len)
        if collective_series:
            collective_label = "Collective Mean"
            fig = self.plotter.create_residue_comparison_plot(
                time_ref,
                collective_label,
                collective_series,
                color_map,
                show_std=self.config.show_std,
            )
            name = "distance_compare_collective"
            self.plotter.save_plot(
                fig,
                plots_dir,
                name,
                fallback={
                    "kind": "compare",
                    "time": time_ref,
                    "series": collective_series,
                    "color_map": color_map,
                    "title": f"{self.config.protein_name} · {collective_label} Distance Comparison",
                    "show_std": self.config.show_std,
                },
            )
            plots_created.append(name)
            data_files.append(self._save_residue_csv(data_dir, time_ref, "collective", collective_series))

        if self.config.include_overall:
            overall_series = self._collect_overall_series(systems_data, min_len)
            if overall_series:
                fig = self.plotter.create_overall_comparison_plot(
                    time_ref,
                    overall_series,
                    color_map,
                    show_std=self.config.show_std
                )
                name = "distance_compare_overall"
                self.plotter.save_plot(
                    fig,
                    plots_dir,
                    name,
                    fallback={
                        "kind": "compare",
                        "time": time_ref,
                        "series": overall_series,
                        "color_map": color_map,
                        "title": f"{self.config.protein_name} · Overall Distance Comparison",
                        "show_std": self.config.show_std,
                    },
                )
                plots_created.append(name)

                csv_path = self._save_residue_csv(data_dir, time_ref, "overall", overall_series)
                data_files.append(csv_path)

        summary_path = self._save_summary(data_dir, systems_data, residue_ids, residue_labels)

        # Heatmaps per ligand (time x residue)
        if self.config.heatmap_enabled:
            heatmap_matrices: Dict[str, np.ndarray] = {}
            heatmap_labels: Dict[str, List[str]] = {}
            for ligand_name, data in systems_data.items():
                residues = data.get("residues", [])
                if not residues:
                    continue
                residue_labels_ordered = [res["label"] for res in residues]
                matrix = self._build_heatmap_matrix_from_residues(residues, min_len)
                if matrix is None:
                    continue
                time_axis = time_ref
                if self.config.heatmap_downsample > 1:
                    step = int(self.config.heatmap_downsample)
                    matrix = matrix[::step]
                    time_axis = time_axis[::step]
                if self.config.heatmap_max_frames:
                    max_frames = int(self.config.heatmap_max_frames)
                    matrix = matrix[:max_frames]
                    time_axis = time_axis[:max_frames]

                fig = self.plotter.create_distance_heatmap(
                    time_axis,
                    residue_labels_ordered,
                    matrix,
                    title=f"{self.config.protein_name} · {ligand_name} Residue Distance Heatmap",
                )
                name = f"distance_heatmap_{ligand_name}"
                self.plotter.save_plot(
                    fig,
                    plots_dir,
                    name,
                    fallback=None,
                )
                plots_created.append(name)
                heatmap_matrices[ligand_name] = matrix
                heatmap_labels[ligand_name] = residue_labels_ordered
                data_files.append(self._save_heatmap_csv(data_dir, ligand_name, time_axis, residue_labels_ordered, matrix))

            # Delta heatmaps vs reference (use common residue set)
            if self.config.delta_enabled and heatmap_matrices:
                common_residue_ids = self._resolve_common_residue_ids(systems_data)
                if common_residue_ids:
                    residue_labels_ordered = [residue_labels.get(rid, f"RES{rid}") for rid in common_residue_ids]
                    common_matrices: Dict[str, np.ndarray] = {}
                    for ligand_name, data in systems_data.items():
                        matrix = self._build_heatmap_matrix(data, common_residue_ids, min_len)
                        if matrix is None:
                            continue
                        common_matrices[ligand_name] = matrix

                    ref_name = self.config.delta_reference or list(common_matrices.keys())[0]
                    ref_matrix = common_matrices.get(ref_name)
                    if ref_matrix is not None:
                        ref_mean = ref_matrix.mean(axis=0)
                        for ligand_name, matrix in common_matrices.items():
                            if ligand_name == ref_name:
                                continue
                            delta = matrix.mean(axis=0) - ref_mean
                            fig = self.plotter.create_delta_heatmap(
                                residue_labels_ordered,
                                delta,
                                title=f"{self.config.protein_name} · ΔDistance {ligand_name} vs {ref_name}",
                            )
                            name = f"distance_delta_{ligand_name}_vs_{ref_name}"
                            self.plotter.save_plot(fig, plots_dir, name, fallback=None)
                            plots_created.append(name)
                            data_files.append(self._save_delta_csv(data_dir, ligand_name, ref_name, residue_labels_ordered, delta))

        logger.info("=" * 60)
        logger.info("Distance Comparison Complete")
        logger.info("=" * 60)
        logger.info(f"Output directory: {self.config.output_dir}")
        logger.info(f"Plots created: {len(plots_created)}")

        return {
            "success": True,
            "output_dir": str(self.config.output_dir),
            "plots_created": plots_created,
            "data_files": [str(p) for p in data_files],
            "summary_file": str(summary_path) if summary_path else None,
        }

    def _create_output_dirs(self) -> None:
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        (self.config.output_dir / "plots").mkdir(parents=True, exist_ok=True)
        (self.config.output_dir / "data").mkdir(parents=True, exist_ok=True)

    def _align_time(self, systems_data: Dict[str, Dict]) -> Tuple[np.ndarray, int]:
        time_arrays = [data["time"] for data in systems_data.values()]
        min_len = min(len(t) for t in time_arrays)
        time_ref = time_arrays[0][:min_len]
        return time_ref, min_len

    def _resolve_residue_ids(self, systems_data: Dict[str, Dict]) -> List[int]:
        residue_sets = []
        for data in systems_data.values():
            residue_sets.append({res["resid"] for res in data["residues"]})

        if not residue_sets:
            return []

        if self.config.compare_mode == "union":
            resids = set().union(*residue_sets)
            return sorted(resids, key=_residue_sort_key)

        resids = set.intersection(*residue_sets)
        if not resids:
            logger.warning(
                "distance_compare intersection produced no common residues; "
                "falling back to union of residue sets."
            )
            resids = set().union(*residue_sets)

        return sorted(resids, key=_residue_sort_key)

    def _build_residue_labels(self, systems_data: Dict[str, Dict]) -> Dict[int, str]:
        labels: Dict[int, str] = {}
        for data in systems_data.values():
            for residue in data["residues"]:
                resid = residue["resid"]
                if resid not in labels:
                    labels[resid] = residue["label"]
        return labels

    def _resolve_common_residue_ids(self, systems_data: Dict[str, Dict]) -> List[int]:
        residue_sets = []
        for data in systems_data.values():
            residue_sets.append({res["resid"] for res in data["residues"]})

        if not residue_sets:
            return []

        resids = set.intersection(*residue_sets)
        return sorted(resids, key=_residue_sort_key)

    def _build_color_map(self, ligand_names: List[str]) -> Dict[str, Tuple[int, int, int]]:
        color_map: Dict[str, Tuple[int, int, int]] = {}
        for idx, name in enumerate(ligand_names):
            color_map[name] = _PALETTE[idx % len(_PALETTE)]
        return color_map

    def _collect_residue_series(self, systems_data: Dict[str, Dict], resid: int, min_len: int) -> Dict[str, Dict[str, np.ndarray]]:
        series: Dict[str, Dict[str, np.ndarray]] = {}
        for ligand_name, data in systems_data.items():
            for residue in data["residues"]:
                if residue["resid"] == resid:
                    series[ligand_name] = {
                        "mean": residue["mean"][:min_len],
                        "std": residue["std"][:min_len],
                    }
                    break
        return series

    def _collect_overall_series(self, systems_data: Dict[str, Dict], min_len: int) -> Dict[str, Dict[str, np.ndarray]]:
        series: Dict[str, Dict[str, np.ndarray]] = {}
        for ligand_name, data in systems_data.items():
            overall = data.get("overall")
            if overall:
                series[ligand_name] = {
                    "mean": overall["mean"][:min_len],
                    "std": overall["std"][:min_len],
                }
        return series

    def _collect_collective_series(
        self,
        systems_data: Dict[str, Dict],
        residue_ids: List[int],
        min_len: int,
    ) -> Dict[str, Dict[str, np.ndarray]]:
        series: Dict[str, Dict[str, np.ndarray]] = {}
        if not residue_ids:
            return series

        residue_id_set = set(residue_ids)
        for ligand_name, data in systems_data.items():
            rows: List[np.ndarray] = []
            for residue in data["residues"]:
                if residue["resid"] in residue_id_set:
                    rows.append(residue["mean"][:min_len])
            if rows:
                matrix = np.vstack(rows)
                series[ligand_name] = {
                    "mean": matrix.mean(axis=0),
                    "std": matrix.std(axis=0),
                }

        return series

    def _save_residue_csv(
        self,
        data_dir: Path,
        time: np.ndarray,
        label: str,
        series: Dict[str, Dict[str, np.ndarray]],
    ) -> Path:
        data: Dict[str, np.ndarray] = {"time_ns": time}
        for ligand_name, stats in series.items():
            safe_name = ligand_name.replace(" ", "_")
            data[f"{safe_name}_mean"] = stats["mean"]
            data[f"{safe_name}_std"] = stats["std"]

        df = pd.DataFrame(data)
        csv_path = data_dir / f"distance_compare_{str(label).replace(':', '_')}.csv"
        df.to_csv(csv_path, index=False)
        return csv_path

    def _save_summary(
        self,
        data_dir: Path,
        systems_data: Dict[str, Dict],
        residue_ids: List[int],
        residue_labels: Dict[int, str],
    ) -> Path:
        rows = []
        for resid in residue_ids:
            label = residue_labels.get(resid, f"RES{resid}")
            for ligand_name, data in systems_data.items():
                for residue in data["residues"]:
                    if residue["resid"] == resid:
                        rows.append({
                            "Residue": label,
                            "Ligand": ligand_name,
                            "MeanDistance": float(np.mean(residue["mean"])),
                            "StdDistance": float(np.mean(residue["std"])),
                        })
                        break
        summary_path = data_dir / "distance_compare_summary.csv"
        pd.DataFrame(rows).to_csv(summary_path, index=False)
        return summary_path

    def _build_heatmap_matrix(self, data: Dict, residue_ids: List[int], min_len: int) -> Optional[np.ndarray]:
        series = []
        for resid in residue_ids:
            found = False
            for residue in data["residues"]:
                if residue["resid"] == resid:
                    series.append(residue["mean"][:min_len])
                    found = True
                    break
            if not found:
                return None
        matrix = np.vstack(series).T
        return matrix

    def _build_heatmap_matrix_from_residues(
        self,
        residues: List[Dict],
        min_len: int,
    ) -> Optional[np.ndarray]:
        if not residues:
            return None
        series = [res["mean"][:min_len] for res in residues]
        if not series:
            return None
        return np.vstack(series).T

    def _save_heatmap_csv(
        self,
        data_dir: Path,
        ligand_name: str,
        time_axis: np.ndarray,
        residue_labels: List[str],
        matrix: np.ndarray,
    ) -> Path:
        df = pd.DataFrame(matrix, columns=residue_labels)
        df.insert(0, "time_ns", time_axis)
        path = data_dir / f"distance_heatmap_{ligand_name}.csv"
        df.to_csv(path, index=False)
        return path

    def _save_delta_csv(
        self,
        data_dir: Path,
        ligand_name: str,
        ref_name: str,
        residue_labels: List[str],
        delta: np.ndarray,
    ) -> Path:
        df = pd.DataFrame({"Residue": residue_labels, "Delta": delta})
        path = data_dir / f"distance_delta_{ligand_name}_vs_{ref_name}.csv"
        df.to_csv(path, index=False)
        return path
