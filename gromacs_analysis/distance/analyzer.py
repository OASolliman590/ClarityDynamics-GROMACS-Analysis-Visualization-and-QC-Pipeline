"""
Distance analysis orchestrator.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional, List

import numpy as np
import pandas as pd

from .config import DistanceConfig
from .processor import DistanceProcessor
from .plotter import DistancePlotter

logger = logging.getLogger(__name__)


def _file_safe(label: str) -> str:
    """Residue labels such as 'A:GLY4' as portable file names ('A_GLY4')."""
    return str(label).replace(":", "_").replace("/", "_")


class DistanceAnalyzer:
    """Run ligand-residue distance analysis."""

    def __init__(self, config: DistanceConfig):
        self.config = config
        self.processor = DistanceProcessor(config)
        self.plotter = DistancePlotter(config)

    def run_analysis(self) -> Dict:
        logger.info("=" * 60)
        logger.info("Starting Distance Analysis")
        logger.info("=" * 60)

        self._create_output_dirs()

        data = self.processor.load_data()
        time = data["time"]
        residues = data["residues"]
        overall = data.get("overall")

        plots_created = []
        plots_dir = self.config.output_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)

        if self.config.per_residue_plots:
            for residue in residues:
                fig = self.plotter.create_per_residue_plot(time, residue)
                name = f"distance_{_file_safe(residue['label'])}"
                self.plotter.save_plot(
                    fig,
                    plots_dir,
                    name,
                    fallback={"kind": "single", "time": time, "residue": residue},
                )
                plots_created.append(name)

        if self.config.combined_plot:
            fig = self.plotter.create_combined_plot(time, residues)
            name = "distance_all_residues"
            self.plotter.save_plot(
                fig,
                plots_dir,
                name,
                fallback={"kind": "combined", "time": time, "residues": residues},
            )
            plots_created.append(name)

        data_files = self._save_data_csv(time, residues, overall)
        summary_file = self._generate_summary(residues)

        logger.info("=" * 60)
        logger.info("Distance Analysis Complete")
        logger.info("=" * 60)
        logger.info(f"Output directory: {self.config.output_dir}")
        logger.info(f"Plots created: {len(plots_created)}")

        return {
            "success": True,
            "output_dir": str(self.config.output_dir),
            "plots_dir": str(plots_dir),
            "plots_created": plots_created,
            "data_files": [str(p) for p in data_files],
            "summary_file": str(summary_file) if summary_file else None,
        }

    def _create_output_dirs(self) -> None:
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        (self.config.output_dir / "plots").mkdir(parents=True, exist_ok=True)
        (self.config.output_dir / "data").mkdir(parents=True, exist_ok=True)

    def _save_data_csv(self, time: np.ndarray, residues: List[Dict], overall: Optional[Dict]) -> List[Path]:
        data_dir = self.config.output_dir / "data"
        data_dir.mkdir(parents=True, exist_ok=True)

        saved = []
        for residue in residues:
            df = pd.DataFrame({
                "time_ns": time,
                "mean_distance": residue["mean"],
                "std_distance": residue["std"],
            })
            file_path = data_dir / f"distance_{_file_safe(residue['label'])}.csv"
            df.to_csv(file_path, index=False)
            saved.append(file_path)

        if overall:
            df = pd.DataFrame({
                "time_ns": time,
                "mean_distance": overall["mean"],
                "std_distance": overall["std"],
            })
            file_path = data_dir / "distance_overall.csv"
            df.to_csv(file_path, index=False)
            saved.append(file_path)

        return saved

    def _generate_summary(self, residues: List[Dict]) -> Optional[Path]:
        try:
            summary_path = self.config.output_dir / "DISTANCE_SUMMARY.md"
            with open(summary_path, "w") as f:
                f.write(f"# Distance Analysis Summary - {self.config.protein_name}\n\n")
                f.write(f"**Ligand:** {self.config.ligand_name}\n\n")
                f.write("## Residue Mean Distances\n\n")
                f.write("| Residue | Mean Distance (A) | Std (A) |\n")
                f.write("|---|---:|---:|\n")
                for residue in residues:
                    mean_val = float(np.mean(residue["mean"]))
                    std_val = float(np.mean(residue["std"]))
                    f.write(f"| {residue['label']} | {mean_val:.3f} | {std_val:.3f} |\n")
            return summary_path
        except Exception as exc:
            logger.error(f"Failed to write summary: {exc}")
            return None
