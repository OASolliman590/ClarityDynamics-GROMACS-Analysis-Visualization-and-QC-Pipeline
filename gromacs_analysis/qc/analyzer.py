"""Preflight quality-control analysis."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from .config import QCConfig

logger = logging.getLogger(__name__)


class QCAnalyzer:
    """Run lightweight and optional trajectory-aware preflight checks."""

    def __init__(self, config: QCConfig):
        self.config = config

    def run_analysis(self) -> Dict:
        self.config.output_dir.mkdir(parents=True, exist_ok=True)

        report = {
            "success": True,
            "errors": [],
            "warnings": [],
            "checks": {},
            "metadata": {},
        }

        self._check_replicates(report)
        self._check_topology(report)
        self._check_index_groups(report)
        self._check_trajectories(report)

        report["success"] = len(report["errors"]) == 0
        report_path = self.config.output_dir / "qc_report.json"
        report_path.write_text(json.dumps(report, indent=2))

        if self.config.fail_on_error and report["errors"]:
            raise RuntimeError("QC preflight failed. See qc_report.json for details.")

        return {
            "success": report["success"],
            "output_dir": str(self.config.output_dir),
            "report_file": str(report_path),
            "errors": report["errors"],
            "warnings": report["warnings"],
        }

    def _check_replicates(self, report: Dict) -> None:
        actual = len(self.config.trajectories)
        expected = self.config.expected_replicates
        payload = {"actual_replicates": actual, "expected_replicates": expected}
        if expected is not None and actual != expected:
            missing = max(expected - actual, 0)
            payload["missing_replicates"] = missing
            report["warnings"].append(
                f"Expected {expected} trajectories but found {actual} (missing={missing})."
            )
        report["checks"]["replicates"] = payload

    def _check_topology(self, report: Dict) -> None:
        topology = self.config.topology
        payload = {"topology": str(topology) if topology else None}
        if topology is None:
            report["warnings"].append("Topology file was not provided.")
            report["checks"]["topology"] = payload
            return

        if not topology.exists():
            report["errors"].append(f"Topology file not found: {topology}")
            report["checks"]["topology"] = payload
            return

        atom_count = self._infer_atom_count(topology)
        payload["atom_count"] = atom_count
        report["checks"]["topology"] = payload

    def _check_index_groups(self, report: Dict) -> None:
        index_file = self.config.index_file
        payload = {"index_file": str(index_file) if index_file else None}
        if index_file is None:
            report["warnings"].append("Index file was not provided; group sanity check skipped.")
            report["checks"]["index"] = payload
            return
        if not index_file.exists():
            report["errors"].append(f"Index file not found: {index_file}")
            report["checks"]["index"] = payload
            return

        groups = self._parse_index_groups(index_file)
        payload["groups"] = groups
        payload["group_count"] = len(groups)
        missing = [group for group in self.config.expected_groups if group not in groups]
        payload["missing_expected_groups"] = missing
        if missing:
            report["warnings"].append(f"Missing expected index groups: {', '.join(missing)}")
        report["checks"]["index"] = payload

    def _check_trajectories(self, report: Dict) -> None:
        entries = []
        for traj in self.config.trajectories:
            entry = {"path": str(traj)}
            if not traj.exists():
                report["errors"].append(f"Trajectory file not found: {traj}")
                entries.append(entry)
                continue

            entry.update(self._extract_trajectory_metadata(traj))
            self._check_sampling(entry, report)
            self._check_ligand_contact(entry, report)
            entries.append(entry)

        report["checks"]["trajectories"] = entries

    def _extract_trajectory_metadata(self, path: Path) -> Dict:
        metadata = {
            "frame_count": None,
            "timestep": None,
            "time_unit": None,
            "box": None,
            "pbc": None,
            "protein_selection_count": None,
            "ligand_selection_count": None,
            "min_ligand_protein_distance": None,
            "max_ligand_protein_min_distance": None,
            "uniform_sampling": None,
        }

        parsed_time = self._parse_text_times(path)
        if parsed_time is not None:
            times, inferred_unit = parsed_time
            metadata["frame_count"] = len(times)
            metadata["time_unit"] = inferred_unit
            if len(times) > 1:
                metadata["timestep"] = float(np.median(np.diff(times)))
                metadata.update(self._sampling_summary(times))
            return metadata

        # Optional MDAnalysis route for binary trajectories/topologies.
        if self.config.topology is None or not self.config.topology.exists():
            return metadata

        try:
            import MDAnalysis as mda  # noqa: WPS433
            from MDAnalysis.analysis import distances  # noqa: WPS433
        except Exception:
            return metadata

        try:
            universe = mda.Universe(str(self.config.topology), str(path))
            metadata["frame_count"] = len(universe.trajectory)
            if self.config.check_frame_times:
                times = np.array([ts.time for ts in universe.trajectory], dtype=float)
                if len(times) > 1:
                    metadata.update(self._sampling_summary(times))
            metadata["timestep"] = float(universe.trajectory.dt) if universe.trajectory.dt else None
            metadata["time_unit"] = "ps"
            dimensions = universe.trajectory.ts.dimensions
            if dimensions is not None and len(dimensions) >= 3:
                box = [float(v) for v in dimensions[:3]]
                metadata["box"] = box
                metadata["pbc"] = any(v > 0 for v in box)

            protein = universe.select_atoms(self.config.protein_selection)
            ligand = universe.select_atoms(self.config.ligand_selection)
            metadata["protein_selection_count"] = int(protein.n_atoms)
            metadata["ligand_selection_count"] = int(ligand.n_atoms)

            if protein.n_atoms > 0 and ligand.n_atoms > 0:
                per_frame = self._sample_min_distances(
                    universe,
                    protein,
                    ligand,
                    distances,
                )
                if per_frame:
                    metadata["min_ligand_protein_distance"] = float(min(per_frame))
                    metadata["max_ligand_protein_min_distance"] = float(max(per_frame))
                    cutoff = self.config.dissociation_cutoff
                    if cutoff is not None:
                        metadata["fraction_frames_ligand_detached"] = float(
                            np.mean(np.asarray(per_frame) > cutoff)
                        )
        except Exception as exc:
            logger.warning(f"Could not parse trajectory metadata for {path}: {exc}")

        return metadata

    def _sample_min_distances(self, universe, protein, ligand, distances_module) -> List[float]:
        """Per-frame minimum ligand-protein heavy-contact distance (Å) on sampled frames."""
        mins: List[float] = []
        sampled = 0
        for frame_idx, _ in enumerate(universe.trajectory):
            if frame_idx % self.config.sample_stride != 0:
                continue
            matrix = distances_module.distance_array(
                ligand.positions,
                protein.positions,
                box=universe.trajectory.ts.dimensions,
            )
            mins.append(float(np.min(matrix)))
            sampled += 1
            if self.config.max_frames and sampled >= self.config.max_frames:
                break
        return mins

    @staticmethod
    def _sampling_summary(times: np.ndarray) -> Dict:
        """Detect non-uniform frame spacing (gaps, overlaps from restarts, duplicates)."""
        diffs = np.diff(np.asarray(times, dtype=float))
        if len(diffs) == 0:
            return {}
        step = float(np.median(diffs))
        tol = max(abs(step) * 1e-3, 1e-6)
        irregular = np.abs(diffs - step) > tol
        return {
            "uniform_sampling": bool(not irregular.any()),
            "irregular_steps": int(irregular.sum()),
            "non_increasing_steps": int(np.sum(diffs <= 0)),
        }

    def _check_sampling(self, entry: Dict, report: Dict) -> None:
        if entry.get("uniform_sampling") is False:
            report["warnings"].append(
                f"{entry['path']}: {entry.get('irregular_steps')} irregular frame interval(s) "
                f"({entry.get('non_increasing_steps')} non-increasing). Check for restarts, "
                "duplicated frames or gaps (e.g. use gmx trjcat -settime) before analysis."
            )

    def _check_ligand_contact(self, entry: Dict, report: Dict) -> None:
        clash = self.config.min_distance_cutoff
        min_dist = entry.get("min_ligand_protein_distance")
        if clash is not None and min_dist is not None and min_dist < clash:
            report["warnings"].append(
                f"{entry['path']}: ligand-protein minimum distance {min_dist:.2f} Å is below "
                f"{clash} Å (possible steric clash or broken PBC)."
            )
        detach = self.config.dissociation_cutoff
        max_min = entry.get("max_ligand_protein_min_distance")
        if detach is not None and max_min is not None and max_min > detach:
            fraction = entry.get("fraction_frames_ligand_detached") or 0.0
            report["warnings"].append(
                f"{entry['path']}: ligand loses contact with the protein (> {detach} Å) in "
                f"{fraction:.0%} of sampled frames. Binding-site metrics and MM/PBSA averages "
                "over these frames are not meaningful, or the trajectory is not PBC-corrected."
            )

    def _infer_atom_count(self, topology: Path) -> Optional[int]:
        suffix = topology.suffix.lower()
        try:
            if suffix == ".gro":
                lines = topology.read_text().splitlines()
                if len(lines) >= 2:
                    return int(lines[1].strip())
            if suffix == ".pdb":
                count = 0
                for line in topology.read_text().splitlines():
                    if line.startswith(("ATOM", "HETATM")):
                        count += 1
                return count
            if suffix == ".psf":
                for line in topology.read_text().splitlines():
                    if "!NATOM" in line:
                        return int(line.split()[0])
        except Exception:
            return None
        return None

    def _parse_index_groups(self, index_file: Path) -> List[str]:
        groups: List[str] = []
        for line in index_file.read_text().splitlines():
            line = line.strip()
            if line.startswith("[") and line.endswith("]"):
                groups.append(line.strip("[] ").strip())
        return groups

    def _parse_text_times(self, path: Path) -> Optional[Tuple[np.ndarray, str]]:
        suffix = path.suffix.lower()
        if suffix not in {".xvg", ".dat", ".txt", ".csv"}:
            return None

        times: List[float] = []
        header_unit: Optional[str] = None
        with open(path, "r") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith("@") and "xaxis" in line and "label" in line:
                    label = line.lower()
                    if "(ps)" in label:
                        header_unit = "ps"
                    elif "(ns)" in label:
                        header_unit = "ns"
                    continue
                if line.startswith(("@", "#", ";")):
                    continue
                parts = line.replace(",", " ").split()
                if len(parts) < 2:
                    continue
                try:
                    times.append(float(parts[0]))
                except ValueError:
                    continue

        if not times:
            return None
        # GROMACS .xvg files declare the x-axis unit in the header; only fall back
        # to a magnitude heuristic (ambiguous for short runs) when it is absent.
        time_unit = header_unit
        if time_unit is None:
            time_unit = "ps" if max(times) > 1000 else "ns"
        return np.array(times, dtype=float), time_unit
