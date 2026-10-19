"""
Distance data processor using MDAnalysis.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from .config import DistanceConfig

logger = logging.getLogger(__name__)


class DistanceProcessor:
    """
    Compute ligand-protein residue distances across trajectories.
    """

    def __init__(self, config: DistanceConfig):
        self.config = config

    def load_data(self) -> Dict:
        """
        Load trajectories and compute distances.

        Returns:
            dict with keys:
              - time: np.ndarray
              - residues: List[dict] with resid, resname, label, mean, std
              - method: str
              - overall: Optional[dict] with mean, std
        """
        try:
            import MDAnalysis as mda  # noqa: WPS433 (runtime import)
            from MDAnalysis.lib.distances import distance_array  # noqa: WPS433
        except ImportError as exc:
            raise ImportError("MDAnalysis is required for distance analysis") from exc

        topology = self._resolve_path(self.config.topology)
        trajectories = [self._resolve_path(p) for p in self.config.trajectories]

        if not topology.exists():
            raise FileNotFoundError(f"Topology not found: {topology}")
        for traj in trajectories:
            if not traj.exists():
                raise FileNotFoundError(f"Trajectory not found: {traj}")

        # Use first universe to determine residue labels
        universe = mda.Universe(str(topology), str(trajectories[0]))
        residue_ids, residue_labels = self._resolve_residues(universe)

        if not residue_ids:
            raise ValueError("No residues resolved for distance analysis")

        all_distances: Dict[int, List[np.ndarray]] = {rid: [] for rid in residue_ids}
        all_times: List[np.ndarray] = []
        overall_series: List[np.ndarray] = []

        for traj_path in trajectories:
            u = mda.Universe(str(topology), str(traj_path))
            ligand = u.select_atoms(self.config.ligand_selection)
            if ligand.n_atoms == 0:
                raise ValueError(f"Ligand selection returned 0 atoms: {self.config.ligand_selection}")

            residue_map = {res.resid: res for res in u.residues}
            residues = [residue_map[rid] for rid in residue_ids if rid in residue_map]
            if not residues:
                raise ValueError("Residue list not found in topology for selected IDs")

            time_ns = []
            distances_this_rep: Dict[int, List[float]] = {rid: [] for rid in residue_ids}
            overall_this_rep: List[float] = []

            for ts in u.trajectory:
                time_ns.append(self._convert_time(ts))
                box = self._box(ts)

                if self.config.method == "com":
                    ligand_com = ligand.center_of_mass()
                    for res in residues:
                        res_com = res.atoms.center_of_mass()
                        dist = self._com_distance(res_com, ligand_com, box) * self.config.distance_scale
                        distances_this_rep[res.resid].append(dist)
                else:
                    lig_positions = ligand.positions
                    for res in residues:
                        res_positions = res.atoms.positions
                        dist = np.min(distance_array(res_positions, lig_positions, box=box)) * self.config.distance_scale
                        distances_this_rep[res.resid].append(dist)

                if self.config.include_overall:
                    protein = u.select_atoms("protein")
                    if protein.n_atoms > 0:
                        overall_dist = self._com_distance(
                            protein.center_of_mass(), ligand.center_of_mass(), box
                        ) * self.config.distance_scale
                        overall_this_rep.append(overall_dist)

            time_array = np.array(time_ns, dtype=float)
            all_times.append(time_array)

            for rid in residue_ids:
                all_distances[rid].append(np.array(distances_this_rep[rid], dtype=float))

            if self.config.include_overall:
                overall_series.append(np.array(overall_this_rep, dtype=float))

        # Align lengths
        min_len = min(len(t) for t in all_times)
        time_ref = all_times[0][:min_len]

        results = []
        for rid in residue_ids:
            reps = [arr[:min_len] for arr in all_distances[rid]]
            rep_array = np.vstack(reps)
            mean = rep_array.mean(axis=0)
            std = rep_array.std(axis=0)
            results.append({
                "resid": rid,
                "resname": residue_labels.get(rid, "UNK"),
                "label": f"{residue_labels.get(rid, 'UNK')}{rid}",
                "mean": mean,
                "std": std,
            })

        overall = None
        if self.config.include_overall and overall_series:
            overall_array = np.vstack([arr[:min_len] for arr in overall_series])
            overall = {
                "mean": overall_array.mean(axis=0),
                "std": overall_array.std(axis=0),
            }

        return {
            "time": time_ref,
            "residues": results,
            "method": self.config.method,
            "overall": overall,
        }

    def _resolve_residues(self, universe) -> Tuple[List[int], Dict[int, str]]:
        residue_ids: List[int] = []
        residue_labels: Dict[int, str] = {}

        if self.config.auto_detect_residues:
            return self._auto_detect_residues(universe)

        if self.config.residue_selection:
            selection = universe.select_atoms(self.config.residue_selection)
            residues = selection.residues
            residue_ids = [int(res.resid) for res in residues]
            residue_labels = {int(res.resid): res.resname for res in residues}
            return residue_ids, residue_labels

        if self.config.residue_ids:
            residue_ids = [int(r) for r in self.config.residue_ids]
            # Map residue names from universe
            for res in universe.residues:
                if res.resid in residue_ids and res.resid not in residue_labels:
                    residue_labels[int(res.resid)] = res.resname
            missing = [rid for rid in residue_ids if rid not in residue_labels]
            if missing:
                logger.warning(f"Residues not found in topology (using UNK): {missing}")
            return residue_ids, residue_labels

        return residue_ids, residue_labels

    def _auto_detect_residues(self, universe) -> Tuple[List[int], Dict[int, str]]:
        try:
            from MDAnalysis.lib.distances import distance_array  # noqa: WPS433
        except ImportError as exc:
            raise ImportError("MDAnalysis is required for distance analysis") from exc

        ligand = universe.select_atoms(self.config.ligand_selection)
        if ligand.n_atoms == 0:
            raise ValueError(f"Ligand selection returned 0 atoms: {self.config.ligand_selection}")

        selection = self.config.auto_detect_selection or "protein"
        candidates = universe.select_atoms(selection)
        residues = candidates.residues
        if not residues:
            raise ValueError(f"Auto-detect selection returned 0 residues: {selection}")

        cutoff = float(self.config.auto_detect_cutoff)
        stride = max(1, int(self.config.auto_detect_stride))
        max_frames = self.config.auto_detect_max_frames
        method = self.config.auto_detect_method

        min_dist: Dict[int, float] = {int(res.resid): np.inf for res in residues}
        frames_used = 0

        for ts in universe.trajectory[::stride]:
            if max_frames is not None and frames_used >= int(max_frames):
                break
            frames_used += 1
            box = self._box(ts)

            if method == "com":
                lig_com = ligand.center_of_mass()
                for res in residues:
                    dist = self._com_distance(res.atoms.center_of_mass(), lig_com, box) * self.config.distance_scale
                    resid = int(res.resid)
                    if dist < min_dist[resid]:
                        min_dist[resid] = dist
            else:
                lig_positions = ligand.positions
                for res in residues:
                    dist = np.min(distance_array(res.atoms.positions, lig_positions, box=box)) * self.config.distance_scale
                    resid = int(res.resid)
                    if dist < min_dist[resid]:
                        min_dist[resid] = dist

        selected = [rid for rid, dist in min_dist.items() if dist <= cutoff]
        selected.sort()

        if not selected:
            raise ValueError(
                "Auto-detect found no residues within cutoff. "
                "Increase auto_detect_cutoff or check ligand selection."
            )

        residue_labels: Dict[int, str] = {}
        for res in residues:
            resid = int(res.resid)
            if resid in selected and resid not in residue_labels:
                residue_labels[resid] = res.resname

        logger.info(
            "Auto-detected %d residues within %.2f (selection=%s, method=%s, frames=%s, stride=%d)",
            len(selected),
            cutoff,
            selection,
            method,
            "all" if max_frames is None else str(max_frames),
            stride,
        )

        return selected, residue_labels

    @staticmethod
    def _box(ts) -> Optional[np.ndarray]:
        """Unit-cell dimensions for minimum-image distances, or None without a valid box."""
        dims = getattr(ts, "dimensions", None)
        if dims is None:
            return None
        dims = np.asarray(dims, dtype=np.float32)
        if len(dims) < 6 or not np.all(np.isfinite(dims[:3])) or np.any(dims[:3] <= 0):
            return None
        return dims

    @staticmethod
    def _com_distance(a: np.ndarray, b: np.ndarray, box: Optional[np.ndarray]) -> float:
        """Distance between two points under the minimum-image convention.

        Centres of mass are only meaningful if each group is whole; trajectories
        should be pre-processed (e.g. ``gmx trjconv -pbc mol -center``).
        """
        delta = np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64)
        if box is not None:
            from MDAnalysis.lib.distances import minimize_vectors  # noqa: WPS433

            delta = minimize_vectors(delta.reshape(1, 3).astype(np.float32), box)[0]
        return float(np.linalg.norm(delta))

    def _convert_time(self, ts) -> float:
        if self.config.time_unit == "frame":
            if self.config.time_step_ns is not None:
                return float(ts.frame) * self.config.time_step_ns
            return float(ts.frame)

        if self.config.time_unit in ("ps", "ns"):
            time_val = float(ts.time) if ts.time is not None else float(ts.frame)
            if self.config.time_unit == "ps":
                return time_val / 1000.0
            return time_val

        # auto
        if ts.time is not None:
            return float(ts.time) / 1000.0
        if self.config.time_step_ns is not None:
            return float(ts.frame) * self.config.time_step_ns
        return float(ts.frame)

    def _resolve_path(self, value: Path) -> Path:
        if value.is_absolute():
            return value
        return (self.config.base_dir / value).resolve()
