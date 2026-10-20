"""ProLIF interaction fingerprint processor."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

from .config import ProlifConfig

logger = logging.getLogger(__name__)


class ProlifProcessor:
    def __init__(self, config: ProlifConfig):
        self.config = config

    def run(self) -> Dict:
        try:
            import MDAnalysis as mda  # noqa: WPS433
            import prolif as plf  # noqa: WPS433
        except ImportError as exc:
            raise ImportError("ProLIF + MDAnalysis are required for ProLIF analysis") from exc

        topology = self._resolve_path(self.config.topology)
        trajectories = [self._resolve_path(p) for p in self.config.trajectories]

        if not topology.exists():
            raise FileNotFoundError(f"Topology not found: {topology}")
        for traj in trajectories:
            if not traj.exists():
                raise FileNotFoundError(f"Trajectory not found: {traj}")

        replicate_outputs: List[Dict] = []
        chemistry_notes: List[Dict] = []
        build_count_fingerprint = "lignetwork" in self.config.outputs and self.config.lignetwork_count_aware

        for rep_idx, trajectory_path in enumerate(trajectories, start=1):
            universe = mda.Universe(str(topology), str(trajectory_path))

            # Ensure elements are present for ProLIF.
            if not hasattr(universe.atoms, "elements"):
                self._guess_elements(universe)

            ligand_selection = universe.select_atoms(self.config.ligand_selection)
            protein_selection = universe.select_atoms(self.config.protein_selection)
            if ligand_selection.n_atoms == 0:
                raise ValueError(
                    f"Ligand selection returned 0 atoms for replicate {rep_idx}: {self.config.ligand_selection}"
                )
            if protein_selection.n_atoms == 0:
                raise ValueError(
                    f"Protein selection returned 0 atoms for replicate {rep_idx}: {self.config.protein_selection}"
                )

            chemistry = self._assess_ligand_chemistry(ligand_selection, rep_idx)
            chemistry_notes.append(chemistry)

            if self.config.interaction_types:
                fp = plf.Fingerprint(self.config.interaction_types, vicinity_cutoff=self.config.vicinity_cutoff)
            else:
                fp = plf.Fingerprint(vicinity_cutoff=self.config.vicinity_cutoff)

            trajectory = self._iter_trajectory(universe)
            fp.run(
                trajectory,
                ligand_selection,
                protein_selection,
                residues=None,
                progress=True,
                n_jobs=self.config.n_jobs,
            )

            fp_count = None
            if build_count_fingerprint:
                if self.config.interaction_types:
                    fp_count = plf.Fingerprint(
                        self.config.interaction_types,
                        vicinity_cutoff=self.config.vicinity_cutoff,
                        count=True,
                    )
                else:
                    fp_count = plf.Fingerprint(vicinity_cutoff=self.config.vicinity_cutoff, count=True)
                fp_count.run(
                    self._iter_trajectory(universe),
                    ligand_selection,
                    protein_selection,
                    residues=None,
                    progress=True,
                    n_jobs=self.config.n_jobs,
                )

            replicate_outputs.append(
                {
                    "replicate": rep_idx,
                    "trajectory": str(trajectory_path),
                    "fingerprint": fp,
                    "fingerprint_count": fp_count,
                    "ligand_mol": plf.Molecule.from_mda(ligand_selection),
                    "protein_mol": plf.Molecule.from_mda(protein_selection),
                    "chemistry": chemistry,
                }
            )

        return {
            "replicates": replicate_outputs,
            "ligand_mol": replicate_outputs[0]["ligand_mol"] if replicate_outputs else None,
            "protein_mol": replicate_outputs[0]["protein_mol"] if replicate_outputs else None,
            "chemistry_requirements": self._chemistry_requirements(),
            "chemistry_summary": self._summarize_chemistry(chemistry_notes),
        }

    @staticmethod
    def _guess_elements(universe) -> None:
        """Add element types (required by ProLIF) across MDAnalysis versions."""
        if hasattr(universe, "guess_TopologyAttrs"):  # MDAnalysis >= 2.8
            universe.guess_TopologyAttrs(to_guess=["elements"])
            return
        from MDAnalysis.topology.guessers import guess_types  # noqa: WPS433

        universe.add_TopologyAttr("elements", guess_types(universe.atoms.names))

    def _iter_trajectory(self, universe):
        if self.config.stride > 1:
            return universe.trajectory[:: self.config.stride]
        return universe.trajectory

    def _resolve_path(self, path: Path) -> Path:
        if path.is_absolute():
            return path
        return (self.config.base_dir / path).resolve()

    def _assess_ligand_chemistry(self, ligand_selection, rep_idx: int) -> Dict:
        elements = list(getattr(ligand_selection, "elements", []))
        if not elements and hasattr(ligand_selection, "names"):
            elements = [str(name)[0] for name in ligand_selection.names]

        heavy_atoms = sum(1 for element in elements if str(element).upper() != "H")
        unknown_elements = sorted(
            {
                str(element)
                for element in elements
                if str(element).upper() in {"", "X", "DU"}
            }
        )
        warnings: List[str] = []
        if heavy_atoms < 3:
            warnings.append("Ligand has very low heavy-atom count; interaction profiling may be unstable.")
        if unknown_elements:
            warnings.append(f"Ligand includes unknown element tags: {', '.join(unknown_elements)}")

        return {
            "replicate": rep_idx,
            "atom_count": int(ligand_selection.n_atoms),
            "heavy_atom_count": int(heavy_atoms),
            "unique_elements": sorted({str(element) for element in elements}),
            "unknown_elements": unknown_elements,
            "warnings": warnings,
        }

    def _summarize_chemistry(self, chemistry_notes: List[Dict]) -> Dict:
        all_warnings: List[str] = []
        min_atoms: Optional[int] = None
        min_heavy: Optional[int] = None
        for note in chemistry_notes:
            atom_count = int(note.get("atom_count", 0))
            heavy = int(note.get("heavy_atom_count", 0))
            min_atoms = atom_count if min_atoms is None else min(min_atoms, atom_count)
            min_heavy = heavy if min_heavy is None else min(min_heavy, heavy)
            all_warnings.extend(note.get("warnings", []))

        return {
            "replicate_count": len(chemistry_notes),
            "min_atom_count": min_atoms,
            "min_heavy_atom_count": min_heavy,
            "warnings": sorted(set(all_warnings)),
        }

    @staticmethod
    def _chemistry_requirements() -> Dict:
        return {
            "requires_ligand_selection_atoms": True,
            "requires_protein_selection_atoms": True,
            "requires_reasonable_heavy_atom_count": ">= 3 recommended",
            "requires_resolved_elements": True,
        }
