"""Configuration for ligand ranking."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional


@dataclass
class RankingConfig:
    output_dir: Path
    protein_name: str
    md_stats_file: Optional[Path] = None
    mmpbsa_data_dir: Optional[Path] = None
    distance_summary_file: Optional[Path] = None
    # Only metrics listed here contribute to rank_score. Defaults are limited to
    # observables with a defensible direction for ligand binding: ligand pose
    # stability, ligand-pocket proximity and MM/PBSA ΔTOTAL. Protein-wide
    # descriptors (rmsd_prot, rmsf, rog, sasa, intra-protein hbonds) have no
    # intrinsic "better" direction; add them explicitly if your hypothesis needs them.
    metric_directions: Dict[str, str] = field(default_factory=lambda: {
        "rmsd_lig": "lower",
        "rmsd_complex": "lower",
        "comcom": "lower",
        "distance_mean": "lower",
        "mmpbsa_total": "lower",
    })

    def __post_init__(self) -> None:
        self.output_dir = Path(self.output_dir)
        if self.md_stats_file:
            self.md_stats_file = Path(self.md_stats_file)
        if self.mmpbsa_data_dir:
            self.mmpbsa_data_dir = Path(self.mmpbsa_data_dir)
        if self.distance_summary_file:
            self.distance_summary_file = Path(self.distance_summary_file)
