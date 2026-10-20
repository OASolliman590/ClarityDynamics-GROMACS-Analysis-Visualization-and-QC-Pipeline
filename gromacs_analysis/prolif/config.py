"""Configuration classes for ProLIF interaction fingerprint analysis."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

from ..md.config import PlotConfig


@dataclass
class ProlifConfig:
    base_dir: Path
    output_dir: Path
    protein_name: str
    ligand_name: str
    topology: Path
    trajectories: List[Path]
    ligand_selection: str = "resname UNK"
    protein_selection: str = "protein"
    interaction_types: List[str] = field(default_factory=lambda: [
        "HBAcceptor",
        "HBDonor",
        "Hydrophobic",
        "PiStacking",
        "Cationic",
        "Anionic",
        "VdWContact",
    ])
    vicinity_cutoff: float = 5.0
    stride: int = 1
    n_jobs: Optional[int] = None
    save_pickle: bool = True
    outputs: List[str] = field(default_factory=lambda: [
        "barcode",
        "lignetwork",
        "occupancy",
        "tanimoto",
    ])
    lignetwork_display_all: bool = True
    lignetwork_count_aware: bool = True
    # "aggregate": interactions present in >= lignetwork_threshold of frames; "frame": first frame only.
    lignetwork_kind: str = "aggregate"
    lignetwork_threshold: float = 0.3
    occupancy_top_n: int = 25
    barcode_n_frame_ticks: int = 10
    barcode_residues_tick_location: str = "top"
    barcode_xlabel: str = "Frame"
    barcode_figsize: Tuple[float, float] = (12.0, 4.0)
    barcode_dpi: int = 150
    barcode_only_interacting_residues: bool = True
    barcode_min_interaction_frames: int = 1
    barcode_max_residues: Optional[int] = None
    plot_config: PlotConfig = field(default_factory=lambda: PlotConfig(
        template="ggplot2",
        style="simple",
        width=1200,
        height=600,
        scale=2,
        font_family="Times New Roman",
        font_size=20,
        save_formats=["html", "svg"],
    ))

    def __post_init__(self) -> None:
        self.base_dir = Path(self.base_dir)
        self.output_dir = Path(self.output_dir)
        self.topology = Path(self.topology)
        self.trajectories = [Path(p) for p in self.trajectories]
        if self.stride < 1:
            raise ValueError("ProlifConfig.stride must be >= 1")
        if self.barcode_n_frame_ticks < 1:
            raise ValueError("ProlifConfig.barcode_n_frame_ticks must be >= 1")
        if self.barcode_dpi < 1:
            raise ValueError("ProlifConfig.barcode_dpi must be >= 1")
        if self.barcode_min_interaction_frames < 1:
            raise ValueError("ProlifConfig.barcode_min_interaction_frames must be >= 1")
        if self.barcode_max_residues is not None and self.barcode_max_residues < 1:
            raise ValueError("ProlifConfig.barcode_max_residues must be >= 1 when provided")
        if self.barcode_residues_tick_location not in {"top", "bottom"}:
            raise ValueError("ProlifConfig.barcode_residues_tick_location must be 'top' or 'bottom'")
        if self.lignetwork_kind not in {"aggregate", "frame"}:
            raise ValueError("ProlifConfig.lignetwork_kind must be 'aggregate' or 'frame'")
        if not 0.0 <= float(self.lignetwork_threshold) <= 1.0:
            raise ValueError("ProlifConfig.lignetwork_threshold must be between 0 and 1")
        if len(self.barcode_figsize) != 2:
            raise ValueError("ProlifConfig.barcode_figsize must contain exactly 2 values")


@dataclass
class ProlifCompareConfig:
    output_dir: Path
    protein_name: str
    systems: List[ProlifConfig]
    compare_mode: str = "intersection"
    plot_config: PlotConfig = field(default_factory=lambda: PlotConfig(
        template="ggplot2",
        style="simple",
        width=1200,
        height=600,
        scale=2,
        font_family="Times New Roman",
        font_size=20,
        save_formats=["html", "svg"],
    ))

    def __post_init__(self) -> None:
        self.output_dir = Path(self.output_dir)
        if self.compare_mode not in ("intersection", "union"):
            raise ValueError("ProlifCompareConfig.compare_mode must be 'intersection' or 'union'")
        if not self.systems:
            raise ValueError("ProlifCompareConfig.systems must not be empty")
