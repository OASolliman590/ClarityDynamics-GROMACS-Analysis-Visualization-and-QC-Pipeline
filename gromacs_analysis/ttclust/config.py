"""Configuration classes for TTClust trajectory clustering."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from ..md.config import PlotConfig


@dataclass
class TTClustConfig:
    """Configuration for one TTClust run."""

    base_dir: Path
    output_dir: Path
    protein_name: str
    system_name: str
    trajectories: List[Path]
    topology: Optional[Path] = None
    executable: Optional[str] = None
    stride: int = 1
    logfile: str = "clustering.log"
    select_traj: str = "all"
    select_alignment: str = "backbone"
    select_rmsd: str = "backbone"
    method: str = "ward"
    cutoff: Optional[float] = None
    n_groups: Optional[str] = None
    autoclust: bool = True
    interactive_matrix: bool = False
    axis: Optional[str] = None
    limit_matrix: Optional[int] = None
    extra_args: List[str] = field(default_factory=list)
    environment: Dict[str, str] = field(default_factory=dict)
    # Styling/size for the cluster summary figures drawn from clustering.log.
    plot_config: PlotConfig = field(default_factory=lambda: PlotConfig(
        template="plotly_white",
        font_family="Arial",
        font_size=16,
        save_formats=["html", "svg"],
    ))

    def __post_init__(self) -> None:
        self.base_dir = Path(self.base_dir)
        self.output_dir = Path(self.output_dir)
        self.trajectories = [Path(p) for p in self.trajectories]
        self.topology = Path(self.topology) if self.topology else None

        if not self.trajectories:
            raise ValueError("TTClustConfig.trajectories must not be empty")
        if self.stride < 1:
            raise ValueError("TTClustConfig.stride must be >= 1")
        if self.method not in ("single", "complete", "average", "weighted", "centroid", "median", "ward"):
            raise ValueError("TTClustConfig.method must be one of scipy linkage methods supported by TTClust")
        if self.limit_matrix is not None and self.limit_matrix < 1:
            raise ValueError("TTClustConfig.limit_matrix must be >= 1")
        if not self.logfile:
            raise ValueError("TTClustConfig.logfile must not be empty")
