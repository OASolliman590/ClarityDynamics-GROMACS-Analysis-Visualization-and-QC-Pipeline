"""
Configuration classes for ligand-protein distance analysis.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Union, List, Optional

from ..md.config import PlotConfig


@dataclass
class DistanceConfig:
    """
    Configuration for ligand-protein distance analysis.

    Attributes:
        base_dir: Base directory for resolving relative paths
        output_dir: Output directory for plots/data
        protein_name: Protein name (for titles)
        ligand_name: Ligand name (for titles)
        topology: Topology file path
        trajectories: List of trajectory file paths (replicates)
        ligand_selection: MDAnalysis selection string for ligand
        residue_ids: List of residue numbers to analyze (optional)
        residue_selection: MDAnalysis selection string to derive residues (optional)
        auto_detect_residues: Auto-detect residues from trajectory if no residues provided
        auto_detect_cutoff: Cutoff (in output units) for auto-detected residues
        auto_detect_stride: Stride for auto-detection frame sampling
        auto_detect_max_frames: Maximum frames to sample for auto-detection
        auto_detect_selection: MDAnalysis selection string for candidate residues
        auto_detect_method: Distance method for auto-detection: 'com' or 'min'
        method: Distance method: 'com' (COM-to-COM) or 'min' (min atom distance)
        distance_scale: Multiply distances by this factor. MDAnalysis coordinates are
            always in Angstrom, so keep 1.0 for Å output and use 0.1 for nm output.
        time_unit: Unit of trajectory time ('ps', 'ns', 'frame', 'auto')
        time_step_ns: Optional timestep (ns) if time unit is 'frame'
        per_residue_plots: Whether to save per-residue plots
        combined_plot: Whether to save combined plot with all residues
        combined_show_std: Whether to show std bands in combined plot
        include_overall: Whether to include overall protein-ligand distance
        plot_config: PlotConfig for styling/output formats
    """
    base_dir: Path
    output_dir: Path
    protein_name: str
    ligand_name: str
    topology: Path
    trajectories: List[Path]
    ligand_selection: str = "resname UNK"
    residue_ids: Optional[List[Union[int, str]]] = None  # ints or chain-qualified "A:45"
    residue_selection: Optional[str] = None
    auto_detect_residues: bool = False
    auto_detect_cutoff: float = 5.0
    auto_detect_stride: int = 1
    auto_detect_max_frames: Optional[int] = None
    auto_detect_selection: str = "protein"
    auto_detect_method: str = "min"
    method: str = "com"
    distance_scale: float = 1.0
    time_unit: str = "auto"
    time_step_ns: Optional[float] = None
    per_residue_plots: bool = True
    combined_plot: bool = True
    combined_show_std: bool = False
    include_overall: bool = False
    plot_config: PlotConfig = field(default_factory=lambda: PlotConfig(
        template="ggplot2",
        style="simple",
        width=1200,
        height=600,
        scale=2,
        font_family="Times New Roman",
        font_size=24,
        save_formats=["html", "svg"],
    ))

    def __post_init__(self) -> None:
        self.base_dir = Path(self.base_dir)
        self.output_dir = Path(self.output_dir)
        self.topology = Path(self.topology)
        self.trajectories = [Path(p) for p in self.trajectories]

        if self.method not in ("com", "min"):
            raise ValueError("DistanceConfig.method must be 'com' or 'min'")

        if self.auto_detect_method not in ("com", "min"):
            raise ValueError("DistanceConfig.auto_detect_method must be 'com' or 'min'")

        if self.auto_detect_stride < 1:
            raise ValueError("DistanceConfig.auto_detect_stride must be >= 1")

        if self.auto_detect_residues and self.auto_detect_cutoff is None:
            raise ValueError("DistanceConfig.auto_detect_cutoff is required when auto_detect_residues is true")

        if self.residue_ids is None and not self.residue_selection and not self.auto_detect_residues:
            raise ValueError("Provide residue_ids/residue_selection or enable auto_detect_residues")


@dataclass
class DistanceCompareConfig:
    """
    Configuration for cross-ligand distance comparison.

    Attributes:
        output_dir: Output directory for comparison plots/data
        protein_name: Protein name (for titles)
        systems: List of DistanceConfig objects (one per ligand/system)
        compare_mode: 'intersection' or 'union' of residue sets
        per_residue_plots: Whether to generate per-residue comparison plots
        show_std: Whether to show std shading per ligand
        include_overall: Whether to include overall protein-ligand distance comparison
        heatmap_enabled: Whether to generate residue/time heatmaps per ligand
        heatmap_downsample: Downsample stride for heatmap time axis
        heatmap_max_frames: Optional cap on frames for heatmap
        delta_enabled: Whether to generate delta heatmaps between ligands
        delta_reference: Reference ligand name for delta maps (defaults to first)
        plot_config: PlotConfig for styling/output formats
    """
    output_dir: Path
    protein_name: str
    systems: List[DistanceConfig]
    compare_mode: str = "intersection"
    per_residue_plots: bool = True
    show_std: bool = False
    include_overall: bool = False
    heatmap_enabled: bool = False
    heatmap_downsample: int = 1
    heatmap_max_frames: Optional[int] = None
    delta_enabled: bool = False
    delta_reference: Optional[str] = None
    plot_config: PlotConfig = field(default_factory=lambda: PlotConfig(
        template="ggplot2",
        style="simple",
        width=1200,
        height=600,
        scale=2,
        font_family="Times New Roman",
        font_size=24,
        save_formats=["html", "svg"],
    ))

    def __post_init__(self) -> None:
        self.output_dir = Path(self.output_dir)

        if self.compare_mode not in ("intersection", "union"):
            raise ValueError("DistanceCompareConfig.compare_mode must be 'intersection' or 'union'")

        if not self.systems:
            raise ValueError("DistanceCompareConfig.systems must not be empty")

        if self.heatmap_downsample < 1:
            raise ValueError("DistanceCompareConfig.heatmap_downsample must be >= 1")
