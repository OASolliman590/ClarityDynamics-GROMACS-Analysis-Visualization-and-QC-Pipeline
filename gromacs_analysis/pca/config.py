"""
Configuration classes for PCA analysis
"""

from dataclasses import dataclass, field
from typing import Dict, Optional, List, Tuple
from pathlib import Path

from ..md.config import PlotConfig


@dataclass
class PCASystemConfig:
    """
    Configuration for a PCA system (APO or HOLO).
    
    Attributes:
        name: System name (e.g., "APO", "LigandA")
        is_apo: Whether this is an APO system
        ligand_name: Ligand name (for HOLO systems, extracted from files if not provided)
    """
    name: str
    is_apo: bool = False
    ligand_name: Optional[str] = None


@dataclass
class FELConfig:
    """Configuration for Free Energy Landscape (FEL) plots."""
    enabled: bool = False
    pairs: List[str] = field(default_factory=lambda: ["12"])
    contours: int = 7
    step: int = 2
    colorscale: str = "Viridis"  # perceptually uniform; avoid rainbow scales such as jet
    overlay_projection: bool = True
    include_probability: bool = False
    include_entropy: bool = False
    include_enthalpy: bool = False


@dataclass
class ClusteringConfig:
    """Configuration for PCA clustering."""
    enabled: bool = False
    method: str = "kmeans"  # kmeans | dbscan
    pair: str = "12"
    n_clusters: int = 3
    eps: float = 1.5
    min_samples: int = 5
    downsample: int = 1
    max_points: Optional[int] = None


@dataclass
class PCATerminalDetectionConfig:
    """Configuration for RMSF-guided terminal trimming before PCA indexing."""
    enabled: bool = True
    rmsf_file: str = "RMSF_.dat"
    pdb_file: str = "step3_input.pdb"
    gro_file: str = "step3_input.gro"
    smoothing_window: int = 3
    mad_multiplier: float = 3.0
    stable_run: int = 3
    max_terminal_fraction: float = 0.25
    core_residue_range: Optional[Tuple[int, int]] = None

    def __post_init__(self):
        self.smoothing_window = max(1, int(self.smoothing_window))
        self.stable_run = max(1, int(self.stable_run))
        self.mad_multiplier = float(self.mad_multiplier)
        self.max_terminal_fraction = max(0.0, min(0.5, float(self.max_terminal_fraction)))


@dataclass
class PCAGromacsRunConfig:
    """Configuration for optional automatic GROMACS PCA output generation."""
    enabled: bool = False
    input_dir: Optional[Path] = None
    gmx: str = "gmx"
    tpr: str = "step5_production.tpr"
    trajectory: str = "step5_production.xtc"
    no_pbc_trajectory: str = "step5_production_noPBC.xtc"
    index: str = "new_index.ndx"
    protein_group: str = "Protein"
    output_group: str = "System"
    core_group_name: str = "core_no_terminals"
    calpha_group_name: str = "core_no_terminals_Calpha"
    n_pcs: int = 10
    window_cosine_enabled: bool = True
    window_max_ns: int = 100
    window_step_ns: int = 10
    extreme_nframes: int = 30
    overwrite: bool = True
    terminal_detection: PCATerminalDetectionConfig = field(default_factory=PCATerminalDetectionConfig)

    def __post_init__(self):
        if self.input_dir is not None:
            self.input_dir = Path(self.input_dir)
        self.n_pcs = max(1, int(self.n_pcs))
        self.window_max_ns = max(1, int(self.window_max_ns))
        self.window_step_ns = max(1, int(self.window_step_ns))
        self.extreme_nframes = max(1, int(self.extreme_nframes))


@dataclass
class PCAConfig:
    """
    Complete configuration for PCA analysis.
    
    Attributes:
        base_dir: Base directory containing PCA files (e.g., '/path/to/2-PCA')
        output_dir: Output directory for plots (default: base_dir/PCA_Vis)
        protein_name: Name of the protein
        ligand_name: Explicit ligand name (overrides auto-detection)
        apo_base_dir: Optional separate directory for APO files (if different from base_dir)
        systems: List of PCASystemConfig objects
        n_pcs: Number of principal components to analyze (default: 10)
        plot_config: PlotConfig object for visualization settings
    """
    base_dir: Path
    output_dir: Optional[Path] = None
    protein_name: str = "Protein"
    ligand_name: Optional[str] = None  # Explicit ligand name
    apo_base_dir: Optional[Path] = None  # Optional separate APO directory
    systems: list = field(default_factory=list)
    n_pcs: int = 10
    # Ligand-bound system colour (hex). None -> first categorical slot; batches assign slots in order.
    color: Optional[str] = None
    plot_config: PlotConfig = field(default_factory=PlotConfig)
    fel: FELConfig = field(default_factory=FELConfig)
    clustering: ClusteringConfig = field(default_factory=ClusteringConfig)
    run_gromacs: PCAGromacsRunConfig = field(default_factory=PCAGromacsRunConfig)
    
    def __post_init__(self):
        """Set default output directory if not provided."""
        self.base_dir = Path(self.base_dir)
        if self.output_dir is None:
            self.output_dir = self.base_dir / 'PCA_Vis'
        else:
            self.output_dir = Path(self.output_dir)
        if self.apo_base_dir is not None:
            self.apo_base_dir = Path(self.apo_base_dir)
        
        # Set default plot config if not customized
        if self.plot_config.template == 'plotly_white' and self.plot_config.font_family == 'Arial':
            # Apply PCA defaults (matching existing code)
            self.plot_config.template = 'plotly_white'
            self.plot_config.font_family = 'Times New Roman'
            self.plot_config.font_size = 24
            self.plot_config.width = 1200
            self.plot_config.height = 600
            self.plot_config.scale = 2


@dataclass
class PCACompareConfig:
    """
    Configuration for cross-ligand PCA comparison.

    Attributes:
        output_dir: Output directory for comparison plots/data
        protein_name: Protein name (for titles)
        systems: List of PCAConfig objects (one per ligand/system)
        compare_pairs: List of 2D projection pairs (e.g., ["12", "13"])
        compare_scree: Whether to compare eigenvalue scree/cumulative variance
        compare_cosine: Whether to compare cosine content across ligands
        downsample: Optional stride for 2D projection points
        max_points: Optional max points per ligand projection
        plot_config: PlotConfig for styling/output formats
    """
    output_dir: Path
    protein_name: str
    systems: List[PCAConfig]
    compare_pairs: List[str] = field(default_factory=lambda: ["12"])
    compare_scree: bool = True
    compare_cosine: bool = False
    downsample: Optional[int] = None
    max_points: Optional[int] = None
    plot_config: PlotConfig = field(default_factory=PlotConfig)

    def __post_init__(self):
        self.output_dir = Path(self.output_dir)
        if not self.systems:
            raise ValueError("PCACompareConfig.systems must not be empty")


@dataclass
class PCAJointSystem:
    """One system of a joint PCA: topology plus one or more replicate trajectories."""
    name: str
    topology: Path
    trajectories: List[Path]
    is_apo: bool = False
    color: Optional[str] = None

    def __post_init__(self):
        self.topology = Path(self.topology)
        self.trajectories = [Path(p) for p in self.trajectories]
        if not self.trajectories:
            raise ValueError(f"PCAJointSystem {self.name}: trajectories must not be empty")


@dataclass
class PCAJointConfig:
    """Joint PCA of several systems on one shared set of eigenvectors."""
    output_dir: Path
    protein_name: str
    systems: List[PCAJointSystem]
    base_dir: Optional[Path] = None
    selection: str = "protein and name CA"
    n_components: int = 10
    n_overlap: int = 10
    stride: int = 1
    max_frames_per_system: Optional[int] = None
    fit_iterations: int = 2
    # Sampling checks: maximum cosine content of PC1 per replicate and minimum RMSIP
    # between the essential subspaces of each system's first and second half of frames.
    convergence: Dict[str, float] = field(default_factory=dict)
    plot_config: PlotConfig = field(default_factory=lambda: PlotConfig(
        template="plotly_white",
        font_family="Arial",
        font_size=16,
        save_formats=["html", "svg"],
    ))

    def __post_init__(self):
        self.output_dir = Path(self.output_dir)
        if self.base_dir is not None:
            self.base_dir = Path(self.base_dir)
        if not self.systems:
            raise ValueError("PCAJointConfig.systems must not be empty")
        names = [s.name for s in self.systems]
        if len(set(names)) != len(names):
            raise ValueError(f"PCAJointConfig system names must be unique: {names}")
        if self.n_components < 1 or self.n_overlap < 1:
            raise ValueError("PCAJointConfig.n_components and n_overlap must be >= 1")
        if self.stride < 1 or self.fit_iterations < 1:
            raise ValueError("PCAJointConfig.stride and fit_iterations must be >= 1")
