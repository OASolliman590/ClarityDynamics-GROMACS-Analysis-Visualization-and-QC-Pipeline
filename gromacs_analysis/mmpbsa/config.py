"""
Configuration classes for MMPBSA analysis

Based on MD analysis config pattern, adapted for MMPBSA needs.
"""

from dataclasses import dataclass, field
from typing import List, Optional
from pathlib import Path

from ..md.config import PlotConfig
from ..utils.publication import PublicationConfig


@dataclass
class MMPBSASystemConfig:
    """
    Configuration for an MMPBSA system (typically HOLO).
    
    Attributes:
        name: System name (e.g., "LigandA")
        dir_pattern: Directory naming pattern (e.g., "3-Results_Holo/r/HOLO_{}")
        replicates: Number of replicates (1 for single, 3+ for replicates)
        results_file_pattern: Pattern for results file (e.g., "FINAL_RESULTS_MMPBSA*.dat")
        decomp_file_pattern: Pattern for decomposition file (e.g., "FINAL_DECOMP_MMPBSA*.dat")
    """
    name: str
    dir_pattern: str
    replicates: int = 1
    results_file_pattern: str = "FINAL_RESULTS_MMPBSA*.dat"
    decomp_file_pattern: str = "FINAL_DECOMP_MMPBSA*.dat"
    
    def __post_init__(self):
        """Validate configuration"""
        if self.replicates < 1:
            raise ValueError(f"Replicates must be >= 1, got {self.replicates}")


@dataclass
class MMPBSAConfig:
    """
    Complete configuration for MMPBSA analysis.
    
    Attributes:
        base_dir: Base directory containing data files
        output_dir: Output directory for results
        protein_name: Name of the protein
        ligand_name: Name of the ligand (for labeling)
        systems: List of MMPBSASystemConfig objects
        file_format: File format ('csv', 'dat', or 'auto' for auto-detection)
        plot_config: PlotConfig object for visualization settings
        amino_acids: Optional list of amino acid numbers for decomposition x-axis
        compare_systems: Whether to generate cross-ligand comparison plots
        compare_binding_key: Binding energy key to compare (TOTAL, GGAS, GSOLV)
        compare_components: Whether to generate component-by-component comparisons
        compare_component_order: Optional list of component names to compare
    """
    base_dir: Path
    output_dir: Path
    protein_name: str
    ligand_name: str
    systems: List[MMPBSASystemConfig]
    file_format: str = 'auto'  # 'csv', 'dat', or 'auto'
    plot_config: PlotConfig = field(default_factory=PlotConfig)
    amino_acids: Optional[List[int]] = None
    compare_systems: bool = False
    compare_binding_key: str = "TOTAL"
    compare_components: bool = False
    compare_component_order: Optional[List[str]] = None
    incomplete_replicates_policy: str = "warn"  # ignore | warn | error
    require_decomp_replicates: bool = False
    # Exclude frames before this "Frame #" (e.g. the equilibration start from md/data/equilibration.json).
    # Applies to frame-level (CSV) results; DAT summaries cannot be re-averaged and are left unchanged.
    start_frame: Optional[int] = None
    publication: PublicationConfig = field(default_factory=PublicationConfig)
    
    def __post_init__(self):
        """Validate configuration"""
        valid_formats = ['csv', 'dat', 'auto']
        if self.file_format not in valid_formats:
            raise ValueError(f"Invalid file_format: {self.file_format}. Choose from {valid_formats}")

        valid_binding_keys = {"TOTAL", "GGAS", "GSOLV"}
        if self.compare_binding_key not in valid_binding_keys:
            raise ValueError(
                f"Invalid compare_binding_key: {self.compare_binding_key}. Choose from {sorted(valid_binding_keys)}"
            )

        if self.compare_component_order is not None and not self.compare_component_order:
            raise ValueError("compare_component_order must not be empty when provided")

        valid_policies = {"ignore", "warn", "error"}
        if self.incomplete_replicates_policy not in valid_policies:
            raise ValueError(
                "Invalid incomplete_replicates_policy: "
                f"{self.incomplete_replicates_policy}. Choose from {sorted(valid_policies)}"
            )
        
        if not self.systems:
            raise ValueError("At least one system must be specified")
        
        # Set default plot config if not customized
        if self.plot_config.template == 'plotly_white' and self.plot_config.font_family == 'Arial':
            # Apply MMPBSA defaults (matching MD analysis)
            self.plot_config.template = 'ggplot2'
            self.plot_config.font_family = 'Times New Roman'
            self.plot_config.font_size = 24
