"""
Configuration classes for MD analysis

Combines best practices from:
- JSON config (holo_md_analysis.py)
- Dataclass config (viz_package/holo_apo_vis.py)
- Dictionary config (mds_parser_replicates.py)
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Union
from pathlib import Path
import json


@dataclass
class AnalysisMetric:
    """
    Definition of an analysis metric.
    
    Attributes:
        name: Internal name (e.g., 'rmsd_prot')
        file_pattern: Pattern to match data files (e.g., 'RMSD_*.dat')
        title: Display title for plots
        ylabel: Y-axis label with units
        is_holo_only: Whether this metric only applies to HOLO systems
        data_format: Format of data ('single_column' or 'two_column')
    """
    name: str
    file_pattern: str
    title: str
    ylabel: str
    is_holo_only: bool = False
    data_format: str = 'two_column'


@dataclass
class SystemConfig:
    """
    Configuration for a molecular system (APO or HOLO).
    
    Attributes:
        name: System name (e.g., 'LigandA', 'APO')
        dir_pattern: Directory naming pattern (e.g., 'ProteinX_HOLO_{}')
        dir_names: Optional explicit list of replicate directories (overrides dir_pattern)
        is_apo: Whether this is an APO system
        replicates: Number of replicates (1 for single trajectory)
        color: RGB color tuple for plots
    """
    name: str
    dir_pattern: str = ""
    is_apo: bool = False
    replicates: int = 1
    color: Optional[Tuple[int, int, int]] = None
    dir_names: Optional[List[str]] = None
    
    def __post_init__(self):
        """Set default colors based on system type."""
        if self.dir_names:
            self.replicates = len(self.dir_names)
        if self.color is None:
            if self.is_apo:
                self.color = (0, 125, 0)  # Green for APO (legacy palette)
            else:
                self.color = (175, 0, 0)  # Red for HOLO


@dataclass
class PlotConfig:
    """
    Configuration for plot appearance.
    
    Attributes:
        template: Plotly template ('plotly_white', 'ggplot2', etc.)
        style: Visualization style ('simple', 'enhanced', 'publication', 'overview', 'comparative')
        width: Plot width in pixels
        height: Plot height in pixels
        scale: Resolution scale factor
        font_family: Font family name
        font_size: Base font size
        colors: Dictionary of colors for different elements
        save_formats: List of formats to save ('html', 'svg', 'png')
    """
    template: str = 'ggplot2'
    style: str = 'simple'  # simple, enhanced, publication, overview, comparative
    width: int = 1200
    height: int = 600
    scale: int = 2
    font_family: str = 'Times New Roman'
    font_size: int = 24
    colors: Dict[str, Tuple[int, int, int]] = field(default_factory=lambda: {
        'holo': (175, 0, 0),
        'apo': (0, 125, 0)
    })
    save_formats: List[str] = field(default_factory=lambda: ['html', 'svg'])
    
    def __post_init__(self):
        """Validate configuration"""
        valid_styles = ['simple', 'enhanced', 'publication', 'overview', 'comparative']
        if self.style not in valid_styles:
            raise ValueError(f"Invalid style: {self.style}. Choose from {valid_styles}")


@dataclass
class StabilityConfig:
    """Configuration for auto-detecting equilibration/stable window."""
    enabled: bool = False
    metric: str = "rmsd_prot"
    window: int = 50
    std_threshold: float = 0.2
    slope_threshold: float = 0.01
    min_points: int = 100
    apply_to_all_metrics: bool = True


@dataclass
class MDConfig:
    """
    Complete configuration for MD trajectory analysis.
    
    Attributes:
        base_dir: Base directory containing data
        output_dir: Output directory for results
        protein_name: Name of the protein
        systems: List of SystemConfig objects
        metrics: List of AnalysisMetric objects to calculate
        plot_config: PlotConfig object
        residue_range: Optional residue range for RMSF (start, end)
    amino_acids: Optional list of residue numbers or labels for RMSF x-axis
    """
    base_dir: Path
    output_dir: Path
    protein_name: str
    systems: List[SystemConfig]
    metrics: List[AnalysisMetric] = field(default_factory=list)
    plot_config: PlotConfig = field(default_factory=PlotConfig)
    residue_range: Optional[Tuple[int, int]] = None
    amino_acids: Optional[List[Union[int, str]]] = None
    stability: StabilityConfig = field(default_factory=StabilityConfig)
    sequence_topology: Optional[Path] = None
    sequence_selection: str = "protein"
    time_unit_input: str = "auto"   # auto | frame | ps | ns
    time_unit_output: str = "ns"    # auto | frame | ps | ns
    time_step_ps: Optional[float] = None
    time_scale: Optional[float] = None
    smoothing_window: int = 1
    # |z| threshold for masking outliers; None/0 disables (default). Large excursions in
    # MD observables are usually real events, so only enable this to drop known artefacts.
    outlier_threshold: Optional[float] = None
    
    def __post_init__(self):
        """Convert paths to Path objects and set default metrics."""
        self.base_dir = Path(self.base_dir)
        self.output_dir = Path(self.output_dir)
        if self.sequence_topology:
            self.sequence_topology = Path(self.sequence_topology)
        valid_units = {"auto", "frame", "ps", "ns"}
        if self.time_unit_input not in valid_units:
            raise ValueError(f"Invalid time_unit_input: {self.time_unit_input}. Choose from {sorted(valid_units)}")
        if self.time_unit_output not in valid_units:
            raise ValueError(f"Invalid time_unit_output: {self.time_unit_output}. Choose from {sorted(valid_units)}")
        if self.time_step_ps is not None:
            self.time_step_ps = float(self.time_step_ps)
            if self.time_step_ps <= 0:
                raise ValueError("time_step_ps must be > 0 when provided")
        if self.time_scale is not None:
            self.time_scale = float(self.time_scale)
        self.smoothing_window = max(1, int(self.smoothing_window))
        if self.outlier_threshold is not None:
            self.outlier_threshold = float(self.outlier_threshold)
            if self.outlier_threshold <= 0:
                self.outlier_threshold = None
        
        # Set default metrics if none provided
        if not self.metrics:
            self.metrics = self.get_default_metrics()
    
    @staticmethod
    def get_default_metrics() -> List[AnalysisMetric]:
        """Get default set of analysis metrics."""
        return [
            AnalysisMetric(
                name='rmsd_prot',
                file_pattern='RMSD_apo*.dat,RMSD_protein*.dat',
                title='Protein RMSD',
                ylabel='Å',
                is_holo_only=False,
                data_format='single_column'
            ),
            AnalysisMetric(
                name='rmsd_complex',
                file_pattern='RMSD_complex*.dat,RMSD_*_complex*.dat',
                title='Protein & Ligands complex RMSD',
                ylabel='Å',
                is_holo_only=True,
                data_format='single_column'
            ),
            AnalysisMetric(
                name='rmsd_lig',
                file_pattern='RMSD_.dat,RMSD_lig*.dat,RMSD_ligand*.dat',
                title='Ligands RMSD',
                ylabel='Å',
                is_holo_only=True,
                data_format='single_column'
            ),
            AnalysisMetric(
                name='rmsf',
                file_pattern='RMSF*.dat',
                title='RMSF',
                ylabel='Å',
                is_holo_only=False,
                data_format='two_column'
            ),
            AnalysisMetric(
                name='rog',
                file_pattern='rog*.dat,RoG*.dat',
                title='Radius of Gyration',
                ylabel='Å',
                is_holo_only=False,
                data_format='single_column'
            ),
            AnalysisMetric(
                name='sasa',
                file_pattern='SASA*.dat,sasa*.dat',
                title='Solvent Accessible Surface Area',
                ylabel='Å²',
                is_holo_only=False,
                data_format='single_column'
            ),
            AnalysisMetric(
                name='hbonds',
                file_pattern='hbonds*.dat,Hbond*.dat',
                title='Number of Hydrogen bonds',
                ylabel='Number of hydrogen bonds',
                is_holo_only=False,
                data_format='two_column'
            ),
            AnalysisMetric(
                name='comcom',
                file_pattern='comcom*.dat,COM*.dat',
                title='COM-COM Distance',
                ylabel='Å',
                is_holo_only=True,
                data_format='two_column'
            ),
        ]
    
    @classmethod
    def from_json(cls, json_file: Path) -> 'MDConfig':
        """
        Load configuration from JSON file.
        
        Args:
            json_file: Path to JSON configuration file
            
        Returns:
            MDConfig object
        """
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        # Convert systems from dicts to SystemConfig objects
        systems = [SystemConfig(**s) for s in data.get('systems', [])]
        
        # Convert metrics from dicts to AnalysisMetric objects
        metrics = [AnalysisMetric(**m) for m in data.get('metrics', [])]
        
        # Convert plot_config from dict to PlotConfig object
        plot_config_data = data.get('plot_config', {})
        plot_config = PlotConfig(**plot_config_data) if plot_config_data else PlotConfig()
        
        return cls(
            base_dir=Path(data['base_dir']),
            output_dir=Path(data.get('output_dir', './output')),
            protein_name=data['protein_name'],
            systems=systems,
            metrics=metrics,
            plot_config=plot_config,
            residue_range=tuple(data['residue_range']) if 'residue_range' in data else None,
            amino_acids=data.get('amino_acids'),
            time_unit_input=data.get('time_unit_input', 'auto'),
            time_unit_output=data.get('time_unit_output', 'ns'),
            time_step_ps=data.get('time_step_ps'),
            time_scale=data.get('time_scale'),
            smoothing_window=data.get('smoothing_window', 1),
            outlier_threshold=data.get('outlier_threshold'),
        )
    
    def to_json(self, json_file: Path):
        """
        Save configuration to JSON file.
        
        Args:
            json_file: Path to output JSON file
        """
        data = {
            'base_dir': str(self.base_dir),
            'output_dir': str(self.output_dir),
            'protein_name': self.protein_name,
            'systems': [
                {
                    'name': s.name,
                    'dir_pattern': s.dir_pattern,
                    'dir_names': s.dir_names,
                    'is_apo': s.is_apo,
                    'replicates': s.replicates,
                    'color': s.color
                }
                for s in self.systems
            ],
            'metrics': [
                {
                    'name': m.name,
                    'file_pattern': m.file_pattern,
                    'title': m.title,
                    'ylabel': m.ylabel,
                    'is_holo_only': m.is_holo_only,
                    'data_format': m.data_format
                }
                for m in self.metrics
            ],
            'plot_config': {
                'template': self.plot_config.template,
                'width': self.plot_config.width,
                'height': self.plot_config.height,
                'scale': self.plot_config.scale,
                'font_family': self.plot_config.font_family,
                'font_size': self.plot_config.font_size,
                'colors': self.plot_config.colors,
                'save_formats': self.plot_config.save_formats
            },
            'residue_range': self.residue_range,
            'amino_acids': self.amino_acids,
            'time_unit_input': self.time_unit_input,
            'time_unit_output': self.time_unit_output,
            'time_step_ps': self.time_step_ps,
            'time_scale': self.time_scale,
            'smoothing_window': self.smoothing_window,
            'outlier_threshold': self.outlier_threshold,
        }
        
        with open(json_file, 'w') as f:
            json.dump(data, f, indent=2)


def create_example_config() -> MDConfig:
    """
    Create an example configuration.
    
    Returns:
        Example MDConfig object
    """
    return MDConfig(
        base_dir=Path('./data'),
        output_dir=Path('./output'),
        protein_name='ProteinX',
        systems=[
            SystemConfig(
                name='LigandA',
                dir_pattern='ProteinX_HOLO_{}',
                is_apo=False,
                replicates=3
            ),
            SystemConfig(
                name='APO',
                dir_pattern='ProteinX_APO_{}',
                is_apo=True,
                replicates=3
            )
        ],
        residue_range=(814, 1166),
        amino_acids=list(range(814, 937)) + list(range(994, 1169))
    )
