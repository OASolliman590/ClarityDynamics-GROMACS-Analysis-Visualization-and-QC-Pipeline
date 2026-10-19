"""
YAML Configuration Loader
=========================

Load and save MDConfig from/to YAML files.
"""

import logging
import re
from pathlib import Path
from typing import List, Optional, Union
import yaml

from ..md.config import MDConfig, SystemConfig, AnalysisMetric, PlotConfig, StabilityConfig
from ..utils.helpers import parse_amino_acid_range

logger = logging.getLogger(__name__)


def parse_amino_acids_config(aa_config: Optional[dict]) -> Optional[List[Union[int, str]]]:
    """
    Parse amino acid configuration for RMSF x-axis labels.

    Supported formats:
      - custom: [814, 815, 816] or ["ALA814", "GLY815"]
      - range: [814, 1166]
      - ranges: [[814, 936], [994, 1168]]
      - labels: ["ALA814", "GLY815"] (used as-is)
      - sequence: "MKT..." or ["ALA", "GLY", "SER"]
        with optional sequence_start to add numbering.
    """
    if not aa_config:
        return None

    if "labels" in aa_config:
        return [str(label) for label in aa_config["labels"]]

    if "sequence" in aa_config:
        sequence = aa_config["sequence"]
        tokens: List[str] = []
        if isinstance(sequence, list):
            tokens = [str(token) for token in sequence if str(token).strip()]
        elif isinstance(sequence, str):
            stripped = sequence.strip()
            if re.search(r"[,\s]", stripped):
                tokens = [t for t in re.split(r"[,\s]+", stripped) if t]
            else:
                tokens = list(stripped)

        tokens = [t.upper() for t in tokens if t]
        start_index = aa_config.get("sequence_start")
        if start_index is None:
            start_index = aa_config.get("start")

        if start_index is not None:
            return [f"{token}{start_index + i}" for i, token in enumerate(tokens)]

        return tokens

    if "custom" in aa_config:
        return aa_config["custom"]
    if "range" in aa_config:
        start, end = aa_config["range"]
        return list(range(start, end + 1))
    if "ranges" in aa_config:
        amino_acids: List[int] = []
        for start, end in aa_config["ranges"]:
            amino_acids.extend(range(start, end + 1))
        return amino_acids

    return None


def load_yaml_config(yaml_file: Path) -> MDConfig:
    """
    Load configuration from YAML file.
    
    Args:
        yaml_file: Path to YAML configuration file
        
    Returns:
        MDConfig object
        
    Raises:
        FileNotFoundError: If YAML file doesn't exist
        ValueError: If YAML structure is invalid
    """
    if not Path(yaml_file).exists():
        raise FileNotFoundError(f"Config file not found: {yaml_file}")
    
    logger.info(f"Loading configuration from: {yaml_file}")
    
    with open(yaml_file, 'r') as f:
        data = yaml.safe_load(f)
    
    # Parse systems
    systems = []
    for sys_data in data.get('systems', []):
        systems.append(SystemConfig(
            name=sys_data['name'],
            dir_pattern=sys_data.get('dir_pattern', ''),
            dir_names=sys_data.get('dir_names'),
            is_apo=sys_data.get('is_apo', False),
            replicates=sys_data.get('replicates', 1),
            color=tuple(sys_data['color']) if 'color' in sys_data else None
        ))
    
    # Parse metrics (use defaults if not specified)
    metrics_data = data.get('metrics', [])
    if metrics_data:
        metrics = []
        for m_data in metrics_data:
            metrics.append(AnalysisMetric(
                name=m_data['name'],
                file_pattern=m_data['file_pattern'],
                title=m_data['title'],
                ylabel=m_data['ylabel'],
                is_holo_only=m_data.get('is_holo_only', False),
                data_format=m_data.get('data_format', 'two_column')
            ))
    else:
        # Use defaults
        metrics = None
    
    # Parse amino acids
    amino_acids = None
    if 'amino_acids' in data:
        amino_acids = parse_amino_acids_config(data['amino_acids'])
    
    # Parse plot config
    plot_data = data.get('plot_config', {})
    plot_config = PlotConfig(
        template=plot_data.get('template', 'plotly_white'),
        style=plot_data.get('style', 'simple'),
        width=plot_data.get('width', 1200),
        height=plot_data.get('height', 600),
        scale=plot_data.get('scale', 2),
        font_family=plot_data.get('font_family', 'Arial'),
        font_size=plot_data.get('font_size', 16),
        colors={k: tuple(v) for k, v in plot_data.get('colors', {}).items()},
        save_formats=plot_data.get('save_formats', ['html', 'svg'])
    )
    
    # Parse residue range
    residue_range = None
    if 'residue_range' in data:
        residue_range = tuple(data['residue_range'])
    
    stability_data = data.get('stability', {})
    stability_config = StabilityConfig(
        enabled=bool(stability_data.get('enabled', False)),
        metric=stability_data.get('metric', 'rmsd_prot'),
        window=stability_data.get('window', 50),
        std_threshold=stability_data.get('std_threshold', 0.2),
        slope_threshold=stability_data.get('slope_threshold', 0.01),
        min_points=stability_data.get('min_points', 100),
        apply_to_all_metrics=bool(stability_data.get('apply_to_all_metrics', True)),
        method=stability_data.get('method', 'chodera'),
        n_candidates=int(stability_data.get('n_candidates', 100)),
    )

    # Create config
    config = MDConfig(
        base_dir=Path(data['base_dir']),
        output_dir=Path(data.get('output_dir', './output')),
        protein_name=data['protein_name'],
        systems=systems,
        metrics=metrics if metrics else [],
        plot_config=plot_config,
        residue_range=residue_range,
        amino_acids=amino_acids,
        stability=stability_config,
        sequence_topology=Path(data["sequence_topology"]) if data.get("sequence_topology") else None,
        sequence_selection=data.get("sequence_selection", "protein"),
        smoothing_window=data.get("smoothing_window", data.get("analysis", {}).get("smoothing_window", 1)),
        outlier_threshold=data.get("outlier_threshold", data.get("analysis", {}).get("outlier_threshold")),
    )
    
    logger.info(f"Configuration loaded successfully")
    logger.info(f"  Protein: {config.protein_name}")
    logger.info(f"  Systems: {len(config.systems)}")
    logger.info(f"  Metrics: {len(config.metrics)}")
    
    return config


def save_config_to_yaml(config: MDConfig, output_file: Path):
    """
    Save MDConfig to YAML file.
    
    Args:
        config: MDConfig object
        output_file: Path to output YAML file
    """
    data = {
        'protein_name': config.protein_name,
        'base_dir': str(config.base_dir),
        'output_dir': str(config.output_dir),
        'systems': [
                {
                    'name': s.name,
                    'dir_pattern': s.dir_pattern,
                    'dir_names': s.dir_names,
                    'is_apo': s.is_apo,
                    'replicates': s.replicates,
                    'color': list(s.color)
                }
            for s in config.systems
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
            for m in config.metrics
        ],
        'plot_config': {
            'template': config.plot_config.template,
            'style': config.plot_config.style,
            'width': config.plot_config.width,
            'height': config.plot_config.height,
            'scale': config.plot_config.scale,
            'font_family': config.plot_config.font_family,
            'font_size': config.plot_config.font_size,
            'colors': {k: list(v) for k, v in config.plot_config.colors.items()},
            'save_formats': config.plot_config.save_formats
        },
        'stability': {
            'enabled': config.stability.enabled,
            'metric': config.stability.metric,
            'window': config.stability.window,
            'std_threshold': config.stability.std_threshold,
            'slope_threshold': config.stability.slope_threshold,
            'min_points': config.stability.min_points,
            'apply_to_all_metrics': config.stability.apply_to_all_metrics,
        },
        'sequence_topology': str(config.sequence_topology) if config.sequence_topology else None,
        'sequence_selection': config.sequence_selection,
        'smoothing_window': config.smoothing_window,
        'outlier_threshold': config.outlier_threshold,
    }
    
    if config.residue_range:
        data['residue_range'] = list(config.residue_range)
    
    if config.amino_acids:
        data['amino_acids'] = {'custom': config.amino_acids}
    
    with open(output_file, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
    
    logger.info(f"Configuration saved to: {output_file}")


def generate_yaml_template(output_file: Path = Path("gromacs_md_config_template.yaml")):
    """
    Generate a template YAML configuration file.
    
    Args:
        output_file: Path to output template file
    """
    template = """# GROMACS MD Analysis Configuration Template
# ============================================

protein_name: "ProteinX"
base_dir: "./data"
output_dir: "./output"

# Systems to analyze
systems:
  - name: "LigandA"                    # System name (e.g., ligand name)
    dir_pattern: "HOLO_LigandA_{}"     # Directory pattern ({} = replicate number)
    # dir_names: ["HOLO_LigandA_1", "HOLO_LigandA_2", "HOLO_LigandA_3"]  # Explicit dirs (override dir_pattern)
    is_apo: false                   # Is this an APO system?
    replicates: 3                   # Number of replicates
    color: [175, 0, 0]              # RGB color for plots

  - name: "APO"
    dir_pattern: "APO_{}"
    is_apo: true
    replicates: 3
    color: [72, 148, 67]            # Green for APO

# Metrics to analyze (optional - will use defaults if not specified)
metrics:
  - name: "rmsd_prot"
    file_pattern: "RMSD_protein*.dat,RMSD_apo*.dat"
    title: "Protein RMSD"
    ylabel: "RMSD (Å)"
    is_holo_only: false
    data_format: "two_column"

  - name: "rmsd_lig"
    file_pattern: "RMSD_[!apo]*.dat,RMSD_ligand*.dat"
    title: "Ligand RMSD"
    ylabel: "RMSD (Å)"
    is_holo_only: true
    data_format: "two_column"

  - name: "rmsf"
    file_pattern: "RMSF*.dat"
    title: "Root Mean Square Fluctuation"
    ylabel: "RMSF (Å)"
    is_holo_only: false
    data_format: "two_column"

  - name: "rog"
    file_pattern: "rog*.dat,RoG*.dat"
    title: "Radius of Gyration"
    ylabel: "Rg (Å)"
    is_holo_only: false
    data_format: "two_column"

  - name: "sasa"
    file_pattern: "SASA*.dat,sasa*.dat"
    title: "Solvent Accessible Surface Area"
    ylabel: "SASA (Å²)"
    is_holo_only: false
    data_format: "two_column"

  - name: "hbonds"
    file_pattern: "hbonds*.dat,Hbond*.dat"
    title: "Hydrogen Bonds"
    ylabel: "Number of H-bonds"
    is_holo_only: false
    data_format: "two_column"

  - name: "comcom"
    file_pattern: "comcom*.dat,COM*.dat"
    title: "COM-COM Distance"
    ylabel: "Distance (Å)"
    is_holo_only: true
    data_format: "two_column"

# Amino acid labeling for RMSF (choose ONE option)
amino_acids:
  # Option 1: Simple range
  range: [814, 1166]
  
  # Option 2: Multiple ranges with gaps (comment out 'range' above)
  # ranges: [[814, 936], [994, 1168]]
  
  # Option 3: Custom list (comment out 'range' above)
  # custom: [814, 815, 816, 994, 995, 996]
  
  # Option 4: Residue labels (use as-is)
  # labels: ["ALA814", "GLY815", "SER816"]
  
  # Option 5: Amino acid sequence (1-letter or 3-letter tokens)
  # sequence: "MKT..."           # contiguous 1-letter string
  # sequence: ["ALA", "GLY", "SER"]  # or token list
  # sequence_start: 814          # optional numbering for labels

# Optional: Residue range (alternative to amino_acids)
# residue_range: [814, 1166]

# Plot configuration
plot_config:
  template: "plotly_white"        # ggplot2, seaborn, simple_white, plotly_white, plotly_dark
  style: "simple"                 # simple, enhanced, publication, overview, comparative
  width: 1200
  height: 600
  scale: 2
  font_family: "Arial"
  font_size: 16
  colors:
    holo: [175, 0, 0]             # Red for HOLO
    apo: [72, 148, 67]            # Green for APO
  save_formats: ["html", "svg"]        # Add "png" only when local Chrome/Kaleido export is verified
"""
    
    with open(output_file, 'w') as f:
        f.write(template)
    
    logger.info(f"Template configuration file created: {output_file}")
    print(f"\n✓ Template YAML configuration file created: {output_file}")
    print(f"  Edit this file with your settings and run: gromacs-md --config {output_file}")
