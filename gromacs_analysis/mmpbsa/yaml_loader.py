"""
YAML Configuration Loader for MMPBSA
====================================

Load and save MMPBSA configurations from YAML files.
"""

import logging
from pathlib import Path
from typing import Dict, List, Any
import yaml

from .config import MMPBSAConfig, MMPBSASystemConfig
from ..md.config import PlotConfig
from ..utils.publication import PublicationConfig

logger = logging.getLogger(__name__)


def load_yaml_config(yaml_path: Path) -> MMPBSAConfig:
    """
    Load MMPBSA configuration from YAML file.
    
    Args:
        yaml_path: Path to YAML configuration file
        
    Returns:
        MMPBSAConfig object
    """
    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)
    
    # Parse systems
    systems = []
    for system_data in data.get('systems', []):
        systems.append(MMPBSASystemConfig(
            name=system_data['name'],
            dir_pattern=system_data['dir_pattern'],
            replicates=system_data.get('replicates', 1),
            results_file_pattern=system_data.get('results_file_pattern', 'FINAL_RESULTS_MMPBSA*.dat'),
            decomp_file_pattern=system_data.get('decomp_file_pattern', 'FINAL_DECOMP_MMPBSA*.dat')
        ))
    
    # Parse plot config
    plot_data = data.get('plot_config', {})
    plot_config = PlotConfig(
        template=plot_data.get('template', 'ggplot2'),
        style=plot_data.get('style', 'publication'),
        width=plot_data.get('width', 1400),
        height=plot_data.get('height', 700),
        scale=plot_data.get('scale', 2),
        font_family=plot_data.get('font_family', 'Times New Roman'),
        font_size=plot_data.get('font_size', 24),
        colors=plot_data.get('colors', {}),
        save_formats=plot_data.get('save_formats', ['html', 'svg'])
    )
    
    # Parse amino acids if provided
    amino_acids = None
    if 'amino_acids' in data:
        aa_data = data['amino_acids']
        if 'range' in aa_data:
            start, end = aa_data['range']
            amino_acids = list(range(start, end + 1))
        elif 'ranges' in aa_data:
            amino_acids = []
            for start, end in aa_data['ranges']:
                amino_acids.extend(list(range(start, end + 1)))
        elif 'custom' in aa_data:
            amino_acids = aa_data['custom']
    
    # Create config
    config = MMPBSAConfig(
        base_dir=Path(data['base_dir']),
        output_dir=Path(data['output_dir']),
        protein_name=data['protein_name'],
        ligand_name=data['ligand_name'],
        systems=systems,
        file_format=data.get('file_format', 'auto'),
        plot_config=plot_config,
        amino_acids=amino_acids,
        compare_systems=bool(data.get('compare_systems', False)),
        compare_binding_key=data.get('compare_binding_key', 'TOTAL'),
        compare_components=bool(data.get('compare_components', False)),
        compare_component_order=data.get('compare_component_order'),
        publication=PublicationConfig(**(data.get('publication') or {})),
        start_frame=data.get('start_frame'),
        incomplete_replicates_policy=data.get('incomplete_replicates_policy', 'warn'),
        require_decomp_replicates=bool(data.get('require_decomp_replicates', False)),
    )
    
    return config


def generate_yaml_template(output_path: Path):
    """
    Generate a YAML configuration template for MMPBSA.
    
    Args:
        output_path: Path where to save the template
    """
    template = """# MMPBSA Analysis Configuration Template
# ============================================

protein_name: "ProteinX"
ligand_name: "LigandA"
base_dir: "/path/to/data"
output_dir: "/path/to/output"

# Systems to analyze
systems:
  - name: "LigandA"                              # System name (ligand name)
    dir_pattern: "3-Results_Holo/r/HOLO_{}"  # Directory pattern ({} = replicate number)
    replicates: 3                             # Number of replicates (1 for single)
    results_file_pattern: "FINAL_RESULTS_MMPBSA*.dat"  # Pattern for results file
    decomp_file_pattern: "FINAL_DECOMP_MMPBSA*.dat"     # Pattern for decomposition file

# File format (auto-detect if not specified)
file_format: "auto"  # Options: "csv", "dat", "auto"

# Cross-ligand comparison (optional)
compare_systems: false
compare_binding_key: "TOTAL"  # TOTAL, GGAS, or GSOLV
compare_components: false
# compare_component_order: ["ΔVDWAALS", "ΔEEL", "ΔEGB", "ΔESURF", "ΔGGAS", "ΔGSOLV"]

# Amino acid numbering for decomposition (optional)
amino_acids:
  # Option 1: Simple range
  range: [814, 1166]
  
  # Option 2: Multiple ranges with gaps (comment out 'range' above)
  # ranges: [[814, 936], [994, 1168]]
  
  # Option 3: Custom list (comment out 'range' above)
  # custom: [814, 815, 816, 994, 995, 996]

# Additional figures and numerical sources below output_dir/publication
publication:
  enabled: true
  subdirectory: publication
  formats: [html, svg, png, pdf]
  scree_components: 10
  moving_window: 11  # positive odd saved-sample count
  top_n: 15         # positive integer or all
  pca_manifest: publication_inputs.json
  sampling_metadata: {}  # validated frame/time mapping; see USAGE.md

# Plot configuration
plot_config:
  template: "ggplot2"               # Grey background (like old codes)
  style: "publication"              # High-quality publication style
  width: 1400
  height: 700
  scale: 2                          # High DPI for publication
  font_family: "Times New Roman"
  font_size: 24
  save_formats: ["html", "svg"]        # Add "png" only when local Chrome/Kaleido export is verified
"""
    
    output_path.write_text(template)
    logger.info(f"Template YAML configuration file created: {output_path}")
    logger.info(f"  Edit this file with your settings and run: gromacs-mmpbsa --config {output_path.name}")
