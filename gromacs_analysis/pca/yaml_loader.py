"""
YAML Configuration Loader for PCA
==================================

Load PCA configuration from YAML files.
"""

import logging
from pathlib import Path
from typing import List
import yaml

from .config import (
    PCAConfig,
    PCASystemConfig,
    FELConfig,
    ClusteringConfig,
    PCAGromacsRunConfig,
    PCATerminalDetectionConfig,
)
from ..md.config import PlotConfig

logger = logging.getLogger(__name__)


def parse_run_gromacs_config(data: dict) -> PCAGromacsRunConfig:
    """Parse optional automatic GROMACS PCA generation config."""
    run_data = data.get('run_gromacs', {}) or {}
    terminal_data = run_data.get('terminal_detection', {}) or {}
    core_range = terminal_data.get('core_residue_range')
    if core_range:
        core_range = tuple(core_range)

    terminal_config = PCATerminalDetectionConfig(
        enabled=terminal_data.get('enabled', True),
        rmsf_file=terminal_data.get('rmsf_file', 'RMSF_.dat'),
        pdb_file=terminal_data.get('pdb_file', 'step3_input.pdb'),
        gro_file=terminal_data.get('gro_file', 'step3_input.gro'),
        smoothing_window=terminal_data.get('smoothing_window', 3),
        mad_multiplier=terminal_data.get('mad_multiplier', 3.0),
        stable_run=terminal_data.get('stable_run', 3),
        max_terminal_fraction=terminal_data.get('max_terminal_fraction', 0.25),
        core_residue_range=core_range,
    )

    window_data = run_data.get('window_cosine', {}) or {}
    return PCAGromacsRunConfig(
        enabled=bool(run_data.get('enabled', False)),
        input_dir=Path(run_data['input_dir']) if run_data.get('input_dir') else None,
        gmx=run_data.get('gmx', 'gmx'),
        tpr=run_data.get('tpr', 'step5_production.tpr'),
        trajectory=run_data.get('trajectory', 'step5_production.xtc'),
        no_pbc_trajectory=run_data.get('no_pbc_trajectory', 'step5_production_noPBC.xtc'),
        index=run_data.get('index', 'new_index.ndx'),
        protein_group=run_data.get('protein_group', 'Protein'),
        output_group=run_data.get('output_group', 'System'),
        core_group_name=run_data.get('core_group_name', 'core_no_terminals'),
        calpha_group_name=run_data.get('calpha_group_name', 'core_no_terminals_Calpha'),
        n_pcs=run_data.get('n_pcs', data.get('n_pcs', 10)),
        window_cosine_enabled=bool(window_data.get('enabled', True)),
        window_max_ns=window_data.get('max_ns', 100),
        window_step_ns=window_data.get('step_ns', 10),
        extreme_nframes=run_data.get('extreme_nframes', 30),
        overwrite=bool(run_data.get('overwrite', True)),
        terminal_detection=terminal_config,
    )


def load_yaml_config(yaml_path: Path) -> PCAConfig:
    """
    Load PCA configuration from YAML file.
    
    Args:
        yaml_path: Path to YAML configuration file
        
    Returns:
        PCAConfig object
    """
    try:
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)
        
        # Parse base directory
        base_dir = Path(data['base_dir'])
        
        # Parse output directory (optional)
        output_dir = data.get('output_dir')
        if output_dir:
            output_dir = Path(output_dir)
        else:
            output_dir = base_dir / 'PCA_Vis'
        
        # Parse protein name
        protein_name = data.get('protein_name', 'Protein')
        
        # Parse ligand name (explicit, overrides auto-detection)
        ligand_name = data.get('ligand_name')
        
        # Parse APO base directory (optional, for shared APO runs)
        apo_base_dir = data.get('apo_base_dir')
        if apo_base_dir:
            apo_base_dir = Path(apo_base_dir)
        
        # Parse number of PCs
        n_pcs = data.get('n_pcs', 10)
        
        # Parse systems (optional, for future use)
        systems = []
        if 'systems' in data:
            for sys_data in data['systems']:
                system = PCASystemConfig(
                    name=sys_data['name'],
                    is_apo=sys_data.get('is_apo', False),
                    ligand_name=sys_data.get('ligand_name')
                )
                systems.append(system)
        
        # Parse plot config
        plot_data = data.get('plot_config', {})
        plot_config = PlotConfig(
            template=plot_data.get('template', 'plotly_white'),
            width=plot_data.get('width', 1200),
            height=plot_data.get('height', 600),
            scale=plot_data.get('scale', 2),
            font_family=plot_data.get('font_family', 'Times New Roman'),
            font_size=plot_data.get('font_size', 24),
            save_formats=plot_data.get('save_formats', ['html', 'svg'])
        )
        
        fel_data = data.get('fel', {})
        fel_config = FELConfig(
            enabled=fel_data.get('enabled', False),
            pairs=fel_data.get('pairs', ["12"]),
            contours=fel_data.get('contours', 7),
            step=fel_data.get('step', 2),
            colorscale=fel_data.get('colorscale', 'Viridis'),
            overlay_projection=fel_data.get('overlay_projection', True),
            include_probability=fel_data.get('include_probability', False),
            include_entropy=fel_data.get('include_entropy', False),
            include_enthalpy=fel_data.get('include_enthalpy', False),
        )

        cluster_data = data.get('clustering', {})
        clustering_config = ClusteringConfig(
            enabled=cluster_data.get('enabled', False),
            method=cluster_data.get('method', 'kmeans'),
            pair=cluster_data.get('pair', '12'),
            n_clusters=cluster_data.get('n_clusters', 3),
            eps=cluster_data.get('eps', 1.5),
            min_samples=cluster_data.get('min_samples', 5),
            downsample=cluster_data.get('downsample', 1),
            max_points=cluster_data.get('max_points'),
        )

        run_gromacs_config = parse_run_gromacs_config(data)

        config = PCAConfig(
            base_dir=base_dir,
            output_dir=output_dir,
            protein_name=protein_name,
            ligand_name=ligand_name,
            apo_base_dir=apo_base_dir,
            systems=systems,
            n_pcs=n_pcs,
            plot_config=plot_config,
            fel=fel_config,
            clustering=clustering_config,
            run_gromacs=run_gromacs_config,
        )
        
        logger.info(f"Loaded PCA configuration from {yaml_path}")
        return config
        
    except Exception as e:
        logger.error(f"Error loading YAML config from {yaml_path}: {e}")
        raise


def generate_yaml_template(output_path: Path):
    """
    Generate a template YAML configuration file.
    
    Args:
        output_path: Path where template will be saved
    """
    template = """# PCA Analysis Configuration Template
# ============================================

# Base directory containing HOLO PCA files (ligand-specific)
# Example: '/path/to/holo_pca'
base_dir: "/path/to/2-PCA"

# Output directory (optional, defaults to base_dir/PCA_Vis)
# output_dir: "/path/to/PCA_Vis"

# Protein name
protein_name: "ProteinX"

# Ligand name (explicit, overrides auto-detection)
# If not specified, will be auto-detected from filenames
ligand_name: "LigandA"

# APO base directory (optional, for shared APO runs)
# If multiple ligands share the same APO data, specify the APO directory here
# If not specified, uses base_dir for APO files
# apo_base_dir: "/path/to/APO/2-PCA"

# Number of principal components to analyze (default: 10)
n_pcs: 10

# Optional: generate PCA files by running GROMACS before plotting.
# run_gromacs:
#   enabled: true
#   input_dir: "/path/to/simulation_folder"
#   gmx: "gmx"
#   tpr: "step5_production.tpr"
#   trajectory: "step5_production.xtc"
#   index: "new_index.ndx"
#   protein_group: "Protein"
#   output_group: "System"
#   core_group_name: "core_no_terminals"
#   calpha_group_name: "core_no_terminals_Calpha"
#   terminal_detection:
#     rmsf_file: "RMSF_.dat"
#     gro_file: "step3_input.gro"
#     pdb_file: "step3_input.pdb"
#     smoothing_window: 3
#     mad_multiplier: 3.0
#     stable_run: 3
#     max_terminal_fraction: 0.25
#     # Optional manual override in RMSF/GRO row numbering:
#     # core_residue_range: [9, 289]
#   window_cosine:
#     enabled: true
#     max_ns: 100
#     step_ns: 10

# Plot configuration
plot_config:
  template: "plotly_white"        # Plotly template
  width: 1200                     # Plot width
  height: 600                     # Plot height
  scale: 2                        # Scale factor for images
  font_family: "Times New Roman"  # Font family
  font_size: 24                   # Font size
  save_formats: ["html", "svg"]   # Output formats

# Free Energy Landscape (FEL) configuration (requires gmx sham outputs)
# fel:
#   enabled: true
#   pairs: ["12", "13", "23"]      # PC pairs to plot
#   contours: 7
#   step: 2                        # Downsample overlay points
#   colorscale: "Viridis"
#   overlay_projection: true
#   include_probability: true
#   include_entropy: false
#   include_enthalpy: false

# Clustering on PCA projections
# clustering:
#   enabled: true
#   method: "kmeans"               # kmeans | dbscan
#   pair: "12"                     # PC pair for clustering
#   n_clusters: 3
#   downsample: 1
#   max_points: 5000

# Systems (optional, for future use)
# systems:
#   - name: "APO"
#     is_apo: true
#   - name: "LigandA"
#     is_apo: false
#     ligand_name: "LigandA"
"""
    
    output_path.write_text(template)
    logger.info(f"Template YAML configuration file created: {output_path}")
    logger.info(f"  Edit this file with your settings and run: gromacs-pca --config {output_path.name}")
