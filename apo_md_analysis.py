#!/usr/bin/env python3
"""
Enhanced Apo MD Analysis with modern visualization capabilities
Processes .dat files from GROMACS simulations for protein-only systems
"""

import sys as _sys
import warnings as _warnings

_DEPRECATION = (
    "apo_md_analysis.py is deprecated and unmaintained: it keeps the old z-score outlier removal and lacks the "
    "replicate statistics, unit handling and other fixes of the gromacs_analysis package. "
    "Use `gromacs-md` or `gromacs-pipeline` instead (see USAGE.md)."
)
_warnings.warn(_DEPRECATION, DeprecationWarning, stacklevel=2)
print("WARNING: " + _DEPRECATION, file=_sys.stderr)


import argparse
import logging
import json
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
from md_analysis_base import MDAnalysisConfig, MDDataProcessor, EnhancedMDVisualizer

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ApoMDAnalyzer(MDDataProcessor, EnhancedMDVisualizer):
    """Enhanced analyzer for apo (protein-only) MD simulations"""
    
    def __init__(self, config_file="md_analysis_config.json"):
        """Initialize with configuration"""
        self.config = MDAnalysisConfig(config_file)
        MDDataProcessor.__init__(self, self.config)
        EnhancedMDVisualizer.__init__(self, self.config)
        
    def setup_analysis(self, data_directory, protein_name, residue_range):
        """Setup analysis parameters and create output directory"""
        self.data_dir = Path(data_directory)
        self.output_dir = Path(data_directory) / "enhanced_analysis_apo"
        self.plots_dir = self.output_dir / "plots"
        
        # Create directories with parents=True
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        
        self.protein_name = protein_name
        self.residue_range = residue_range
        
        # Setup analysis types for apo systems
        self.analysis_types = {
            'rmsd': {'pattern': '*apo*.dat', 'ylabel': 'RMSD (nm)'},
            'hbonds': {'pattern': 'hbonds*apo*.dat', 'ylabel': 'H-bonds'},
            'rog': {'pattern': 'rog*apo*.dat', 'ylabel': 'Radius of Gyration (nm)'},
            'sasa': {'pattern': 'SASA*apo*.dat', 'ylabel': 'SASA (nm²)'},
            'rmsf': {'pattern': 'RMSF*apo*.dat', 'ylabel': 'RMSF (nm)'}
        }
        
    def find_data_files(self):
        """Find all relevant .dat files for apo analysis"""
        found_files = {}
        
        for analysis_type, info in self.analysis_types.items():
            pattern = info['pattern']
            files = list(self.data_dir.glob(pattern))
            
            if files:
                found_files[analysis_type] = files[0]  # Take first match
                logger.info(f"Found {analysis_type} file: {files[0]}")
            else:
                logger.warning(f"Missing {analysis_type} file with pattern: {pattern}")
        
        return found_files
    
    def process_all_data(self, found_files):
        """Process all found data files"""
        processed_data = {}
        all_statistics = {}
        
        for analysis_type, file_path in found_files.items():
            ylabel = self.analysis_types[analysis_type]['ylabel']
            
            try:
                # Process data
                time_data, value_data, stats = self.process_data_file(file_path)
                
                if len(time_data) > 0 and len(value_data) > 0:
                    processed_data[analysis_type] = {
                        'apo': {
                            'time': time_data,
                            'values': value_data,
                            'label': f"{self.protein_name} (Apo)"
                        }
                    }
                    all_statistics[analysis_type] = {'apo': stats}
                    logger.info(f"Processed {analysis_type} data: {len(time_data)} points")
                else:
                    logger.warning(f"No valid data in {file_path}")
                    
            except Exception as e:
                logger.error(f"Error processing {file_path}: {e}")
                continue
        
        return processed_data, all_statistics
    
    def create_visualizations(self, processed_data, all_statistics):
        """Create enhanced visualizations for apo analysis"""
        plots_created = []
        
        for analysis_type, data_dict in processed_data.items():
            if not data_dict:
                continue
                
            ylabel = self.analysis_types[analysis_type]['ylabel']
            title = f"{self.protein_name} (Apo) - {analysis_type.upper()} Analysis"
            
            # Create the plot
            fig = self.create_enhanced_plot(
                data_dict, 
                title=title,
                ylabel=ylabel,
                analysis_type=analysis_type
            )
            
            if fig:
                # Save in multiple formats
                output_base = self.plots_dir / f"{analysis_type}_apo_analysis"
                
                # Save HTML (interactive)
                fig.write_html(f"{output_base}.html")
                
                # Save static formats
                try:
                    fig.write_image(f"{output_base}.png", width=1200, height=800, scale=2)
                    fig.write_image(f"{output_base}.svg", width=1200, height=800)
                except Exception as e:
                    logger.warning(f"Could not save static images: {e}")
                
                plots_created.append(f"{analysis_type}_apo_analysis")
                logger.info(f"Created visualization: {analysis_type}")
        
        # Create overview plot if we have key data
        self.create_overview_plot(processed_data)
        
        return plots_created
    
    def create_overview_plot(self, processed_data):
        """Create multi-panel overview plot with key metrics"""
        key_analyses = ['rmsd', 'hbonds', 'rog', 'sasa']
        available_analyses = [anal for anal in key_analyses if anal in processed_data and processed_data[anal]]
        
        if len(available_analyses) < 2:
            logger.info("Not enough data for overview plot")
            return None
        
        n_plots = len(available_analyses)
        rows = 2
        cols = 2 if n_plots <= 4 else 3
        
        fig = make_subplots(
            rows=rows, cols=cols,
            subplot_titles=[f"{anal.upper()}" for anal in available_analyses[:4]],
            vertical_spacing=0.1,
            horizontal_spacing=0.1
        )
        
        colors = self.config.get_colors()
        
        for i, analysis_type in enumerate(available_analyses[:4]):
            row = (i // cols) + 1
            col = (i % cols) + 1
            
            data_dict = processed_data[analysis_type]
            ylabel = self.analysis_types[analysis_type]['ylabel']
            
            # Add apo data
            apo_data = data_dict['apo']
            fig.add_trace(
                go.Scatter(
                    x=apo_data['time'],
                    y=apo_data['values'],
                    name=f"{analysis_type} - Apo",
                    line=dict(color=colors[0], width=2),
                    showlegend=(i == 0)  # Only show legend for first subplot
                ),
                row=row, col=col
            )
            
            # Update axes labels
            fig.update_xaxes(title_text="Time (ps)" if row == rows else "", row=row, col=col)
            fig.update_yaxes(title_text=ylabel, row=row, col=col)
        
        # Update layout
        fig.update_layout(
            title=f"{self.protein_name} (Apo) - MD Analysis Overview",
            height=800,
            template=self.config.get_template(),
            font=dict(size=self.config.get_font_size())
        )
        
        # Save overview plot
        output_base = self.plots_dir / "overview_apo_analysis"
        fig.write_html(f"{output_base}.html")
        
        try:
            fig.write_image(f"{output_base}.png", width=1400, height=800, scale=2)
            fig.write_image(f"{output_base}.svg", width=1400, height=800)
        except Exception as e:
            logger.warning(f"Could not save overview static images: {e}")
        
        logger.info("Created apo overview visualization")
        return fig
    
    def save_statistics(self, all_statistics):
        """Save comprehensive statistics to CSV"""
        stats_rows = []
        
        for analysis_type, stats_dict in all_statistics.items():
            for condition, stats in stats_dict.items():
                row = {
                    'Analysis_Type': analysis_type,
                    'Condition': condition,
                    'System': f"{self.protein_name} ({condition})",
                    'Mean': stats.get('mean', 0),
                    'Std': stats.get('std', 0),
                    'Min': stats.get('min', 0),
                    'Max': stats.get('max', 0),
                    'Median': stats.get('median', 0),
                    'Q1': stats.get('q1', 0),
                    'Q3': stats.get('q3', 0),
                    'Count': stats.get('count', 0)
                }
                stats_rows.append(row)
        
        if stats_rows:
            df_stats = pd.DataFrame(stats_rows)
            stats_file = self.output_dir / "summary_statistics_apo.csv"
            df_stats.to_csv(stats_file, index=False)
            logger.info(f"Apo statistics saved to {stats_file}")
            return stats_file
        
        return None
    
    def run_analysis(self, data_directory, protein_name, residue_range=None):
        """Run complete apo MD analysis"""
        logger.info("Starting Apo MD Analysis...")
        
        # Setup
        self.setup_analysis(data_directory, protein_name, residue_range)
        
        # Find data files
        found_files = self.find_data_files()
        if not found_files:
            raise ValueError("No apo .dat files found matching expected patterns")
        
        # Process data
        processed_data, all_statistics = self.process_all_data(found_files)
        
        # Create visualizations
        plots_created = self.create_visualizations(processed_data, all_statistics)
        
        # Save statistics
        stats_file = self.save_statistics(all_statistics)
        
        results = {
            'output_directory': self.output_dir,
            'plots_directory': self.plots_dir,
            'plots_created': plots_created,
            'statistics_file': stats_file,
            'processed_files': len(found_files)
        }
        
        logger.info(f"Apo analysis complete! Results saved to {self.output_dir}")
        logger.info(f"Created {len(plots_created)} visualizations")
        
        return results, all_statistics

def main():
    """Main function with command line interface"""
    parser = argparse.ArgumentParser(description="Enhanced Apo MD Analysis")
    parser.add_argument("--config", default="md_analysis_config.json", help="Configuration file")
    parser.add_argument("--data_dir", default=".", help="Directory containing .dat files")
    parser.add_argument("--protein", required=True, help="Protein name")
    parser.add_argument("--res_range", nargs=2, type=int, help="Residue range (start end)")
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = ApoMDAnalyzer(args.config)
    
    # Run analysis
    results, statistics = analyzer.run_analysis(
        args.data_dir, 
        args.protein, 
        args.res_range
    )
    
    print(f"\n{'='*50}")
    print("APO MD ANALYSIS COMPLETE")
    print(f"{'='*50}")
    print(f"Output directory: {results['output_directory']}")
    print(f"Plots created: {len(results['plots_created'])}")
    print(f"Files processed: {results['processed_files']}")
    print(f"Statistics saved: {results['statistics_file']}")
    print("\nCheck the enhanced_analysis_apo/plots/ directory for visualizations!")

if __name__ == "__main__":
    main()
