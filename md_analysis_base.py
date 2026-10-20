#!/usr/bin/env python3
"""
Base classes for enhanced MD analysis with modern visualization capabilities
Provides common functionality for both holo and apo analysis workflows
"""

import sys as _sys
import warnings as _warnings

_DEPRECATION = (
    "md_analysis_base.py is deprecated and unmaintained: it keeps the old z-score outlier removal and lacks the "
    "replicate statistics, unit handling and other fixes of the gromacs_analysis package. "
    "Use `gromacs_analysis.md` (CLI: `gromacs-md`) instead (see USAGE.md)."
)
_warnings.warn(_DEPRECATION, DeprecationWarning, stacklevel=2)
print("WARNING: " + _DEPRECATION, file=_sys.stderr)


import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
from pathlib import Path
import logging
from scipy import stats

logger = logging.getLogger(__name__)

class MDAnalysisConfig:
    """Configuration manager for MD analysis parameters"""
    
    def __init__(self, config_file="md_analysis_config.json"):
        """Load configuration from JSON file with defaults"""
        self.config_file = config_file
        self.config = self._load_config()
        
    def _load_config(self):
        """Load configuration with fallback to defaults"""
        defaults = {
            'visualization': {
                'template': 'plotly_dark',
                'colors': ['#009e73', '#e69f00', '#56b4e9', '#cc79a7', '#f0e442'],
                'title_font_size': 22,
                'axis_font_size': 14,
                'font_family': 'Arial'
            },
            'analysis': {
                'smoothing_window': 8,
                'outlier_threshold': 2.5,
                'time_unit': 'ps'
            },
            'output': {
                'formats': ['html', 'svg', 'png'],
                'dpi': 300,
                'width': 1200,
                'height': 800
            }
        }
        
        try:
            if Path(self.config_file).exists():
                with open(self.config_file, 'r') as f:
                    user_config = json.load(f)
                # Merge user config with defaults
                config = defaults.copy()
                for key, value in user_config.items():
                    if key in config:
                        config[key].update(value)
                    else:
                        config[key] = value
                logger.info(f"Configuration loaded from {self.config_file}")
                return config
            else:
                logger.info(f"Config file {self.config_file} not found, using defaults")
                return defaults
        except Exception as e:
            logger.warning(f"Error loading config: {e}, using defaults")
            return defaults
    
    def get_template(self):
        return self.config['visualization']['template']
    
    def get_colors(self):
        return self.config['visualization']['colors']
    
    def get_font_size(self):
        return self.config['visualization']['title_font_size']
    
    def get_smoothing_window(self):
        return self.config['analysis']['smoothing_window']
    
    def get_outlier_threshold(self):
        return self.config['analysis']['outlier_threshold']

class MDDataProcessor:
    """Enhanced data processing for MD simulation results"""
    
    def __init__(self, config):
        """Initialize with configuration"""
        self.config = config
        
    def process_data_file(self, file_path):
        """Process a single .dat file with enhanced cleaning"""
        try:
            # Read file, skipping Grace comments
            data = []
            time_step = 0  # Auto-generate time if not present
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith(('@', '#')):
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                time_val = float(parts[0])
                                value_val = float(parts[1])
                                data.append([time_val, value_val])
                            except ValueError:
                                continue
                        elif len(parts) == 1:
                            try:
                                value_val = float(parts[0])
                                data.append([time_step, value_val])
                                time_step += 1
                            except ValueError:
                                continue
            
            if not data:
                logger.warning(f"No valid data found in {file_path}")
                return np.array([]), np.array([]), {}
            
            # Convert to numpy arrays
            data_array = np.array(data)
            time_data = data_array[:, 0]
            value_data = data_array[:, 1]
            
            # Remove outliers
            time_clean, value_clean = self._remove_outliers(time_data, value_data)
            
            # Apply smoothing if requested
            smoothing_window = self.config.get_smoothing_window()
            if smoothing_window > 1 and len(value_clean) > smoothing_window:
                value_smooth = self._apply_smoothing(value_clean, smoothing_window)
            else:
                value_smooth = value_clean
            
            # Calculate statistics
            statistics = self._calculate_statistics(value_smooth)
            
            return time_clean, value_smooth, statistics
            
        except Exception as e:
            logger.error(f"Error processing {file_path}: {e}")
            return np.array([]), np.array([]), {}
    
    def _remove_outliers(self, time_data, value_data):
        """Remove outliers using Z-score method"""
        threshold = self.config.get_outlier_threshold()
        
        if len(value_data) < 3:  # Need at least 3 points for Z-score
            return time_data, value_data
        
        z_scores = np.abs(stats.zscore(value_data))
        mask = z_scores < threshold
        
        n_outliers = len(value_data) - np.sum(mask)
        if n_outliers > 0:
            logger.info(f"Removed {n_outliers} outliers (threshold: {threshold})")
        
        return time_data[mask], value_data[mask]
    
    def _apply_smoothing(self, data, window):
        """Apply moving average smoothing"""
        if len(data) < window:
            return data
        
        # Pandas rolling mean
        df = pd.DataFrame({'values': data})
        smoothed = df['values'].rolling(window=window, center=True).mean()
        
        # Fill NaN values at edges
        smoothed = smoothed.bfill().ffill()
        
        return smoothed.values
    
    def _calculate_statistics(self, data):
        """Calculate comprehensive statistics"""
        if len(data) == 0:
            return {}
        
        return {
            'mean': np.mean(data),
            'std': np.std(data),
            'min': np.min(data),
            'max': np.max(data),
            'median': np.median(data),
            'q1': np.percentile(data, 25),
            'q3': np.percentile(data, 75),
            'count': len(data)
        }

class EnhancedMDVisualizer:
    """Enhanced visualization with modern Plotly plots"""
    
    def __init__(self, config):
        """Initialize with configuration"""
        self.config = config
        
    def create_enhanced_plot(self, data_dict, title, ylabel, analysis_type):
        """Create an enhanced interactive plot"""
        if not data_dict:
            return None
        
        fig = go.Figure()
        colors = self.config.get_colors()
        
        color_idx = 0
        for system_name, data in data_dict.items():
            time_data = data['time']
            value_data = data['values']
            label = data['label']
            
            # Main trace
            fig.add_trace(go.Scatter(
                x=time_data,
                y=value_data,
                mode='lines',
                name=label,
                line=dict(
                    color=colors[color_idx % len(colors)],
                    width=2
                ),
                hovertemplate=f'<b>{label}</b><br>' +
                            'Time: %{x:.1f} ps<br>' +
                            f'{ylabel}: %{{y:.3f}}<br>' +
                            '<extra></extra>'
            ))
            
            # Add mean line
            mean_val = np.mean(value_data)
            fig.add_hline(
                y=mean_val,
                line_dash="dash",
                line_color=colors[color_idx % len(colors)],
                opacity=0.7,
                annotation_text=f"Mean: {mean_val:.3f}",
                annotation_position="top right"
            )
            
            color_idx += 1
        
        # Update layout
        fig.update_layout(
            title=dict(
                text=title,
                font=dict(size=self.config.get_font_size())
            ),
            xaxis_title="Time (ps)",
            yaxis_title=ylabel,
            template=self.config.get_template(),
            hovermode='x unified',
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99
            )
        )
        
        return fig

# Print success message when base is imported
print("Enhanced MD Analysis base classes created successfully!")
print("Key features:")
print("- Configurable visualization parameters")
print("- Enhanced data cleaning with outlier detection") 
print("- Statistical analysis capabilities")
print("- Modern plotly visualizations")
print("- Modular design for easy customization")
