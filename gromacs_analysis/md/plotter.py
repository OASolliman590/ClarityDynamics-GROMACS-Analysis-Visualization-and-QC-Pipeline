"""
MD Plotter Module
=================

Consolidated plotting functionality with 5 visualization styles:
1. Simple: Lines with error bands (default)
2. Enhanced: With mean lines and annotations
3. Publication: High-quality with professional styling
4. Overview: Multi-panel (2x2 grid)
5. Comparative: With statistical comparison boxes

Combines features from:
- holo_apo_vis.py (simple style, error bands)
- md_analysis_base.py (enhanced style, mean lines)
- MMPBSA_per_res.py (publication style, colors)
- holo_md_analysis.py (overview style, multi-panel)
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
from scipy import stats as scipy_stats

from .config import MDConfig, AnalysisMetric
from ..utils.mpl_fallback import init_matplotlib, rgb_to_mpl
from ..utils.figure_export import write_static
from ..utils.palette import MUTED_HEX, TEXT_HEX, assign_system_colors, dash_for, rgba

logger = logging.getLogger(__name__)


class MDPlotter:
    """
    Creates publication-quality plots for MD trajectory analysis.
    
    Supports 5 visualization styles:
    - simple: Clean lines with error bands
    - enhanced: Adds mean lines and hover info
    - publication: Professional styling for papers
    - overview: Multi-panel summary (2x2 grid)
    - comparative: Statistical comparison between systems
    
    Attributes:
        config: MDConfig object with plot settings
        style: Selected visualization style
    """
    
    def __init__(self, config: MDConfig):
        """
        Initialize plotter with configuration.
        
        Args:
            config: MDConfig object
        """
        self.config = config
        self.style = config.plot_config.style
        names = [system.name for system in config.systems]
        self.system_color_map = assign_system_colors(
            names,
            is_apo={system.name: system.is_apo for system in config.systems},
            explicit={system.name: system.color for system in config.systems if system.color is not None},
        )
        # Beyond eight ligand hues, a line dash is the secondary encoding (hues are never cycled alone).
        holo = [system.name for system in config.systems if not system.is_apo]
        self._system_dash = {name: dash_for(i) for i, name in enumerate(holo)}
        
        # Configure plotly defaults
        pio.defaults.default_format = "svg"
        pio.defaults.default_width = config.plot_config.width
        pio.defaults.default_height = config.plot_config.height
        pio.defaults.default_scale = config.plot_config.scale
        pio.templates.default = config.plot_config.template
        
        logger.info(f"MDPlotter initialized with style: {self.style}")

    @staticmethod
    def _has_values(values) -> bool:
        """Return True when a list-like or array-like object contains values."""
        return values is not None and len(values) > 0

    def _time_axis_label(self) -> str:
        unit = self.config.time_unit_output
        if unit == "auto":
            unit = self.config.time_unit_input if self.config.time_unit_input != "auto" else "time"
        if unit == "frame":
            return "Frame"
        if unit in {"ps", "ns"}:
            return f"Time ({unit})"
        return "Time"
    
    # ==================== MAIN PLOTTING DISPATCHER ====================
    
    def create_plot(self, data: Dict, metric: AnalysisMetric, 
                    x_values: Optional[List] = None) -> go.Figure:
        """
        Create plot based on selected style.
        
        Args:
            data: Dict of {system_name: processed_data_dict}
            metric: AnalysisMetric object
            x_values: Optional custom x-axis values (for RMSF)
            
        Returns:
            Plotly Figure object
        """
        if self.style == 'simple':
            return self._plot_simple(data, metric, x_values)
        elif self.style == 'enhanced':
            return self._plot_enhanced(data, metric, x_values)
        elif self.style == 'publication':
            return self._plot_publication(data, metric, x_values)
        elif self.style == 'comparative':
            return self._plot_comparative(data, metric, x_values)
        else:
            # Default to simple
            return self._plot_simple(data, metric, x_values)
    
    def create_overview_plot(self, all_data: Dict) -> Optional[go.Figure]:
        """
        Create multi-panel overview plot (2x2 grid).
        
        Args:
            all_data: Complete data structure from processor
            
        Returns:
            Plotly Figure or None
        """
        return self._plot_overview(all_data)
    
    # ==================== STYLE 1: SIMPLE (Lines + Error Bands) ====================
    
    def _plot_simple(self, data: Dict, metric: AnalysisMetric,
                     x_values: Optional[List] = None) -> go.Figure:
        """
        Line plots with ±SD bands for replicates.

        Every system is drawn on its own x-axis (time series) so systems with
        different lengths or sampling are not forced onto the first system's
        axis. RMSF profiles share the residue axis and, when enabled, flexible
        terminal overhangs are drawn faded between two boundary lines.
        """
        fig = go.Figure()
        is_profile = metric.base_name == 'rmsf'
        residue_axis = self._residue_axis(data, x_values) if is_profile else None
        terminal = self._rmsf_terminal_region(data) if is_profile else None

        for idx, (system_name, system_data) in enumerate(data.items()):
            color = self._get_system_color(system_name)
            dash = self._system_dash.get(system_name, 'solid')
            display_name = self._format_legend_name(system_name, metric.name)
            opacity = 0.75 if metric.base_name == 'hbonds' else 1.0

            if system_data['mode'] == 'single':
                y_vals = np.asarray(system_data['values'], dtype=float)
                upper = lower = None
            else:
                y_vals = np.asarray(system_data['mean'], dtype=float)
                upper = np.asarray(system_data['upper'], dtype=float)
                lower = np.asarray(system_data['lower'], dtype=float)

            if is_profile:
                x_plot = list(residue_axis[:len(y_vals)])
                if len(x_plot) < len(y_vals):
                    x_plot += [str(i + 1) for i in range(len(x_plot), len(y_vals))]
            else:
                x_plot = np.asarray(system_data.get('time', np.arange(len(y_vals))))
            n = min(len(x_plot), len(y_vals))
            x_plot = list(x_plot[:n])
            y_vals = y_vals[:n]

            if upper is not None and lower is not None and len(upper) and len(lower):
                m = min(n, len(upper), len(lower))
                fig.add_trace(go.Scatter(
                    x=x_plot[:m] + x_plot[:m][::-1],
                    y=list(upper[:m]) + list(lower[:m])[::-1],
                    fill='toself',
                    fillcolor=rgba(color, 0.18),
                    line=dict(color=rgba(color, 0)),
                    showlegend=False,
                    hoverinfo='skip',
                    legendgroup=system_name,
                ))

            line = dict(color=rgba(color, 1.0), width=2, dash=dash)
            if terminal is not None:
                core_start, core_end = terminal['core_start'], terminal['core_end']
                core_y = [y if core_start <= i <= core_end else None for i, y in enumerate(y_vals)]
                fig.add_trace(go.Scatter(
                    x=x_plot, y=core_y, name=display_name, mode='lines', opacity=opacity,
                    line=line, legendgroup=system_name, connectgaps=False,
                ))
                # Terminal overhangs: same colour, faded, joined to the core line.
                for lo, hi in ((0, core_start), (core_end, n - 1)):
                    if hi <= lo:
                        continue
                    fig.add_trace(go.Scatter(
                        x=x_plot[lo:hi + 1], y=list(y_vals[lo:hi + 1]), mode='lines',
                        line=dict(color=rgba(color, 0.3), width=2, dash=dash),
                        showlegend=False, legendgroup=system_name,
                        hovertemplate='%{x}: %{y:.2f} (flexible terminus)<extra>' + display_name + '</extra>',
                    ))
            else:
                fig.add_trace(go.Scatter(
                    x=x_plot, y=list(y_vals), name=display_name, mode='lines',
                    opacity=opacity, line=line, legendgroup=system_name,
                ))

        if terminal is not None:
            self._draw_terminal_markers(fig, terminal, residue_axis)

        self._update_layout_simple(fig, metric, x_values)
        self._layout_legend(fig, len(data))
        return fig

    def _residue_axis(self, data: Dict, x_values: Optional[List]) -> List[str]:
        """Residue labels for the shared RMSF axis (falls back to the data's own residue index)."""
        if self._has_values(x_values):
            return [str(x) for x in (x_values.tolist() if isinstance(x_values, np.ndarray) else x_values)]
        first = next(iter(data.values()))
        idx = first.get('time')
        if idx is None:
            n = len(first.get('values', first.get('mean', [])))
            idx = range(1, n + 1)
        return [str(int(i)) if float(i).is_integer() else str(i) for i in idx]

    def _rmsf_terminal_region(self, data: Dict) -> Optional[Dict]:
        """Detect flexible termini on the per-residue maximum RMSF across systems.

        Uses the same robust rule as the PCA core trimming (median + k·MAD on the
        smoothed profile). Using the maximum over systems gives one conservative
        region that is valid for every system drawn on the shared axis.
        """
        cfg = getattr(self.config, 'rmsf_terminals', None)
        if cfg is None or not cfg.enabled:
            return None
        from ..pca.prep import robust_threshold, smooth_values, _terminal_trim

        profiles = []
        for system_data in data.values():
            vals = system_data.get('values') if system_data.get('mode') == 'single' else system_data.get('mean')
            if vals is not None and len(vals):
                profiles.append(np.asarray(vals, dtype=float))
        if not profiles:
            return None
        n = min(len(p) for p in profiles)
        if n < 4 * max(1, cfg.stable_run):
            return None
        envelope = np.nanmax(np.vstack([p[:n] for p in profiles]), axis=0)
        smoothed = smooth_values(list(envelope), cfg.smoothing_window)
        threshold, median, mad = robust_threshold(smoothed, cfg.mad_multiplier)
        limit = max(1, int(n * cfg.max_terminal_fraction))
        n_trim = _terminal_trim(smoothed, threshold, cfg.stable_run, from_end=False, limit=limit)
        c_trim = _terminal_trim(smoothed, threshold, cfg.stable_run, from_end=True, limit=limit)
        if n_trim == 0 and c_trim == 0:
            return None
        if n_trim + c_trim >= n:
            return None
        return {
            'core_start': n_trim,
            'core_end': n - 1 - c_trim,
            'n_trimmed': n_trim,
            'c_trimmed': c_trim,
            'n_residues': n,
            'threshold': float(threshold),
            'median': float(median),
            'mad': float(mad),
            'mad_multiplier': cfg.mad_multiplier,
            'show_threshold': cfg.show_threshold,
        }

    def _draw_terminal_markers(self, fig: go.Figure, terminal: Dict, residue_axis: List[str]) -> None:
        """Shade terminal overhangs, draw the two core boundaries and explain them."""
        n = terminal['n_residues']
        muted = rgba(MUTED_HEX, 0.9)
        notes = []
        if terminal['n_trimmed']:
            boundary = terminal['core_start'] - 0.5
            fig.add_vrect(x0=-0.5, x1=boundary, fillcolor=rgba(MUTED_HEX, 0.08), line_width=0, layer='below')
            fig.add_vline(x=boundary, line=dict(color=muted, width=1.5, dash='dash'))
            notes.append(f"N-terminus: {terminal['n_trimmed']} residue(s) before {residue_axis[terminal['core_start']]}")
        if terminal['c_trimmed']:
            boundary = terminal['core_end'] + 0.5
            fig.add_vrect(x0=boundary, x1=n - 0.5, fillcolor=rgba(MUTED_HEX, 0.08), line_width=0, layer='below')
            fig.add_vline(x=boundary, line=dict(color=muted, width=1.5, dash='dash'))
            notes.append(f"C-terminus: {terminal['c_trimmed']} residue(s) after {residue_axis[terminal['core_end']]}")
        if terminal['show_threshold']:
            fig.add_hline(
                y=terminal['threshold'],
                line=dict(color=rgba(MUTED_HEX, 0.6), width=1, dash='dot'),
            )
        fig.add_annotation(
            text=(
                "Grey bands, faded lines: flexible terminal overhangs (smoothed RMSF above "
                f"median + {terminal['mad_multiplier']:g}·MAD = {terminal['threshold']:.2f}, dotted line). "
                "Dashed lines mark where the rigid core used for PCA starts and ends.<br>"
                + " · ".join(notes)
            ),
            xref='paper', yref='paper', x=0.0, y=1.0, xanchor='left', yanchor='bottom',
            showarrow=False, align='left',
            font=dict(size=max(10, self.config.plot_config.font_size - 6), color=TEXT_HEX),
        )

    def _layout_legend(self, fig: go.Figure, n_series: int) -> None:
        """Horizontal legend under the plot for few series; a right-hand column for many."""
        if n_series > 4:
            fig.update_layout(legend=dict(orientation='v', x=1.02, xanchor='left', y=1.0, yanchor='top'))
        else:
            fig.update_layout(legend=dict(orientation='h', x=0.5, xanchor='center', y=-0.25, yanchor='top'))

    # ==================== STYLE 2: ENHANCED (Mean Lines + Annotations) ====================
    
    def _plot_enhanced(self, data: Dict, metric: AnalysisMetric,
                       x_values: Optional[List] = None) -> go.Figure:
        """
        Enhanced plots with horizontal mean lines and statistics.
        From md_analysis_base.py create_enhanced_plot().
        """
        fig = self._plot_simple(data, metric, x_values)
        
        # Add horizontal mean lines for each system
        colors = [self._get_system_color(name) for name in data.keys()]
        
        for idx, (system_name, system_data) in enumerate(data.items()):
            color = colors[idx]
            display_name = self._format_legend_name(system_name, metric.name)
            
            # Calculate mean
            if system_data['mode'] == 'single':
                mean_val = np.mean(system_data['values'])
            else:
                mean_val = np.mean(system_data['mean'])
            
            # Add horizontal line
            fig.add_hline(
                y=mean_val,
                line_dash="dash",
                line_color=f'rgba({color[0]},{color[1]},{color[2]},0.7)',
                opacity=0.7,
                annotation_text=f"{display_name} mean: {mean_val:.3f}",
                annotation_position="top right",
                annotation_font_size=10
            )
        
        return fig
    
    # ==================== STYLE 3: PUBLICATION (High Quality) ====================
    
    def _plot_publication(self, data: Dict, metric: AnalysisMetric,
                          x_values: Optional[List] = None) -> go.Figure:
        """
        Publication-quality plots with professional styling.
        Based on MMPBSA_per_res.py styling.
        """
        fig = self._plot_simple(data, metric, x_values)
        
        # Apply publication styling with grey background (matching ggplot2)
        fig.update_layout(
            font=dict(
                family="Times New Roman, serif" if 'Times' in self.config.plot_config.font_family else self.config.plot_config.font_family,
                size=self.config.plot_config.font_size,
                color='black'
            ),
            title_font=dict(
                size=self.config.plot_config.font_size + 4,
                color='black'
            ),
            paper_bgcolor='white',
            plot_bgcolor='#E5ECF6',  # Grey background (ggplot2 style)
            legend=dict(
                bgcolor='rgba(255, 255, 255, 0.9)',
                bordercolor='rgba(0, 0, 0, 0.2)',
                borderwidth=1
            )
        )
        
        # Clean axis styling
        fig.update_xaxes(
            showgrid=True,
            gridcolor='white',  # White grid lines on grey background
            showline=True,
            linecolor='black',
            linewidth=1.5,
            mirror=True
        )
        
        fig.update_yaxes(
            showgrid=True,
            gridcolor='white',  # White grid lines on grey background
            showline=True,
            linecolor='black',
            linewidth=1.5,
            mirror=True,
            zeroline=True,
            zerolinecolor='rgba(0, 0, 0, 0.2)'
        )
        
        return fig
    
    # ==================== STYLE 4: OVERVIEW (Multi-Panel) ====================
    
    def _plot_overview(self, all_data: Dict) -> Optional[go.Figure]:
        """
        Create 2x2 multi-panel overview plot.
        From holo_md_analysis.py create_overview_plot().
        
        Shows key metrics: RMSD, COM-COM, H-bonds, Rg
        """
        # Select key metrics for overview
        key_metrics = ['rmsd_prot', 'comcom', 'hbonds', 'rog']
        
        # Find available metrics
        available_metrics = []
        metric_objects = {}
        
        # Union over systems: the APO system has no ligand metrics, ligands may lack others.
        present = set()
        for system_data in all_data.values():
            present.update(system_data.keys())

        for metric_name in key_metrics:
            if metric_name in present:
                available_metrics.append(metric_name)
                # Find metric object
                for m in self.config.metrics:
                    if m.name == metric_name:
                        metric_objects[metric_name] = m
                        break
        
        if len(available_metrics) < 2:
            logger.warning("Not enough metrics for overview plot")
            return None
        
        # Create subplots
        n_plots = min(len(available_metrics), 4)
        cols = 2
        rows = (n_plots + cols - 1) // cols
        
        fig = make_subplots(
            rows=rows, cols=cols,
            subplot_titles=[metric_objects[m].title for m in available_metrics[:4]],
            vertical_spacing=0.12,
            horizontal_spacing=0.10
        )
        
        # Add traces for each metric; each system gets one legend entry.
        legend_shown = set()
        for i, metric_name in enumerate(available_metrics[:4]):
            if i >= 4:
                break
            
            row = (i // cols) + 1
            col = (i % cols) + 1
            
            metric = metric_objects[metric_name]
            
            # Plot each system
            for system_name, system_data in all_data.items():
                if metric_name not in system_data:
                    continue
                
                data = system_data[metric_name]
                color = self._get_system_color(system_name)
                display_name = self._format_legend_name(system_name, metric_name)
                
                # Get y-values
                if data['mode'] == 'single':
                    y_values = data['values']
                else:
                    y_values = data['mean']
                
                # Generate x-values
                x_values = data.get('time', np.arange(len(y_values)))
                
                fig.add_trace(
                    go.Scatter(
                        x=x_values,
                        y=y_values,
                        name=system_name,
                        legendgroup=system_name,
                        line=dict(
                            color=f'rgb({color[0]},{color[1]},{color[2]})',
                            width=2,
                            dash=self._system_dash.get(system_name, 'solid'),
                        ),
                        showlegend=system_name not in legend_shown,
                    ),
                    row=row, col=col
                )
                legend_shown.add(system_name)
            
            # Update axes
            fig.update_xaxes(title_text=self._time_axis_label() if row == rows else "", row=row, col=col)
            fig.update_yaxes(title_text=metric.ylabel, row=row, col=col)
        
        # Update layout
        fig.update_layout(
            title=f"{self.config.protein_name} - MD Analysis Overview",
            height=420 * rows + 120,
            template=self.config.plot_config.template,
            font=dict(size=self.config.plot_config.font_size - 2)
        )
        
        return fig
    
    # ==================== STYLE 5: COMPARATIVE (Statistical Comparison) ====================
    
    def _plot_comparative(self, data: Dict, metric: AnalysisMetric,
                          x_values: Optional[List] = None) -> go.Figure:
        """
        Comparative plots with statistical comparison boxes.
        Shows mean, std, and if possible, statistical significance.
        """
        fig = self._plot_simple(data, metric, x_values)
        
        # Add statistics annotation box
        stats_text = self._generate_statistics_text(data)
        
        fig.add_annotation(
            text=stats_text,
            xref="paper", yref="paper",
            x=0.98, y=0.98,
            xanchor='right', yanchor='top',
            showarrow=False,
            bordercolor='black',
            borderwidth=1,
            borderpad=10,
            bgcolor='rgba(255, 255, 255, 0.9)',
            font=dict(size=10, family='monospace')
        )
        
        return fig
    
    # ==================== HELPER METHODS ====================
    
    def _get_system_color(self, system_name: str) -> Tuple[int, int, int]:
        """
        Get RGB color for a system (matching old code).
        
        Old code uses: [[175,0,0], [0,175,0]]
        - HOLO: [175, 0, 0] (red)
        - APO: [0, 175, 0] (green)
        
        Args:
            system_name: System name
            
        Returns:
            RGB tuple
        """
        if system_name not in self.system_color_map:
            # Systems not declared in the config (e.g. programmatic use) take the next free slot.
            names = list(self.system_color_map) + [system_name]
            self.system_color_map = assign_system_colors(
                names,
                explicit=dict(self.system_color_map),
            )
        return self.system_color_map[system_name]
    
    def _format_legend_name(self, system_name: str, metric_name: str) -> str:
        """
        Format system name for legend display based on metric type.
        
        Rules:
        - COM-COM: ligand name only (ligand-only figure)
        - Ligand RMSD: ligand name only (ligand-only figure)
        - H-bonds: Protein_name_Ligand_name_Complex and Apo
        - RMSD Complex: Protein_name_Ligand_name_Complex
        - RMSD Protein: Protein_name Backbone (for both HOLO and APO)
        - RMSF: Protein_name_Ligand_name_Complex
        - Rg: Protein_name_Ligand_name_Complex
        - SASA: Protein_name_Ligand_name_Complex
        
        Args:
            system_name: System name from config
            metric_name: Name of the metric being plotted
            
        Returns:
            Formatted display name
        """
        base = self._metric_base_name(metric_name)
        system = next((s for s in self.config.systems if s.name == system_name), None)
        is_apo = system.is_apo if system is not None else 'APO' in system_name.upper()
        if is_apo:
            if base == 'rmsd_prot':
                return f"{self.config.protein_name} backbone · Apo"
            return 'Apo'

        if base in ['comcom', 'rmsd_lig']:
            # Ligand-only figures: just the system (ligand) name
            return system_name
        if base == 'rmsd_prot':
            # One backbone curve per complex, so name the ligand to tell them apart.
            return f"{self.config.protein_name} backbone · {system_name}"
        return f"{self.config.protein_name}–{system_name} complex"

    def _metric_base_name(self, metric_name: str) -> str:
        for metric in self.config.metrics:
            if metric.name == metric_name:
                return metric.base_name
        return metric_name
    
    def _update_layout_simple(self, fig: go.Figure, metric: AnalysisMetric, 
                              x_values: List):
        """Update layout for simple style plots"""
        # Determine x-axis title
        if metric.base_name == 'rmsf':
            xaxis_title = "Residue"
        else:
            xaxis_title = self._time_axis_label()

        title_text = metric.title
        if metric.name == 'comcom' and metric.title == 'COM-COM Distance':
            title_text = (
                f"Distance from the center of mass of {self.config.protein_name} "
                "to the center of mass of the ligands"
            )

        fig.update_layout(
            title=dict(
                text=title_text,
                y=0.95,
                x=0.5,
                xanchor='center',
                yanchor='top'
            ),
            xaxis_title=xaxis_title,
            yaxis_title=metric.ylabel,
            legend=dict(
                y=-0.4,  # Matching old code: y=-0.4
                xanchor="center",
                x=0.5,
                orientation='h'
            ),
            font_family=self.config.plot_config.font_family,
            font_color="black",
            font_size=self.config.plot_config.font_size,
            title_font_family=self.config.plot_config.font_family,  # Matching old code
            title_font_color="black",  # Matching old code
            template=self.config.plot_config.template
        )
        
        # Apply metric-specific y-axis limits
        self._apply_axis_limits(fig, metric, list(fig.data))
    
    def _apply_axis_limits(self, fig: go.Figure, metric: AnalysisMetric, traces: List):
        """
        Apply metric-specific y-axis limits.
        From holo_apo_vis.py modify_limits().
        """
        if not traces:
            return
        
        # Extract all y-values
        all_y_values = []
        for trace in traces:
            if hasattr(trace, 'y') and trace.y is not None:
                all_y_values.extend([y for y in trace.y if y is not None and np.isfinite(y)])
        
        if not all_y_values:
            return
        
        metric_name = metric.base_name

        try:
            if metric_name == 'comcom':
                upper_limit = max(all_y_values) + 1
                fig.update_yaxes(range=[0, upper_limit])
            
            elif metric_name == 'hbonds':
                upper_limit = max(all_y_values) + 2
                fig.update_yaxes(range=[0, upper_limit])
            
            elif metric_name == 'rog':
                lower_limit = min(all_y_values) - 0.5
                upper_limit = max(all_y_values) + 0.5
                fig.update_yaxes(range=[lower_limit, upper_limit])
            
            elif metric_name == 'sasa':
                lower_limit = min(all_y_values) - 1000
                upper_limit = max(all_y_values) + 1000
                fig.update_yaxes(range=[lower_limit, upper_limit])
            
            elif metric_name in ['rmsd_prot', 'rmsd_lig', 'rmsd_complex', 'rmsd_combined']:
                upper_limit = max(all_y_values) + 1
                fig.update_yaxes(range=[0, upper_limit])
            
            elif metric_name == 'rmsf':
                upper_limit = max(all_y_values) + 1
                fig.update_yaxes(range=[0, upper_limit])
        
        except (ValueError, TypeError) as e:
            logger.warning(f"Could not set axis limits for {metric_name}: {e}")
    
    def _generate_statistics_text(self, data: Dict) -> str:
        """
        Generate statistics text box for comparative plots.
        
        Args:
            data: Dict of {system_name: processed_data_dict}
            
        Returns:
            Formatted statistics string
        """
        lines = ["<b>Statistics</b>", ""]
        
        for system_name, system_data in data.items():
            stats = system_data.get('stats', {})
            
            lines.append(f"<b>{system_name}:</b>")
            lines.append(f"  Mean: {stats.get('mean', 0):.3f}")
            lines.append(f"  Std:  {stats.get('std', 0):.3f}")
            lines.append(f"  Min:  {stats.get('min', 0):.3f}")
            lines.append(f"  Max:  {stats.get('max', 0):.3f}")
            lines.append("")
        
        return "<br>".join(lines)
    
    # ==================== SAVE PLOTS ====================
    
    def save_plot(self, fig: go.Figure, output_path: Path, base_name: str, fallback: Optional[Dict] = None):
        """
        Save plot in multiple formats.
        
        Args:
            fig: Plotly Figure object
            output_path: Output directory
            base_name: Base filename (without extension)
            fallback: Optional fallback payload for matplotlib SVG export
        """
        output_path.mkdir(parents=True, exist_ok=True)
        
        for fmt in self.config.plot_config.save_formats:
            filepath = output_path / f"{base_name}.{fmt}"
            
            try:
                if fmt == 'html':
                    fig.write_html(filepath)
                elif fmt in ['svg', 'png', 'pdf']:
                    write_static(fig, filepath, fmt, self.config.plot_config.width, self.config.plot_config.height, self.config.plot_config.scale)
                logger.info(f"Saved plot: {filepath}")
            except Exception as e:
                if fmt == 'svg' and fallback:
                    if self._save_svg_fallback(filepath, fallback):
                        logger.info(f"Saved SVG via matplotlib fallback: {filepath}")
                        continue
                logger.error(f"Could not save {filepath}: {e}")
    
    def create_all_plots(self, all_data: Dict) -> List[str]:
        """
        Create all plots for all metrics.
        
        Args:
            all_data: Complete data structure from processor
            
        Returns:
            List of created plot names
        """
        plots_created = []
        output_dir = self.config.output_dir / "plots"
        
        # Get all metrics
        metrics_to_plot = set()
        for system_data in all_data.values():
            metrics_to_plot.update(system_data.keys())
        
        # Create plot for each metric
        for metric_name in metrics_to_plot:
            # Find metric object
            metric = None
            for m in self.config.metrics:
                if m.name == metric_name:
                    metric = m
                    break
            
            if not metric:
                logger.warning(f"Metric object not found for: {metric_name}")
                continue
            
            # Collect data for this metric from all systems
            metric_data = {}
            for system_name, system_data in all_data.items():
                if metric_name in system_data:
                    metric_data[system_name] = system_data[metric_name]
            
            if not metric_data:
                continue
            
            # Prepare x-values for RMSF
            x_values = None
            if metric.name == 'rmsf':
                if self.config.amino_acids:
                    x_values = [str(aa) for aa in self.config.amino_acids]
                else:
                    try:
                        from .data_processor import MDDataProcessor
                        x_values = MDDataProcessor(self.config).prepare_amino_acid_axis()
                    except Exception:
                        x_values = None
            
            # Create plot
            fig = self.create_plot(metric_data, metric, x_values)
            
            # Save plot
            self.save_plot(
                fig,
                output_dir,
                metric_name,
                fallback={
                    'kind': 'metric',
                    'metric': metric,
                    'metric_data': metric_data,
                    'x_values': x_values
                }
            )
            plots_created.append(metric_name)

            if metric.base_name == 'rmsf':
                self._write_terminal_detection(metric_data, x_values)

        # Create overview plot if not in overview style mode
        if self.style != 'overview':
            overview_fig = self.create_overview_plot(all_data)
            if overview_fig:
                self.save_plot(
                    overview_fig,
                    output_dir,
                    'overview',
                    fallback={'kind': 'overview', 'all_data': all_data}
                )
                plots_created.append('overview')

        # Create combined RMSD plot (APO protein vs HOLO complex)
        combined_fig = self.create_combined_rmsd_plot(all_data)
        if combined_fig:
            self.save_plot(
                combined_fig,
                output_dir,
                'rmsd_combined',
                fallback={'kind': 'combined', 'all_data': all_data}
            )
            plots_created.append('rmsd_combined')
        
        return plots_created

    def _write_terminal_detection(self, metric_data: Dict, x_values: Optional[List]) -> Optional[Path]:
        """Record the RMSF terminal region shown on the plot, so it can be cited or reused."""
        terminal = self._rmsf_terminal_region(metric_data)
        if terminal is None:
            return None
        axis = self._residue_axis(metric_data, x_values)
        data_dir = self.config.output_dir / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        path = data_dir / "rmsf_terminal_detection.json"
        payload = dict(terminal)
        payload["core_first_residue"] = axis[terminal["core_start"]] if terminal["core_start"] < len(axis) else None
        payload["core_last_residue"] = axis[terminal["core_end"]] if terminal["core_end"] < len(axis) else None
        payload["systems"] = list(metric_data.keys())
        path.write_text(json.dumps(payload, indent=2) + "\n")
        return path

    def _save_svg_fallback(self, filepath: Path, fallback: Dict) -> bool:
        kind = fallback.get('kind')
        try:
            if kind == 'metric':
                return self._save_metric_svg(
                    filepath,
                    fallback.get('metric'),
                    fallback.get('metric_data', {}),
                    fallback.get('x_values')
                )
            if kind == 'overview':
                return self._save_overview_svg(filepath, fallback.get('all_data', {}))
            if kind == 'combined':
                return self._save_combined_svg(filepath, fallback.get('all_data', {}))
        except Exception as exc:
            logger.error(f"Matplotlib fallback failed for {filepath}: {exc}")
        return False

    def _save_metric_svg(self, filepath: Path, metric: AnalysisMetric,
                         metric_data: Dict, x_values: Optional[List]) -> bool:
        if not metric or not metric_data:
            return False

        _, plt = init_matplotlib()
        fig, ax = plt.subplots(
            figsize=(self.config.plot_config.width / 100, self.config.plot_config.height / 100),
            dpi=100 * self.config.plot_config.scale
        )

        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size

        all_y_values = []
        label_map = None
        terminal = self._rmsf_terminal_region(metric_data) if metric.base_name == 'rmsf' else None
        terminal_x = None

        for system_name, system_data in metric_data.items():
            color = rgb_to_mpl(self._get_system_color(system_name))
            display_name = self._format_legend_name(system_name, metric.name)
            opacity = 0.75 if metric.base_name == 'hbonds' else 1.0

            if system_data['mode'] == 'single':
                y_vals = list(system_data.get('values', []))
                x_plot, label_map = self._get_x_values(metric, system_data, y_vals, x_values)
                if terminal is not None:
                    lo, hi = terminal['core_start'], terminal['core_end']
                    terminal_x = x_plot
                    ax.plot(x_plot, y_vals[:len(x_plot)], color=color, linewidth=2, alpha=0.3)
                    ax.plot(x_plot[lo:hi + 1], y_vals[lo:hi + 1], color=color, linewidth=2, label=display_name)
                else:
                    ax.plot(x_plot, y_vals[:len(x_plot)], color=color, linewidth=2, alpha=opacity, label=display_name)
                all_y_values.extend(y_vals)
            else:
                mean_vals = list(system_data.get('mean', []))
                upper_vals = list(system_data.get('upper', []))
                lower_vals = list(system_data.get('lower', []))
                x_plot, label_map = self._get_x_values(metric, system_data, mean_vals, x_values)
                min_len = min(len(x_plot), len(mean_vals), len(upper_vals), len(lower_vals))
                x_plot = x_plot[:min_len]
                mean_vals = mean_vals[:min_len]
                upper_vals = upper_vals[:min_len]
                lower_vals = lower_vals[:min_len]
                if terminal is not None:
                    lo, hi = terminal['core_start'], terminal['core_end']
                    terminal_x = x_plot
                    ax.plot(x_plot, mean_vals, color=color, linewidth=2, alpha=0.3)
                    ax.plot(x_plot[lo:hi + 1], mean_vals[lo:hi + 1], color=color, linewidth=2, label=display_name)
                else:
                    ax.plot(x_plot, mean_vals, color=color, linewidth=2, alpha=opacity, label=display_name)
                ax.fill_between(x_plot, lower_vals, upper_vals, color=color, alpha=0.18)
                all_y_values.extend(mean_vals)
                all_y_values.extend(upper_vals)
                all_y_values.extend(lower_vals)

        if terminal is not None and terminal_x:
            muted = MUTED_HEX
            if terminal['n_trimmed']:
                ax.axvspan(terminal_x[0], terminal_x[terminal['core_start']], color=muted, alpha=0.08, lw=0)
                ax.axvline(terminal_x[terminal['core_start']], color=muted, ls='--', lw=1.2)
            if terminal['c_trimmed']:
                ax.axvspan(terminal_x[terminal['core_end']], terminal_x[-1], color=muted, alpha=0.08, lw=0)
                ax.axvline(terminal_x[terminal['core_end']], color=muted, ls='--', lw=1.2)
            if terminal['show_threshold']:
                ax.axhline(terminal['threshold'], color=muted, ls=':', lw=1)
            ax.text(
                0.0, 1.02,
                "Faded: flexible termini (smoothed RMSF > median + "
                f"{terminal['mad_multiplier']:g}·MAD, dotted). Dashed: PCA core boundaries.",
                transform=ax.transAxes, fontsize=max(8, self.config.plot_config.font_size - 8),
            )

        ax.set_title(metric.title, pad=28 if terminal is not None else 6)
        ax.set_xlabel("Residue" if metric.base_name == 'rmsf' else self._time_axis_label())
        ax.set_ylabel(metric.ylabel)
        ax.grid(True, alpha=0.2)

        y_limits = self._compute_axis_limits(metric.name, all_y_values)
        if y_limits:
            ax.set_ylim(y_limits)

        if label_map:
            step = max(1, len(label_map) // 30)
            ax.set_xticks(list(range(0, len(label_map), step)))
            ax.set_xticklabels([label_map[i] for i in range(0, len(label_map), step)], rotation=90)

        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False)
        fig.subplots_adjust(bottom=0.25)
        fig.savefig(filepath, format='svg', bbox_inches='tight')
        plt.close(fig)
        return True

    def _save_overview_svg(self, filepath: Path, all_data: Dict) -> bool:
        if not all_data:
            return False

        key_metrics = ['rmsd_prot', 'comcom', 'hbonds', 'rog']
        available_metrics = []
        metric_objects = {}
        first_system = next(iter(all_data.values()))

        for metric_name in key_metrics:
            if metric_name in first_system:
                available_metrics.append(metric_name)
                for m in self.config.metrics:
                    if m.name == metric_name:
                        metric_objects[metric_name] = m
                        break

        if len(available_metrics) < 2:
            return False

        _, plt = init_matplotlib()
        fig, axes = plt.subplots(
            2, 2,
            figsize=(self.config.plot_config.width / 100, self.config.plot_config.height / 100),
            dpi=100 * self.config.plot_config.scale
        )
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size - 2

        handles = []
        labels = []

        for idx, metric_name in enumerate(available_metrics[:4]):
            row = idx // 2
            col = idx % 2
            ax = axes[row][col]
            metric = metric_objects[metric_name]

            all_y_values = []
            for system_name, system_data in all_data.items():
                if metric_name not in system_data:
                    continue
                data = system_data[metric_name]
                color = rgb_to_mpl(self._get_system_color(system_name))
                display_name = self._format_legend_name(system_name, metric_name)
                y_values = data['values'] if data['mode'] == 'single' else data['mean']
                x_values = data.get('time', list(range(len(y_values))))
                ax.plot(x_values, y_values, color=color, linewidth=2, label=display_name)
                all_y_values.extend(y_values)

            ax.set_title(metric.title)
            ax.set_xlabel(self._time_axis_label() if row == 1 else "")
            ax.set_ylabel(metric.ylabel)
            ax.grid(True, alpha=0.2)

            y_limits = self._compute_axis_limits(metric.name, all_y_values)
            if y_limits:
                ax.set_ylim(y_limits)

            if idx == 0:
                handles, labels = ax.get_legend_handles_labels()

        if handles:
            fig.legend(handles, labels, loc='lower center', ncol=2, frameon=False)

        fig.suptitle(f"{self.config.protein_name} - MD Analysis Overview")
        fig.subplots_adjust(bottom=0.2, top=0.9, hspace=0.4, wspace=0.3)
        fig.savefig(filepath, format='svg', bbox_inches='tight')
        plt.close(fig)
        return True

    def _save_combined_svg(self, filepath: Path, all_data: Dict) -> bool:
        if not all_data:
            return False

        apo_system = next((s for s in self.config.systems if s.is_apo), None)
        holo_series = [
            (s.name, all_data.get(s.name, {}).get('rmsd_complex'))
            for s in self.config.systems if not s.is_apo
        ]
        holo_series = [(name, data) for name, data in holo_series if data is not None]
        if not apo_system or not holo_series:
            return False

        apo_data = all_data.get(apo_system.name, {}).get('rmsd_prot')
        if apo_data is None:
            return False

        _, plt = init_matplotlib()
        fig, ax = plt.subplots(
            figsize=(self.config.plot_config.width / 100, self.config.plot_config.height / 100),
            dpi=100 * self.config.plot_config.scale
        )
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size

        def add_series(label: str, data: Dict, color: Tuple[int, int, int]):
            color = rgb_to_mpl(color)
            if data.get('mode') == 'single':
                y_values = list(data.get('values', []))
                x_values = data.get('time', list(range(len(y_values))))
            else:
                y_values = list(data.get('mean', []))
                x_values = data.get('time', list(range(len(y_values))))
                upper = data.get('upper')
                lower = data.get('lower')
                if upper is not None and lower is not None:
                    min_len = min(len(x_values), len(upper), len(lower))
                    ax.fill_between(
                        x_values[:min_len],
                        list(lower)[:min_len],
                        list(upper)[:min_len],
                        color=color,
                        alpha=0.25
                    )

            ax.plot(x_values, y_values, color=color, linewidth=2, label=label)
            return y_values

        all_y = []
        all_y.extend(add_series(f"{self.config.protein_name} · Apo", apo_data, self._get_system_color(apo_system.name)))
        for name, data in holo_series:
            all_y.extend(add_series(f"{self.config.protein_name}–{name} complex", data, self._get_system_color(name)))

        ax.set_title("Protein RMSD: apo vs complexes")
        ax.set_xlabel(self._time_axis_label())
        ax.set_ylabel("RMSD (Å)")
        ax.grid(True, alpha=0.2)

        y_limits = self._compute_axis_limits('rmsd_combined', all_y)
        if y_limits:
            ax.set_ylim(y_limits)

        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False)
        fig.subplots_adjust(bottom=0.25)
        fig.savefig(filepath, format='svg', bbox_inches='tight')
        plt.close(fig)
        return True

    def _get_x_values(self, metric: AnalysisMetric, system_data: Dict,
                      y_vals: List, x_values: Optional[List]):
        if self._has_values(x_values):
            x_vals = list(x_values)
        else:
            x_vals = system_data.get('time')
            if x_vals is None or len(x_vals) == 0:
                x_vals = list(range(len(y_vals)))
        if isinstance(x_vals, np.ndarray):
            x_vals = x_vals.tolist()

        x_vals = list(x_vals)[:len(y_vals)]

        if metric.name == 'rmsf':
            try:
                x_numeric = [int(float(x)) for x in x_vals]
                return x_numeric, None
            except Exception:
                return list(range(len(y_vals))), [str(x) for x in x_vals]

        return x_vals, None

    def _compute_axis_limits(self, metric_name: str, all_y_values: List):
        cleaned = [y for y in all_y_values if y is not None and not np.isnan(y)]
        if not cleaned:
            return None

        try:
            if metric_name == 'comcom':
                return (0, max(cleaned) + 1)
            if metric_name == 'hbonds':
                return (0, max(cleaned) + 2)
            if metric_name == 'rog':
                return (min(cleaned) - 0.5, max(cleaned) + 0.5)
            if metric_name == 'sasa':
                return (min(cleaned) - 1000, max(cleaned) + 1000)
            if metric_name in ['rmsd_prot', 'rmsd_lig', 'rmsd_complex', 'rmsd_combined', 'rmsf']:
                return (0, max(cleaned) + 1)
        except (ValueError, TypeError):
            return None

        return None

    def create_combined_rmsd_plot(self, all_data: Dict) -> Optional[go.Figure]:
        """
        Combined RMSD plot: APO protein backbone vs every ligand-bound complex.

        Uses the standard combined RMSD output name: rmsd_combined.
        """
        apo_system = next((s for s in self.config.systems if s.is_apo), None)
        holo_systems = [s for s in self.config.systems if not s.is_apo]
        if not apo_system or not holo_systems:
            return None

        apo_data = all_data.get(apo_system.name, {}).get('rmsd_prot')
        holo_series = [
            (s.name, all_data.get(s.name, {}).get('rmsd_complex'))
            for s in holo_systems
        ]
        holo_series = [(name, data) for name, data in holo_series if data is not None]
        if apo_data is None or not holo_series:
            return None

        fig = go.Figure()

        def add_series(label: str, system_name: str, data: Dict):
            color = self._get_system_color(system_name)
            if data.get('mode') == 'single':
                y_values = data.get('values', [])
            else:
                y_values = data.get('mean', [])
            x_values = data.get('time', np.arange(len(y_values)))
            upper = data.get('upper') if data.get('mode') != 'single' else None
            lower = data.get('lower') if data.get('mode') != 'single' else None
            if upper is not None and lower is not None:
                x_list = list(x_values)
                fig.add_trace(go.Scatter(
                    x=x_list + x_list[::-1],
                    y=list(upper) + list(lower)[::-1],
                    fill='toself',
                    fillcolor=rgba(color, 0.18),
                    line=dict(width=0),
                    showlegend=False,
                    hoverinfo='skip',
                    legendgroup=label,
                ))
            fig.add_trace(go.Scatter(
                x=x_values,
                y=y_values,
                mode='lines',
                name=label,
                legendgroup=label,
                line=dict(color=rgba(color, 1.0), width=2, dash=self._system_dash.get(system_name, 'solid')),
            ))

        add_series(f"{self.config.protein_name} · Apo", apo_system.name, apo_data)
        for name, data in holo_series:
            add_series(f"{self.config.protein_name}–{name} complex", name, data)

        metric = AnalysisMetric(
            name='rmsd_combined',
            file_pattern='',
            title='Protein RMSD: apo vs complexes',
            ylabel='RMSD (Å)',
            data_format='single_column'
        )
        self._update_layout_simple(fig, metric, [])
        self._layout_legend(fig, 1 + len(holo_series))

        return fig
