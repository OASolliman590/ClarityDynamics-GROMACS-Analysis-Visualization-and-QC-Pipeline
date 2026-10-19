"""
MMPBSA Plotter Module
=====================

Publication-quality visualizations for MMPBSA analysis.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, List, Tuple
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .config import MMPBSAConfig
from ..utils.mpl_fallback import init_matplotlib
from ..utils.figure_export import write_static
from ..utils.palette import categorical

logger = logging.getLogger(__name__)


class MMPBSAPlotter:
    """
    Generate publication-quality MMPBSA visualizations.
    
    Features:
    - Component plots with entropy annotations
    - Decomposition plots with binding affinity display
    - Error bars for replicates
    - Publication styling
    """
    
    def __init__(self, config: MMPBSAConfig):
        """
        Initialize plotter with configuration.
        
        Args:
            config: MMPBSAConfig object
        """
        self.config = config
        # Ligands keep the same palette slot in every MM/PBSA figure (config order).
        self._slots = {system.name: idx for idx, system in enumerate(getattr(config, 'systems', []) or [])}
        logger.info("MMPBSAPlotter initialized")

    def _system_color(self, system_name: str) -> str:
        if system_name not in self._slots:
            self._slots[system_name] = len(self._slots)
        return categorical(self._slots[system_name])
    
    def create_component_plot(self, results_data: Dict, system_name: str) -> go.Figure:
        """
        Create energy components bar chart with entropy annotations.
        
        Args:
            results_data: Processed results data from processor
            system_name: Name of the system
            
        Returns:
            Plotly figure
        """
        delta_components = results_data.get('delta_components', {})
        entropy = results_data.get('entropy', {})
        mode = results_data.get('mode', 'single')
        
        # Extract component values
        components = []
        values = []
        errors = []
        
        # Order: VDWAALS, EEL, EGB, ESURF, GGAS, GSOLV
        component_order = ['ΔVDWAALS', 'ΔEEL', 'ΔEGB', 'ΔESURF', 'ΔGGAS', 'ΔGSOLV']
        
        for comp in component_order:
            if comp in delta_components:
                comp_data = delta_components[comp]
                components.append(comp.replace('Δ', ''))
                values.append(comp_data['mean'])
                errors.append(comp_data['std'] if mode == 'replicate' else 0.0)
        
        # Create figure
        fig = go.Figure()
        
        # Add component bars (no legend entry for individual components)
        fig.add_trace(go.Bar(
            x=components,
            y=values,
            error_y=dict(
                type='data',
                array=errors,
                visible=True if mode == 'replicate' else False,
                thickness=1.5,
                color='rgba(0, 0, 0, 0.4)'
            ),
            marker_color='#2f6690',
            marker_line=dict(color='rgba(0, 0, 0, 0.15)', width=1),
            showlegend=False,  # Hide from legend to reduce clutter
            hovertemplate='<b>%{x}</b><br>Δ Energy: %{y:.2f} kcal/mol<extra></extra>'
        ))
        
        # Add TOTAL bar separately (highlighted)
        if 'ΔTOTAL' in delta_components or 'ΔGTOTAL' in delta_components:
            total_comp = delta_components.get('ΔTOTAL') or delta_components.get('ΔGTOTAL')
            total_value = total_comp['mean']
            total_error = total_comp['std'] if mode == 'replicate' else 0.0
            
            fig.add_trace(go.Bar(
                x=['TOTAL'],
                y=[total_value],
                error_y=dict(
                    type='data',
                    array=[total_error],
                    visible=True if mode == 'replicate' else False,
                    thickness=1.5,
                    color='rgba(0, 0, 0, 0.45)'
                ),
                marker_color='#ef6351',  # Red for TOTAL
                marker_line=dict(color='rgba(0, 0, 0, 0.18)', width=1),
                showlegend=False,  # Hide from legend
                hovertemplate='<b>ΔTOTAL</b><br>Δ Energy: %{y:.2f} kcal/mol<extra></extra>'
            ))
        
        # Add zero line
        fig.add_hline(
            y=0,
            line=dict(color='rgba(0, 0, 0, 0.25)', dash='dot', width=1)
        )
        
        # Add entropy annotations (positioned lower to avoid overlap with title)
        annotations = []
        y_pos = 0.98  # Lower position to avoid title overlap
        
        if 'interaction_entropy' in entropy:
            ie_data = entropy['interaction_entropy']
            ie_value = ie_data['mean']
            ie_std = ie_data['std'] if mode == 'replicate' else 0.0
            
            color = '#2a9d8f' if ie_value >= 0 else '#b9375e'
            text = f"Interaction Entropy: {ie_value:.2f}"
            if mode == 'replicate' and ie_std > 0:
                text += f" ± {ie_std:.2f}"
            text += " kcal/mol"
            
            annotations.append(dict(
                x=0.02,  # Left side to avoid center overlap
                y=y_pos,
                xref='paper',
                yref='paper',
                xanchor='left',
                yanchor='top',
                text=text,
                showarrow=False,
                font=dict(color=color, size=12, family=self.config.plot_config.font_family),
                bgcolor='rgba(255, 255, 255, 0.92)',
                bordercolor='#4f4f4f',
                borderwidth=1,
                borderpad=6
            ))
            y_pos -= 0.06
        
        if 'c2_entropy' in entropy:
            c2_data = entropy['c2_entropy']
            c2_value = c2_data['mean']
            c2_std = c2_data['std'] if mode == 'replicate' else 0.0
            
            color = '#2a9d8f' if c2_value >= 0 else '#b9375e'
            text = f"C2 Entropy: {c2_value:.2f}"
            if mode == 'replicate' and c2_std > 0:
                text += f" ± {c2_std:.2f}"
            text += " kcal/mol"
            
            annotations.append(dict(
                x=0.02,  # Left side to avoid center overlap
                y=y_pos,
                xref='paper',
                yref='paper',
                xanchor='left',
                yanchor='top',
                text=text,
                showarrow=False,
                font=dict(color=color, size=12, family=self.config.plot_config.font_family),
                bgcolor='rgba(255, 255, 255, 0.92)',
                bordercolor='#4f4f4f',
                borderwidth=1,
                borderpad=6
            ))
        
        # Update layout
        # Format title with protein and ligand names (remove "TOTAL" from title)
        title_text = f"Energy Components · {self.config.protein_name}-{self.config.ligand_name}"
        fig.update_layout(
            title=dict(
                text=title_text,
                font=dict(size=self.config.plot_config.font_size)
            ),
            xaxis=dict(
                title='Energy Terms',
                showgrid=False
            ),
            yaxis=dict(
                title='ΔEnergy (kcal/mol)',
                gridcolor='rgba(0, 0, 0, 0.08)',
                zeroline=False
            ),
            font=dict(
                family=self.config.plot_config.font_family,
                size=self.config.plot_config.font_size,
                color='#1f1f1f'
            ),
            showlegend=False,  # Hide legend - TOTAL bar is visually distinct (red color)
            margin=dict(l=70, r=140, t=120, b=60),  # Reduced top margin since legend is lower
            bargap=0.25,
            width=self.config.plot_config.width,
            height=self.config.plot_config.height,
            template=self.config.plot_config.template,
            annotations=annotations
        )
        
        return fig
    
    def create_decomp_plot(self, decomp_data: pd.DataFrame, results_data: Dict, 
                          system_name: str) -> go.Figure:
        """
        Create per-residue decomposition plot with binding affinity annotation.
        
        Args:
            decomp_data: DataFrame with per-residue contributions
            results_data: Processed results data (for binding affinity)
            system_name: Name of the system
            
        Returns:
            Plotly figure
        """
        if decomp_data is None or decomp_data.empty:
            logger.warning("No decomposition data provided")
            return go.Figure()
        
        # Get binding affinity from results
        binding_energy = results_data.get('binding_energy', {})
        mode = results_data.get('mode', 'single')
        
        total_binding = binding_energy.get('TOTAL', {})
        bind_value = total_binding.get('mean', 0.0)
        bind_std = total_binding.get('std', 0.0) if mode == 'replicate' else 0.0
        
        # Prepare data
        if 'TOTAL Avg.' in decomp_data.columns:
            residues = decomp_data['Residue'].values
            energies = decomp_data['TOTAL Avg.'].values
            errors = decomp_data.get('TOTAL Avg. SD', pd.Series([0.0] * len(decomp_data))).values if mode == 'replicate' else np.zeros(len(decomp_data))
        else:
            logger.error("TOTAL Avg. column not found in decomposition data")
            return go.Figure()
        
        # Format residue labels
        residue_labels = [self._format_residue_label(r) for r in residues]
        
        # Create figure
        fig = go.Figure()
        
        # Add bars
        fig.add_trace(go.Bar(
            x=residue_labels,
            y=energies,
            error_y=dict(
                type='data',
                array=errors,
                visible=True if mode == 'replicate' else False,
                thickness=1.5,
                color='rgba(47, 102, 144, 0.6)'
            ),
            marker_color='#2f6690',
            marker_line=dict(color='rgba(0, 0, 0, 0.15)', width=1),
            hovertemplate='<b>%{x}</b><br>Total Energy: %{y:.2f} kcal/mol<extra></extra>'
        ))
        
        # Add zero line
        fig.add_hline(
            y=0,
            line=dict(color='rgba(0, 0, 0, 0.3)', dash='dot', width=1)
        )
        
        # Calculate y-axis range
        if len(energies) > 0:
            min_val = float(np.min(energies - errors if mode == 'replicate' else energies))
            max_val = float(np.max(energies + errors if mode == 'replicate' else energies))
            span = max_val - min_val
            padding = max(span * 0.08, 0.5) if span else max(abs(max_val) * 0.1, 0.5)
            
            y_min = min_val - padding
            y_max = max_val + padding
            
            # Ensure zero is visible
            if y_min > 0:
                y_min = -0.5
            if y_max < 0:
                y_max = 0.5
            
            fig.update_yaxes(range=[y_min, y_max])
        
        # Format binding affinity for title or legend
        bind_text = f"ΔG<sub>bind</sub> = {bind_value:.2f}"
        if mode == 'replicate' and bind_std > 0:
            bind_text += f" ± {bind_std:.2f}"
        bind_text += " kcal/mol"
        
        # Update layout
        # Format title with protein, ligand names, and binding affinity
        title_text = f"Per-Residue Energy Contribution · {self.config.protein_name}-{self.config.ligand_name}<br><span style='font-size: {self.config.plot_config.font_size - 4}px; color: #d32f2f;'>Binding Affinity: {bind_text}</span>"
        fig.update_layout(
            title=dict(
                text=title_text,
                font=dict(size=self.config.plot_config.font_size)
            ),
            xaxis=dict(
                title='Residue',
                tickangle=-75,
                tickfont=dict(size=self.config.plot_config.font_size - 4),  # Smaller font for dense labels
                automargin=True,
                # Use tickmode='array' with selected indices to show every Nth residue
                tickmode='array',
                tickvals=list(range(0, len(residue_labels), max(1, len(residue_labels) // 30))),  # Show ~30 labels max
                ticktext=[residue_labels[i] for i in range(0, len(residue_labels), max(1, len(residue_labels) // 30))]
            ),
            yaxis=dict(
                title='Total Energy (kcal/mol)',
                zeroline=False,
                gridcolor='rgba(0, 0, 0, 0.1)'
            ),
            font=dict(
                family=self.config.plot_config.font_family,
                size=self.config.plot_config.font_size,
                color='#1f1f1f'
            ),
            margin=dict(l=80, r=30, t=100, b=140),  # Increased top margin for title with binding affinity
            bargap=0.3,  # Slightly more spacing between bars
            width=self.config.plot_config.width + 200,  # Wider for better x-axis label spacing
            height=self.config.plot_config.height,
            template=self.config.plot_config.template
        )
        
        return fig

    def create_binding_energy_comparison(self, all_data: Dict, binding_key: str = "TOTAL") -> go.Figure:
        """
        Create cross-ligand comparison plot for binding energy.
        """
        systems = []
        means = []
        stds = []

        for system_name, system_data in all_data.items():
            results_data = system_data.get('results')
            if not results_data:
                continue
            mean_val, std_val = self._extract_binding_value(results_data, binding_key)
            if mean_val is None:
                continue
            systems.append(system_name)
            means.append(mean_val)
            stds.append(std_val)

        if not systems:
            raise ValueError("No binding energy data found for comparison")


        colors = [self._system_color(name) for name in systems]

        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=systems,
            y=means,
            error_y=dict(
                type='data',
                array=stds,
                visible=True if any(s > 0 for s in stds) else False,
                thickness=1.5,
                color='rgba(0, 0, 0, 0.4)'
            ),
            marker_color=colors,
            marker_line=dict(color='rgba(0, 0, 0, 0.15)', width=1),
            hovertemplate='<b>%{x}</b><br>%{y:.2f} kcal/mol<extra></extra>'
        ))

        fig.add_hline(
            y=0,
            line=dict(color='rgba(0, 0, 0, 0.25)', dash='dot', width=1)
        )

        fig.update_layout(
            title=dict(
                text=f"Binding Energy Comparison · {self.config.protein_name}",
                font=dict(size=self.config.plot_config.font_size)
            ),
            xaxis=dict(title='Ligand/System', showgrid=False),
            yaxis=dict(title=f'Δ{binding_key} (kcal/mol)', gridcolor='rgba(0, 0, 0, 0.08)', zeroline=False),
            font=dict(
                family=self.config.plot_config.font_family,
                size=self.config.plot_config.font_size,
                color='#1f1f1f'
            ),
            showlegend=False,
            margin=dict(l=70, r=40, t=100, b=80),
            width=self.config.plot_config.width,
            height=self.config.plot_config.height,
            template=self.config.plot_config.template,
        )

        return fig

    def create_components_comparison(self, all_data: Dict, component_order: Optional[List[str]] = None) -> go.Figure:
        """
        Create cross-ligand comparison plot for energy components.

        Uses per-component subplots so each term has its own y-scale.
        """
        systems = []
        system_results: Dict[str, Dict] = {}
        for system_name, system_data in all_data.items():
            results_data = system_data.get('results')
            if not results_data:
                continue
            systems.append(system_name)
            system_results[system_name] = results_data

        if not systems:
            raise ValueError("No results available for component comparison")

        if component_order is None:
            component_order = ['ΔVDWAALS', 'ΔEEL', 'ΔEGB', 'ΔESURF', 'ΔGGAS', 'ΔGSOLV', 'ΔTOTAL']

        components = []
        for comp in component_order:
            if any(comp in system_results[s].get('delta_components', {}) for s in systems):
                components.append(comp)

        if not components:
            raise ValueError("No matching delta components found for comparison")


        cols = min(3, len(components))
        rows = int(np.ceil(len(components) / cols))
        subplot_titles = [comp.replace('Δ', '') for comp in components]
        fig = make_subplots(
            rows=rows,
            cols=cols,
            subplot_titles=subplot_titles,
            horizontal_spacing=0.08,
            vertical_spacing=0.16
        )

        for comp_idx, comp in enumerate(components):
            row = comp_idx // cols + 1
            col = comp_idx % cols + 1
            for sys_idx, system_name in enumerate(systems):
                results_data = system_results[system_name]
                delta_components = results_data.get('delta_components', {})
                comp_data = delta_components.get(comp, {})
                value = comp_data.get('mean')
                error = comp_data.get('std', 0.0)
                fig.add_trace(
                    go.Bar(
                        x=[system_name],
                        y=[value],
                        name=system_name,
                        marker_color=self._system_color(system_name),
                        error_y=dict(
                            type='data',
                            array=[error],
                            visible=True if error and error > 0 else False,
                            thickness=1.2,
                            color='rgba(0, 0, 0, 0.4)'
                        ),
                        hovertemplate='<b>%{x}</b><br>Δ Energy: %{y:.2f} kcal/mol<extra></extra>',
                        showlegend=True if comp_idx == 0 else False
                    ),
                    row=row,
                    col=col
                )
            fig.add_hline(
                y=0,
                line=dict(color='rgba(0, 0, 0, 0.25)', dash='dot', width=1),
                row=row,
                col=col
            )
            fig.update_yaxes(
                title_text='ΔEnergy (kcal/mol)',
                gridcolor='rgba(0, 0, 0, 0.08)',
                zeroline=False,
                row=row,
                col=col
            )

        min_row_height = 280
        dynamic_height = max(self.config.plot_config.height, rows * min_row_height + 140)
        fig.update_layout(
            title=dict(
                text=f"Energy Components Comparison · {self.config.protein_name}",
                font=dict(size=self.config.plot_config.font_size)
            ),
            font=dict(
                family=self.config.plot_config.font_family,
                size=self.config.plot_config.font_size,
                color='#1f1f1f'
            ),
            barmode='group',
            margin=dict(l=70, r=40, t=120, b=110),
            width=self.config.plot_config.width,
            height=dynamic_height,
            template=self.config.plot_config.template,
            legend=dict(orientation="h", yanchor="bottom", y=-0.2, xanchor="center", x=0.5)
        )
        fig.update_xaxes(tickangle=30, automargin=True)
        fig.update_annotations(font_size=max(10, self.config.plot_config.font_size - 2))

        return fig

    # ==================== ENTROPY ====================

    ENTROPY_LABELS = {
        'interaction_entropy': 'IE',
        'c2_entropy': 'C2',
    }

    def entropy_terms(self, results_data: Dict) -> List[Dict]:
        """ΔTOTAL, −TΔS and ΔG = ΔTOTAL + (−TΔS) for every entropy method present.

        Uncertainties combine in quadrature (independent estimates).
        """
        total_mean, total_err = self._binding_with_error(results_data)
        if total_mean is None:
            return []
        rows = []
        for key, label in self.ENTROPY_LABELS.items():
            ent = (results_data.get('entropy') or {}).get(key)
            if not ent or ent.get('mean') is None:
                continue
            ent_err = self._error_of(ent)
            rows.append({
                'method': label,
                'key': key,
                'total': total_mean,
                'total_err': total_err,
                'minus_tds': float(ent['mean']),
                'minus_tds_err': ent_err,
                'dg': total_mean + float(ent['mean']),
                'dg_err': float(np.hypot(total_err, ent_err)),
                'sigma': ent.get('sigma'),
                'reliable': ent.get('reliable'),
            })
        return rows

    @staticmethod
    def _error_of(entry: Dict) -> float:
        """SEM over replicates when available, otherwise the reported SD."""
        if entry.get('uncertainty') == 'between_replicates' and entry.get('sem') is not None:
            return float(entry['sem'])
        return float(entry.get('std') or 0.0)

    def _binding_with_error(self, results_data: Dict) -> Tuple[Optional[float], float]:
        binding = results_data.get('binding_energy', {}) or {}
        entry = binding.get('TOTAL')
        if entry is None:
            for key in ("ΔTOTAL", "ΔGTOTAL"):
                entry = (results_data.get('delta_components') or {}).get(key)
                if entry:
                    break
        if not entry:
            return None, 0.0
        return entry.get('mean'), self._error_of(entry)

    def create_entropy_plot(self, results_data: Dict, system_name: str) -> Optional[go.Figure]:
        """Bars for ΔTOTAL, −TΔS and the entropy-corrected ΔG per available method.

        Entropy terms flagged unreliable (IE with σ(ΔE_int) > 3.6 kcal/mol) are hatched.
        """
        rows = self.entropy_terms(results_data)
        if not rows:
            return None

        total_color, entropy_color, dg_color = categorical(0), categorical(1), categorical(2)
        x = ['ΔTOTAL']
        y = [rows[0]['total']]
        err = [rows[0]['total_err']]
        colors = [total_color]
        patterns = ['']
        hover = ['Effective binding energy (no entropy)']
        for row in rows:
            unreliable = row['reliable'] is False
            sigma_txt = f"σ(ΔE_int) = {row['sigma']:.2f} kcal/mol" if row['sigma'] is not None else ''
            x += [f"−TΔS ({row['method']})", f"ΔG ({row['method']})"]
            y += [row['minus_tds'], row['dg']]
            err += [row['minus_tds_err'], row['dg_err']]
            colors += [entropy_color, dg_color]
            patterns += ['/' if unreliable else '', '/' if unreliable else '']
            note = ' · unreliable (σ > 3.6)' if unreliable else ''
            hover += [f"Entropy penalty {sigma_txt}{note}", f"ΔTOTAL + (−TΔS) {note}"]

        fig = go.Figure(go.Bar(
            x=x, y=y,
            error_y=dict(type='data', array=err, visible=any(e > 0 for e in err), thickness=1.5,
                         color='rgba(0,0,0,0.45)'),
            marker=dict(color=colors, pattern=dict(shape=patterns, fgcolor='white', size=8),
                        line=dict(color='white', width=2)),
            text=[f"{v:.2f}" for v in y], textposition='outside',
            customdata=hover,
            hovertemplate='<b>%{x}</b><br>%{y:.2f} kcal/mol<br>%{customdata}<extra></extra>',
        ))
        fig.add_hline(y=0, line=dict(color='rgba(0,0,0,0.3)', width=1))
        notes = [
            f"{r['method']}: σ(ΔE_int) = {r['sigma']:.2f} kcal/mol"
            + (" — above 3.6, entropy estimate unreliable (hatched)" if r['reliable'] is False else "")
            for r in rows if r['sigma'] is not None
        ]
        if notes:
            fig.add_annotation(text="<br>".join(notes), xref='paper', yref='paper', x=0, y=1.02,
                               xanchor='left', yanchor='bottom', showarrow=False, align='left',
                               font=dict(size=max(10, self.config.plot_config.font_size - 6)))
        fig.update_layout(
            title=dict(text=f"Entropy-corrected binding free energy · {system_name}", y=0.97),
            yaxis=dict(title='kcal/mol', zeroline=False, gridcolor='rgba(0,0,0,0.08)'),
            xaxis=dict(showgrid=False),
            showlegend=False,
            font=dict(family=self.config.plot_config.font_family, size=self.config.plot_config.font_size),
            template=self.config.plot_config.template,
            margin=dict(t=140),
            width=max(self.config.plot_config.width, 220 * len(x)),
            height=self.config.plot_config.height,
        )
        return fig

    def create_entropy_comparison(self, all_data: Dict) -> Optional[go.Figure]:
        """ΔTOTAL vs entropy-corrected ΔG for every system, grouped by system."""
        per_system = {}
        for system_name, system_data in all_data.items():
            results_data = system_data.get('results')
            if results_data:
                rows = self.entropy_terms(results_data)
                if rows:
                    per_system[system_name] = rows
        if not per_system:
            return None

        systems = list(per_system)
        methods = sorted({r['method'] for rows in per_system.values() for r in rows})
        series = [('ΔTOTAL', None)] + [(f"ΔG (ΔTOTAL + −TΔS {m})", m) for m in methods]
        fig = go.Figure()
        for idx, (label, method) in enumerate(series):
            ys, errs, patterns = [], [], []
            for system in systems:
                rows = per_system[system]
                if method is None:
                    ys.append(rows[0]['total']); errs.append(rows[0]['total_err']); patterns.append('')
                    continue
                row = next((r for r in rows if r['method'] == method), None)
                ys.append(row['dg'] if row else None)
                errs.append(row['dg_err'] if row else 0.0)
                patterns.append('/' if row and row['reliable'] is False else '')
            fig.add_trace(go.Bar(
                name=label, x=systems, y=ys,
                error_y=dict(type='data', array=errs, visible=any(e for e in errs), thickness=1.5,
                             color='rgba(0,0,0,0.45)'),
                marker=dict(color=categorical(idx), pattern=dict(shape=patterns, fgcolor='white', size=8),
                            line=dict(color='white', width=2)),
                hovertemplate='<b>%{x}</b><br>' + label + ': %{y:.2f} kcal/mol<extra></extra>',
            ))
        fig.add_hline(y=0, line=dict(color='rgba(0,0,0,0.3)', width=1))
        fig.update_layout(
            barmode='group', bargap=0.25,
            title=dict(text=f"Binding energy with and without entropy · {self.config.protein_name}"),
            yaxis=dict(title='kcal/mol', zeroline=False, gridcolor='rgba(0,0,0,0.08)'),
            xaxis=dict(title='Ligand/System', showgrid=False),
            legend=dict(orientation='h', x=0.5, xanchor='center', y=-0.2, yanchor='top'),
            font=dict(family=self.config.plot_config.font_family, size=self.config.plot_config.font_size),
            template=self.config.plot_config.template,
            width=max(self.config.plot_config.width, 160 * len(systems) * len(series) // 2 + 300),
            height=self.config.plot_config.height,
        )
        fig.add_annotation(
            text="Hatched: interaction entropy with σ(ΔE_int) > 3.6 kcal/mol (unreliable).",
            xref='paper', yref='paper', x=0, y=1.02, xanchor='left', yanchor='bottom', showarrow=False,
            font=dict(size=max(10, self.config.plot_config.font_size - 6)),
        )
        return fig

    def _extract_binding_value(self, results_data: Dict, binding_key: str) -> Tuple[Optional[float], float]:
        binding = results_data.get('binding_energy', {})
        if binding_key in binding:
            entry = binding[binding_key]
            return entry.get('mean'), entry.get('std', 0.0)

        if binding_key == "TOTAL":
            delta_components = results_data.get('delta_components', {})
            for key in ("ΔTOTAL", "ΔGTOTAL"):
                if key in delta_components:
                    entry = delta_components[key]
                    return entry.get('mean'), entry.get('std', 0.0)

        return None, 0.0
    
    def _format_residue_label(self, residue: str) -> str:
        """
        Format residue identifier for cleaner axis labels.
        
        Args:
            residue: Residue identifier (e.g., "R:A:TRP:106")
            
        Returns:
            Formatted label (e.g., "TRP106")
        """
        if not residue:
            return residue
        
        parts = [p for p in residue.split(':') if p]
        if not parts:
            return residue
        
        # Skip ligand residues (L:)
        if parts[0] == 'L':
            return f"LIG{parts[-1]}" if len(parts) > 1 else residue
        
        # Extract residue name and number
        if len(parts) >= 4:
            return f"{parts[2]}{parts[3]}"
        elif len(parts) >= 3:
            return f"{parts[1]}{parts[2]}"
        elif len(parts) >= 2:
            return ''.join(parts[-2:])
        
        return residue.replace('UNK', '')
    
    def save_plot(self, fig: go.Figure, output_path: Path, name: str,
                  fallback: Optional[Dict] = None):
        """
        Save plot in multiple formats.
        
        Args:
            fig: Plotly figure
            output_path: Output directory
            name: Base name for files
            fallback: Optional fallback payload for matplotlib SVG export
        """
        output_path.mkdir(parents=True, exist_ok=True)
        
        for fmt in self.config.plot_config.save_formats:
            file_path = output_path / f"{name}.{fmt}"
            
            try:
                if fmt == 'html':
                    fig.write_html(str(file_path))
                elif fmt in ('svg', 'png', 'pdf'):
                    write_static(fig, file_path, fmt, self.config.plot_config.width, self.config.plot_config.height, self.config.plot_config.scale)
                
                logger.info(f"Saved {fmt} plot: {file_path}")
            except Exception as e:
                if fmt == 'svg' and fallback:
                    if self._save_svg_fallback(file_path, fallback):
                        logger.info(f"Saved SVG via matplotlib fallback: {file_path}")
                        continue
                logger.error(f"Error saving {fmt} plot {file_path}: {e}")

    def _save_svg_fallback(self, file_path: Path, fallback: Dict) -> bool:
        kind = fallback.get('kind')
        try:
            if kind == 'components':
                return self._save_components_svg(
                    file_path,
                    fallback.get('results_data', {}),
                    fallback.get('system_name', '')
                )
            if kind == 'decomp':
                return self._save_decomp_svg(
                    file_path,
                    fallback.get('decomp_data'),
                    fallback.get('results_data', {}),
                    fallback.get('system_name', '')
                )
            if kind == 'binding_compare':
                return self._save_binding_compare_svg(
                    file_path,
                    fallback.get('systems', []),
                    fallback.get('means', []),
                    fallback.get('stds', [])
                )
            if kind == 'components_compare':
                return self._save_components_compare_svg(
                    file_path,
                    fallback.get('components', []),
                    fallback.get('series', {})
                )
        except Exception as exc:
            logger.error(f"Matplotlib fallback failed for {file_path}: {exc}")
        return False

    def _save_binding_compare_svg(self, file_path: Path, systems: List[str],
                                  means: List[float], stds: List[float]) -> bool:
        if not systems:
            return False

        _, plt = init_matplotlib()
        fig, ax = plt.subplots(
            figsize=(self.config.plot_config.width / 100, self.config.plot_config.height / 100),
            dpi=100 * self.config.plot_config.scale
        )
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size

        x = np.arange(len(systems))
        ax.bar(x, means, yerr=stds if any(s > 0 for s in stds) else None, capsize=4)
        ax.axhline(0, color='gray', linewidth=1, linestyle='--', alpha=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(systems, rotation=45, ha='right')
        ax.set_ylabel('ΔG (kcal/mol)')
        ax.set_title(f"Binding Energy Comparison · {self.config.protein_name}")
        fig.tight_layout()
        fig.savefig(file_path, format="svg", bbox_inches="tight")
        plt.close(fig)
        return True

    def _save_components_compare_svg(self, file_path: Path, components: List[str],
                                     series: Dict[str, Dict[str, List[float]]]) -> bool:
        if not components or not series:
            return False

        _, plt = init_matplotlib()
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size

        systems = list(series.keys())
        cols = min(3, len(components))
        rows = int(np.ceil(len(components) / cols))
        fig, axes = plt.subplots(
            rows,
            cols,
            figsize=(self.config.plot_config.width / 100, self.config.plot_config.height / 100),
            dpi=100 * self.config.plot_config.scale
        )
        axes = np.atleast_1d(axes).flatten()

        for idx, comp in enumerate(components):
            ax = axes[idx]
            values = []
            stds = []
            for system_name in systems:
                means = series[system_name]['means']
                errs = series[system_name]['stds']
                comp_idx = components.index(comp)
                values.append(means[comp_idx])
                stds.append(errs[comp_idx])

            x = np.arange(len(systems))
            ax.bar(x, values, yerr=stds if any(s > 0 for s in stds) else None, capsize=3)
            ax.axhline(0, color='gray', linewidth=1, linestyle='--', alpha=0.5)
            ax.set_title(comp.replace('Δ', ''))
            ax.set_xticks(x)
            ax.set_xticklabels(systems, rotation=45, ha='right')
            ax.set_ylabel('ΔEnergy (kcal/mol)')

        for j in range(len(components), len(axes)):
            fig.delaxes(axes[j])

        fig.suptitle(f"Energy Components Comparison · {self.config.protein_name}")
        fig.tight_layout(rect=[0, 0.05, 1, 0.95])
        fig.savefig(file_path, format="svg", bbox_inches="tight")
        plt.close(fig)
        return True

    def _save_components_svg(self, file_path: Path, results_data: Dict, system_name: str) -> bool:
        if not results_data:
            return False

        _, plt = init_matplotlib()
        fig, ax = plt.subplots(
            figsize=(self.config.plot_config.width / 100, self.config.plot_config.height / 100),
            dpi=100 * self.config.plot_config.scale
        )
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size

        delta_components = results_data.get('delta_components', {})
        entropy = results_data.get('entropy', {})
        mode = results_data.get('mode', 'single')

        component_order = ['ΔVDWAALS', 'ΔEEL', 'ΔEGB', 'ΔESURF', 'ΔGGAS', 'ΔGSOLV']
        components = []
        values = []
        errors = []

        for comp in component_order:
            if comp in delta_components:
                comp_data = delta_components[comp]
                components.append(comp.replace('Δ', ''))
                values.append(comp_data['mean'])
                errors.append(comp_data['std'] if mode == 'replicate' else 0.0)

        x_positions = np.arange(len(components))
        ax.bar(
            x_positions,
            values,
            yerr=errors if mode == 'replicate' else None,
            color='#2f6690',
            edgecolor='black',
            linewidth=0.5
        )

        if 'ΔTOTAL' in delta_components or 'ΔGTOTAL' in delta_components:
            total_comp = delta_components.get('ΔTOTAL') or delta_components.get('ΔGTOTAL')
            total_value = total_comp['mean']
            total_error = total_comp['std'] if mode == 'replicate' else 0.0
            ax.bar(
                [len(components)],
                [total_value],
                yerr=[total_error] if mode == 'replicate' else None,
                color='#ef6351',
                edgecolor='black',
                linewidth=0.5
            )
            components.append('TOTAL')

        ax.axhline(0, color='black', linewidth=0.8, linestyle='--', alpha=0.4)
        ax.set_xticks(np.arange(len(components)))
        ax.set_xticklabels(components)

        ax.set_xlabel('Energy Terms')
        ax.set_ylabel('Delta Energy (kcal/mol)')

        title_text = f"Energy Components - {self.config.protein_name}-{self.config.ligand_name}"
        ax.set_title(title_text)

        annotation_lines = []
        if 'interaction_entropy' in entropy:
            ie_data = entropy['interaction_entropy']
            ie_value = ie_data['mean']
            ie_std = ie_data['std'] if mode == 'replicate' else 0.0
            line = f"Interaction Entropy: {ie_value:.2f}"
            if mode == 'replicate' and ie_std > 0:
                line += f" ± {ie_std:.2f}"
            annotation_lines.append(line)

        if 'c2_entropy' in entropy:
            c2_data = entropy['c2_entropy']
            c2_value = c2_data['mean']
            c2_std = c2_data['std'] if mode == 'replicate' else 0.0
            line = f"C2 Entropy: {c2_value:.2f}"
            if mode == 'replicate' and c2_std > 0:
                line += f" ± {c2_std:.2f}"
            annotation_lines.append(line)

        if annotation_lines:
            fig.text(0.02, 0.95, "\n".join(annotation_lines), ha='left', va='top', fontsize=10)

        fig.tight_layout()
        fig.savefig(file_path, format='svg', bbox_inches='tight')
        plt.close(fig)
        return True

    def _save_decomp_svg(self, file_path: Path, decomp_data: Optional[pd.DataFrame],
                         results_data: Dict, system_name: str) -> bool:
        if decomp_data is None or decomp_data.empty:
            return False

        _, plt = init_matplotlib()
        fig, ax = plt.subplots(
            figsize=((self.config.plot_config.width + 200) / 100, self.config.plot_config.height / 100),
            dpi=100 * self.config.plot_config.scale
        )
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size

        mode = results_data.get('mode', 'single')
        binding_energy = results_data.get('binding_energy', {})
        total_binding = binding_energy.get('TOTAL', {})
        bind_value = total_binding.get('mean', 0.0)
        bind_std = total_binding.get('std', 0.0) if mode == 'replicate' else 0.0

        if 'TOTAL Avg.' not in decomp_data.columns:
            return False

        residues = decomp_data['Residue'].values
        energies = decomp_data['TOTAL Avg.'].values
        errors = (decomp_data.get('TOTAL Avg. SD', pd.Series([0.0] * len(decomp_data))).values
                  if mode == 'replicate' else np.zeros(len(decomp_data)))

        labels = [self._format_residue_label(r) for r in residues]
        x_positions = np.arange(len(labels))

        ax.bar(
            x_positions,
            energies,
            yerr=errors if mode == 'replicate' else None,
            color='#2f6690',
            edgecolor='black',
            linewidth=0.3
        )
        ax.axhline(0, color='black', linewidth=0.8, linestyle='--', alpha=0.4)

        ax.set_xlabel('Residue')
        ax.set_ylabel('Total Energy (kcal/mol)')

        bind_text = f"DeltaGbind = {bind_value:.2f}"
        if mode == 'replicate' and bind_std > 0:
            bind_text += f" ± {bind_std:.2f}"
        bind_text += " kcal/mol"

        title_text = f"Per-Residue Energy Contribution - {self.config.protein_name}-{self.config.ligand_name}\nBinding Affinity: {bind_text}"
        ax.set_title(title_text)

        step = max(1, len(labels) // 30)
        ax.set_xticks(list(range(0, len(labels), step)))
        ax.set_xticklabels([labels[i] for i in range(0, len(labels), step)], rotation=75)

        fig.tight_layout()
        fig.savefig(file_path, format='svg', bbox_inches='tight')
        plt.close(fig)
        return True
