"""
PCA Plotter Module
==================

Create all PCA visualizations.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import plotly.io as pio

from .config import PCAConfig
from ..utils.mpl_fallback import init_matplotlib
from ..utils.figure_export import StaticExportUnavailable, write_static
from ..utils.palette import APO_HEX, MUTED_HEX, categorical, mix, rgba, single_hue_scale

# Set default Plotly export settings.
pio.defaults.default_format = "svg"
pio.defaults.default_width = 1200
pio.defaults.default_height = 600
pio.defaults.default_scale = 2

logger = logging.getLogger(__name__)


class PCAPlotter:
    """
    Generate all PCA visualizations.
    
    Features:
    - Scree plot
    - Cosine content
    - PC histograms
    - PC time series
    - 2D projections (3 combinations)
    - 3D projection
    - RMSIP heatmap
    """
    
    def __init__(self, config: PCAConfig):
        """
        Initialize plotter with configuration.
        
        Args:
            config: PCAConfig object
        """
        self.config = config
        self.holo_color = getattr(config, 'color', None) or categorical(0)
        logger.info("PCAPlotter initialized")
    
    def calculate_variance(self, eigenvals: List[float]) -> List[float]:
        """Calculate cumulative variance from eigenvalues."""
        all_vals = sum(eigenvals)
        variances = []
        for idx, val in enumerate(eigenvals):
            variances.append(sum(eigenvals[0:idx+1]) / all_vals * 100)
        return variances
    
    def update_fig(self, fig: go.Figure, titlename: str, xname: Optional[str], 
                   yname: Optional[str], y: float = 0.9) -> go.Figure:
        """Update figure layout with consistent styling."""
        fig.update_layout(
            title={
                'text': titlename,
                'y': y,
                'x': 0.5,
                'xanchor': 'center',
                'yanchor': 'top'
            },
            xaxis_title=xname,
            yaxis_title=yname,
            font_family=self.config.plot_config.font_family,
            font_color="black",
            font_size=self.config.plot_config.font_size,
            title_font_family=self.config.plot_config.font_family,
            title_font_color="black",
        )
        return fig
    
    def create_scree_plot(self, eigenvals_data: Dict) -> go.Figure:
        """
        Scree plot: eigenvalues (top) and cumulative explained variance (bottom).

        Two stacked panels sharing the eigenvector axis instead of one plot with
        two y-scales, so neither quantity is read against the wrong axis.
        """
        n = min(self.config.n_pcs, len(eigenvals_data['y']))
        xs = list(eigenvals_data['x'][:n])
        ys = list(eigenvals_data['y'][:n])
        variances = self.calculate_variance(eigenvals_data['y'])[:n]

        fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
                            row_heights=[0.55, 0.45])
        fig.add_trace(
            go.Scatter(x=xs, y=ys, mode='lines+markers', name='Eigenvalue',
                       line=dict(color=self.holo_color, width=2), marker=dict(size=9)),
            row=1, col=1,
        )
        fig.add_trace(
            go.Scatter(x=xs, y=variances, mode='lines+markers', name='Cumulative variance',
                       line=dict(color=self.holo_color, width=2, dash='dot'), marker=dict(size=9),
                       text=[f"{v:.1f}%" for v in variances], textposition='top center'),
            row=2, col=1,
        )
        self.update_fig(fig, 'Eigenvalue spectrum and cumulative variance', None, None, y=0.97)
        fig.update_xaxes(dtick=1, range=(0.5, n + 0.5))
        fig.update_xaxes(title_text='Eigenvector', row=2, col=1)
        fig.update_yaxes(title_text='Eigenvalue (nm²)', rangemode='tozero', row=1, col=1)
        fig.update_yaxes(title_text='Cumulative (%)', range=(0, 100), row=2, col=1)
        fig.update_layout(showlegend=False, height=max(self.config.plot_config.height, 750))
        return fig

    def create_cosine_content_plot(self, cc_apo: Dict, cc_holo: Dict,
                                   apo_name: str, holo_name: str) -> go.Figure:
        """
        Cosine content per PC (Hess 2002). Values approaching 1 indicate random-diffusion-like,
        unconverged sampling along that PC, so a reference line is drawn at 0.5.
        """
        fig = go.Figure(data=[
            go.Bar(x=cc_holo['x'], y=cc_holo['y'], name=holo_name, marker_color=self.holo_color),
            go.Bar(x=cc_apo['x'], y=cc_apo['y'], name=apo_name, marker_color=APO_HEX),
        ])
        fig.update_layout(barmode='group', bargap=0.25, bargroupgap=0.08)
        fig.add_hline(y=0.5, line=dict(color=MUTED_HEX, width=1, dash='dot'),
                      annotation_text='0.5: sampling resembles random diffusion above this',
                      annotation_position='top left', annotation_font_size=12)
        fig.update_xaxes(dtick=1, range=(0.5, self.config.n_pcs + 0.5))
        fig.update_yaxes(range=(0, 1))

        title = cc_holo.get('title') or "Cosine content"
        ylabel = cc_holo.get('y_title') or "Cosine content"
        self.update_fig(fig, title, 'Eigenvector', ylabel)
        fig.update_layout(legend=dict(orientation='h', x=0.5, xanchor='center', y=-0.2, yanchor='top'))
        return fig

    def create_pc_histograms(self, proj_data: Dict) -> go.Figure:
        """Distribution of the projection on each of the first N PCs (auto-binned per PC)."""
        n = min(self.config.n_pcs, len(proj_data['ys']))
        ys = proj_data['ys'][:n]
        cols = 5
        rows = (n + cols - 1) // cols
        subplots = make_subplots(rows=rows, cols=cols, subplot_titles=[f'PC {i + 1}' for i in range(n)])

        for pc_idx, y in enumerate(ys):
            row, col = pc_idx // cols + 1, pc_idx % cols + 1
            values = np.asarray(y, dtype=float)
            subplots.add_trace(
                go.Histogram(x=values, nbinsx=40, histnorm='probability density',
                             marker=dict(color=rgba(self.holo_color, 0.55), line=dict(color='white', width=0.5)),
                             name=f'PC {pc_idx + 1}'),
                row=row, col=col,
            )
            kde = self._kde(values)
            if kde is not None:
                subplots.add_trace(
                    go.Scatter(x=kde[0], y=kde[1], mode='lines',
                               line=dict(color=self.holo_color, width=2), name=f'PC {pc_idx + 1} KDE'),
                    row=row, col=col,
                )
            subplots.update_xaxes(title_text='nm', row=row, col=col)

        subplots.update_layout(showlegend=False, height=max(self.config.plot_config.height, 360 * rows + 120))
        self.update_fig(subplots, f'Distribution along the first {n} PCs', None, None, y=0.98)
        return subplots

    @staticmethod
    def _kde(values: np.ndarray, points: int = 200):
        values = values[np.isfinite(values)]
        if len(values) < 5 or np.std(values) == 0:
            return None
        try:
            from scipy.stats import gaussian_kde
        except ImportError:
            return None
        grid = np.linspace(values.min(), values.max(), points)
        return grid, gaussian_kde(values)(grid)

    def create_pc_timeseries(self, proj_data: Dict) -> go.Figure:
        """Projection on each of the first N PCs versus simulation time (x from proj.xvg, ps -> ns)."""
        n = min(self.config.n_pcs, len(proj_data['ys']))
        xs = proj_data['xs'][:n]
        ys = proj_data['ys'][:n]
        x_titles = proj_data.get('x_titles') or []
        x_unit = str(x_titles[0]).lower() if x_titles and x_titles[0] else 'time (ps)'
        to_ns = 1e-3 if 'ps' in x_unit else 1.0

        cols = 2
        rows = (n + cols - 1) // cols
        subplots = make_subplots(rows=rows, cols=cols, subplot_titles=[f'PC {i + 1}' for i in range(n)],
                                 x_title='Time (ns)', vertical_spacing=min(0.1, 0.5 / max(rows, 1)))
        for pc_idx, (x, y) in enumerate(zip(xs, ys)):
            subplots.add_trace(
                go.Scatter(x=np.asarray(x, dtype=float) * to_ns, y=y, mode='lines',
                           line=dict(color=self.holo_color, width=1.2)),
                row=pc_idx // cols + 1, col=pc_idx % cols + 1,
            )
        subplots.update_yaxes(title_text='nm')
        subplots.update_layout(showlegend=False, height=max(self.config.plot_config.height, 200 * rows + 160))
        self.update_fig(subplots, f'Projection of the trajectory on the first {n} PCs', None, None, y=0.98)
        return subplots

    def _time_colored_markers(self, n: int, color: str, title: str, bar_x: float, size: int = 6) -> Dict:
        return dict(
            size=size,
            color=list(range(n)),
            colorscale=single_hue_scale(color),
            opacity=0.75,
            line=dict(width=0),
            showscale=True,
            colorbar=dict(title=dict(text=title, side='right'), x=bar_x, len=0.6, thickness=12,
                          tickvals=[0, max(n - 1, 0)], ticktext=['start', 'end']),
        )

    def _mean_marker(self, color: str, name: str, x: float, y: float, z: Optional[float] = None):
        marker = dict(size=16, symbol='diamond', color=mix(color, '#000000', 0.25),
                      line=dict(color='white', width=2))
        if z is None:
            return go.Scatter(x=[x], y=[y], mode='markers', name=name, marker=marker)
        marker['size'] = 8
        return go.Scatter3d(x=[x], y=[y], z=[z], mode='markers', name=name, marker=marker)

    def create_2d_projection(self, proj_apo: Dict, proj_holo: Dict,
                             pc_pair: Tuple[int, int], apo_name: str, holo_name: str) -> go.Figure:
        """
        2D projection, one hue per system (APO grey, ligand in its own colour).

        Within a system, light-to-dark encodes trajectory progress so early and late
        sampling can be told apart; a diamond marks each system's mean.
        """
        x1, y1 = proj_apo['x'], proj_apo['y']
        x2, y2 = proj_holo['x'], proj_holo['y']
        pc1, pc2 = pc_pair

        fig = go.Figure(data=[
            go.Scatter(x=x1, y=y1, mode='markers', name=apo_name,
                       text=[f'Frame {i}' for i in range(len(x1))],
                       marker=self._time_colored_markers(len(x1), APO_HEX, apo_name, 1.02)),
            go.Scatter(x=x2, y=y2, mode='markers', name=holo_name,
                       text=[f'Frame {i}' for i in range(len(x2))],
                       marker=self._time_colored_markers(len(x2), self.holo_color, holo_name, 1.14)),
            self._mean_marker(APO_HEX, f'{apo_name} mean', float(np.mean(x1)), float(np.mean(y1))),
            self._mean_marker(self.holo_color, f'{holo_name} mean', float(np.mean(x2)), float(np.mean(y2))),
        ])
        self.update_fig(fig, f'Projection on PC {pc1} and PC {pc2}', f'PC{pc1} (nm)', f'PC{pc2} (nm)', y=0.95)
        fig.update_layout(legend=dict(orientation='h', x=0.5, xanchor='center', y=-0.18, yanchor='top'),
                          margin=dict(r=180))
        return fig

    def create_3d_projection(self, proj_apo: Tuple[List[float], List[float], List[float]],
                             proj_holo: Tuple[List[float], List[float], List[float]],
                             apo_name: str, holo_name: str) -> go.Figure:
        """3D projection on PC1-3 with the same colour encoding as the 2D projections."""
        x1, y1, z1 = proj_apo
        x2, y2, z2 = proj_holo

        fig = go.Figure(data=[
            go.Scatter3d(x=x1, y=y1, z=z1, mode='markers', name=apo_name,
                         text=[f'Frame {i}' for i in range(len(x1))],
                         marker=self._time_colored_markers(len(x1), APO_HEX, apo_name, 1.0, size=3)),
            go.Scatter3d(x=x2, y=y2, z=z2, mode='markers', name=holo_name,
                         text=[f'Frame {i}' for i in range(len(x2))],
                         marker=self._time_colored_markers(len(x2), self.holo_color, holo_name, 1.1, size=3)),
            self._mean_marker(APO_HEX, f'{apo_name} mean',
                              float(np.mean(x1)), float(np.mean(y1)), float(np.mean(z1))),
            self._mean_marker(self.holo_color, f'{holo_name} mean',
                              float(np.mean(x2)), float(np.mean(y2)), float(np.mean(z2))),
        ])
        fig.update_layout(
            scene=dict(xaxis_title='PC1 (nm)', yaxis_title='PC2 (nm)', zaxis_title='PC3 (nm)'),
            legend=dict(orientation='h', x=0.5, xanchor='center', y=-0.05, yanchor='top'),
            font_family=self.config.plot_config.font_family,
            font_size=self.config.plot_config.font_size,
            height=max(self.config.plot_config.height, 800),
        )
        return fig

    def create_rmsip_heatmap(self, xpm_data: np.ndarray) -> go.Figure:
        """
        Create RMSIP heatmap from XPM data.
        
        Args:
            xpm_data: NumPy array with matrix values
            
        Returns:
            Plotly figure
        """
        fig = go.Figure(data=[go.Heatmap(
            z=xpm_data,
            colorscale='Greys',
            x=list(range(1, xpm_data.shape[1] + 1)),
            y=list(range(1, xpm_data.shape[0] + 1))
        )])
        
        self.update_fig(fig, 'Root Mean Square Inner Product', 'run 1', 'run 2')
        
        return fig

    def create_fel_plot(
        self,
        matrix_vals: List[List[float]],
        x_coords: List[float],
        y_coords: List[float],
        pc_pair: Tuple[int, int],
        overlay: Optional[Tuple[List[float], List[float]]] = None,
        colorscale: str = "Viridis",
        contours: int = 7,
        title: Optional[str] = None,
        z_title: str = "G (kJ/mol)",
    ) -> go.Figure:
        """Create FEL heatmap + surface subplot with optional projection overlay."""
        fig = make_subplots(
            rows=1,
            cols=2,
            specs=[[{'type': 'heatmap'}, {'type': 'surface'}]],
            column_widths=[1, 1]
        )

        z_vals = np.array(matrix_vals)
        if z_vals.size == 0:
            return fig

        zmin = float(np.nanmin(z_vals))
        zmax = float(np.nanmax(z_vals))

        fig.add_trace(
            go.Heatmap(
                z=z_vals[::-1],
                x=x_coords,
                y=y_coords,
                colorscale=colorscale,
                coloraxis='coloraxis',
                zsmooth='best'
            ),
            row=1,
            col=1
        )

        fig.add_trace(
            go.Surface(
                z=z_vals[::-1],
                x=x_coords,
                y=y_coords,
                colorscale=colorscale,
                coloraxis='coloraxis',
                showscale=False
            ),
            row=1,
            col=2
        )

        if overlay:
            overlay_x, overlay_y = overlay
            fig.add_trace(
                go.Scatter(
                    x=overlay_x,
                    y=overlay_y,
                    mode='lines+markers',
                    marker_size=5,
                    line_width=0.6,
                    line_color='rgba(255,255,255,0.7)',
                    marker_colorscale='greys',
                    marker_showscale=False,
                    showlegend=False,
                ),
                row=1,
                col=1
            )

        fig.update_traces(
            contours_z=dict(
                project_z=True,
                show=True,
                usecolormap=True,
                size=round((zmax - zmin) / max(contours, 1), 2),
                start=zmin,
                end=zmax,
            ),
            selector=dict(type='surface')
        )

        pc1, pc2 = pc_pair
        fig.update_layout(
            scene_camera_eye=dict(x=1.9, y=1.9, z=1.5),
            coloraxis_colorbar_title=z_title,
            coloraxis_colorscale=colorscale,
            height=self.config.plot_config.height,
            width=self.config.plot_config.width,
            scene=dict(
                xaxis_title=f'PC{pc1}',
                yaxis_title=f'PC{pc2}',
                zaxis_title=z_title,
                xaxis_nticks=5,
                yaxis_nticks=5,
                zaxis_nticks=2,
            ),
        )

        title_text = title or f"FEL on PC{pc1} and PC{pc2}"
        self.update_fig(fig, title_text, f'PC{pc1}', f'PC{pc2}', y=0.9)

        return fig

    def create_xpm_heatmap(
        self,
        matrix_vals: List[List[float]],
        x_coords: List[float],
        y_coords: List[float],
        title: str,
        x_label: str,
        y_label: str,
        colorscale: str = "Viridis",
    ) -> go.Figure:
        """Create heatmap from XPM data."""
        fig = go.Figure(
            data=[go.Heatmap(
                z=np.array(matrix_vals)[::-1],
                x=x_coords,
                y=y_coords,
                colorscale=colorscale,
                coloraxis='coloraxis',
                zsmooth='best'
            )]
        )
        self.update_fig(fig, title, x_label, y_label, y=0.9)
        return fig
    
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
                if isinstance(e, StaticExportUnavailable):
                    logger.debug("Static export unavailable; skipped %s", file_path)
                    continue
                logger.error(f"Error saving {fmt} plot {file_path}: {e}")

    def _save_svg_fallback(self, file_path: Path, fallback: Dict) -> bool:
        kind = fallback.get('kind')
        try:
            if kind == 'scree':
                return self._save_scree_svg(file_path, fallback.get('eigenvals'))
            if kind == 'cosine':
                return self._save_cosine_svg(
                    file_path,
                    fallback.get('cc_apo'),
                    fallback.get('cc_holo'),
                    fallback.get('apo_name'),
                    fallback.get('holo_name')
                )
            if kind == 'hist':
                return self._save_hist_svg(file_path, fallback.get('proj_data'))
            if kind == 'line':
                return self._save_line_svg(file_path, fallback.get('proj_data'))
            if kind == 'proj2d':
                return self._save_proj2d_svg(
                    file_path,
                    fallback.get('proj_apo'),
                    fallback.get('proj_holo'),
                    fallback.get('pc_pair'),
                    fallback.get('apo_name'),
                    fallback.get('holo_name')
                )
            if kind == 'proj3d':
                return self._save_proj3d_svg(
                    file_path,
                    fallback.get('proj_apo'),
                    fallback.get('proj_holo'),
                    fallback.get('apo_name'),
                    fallback.get('holo_name')
                )
            if kind == 'rmsip':
                return self._save_rmsip_svg(file_path, fallback.get('xpm_data'))
        except Exception as exc:
            logger.error(f"Matplotlib fallback failed for {file_path}: {exc}")
        return False

    def _base_fig(self, width: Optional[int] = None, height: Optional[int] = None):
        _, plt = init_matplotlib()
        fig, ax = plt.subplots(
            figsize=((width or self.config.plot_config.width) / 100,
                     (height or self.config.plot_config.height) / 100),
            dpi=100 * self.config.plot_config.scale
        )
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size
        return plt, fig, ax

    def _save_scree_svg(self, file_path: Path, eigenvals_data: Optional[Dict]) -> bool:
        if not eigenvals_data:
            return False

        xs = eigenvals_data['x'][:self.config.n_pcs]
        ys = eigenvals_data['y'][:self.config.n_pcs]
        variances = self.calculate_variance(eigenvals_data['y'])[:self.config.n_pcs]

        _, plt, ax1 = self._base_fig()
        ax2 = ax1.twinx()

        ax1.plot(xs, ys, marker='o', label='Individual Eigenvalues')
        ax2.plot(xs, variances, marker='o', color=self.holo_color, linestyle=':',
                 label='Cumulative Contribution to Total Variance')

        ax1.set_xlabel('Eigenvector Number')
        ax1.set_ylabel('Eigenvalue (nm^2)')
        ax2.set_ylabel('Proportion (%)')

        ax1.set_xlim(1, self.config.n_pcs)
        ax2.set_ylim(0, 100)
        if ys:
            ax1.set_ylim(0, ys[0] + 0.1 * ys[0])

        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', frameon=False)

        ax1.set_title('Change of Eigenvalues with Eigenvectors')
        plt.tight_layout()
        plt.savefig(file_path, format='svg', bbox_inches='tight')
        plt.close()
        return True

    def _save_cosine_svg(self, file_path: Path, cc_apo: Optional[Dict],
                         cc_holo: Optional[Dict], apo_name: str, holo_name: str) -> bool:
        if not cc_apo or not cc_holo:
            return False

        xs_apo = cc_apo['x']
        ys_apo = cc_apo['y']
        xs_holo = cc_holo['x']
        ys_holo = cc_holo['y']

        _, plt, ax = self._base_fig()
        width = 0.4
        ax.bar([x - width / 2 for x in xs_holo], ys_holo, width=width, label=holo_name)
        ax.bar([x + width / 2 for x in xs_apo], ys_apo, width=width, label=apo_name)
        ax.set_xlim(0.5, self.config.n_pcs + 0.5)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Eigenvector Number')
        ax.set_ylabel('Cosine content')
        ax.set_title(cc_holo.get('title') or 'Cosine content')
        ax.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(file_path, format='svg', bbox_inches='tight')
        plt.close()
        return True

    def _save_hist_svg(self, file_path: Path, proj_data: Optional[Dict]) -> bool:
        if not proj_data:
            return False

        titles = proj_data['titles'][:self.config.n_pcs]
        ys = proj_data['ys'][:self.config.n_pcs]

        _, plt = init_matplotlib()
        fig, axes = plt.subplots(
            2, 5,
            figsize=(self.config.plot_config.width / 100, self.config.plot_config.height / 100),
            dpi=100 * self.config.plot_config.scale
        )
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size - 2

        for idx, y_vals in enumerate(ys[:self.config.n_pcs]):
            row = idx // 5
            col = idx % 5
            ax = axes[row][col]
            ax.hist(y_vals, bins=30, density=True, color=self.holo_color, alpha=0.6, edgecolor='white')
            ax.set_title(f"PC {idx + 1}")
            ax.tick_params(axis='both', labelsize=8)

        fig.suptitle('Histogram of the first 10 PCs')
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        fig.savefig(file_path, format='svg', bbox_inches='tight')
        plt.close(fig)
        return True

    def _save_line_svg(self, file_path: Path, proj_data: Optional[Dict]) -> bool:
        if not proj_data:
            return False

        xs = proj_data['xs'][:self.config.n_pcs]
        ys = proj_data['ys'][:self.config.n_pcs]

        _, plt = init_matplotlib()
        fig, axes = plt.subplots(
            5, 2,
            figsize=(self.config.plot_config.width / 100, self.config.plot_config.height / 100),
            dpi=100 * self.config.plot_config.scale
        )
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size - 2

        for idx, (x_vals, y_vals) in enumerate(zip(xs, ys)):
            if idx >= self.config.n_pcs:
                break
            row = idx // 2
            col = idx % 2
            ax = axes[row][col]
            x_plot = np.arange(len(x_vals)) / 10
            ax.plot(x_plot, y_vals, linewidth=1)
            ax.set_title(f"PC {idx + 1}")
            ax.tick_params(axis='both', labelsize=8)

        fig.suptitle('Projection of the trajectory on the first 10 PCs')
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        fig.savefig(file_path, format='svg', bbox_inches='tight')
        plt.close(fig)
        return True

    def _save_proj2d_svg(self, file_path: Path, proj_apo: Optional[Dict],
                         proj_holo: Optional[Dict], pc_pair: Optional[Tuple[int, int]],
                         apo_name: str, holo_name: str) -> bool:
        if not proj_apo or not proj_holo or not pc_pair:
            return False

        x1 = proj_apo['x']
        y1 = proj_apo['y']
        x2 = proj_holo['x']
        y2 = proj_holo['y']
        pc1, pc2 = pc_pair

        _, plt, ax = self._base_fig()
        ax.scatter(x1, y1, s=8, alpha=0.6, label=apo_name, color=APO_HEX)
        ax.scatter(x2, y2, s=8, alpha=0.6, label=holo_name, color=self.holo_color)

        ax.scatter([sum(x1) / len(x1)], [sum(y1) / len(y1)], s=90, marker='D', color=mix(APO_HEX, '#000000', 0.25), edgecolors='white', label=f'{apo_name} mean')
        ax.scatter([sum(x2) / len(x2)], [sum(y2) / len(y2)], s=90, marker='D', color=mix(self.holo_color, '#000000', 0.25), edgecolors='white', label=f'{holo_name} mean')

        ax.set_title(f'Bidimensional projection of the trajectory along PC {pc1} and PC {pc2}')
        ax.set_xlabel(f'PC{pc1}')
        ax.set_ylabel(f'PC{pc2}')
        ax.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(file_path, format='svg', bbox_inches='tight')
        plt.close()
        return True

    def _save_proj3d_svg(self, file_path: Path, proj_apo: Optional[Tuple],
                         proj_holo: Optional[Tuple], apo_name: str, holo_name: str) -> bool:
        if not proj_apo or not proj_holo:
            return False

        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        x1, y1, z1 = proj_apo
        x2, y2, z2 = proj_holo

        _, plt = init_matplotlib()
        fig = plt.figure(figsize=(self.config.plot_config.width / 100,
                                  self.config.plot_config.height / 100),
                         dpi=100 * self.config.plot_config.scale)
        plt.rcParams['font.family'] = self.config.plot_config.font_family
        plt.rcParams['font.size'] = self.config.plot_config.font_size - 2
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(x1, y1, z1, s=6, alpha=0.6, label=apo_name, color=APO_HEX)
        ax.scatter(x2, y2, z2, s=6, alpha=0.6, label=holo_name, color=self.holo_color)
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_zlabel('PC3')
        ax.set_title('3D Projection of Trajectory')
        ax.legend(frameon=False)

        fig.tight_layout()
        fig.savefig(file_path, format='svg', bbox_inches='tight')
        plt.close(fig)
        return True

    def _save_rmsip_svg(self, file_path: Path, xpm_data: Optional[np.ndarray]) -> bool:
        if xpm_data is None:
            return False

        _, plt, ax = self._base_fig()
        ax.imshow(xpm_data, cmap='Greys', origin='lower')
        ax.set_title('Root Mean Square Inner Product')
        ax.set_xlabel('run 2')
        ax.set_ylabel('run 1')
        plt.tight_layout()
        plt.savefig(file_path, format='svg', bbox_inches='tight')
        plt.close()
        return True
