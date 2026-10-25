"""
PCA Analyzer Module
===================

Main orchestrator for PCA analysis.
Auto-detects files and creates all visualizations.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, List

from .config import PCAConfig
from .processor import PCAProcessor
from .plotter import PCAPlotter
from .runner import PCAGromacsRunner

logger = logging.getLogger(__name__)


class PCAAnalyzer:
    """
    Main PCA analysis orchestrator.
    
    Automatically:
    1. Detects PCA files in directory structure
    2. Parses all data files
    3. Creates all visualizations
    4. Saves plots to output directory
    """
    
    def __init__(self, config: PCAConfig):
        """
        Initialize analyzer.
        
        Args:
            config: PCAConfig object
        """
        self.config = config
        self.processor = PCAProcessor(config)
        self.plotter = PCAPlotter(config)
        logger.info("PCAAnalyzer initialized")
    
    def run_analysis(self) -> Dict:
        """
        Run complete PCA analysis.
        
        Returns:
            Dict with analysis results:
            {
                'success': bool,
                'output_dir': str,
                'plots_created': List[str],
                'ligand_name': str,
                'gromacs_run': Dict
            }
        """
        self._create_output_dirs()
        
        logger.info("="*60)
        logger.info("Starting PCA Analysis")
        logger.info("="*60)
        logger.info(f"Base directory: {self.config.base_dir}")
        logger.info(f"Output directory: {self.config.output_dir}")
        logger.info(f"Protein: {self.config.protein_name}")

        gromacs_result = None
        if self.config.run_gromacs.enabled:
            logger.info("\n[0/4] Running GROMACS PCA preparation...")
            gromacs_result = PCAGromacsRunner(self.config).run()
        
        # Auto-detect files
        logger.info("\n[1/4] Auto-detecting PCA files...")
        detected_files = self.processor.auto_detect_files()
        
        # Extract ligand name if found
        ligand_name = detected_files.get('ligand_name')
        if ligand_name:
            logger.info(f"Detected ligand name: {ligand_name}")
        
        # Load and parse data
        logger.info("\n[2/4] Loading and parsing data...")
        all_data = self.processor.load_all_data(detected_files)
        
        # Create visualizations
        logger.info("\n[3/4] Creating visualizations...")
        plots_created = []
        
        # 1. Scree plot
        if 'eigenvals' in all_data:
            try:
                fig = self.plotter.create_scree_plot(all_data['eigenvals'])
                self.plotter.save_plot(
                    fig,
                    self.config.output_dir,
                    "scree_plot",
                    fallback={
                        'kind': 'scree',
                        'eigenvals': all_data['eigenvals']
                    }
                )
                plots_created.append("scree_plot")
                logger.info("  ✓ Created scree plot")
            except Exception as e:
                logger.error(f"  ✗ Error creating scree plot: {e}")
        
        # 2. Cosine content
        if 'cc_holo' in all_data and 'cc_apo' in all_data:
            try:
                apo_name = f"{self.config.protein_name} Apo"
                holo_name = f"{self.config.protein_name} Holo"
                if ligand_name:
                    holo_name = f"{self.config.protein_name}-{ligand_name} Holo"
                
                fig = self.plotter.create_cosine_content_plot(
                    all_data['cc_apo'],
                    all_data['cc_holo'],
                    apo_name,
                    holo_name
                )
                self.plotter.save_plot(
                    fig,
                    self.config.output_dir,
                    "CC_all",
                    fallback={
                        'kind': 'cosine',
                        'cc_apo': all_data['cc_apo'],
                        'cc_holo': all_data['cc_holo'],
                        'apo_name': apo_name,
                        'holo_name': holo_name
                    }
                )
                plots_created.append("CC_all")
                logger.info("  ✓ Created cosine content plot")
            except Exception as e:
                logger.error(f"  ✗ Error creating cosine content plot: {e}")
        
        # 3. PC histograms
        if 'proj_1d' in all_data:
            try:
                fig = self.plotter.create_pc_histograms(all_data['proj_1d'])
                self.plotter.save_plot(
                    fig,
                    self.config.output_dir,
                    "hist",
                    fallback={
                        'kind': 'hist',
                        'proj_data': all_data['proj_1d']
                    }
                )
                plots_created.append("hist")
                logger.info("  ✓ Created PC histograms")
            except Exception as e:
                logger.error(f"  ✗ Error creating PC histograms: {e}")
        
        # 4. PC time series
        if 'proj_1d' in all_data:
            try:
                fig = self.plotter.create_pc_timeseries(all_data['proj_1d'])
                self.plotter.save_plot(
                    fig,
                    self.config.output_dir,
                    "line",
                    fallback={
                        'kind': 'line',
                        'proj_data': all_data['proj_1d']
                    }
                )
                plots_created.append("line")
                logger.info("  ✓ Created PC time series")
            except Exception as e:
                logger.error(f"  ✗ Error creating PC time series: {e}")
        
        # 5. 2D projections
        apo_name = f"{self.config.protein_name} Apo"
        holo_name = f"{self.config.protein_name} Holo"
        if ligand_name:
            holo_name = f"{self.config.protein_name}-{ligand_name} Holo"
        
        for pair_key, pc_pair in [('12', (1, 2)), ('13', (1, 3)), ('23', (2, 3))]:
            if pair_key in all_data.get('proj_2d_apo', {}) and pair_key in all_data.get('proj_2d_holo', {}):
                try:
                    fig = self.plotter.create_2d_projection(
                        all_data['proj_2d_apo'][pair_key],
                        all_data['proj_2d_holo'][pair_key],
                        pc_pair,
                        apo_name,
                        holo_name
                    )
                    self.plotter.save_plot(
                        fig,
                        self.config.output_dir,
                        f"proj_{pair_key}",
                        fallback={
                            'kind': 'proj2d',
                            'proj_apo': all_data['proj_2d_apo'][pair_key],
                            'proj_holo': all_data['proj_2d_holo'][pair_key],
                            'pc_pair': pc_pair,
                            'apo_name': apo_name,
                            'holo_name': holo_name
                        }
                    )
                    plots_created.append(f"proj_{pair_key}")
                    logger.info(f"  ✓ Created 2D projection PC{pc_pair[0]}-PC{pc_pair[1]}")
                except Exception as e:
                    logger.error(f"  ✗ Error creating 2D projection PC{pc_pair[0]}-PC{pc_pair[1]}: {e}")
        
        # 6. 3D projection.
        if 'proj_3d_apo' in all_data and 'proj_3d_holo' in all_data:
            try:
                fig = self.plotter.create_3d_projection(
                    all_data['proj_3d_apo'],
                    all_data['proj_3d_holo'],
                    apo_name,
                    holo_name
                )
                # Use the standard projection output name.
                self.plotter.save_plot(
                    fig,
                    self.config.output_dir,
                    "proj_123",
                    fallback={
                        'kind': 'proj3d',
                        'proj_apo': all_data['proj_3d_apo'],
                        'proj_holo': all_data['proj_3d_holo'],
                        'apo_name': apo_name,
                        'holo_name': holo_name
                    }
                )
                plots_created.append("proj_123")
                logger.info("  ✓ Created 3D projection")
            except Exception as e:
                logger.error(f"  ✗ Error creating 3D projection: {e}")

        # 7. Free Energy Landscape (FEL) plots
        if self.config.fel.enabled and 'fel' in all_data:
            for pair_key in self.config.fel.pairs:
                pair_key = str(pair_key)
                fel_entry = all_data['fel'].get(pair_key, {})
                if not fel_entry:
                    continue
                pc_pair = (int(pair_key[0]), int(pair_key[1]))

                overlay = None
                if self.config.fel.overlay_projection:
                    proj = all_data.get('proj_eig', {}).get(pair_key)
                    if proj and len(proj.get('ys', [])) >= 2:
                        step = max(1, int(self.config.fel.step))
                        overlay = (
                            proj['ys'][0][::step],
                            proj['ys'][1][::step],
                        )

                if 'gibbs' in fel_entry:
                    try:
                        fel_data = fel_entry['gibbs']
                        fig = self.plotter.create_fel_plot(
                            fel_data['matrix'],
                            fel_data['x_coords'],
                            fel_data['y_coords'],
                            pc_pair,
                            overlay=overlay,
                            colorscale=self.config.fel.colorscale,
                            contours=self.config.fel.contours,
                            title=f"2D and 3D projection of FEL on PC{pc_pair[0]} and PC{pc_pair[1]}",
                            z_title="G (kJ/mol)"
                        )
                        name = f"FEL_{pair_key}"
                        self.plotter.save_plot(fig, self.config.output_dir, name)
                        plots_created.append(name)
                        logger.info(f"  ✓ Created FEL plot PC{pc_pair[0]}-PC{pc_pair[1]}")
                    except Exception as e:
                        logger.error(f"  ✗ Error creating FEL plot {pair_key}: {e}")

                if self.config.fel.include_probability and 'probability' in fel_entry:
                    try:
                        prob_data = fel_entry['probability']
                        fig = self.plotter.create_xpm_heatmap(
                            prob_data['matrix'],
                            prob_data['x_coords'],
                            prob_data['y_coords'],
                            f"Probability distribution projected on PC{pc_pair[0]} and PC{pc_pair[1]}",
                            f"PC{pc_pair[0]}",
                            f"PC{pc_pair[1]}",
                            colorscale=self.config.fel.colorscale
                        )
                        name = f"FEL_prob_{pair_key}"
                        self.plotter.save_plot(fig, self.config.output_dir, name)
                        plots_created.append(name)
                    except Exception as e:
                        logger.error(f"  ✗ Error creating FEL probability plot {pair_key}: {e}")

                if self.config.fel.include_entropy and 'entropy' in fel_entry:
                    try:
                        ent_data = fel_entry['entropy']
                        fig = self.plotter.create_xpm_heatmap(
                            ent_data['matrix'],
                            ent_data['x_coords'],
                            ent_data['y_coords'],
                            f"Entropy projected on PC{pc_pair[0]} and PC{pc_pair[1]}",
                            f"PC{pc_pair[0]}",
                            f"PC{pc_pair[1]}",
                            colorscale=self.config.fel.colorscale
                        )
                        name = f"FEL_entropy_{pair_key}"
                        self.plotter.save_plot(fig, self.config.output_dir, name)
                        plots_created.append(name)
                    except Exception as e:
                        logger.error(f"  ✗ Error creating FEL entropy plot {pair_key}: {e}")

                if self.config.fel.include_enthalpy and 'enthalpy' in fel_entry:
                    try:
                        enth_data = fel_entry['enthalpy']
                        fig = self.plotter.create_xpm_heatmap(
                            enth_data['matrix'],
                            enth_data['x_coords'],
                            enth_data['y_coords'],
                            f"Enthalpy projected on PC{pc_pair[0]} and PC{pc_pair[1]}",
                            f"PC{pc_pair[0]}",
                            f"PC{pc_pair[1]}",
                            colorscale=self.config.fel.colorscale
                        )
                        name = f"FEL_enthalpy_{pair_key}"
                        self.plotter.save_plot(fig, self.config.output_dir, name)
                        plots_created.append(name)
                    except Exception as e:
                        logger.error(f"  ✗ Error creating FEL enthalpy plot {pair_key}: {e}")

        # 8. Clustering on PCA projections
        if self.config.clustering.enabled:
            pair_key = str(self.config.clustering.pair)
            proj_source = None
            if pair_key in all_data.get('proj_2d_holo', {}):
                proj_source = all_data['proj_2d_holo'][pair_key]
                cluster_label = "Holo"
            elif pair_key in all_data.get('proj_2d_apo', {}):
                proj_source = all_data['proj_2d_apo'][pair_key]
                cluster_label = "Apo"

            if proj_source:
                try:
                    points = list(zip(proj_source['x'], proj_source['y']))
                    step = max(1, int(self.config.clustering.downsample))
                    if step > 1:
                        points = points[::step]
                    if self.config.clustering.max_points:
                        points = points[: int(self.config.clustering.max_points)]
                    labels, centroids = self._cluster_points(points)
                    if labels is not None:
                        fig = self._cluster_plot(points, labels, pair_key, cluster_label)
                        name = f"clusters_{pair_key}"
                        self.plotter.save_plot(fig, self.config.output_dir, name)
                        plots_created.append(name)
                        self._save_cluster_csv(points, labels, pair_key)
                        self._save_cluster_representatives(points, labels, centroids, pair_key)
                        logger.info(f"  ✓ Created clustering plot for PC{pair_key[0]}-PC{pair_key[1]}")
                except Exception as e:
                    logger.error(f"  ✗ Error creating clustering for {pair_key}: {e}")
        
        publication = None
        if self.config.publication.enabled:
            from .publication import save_publication
            payload = self.processor.build_publication_payload(all_data, detected_files)
            publication = save_publication(payload, self.config.output_dir, self.config.publication, self.plotter)

        # Generate summary
        logger.info("\n[4/4] Generating summary...")
        summary_file = self._generate_summary_report(plots_created, ligand_name)
        
        logger.info("\n" + "="*60)
        logger.info("PCA Analysis Complete!")
        logger.info("="*60)
        logger.info(f"Output directory: {self.config.output_dir}")
        logger.info(f"Plots created: {len(plots_created)}")
        if ligand_name:
            logger.info(f"Ligand name: {ligand_name}")
        
        return {
            'success': True,
            'output_dir': str(self.config.output_dir),
            'plots_created': plots_created,
            'ligand_name': ligand_name,
            'summary_file': str(summary_file) if summary_file else None,
            'gromacs_run': gromacs_result,
            'publication': publication,
        }
    
    def _create_output_dirs(self):
        """Create output directories if they don't exist."""
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Output directory: {self.config.output_dir}")

    def _cluster_points(self, points: List[tuple]):
        """Cluster 2D points using configured method."""
        if not points:
            return None, None

        method = self.config.clustering.method.lower()
        if method == "dbscan":
            try:
                from sklearn.cluster import DBSCAN  # noqa: WPS433
                import numpy as np

                arr = np.array(points)
                model = DBSCAN(
                    eps=self.config.clustering.eps,
                    min_samples=self.config.clustering.min_samples,
                ).fit(arr)
                labels = model.labels_.tolist()
                centroids = self._compute_centroids(arr, labels)
                return labels, centroids
            except Exception as exc:
                logger.warning(f"DBSCAN unavailable or failed ({exc}); falling back to kmeans.")

        return self._kmeans_cluster(points)

    def _kmeans_cluster(self, points: List[tuple]):
        import numpy as np

        arr = np.array(points)
        n_clusters = max(1, int(self.config.clustering.n_clusters))
        if len(arr) < n_clusters:
            n_clusters = max(1, len(arr))

        # Simple k-means
        rng = np.random.default_rng(42)
        centroids = arr[rng.choice(len(arr), size=n_clusters, replace=False)]

        for _ in range(50):
            distances = np.linalg.norm(arr[:, None, :] - centroids[None, :, :], axis=2)
            labels = distances.argmin(axis=1)
            new_centroids = []
            for k in range(n_clusters):
                members = arr[labels == k]
                if len(members) == 0:
                    new_centroids.append(centroids[k])
                else:
                    new_centroids.append(members.mean(axis=0))
            new_centroids = np.array(new_centroids)
            if np.allclose(new_centroids, centroids):
                break
            centroids = new_centroids

        return labels.tolist(), centroids

    def _compute_centroids(self, arr, labels):
        import numpy as np

        centroids = []
        for label in sorted(set(labels)):
            if label == -1:
                continue
            members = arr[np.array(labels) == label]
            if len(members) == 0:
                centroids.append(None)
            else:
                centroids.append(members.mean(axis=0))
        return centroids

    def _cluster_plot(self, points: List[tuple], labels: List[int], pair_key: str, cluster_label: str):
        import plotly.graph_objects as go

        pc1 = int(pair_key[0])
        pc2 = int(pair_key[1])
        fig = go.Figure()

        unique_labels = sorted(set(labels))
        for lab in unique_labels:
            xs = [p[0] for i, p in enumerate(points) if labels[i] == lab]
            ys = [p[1] for i, p in enumerate(points) if labels[i] == lab]
            name = f"Cluster {lab}" if lab != -1 else "Noise"
            fig.add_trace(go.Scatter(
                x=xs,
                y=ys,
                mode='markers',
                name=name
            ))

        fig.update_layout(
            title=f"{self.config.protein_name} {cluster_label} · PCA Clusters (PC{pc1}-PC{pc2})",
            xaxis_title=f"PC{pc1}",
            yaxis_title=f"PC{pc2}",
            font_family=self.config.plot_config.font_family,
            font_size=self.config.plot_config.font_size,
            template=self.config.plot_config.template
        )
        return fig

    def _save_cluster_csv(self, points: List[tuple], labels: List[int], pair_key: str) -> None:
        import pandas as pd

        data = {
            "frame": list(range(len(points))),
            "pc1": [p[0] for p in points],
            "pc2": [p[1] for p in points],
            "cluster": labels,
        }
        df = pd.DataFrame(data)
        out_path = self.config.output_dir / f"clusters_{pair_key}.csv"
        df.to_csv(out_path, index=False)

    def _save_cluster_representatives(self, points: List[tuple], labels: List[int], centroids, pair_key: str) -> None:
        import numpy as np
        import pandas as pd

        arr = np.array(points)
        reps = []
        for idx, centroid in enumerate(centroids):
            if centroid is None:
                continue
            distances = np.linalg.norm(arr - centroid, axis=1)
            rep_idx = int(np.argmin(distances))
            reps.append({
                "cluster": idx,
                "frame": rep_idx,
                "pc1": points[rep_idx][0],
                "pc2": points[rep_idx][1],
            })

        if reps:
            out_path = self.config.output_dir / f"cluster_representatives_{pair_key}.csv"
            pd.DataFrame(reps).to_csv(out_path, index=False)
    
    def _generate_summary_report(self, plots_created: List[str], ligand_name: Optional[str]) -> Optional[Path]:
        """
        Generate summary report.
        
        Args:
            plots_created: List of plot names created
            ligand_name: Detected ligand name
            
        Returns:
            Path to summary file or None
        """
        try:
            summary_file = self.config.output_dir / "PCA_analysis_summary.md"
            
            with open(summary_file, 'w') as f:
                f.write("# PCA Analysis Summary\n\n")
                f.write(f"**Protein:** {self.config.protein_name}\n")
                if ligand_name:
                    f.write(f"**Ligand:** {ligand_name}\n")
                f.write(f"**Base Directory:** {self.config.base_dir}\n")
                f.write(f"**Output Directory:** {self.config.output_dir}\n\n")
                
                f.write("## Plots Created\n\n")
                for plot in plots_created:
                    f.write(f"- {plot}\n")
                
                f.write("\n## Files\n\n")
                f.write("All plots are saved in multiple formats (HTML, SVG, PNG) as configured.\n")
            
            logger.info(f"  ✓ Generated summary report: {summary_file}")
            return summary_file
            
        except Exception as e:
            logger.error(f"Error generating summary report: {e}")
            return None
