"""
MMPBSA Analyzer Module
======================

Main orchestrator for MMPBSA analysis.
Coordinates data processing and plotting.
"""

import logging
import re
from pathlib import Path
from typing import Dict, Optional, List
import pandas as pd

from .config import MMPBSAConfig
from .processor import MMPBSAProcessor
from .plotter import MMPBSAPlotter
from ..utils.publication import export_publication_figure, write_metadata, write_source_table

logger = logging.getLogger(__name__)


def _safe_publication_token(value: str) -> str:
    token = re.sub(r'[^A-Za-z0-9_.-]+', '_', str(value)).strip('._')
    return token or 'system'


class MMPBSAAnalyzer:
    """
    Main orchestrator for MMPBSA analysis.
    
    Coordinates the complete workflow:
    - Data processing (via MMPBSAProcessor)
    - Visualization (via MMPBSAPlotter)
    - Statistics export
    - Summary reports
    
    Attributes:
        config: MMPBSAConfig object
        processor: MMPBSAProcessor instance
        plotter: MMPBSAPlotter instance
    """
    
    def __init__(self, config: MMPBSAConfig):
        """
        Initialize analyzer with configuration.
        
        Args:
            config: MMPBSAConfig object
        """
        self.config = config
        self.processor = MMPBSAProcessor(config)
        self.plotter = MMPBSAPlotter(config)
        
        logger.info("MMPBSAAnalyzer initialized")
    
    def run_analysis(self) -> Dict:
        """
        Run complete MMPBSA analysis pipeline.
        
        Workflow:
        1. Create output directories
        2. Find and process all data
        3. Create all plots
        4. Save statistics CSV
        5. Generate summary report
        
        Returns:
            Dict with results:
            - output_dir: Path to output directory
            - plots_dir: Path to plots directory
            - plots_created: List of created plot names
            - statistics_file: Path to statistics CSV
            - n_systems: Number of systems analyzed
        """
        logger.info("="*60)
        logger.info("Starting MMPBSA Analysis")
        logger.info("="*60)
        
        # Step 1: Create output directories
        logger.info("\n[1/5] Creating output directories...")
        self._create_output_dirs()
        
        # Step 2: Find and process all data
        logger.info("\n[2/5] Finding and processing data files...")
        all_data = self.processor.load_all_data()
        
        if not all_data:
            logger.error("No data found! Check your configuration and data directories.")
            return {
                'success': False,
                'error': 'No data found'
            }
        
        logger.info(f"  Processed {len(all_data)} systems")
        for system_name, system_data in all_data.items():
            has_results = system_data.get('results') is not None
            has_decomp = system_data.get('decomp') is not None
            logger.info(f"    {system_name}: Results={has_results}, Decomp={has_decomp}")
        
        # Step 3: Create all plots
        logger.info("\n[3/5] Creating visualizations...")
        plots_created = []
        
        plots_dir = self.config.output_dir / 'plots'
        plots_dir.mkdir(parents=True, exist_ok=True)
        
        for system_name, system_data in all_data.items():
            results_data = system_data.get('results')
            decomp_data = system_data.get('decomp')
            
            # Component plot
            if results_data:
                try:
                    fig = self.plotter.create_component_plot(results_data, system_name)
                    self.plotter.save_plot(
                        fig,
                        plots_dir,
                        f"components_{system_name}",
                        fallback={
                            'kind': 'components',
                            'results_data': results_data,
                            'system_name': system_name
                        }
                    )
                    plots_created.append(f"components_{system_name}")
                    logger.info(f"  ✓ Created component plot for {system_name}")
                except Exception as e:
                    logger.error(f"  ✗ Error creating component plot for {system_name}: {e}")
            
            # Entropy-corrected binding energy (only when IE/C2 entropy was computed)
            if results_data and results_data.get('entropy'):
                try:
                    fig = self.plotter.create_entropy_plot(results_data, system_name)
                    if fig is not None:
                        self.plotter.save_plot(fig, plots_dir, f"entropy_{system_name}")
                        plots_created.append(f"entropy_{system_name}")
                        logger.info(f"  ✓ Created entropy plot for {system_name}")
                except Exception as e:
                    logger.error(f"  ✗ Error creating entropy plot for {system_name}: {e}")
            elif results_data:
                logger.info(
                    f"  No interaction/C2 entropy found for {system_name} "
                    "(enable entropy in mmpbsa.in or keep FINAL_RESULTS*.dat next to CSV results)."
                )

            # Decomposition plot
            if decomp_data is not None and not decomp_data.empty and results_data:
                try:
                    fig = self.plotter.create_decomp_plot(decomp_data, results_data, system_name)
                    self.plotter.save_plot(
                        fig,
                        plots_dir,
                        f"decomp_{system_name}",
                        fallback={
                            'kind': 'decomp',
                            'decomp_data': decomp_data,
                            'results_data': results_data,
                            'system_name': system_name
                        }
                    )
                    plots_created.append(f"decomp_{system_name}")
                    logger.info(f"  ✓ Created decomposition plot for {system_name}")
                except Exception as e:
                    logger.error(f"  ✗ Error creating decomposition plot for {system_name}: {e}")

        # Cross-ligand comparison
        if self.config.compare_systems and len(all_data) > 1:
            try:
                fig = self.plotter.create_binding_energy_comparison(
                    all_data,
                    binding_key=self.config.compare_binding_key
                )

                systems = []
                means = []
                stds = []
                for system_name, system_data in all_data.items():
                    results_data = system_data.get('results')
                    if not results_data:
                        continue
                    mean_val, std_val = self.plotter._extract_binding_value(
                        results_data,
                        self.config.compare_binding_key
                    )
                    if mean_val is None:
                        continue
                    systems.append(system_name)
                    means.append(mean_val)
                    stds.append(std_val)

                self.plotter.save_plot(
                    fig,
                    plots_dir,
                    "binding_energy_compare",
                    fallback={
                        'kind': 'binding_compare',
                        'systems': systems,
                        'means': means,
                        'stds': stds
                    }
                )
                plots_created.append("binding_energy_compare")
                logger.info("  ✓ Created binding energy comparison plot")
            except Exception as e:
                logger.error(f"  ✗ Error creating binding energy comparison plot: {e}")

        # Entropy comparison across systems
        if self.config.compare_systems and len(all_data) > 1:
            try:
                fig = self.plotter.create_entropy_comparison(all_data)
                if fig is not None:
                    self.plotter.save_plot(fig, plots_dir, "entropy_compare")
                    plots_created.append("entropy_compare")
                    logger.info("  ✓ Created entropy comparison plot")
            except Exception as e:
                logger.error(f"  ✗ Error creating entropy comparison plot: {e}")

        entropy_rows = []
        for system_name, system_data in all_data.items():
            results_data = system_data.get('results')
            if results_data:
                for row in self.plotter.entropy_terms(results_data):
                    entropy_rows.append({'System': system_name, **{k: v for k, v in row.items() if k != 'key'}})
        if entropy_rows:
            entropy_file = self.config.output_dir / 'data' / 'entropy_summary.csv'
            entropy_file.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(entropy_rows).rename(columns={
                'method': 'Method', 'total': 'dTOTAL', 'total_err': 'dTOTAL_err',
                'minus_tds': 'minus_TdS', 'minus_tds_err': 'minus_TdS_err',
                'dg': 'dG_bind', 'dg_err': 'dG_bind_err', 'sigma': 'sigma_Eint', 'reliable': 'Reliable',
            }).to_csv(entropy_file, index=False)

        # Component-by-component comparison
        if self.config.compare_components and len(all_data) > 1:
            try:
                fig = self.plotter.create_components_comparison(
                    all_data,
                    component_order=self.config.compare_component_order
                )

                components = []
                systems_series: Dict[str, Dict[str, List[float]]] = {}

                # Resolve components list for fallback data
                order = self.config.compare_component_order or [
                    'ΔVDWAALS', 'ΔEEL', 'ΔEGB', 'ΔESURF', 'ΔGGAS', 'ΔGSOLV', 'ΔTOTAL'
                ]
                for comp in order:
                    if any(comp in (data.get('results') or {}).get('delta_components', {}) for data in all_data.values()):
                        components.append(comp)

                for system_name, system_data in all_data.items():
                    results_data = system_data.get('results')
                    if not results_data:
                        continue
                    delta_components = results_data.get('delta_components', {})
                    means = []
                    stds = []
                    for comp in components:
                        comp_data = delta_components.get(comp)
                        if comp_data:
                            means.append(comp_data.get('mean'))
                            stds.append(comp_data.get('std', 0.0))
                        else:
                            means.append(None)
                            stds.append(0.0)
                    systems_series[system_name] = {"means": means, "stds": stds}

                self.plotter.save_plot(
                    fig,
                    plots_dir,
                    "components_compare",
                    fallback={
                        'kind': 'components_compare',
                        'components': components,
                        'series': systems_series,
                    }
                )
                plots_created.append("components_compare")
                logger.info("  ✓ Created components comparison plot")
            except Exception as e:
                logger.error(f"  ✗ Error creating components comparison plot: {e}")
        
        publication_outputs = self._write_publication_outputs(all_data)

        # Step 4: Save parsed data as CSV
        logger.info("\n[4/6] Exporting parsed data to CSV...")
        data_files = self._save_parsed_data_csv(all_data)
        
        # Step 5: Save statistics
        logger.info("\n[5/6] Exporting statistics...")
        statistics_file = self._save_statistics_csv(all_data)
        
        # Step 6: Generate summary report
        logger.info("\n[6/6] Generating summary report...")
        summary_file = self._generate_summary_report(all_data, plots_created)
        
        logger.info("\n" + "="*60)
        logger.info("MMPBSA Analysis Complete!")
        logger.info("="*60)
        logger.info(f"Output directory: {self.config.output_dir}")
        logger.info(f"Plots created: {len(plots_created)}")
        
        return {
            'success': True,
            'output_dir': str(self.config.output_dir),
            'plots_dir': str(plots_dir),
            'plots_created': plots_created,
            'statistics_file': str(statistics_file) if statistics_file else None,
            'summary_file': str(summary_file) if summary_file else None,
            'n_systems': len(all_data),
            'publication_outputs': publication_outputs,
        }

    def _write_publication_outputs(self, all_data: Dict) -> Dict:
        """Write additive publication figures, source tables, and availability metadata."""
        publication = getattr(self.config, 'publication', None)
        if publication is None or not publication.enabled:
            return {}
        output_dir = self.config.output_dir / publication.subdirectory
        output_dir.mkdir(parents=True, exist_ok=True)
        outputs: Dict[str, object] = {}
        availability: Dict[str, object] = {}

        component_payload = self.processor.publication_component_payload(all_data)
        if component_payload.empty:
            availability['components'] = {
                'status': 'unavailable', 'reason': 'No supported component summaries were parsed.'
            }
        else:
            missing_components = {
                system: component_payload.loc[
                    (component_payload['system'] == system) & (component_payload['status'] == 'unavailable'),
                    ['component', 'reason'],
                ].to_dict(orient='records')
                for system in component_payload['system'].unique()
            }
            source = write_source_table(component_payload, output_dir, 'components_source')
            figure = self.plotter.create_publication_components_figure(component_payload)
            exports = export_publication_figure(figure, output_dir, 'components', publication)
            outputs['components'] = {'source': source, 'formats': exports}
            availability['components'] = {
                'status': 'available', 'source': source, 'formats': exports,
                'missing_components': missing_components,
                'provenance': component_payload.attrs.get('provenance', {}),
            }

        residue_payloads = {}
        for system_index, (system_name, system_data) in enumerate(all_data.items(), start=1):
            provenance = (system_data.get('qc') or {}).get('results_files', [])
            results = system_data.get('results') or {}
            sampling = (publication.sampling_metadata or {}).get(system_name)
            energy = self.processor.publication_energy_payload(
                results, window=publication.moving_window, sampling_metadata=sampling,
            )
            if energy.get('status') != 'available':
                availability[f'saved_total:{system_name}'] = energy
            else:
                sample_data = energy['samples'].copy()
                safe_name = _safe_publication_token(system_name)
                artifact_stem = f'saved_total_{system_index}_{safe_name}'
                energy_source = write_source_table(sample_data, output_dir, f'{artifact_stem}_source')
                figure = self.plotter.create_publication_energy_figure(
                    sample_data, system_name, window=publication.moving_window, x_label=energy['x_label'],
                )
                exports = export_publication_figure(figure, output_dir, artifact_stem, publication)
                outputs[f'saved_total:{system_name}'] = {'source': energy_source, 'formats': exports}
                availability[f'saved_total:{system_name}'] = {
                    'status': 'available', 'source': energy_source, 'formats': exports,
                    'x_axis': energy['x_label'], 'moving_window': publication.moving_window,
                    'physical_time': energy['physical_time'],
                    'invalid_sample_count': energy.get('invalid_sample_count', 0),
                    'invalid_samples': energy.get('invalid_samples', []),
                    'overall_mean': energy['overall_mean'], 'provenance': provenance,
                    'replicate_provenance': results.get('replicate_sources', {})
                    if isinstance(energy['overall_mean'], dict) else {},
                }

            decomp = system_data.get('decomp')
            if decomp is None or decomp.empty:
                availability[f'residues:{system_name}'] = {
                    'status': 'unavailable', 'reason': 'No residue decomposition input was parsed.'
                }
                continue
            selected = self.processor.publication_residue_payload(decomp, top_n=publication.top_n)
            if selected.empty:
                availability[f'residues:{system_name}'] = {
                    'status': 'unavailable', 'reason': 'Parsed decomposition has no numeric TOTAL residue contributions.',
                    'provenance': (system_data.get('qc') or {}).get('decomp_files', []),
                }
                continue
            selected.insert(0, 'System', system_name)
            residue_payloads[system_name] = selected
            availability[f'residues:{system_name}'] = {
                'status': 'available', 'top_n': publication.top_n,
                'selection': 'all' if publication.top_n == 'all' else f'most_negative_{publication.top_n}',
                'provenance': (system_data.get('qc') or {}).get('decomp_files', []),
            }

        if residue_payloads:
            residue_source_data = pd.concat(residue_payloads.values(), ignore_index=True)
            source = write_source_table(residue_source_data, output_dir, 'residues_source')
            plots_payloads = {name: df.drop(columns=['System']) for name, df in residue_payloads.items()}
            figure = self.plotter.create_publication_residue_figure(plots_payloads)
            exports = export_publication_figure(figure, output_dir, 'residues', publication)
            outputs['residues'] = {'source': source, 'formats': exports}
            availability['residues'] = {'status': 'available', 'source': source, 'formats': exports}
        else:
            availability['residues'] = {
                'status': 'unavailable', 'reason': 'No residue decomposition inputs were available.'
            }

        metadata = {
            'schema_version': 1,
            'description': 'MM/GBSA endpoint decomposition and saved-sample estimates; not experimental binding free energies.',
            'units': 'kcal/mol',
            'moving_window': publication.moving_window,
            'top_n': publication.top_n,
            'sampling_metadata': publication.sampling_metadata,
            'systems': {
                system_name: {
                    'source_files': (system_data.get('results') or {}).get('source_files', []),
                    'replicate_sources': (system_data.get('results') or {}).get('replicate_sources', {}),
                }
                for system_name, system_data in all_data.items()
            },
            'availability': availability,
        }
        metadata_path = write_metadata(metadata, output_dir, 'metadata')
        outputs['metadata'] = metadata_path
        outputs['availability'] = availability
        return outputs
    
    def _create_output_dirs(self):
        """Create output directory structure."""
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        (self.config.output_dir / 'plots').mkdir(parents=True, exist_ok=True)
        (self.config.output_dir / 'data').mkdir(parents=True, exist_ok=True)
    
    def _save_parsed_data_csv(self, all_data: Dict) -> List[Path]:
        """
        Save parsed data to CSV files in the data folder.
        
        Args:
            all_data: Processed data from processor
            
        Returns:
            List of saved file paths
        """
        data_dir = self.config.output_dir / 'data'
        data_dir.mkdir(parents=True, exist_ok=True)
        saved_files = []
        
        try:
            for system_name, system_data in all_data.items():
                # Save results data
                results_data = system_data.get('results')
                qc_data = system_data.get('qc') or {}
                if results_data:
                    # Create DataFrame from delta components
                    rows = []
                    for comp_name, comp_data in results_data.get('delta_components', {}).items():
                        rows.append({
                            'Component': comp_name,
                            'Mean': comp_data['mean'],
                            'Std': comp_data['std'],
                            'SEM': comp_data.get('sem'),
                            'Uncertainty': comp_data.get('uncertainty'),
                            'N': comp_data['n']
                        })
                    
                    # Add binding energy
                    for key, data in results_data.get('binding_energy', {}).items():
                        rows.append({
                            'Component': f'Binding_{key}',
                            'Mean': data['mean'],
                            'Std': data['std'],
                            'SEM': data.get('sem'),
                            'Uncertainty': data.get('uncertainty'),
                            'N': data['n']
                        })
                    
                    # Add entropy
                    for key, data in results_data.get('entropy', {}).items():
                        rows.append({
                            'Component': f'Entropy_{key}',
                            'Mean': data['mean'],
                            'Std': data['std'],
                            'SEM': data.get('sem'),
                            'Uncertainty': data.get('uncertainty'),
                            'N': data['n']
                        })

                    if rows and qc_data:
                        for row in rows:
                            row.update({
                                'QC_ExpectedReplicates': qc_data.get('expected_replicates'),
                                'QC_FoundResultsReplicates': qc_data.get('found_results_replicates'),
                                'QC_MissingResultsReplicates': qc_data.get('missing_results_replicates'),
                                'QC_ReplicateComplete': qc_data.get('replicate_complete'),
                                'QC_GroupIdentityOK': qc_data.get('group_identity_ok'),
                            })

                    if rows:
                        df_results = pd.DataFrame(rows)
                        results_file = data_dir / f"{system_name}_results.csv"
                        df_results.to_csv(results_file, index=False)
                        saved_files.append(results_file)
                        logger.info(f"  ✓ Saved results CSV: {results_file}")

                    frame_data = results_data.get('frame_data')
                    if frame_data is not None and not frame_data.empty:
                        frame_file = data_dir / f"{system_name}_frame_results.csv"
                        frame_data.to_csv(frame_file, index=False)
                        saved_files.append(frame_file)
                        logger.info(f"  ✓ Saved per-frame results CSV: {frame_file}")

                if qc_data:
                    qc_file = data_dir / f"{system_name}_qc.csv"
                    pd.DataFrame([qc_data]).to_csv(qc_file, index=False)
                    saved_files.append(qc_file)
                    logger.info(f"  ✓ Saved QC metadata CSV: {qc_file}")
                
                # Save decomposition data
                decomp_data = system_data.get('decomp')
                if decomp_data is not None and not decomp_data.empty:
                    decomp_file = data_dir / f"{system_name}_decomposition.csv"
                    decomp_data.to_csv(decomp_file, index=False)
                    saved_files.append(decomp_file)
                    logger.info(f"  ✓ Saved decomposition CSV: {decomp_file}")
            
            return saved_files
            
        except Exception as e:
            logger.error(f"Error saving parsed data to CSV: {e}")
            return saved_files
    
    def _save_statistics_csv(self, all_data: Dict) -> Optional[Path]:
        """
        Save statistics to CSV file.
        
        Args:
            all_data: Processed data from processor
            
        Returns:
            Path to statistics file
        """
        try:
            stats_data = []
            
            for system_name, system_data in all_data.items():
                results_data = system_data.get('results')
                if not results_data:
                    continue
                
                # Delta components
                for comp_name, comp_data in results_data.get('delta_components', {}).items():
                    stats_data.append({
                        'System': system_name,
                        'Component': comp_name,
                        'Mean': comp_data['mean'],
                        'Std': comp_data['std'],
                        'SEM': comp_data.get('sem'),
                        'Uncertainty': comp_data.get('uncertainty'),
                        'N': comp_data['n']
                    })
                
                # Binding energy
                for key, data in results_data.get('binding_energy', {}).items():
                    stats_data.append({
                        'System': system_name,
                        'Component': f'Binding_{key}',
                        'Mean': data['mean'],
                        'Std': data['std'],
                        'SEM': data.get('sem'),
                        'Uncertainty': data.get('uncertainty'),
                        'N': data['n']
                    })
                
                # Entropy
                for key, data in results_data.get('entropy', {}).items():
                    stats_data.append({
                        'System': system_name,
                        'Component': f'Entropy_{key}',
                        'Mean': data['mean'],
                        'Std': data['std'],
                        'SEM': data.get('sem'),
                        'Uncertainty': data.get('uncertainty'),
                        'N': data['n']
                    })
            
            if stats_data:
                df = pd.DataFrame(stats_data)
                stats_file = self.config.output_dir / 'statistics_summary.csv'
                df.to_csv(stats_file, index=False)
                logger.info(f"  ✓ Statistics saved to: {stats_file}")
                return stats_file
            
            return None
            
        except Exception as e:
            logger.error(f"Error saving statistics: {e}")
            return None
    
    def _generate_summary_report(self, all_data: Dict, plots_created: List) -> Optional[Path]:
        """
        Generate summary report in Markdown format.
        
        Args:
            all_data: Processed data from processor
            plots_created: List of created plot names
            
        Returns:
            Path to summary file
        """
        try:
            report_lines = [
                f"# MMPBSA Analysis Summary - {self.config.protein_name}",
                "",
                "---",
                "",
                "## Configuration",
                "",
                f"- **Protein:** {self.config.protein_name}",
                f"- **Ligand:** {self.config.ligand_name}",
                f"- **Base Directory:** {self.config.base_dir}",
                f"- **Output Directory:** {self.config.output_dir}",
                "",
                "## Systems Analyzed",
                ""
            ]
            
            for system in self.config.systems:
                mode = "replicate" if system.replicates > 1 else "single"
                report_lines.append(f"- **{system.name}**: {system.replicates} replicate(s) ({mode})")
            
            report_lines.extend([
                "",
                "## Results Summary",
                ""
            ])
            
            for system_name, system_data in all_data.items():
                results_data = system_data.get('results')
                qc_data = system_data.get('qc') or {}
                if results_data:
                    binding_energy = results_data.get('binding_energy', {})
                    total = binding_energy.get('TOTAL', {})
                    
                    if total:
                        mean = total.get('mean', 0.0)
                        std = total.get('std', 0.0)
                        report_lines.append(f"### {system_name}")
                        sem = total.get('sem')
                        spread = (
                            f"± {sem:.2f} (SEM over {total.get('n')} replicates)"
                            if total.get('uncertainty') == 'between_replicates' and sem is not None
                            else f"(frame SD {std:.2f}; single run, no replicate uncertainty)"
                        )
                        report_lines.append(
                            f"- **ΔTOTAL (effective binding energy, no entropy):** {mean:.2f} {spread} kcal/mol"
                        )
                        for key, label in (('interaction_entropy', 'IE'), ('c2_entropy', 'C2')):
                            ent = results_data.get('entropy', {}).get(key)
                            if ent:
                                report_lines.append(
                                    f"- **ΔG_bind ≈ ΔTOTAL + (−TΔS, {label}):** "
                                    f"{mean + ent.get('mean', 0.0):.2f} kcal/mol"
                                )
                        if qc_data:
                            found = qc_data.get('found_results_replicates')
                            expected = qc_data.get('expected_replicates')
                            missing = qc_data.get('missing_results_replicates')
                            report_lines.append(
                                f"- **Replicate completeness:** {found}/{expected} (missing: {missing})"
                            )
                            report_lines.append(
                                f"- **Group identity check:** {qc_data.get('group_identity_ok')}"
                            )
                        report_lines.append("")
            
            report_lines.extend([
                "> MM/PBSA and MM/GBSA ΔTOTAL values are end-point estimates suited to ranking "
                "related ligands; they are not absolute binding affinities and typically "
                "overestimate their magnitude. Entropy corrections (IE/C2) are noisy and "
                "should be interpreted with care.",
                "",
                "## Generated Plots",
                ""
            ])
            
            for plot_name in plots_created:
                report_lines.append(f"- `{plot_name}.html` / `{plot_name}.svg`")
            
            report_lines.extend([
                "",
                "---",
                "",
                f"Generated by GROMACS Analysis Toolkit"
            ])
            
            report_file = self.config.output_dir / 'ANALYSIS_SUMMARY.md'
            report_file.write_text('\n'.join(report_lines))
            logger.info(f"  ✓ Summary report saved to: {report_file}")
            return report_file
            
        except Exception as e:
            logger.error(f"Error generating summary report: {e}")
            return None
    
    def get_data(self) -> Dict:
        """
        Get processed data (for programmatic access).
        
        Returns:
            Dictionary with all processed data
        """
        return self.processor.load_all_data()
