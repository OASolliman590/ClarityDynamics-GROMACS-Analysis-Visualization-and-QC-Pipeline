"""
MD Analyzer Module
==================

Main orchestrator for MD trajectory analysis.
Coordinates data processing and plotting.

Workflow:
1. Create output directories
2. Find and process data files
3. Create visualizations
4. Save statistics
5. Generate summary report
"""

import json
import logging
from pathlib import Path
from typing import Dict, Optional, List
import pandas as pd

from .config import MDConfig
from .data_processor import MDDataProcessor
from .plotter import MDPlotter

logger = logging.getLogger(__name__)


class MDAnalyzer:
    """
    Main orchestrator for MD trajectory analysis.
    
    Coordinates the complete workflow:
    - Data processing (via MDDataProcessor)
    - Visualization (via MDPlotter)
    - Statistics export
    - Summary reports
    
    Attributes:
        config: MDConfig object
        processor: MDDataProcessor instance
        plotter: MDPlotter instance
    """
    
    def __init__(self, config: MDConfig):
        """
        Initialize analyzer with configuration.
        
        Args:
            config: MDConfig object
        """
        self.config = config
        self.processor = MDDataProcessor(config)
        self.plotter = MDPlotter(config)
        
        logger.info("MDAnalyzer initialized")
    
    def run_analysis(self) -> Dict:
        """
        Run complete analysis pipeline.
        
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
            - n_metrics: Number of metrics processed
        """
        logger.info("="*60)
        logger.info("Starting MD Trajectory Analysis")
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
            logger.info(f"    {system_name}: {len(system_data)} metrics")
        
        # Step 3: Create all plots
        logger.info("\n[3/5] Creating visualizations...")
        plots_created = self.plotter.create_all_plots(all_data)
        logger.info(f"  Created {len(plots_created)} plots")

        # Step 4: Export processed metric data
        logger.info("\n[4/6] Exporting processed metric data...")
        data_exports = self._export_metric_data_csv(all_data)
        logger.info(f"  Exported {len(data_exports)} data file(s)")

        if self.config.stability.enabled:
            self._write_equilibration(self.processor.equilibration_records)

        convergence_file = self._write_convergence(all_data)

        # Step 5: Save statistics
        logger.info("\n[5/6] Saving statistics...")
        statistics_file = self._save_statistics_csv(all_data)
        if statistics_file:
            logger.info(f"  Statistics saved to: {statistics_file}")

        # Step 6: Generate summary report
        logger.info("\n[6/6] Generating summary report...")
        summary_file = self._generate_summary_report(all_data, plots_created)
        if summary_file:
            logger.info(f"  Summary report saved to: {summary_file}")
        
        # Compile results
        results = {
            'success': True,
            'output_dir': self.config.output_dir,
            'plots_dir': self.config.output_dir / 'plots',
            'plots_created': plots_created,
            'data_exports': data_exports,
            'statistics_file': statistics_file,
            'convergence_file': convergence_file,
            'summary_file': summary_file,
            'n_systems': len(all_data),
            'n_metrics': sum(len(system_data) for system_data in all_data.values())
        }
        
        logger.info("\n" + "="*60)
        logger.info("Analysis Complete!")
        logger.info("="*60)
        logger.info(f"Output directory: {results['output_dir']}")
        logger.info(f"Plots created: {len(plots_created)}")
        logger.info(f"Systems analyzed: {results['n_systems']}")
        logger.info(f"Data files exported: {len(data_exports)}")
        logger.info("="*60 + "\n")
        
        return results
    
    def _create_output_dirs(self):
        """Create output directory structure"""
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        (self.config.output_dir / 'plots').mkdir(parents=True, exist_ok=True)
        (self.config.output_dir / 'data').mkdir(parents=True, exist_ok=True)
        
        logger.debug(f"  Output directory: {self.config.output_dir}")
        logger.debug(f"  Plots directory: {self.config.output_dir / 'plots'}")
    
    def _write_equilibration(self, records: List[Dict]) -> Optional[Path]:
        """Per-system production start (time, index, N_eff) for reuse by RMSF/MM/PBSA."""
        path = self.config.output_dir / "data" / "equilibration.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(records, indent=2) + "\n")
        for record in records:
            if record.get("status") == "detected":
                logger.info(
                    "  %s: production starts at %s %s (index %s, N_eff=%s)",
                    record["system"], record["start_time"], record.get("time_unit") or "",
                    record["start_index"],
                    f"{record['n_effective']:.0f}" if record.get("n_effective") else "n/a",
                )
        return path

    def _write_convergence(self, all_data: Dict) -> Optional[Path]:
        """Per system/metric/replicate sampling checks over the analysed (production) window."""
        from ..utils.convergence import series_convergence  # noqa: WPS433

        rows = []
        for system_name, system_data in all_data.items():
            for metric_name, metric_data in system_data.items():
                if self.processor.is_profile_metric(metric_name):
                    continue  # per-residue profiles are not time series
                if metric_data.get('mode') == 'replicate':
                    series = metric_data.get('replicates_raw') or metric_data.get('replicates') or []
                else:
                    series = [metric_data.get('values')] if metric_data.get('values') is not None else []
                for index, values in enumerate(series, start=1):
                    row = {'System': system_name, 'Metric': metric_name, 'Replicate': index}
                    row.update(series_convergence(values, self.config.convergence))
                    rows.append(row)
        if not rows:
            return None
        frame = pd.DataFrame(rows)
        path = self.config.output_dir / 'data' / 'convergence.csv'
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, index=False, float_format='%.6g')
        flagged = frame[frame['Flags'].astype(str).str.len() > 0]
        if not flagged.empty:
            logger.warning("  Sampling checks flagged %d of %d replicate series (see %s)",
                           len(flagged), len(frame), path)
        return path

    def _save_statistics_csv(self, all_data: Dict) -> Optional[Path]:
        """
        Save comprehensive statistics to CSV.
        
        Args:
            all_data: Complete data structure from processor
            
        Returns:
            Path to statistics file or None
        """
        rows = []
        
        for system_name, system_data in all_data.items():
            for metric_name, metric_data in system_data.items():
                stats = metric_data.get('stats', {})
                
                row = {
                    'System': system_name,
                    'Metric': metric_name,
                    'Mode': metric_data.get('mode', 'unknown'),
                    'N_Points': metric_data.get('n_points', 0),
                    'Mean': stats.get('mean', 0),
                    'Std': stats.get('std', 0),
                    'Min': stats.get('min', 0),
                    'Max': stats.get('max', 0),
                    'Median': stats.get('median', 0),
                    'Q1': stats.get('q1', 0),
                    'Q3': stats.get('q3', 0)
                }

                time_axis = metric_data.get("time")
                if time_axis is not None and len(time_axis) > 0:
                    row["Time_Start"] = float(time_axis[0])
                    row["Time_End"] = float(time_axis[-1])

                time_meta = metric_data.get("time_metadata") or {}
                row["Time_Input_Unit"] = time_meta.get("input_unit")
                row["Time_Inferred_Unit"] = time_meta.get("inferred_input_unit")
                row["Time_Output_Unit"] = time_meta.get("output_unit")
                row["Time_Conversion_Scale"] = time_meta.get("applied_scale")
                row["Time_Step_ps"] = time_meta.get("time_step_ps")
                row["Time_Conversion_Mode"] = time_meta.get("conversion_mode")
                
                # Add replicate info if available
                if 'n_replicates' in metric_data:
                    row['N_Replicates'] = metric_data['n_replicates']
                if 'std_mean' in stats:
                    row['Std_Mean'] = stats['std_mean']
                # Uncertainty of the mean: between-replicate SEM for replicates,
                # block-averaged SEM (autocorrelation-aware) for single runs.
                rep_means = stats.get('replicate_means')
                row['Replicate_Means'] = (
                    ";".join(f"{float(v):.6g}" for v in rep_means) if rep_means else None
                )
                row['SEM'] = stats.get('sem')
                row['SEM_Method'] = stats.get('sem_method')
                row['Between_Replicate_SD'] = stats.get('between_replicate_sd')
                row['N_Effective'] = stats.get('n_effective')
                row['Stable_Start_Time'] = metric_data.get('stable_start_time')
                unit_meta = metric_data.get('unit_metadata') or {}
                row['Value_Unit_Declared'] = unit_meta.get('declared')
                row['Value_Unit'] = unit_meta.get('display')
                row['Value_Unit_Factor'] = unit_meta.get('factor')
                row['Value_Unit_Status'] = unit_meta.get('status')

                rows.append(row)
        
        if not rows:
            logger.warning("  No statistics to save")
            return None
        
        # Create DataFrame and save
        df = pd.DataFrame(rows)
        stats_file = self.config.output_dir / 'statistics_summary.csv'
        df.to_csv(stats_file, index=False, float_format='%.4f')
        
        return stats_file

    def _export_metric_data_csv(self, all_data: Dict) -> List[Path]:
        """
        Export per-system, per-metric processed arrays into `output_dir/data`.

        Returns:
            List of exported CSV paths.
        """
        data_dir = self.config.output_dir / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        exported: List[Path] = []

        rmsf_labels: Optional[List[str]] = None

        for system_name, system_data in all_data.items():
            for metric_name, metric_data in system_data.items():
                mode = metric_data.get("mode")
                filename = f"{system_name}_{metric_name}.csv"
                output_path = data_dir / filename

                if metric_name == "rmsf":
                    if rmsf_labels is None:
                        try:
                            rmsf_labels = self.processor.prepare_amino_acid_axis()
                        except Exception:
                            rmsf_labels = None

                if mode == "single":
                    x_vals = metric_data.get("time")
                    y_vals = metric_data.get("values")
                    if x_vals is None or y_vals is None:
                        continue

                    frame_df = pd.DataFrame({"x": x_vals, "value": y_vals})
                    if metric_name == "rmsf":
                        frame_df.rename(columns={"x": "residue_index", "value": "rmsf"}, inplace=True)
                        if rmsf_labels:
                            frame_df["residue_label"] = self._align_labels(rmsf_labels, len(frame_df))
                    frame_df.to_csv(output_path, index=False)
                    exported.append(output_path)
                    continue

                if mode == "replicate":
                    x_vals = metric_data.get("time")
                    mean_vals = metric_data.get("mean")
                    upper_vals = metric_data.get("upper")
                    lower_vals = metric_data.get("lower")
                    if x_vals is None or mean_vals is None:
                        continue

                    frame_df = pd.DataFrame(
                        {
                            "x": x_vals,
                            "mean": mean_vals,
                            "upper": upper_vals,
                            "lower": lower_vals,
                        }
                    )
                    replicates = metric_data.get("replicates") or []
                    for idx, rep_values in enumerate(replicates, start=1):
                        frame_df[f"rep{idx}"] = rep_values[: len(frame_df)]

                    if metric_name == "rmsf":
                        frame_df.rename(columns={"x": "residue_index"}, inplace=True)
                        if rmsf_labels:
                            frame_df["residue_label"] = self._align_labels(rmsf_labels, len(frame_df))
                    frame_df.to_csv(output_path, index=False)
                    exported.append(output_path)

        return exported

    @staticmethod
    def _align_labels(labels: List[str], target_len: int) -> List[str]:
        """Align residue labels to data length without crashing on mismatch."""
        if not labels:
            return []
        if len(labels) == target_len:
            return labels
        if len(labels) > target_len:
            return labels[:target_len]
        # Pad with numeric strings if labels are shorter than RMSF rows.
        padded = list(labels)
        start = len(labels)
        for idx in range(start, target_len):
            padded.append(str(idx))
        return padded
    
    def _generate_summary_report(self, all_data: Dict, plots_created: List) -> Optional[Path]:
        """
        Generate HTML summary report.
        
        Args:
            all_data: Complete data structure
            plots_created: List of created plot names
            
        Returns:
            Path to summary report or None
        """
        try:
            report_file = self.config.output_dir / 'ANALYSIS_SUMMARY.md'
            
            with open(report_file, 'w') as f:
                # Header
                f.write(f"# MD Analysis Summary - {self.config.protein_name}\n\n")
                f.write("---\n\n")
                
                # Configuration
                f.write("## Configuration\n\n")
                f.write(f"- **Protein:** {self.config.protein_name}\n")
                f.write(f"- **Base Directory:** {self.config.base_dir}\n")
                f.write(f"- **Output Directory:** {self.config.output_dir}\n")
                f.write(f"- **Visualization Style:** {self.config.plot_config.style}\n")
                f.write(f"- **Plot Template:** {self.config.plot_config.template}\n\n")
                f.write(f"- **Time Input Unit:** {self.config.time_unit_input}\n")
                f.write(f"- **Time Output Unit:** {self.config.time_unit_output}\n")
                if self.config.time_step_ps is not None:
                    f.write(f"- **Time Step (ps):** {self.config.time_step_ps}\n")
                if self.config.time_scale is not None:
                    f.write(f"- **Explicit Time Scale:** {self.config.time_scale}\n")
                f.write("\n")
                
                # Systems
                f.write("## Systems Analyzed\n\n")
                for system in self.config.systems:
                    apo_status = "APO" if system.is_apo else "HOLO"
                    f.write(f"- **{system.name}** ({apo_status}): {system.replicates} replicate(s)\n")
                f.write("\n")
                
                # Metrics
                f.write("## Metrics Processed\n\n")
                for metric in self.config.metrics:
                    f.write(f"- **{metric.title}** ({metric.name})\n")
                f.write("\n")
                
                # Results
                f.write("## Results\n\n")
                f.write(f"- **Total Systems:** {len(all_data)}\n")
                f.write(f"- **Total Plots:** {len(plots_created)}\n")
                f.write(f"- **Statistics File:** `statistics_summary.csv`\n\n")
                
                # Plots
                f.write("## Generated Plots\n\n")
                for plot_name in sorted(plots_created):
                    f.write(f"- `{plot_name}.html` / `{plot_name}.svg`\n")
                f.write("\n")
                
                # Statistics Summary
                f.write("## Statistics Summary\n\n")
                f.write("| System | Metric | Mean | SEM | Std | Min | Max |\n")
                f.write("|--------|--------|------|-----|-----|-----|-----|\n")

                for system_name, system_data in all_data.items():
                    for metric_name, metric_data in system_data.items():
                        stats = metric_data.get('stats', {})
                        sem = stats.get('sem')
                        sem_text = f"{sem:.3f}" if sem is not None else "n/a"
                        f.write(f"| {system_name} | {metric_name} | "
                               f"{stats.get('mean', 0):.3f} | "
                               f"{sem_text} | "
                               f"{stats.get('std', 0):.3f} | "
                               f"{stats.get('min', 0):.3f} | "
                               f"{stats.get('max', 0):.3f} |\n")

                f.write("\n")
                f.write(
                    "Std is the spread of the observable over frames (pooled over replicates). "
                    "SEM is the uncertainty of the mean: between-replicate SD/sqrt(n) when "
                    "replicates exist, otherwise a block-averaged estimate that accounts for "
                    "frame-to-frame correlation. Compare systems using SEM, not Std.\n\n"
                )
                
                # Footer
                f.write("---\n\n")
                f.write("Generated by GROMACS Analysis Toolkit\n")
            
            return report_file
            
        except Exception as e:
            logger.error(f"  Could not generate summary report: {e}")
            return None
    
    def get_data(self) -> Dict:
        """
        Get processed data without running full analysis.
        Useful for programmatic access.
        
        Returns:
            Complete processed data structure
        """
        return self.processor.load_all_data()
    
    def create_custom_plot(self, metric_name: str, systems: Optional[List[str]] = None):
        """
        Create a custom plot for specific metric and systems.
        
        Args:
            metric_name: Name of metric to plot
            systems: Optional list of system names (None = all)
        
        Returns:
            Plotly Figure object
        """
        # Load data
        all_data = self.processor.load_all_data()
        
        # Filter systems if specified
        if systems:
            filtered_data = {k: v for k, v in all_data.items() if k in systems}
        else:
            filtered_data = all_data
        
        # Extract metric data
        metric_data = {}
        for system_name, system_data in filtered_data.items():
            if metric_name in system_data:
                metric_data[system_name] = system_data[metric_name]
        
        if not metric_data:
            raise ValueError(f"No data found for metric: {metric_name}")
        
        # Find metric object
        metric = None
        for m in self.config.metrics:
            if m.name == metric_name:
                metric = m
                break
        
        if not metric:
            raise ValueError(f"Metric configuration not found: {metric_name}")
        
        # Create plot
        x_values = None
        if metric.name == 'rmsf':
            if self.config.amino_acids:
                x_values = [str(aa) for aa in self.config.amino_acids]
            else:
                try:
                    x_values = self.processor.prepare_amino_acid_axis()
                except Exception:
                    x_values = None
        
        return self.plotter.create_plot(metric_data, metric, x_values)
