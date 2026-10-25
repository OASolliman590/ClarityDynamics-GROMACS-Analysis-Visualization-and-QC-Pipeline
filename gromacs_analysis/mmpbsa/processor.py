"""
MMPBSA Processor Module
=======================

Handles data processing, replicate averaging, and SD calculation.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np

from .config import MMPBSAConfig, MMPBSASystemConfig
from .parser import MMPBSAParser, _canonical_component_name, _canonical_component_payload

logger = logging.getLogger(__name__)


class MMPBSAProcessor:
    """
    Process MMPBSA data files and calculate statistics.
    
    Handles:
    - Finding data files
    - Parsing CSV/DAT files
    - Replicate averaging
    - SD calculation
    """
    
    def __init__(self, config: MMPBSAConfig):
        """
        Initialize processor with configuration.
        
        Args:
            config: MMPBSAConfig object
        """
        self.config = config
        self.parser = MMPBSAParser(debug=False)
        self._frames_used: List[Optional[int]] = []
        logger.info("MMPBSAProcessor initialized")

    @staticmethod
    def publication_component_payload(all_data: Dict[str, Dict]) -> pd.DataFrame:
        """Return long-form component means and dispersion with scientific grouping."""
        component_order = ('VDWAALS', 'EEL', 'EGB', 'ESURF', 'GGAS', 'GSOLV', 'TOTAL')
        legacy_names = {name: f'Δ{name}' for name in component_order}
        rows = []
        for system_name, system_data in all_data.items():
            results = system_data.get('results') or {}
            components = results.get('canonical_components') or results.get('delta_components') or {}
            for name in component_order:
                value = components.get(name, components.get(legacy_names[name]))
                group = 'raw' if name in component_order[:4] else 'grouped' if name in ('GGAS', 'GSOLV') else 'total'
                if not value:
                    rows.append({
                        'system': system_name, 'component': name, 'group': group, 'mean': None,
                        'sd': None, 'uncertainty': None, 'status': 'unavailable',
                        'n': 0, 'reason': f'{name} was not present in the parsed results.',
                    })
                    continue
                uncertainty = value.get('uncertainty')
                if uncertainty == 'between_replicates' and value.get('n', 0) > 1:
                    uncertainty = 'replicate_sd'
                elif results.get('mode') == 'replicate' and value.get('n', 0) <= 1:
                    uncertainty = 'unavailable_single_replicate'
                elif uncertainty == 'frame_sd':
                    uncertainty = 'saved_sample_sd'
                elif not uncertainty:
                    uncertainty = 'saved_sample_sd' if results.get('mode') == 'single' else 'unavailable_single_replicate'
                spread = value.get('std', value.get('sample_sd', value.get('std_dev', 0.0)))
                if uncertainty == 'unavailable_single_replicate':
                    spread = None
                mean = value.get('mean', value.get('value'))
                if mean is None or not np.isfinite(float(mean)) or (
                    spread is not None and not np.isfinite(float(spread))
                ):
                    rows.append({
                        'system': system_name, 'component': name, 'group': group, 'mean': None,
                        'sd': None, 'n': value.get('n', 0), 'uncertainty': None,
                        'status': 'unavailable', 'reason': 'Component mean or dispersion is not finite.',
                    })
                    continue
                rows.append({
                    'system': system_name, 'component': name, 'group': group,
                    'mean': mean,
                    'sd': spread,
                    'n': value.get('n', 0),
                    'uncertainty': uncertainty,
                    'status': 'available', 'reason': None,
                })
        return pd.DataFrame(rows, columns=['system', 'component', 'group', 'mean', 'sd', 'n', 'uncertainty', 'status', 'reason'])

    @staticmethod
    def publication_energy_payload(
        results: Dict, window: int = 11, sampling_metadata: Optional[Dict] = None,
    ) -> Dict[str, object]:
        """Prepare saved TOTAL samples and centered mean without inferring time."""
        if isinstance(window, bool) or not isinstance(window, int) or window < 1 or window % 2 == 0:
            raise ValueError('moving window must be a positive odd integer')
        frame_data = results.get('frame_data')
        if frame_data is None or frame_data.empty:
            return {'status': 'unavailable', 'reason': 'No saved frame-level results are available.'}
        total_column = next((c for c in frame_data.columns if _canonical_component_name(c) == 'TOTAL'), None)
        if total_column is None:
            return {'status': 'unavailable', 'reason': 'Saved frame-level results do not contain TOTAL.'}
        all_total_values = pd.to_numeric(frame_data[total_column], errors='coerce')
        invalid_sample_count = int((~np.isfinite(all_total_values.to_numpy())).sum())
        if np.isinf(all_total_values.to_numpy()).any():
            return {
                'status': 'unavailable',
                'reason': 'Saved TOTAL samples include non-finite values; no trace was produced.',
                'invalid_sample_count': invalid_sample_count,
            }
        groups = list(frame_data.groupby('Replicate', sort=False)) if 'Replicate' in frame_data.columns else [(None, frame_data)]
        sample_frames = []
        for replicate, group in groups:
            total_values = pd.to_numeric(group[total_column], errors='coerce')
            valid_rows = total_values.notna()
            total = total_values.loc[valid_rows].reset_index(drop=True)
            if total.empty:
                continue
            samples = pd.DataFrame({
                'sample_index': np.arange(1, len(total) + 1), 'total': total,
                'moving_mean': total.rolling(window=window, center=True, min_periods=1).mean().values,
                'overall_mean': float(total.mean()),
            })
            if 'Frame #' in group.columns:
                samples['frame'] = group.loc[valid_rows, 'Frame #'].reset_index(drop=True).to_numpy()
            if replicate is not None:
                samples['replicate'] = replicate
            sample_frames.append(samples)
        if not sample_frames:
            return {
                'status': 'unavailable', 'reason': 'Saved TOTAL samples contain no finite values.',
                'invalid_sample_count': invalid_sample_count,
            }
        samples = pd.concat(sample_frames, ignore_index=True)
        axis = 'Saved sample index'
        metadata = sampling_metadata or {}
        physical_time = {
            'status': 'unavailable',
            'reason': 'Validated sampling metadata was not supplied.',
        }
        if metadata.get('validated') is True:
            times = metadata.get('time')
            frames = metadata.get('frames')
            source = metadata.get('source')
            if not isinstance(source, str) or not source.strip():
                physical_time['reason'] = 'Validated sampling metadata must identify its source.'
            elif frames is None or times is None or not metadata.get('unit') or len(times) != len(frames):
                physical_time['reason'] = 'Validated sampling metadata requires matching frames, time values, and a unit.'
            elif not pd.Series(frames).is_unique:
                physical_time['reason'] = 'Frame identifiers must be unique for a one-to-one time mapping.'
            else:
                time_values = pd.to_numeric(pd.Series(times), errors='coerce')
                if not time_values.notna().all() or not np.isfinite(time_values).all() or not (np.diff(time_values.to_numpy()) > 0).all():
                    physical_time['reason'] = 'Mapped time values must be finite and strictly increasing.'
                else:
                    valid_mapping = True
                    for replicate, group in groups:
                        if 'Frame #' not in group.columns or group['Frame #'].tolist() != list(frames):
                            valid_mapping = False
                            break
                    if valid_mapping and 'frame' in samples:
                        time_lookup = dict(zip(frames, time_values))
                        samples['time'] = samples['frame'].map(time_lookup)
                        axis = f"Time ({metadata['unit']})"
                        physical_time = {
                            'status': 'available', 'unit': metadata['unit'], 'source': source,
                        }
                    else:
                        physical_time['reason'] = 'Sampling frame identifiers do not match saved result frames in order.'
        return {
            'status': 'available', 'samples': samples, 'x': 'time' if 'time' in samples else 'sample_index',
            'x_label': axis, 'window': window,
            'physical_time': physical_time,
            'invalid_sample_count': invalid_sample_count,
            'overall_mean': {str(row['replicate']): float(row['overall_mean']) for _, row in samples.drop_duplicates('replicate').iterrows()}
                if 'replicate' in samples else float(samples['overall_mean'].iloc[0]),
            'source': results.get('source'),
        }

    @staticmethod
    def publication_residue_payload(decomp_data: pd.DataFrame, top_n=15) -> pd.DataFrame:
        """Select favorable TOTAL residue contributions with deterministic tie order."""
        if top_n != 'all' and (isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1):
            raise ValueError("top_n must be a positive integer or 'all'")
        if decomp_data is None or decomp_data.empty or 'TOTAL Avg.' not in decomp_data.columns:
            return pd.DataFrame(columns=[] if decomp_data is None else decomp_data.columns)
        selected = decomp_data.copy()
        selected['TOTAL Avg.'] = pd.to_numeric(selected['TOTAL Avg.'], errors='coerce')
        selected = selected.dropna(subset=['TOTAL Avg.']).sort_values(
            ['TOTAL Avg.'] + (['Residue'] if 'Residue' in selected else []), kind='mergesort'
        )
        if top_n != 'all':
            selected = selected.head(top_n)
        return selected.reset_index(drop=True)
    
    def load_all_data(self) -> Dict[str, Dict]:
        """
        Load all data for all systems.
        
        Returns:
            Dictionary with structure:
            {
                'system_name': {
                    'results': {...},  # Parsed results data
                    'decomp': DataFrame  # Parsed decomposition data
                }
            }
        """
        all_data = {}
        
        for system in self.config.systems:
            logger.info(f"Processing system: {system.name}")
            
            # Find files for this system
            results_files, decomp_files = self._find_system_files(system)
            qc_metadata = self._build_qc_metadata(system, results_files, decomp_files)
            self._enforce_replicate_policy(system, qc_metadata)
            
            if not results_files:
                logger.warning(f"No results files found for system {system.name}")
                continue
            
            if not decomp_files:
                logger.warning(f"No decomposition files found for system {system.name}")
            
            # Parse and process
            results_data = self._process_results_files(results_files, system.replicates)
            if self.config.start_frame is not None:
                qc_metadata["start_frame"] = self.config.start_frame
                qc_metadata["frames_used"] = list(self._frames_used)
            decomp_data = self._process_decomp_files(decomp_files, system.replicates) if decomp_files else None
            
            all_data[system.name] = {
                'results': results_data,
                'decomp': decomp_data,
                'qc': qc_metadata,
            }
        
        return all_data

    def _build_qc_metadata(
        self,
        system: MMPBSASystemConfig,
        results_files: List[Path],
        decomp_files: List[Path],
    ) -> Dict[str, object]:
        expected = int(system.replicates)
        found_results = len(results_files)
        found_decomp = len(decomp_files)

        missing_results = max(expected - found_results, 0)
        missing_decomp = max(expected - found_decomp, 0)

        group_identity_ok, group_identity_note = self._validate_group_identity(system, results_files)
        return {
            "expected_replicates": expected,
            "found_results_replicates": found_results,
            "found_decomp_replicates": found_decomp,
            "missing_results_replicates": missing_results,
            "missing_decomp_replicates": missing_decomp,
            "replicate_complete": missing_results == 0,
            "group_identity_ok": group_identity_ok,
            "group_identity_note": group_identity_note,
            "results_files": [str(path) for path in results_files],
            "decomp_files": [str(path) for path in decomp_files],
        }

    def _enforce_replicate_policy(self, system: MMPBSASystemConfig, qc: Dict[str, object]) -> None:
        policy = self.config.incomplete_replicates_policy
        missing_results = int(qc.get("missing_results_replicates", 0) or 0)
        missing_decomp = int(qc.get("missing_decomp_replicates", 0) or 0)

        if policy == "ignore":
            return

        if missing_results > 0:
            msg = (
                f"System {system.name} missing {missing_results} result replicate(s); "
                f"found {qc.get('found_results_replicates')}/{qc.get('expected_replicates')}."
            )
            if policy == "error":
                raise ValueError(msg)
            logger.warning(msg)

        if self.config.require_decomp_replicates and missing_decomp > 0:
            msg = (
                f"System {system.name} missing {missing_decomp} decomposition replicate(s); "
                f"found {qc.get('found_decomp_replicates')}/{qc.get('expected_replicates')}."
            )
            if policy == "error":
                raise ValueError(msg)
            logger.warning(msg)

    def _validate_group_identity(self, system: MMPBSASystemConfig, files: List[Path]) -> Tuple[bool, str]:
        if system.replicates <= 1 or '{}' not in system.dir_pattern:
            return True, "Not applicable for non-patterned replicate layout."
        if not files:
            return False, "No files found to validate directory identity."

        base_path = Path(self.config.base_dir)
        expected_dirs = {
            (base_path / system.dir_pattern.format(rep)).resolve()
            for rep in range(1, system.replicates + 1)
        }
        found_dirs = {path.parent.resolve() for path in files}
        unexpected = sorted(str(path) for path in found_dirs - expected_dirs)
        if unexpected:
            return False, f"Unexpected replicate directories: {', '.join(unexpected)}"
        return True, "All replicate files matched expected directory identities."
    
    def _find_system_files(self, system: MMPBSASystemConfig) -> Tuple[List[Path], List[Path]]:
        """
        Find results and decomposition files for a system.
        
        Args:
            system: System configuration
            
        Returns:
            Tuple of (results_files, decomp_files)
        """
        base_path = Path(self.config.base_dir)
        results_files = []
        decomp_files = []
        
        # Determine file format
        file_format = self.config.file_format
        if file_format == 'auto':
            # Try to detect from first file pattern
            if '.csv' in system.results_file_pattern.lower():
                file_format = 'csv'
            elif '.dat' in system.results_file_pattern.lower():
                file_format = 'dat'
            else:
                file_format = 'dat'  # Default
        
        # Determine file extension
        ext = '.csv' if file_format == 'csv' else '.dat'
        results_pattern = self._normalize_file_pattern(system.results_file_pattern, ext)
        decomp_pattern = self._normalize_file_pattern(system.decomp_file_pattern, ext)
        
        # Find files based on replicates
        if system.replicates > 1:
            # Check if dir_pattern contains {} (separate directories) or not (same directory with numbered files)
            if '{}' in system.dir_pattern:
                # Multiple replicates: look in separate directories
                for i in range(1, system.replicates + 1):
                    dir_path = base_path / system.dir_pattern.format(i)
                    
                    # Find results file
                    results_glob = dir_path.glob(results_pattern)
                    results_files.extend([f for f in results_glob if f.is_file()])
                    
                    # Find decomp file
                    decomp_glob = dir_path.glob(decomp_pattern)
                    decomp_files.extend([f for f in decomp_glob if f.is_file()])
            else:
                # Same directory: look for numbered files (e.g., FINAL_RESULTS_MMPBSA_1.dat, _2.dat, _3.dat)
                dir_path = base_path / system.dir_pattern
                
                # Use glob pattern with wildcard to find all numbered files
                results_glob = dir_path.glob(results_pattern)
                results_files = [f for f in results_glob if f.is_file()]
                # Sort to ensure correct order
                results_files.sort()
                
                # Limit to number of replicates if more files found
                if len(results_files) > system.replicates:
                    results_files = results_files[:system.replicates]
                
                # Same for decomp files
                decomp_glob = dir_path.glob(decomp_pattern)
                decomp_files = [f for f in decomp_glob if f.is_file()]
                decomp_files.sort()
                
                if len(decomp_files) > system.replicates:
                    decomp_files = decomp_files[:system.replicates]
        else:
            # Single trajectory: look in dir_pattern without replacement
            dir_path = base_path / system.dir_pattern.replace('{}', '')
            
            # Find results file (use wildcard if present)
            results_glob = dir_path.glob(results_pattern)
            results_files.extend([f for f in results_glob if f.is_file()])
            
            # Find decomp file
            decomp_glob = dir_path.glob(decomp_pattern)
            decomp_files.extend([f for f in decomp_glob if f.is_file()])
        
        # Sort files
        results_files.sort()
        decomp_files.sort()
        
        logger.info(f"Found {len(results_files)} results files and {len(decomp_files)} decomp files for {system.name}")
        
        return results_files, decomp_files

    @staticmethod
    def _normalize_file_pattern(pattern: str, ext: str) -> str:
        """Switch CSV/DAT extensions without removing wildcards."""
        if ext == ".csv":
            return pattern.replace(".dat", ".csv")
        return pattern.replace(".csv", ".dat")
    
    def _process_results_files(self, files: List[Path], n_replicates: int) -> Optional[Dict]:
        """
        Process results files (single or replicates).
        
        Args:
            files: List of file paths
            n_replicates: Number of replicates (1 for single)
            
        Returns:
            Dictionary with processed results data
        """
        if not files:
            return None
        
        # Parse all files
        parsed_data = []
        self._frames_used = []
        for file_path in files:
            data = self.parser.parse_results(file_path, file_format='auto')
            if data:
                data['source'] = str(file_path)
                if not data.get('entropy'):
                    data['entropy'] = self._entropy_from_sibling_dat(file_path)
                if self.config.start_frame is not None:
                    data = self._apply_start_frame(data, file_path)
                parsed_data.append(data)
        
        if not parsed_data:
            logger.error("No valid results files parsed")
            return None
        
        # Process based on number of replicates
        if n_replicates == 1 or len(parsed_data) == 1:
            # Single trajectory
            return self._format_single_results(parsed_data[0])
        else:
            # Multiple replicates: average and calculate SD
            return self._process_replicate_results(parsed_data)
    
    def _apply_start_frame(self, data: Dict, file_path: Path) -> Dict:
        """Re-average frame-level results using only frames >= start_frame."""
        frames = data.get('frame_data')
        start = int(self.config.start_frame)
        if frames is None or frames.empty or 'Frame #' not in frames.columns:
            logger.warning(
                f"start_frame={start} ignored for {file_path.name}: no frame-level data "
                "(use CSV output, or set startframe in the gmx_MMPBSA input)."
            )
            self._frames_used.append(None)
            return data
        frame_numbers = pd.to_numeric(frames['Frame #'], errors='coerce')
        kept = frames[frame_numbers >= start]
        self._frames_used.append(int(len(kept)))
        if kept.empty:
            raise ValueError(f"start_frame={start} leaves no frames in {file_path}")
        numeric = [c for c in kept.columns if c != 'Frame #']
        delta = {}
        for col in numeric:
            if col in ('GGAS', 'GSOLV', 'TOTAL'):
                continue
            values = pd.to_numeric(kept[col], errors='coerce').dropna()
            std = float(values.std()) if len(values) > 1 else 0.0
            delta[col] = {
                'value': float(values.mean()),
                'std_dev': std,
                'std_err': std / np.sqrt(len(values)) if len(values) > 1 else 0.0,
            }
        binding = {
            key: float(pd.to_numeric(kept[key], errors='coerce').mean())
            for key in ('GGAS', 'GSOLV', 'TOTAL') if key in kept.columns
        }
        updated = dict(data)
        updated['delta_components'] = delta
        updated['binding_energy'] = binding
        updated['frame_data'] = kept.copy()
        updated['canonical_components'] = _canonical_component_payload(delta, kept, binding)
        updated.setdefault('metadata', {})['frames_used'] = int(len(kept))
        logger.info(f"{file_path.name}: using {len(kept)} frames with Frame # >= {start}")
        return updated

    def _entropy_from_sibling_dat(self, results_file: Path) -> Dict:
        """CSV exports carry no entropy; read IE/C2 from the FINAL_RESULTS DAT next to them."""
        if results_file.suffix.lower() == '.dat':
            return {}
        candidates = sorted(results_file.parent.glob('FINAL_RESULTS*.dat'))
        stem_match = [c for c in candidates if c.stem in results_file.stem or results_file.stem in c.stem]
        for candidate in stem_match + [c for c in candidates if c not in stem_match]:
            entropy = self.parser.extract_entropy_from_file(candidate)
            if entropy:
                logger.info(f"Entropy terms for {results_file.name} read from {candidate.name}")
                return entropy
        return {}

    def _process_replicate_results(self, parsed_data: List[Dict]) -> Dict:
        """
        Process multiple replicate results files.
        
        Calculate mean and SD across replicates.
        
        Args:
            parsed_data: List of parsed results dictionaries
            
        Returns:
            Dictionary with averaged data and SD
        """
        n_replicates = len(parsed_data)
        
        # Collect all delta component values
        delta_components = {}
        canonical_components = {}
        binding_energy = {'GGAS': [], 'GSOLV': [], 'TOTAL': []}
        entropy = {'interaction_entropy': [], 'c2_entropy': []}
        
        for data in parsed_data:
            # Delta components
            for comp_name, comp_data in data.get('delta_components', {}).items():
                if comp_name not in delta_components:
                    delta_components[comp_name] = []
                delta_components[comp_name].append(comp_data['value'])

            for comp_name, comp_data in data.get('canonical_components', {}).items():
                canonical_components.setdefault(comp_name, []).append(
                    comp_data.get('mean', comp_data.get('value'))
                )
            
            # Binding energy
            for key in ['GGAS', 'GSOLV', 'TOTAL']:
                if key in data.get('binding_energy', {}):
                    binding_energy[key].append(data['binding_energy'][key])
            
            # Entropy
            for key in ['interaction_entropy', 'c2_entropy']:
                if key in data.get('entropy', {}):
                    entropy[key].append(data['entropy'][key].get('value', 0.0))

        frame_tables = []
        for rep_idx, data in enumerate(parsed_data, start=1):
            frame_data = data.get("frame_data")
            if frame_data is not None and not frame_data.empty:
                frame_df = frame_data.copy()
                frame_df.insert(0, "Replicate", rep_idx)
                frame_tables.append(frame_df)
        
        # Calculate statistics
        # Replicates are the independent samples: report the sample SD (ddof=1)
        # across replicate averages and the corresponding standard error.
        delta_result = {
            comp_name: self._replicate_summary(values)
            for comp_name, values in delta_components.items()
            if values
        }
        canonical_result = {
            comp_name: self._replicate_summary(values)
            for comp_name, values in canonical_components.items() if values
        }
        binding_result = {
            key: self._replicate_summary(values)
            for key, values in binding_energy.items()
            if values
        }
        entropy_result = {
            key: self._replicate_summary(values)
            for key, values in entropy.items()
            if values
        }
        for key, summary in entropy_result.items():
            sigmas = [d['entropy'][key].get('sigma') for d in parsed_data
                      if key in d.get('entropy', {}) and d['entropy'][key].get('sigma') is not None]
            if sigmas:
                summary['sigma'] = float(np.mean(sigmas))
                if key == 'interaction_entropy':
                    summary['reliable'] = bool(max(sigmas) <= self.parser.IE_SIGMA_LIMIT)
        
        return {
            'delta_components': delta_result,
            'canonical_components': canonical_result,
            'binding_energy': binding_result,
            'replicate_sources': {
                str(index): data['source']
                for index, data in enumerate(parsed_data, start=1) if data.get('source')
            },
            'entropy': entropy_result,
            'n_replicates': n_replicates,
            'mode': 'replicate',
            'frame_data': pd.concat(frame_tables, ignore_index=True) if frame_tables else None,
        }
    
    def _format_single_results(self, data: Dict) -> Dict:
        """
        Format single trajectory results data.
        
        Args:
            data: Parsed results dictionary
            
        Returns:
            Formatted dictionary
        """
        # Single run: 'std' is the frame-to-frame SD and 'sem' the SEM reported by
        # gmx_MMPBSA. That SEM assumes uncorrelated frames and is optimistic;
        # independent replicates are needed for a defensible uncertainty.
        delta_components = data.get('delta_components', {})
        delta_result = {}
        for comp_name, comp_data in delta_components.items():
            delta_result[comp_name] = {
                'mean': comp_data['value'],
                'std': comp_data.get('std_dev', 0.0),
                'sem': comp_data.get('std_err', 0.0),
                'n': 1,
                'uncertainty': 'frame_sd',
            }

        # Convert binding energy format, reusing the matching Delta term's spread.
        binding_sources = {
            'GGAS': ['ΔGGAS'],
            'GSOLV': ['ΔGSOLV'],
            'TOTAL': ['ΔTOTAL', 'ΔGTOTAL'],
        }
        binding_result = {}
        for key, value in data.get('binding_energy', {}).items():
            source = next(
                (delta_components[name] for name in binding_sources.get(key, []) if name in delta_components),
                {},
            )
            binding_result[key] = {
                'mean': value,
                'std': source.get('std_dev', 0.0),
                'sem': source.get('std_err', 0.0),
                'n': 1,
                'uncertainty': 'frame_sd',
            }

        # Convert entropy format
        entropy_result = {}
        for key, entropy_data in data.get('entropy', {}).items():
            entropy_result[key] = {
                'mean': entropy_data.get('value', 0.0),
                'std': entropy_data.get('std_dev', 0.0),
                'sem': entropy_data.get('std_err', 0.0),
                'n': 1,
                'uncertainty': 'frame_sd',
                'sigma': entropy_data.get('sigma'),
                'reliable': entropy_data.get('reliable'),
            }
        
        return {
            'delta_components': delta_result,
            'binding_energy': binding_result,
            'canonical_components': data.get('canonical_components', {}),
            'source': data.get('source'),
            'entropy': entropy_result,
            'n_replicates': 1,
            'mode': 'single',
            'frame_data': data.get('frame_data'),
        }
    
    def _process_decomp_files(self, files: List[Path], n_replicates: int) -> Optional[pd.DataFrame]:
        """
        Process decomposition files (single or replicates).
        
        Args:
            files: List of file paths
            n_replicates: Number of replicates (1 for single)
            
        Returns:
            DataFrame with processed decomposition data
        """
        if not files:
            return None
        
        # Parse all files
        parsed_dfs = []
        for file_path in files:
            df = self.parser.parse_decomp(file_path, file_format='auto')
            if df is not None and not df.empty:
                parsed_dfs.append(df)
        
        if not parsed_dfs:
            logger.error("No valid decomposition files parsed")
            return None
        
        # Process based on number of replicates
        if n_replicates == 1 or len(parsed_dfs) == 1:
            # Single trajectory
            return parsed_dfs[0]
        else:
            # Multiple replicates: average and calculate SD
            return self._process_replicate_decomp(parsed_dfs)
    
    def _process_replicate_decomp(self, dfs: List[pd.DataFrame]) -> pd.DataFrame:
        """
        Process multiple replicate decomposition DataFrames.
        
        Calculate mean and SD per residue across replicates.
        
        Args:
            dfs: List of decomposition DataFrames
            
        Returns:
            Combined DataFrame with averaged values and SD
        """
        # Get all unique residues
        all_residues = set()
        for df in dfs:
            if 'Residue' in df.columns:
                all_residues.update(df['Residue'].values)
        
        # Get all value columns (exclude SD columns)
        value_columns = []
        for df in dfs:
            for col in df.columns:
                if col != 'Residue' and 'SD' not in col and col not in value_columns:
                    value_columns.append(col)
        
        # Process each residue
        combined_data = []
        for residue in sorted(all_residues):
            row = {'Residue': residue}
            
            for col in value_columns:
                values = []
                for df in dfs:
                    residue_rows = df[df['Residue'] == residue]
                    if not residue_rows.empty and col in residue_rows.columns:
                        val = residue_rows[col].iloc[0]
                        if pd.notna(val):
                            values.append(float(val))
                
                if values:
                    row[col] = np.mean(values)
                    row[f'{col} SD'] = np.std(values, ddof=1) if len(values) > 1 else 0.0
            
            combined_data.append(row)
        
        return pd.DataFrame(combined_data)

    @staticmethod
    def _replicate_summary(values: List[float]) -> Dict[str, float]:
        """Mean, sample SD (ddof=1) and SEM over replicate averages."""
        arr = np.asarray(values, dtype=float)
        n = len(arr)
        std = float(np.std(arr, ddof=1)) if n > 1 else 0.0
        return {
            'mean': float(np.mean(arr)),
            'std': std,
            'sem': std / np.sqrt(n) if n > 1 else 0.0,
            'n': n,
            'uncertainty': 'between_replicates' if n > 1 else 'none',
        }
