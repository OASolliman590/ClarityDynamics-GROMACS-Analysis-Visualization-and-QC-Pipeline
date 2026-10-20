"""
MMPBSA Parser Module
====================

Unified parser for MMPBSA CSV and DAT files.
Based on fixed_dat_file_parsing.py with enhancements.
"""

import re
import logging
from pathlib import Path
from typing import Dict, Optional, Tuple, List
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class MMPBSAParser:
    """
    Unified parser for MMPBSA files (CSV and DAT formats).
    
    Handles:
    - Results files (energy components, Delta values, entropy)
    - Decomposition files (per-residue contributions)
    - Auto-detection of file format
    - UNK residue filtering
    """
    
    def __init__(self, debug: bool = False):
        """
        Initialize parser.
        
        Args:
            debug: Enable debug logging
        """
        self.debug = debug
    
    def parse_results(self, file_path: Path, file_format: str = 'auto') -> Optional[Dict]:
        """
        Parse MMPBSA results file (CSV or DAT).
        
        Args:
            file_path: Path to results file
            file_format: 'csv', 'dat', or 'auto' for auto-detection
            
        Returns:
            Dictionary with:
            - delta_components: Dict of Delta energy components
            - binding_energy: Dict with GGAS, GSOLV, TOTAL
            - entropy: Dict with interaction_entropy, c2_entropy
            - metadata: Dict with file metadata
        """
        # Auto-detect format
        if file_format == 'auto':
            if str(file_path).lower().endswith('.csv'):
                file_format = 'csv'
            elif str(file_path).lower().endswith('.dat'):
                file_format = 'dat'
            else:
                logger.warning(f"Unknown file format for {file_path}, trying DAT format")
                file_format = 'dat'
        
        if file_format == 'dat':
            return self._parse_dat_results(file_path)
        elif file_format == 'csv':
            return self._parse_csv_results(file_path)
        else:
            logger.error(f"Unknown file format: {file_format}")
            return None
    
    def parse_decomp(self, file_path: Path, file_format: str = 'auto') -> Optional[pd.DataFrame]:
        """
        Parse MMPBSA decomposition file (CSV or DAT).
        
        Args:
            file_path: Path to decomposition file
            file_format: 'csv', 'dat', or 'auto' for auto-detection
            
        Returns:
            DataFrame with per-residue contributions (UNK filtered)
        """
        # Auto-detect format
        if file_format == 'auto':
            if str(file_path).lower().endswith('.csv'):
                file_format = 'csv'
            elif str(file_path).lower().endswith('.dat'):
                file_format = 'dat'
            else:
                logger.warning(f"Unknown file format for {file_path}, trying DAT format")
                file_format = 'dat'
        
        if file_format == 'dat':
            df = self._parse_dat_decomp(file_path)
        elif file_format == 'csv':
            df = self._parse_csv_decomp(file_path)
        else:
            logger.error(f"Unknown file format: {file_format}")
            return None
        
        # Filter UNK residues
        if df is not None:
            df = self.filter_unk_residues(df)
        
        return df
    
    def filter_unk_residues(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Filter out UNK residues from decomposition data.
        
        Args:
            data: DataFrame with Residue column
            
        Returns:
            DataFrame with UNK residues removed
        """
        if 'Residue' not in data.columns:
            return data
        
        # Filter rows where Residue contains UNK
        mask = ~data['Residue'].astype(str).str.contains('UNK', case=False, na=False)
        filtered = data[mask].copy()
        
        if self.debug:
            removed = len(data) - len(filtered)
            if removed > 0:
                logger.info(f"Filtered out {removed} UNK residues")
        
        return filtered
    
    # ==================== DAT FILE PARSING ====================
    
    def _parse_dat_results(self, file_path: Path) -> Optional[Dict]:
        """
        Parse DAT results file (based on fixed_dat_file_parsing.py).
        
        Extracts:
        - Delta components (ΔVDWAALS, ΔEEL, ΔEGB, ΔESURF, ΔGGAS, ΔGSOLV, ΔTOTAL)
        - Entropy (Interaction Entropy, C2 Entropy)
        - Metadata
        """
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            
            if self.debug:
                logger.info(f"Parsing DAT results file: {file_path}")
            
            # Extract metadata
            metadata = self._extract_metadata(content)
            
            # Extract Delta components (this is the key part)
            delta_components = self._extract_delta_components(content)
            
            # Extract binding energy from delta components
            binding_energy = self._extract_binding_energy(delta_components)
            
            # Extract entropy
            entropy = self._extract_entropy(content)
            
            result = {
                'delta_components': delta_components,
                'binding_energy': binding_energy,
                'entropy': entropy,
                'metadata': metadata
            }
            
            if self.debug:
                logger.info(f"Found {len(delta_components)} delta components")
                logger.info(f"Found {len(entropy)} entropy values")
            
            return result
            
        except Exception as e:
            logger.error(f"Error parsing DAT results file {file_path}: {e}")
            if self.debug:
                import traceback
                traceback.print_exc()
            return None
    
    def _extract_metadata(self, content: str) -> Dict:
        """Extract metadata from DAT file."""
        metadata = {}
        
        # Run date
        match = re.search(r'Run on (.+)', content)
        if match:
            metadata['run_date'] = match.group(1).strip()
        
        # Version
        match = re.search(r'gmx_MMPBSA Version=(.+)', content)
        if match:
            metadata['version'] = match.group(1).strip()
        
        # Complex structure file
        match = re.search(r'Complex Structure file:\s+(.+)', content)
        if match:
            metadata['complex_file'] = match.group(1).strip()
        
        # Number of frames
        match = re.search(r'Calculations performed using (\d+) complex frames', content)
        if match:
            metadata['num_frames'] = int(match.group(1))
        
        # Temperature
        match = re.search(r'Using temperature = (\d+\.\d+) K', content)
        if match:
            metadata['temperature'] = float(match.group(1))
        
        return metadata
    
    def _extract_delta_components(self, content: str) -> Dict:
        """
        Extract Delta components from DAT file.
        
        Looks for "Delta (Complex - Receptor - Ligand):" section.
        Extracts: ΔVDWAALS, ΔEEL, ΔEGB, ΔESURF, ΔGGAS, ΔGSOLV, ΔTOTAL
        """
        delta_components = {}
        
        # Primary method: Look for "Delta (Complex - Receptor - Ligand):" section.
        # gmx_MMPBSA >= 1.5 columns: Average, SD(Prop.), SD, SEM(Prop.), SEM.
        # Older releases: Average, Std. Dev., Std. Err. of Mean.
        delta_section_match = re.search(
            r"Delta \(Complex - Receptor - Ligand\):\s*\n([\s\S]*?)\n-+\s*\n-+\s*\n",
            content
        ) or re.search(
            r"Delta \(Complex - Receptor - Ligand\):\s*\n([\s\S]*?)(?:\n\s*\n|\Z)",
            content
        )

        if delta_section_match:
            delta_section = delta_section_match.group(1)
            components = [
                'VDWAALS', 'EEL', 'EGB', 'EPB', 'ESURF', 'ENPOLAR', 'EDISPER',
                'GGAS', 'GSOLV', 'TOTAL',
            ]
            number = r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"
            for comp in components:
                match = re.search(
                    rf"^[ \t]*Δ{comp}[ \t]*=?[ \t]+((?:{number}[ \t()]*)+)$",
                    delta_section,
                    re.MULTILINE,
                )
                if not match:
                    continue
                values = [float(v) for v in re.findall(number, match.group(1))]
                if not values:
                    continue
                if len(values) >= 5:
                    std_dev, std_err = values[2], values[4]
                elif len(values) >= 3:
                    std_dev, std_err = values[1], values[2]
                elif len(values) == 2:
                    std_dev, std_err = values[1], 0.0
                else:
                    std_dev, std_err = 0.0, 0.0
                delta_components[f'Δ{comp}'] = {
                    'value': values[0],
                    'std_dev': std_dev,
                    'std_err': std_err,
                }

        # Alternative method: Look for FINAL RESULTS section
        if not delta_components:
            final_results_match = re.search(r'FINAL RESULTS:\s+\n([\s\S]+?)(?:\n\n|\Z)', content)
            if final_results_match:
                final_results = final_results_match.group(1)
                
                components = ['VDWAALS', 'EEL', 'EGB', 'ESURF', 'GGAS', 'GSOLV', 'TOTAL']
                for comp in components:
                    pattern = rf"Δ{comp}\s+=\s+([\d.-]+)"
                    match = re.search(pattern, final_results)
                    if match:
                        delta_components[f'Δ{comp}'] = {
                            'value': float(match.group(1)),
                            'std_dev': 0.0,
                            'std_err': 0.0
                        }
        
        # Also try format with equals and parentheses
        if not delta_components:
            pattern = r'(ΔGGAS|ΔGSOLV|ΔGTOTAL|ΔVDWAALS|ΔEEL|ΔESURF|ΔEGB)\s+=\s+([\-\d\.]+)\s+\(([\-\d\.]+)\)\s+\(([\-\d\.]+)\)'
            matches = re.findall(pattern, content)
            
            for match in matches:
                component = match[0]
                value = float(match[1])
                std_dev = float(match[2])
                std_err = float(match[3])
                
                delta_components[component] = {
                    'value': value,
                    'std_dev': std_dev,
                    'std_err': std_err
                }
        
        return delta_components
    
    def _extract_binding_energy(self, delta_components: Dict) -> Dict:
        """Extract binding energy components from delta components."""
        binding_energy = {}
        
        # Direct mapping
        if 'ΔGGAS' in delta_components:
            binding_energy['GGAS'] = delta_components['ΔGGAS']['value']
        if 'ΔGSOLV' in delta_components:
            binding_energy['GSOLV'] = delta_components['ΔGSOLV']['value']
        if 'ΔGTOTAL' in delta_components:
            binding_energy['TOTAL'] = delta_components['ΔGTOTAL']['value']
        elif 'ΔTOTAL' in delta_components:
            binding_energy['TOTAL'] = delta_components['ΔTOTAL']['value']
        
        # Calculate if not found
        if not binding_energy and 'ΔVDWAALS' in delta_components and 'ΔEEL' in delta_components:
            vdw = delta_components['ΔVDWAALS']['value']
            eel = delta_components['ΔEEL']['value']
            binding_energy['GGAS'] = vdw + eel
            
            if 'ΔESURF' in delta_components and 'ΔEGB' in delta_components:
                esurf = delta_components['ΔESURF']['value']
                egb = delta_components['ΔEGB']['value']
                binding_energy['GSOLV'] = esurf + egb
                binding_energy['TOTAL'] = binding_energy['GGAS'] + binding_energy['GSOLV']
        
        return binding_energy
    
    # σ(ΔE_int) above which interaction entropy is unreliable (Duan et al. 2016;
    # Sun et al. 2018 recommend ≲ 3.6 kcal/mol).
    IE_SIGMA_LIMIT = 3.6

    def _extract_entropy(self, content: str) -> Dict:
        """Extract interaction-entropy (IE) and C2 entropy results.

        Values are −TΔS in kcal/mol, as printed by gmx_MMPBSA. Handles the
        tabular layout of gmx_MMPBSA >= 1.5::

            ENTROPY RESULTS (INTERACTION ENTROPY):
            Model     σ(Int. Energy)    Average    SD    SEM
            GB             2.13           2.18    0.06   0.01

        and the older single-line summaries.
        """
        entropy = {}
        number = r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"

        for key, header in (('interaction_entropy', 'INTERACTION ENTROPY'), ('c2_entropy', 'C2 ENTROPY')):
            section = re.search(
                rf"ENTROPY RESULTS \({header}\):\s*\n([\s\S]*?)(?:\n\s*\n|\Z)", content
            )
            if not section:
                continue
            for line in section.group(1).splitlines():
                line = line.strip()
                if not line or line.startswith(('-', 'Model', 'Energy')):
                    continue
                parts = line.split()
                model = parts[0] if parts and not re.fullmatch(number, parts[0]) else None
                values = re.findall(number, line[len(model):] if model else line)
                # C2 confidence interval "2.05-2.53" splits into two numbers; keep raw text too.
                ci_match = re.search(rf"({number})\s*-\s*({number})\s*$", line) if key == 'c2_entropy' else None
                if len(values) < 2:
                    continue
                entry = {
                    'model': model,
                    'sigma': float(values[0]),
                    'value': float(values[1]),
                    'std_dev': float(values[2]) if len(values) > 2 else 0.0,
                }
                if key == 'interaction_entropy':
                    entry['std_err'] = float(values[3]) if len(values) > 3 else 0.0
                elif ci_match:
                    entry['confidence_interval'] = f"{ci_match.group(1)}-{ci_match.group(2)}"
                # The σ criterion is established for IE; C2 is reported without a verdict.
                entry['reliable'] = entry['sigma'] <= self.IE_SIGMA_LIMIT if key == 'interaction_entropy' else None
                entropy[key] = entry
                break

        # Older gmx_MMPBSA summaries: "... INTERACTION ENTROPY): ... gb  <sigma>  <value>"
        for key, header in (('interaction_entropy', 'INTERACTION ENTROPY'), ('c2_entropy', 'C2 ENTROPY')):
            if key in entropy:
                continue
            match = re.search(rf"{header}\):[\s\S]*?\b(?:gb|pb|GB|PB)\s+({number})\s+({number})", content)
            if match:
                sigma = float(match.group(1))
                entropy[key] = {
                    'model': None,
                    'sigma': sigma,
                    'value': float(match.group(2)),
                    'std_dev': 0.0,
                    'std_err': 0.0,
                    'reliable': (sigma <= self.IE_SIGMA_LIMIT) if key == 'interaction_entropy' else None,
                }

        return entropy

    def extract_entropy_from_file(self, file_path: Path) -> Dict:
        """Entropy terms from a gmx_MMPBSA DAT file (used alongside CSV results)."""
        try:
            return self._extract_entropy(Path(file_path).read_text(errors="replace"))
        except OSError as exc:
            logger.warning(f"Could not read entropy from {file_path}: {exc}")
            return {}

    def _parse_dat_decomp(self, file_path: Path) -> Optional[pd.DataFrame]:
        """
        Parse DAT decomposition file (based on fixed_dat_file_parsing.py).
        
        Extracts per-residue contributions, filters UNK.
        """
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            if self.debug:
                logger.info(f"Parsing DAT decomposition file: {file_path}")
            
            # Find the Total Energy Decomposition section
            decomp_section_start = -1
            for i, line in enumerate(lines):
                if "Total Energy Decomposition:" in line:
                    decomp_section_start = i
                    break
            
            if decomp_section_start == -1:
                logger.warning(f"Could not find Total Energy Decomposition section in {file_path}")
                return None
            
            # Find header line
            header_line = None
            for i in range(decomp_section_start, min(decomp_section_start + 10, len(lines))):
                if "Residue," in lines[i]:
                    header_line = i
                    break
            
            if header_line is None:
                logger.warning(f"Could not find header line in decomposition section of {file_path}")
                return None
            
            # Skip header and subheader
            data_start = header_line + 2
            processed_data = []
            
            # Process each residue line
            for i in range(data_start, len(lines)):
                line = lines[i].strip()
                if not line or line.startswith('--') or line.startswith('Complex:'):
                    continue
                
                parts = line.split(',')
                if len(parts) < 18:
                    continue
                
                residue = parts[0].strip()
                
                # Skip UNK residues
                if 'UNK' in residue.upper():
                    continue
                
                try:
                    row_data = {
                        'Residue': residue,
                        'van der Waals Avg.': float(parts[4]),
                        'van der Waals Avg. SD': float(parts[5]),
                        'Electrostatic Avg.': float(parts[7]),
                        'Electrostatic Avg. SD': float(parts[8]),
                        'Polar Solvation Avg.': float(parts[10]),
                        'Polar Solvation Avg. SD': float(parts[11]),
                        'Non-Polar Solv. Avg.': float(parts[13]),
                        'Non-Polar Solv. Avg. SD': float(parts[14]),
                        'TOTAL Avg.': float(parts[16]),
                        'TOTAL Avg. SD': float(parts[17])
                    }
                    processed_data.append(row_data)
                except (IndexError, ValueError) as e:
                    if self.debug:
                        logger.debug(f"Error processing line {i}: {e}")
                    continue
            
            df = pd.DataFrame(processed_data)
            
            if self.debug:
                logger.info(f"Parsed {len(df)} residues from decomposition file")
            
            return df
            
        except Exception as e:
            logger.error(f"Error parsing DAT decomposition file {file_path}: {e}")
            if self.debug:
                import traceback
                traceback.print_exc()
            return None
    
    # ==================== CSV FILE PARSING ====================
    
    def _parse_csv_results(self, file_path: Path, chunksize: int = 10000) -> Optional[Dict]:
        """
        Parse CSV results file.
        
        Extracts Delta components and binding energy.
        Supports chunked reading for large files.
        """
        try:
            if self.debug:
                logger.info(f"Parsing CSV results file: {file_path}")
            
            # Check file size to decide on chunking
            file_size = file_path.stat().st_size
            use_chunks = file_size > 10 * 1024 * 1024  # 10 MB
            
            if use_chunks:
                return self._parse_csv_results_chunked(file_path, chunksize)
            
            # Read entire file
            df = pd.read_csv(file_path, dtype=object, low_memory=False)
            
            if 'Frame #' not in df.columns:
                logger.error(f"CSV file {file_path} does not contain 'Frame #' column")
                return None
            
            # Convert numeric columns
            numeric_cols = [col for col in df.columns if col != 'Frame #']
            for col in numeric_cols:
                df[col] = pd.to_numeric(df[col], errors='coerce')
            
            # Extract Delta components (all columns except Frame #, GGAS, GSOLV, TOTAL)
            delta_components = {}
            for col in numeric_cols:
                if col not in ['GGAS', 'GSOLV', 'TOTAL']:
                    mean_val = df[col].mean()
                    std_val = df[col].std() if len(df) > 1 else 0.0
                    delta_components[col] = {
                        'value': mean_val,
                        'std_dev': std_val,
                        'std_err': std_val / np.sqrt(len(df)) if len(df) > 1 else 0.0
                    }
            
            # Extract binding energy
            binding_energy = {}
            if 'GGAS' in df.columns:
                binding_energy['GGAS'] = df['GGAS'].mean()
            if 'GSOLV' in df.columns:
                binding_energy['GSOLV'] = df['GSOLV'].mean()
            if 'TOTAL' in df.columns:
                binding_energy['TOTAL'] = df['TOTAL'].mean()
            
            return {
                'delta_components': delta_components,
                'binding_energy': binding_energy,
                'entropy': {},  # CSV typically doesn't have entropy
                'metadata': {'num_frames': len(df)},
                'frame_data': df.copy(),
            }
            
        except Exception as e:
            logger.error(f"Error parsing CSV results file {file_path}: {e}")
            if self.debug:
                import traceback
                traceback.print_exc()
            return None
    
    def _parse_csv_results_chunked(self, file_path: Path, chunksize: int) -> Optional[Dict]:
        """Parse large CSV results file using chunking."""
        try:
            chunks = pd.read_csv(file_path, dtype=object, chunksize=chunksize)
            
            delta_components = {}
            binding_energy = {'GGAS': [], 'GSOLV': [], 'TOTAL': []}
            frame_count = 0
            frame_chunks = []
            
            for chunk in chunks:
                numeric_cols = [col for col in chunk.columns if col != 'Frame #']
                for col in numeric_cols:
                    chunk[col] = pd.to_numeric(chunk[col], errors='coerce')
                frame_chunks.append(chunk.copy())
                
                # Accumulate values
                for col in numeric_cols:
                    if col not in ['GGAS', 'GSOLV', 'TOTAL']:
                        if col not in delta_components:
                            delta_components[col] = []
                        delta_components[col].extend(chunk[col].dropna().tolist())
                
                if 'GGAS' in chunk.columns:
                    binding_energy['GGAS'].extend(chunk['GGAS'].dropna().tolist())
                if 'GSOLV' in chunk.columns:
                    binding_energy['GSOLV'].extend(chunk['GSOLV'].dropna().tolist())
                if 'TOTAL' in chunk.columns:
                    binding_energy['TOTAL'].extend(chunk['TOTAL'].dropna().tolist())
                
                frame_count += len(chunk)
            
            # Calculate means and SDs
            delta_result = {}
            for col, values in delta_components.items():
                if values:
                    delta_result[col] = {
                        'value': np.mean(values),
                        'std_dev': np.std(values) if len(values) > 1 else 0.0,
                        'std_err': np.std(values) / np.sqrt(len(values)) if len(values) > 1 else 0.0
                    }
            
            binding_result = {}
            for key, values in binding_energy.items():
                if values:
                    binding_result[key] = np.mean(values)
            
            return {
                'delta_components': delta_result,
                'binding_energy': binding_result,
                'entropy': {},
                'metadata': {'num_frames': frame_count},
                'frame_data': pd.concat(frame_chunks, ignore_index=True) if frame_chunks else None,
            }
            
        except Exception as e:
            logger.error(f"Error parsing chunked CSV results file {file_path}: {e}")
            return None
    
    def _parse_csv_decomp(self, file_path: Path, chunksize: int = 10000) -> Optional[pd.DataFrame]:
        """
        Parse CSV decomposition file.
        
        Filters UNK residues and calculates per-residue averages.
        Supports chunked reading for large files.
        """
        try:
            if self.debug:
                logger.info(f"Parsing CSV decomposition file: {file_path}")
            
            # Check file size
            file_size = file_path.stat().st_size
            use_chunks = file_size > 10 * 1024 * 1024  # 10 MB
            
            if use_chunks:
                return self._parse_csv_decomp_chunked(file_path, chunksize)
            
            # Read entire file
            df = pd.read_csv(file_path, dtype=object, low_memory=False)
            
            if 'Residue' not in df.columns:
                logger.error(f"CSV file {file_path} does not contain 'Residue' column")
                return None
            
            # Filter UNK
            df = df[~df['Residue'].astype(str).str.contains('UNK', case=False, na=False)]
            
            # Convert numeric columns
            energy_cols = [col for col in df.columns if col not in ['Frame #', 'Residue']]
            for col in energy_cols:
                df[col] = pd.to_numeric(df[col], errors='coerce')
            
            # If Frame # exists, calculate averages per residue
            if 'Frame #' in df.columns:
                grouped = df.groupby('Residue')
                avg_data = []
                
                for residue, group in grouped:
                    row = {'Residue': residue}
                    for col in energy_cols:
                        row[f'{col} Avg.'] = group[col].mean()
                        row[f'{col} Avg. SD'] = group[col].std() if len(group) > 1 else 0.0
                    avg_data.append(row)
                
                df = pd.DataFrame(avg_data)
            
            return df
            
        except Exception as e:
            logger.error(f"Error parsing CSV decomposition file {file_path}: {e}")
            if self.debug:
                import traceback
                traceback.print_exc()
            return None
    
    def _parse_csv_decomp_chunked(self, file_path: Path, chunksize: int) -> Optional[pd.DataFrame]:
        """Parse large CSV decomposition file using chunking."""
        try:
            chunks = pd.read_csv(file_path, dtype=object, chunksize=chunksize)
            
            residue_data = {}
            
            for chunk in chunks:
                if 'Residue' not in chunk.columns:
                    continue
                
                # Filter UNK
                chunk = chunk[~chunk['Residue'].astype(str).str.contains('UNK', case=False, na=False)]
                
                energy_cols = [col for col in chunk.columns if col not in ['Frame #', 'Residue']]
                for col in energy_cols:
                    chunk[col] = pd.to_numeric(chunk[col], errors='coerce')
                
                # Group by residue and accumulate
                if 'Frame #' in chunk.columns:
                    grouped = chunk.groupby('Residue')
                    for residue, group in grouped:
                        if residue not in residue_data:
                            residue_data[residue] = {col: [] for col in energy_cols}
                        
                        for col in energy_cols:
                            residue_data[residue][col].extend(group[col].dropna().tolist())
            
            # Calculate averages and SDs
            avg_data = []
            for residue, data in residue_data.items():
                row = {'Residue': residue}
                for col, values in data.items():
                    if values:
                        row[f'{col} Avg.'] = np.mean(values)
                        row[f'{col} Avg. SD'] = np.std(values) if len(values) > 1 else 0.0
                avg_data.append(row)
            
            return pd.DataFrame(avg_data)
            
        except Exception as e:
            logger.error(f"Error parsing chunked CSV decomposition file {file_path}: {e}")
            return None

