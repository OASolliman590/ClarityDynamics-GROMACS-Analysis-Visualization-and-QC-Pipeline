"""
PCA Parser Module
=================

Parse XVG, PDB, and XPM files for PCA analysis.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np

logger = logging.getLogger(__name__)


class PCAParser:
    """
    Parse PCA-related files (XVG, PDB, XPM).
    
    Handles:
    - XVG files (eigenvalues, cosine content, projections)
    - PDB files (3D projections)
    - XPM files (RMSIP heatmaps)
    """
    
    def __init__(self, debug: bool = False):
        """
        Initialize parser.
        
        Args:
            debug: Enable debug logging
        """
        self.debug = debug
    
    def parse_xvg_single(self, file_path: Path) -> Tuple[Optional[str], List[float], List[float], Optional[str], Optional[str]]:
        """
        Parse single dataset from XVG file.
        
        Args:
            file_path: Path to XVG file
            
        Returns:
            Tuple of (title, x_values, y_values, x_title, y_title)
        """
        try:
            with open(file_path, 'r') as f:
                data = f.read()
            
            return self._parse_xvg_content(data)
            
        except FileNotFoundError:
            logger.error(f"File not found: {file_path}")
            return None, [], [], None, None
        except Exception as e:
            logger.error(f"Error reading {file_path}: {e}")
            return None, [], [], None, None
    
    def parse_xvg_multiple(self, file_path: Path) -> Tuple[List[str], List[List[float]], List[List[float]], List[Optional[str]], List[Optional[str]]]:
        """
        Parse multiple datasets from XVG file.
        
        Handles:
        - Files with '&' separator
        - Multi-column files
        - Individual PC files (eigenvec_1.xvg, etc.)
        
        Args:
            file_path: Path to XVG file
            
        Returns:
            Tuple of (titles, x_axes_list, y_axes_list, x_titles, y_titles)
        """
        try:
            if not file_path.exists():
                logger.warning(f"File not found: {file_path}")
                return [], [], [], [], []
            
            # Try individual PC files first
            base_dir = file_path.parent
            individual_files = []
            for i in range(1, 11):  # Check PC 1-10
                pc_file = base_dir / f"eigenvec_{i}.xvg"
                if pc_file.exists():
                    individual_files.append(pc_file)
            
            if individual_files:
                logger.info(f"Found {len(individual_files)} individual PC files")
                titles = []
                x_axes = []
                y_axes = []
                x_titles = []
                y_titles = []
                
                for pc_file in sorted(individual_files):
                    title, x, y, x_t, y_t = self.parse_xvg_single(pc_file)
                    if x and y:
                        pc_num = pc_file.stem.split('_')[-1]
                        titles.append(f"PC {pc_num}")
                        x_axes.append(x)
                        y_axes.append(y)
                        x_titles.append(x_t or "Frame")
                        y_titles.append(y_t or f"PC {pc_num}")
                
                return titles, x_axes, y_axes, x_titles, y_titles
            
            # Otherwise parse combined file
            with open(file_path, 'r') as f:
                content = f.read()
            
            # Check for '&' separator
            if '&' in content:
                sets = content.split('&')
                titles = []
                x_axes = []
                y_axes = []
                x_titles = []
                y_titles = []
                
                for data_set in sets:
                    if len(data_set.strip()) > 10:
                        title, x, y, x_t, y_t = self._parse_xvg_content(data_set)
                        if x and y:
                            titles.append(title or "PC")
                            x_axes.append(x)
                            y_axes.append(y)
                            x_titles.append(x_t or "Frame")
                            y_titles.append(y_t or "Value")
                
                return titles, x_axes, y_axes, x_titles, y_titles
            
            # Try multi-column format
            return self._parse_xvg_multicolumn(content)
            
        except Exception as e:
            logger.error(f"Error parsing multiple datasets from {file_path}: {e}")
            return [], [], [], [], []
    
    def _parse_xvg_content(self, data: str) -> Tuple[Optional[str], List[float], List[float], Optional[str], Optional[str]]:
        """Parse XVG content and extract data."""
        x_axis = []
        y_axis = []
        x_title = None
        y_title = None
        title = None
        
        for line in data.splitlines():
            # Parse title
            if '@' in line and 'title' in line:
                try:
                    title = line.split('"')[1]
                except:
                    title = "PCA Analysis"
                continue
            
            # Parse x-axis label
            if '@' in line and 'xaxis' in line and 'label' in line:
                try:
                    x_title = line.split('"')[1]
                except:
                    x_title = "PC"
                continue
            
            # Parse y-axis label
            if '@' in line and 'yaxis' in line and 'label' in line:
                try:
                    y_title = line.split('"')[1]
                except:
                    y_title = "Value"
                continue
            
            # Skip comment and metadata lines
            if line.startswith('#') or line.startswith('@'):
                continue
            
            # Parse data line
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    x_axis.append(float(parts[0]))
                    y_axis.append(float(parts[1]))
                except:
                    continue
        
        return title, x_axis, y_axis, x_title, y_title
    
    def _parse_xvg_multicolumn(self, content: str) -> Tuple[List[str], List[List[float]], List[List[float]], List[Optional[str]], List[Optional[str]]]:
        """Parse multi-column XVG file."""
        cols = []
        data_title = None
        x_data = []
        
        for line in content.splitlines():
            if '@' in line and 'title' in line:
                try:
                    data_title = line.split('"')[1]
                except:
                    data_title = "PCA Analysis"
                continue
            
            if line.startswith('#') or line.startswith('@'):
                continue
            
            parts = line.strip().split()
            if len(parts) > 1:
                try:
                    if not x_data:
                        cols = [[] for _ in range(len(parts) - 1)]
                    
                    x_val = float(parts[0])
                    x_data.append(x_val)
                    
                    for i in range(1, len(parts)):
                        if i-1 < len(cols):
                            cols[i-1].append(float(parts[i]))
                except ValueError:
                    continue
        
        if x_data and cols:
            titles = []
            x_axes = []
            y_axes = []
            x_titles = []
            y_titles = []
            
            for i, col in enumerate(cols):
                if len(col) == len(x_data):
                    titles.append(f"PC {i+1}")
                    x_axes.append(x_data)
                    y_axes.append(col)
                    x_titles.append("Frame")
                    y_titles.append(f"PC {i+1}")
            
            return titles, x_axes, y_axes, x_titles, y_titles
        
        return [], [], [], [], []

    def parse_xpm(self, file_path: Path) -> Tuple[List[List[float]], List[float], List[float]]:
        """
        Parse an XPM file into matrix values and axis coordinates.

        Returns:
            matrix_vals (list of rows), x_coords, y_coords
        """
        letter_val: Dict[str, float] = {}
        matrix_vals: List[List[float]] = []
        x_coords: List[float] = []
        y_coords: List[float] = []

        if not file_path.exists():
            logger.warning(f"XPM file not found: {file_path}")
            return matrix_vals, x_coords, y_coords

        try:
            with open(file_path, "r") as f:
                enter = False
                matrix = False
                for line in f:
                    if line.startswith('"A '):
                        enter = True
                    if 'x-axis' in line:
                        enter = False
                        try:
                            x_coords = list(map(float, line.split()[2:-1]))
                        except Exception:
                            x_coords = []
                    if enter:
                        start = line.find('/*')
                        end = line.find('*/')
                        if start != -1 and end != -1:
                            val = line[start + 2:end].strip()[1:-1]
                            letter = line[1]
                            try:
                                letter_val[letter] = float(val)
                            except ValueError:
                                continue
                    if 'y-axis' in line:
                        matrix = True
                        try:
                            y_coords = list(map(float, line.split()[2:-1]))
                        except Exception:
                            y_coords = []
                        continue
                    if matrix:
                        if ',' in line:
                            row = line.split(',')[0][1:-1]
                        else:
                            row = line[1:-2]
                        row_val = []
                        for letter in row:
                            if letter in letter_val:
                                row_val.append(letter_val[letter])
                        if row_val:
                            matrix_vals.append(row_val)

            return matrix_vals, x_coords, y_coords
        except Exception as exc:
            logger.error(f"Error parsing XPM file {file_path}: {exc}")
            return [], [], []
    
    def parse_pdb_coords(self, file_path: Path) -> Tuple[List[float], List[float], List[float]]:
        """
        Extract x, y, z coordinates from PDB file.
        
        Args:
            file_path: Path to PDB file
            
        Returns:
            Tuple of (x_coords, y_coords, z_coords)
        """
        try:
            x = []
            y = []
            z = []
            
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        try:
                            coords = line[30:56]
                            x.append(float(coords[:8].strip()))
                            y.append(float(coords[8:16].strip()))
                            z.append(float(coords[16:].strip()))
                        except:
                            continue
            
            return x, y, z
            
        except FileNotFoundError:
            logger.error(f"File not found: {file_path}")
            return [], [], []
        except Exception as e:
            logger.error(f"Error reading PDB file {file_path}: {e}")
            return [], [], []
    

    def parse_publication_manifest(self, file_path: Path) -> Dict:
        """Read explicit imported projection/basis and diagnostic contracts.

        Each state names a projections CSV (PC1, PC2, optional PC3), an
        eigenvalues XVG, basis_id, basis_verified, basis_source and unit.
        Diagnostic CSV values require explicit compatible=True and provenance;
        window/split-half values are never inferred from filenames.
        """
        import json
        import pandas as pd
        from .publication import diagnostic, DIAGNOSTICS

        states, records = [], []
        try:
            manifest = json.loads(Path(file_path).read_text())
            names = set()
            for item in manifest.get('states', []):
                name = str(item['name'])
                if name in names:
                    raise ValueError(f'Duplicate PCA state name: {name}')
                names.add(name)
                source = Path(file_path).parent / item['projections']
                frame = pd.read_csv(source)
                columns = [c for c in ('PC1', 'PC2', 'PC3') if c in frame]
                if columns[:2] != ['PC1', 'PC2'] or frame.empty:
                    raise ValueError(f'{name}: PC1 and PC2 samples are required')
                values = frame[columns].to_numpy(dtype=float)
                if not np.isfinite(values).all():
                    raise ValueError(f'{name}: projection samples must be finite')
                spectrum_path = Path(file_path).parent / item['eigenvalues'] if item.get('eigenvalues') else None
                eigenvalues = self.parse_xvg_single(spectrum_path)[2] if spectrum_path else []
                states.append(dict(name=name, projections={c: frame[c].astype(float).tolist() for c in columns},
                    eigenvalues=eigenvalues, basis_id=item.get('basis_id'),
                    basis_verified=item.get('basis_verified') is True, basis_source=item.get('basis_source', ''),
                    source=str(source), eigenvalue_source=str(spectrum_path or ''), unit=item.get('unit', 'source coordinate units')))
            if manifest.get('diagnostics'):
                source = Path(file_path).parent / manifest['diagnostics']
                try:
                    frame = pd.read_csv(source).fillna('')
                except (OSError, ValueError) as exc:
                    records.append(diagnostic('diagnostic_inputs', reason=f'Optional diagnostic table unavailable: {exc}', source=str(source)))
                    frame = pd.DataFrame()
                for row in frame.to_dict('records'):
                    metric = row.get('metric')
                    if metric not in DIAGNOSTICS:
                        continue
                    value = pd.to_numeric(row.get('value'), errors='coerce')
                    valid = (str(row.get('compatible')).lower() == 'true' and bool(row.get('source'))
                             and np.isfinite(value) and 0 <= value <= 1)
                    if metric == 'window_sensitivity' and not row.get('window'):
                        valid = False
                    records.append(diagnostic(metric, float(value) if valid else None,
                        reason='' if valid else 'Requires a finite fraction, compatible subspace, provenance and dedicated window when applicable',
                        state_pair=row.get('state_pair', ''), window=row.get('window', ''), source=row.get('source', str(source))))
        except (OSError, ValueError, KeyError, TypeError) as exc:
            records.append(diagnostic('publication_inputs', reason=f'Invalid publication manifest: {exc}', source=str(file_path)))
            # An invalid manifest cannot authorize a partial multi-state overlay.
            states = []
        return {'states': states, 'diagnostics': records}
