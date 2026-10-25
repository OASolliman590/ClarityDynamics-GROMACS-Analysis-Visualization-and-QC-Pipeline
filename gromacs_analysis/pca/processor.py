"""
PCA Processor Module
====================

Auto-detect and process PCA files.
Handles directory structure:
- base_dir: HOLO files + shared files (eigenvals, proj.xvg)
- base_dir/PCA_Vis/ccall/: APO cosine content
- base_dir/PCA_Vis/: APO 2D/3D projections + HOLO 2D/3D projections
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re

from .config import PCAConfig, PCASystemConfig
from .parser import PCAParser

logger = logging.getLogger(__name__)


class PCAProcessor:
    """
    Auto-detect and process PCA files.
    
    Automatically finds:
    - Eigenvalues files
    - Cosine content files (APO and HOLO)
    - Projection files (1D, 2D, 3D)
    - Ligand names from filenames
    """
    
    def __init__(self, config: PCAConfig):
        """
        Initialize processor.
        
        Args:
            config: PCAConfig object
        """
        self.config = config
        self.parser = PCAParser()
        logger.info("PCAProcessor initialized")
    
    def auto_detect_files(self) -> Dict[str, Dict]:
        """
        Auto-detect all PCA files based on directory structure.
        
        Uses explicit ligand_name from config if provided, otherwise auto-detects.
        Uses apo_base_dir if specified, otherwise uses base_dir for APO files.
        
        Returns:
            Dict with detected file paths:
            {
                'eigenvals': Path,
                'proj_1d': Path,
                'cc_holo': Path,
                'cc_apo': Path,
                'proj_2d_apo': {12: Path, 13: Path, 23: Path},
                'proj_2d_holo': {12: Path, 13: Path, 23: Path},
                'proj_3d_apo': Path,
                'proj_3d_holo': Path,
                'ligand_name': str
            }
        """
        base_dir = Path(self.config.base_dir)
        apo_dir = Path(self.config.apo_base_dir) if self.config.apo_base_dir else base_dir
        vis_dir = base_dir / 'PCA_Vis'
        apo_vis_dir = apo_dir / 'PCA_Vis' if self.config.apo_base_dir else vis_dir
        # Check for ccall in base_dir first (when files are in same directory)
        # Then check in PCA_Vis subdirectory as fallback
        ccall_dir = apo_dir / 'ccall'  # Direct ccall folder in base directory
        ccall_dir_vis = apo_vis_dir / 'ccall'  # ccall in PCA_Vis subdirectory
        
        detected = {
            'eigenvals': None,
            'proj_1d': None,
            'cc_holo': None,
            'cc_apo': None,
            'proj_2d_apo': {},
            'proj_2d_holo': {},
            'proj_3d_apo': None,
            'proj_3d_holo': None,
            'ligand_name': None,
            'proj_eig': {},
            'fel': {}
        }
        
        # Find eigenvalues (in base_dir)
        eigenvals_patterns = [
            'eigenvals_all_ref.xvg',
            'eigenvals*.xvg'
        ]
        for pattern in eigenvals_patterns:
            matches = list(base_dir.glob(pattern))
            if matches:
                detected['eigenvals'] = sorted(matches)[0]
                logger.info(f"Found eigenvalues: {detected['eigenvals']}")
                break
        
        # Find 1D projection (in base_dir)
        proj_patterns = ['proj.xvg']
        for pattern in proj_patterns:
            matches = list(base_dir.glob(pattern))
            if matches:
                detected['proj_1d'] = matches[0]
                logger.info(f"Found 1D projection: {detected['proj_1d']}")
                break
        
        # Find HOLO cosine content (in base_dir)
        cc_holo_patterns = ['cc_all.xvg']
        for pattern in cc_holo_patterns:
            matches = list(base_dir.glob(pattern))
            if matches:
                detected['cc_holo'] = matches[0]
                logger.info(f"Found HOLO cosine content: {detected['cc_holo']}")
                break
        
        # Find APO cosine content in the common PCA directory layouts.
        # Priority: 1) base_dir/ccall, 2) base_dir/PCA_Vis/ccall, 3) apo_base_dir/ccall
        # IMPORTANT: Do NOT check apo_dir directly to avoid picking up HOLO file
        ccall_dirs = []
        # First check: direct ccall folder in base directory (most common)
        if ccall_dir.exists():
            ccall_dirs.append(ccall_dir)
        # Second check: ccall in PCA_Vis subdirectory
        if ccall_dir_vis.exists() and ccall_dir_vis != ccall_dir:
            ccall_dirs.append(ccall_dir_vis)
        # Third check: if apo_base_dir is different, check its ccall folder
        if self.config.apo_base_dir and apo_dir != base_dir:
            apo_ccall = apo_dir / 'ccall'
            if apo_ccall.exists() and apo_ccall not in ccall_dirs:
                ccall_dirs.append(apo_ccall)
            apo_ccall_vis = apo_dir / 'PCA_Vis' / 'ccall'
            if apo_ccall_vis.exists() and apo_ccall_vis not in ccall_dirs:
                ccall_dirs.append(apo_ccall_vis)
        
        cc_apo_patterns = ['cc_all.xvg']
        for ccall_dir_to_check in ccall_dirs:
            for pattern in cc_apo_patterns:
                matches = list(ccall_dir_to_check.glob(pattern))
                if matches:
                    detected['cc_apo'] = matches[0]
                    logger.info(f"Found APO cosine content: {detected['cc_apo']}")
                    break
            if detected['cc_apo']:
                break
        
        # Use explicit ligand name from config if provided, otherwise auto-detect
        ligand_name = self.config.ligand_name
        if not ligand_name and vis_dir.exists():
            ligand_name = self._extract_ligand_name(vis_dir)
        if ligand_name:
            detected['ligand_name'] = ligand_name
            logger.info(f"Using ligand name: {ligand_name} {'(from config)' if self.config.ligand_name else '(auto-detected)'}")
        
        # Find 2D and 3D projections
        # Some workflows write projections in base_dir, while others write them
        # under PCA_Vis. Prefer base_dir and keep PCA_Vis as a fallback.
        
        # Determine search directories
        # For APO: use apo_base_dir if specified, otherwise base_dir
        apo_search_base = apo_dir if self.config.apo_base_dir else base_dir
        # For HOLO: use base_dir
        holo_search_base = base_dir
        
        # Try base_dir first, then PCA_Vis as fallback.
        search_dirs_apo = [apo_search_base, apo_search_base / 'PCA_Vis'] if apo_search_base.exists() else []
        search_dirs_holo = [holo_search_base, vis_dir] if holo_search_base.exists() else []
        
        # Find 2D projections
        # Handle both naming conventions:
        # 1. With suffixes: 2dproj_12_Apo.xvg / 2dproj_12_{ligand}.xvg
        # 2. Without suffixes: 2dproj_12.xvg (distinguished by directory)
        for pair in [(1, 2), (1, 3), (2, 3)]:
            pc1, pc2 = pair
            # APO 2D projections - try patterns with suffix first, then without
            patterns_apo = [
                f'2dproj_{pc1}{pc2}_Apo.xvg',
                f'2dproj_{pc1}{pc2}_APO.xvg',
                f'2dproj_{pc1}{pc2}_apo.xvg',
                f'2dproj_{pc1}{pc2}.xvg'  # Fallback: no suffix (distinguished by directory)
            ]
            for search_dir in search_dirs_apo:
                for pattern in patterns_apo:
                    matches = list(search_dir.glob(pattern))
                    if matches:
                        detected['proj_2d_apo'][f'{pc1}{pc2}'] = matches[0]
                        logger.info(f"Found APO 2D projection PC{pc1}-PC{pc2}: {matches[0]}")
                        break
                if f'{pc1}{pc2}' in detected['proj_2d_apo']:
                    break
            
            # HOLO 2D projections - try patterns with ligand name first, then without
            if ligand_name:
                patterns_holo = [
                    f'2dproj_{pc1}{pc2}_{ligand_name}.xvg',
                    f'2dproj_{pc1}{pc2}_{ligand_name.upper()}.xvg',
                    f'2dproj_{pc1}{pc2}_{ligand_name.lower()}.xvg',
                    f'2dproj_{pc1}{pc2}.xvg'  # Fallback: no suffix (in HOLO directory)
                ]
            else:
                patterns_holo = [
                    f'2dproj_{pc1}{pc2}_*.xvg',
                    f'2dproj_{pc1}{pc2}.xvg'  # Fallback: no suffix
                ]
            
            for search_dir in search_dirs_holo:
                for pattern in patterns_holo:
                    matches = list(search_dir.glob(pattern))
                    # Filter out APO files when looking for HOLO projections.
                    matches = [m for m in matches if 'Apo' not in m.name and 'APO' not in m.name and 'apo' not in m.name]
                    if matches:
                        detected['proj_2d_holo'][f'{pc1}{pc2}'] = matches[0]
                        logger.info(f"Found HOLO 2D projection PC{pc1}-PC{pc2}: {matches[0]}")
                        break
                if f'{pc1}{pc2}' in detected['proj_2d_holo']:
                    break
        
        # Find 3D projections - handle both naming conventions
        apo_3d_patterns = [
            '3dproj_123_Apo.pdb',
            '3dproj_123_APO.pdb',
            '3dproj_123_apo.pdb',
            '3dproj_123.pdb'  # Fallback: no suffix (distinguished by directory)
        ]
        for search_dir in search_dirs_apo:
            for pattern in apo_3d_patterns:
                matches = list(search_dir.glob(pattern))
                if matches:
                    detected['proj_3d_apo'] = matches[0]
                    logger.info(f"Found APO 3D projection: {detected['proj_3d_apo']}")
                    break
            if detected['proj_3d_apo']:
                break
        
        # HOLO 3D projections
        if ligand_name:
            holo_3d_patterns = [
                f'3dproj_123_{ligand_name}.pdb',
                f'3dproj_123_{ligand_name.upper()}.pdb',
                f'3dproj_123_{ligand_name.lower()}.pdb',
                '3dproj_123.pdb'  # Fallback: no suffix (in HOLO directory)
            ]
        else:
            holo_3d_patterns = [
                '3dproj_123_*.pdb',
                '3dproj_123.pdb'  # Fallback: no suffix
            ]
        
        for search_dir in search_dirs_holo:
            for pattern in holo_3d_patterns:
                matches = list(search_dir.glob(pattern))
                # Filter out APO files when looking for HOLO projections.
                matches = [m for m in matches if 'Apo' not in m.name and 'APO' not in m.name and 'apo' not in m.name]
                if matches:
                    detected['proj_3d_holo'] = matches[0]
                    logger.info(f"Found HOLO 3D projection: {detected['proj_3d_holo']}")
                    break
            if detected['proj_3d_holo']:
                break

        # Detect proj_eig files for FEL overlay
        for pair_key in ["12", "13", "23"]:
            proj_eig_name = f"proj_eig_{pair_key}.xvg"
            proj_eig_path = base_dir / proj_eig_name
            if proj_eig_path.exists():
                detected['proj_eig'][pair_key] = proj_eig_path
                logger.info(f"Found projection eigenvector file: {proj_eig_path}")

        # Detect FEL XPM files (gmx sham outputs)
        if self.config.fel.enabled:
            for pair_key in self.config.fel.pairs:
                pair_key = str(pair_key)
                fel_entry: Dict[str, Path] = {}
                # Gibbs/free energy surface
                gibbs_patterns = [
                    f"gibbs_{pair_key}.xpm",
                    "gibbs.xpm" if pair_key == "12" else None,
                ]
                for pattern in [p for p in gibbs_patterns if p]:
                    matches = list(base_dir.glob(pattern))
                    if matches:
                        fel_entry["gibbs"] = matches[0]
                        logger.info(f"Found FEL Gibbs surface ({pair_key}): {matches[0]}")
                        break

                # Probability
                prob_patterns = [
                    f"prob_gibbs_{pair_key}.xpm",
                    "prob.xpm" if pair_key == "12" else None,
                ]
                for pattern in [p for p in prob_patterns if p]:
                    matches = list(base_dir.glob(pattern))
                    if matches:
                        fel_entry["probability"] = matches[0]
                        logger.info(f"Found FEL probability ({pair_key}): {matches[0]}")
                        break

                # Entropy
                entropy_patterns = [
                    f"entropy_gibbs_{pair_key}.xpm",
                    "entropy.xpm" if pair_key == "12" else None,
                ]
                for pattern in [p for p in entropy_patterns if p]:
                    matches = list(base_dir.glob(pattern))
                    if matches:
                        fel_entry["entropy"] = matches[0]
                        logger.info(f"Found FEL entropy ({pair_key}): {matches[0]}")
                        break

                # Enthalpy
                enthalpy_patterns = [
                    f"enthalpy_gibbs_{pair_key}.xpm",
                    "enthalpy.xpm" if pair_key == "12" else None,
                ]
                for pattern in [p for p in enthalpy_patterns if p]:
                    matches = list(base_dir.glob(pattern))
                    if matches:
                        fel_entry["enthalpy"] = matches[0]
                        logger.info(f"Found FEL enthalpy ({pair_key}): {matches[0]}")
                        break

                if fel_entry:
                    detected["fel"][pair_key] = fel_entry
        
        return detected
    
    def _extract_ligand_name(self, vis_dir: Path) -> Optional[str]:
        """
        Extract ligand name from 2D projection filenames.
        
        Looks for patterns like: 2dproj_12_LIGAND.xvg
        
        Args:
            vis_dir: PCA_Vis directory
            
        Returns:
            Ligand name or None
        """
        # Look for 2D projection files
        proj_files = list(vis_dir.glob('2dproj_*_*.xvg'))
        
        for file in proj_files:
            # Skip APO files
            if 'Apo' in file.name or 'APO' in file.name or 'apo' in file.name:
                continue
            
            # Extract ligand name from pattern: 2dproj_12_LIGAND.xvg
            match = re.search(r'2dproj_\d+_(.+)\.xvg', file.name)
            if match:
                ligand = match.group(1)
                # Clean up ligand name
                ligand = ligand.replace('_', '').strip()
                if ligand:
                    return ligand
        
        return None
    
    def load_all_data(self, detected_files: Dict) -> Dict[str, Dict]:
        """
        Load and parse all detected files.
        
        Args:
            detected_files: Dict from auto_detect_files()
            
        Returns:
            Dict with parsed data:
            {
                'eigenvals': {...},
                'proj_1d': {...},
                'cc_holo': {...},
                'cc_apo': {...},
                'proj_2d_apo': {12: {...}, 13: {...}, 23: {...}},
                'proj_2d_holo': {12: {...}, 13: {...}, 23: {...}},
                'proj_3d_apo': (x, y, z),
                'proj_3d_holo': (x, y, z)
            }
        """
        all_data = {}
        
        # Parse eigenvalues
        if detected_files['eigenvals']:
            title, x, y, x_t, y_t = self.parser.parse_xvg_single(detected_files['eigenvals'])
            if x and y:
                all_data['eigenvals'] = {
                    'title': title,
                    'x': x,
                    'y': y,
                    'x_title': x_t,
                    'y_title': y_t
                }
                logger.info(f"Parsed eigenvalues: {len(x)} points")
        
        # Parse 1D projection
        if detected_files['proj_1d']:
            titles, xs, ys, x_ts, y_ts = self.parser.parse_xvg_multiple(detected_files['proj_1d'])
            if titles and xs and ys:
                all_data['proj_1d'] = {
                    'titles': titles,
                    'xs': xs,
                    'ys': ys,
                    'x_titles': x_ts,
                    'y_titles': y_ts
                }
                logger.info(f"Parsed 1D projection: {len(titles)} PCs")
        
        # Parse cosine content
        if detected_files['cc_holo']:
            title, x, y, x_t, y_t = self.parser.parse_xvg_single(detected_files['cc_holo'])
            if x and y:
                all_data['cc_holo'] = {
                    'title': title,
                    'x': x,
                    'y': y,
                    'x_title': x_t,
                    'y_title': y_t
                }
                logger.info(f"Parsed HOLO cosine content: {len(x)} points")
        
        if detected_files['cc_apo']:
            title, x, y, x_t, y_t = self.parser.parse_xvg_single(detected_files['cc_apo'])
            if x and y:
                all_data['cc_apo'] = {
                    'title': title,
                    'x': x,
                    'y': y,
                    'x_title': x_t,
                    'y_title': y_t
                }
                logger.info(f"Parsed APO cosine content: {len(x)} points")
        
        # Parse 2D projections
        all_data['proj_2d_apo'] = {}
        all_data['proj_2d_holo'] = {}
        
        for pair_key in ['12', '13', '23']:
            # APO
            if pair_key in detected_files['proj_2d_apo']:
                title, x, y, x_t, y_t = self.parser.parse_xvg_single(
                    detected_files['proj_2d_apo'][pair_key]
                )
                if x and y:
                    all_data['proj_2d_apo'][pair_key] = {
                        'title': title,
                        'x': x,
                        'y': y,
                        'x_title': x_t,
                        'y_title': y_t
                    }
                    logger.info(f"Parsed APO 2D projection PC{pair_key[0]}-PC{pair_key[1]}: {len(x)} points")
            
            # HOLO
            if pair_key in detected_files['proj_2d_holo']:
                title, x, y, x_t, y_t = self.parser.parse_xvg_single(
                    detected_files['proj_2d_holo'][pair_key]
                )
                if x and y:
                    all_data['proj_2d_holo'][pair_key] = {
                        'title': title,
                        'x': x,
                        'y': y,
                        'x_title': x_t,
                        'y_title': y_t
                    }
                    logger.info(f"Parsed HOLO 2D projection PC{pair_key[0]}-PC{pair_key[1]}: {len(x)} points")
        
        # Parse 3D projections
        if detected_files['proj_3d_apo']:
            x, y, z = self.parser.parse_pdb_coords(detected_files['proj_3d_apo'])
            if x and y and z:
                all_data['proj_3d_apo'] = (x, y, z)
                logger.info(f"Parsed APO 3D projection: {len(x)} points")

        if detected_files['proj_3d_holo']:
            x, y, z = self.parser.parse_pdb_coords(detected_files['proj_3d_holo'])
            if x and y and z:
                all_data['proj_3d_holo'] = (x, y, z)
                logger.info(f"Parsed HOLO 3D projection: {len(x)} points")

        # Parse proj_eig files for FEL overlays
        if detected_files.get('proj_eig'):
            all_data['proj_eig'] = {}
            for pair_key, file_path in detected_files['proj_eig'].items():
                titles, xs, ys, x_ts, y_ts = self.parser.parse_xvg_multiple(file_path)
                if ys:
                    all_data['proj_eig'][pair_key] = {
                        'titles': titles,
                        'xs': xs,
                        'ys': ys,
                        'x_titles': x_ts,
                        'y_titles': y_ts
                    }
                    logger.info(f"Parsed proj_eig {pair_key}: {len(ys)} sets")

        # Parse FEL XPM files
        if detected_files.get('fel'):
            all_data['fel'] = {}
            for pair_key, fel_entry in detected_files['fel'].items():
                fel_data: Dict[str, Dict] = {}
                for key, path in fel_entry.items():
                    matrix_vals, x_coords, y_coords = self.parser.parse_xpm(path)
                    if matrix_vals:
                        fel_data[key] = {
                            'matrix': matrix_vals,
                            'x_coords': x_coords,
                            'y_coords': y_coords,
                            'path': path
                        }
                        logger.info(f"Parsed FEL {key} ({pair_key}): {len(matrix_vals)} rows")
                if fel_data:
                    all_data['fel'][pair_key] = fel_data
        
        return all_data

    def build_publication_payload(self, all_data: Dict, detected_files: Dict) -> Dict:
        """Prepare numerical sources and explicit scientific availability records."""
        from .publication import imported_payload
        return imported_payload(self, all_data, detected_files)
