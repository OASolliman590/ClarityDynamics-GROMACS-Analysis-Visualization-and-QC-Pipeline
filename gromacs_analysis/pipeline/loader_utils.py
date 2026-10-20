"""Path and default-directory helpers shared by the pipeline YAML builders."""

from pathlib import Path
from typing import Optional


def _resolve_path(value: Optional[str], base_path: Path) -> Optional[Path]:
    if not value:
        return None
    path = Path(value)
    if not path.is_absolute():
        path = (base_path / path).resolve()
    return path


def _path_parts(path: Path) -> tuple:
    """Return normalized path parts without a leading './'."""
    return tuple(part for part in path.parts if part not in ("", "."))


def _is_prefixed_by_base(value: str, raw_base_dir: Optional[str]) -> bool:
    """Return True when a relative value already includes the stage base_dir."""
    if not raw_base_dir:
        return False
    value_path = Path(value)
    base_path = Path(str(raw_base_dir))
    if value_path.is_absolute() or base_path.is_absolute():
        return False

    value_parts = _path_parts(value_path)
    base_parts = _path_parts(base_path)
    if not base_parts or len(value_parts) < len(base_parts):
        return False
    return value_parts[:len(base_parts)] == base_parts


def _resolve_stage_input_path(
    value: str,
    base_path: Path,
    base_dir: Optional[Path],
    raw_base_dir: Optional[str],
) -> Path:
    """
    Resolve stage input files consistently.

    Relative paths that already include the configured base_dir are resolved
    from the YAML file location. Short paths such as "rep1/traj.xtc" are
    resolved from the stage base_dir.
    """
    path = Path(value)
    if path.is_absolute():
        return path
    if base_dir and not _is_prefixed_by_base(value, raw_base_dir):
        return (base_dir / path).resolve()
    return (base_path / path).resolve()


def _default_output_dir(stage: str, output_root: Optional[Path], base_dir: Optional[Path]) -> Optional[Path]:
    if output_root:
        return output_root / stage
    if stage == "qc":
        return Path("./output_qc")
    if stage == "md":
        return Path("./output")
    if stage == "md_compare":
        return Path("./output_md_compare")
    if stage == "mmpbsa":
        return Path("./output_mmpbsa")
    if stage == "pca":
        return base_dir / "PCA_Vis" if base_dir else None
    if stage == "pca_compare":
        return Path("./output_pca_compare")
    if stage == "distance":
        return Path("./output_distance")
    if stage == "distance_compare":
        return Path("./output_distance_compare")
    return None
