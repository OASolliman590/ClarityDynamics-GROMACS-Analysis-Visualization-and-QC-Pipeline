"""
YAML loader/saver for pipeline configurations.
"""

import logging
from pathlib import Path
from typing import Any, Dict, Optional, List

import yaml

from .config import PipelineConfig, PipelineStage
from ..utils.palette import categorical
from .loader_utils import (  # noqa: F401 - re-exported for backward compatibility
    _default_output_dir,
    _is_prefixed_by_base,
    _path_parts,
    _resolve_path,
    _resolve_stage_input_path,
)
from .stage_builders import (  # noqa: F401 - re-exported for backward compatibility
    _batch_systems,
    _build_pca_joint_config,
    _build_prolif_batch_configs,
    _build_prolif_compare_config,
    _build_prolif_config,
    _build_ttclust_batch_configs,
    _build_ttclust_config,
)
from .templates import PIPELINE_TEMPLATE
from ..md.config import MDConfig, SystemConfig, AnalysisMetric, PlotConfig, StabilityConfig, RMSFTerminalConfig
from ..mmpbsa.config import MMPBSAConfig, MMPBSASystemConfig
from ..pca.config import PCAConfig, PCASystemConfig, PCACompareConfig, FELConfig, ClusteringConfig
from ..distance.config import DistanceConfig, DistanceCompareConfig
from ..interactions.config import InteractionConfig, InteractionCompareConfig
from ..ranking.config import RankingConfig
from ..prolif.config import ProlifConfig, ProlifCompareConfig
from ..ttclust.config import TTClustConfig
from ..qc.config import QCConfig
from ..config.yaml_loader import load_yaml_config as load_md_yaml, parse_amino_acids_config
from ..mmpbsa.yaml_loader import load_yaml_config as load_mmpbsa_yaml
from ..pca.yaml_loader import load_yaml_config as load_pca_yaml, parse_run_gromacs_config
from ..distance.yaml_loader import load_yaml_config as load_distance_yaml

logger = logging.getLogger(__name__)



def _apply_relative_paths(config: Any, base_path: Path) -> None:
    """Resolve relative paths inside config objects based on base_path."""
    for attr in ("base_dir", "output_dir", "apo_base_dir"):
        if hasattr(config, attr):
            current = getattr(config, attr)
            if isinstance(current, Path) and not current.is_absolute():
                setattr(config, attr, (base_path / current).resolve())


def _clone_system(system: SystemConfig) -> SystemConfig:
    return SystemConfig(
        name=system.name,
        dir_pattern=system.dir_pattern,
        dir_names=list(system.dir_names) if system.dir_names else None,
        is_apo=system.is_apo,
        replicates=system.replicates,
        color=tuple(system.color) if system.color else None,
    )


def _assign_palette_if_needed(systems: List[SystemConfig], palette: Optional[List[tuple]] = None) -> None:
    """Apply a user-supplied comparison palette to ligand systems without explicit colours.

    Without a palette the plotters assign the shared validated palette, so nothing is set here.
    """
    if not palette:
        return
    holo_systems = [s for s in systems if not s.is_apo]
    if len(holo_systems) <= 1:
        return
    if all(s.color is None for s in holo_systems):
        for idx, sys in enumerate(holo_systems):
            sys.color = tuple(palette[idx % len(palette)])


def _build_md_config(data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]) -> MDConfig:
    base_dir = data.get("base_dir") or defaults.get("data_root")
    if not base_dir:
        raise ValueError("MD config requires base_dir (or paths.data_root)")
    base_dir = _resolve_path(base_dir, base_path)

    output_dir = data.get("output_dir") or _default_output_dir("md", defaults.get("output_root"), base_dir)
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None

    protein_name = data.get("protein_name") or defaults.get("protein_name")
    if not protein_name:
        raise ValueError("MD config requires protein_name (or project.protein_name)")

    systems_data = data.get("systems", [])
    if not systems_data:
        raise ValueError("MD config requires at least one system entry")
    systems = [
        SystemConfig(
            name=s["name"],
            dir_pattern=s.get("dir_pattern", ""),
            dir_names=s.get("dir_names"),
            is_apo=s.get("is_apo", False),
            replicates=s.get("replicates", 1),
            color=tuple(s["color"]) if "color" in s else None,
        )
        for s in systems_data
    ]

    metrics_data = data.get("metrics", [])
    metrics = []
    if metrics_data:
        for m in metrics_data:
            metrics.append(
                AnalysisMetric(
                    name=m["name"],
                    file_pattern=m["file_pattern"],
                    title=m["title"],
                    ylabel=m["ylabel"],
                    is_holo_only=m.get("is_holo_only", False),
                    data_format=m.get("data_format", "two_column"),
                    ligand=m.get("ligand"),
                    unit=m.get("unit"),
                )
            )

    plot_data = data.get("plot_config", {})
    plot_config = PlotConfig(
        template=plot_data.get("template", "plotly_white"),
        style=plot_data.get("style", "simple"),
        width=plot_data.get("width", 1200),
        height=plot_data.get("height", 600),
        scale=plot_data.get("scale", 2),
        font_family=plot_data.get("font_family", "Arial"),
        font_size=plot_data.get("font_size", 16),
        colors={k: tuple(v) for k, v in plot_data.get("colors", {}).items()},
        save_formats=plot_data.get("save_formats", ["html", "svg"]),
    )

    residue_range = None
    if "residue_range" in data:
        residue_range = tuple(data["residue_range"])

    amino_acids = None
    if "amino_acids" in data:
        amino_acids = parse_amino_acids_config(data["amino_acids"])

    stability_data = data.get("stability", {})
    stability_config = StabilityConfig(
        enabled=bool(stability_data.get("enabled", False)),
        metric=stability_data.get("metric", "rmsd_prot"),
        window=stability_data.get("window", 50),
        std_threshold=stability_data.get("std_threshold", 0.2),
        slope_threshold=stability_data.get("slope_threshold", 0.01),
        min_points=stability_data.get("min_points", 100),
        apply_to_all_metrics=bool(stability_data.get("apply_to_all_metrics", True)),
        method=stability_data.get("method", "chodera"),
        n_candidates=int(stability_data.get("n_candidates", 100)),
    )

    sequence_topology = data.get("sequence_topology")
    sequence_topology = _resolve_path(sequence_topology, base_path) if sequence_topology else None
    sequence_selection = data.get("sequence_selection", "protein")

    return MDConfig(
        base_dir=base_dir,
        output_dir=output_dir,
        protein_name=protein_name,
        systems=systems,
        metrics=metrics,
        plot_config=plot_config,
        residue_range=residue_range,
        amino_acids=amino_acids,
        stability=stability_config,
        sequence_topology=sequence_topology,
        sequence_selection=sequence_selection,
        time_unit_input=data.get("time_unit_input", "auto"),
        time_unit_output=data.get("time_unit_output", "ns"),
        time_step_ps=data.get("time_step_ps"),
        time_scale=data.get("time_scale"),
        smoothing_window=data.get("smoothing_window", data.get("analysis", {}).get("smoothing_window", 1)),
        outlier_threshold=data.get("outlier_threshold", data.get("analysis", {}).get("outlier_threshold")),
        ligands=data.get("ligands") or [],
        rmsf_terminals=RMSFTerminalConfig(**(data.get("rmsf_terminals") or {})),
        convert_units=bool(data.get("convert_units", True)),
        convergence={k: float(v) for k, v in (data.get("convergence") or {}).items()},
    )


def _build_mmpbsa_config(data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]) -> MMPBSAConfig:
    base_dir = data.get("base_dir") or defaults.get("data_root")
    if not base_dir:
        raise ValueError("MMPBSA config requires base_dir (or paths.data_root)")
    base_dir = _resolve_path(base_dir, base_path)

    output_dir = data.get("output_dir") or _default_output_dir("mmpbsa", defaults.get("output_root"), base_dir)
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None

    protein_name = data.get("protein_name") or defaults.get("protein_name")
    if not protein_name:
        raise ValueError("MMPBSA config requires protein_name (or project.protein_name)")

    ligand_name = data.get("ligand_name")
    if not ligand_name:
        raise ValueError("MMPBSA config requires ligand_name")

    systems_data = data.get("systems", [])
    if not systems_data:
        raise ValueError("MMPBSA config requires at least one system entry")
    systems = [
        MMPBSASystemConfig(
            name=s["name"],
            dir_pattern=s["dir_pattern"],
            replicates=s.get("replicates", 1),
            results_file_pattern=s.get("results_file_pattern", "FINAL_RESULTS_MMPBSA*.dat"),
            decomp_file_pattern=s.get("decomp_file_pattern", "FINAL_DECOMP_MMPBSA*.dat"),
        )
        for s in systems_data
    ]

    plot_data = data.get("plot_config", {})
    plot_config = PlotConfig(
        template=plot_data.get("template", "ggplot2"),
        style=plot_data.get("style", "publication"),
        width=plot_data.get("width", 1400),
        height=plot_data.get("height", 700),
        scale=plot_data.get("scale", 2),
        font_family=plot_data.get("font_family", "Times New Roman"),
        font_size=plot_data.get("font_size", 24),
        colors=plot_data.get("colors", {}),
        save_formats=plot_data.get("save_formats", ["html", "svg"]),
    )

    amino_acids = None
    if "amino_acids" in data:
        aa_cfg = data["amino_acids"]
        if "custom" in aa_cfg:
            amino_acids = aa_cfg["custom"]
        elif "range" in aa_cfg:
            start, end = aa_cfg["range"]
            amino_acids = list(range(start, end + 1))
        elif "ranges" in aa_cfg:
            amino_acids = []
            for start, end in aa_cfg["ranges"]:
                amino_acids.extend(range(start, end + 1))

    return MMPBSAConfig(
        base_dir=base_dir,
        output_dir=output_dir,
        protein_name=protein_name,
        ligand_name=ligand_name,
        systems=systems,
        file_format=data.get("file_format", "auto"),
        plot_config=plot_config,
        amino_acids=amino_acids,
        compare_systems=bool(data.get("compare_systems", False)),
        compare_binding_key=data.get("compare_binding_key", "TOTAL"),
        compare_components=bool(data.get("compare_components", False)),
        compare_component_order=data.get("compare_component_order"),
        incomplete_replicates_policy=data.get("incomplete_replicates_policy", "warn"),
        start_frame=data.get("start_frame"),
        require_decomp_replicates=bool(data.get("require_decomp_replicates", False)),
    )


def _build_qc_config(data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]) -> QCConfig:
    raw_base_dir = data.get("base_dir") or defaults.get("data_root") or "."
    base_dir = _resolve_path(str(raw_base_dir), base_path)

    output_dir = data.get("output_dir") or _default_output_dir("qc", defaults.get("output_root"), base_dir)
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None
    if output_dir is None:
        raise ValueError("qc config requires output_dir or output_root")

    topology = data.get("topology")
    topology_path = _resolve_stage_input_path(topology, base_path, base_dir, str(raw_base_dir)) if topology else None

    index_file = data.get("index_file")
    index_path = _resolve_stage_input_path(index_file, base_path, base_dir, str(raw_base_dir)) if index_file else None

    trajectories = [
        _resolve_stage_input_path(path, base_path, base_dir, str(raw_base_dir))
        for path in data.get("trajectories", [])
    ]

    return QCConfig(
        base_dir=base_dir,
        output_dir=output_dir,
        topology=topology_path,
        trajectories=trajectories,
        index_file=index_path,
        expected_groups=data.get("expected_groups", []),
        protein_selection=data.get("protein_selection", "protein"),
        ligand_selection=data.get("ligand_selection", "resname UNK"),
        expected_replicates=data.get("expected_replicates"),
        min_distance_cutoff=data.get("min_distance_cutoff"),
        dissociation_cutoff=data.get("dissociation_cutoff", 5.0),
        check_frame_times=bool(data.get("check_frame_times", True)),
        sample_stride=data.get("sample_stride", 10),
        max_frames=data.get("max_frames", 500),
        fail_on_error=bool(data.get("fail_on_error", False)),
    )


def _build_pca_config(data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]) -> PCAConfig:
    base_dir = data.get("base_dir") or defaults.get("data_root")
    if not base_dir:
        raise ValueError("PCA config requires base_dir (or paths.data_root)")
    base_dir = _resolve_path(base_dir, base_path)

    output_dir = data.get("output_dir") or _default_output_dir("pca", defaults.get("output_root"), base_dir)
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None

    protein_name = data.get("protein_name") or defaults.get("protein_name") or "Protein"
    ligand_name = data.get("ligand_name")

    apo_base_dir = data.get("apo_base_dir")
    apo_base_dir = _resolve_path(apo_base_dir, base_path) if apo_base_dir else None

    n_pcs = data.get("n_pcs", 10)

    systems = []
    for sys_data in data.get("systems", []):
        systems.append(
            PCASystemConfig(
                name=sys_data["name"],
                is_apo=sys_data.get("is_apo", False),
                ligand_name=sys_data.get("ligand_name"),
            )
        )

    plot_data = data.get("plot_config", {})
    plot_config = PlotConfig(
        template=plot_data.get("template", "plotly_white"),
        width=plot_data.get("width", 1200),
        height=plot_data.get("height", 600),
        scale=plot_data.get("scale", 2),
        font_family=plot_data.get("font_family", "Times New Roman"),
        font_size=plot_data.get("font_size", 24),
        save_formats=plot_data.get("save_formats", ["html", "svg"]),
    )

    fel_data = data.get("fel", {})
    fel_config = FELConfig(
        enabled=fel_data.get("enabled", False),
        pairs=fel_data.get("pairs", ["12"]),
        contours=fel_data.get("contours", 7),
        step=fel_data.get("step", 2),
        colorscale=fel_data.get("colorscale", "Viridis"),
        overlay_projection=fel_data.get("overlay_projection", True),
        include_probability=fel_data.get("include_probability", False),
        include_entropy=fel_data.get("include_entropy", False),
        include_enthalpy=fel_data.get("include_enthalpy", False),
    )

    cluster_data = data.get("clustering", {})
    clustering_config = ClusteringConfig(
        enabled=cluster_data.get("enabled", False),
        method=cluster_data.get("method", "kmeans"),
        pair=cluster_data.get("pair", "12"),
        n_clusters=cluster_data.get("n_clusters", 3),
        eps=cluster_data.get("eps", 1.5),
        min_samples=cluster_data.get("min_samples", 5),
        downsample=cluster_data.get("downsample", 1),
        max_points=cluster_data.get("max_points"),
    )

    run_gromacs_config = parse_run_gromacs_config(data)
    if run_gromacs_config.input_dir and not run_gromacs_config.input_dir.is_absolute():
        run_gromacs_config.input_dir = (base_path / run_gromacs_config.input_dir).resolve()

    return PCAConfig(
        base_dir=base_dir,
        output_dir=output_dir,
        protein_name=protein_name,
        ligand_name=ligand_name,
        apo_base_dir=apo_base_dir,
        systems=systems,
        n_pcs=n_pcs,
        color=data.get("color"),
        plot_config=plot_config,
        fel=fel_config,
        clustering=clustering_config,
        run_gromacs=run_gromacs_config,
    )


def _parse_residue_config(residue_data: Any) -> Dict[str, Any]:
    residue_ids = None
    residue_selection = None
    auto_detect_residues = None
    auto_detect_cutoff = None
    auto_detect_stride = None
    auto_detect_max_frames = None
    auto_detect_selection = None
    auto_detect_method = None

    if isinstance(residue_data, dict):
        if "auto" in residue_data or "auto_detect" in residue_data:
            auto_flag = residue_data.get("auto")
            if auto_flag is None:
                auto_flag = residue_data.get("auto_detect")
            auto_detect_residues = bool(auto_flag)
            if "cutoff" in residue_data:
                auto_detect_cutoff = residue_data.get("cutoff")
            if "stride" in residue_data:
                auto_detect_stride = residue_data.get("stride")
            if "max_frames" in residue_data:
                auto_detect_max_frames = residue_data.get("max_frames")
            if "selection" in residue_data:
                auto_detect_selection = residue_data.get("selection")
            if "method" in residue_data:
                auto_detect_method = residue_data.get("method")

        if auto_detect_residues is not True:
            if "selection" in residue_data:
                residue_selection = residue_data["selection"]
            else:
                residue_ids = parse_amino_acids_config(residue_data)
    elif isinstance(residue_data, list):
        residue_ids = []
        for item in residue_data:
            if isinstance(item, int):
                residue_ids.append(item)
            elif isinstance(item, str) and ":" in item:
                residue_ids.append(item.strip())  # chain-qualified, e.g. "A:45"
            elif isinstance(item, str):
                digits = "".join(ch for ch in item if ch.isdigit())
                if digits:
                    residue_ids.append(int(digits))
    elif isinstance(residue_data, str):
        residue_selection = residue_data

    return {
        "residue_ids": residue_ids,
        "residue_selection": residue_selection,
        "auto_detect_residues": auto_detect_residues,
        "auto_detect_cutoff": auto_detect_cutoff,
        "auto_detect_stride": auto_detect_stride,
        "auto_detect_max_frames": auto_detect_max_frames,
        "auto_detect_selection": auto_detect_selection,
        "auto_detect_method": auto_detect_method,
    }


def _expand_trajectories(
    data: Dict[str, Any],
    base_path: Path,
    base_dir: Optional[Path],
    raw_base_dir: Optional[str] = None,
) -> List[Path]:
    trajectories = data.get("trajectories") or []
    if trajectories:
        return [
            _resolve_stage_input_path(p, base_path, base_dir, raw_base_dir)
            for p in trajectories
        ]

    pattern = data.get("trajectory_pattern")
    replicates = data.get("replicates")
    if pattern and replicates:
        trajs = []
        for i in range(1, int(replicates) + 1):
            if "{}" in pattern:
                value = pattern.format(i)
            else:
                value = pattern.format(rep=i)
            trajs.append(_resolve_stage_input_path(value, base_path, base_dir, raw_base_dir))
        return trajs

    return []


def _build_distance_config(data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]) -> DistanceConfig:
    raw_base_dir = data.get("base_dir") or defaults.get("data_root")
    if not raw_base_dir:
        raise ValueError("distance config requires base_dir (or paths.data_root)")
    base_dir = _resolve_path(raw_base_dir, base_path)

    output_dir = data.get("output_dir") or _default_output_dir("distance", defaults.get("output_root"), base_dir)
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None

    protein_name = data.get("protein_name") or defaults.get("protein_name") or "Protein"
    ligand_name = data.get("ligand_name") or data.get("name") or "Ligand"

    topology = data.get("topology")
    if not topology:
        raise ValueError("distance config requires topology")
    topology = _resolve_stage_input_path(topology, base_path, base_dir, raw_base_dir)

    trajectories = _expand_trajectories(data, base_path, base_dir, raw_base_dir)
    if not trajectories:
        raise ValueError("distance config requires trajectories (or trajectory_pattern + replicates)")

    residue_cfg = _parse_residue_config(data.get("residues") or data.get("residue_ids"))
    auto_detect_residues = residue_cfg.get("auto_detect_residues")
    if "auto_detect_residues" in data:
        auto_detect_residues = bool(data.get("auto_detect_residues"))
    if auto_detect_residues is None:
        auto_detect_residues = False

    auto_detect_cutoff = data.get("auto_detect_cutoff")
    if auto_detect_cutoff is None:
        auto_detect_cutoff = residue_cfg.get("auto_detect_cutoff")
    if auto_detect_cutoff is None:
        auto_detect_cutoff = 5.0

    auto_detect_stride = data.get("auto_detect_stride")
    if auto_detect_stride is None:
        auto_detect_stride = residue_cfg.get("auto_detect_stride")
    if auto_detect_stride is None:
        auto_detect_stride = 1

    auto_detect_max_frames = data.get("auto_detect_max_frames")
    if auto_detect_max_frames is None:
        auto_detect_max_frames = residue_cfg.get("auto_detect_max_frames")

    auto_detect_selection = data.get("auto_detect_selection")
    if auto_detect_selection is None:
        auto_detect_selection = residue_cfg.get("auto_detect_selection")
    if auto_detect_selection is None:
        auto_detect_selection = "protein"

    auto_detect_method = data.get("auto_detect_method")
    if auto_detect_method is None:
        auto_detect_method = residue_cfg.get("auto_detect_method")
    if auto_detect_method is None:
        auto_detect_method = "min"

    plot_data = data.get("plot_config", {})
    plot_config = PlotConfig(
        template=plot_data.get("template", "ggplot2"),
        style=plot_data.get("style", "simple"),
        width=plot_data.get("width", 1200),
        height=plot_data.get("height", 600),
        scale=plot_data.get("scale", 2),
        font_family=plot_data.get("font_family", "Times New Roman"),
        font_size=plot_data.get("font_size", 24),
        colors={k: tuple(v) for k, v in plot_data.get("colors", {}).items()},
        save_formats=plot_data.get("save_formats", ["html", "svg"]),
    )

    return DistanceConfig(
        base_dir=base_dir,
        output_dir=output_dir,
        protein_name=protein_name,
        ligand_name=ligand_name,
        topology=topology,
        trajectories=trajectories,
        ligand_selection=data.get("ligand_selection", "resname UNK"),
        residue_ids=residue_cfg["residue_ids"],
        residue_selection=residue_cfg["residue_selection"],
        auto_detect_residues=auto_detect_residues,
        auto_detect_cutoff=auto_detect_cutoff,
        auto_detect_stride=auto_detect_stride,
        auto_detect_max_frames=auto_detect_max_frames,
        auto_detect_selection=auto_detect_selection,
        auto_detect_method=auto_detect_method,
        method=data.get("method", "com"),
        distance_scale=data.get("distance_scale", 1.0),
        time_unit=data.get("time_unit", "auto"),
        time_step_ns=data.get("time_step_ns"),
        per_residue_plots=data.get("per_residue_plots", True),
        combined_plot=data.get("combined_plot", True),
        combined_show_std=data.get("combined_show_std", False),
        include_overall=data.get("include_overall", False),
        plot_config=plot_config,
    )


def _build_distance_batch_configs(stage_data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]) -> List[DistanceConfig]:
    raw_base_dir = stage_data.get("base_dir") or defaults.get("data_root")
    if not raw_base_dir:
        raise ValueError("distance_batch requires base_dir (or paths.data_root)")
    base_dir = _resolve_path(raw_base_dir, base_path)

    output_dir = stage_data.get("output_dir") or _default_output_dir("distance", defaults.get("output_root"), base_dir)
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None

    protein_name = stage_data.get("protein_name") or defaults.get("protein_name") or "Protein"

    systems_data = stage_data.get("systems", [])
    if not systems_data:
        raise ValueError("distance_batch requires systems list")

    residue_cfg = _parse_residue_config(stage_data.get("residues") or stage_data.get("residue_ids"))
    stage_auto_detect_residues = residue_cfg.get("auto_detect_residues")
    if "auto_detect_residues" in stage_data:
        stage_auto_detect_residues = bool(stage_data.get("auto_detect_residues"))
    if stage_auto_detect_residues is None:
        stage_auto_detect_residues = False

    stage_auto_detect_cutoff = stage_data.get("auto_detect_cutoff")
    if stage_auto_detect_cutoff is None:
        stage_auto_detect_cutoff = residue_cfg.get("auto_detect_cutoff")
    if stage_auto_detect_cutoff is None:
        stage_auto_detect_cutoff = 5.0

    stage_auto_detect_stride = stage_data.get("auto_detect_stride")
    if stage_auto_detect_stride is None:
        stage_auto_detect_stride = residue_cfg.get("auto_detect_stride")
    if stage_auto_detect_stride is None:
        stage_auto_detect_stride = 1

    stage_auto_detect_max_frames = stage_data.get("auto_detect_max_frames")
    if stage_auto_detect_max_frames is None:
        stage_auto_detect_max_frames = residue_cfg.get("auto_detect_max_frames")

    stage_auto_detect_selection = stage_data.get("auto_detect_selection")
    if stage_auto_detect_selection is None:
        stage_auto_detect_selection = residue_cfg.get("auto_detect_selection")
    if stage_auto_detect_selection is None:
        stage_auto_detect_selection = "protein"

    stage_auto_detect_method = stage_data.get("auto_detect_method")
    if stage_auto_detect_method is None:
        stage_auto_detect_method = residue_cfg.get("auto_detect_method")
    if stage_auto_detect_method is None:
        stage_auto_detect_method = "min"

    output_pattern = stage_data.get("output_subdir_pattern", "{name}")
    configs: List[DistanceConfig] = []

    for sys_data in systems_data:
        sys_name = sys_data.get("name")
        if not sys_name:
            raise ValueError("distance_batch systems entries require name")

        sys_raw_base_dir = sys_data.get("base_dir") or raw_base_dir
        sys_base_dir = _resolve_path(sys_raw_base_dir, base_path)
        sys_residue_cfg = _parse_residue_config(sys_data.get("residues") or sys_data.get("residue_ids"))

        topology = sys_data.get("topology") or stage_data.get("topology")
        if not topology:
            raise ValueError(f"distance_batch system '{sys_name}' requires topology (or stage-level topology)")
        topology = _resolve_stage_input_path(topology, base_path, sys_base_dir, sys_raw_base_dir)

        trajectories = _expand_trajectories({**stage_data, **sys_data}, base_path, sys_base_dir, sys_raw_base_dir)
        if not trajectories:
            raise ValueError(f"distance_batch system '{sys_name}' requires trajectories")

        sys_output = output_dir / output_pattern.format(name=sys_name) if output_dir else None

        auto_detect_residues = sys_residue_cfg.get("auto_detect_residues")
        if "auto_detect_residues" in sys_data:
            auto_detect_residues = bool(sys_data.get("auto_detect_residues"))
        if auto_detect_residues is None:
            if sys_residue_cfg.get("residue_ids") or sys_residue_cfg.get("residue_selection"):
                auto_detect_residues = False
            else:
                auto_detect_residues = stage_auto_detect_residues

        auto_detect_cutoff = sys_data.get("auto_detect_cutoff")
        if auto_detect_cutoff is None:
            auto_detect_cutoff = sys_residue_cfg.get("auto_detect_cutoff")
        if auto_detect_cutoff is None:
            auto_detect_cutoff = stage_auto_detect_cutoff
        if auto_detect_cutoff is None:
            auto_detect_cutoff = 5.0

        auto_detect_stride = sys_data.get("auto_detect_stride")
        if auto_detect_stride is None:
            auto_detect_stride = sys_residue_cfg.get("auto_detect_stride")
        if auto_detect_stride is None:
            auto_detect_stride = stage_auto_detect_stride
        if auto_detect_stride is None:
            auto_detect_stride = 1

        auto_detect_max_frames = sys_data.get("auto_detect_max_frames")
        if auto_detect_max_frames is None:
            auto_detect_max_frames = sys_residue_cfg.get("auto_detect_max_frames")
        if auto_detect_max_frames is None:
            auto_detect_max_frames = stage_auto_detect_max_frames

        auto_detect_selection = sys_data.get("auto_detect_selection")
        if auto_detect_selection is None:
            auto_detect_selection = sys_residue_cfg.get("auto_detect_selection")
        if auto_detect_selection is None:
            auto_detect_selection = stage_auto_detect_selection
        if auto_detect_selection is None:
            auto_detect_selection = "protein"

        auto_detect_method = sys_data.get("auto_detect_method")
        if auto_detect_method is None:
            auto_detect_method = sys_residue_cfg.get("auto_detect_method")
        if auto_detect_method is None:
            auto_detect_method = stage_auto_detect_method
        if auto_detect_method is None:
            auto_detect_method = "min"

        plot_data = sys_data.get("plot_config", stage_data.get("plot_config", {}))
        plot_config = PlotConfig(
            template=plot_data.get("template", "ggplot2"),
            style=plot_data.get("style", "simple"),
            width=plot_data.get("width", 1200),
            height=plot_data.get("height", 600),
            scale=plot_data.get("scale", 2),
            font_family=plot_data.get("font_family", "Times New Roman"),
            font_size=plot_data.get("font_size", 24),
            colors={k: tuple(v) for k, v in plot_data.get("colors", {}).items()},
            save_formats=plot_data.get("save_formats", ["html", "svg"]),
        )

        configs.append(
            DistanceConfig(
                base_dir=sys_base_dir,
                output_dir=sys_output,
                protein_name=protein_name,
                ligand_name=sys_data.get("ligand_name", sys_name),
                topology=topology,
                trajectories=trajectories,
                ligand_selection=sys_data.get(
                    "ligand_selection",
                    stage_data.get("ligand_selection", "resname UNK")
                ),
                residue_ids=(
                    sys_residue_cfg["residue_ids"]
                    if sys_residue_cfg.get("residue_ids") is not None
                    else residue_cfg["residue_ids"]
                ),
                residue_selection=(
                    sys_residue_cfg["residue_selection"]
                    if sys_residue_cfg.get("residue_selection") is not None
                    else residue_cfg["residue_selection"]
                ),
                auto_detect_residues=auto_detect_residues,
                auto_detect_cutoff=auto_detect_cutoff,
                auto_detect_stride=auto_detect_stride,
                auto_detect_max_frames=auto_detect_max_frames,
                auto_detect_selection=auto_detect_selection,
                auto_detect_method=auto_detect_method,
                method=sys_data.get("method", stage_data.get("method", "com")),
                distance_scale=sys_data.get("distance_scale", stage_data.get("distance_scale", 1.0)),
                time_unit=sys_data.get("time_unit", stage_data.get("time_unit", "auto")),
                time_step_ns=sys_data.get("time_step_ns", stage_data.get("time_step_ns")),
                per_residue_plots=sys_data.get("per_residue_plots", stage_data.get("per_residue_plots", True)),
                combined_plot=sys_data.get("combined_plot", stage_data.get("combined_plot", True)),
                combined_show_std=sys_data.get("combined_show_std", stage_data.get("combined_show_std", False)),
                include_overall=sys_data.get("include_overall", stage_data.get("include_overall", False)),
                plot_config=plot_config,
            )
        )

    return configs


def _build_distance_compare_config(
    stage_data: Dict[str, Any],
    base_path: Path,
    defaults: Dict[str, Any],
    distance_batch_configs: Optional[List[DistanceConfig]],
) -> DistanceCompareConfig:
    from_batch = bool(stage_data.get("from_distance_batch", False))

    if from_batch:
        if not distance_batch_configs:
            raise ValueError("distance_compare from_distance_batch requires distance_batch stage")
        systems = distance_batch_configs
    elif stage_data.get("systems"):
        systems = _build_distance_batch_configs(stage_data, base_path, defaults)
    else:
        raise ValueError("distance_compare requires from_distance_batch or systems list")

    base_dir = stage_data.get("base_dir") or defaults.get("data_root")
    if base_dir:
        base_dir = _resolve_path(base_dir, base_path)
    else:
        base_dir = systems[0].base_dir

    output_dir = stage_data.get("output_dir") or _default_output_dir(
        "distance_compare",
        defaults.get("output_root"),
        base_dir,
    )
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None
    if output_dir is None:
        raise ValueError("distance_compare requires output_dir or output_root")

    protein_name = stage_data.get("protein_name") or defaults.get("protein_name") or systems[0].protein_name

    plot_data = stage_data.get("plot_config", {})
    base_plot = systems[0].plot_config
    plot_config = PlotConfig(
        template=plot_data.get("template", base_plot.template),
        style=plot_data.get("style", base_plot.style),
        width=plot_data.get("width", base_plot.width),
        height=plot_data.get("height", base_plot.height),
        scale=plot_data.get("scale", base_plot.scale),
        font_family=plot_data.get("font_family", base_plot.font_family),
        font_size=plot_data.get("font_size", base_plot.font_size),
        colors={k: tuple(v) for k, v in plot_data.get("colors", {}).items()},
        save_formats=plot_data.get("save_formats", base_plot.save_formats),
    )

    return DistanceCompareConfig(
        output_dir=output_dir,
        protein_name=protein_name,
        systems=systems,
        compare_mode=stage_data.get("compare_mode", "intersection"),
        per_residue_plots=bool(stage_data.get("per_residue_plots", True)),
        show_std=bool(stage_data.get("show_std", False)),
        include_overall=bool(stage_data.get("include_overall", False)),
        heatmap_enabled=bool(stage_data.get("heatmap_enabled", False)),
        heatmap_downsample=stage_data.get("heatmap_downsample", 1),
        heatmap_max_frames=stage_data.get("heatmap_max_frames"),
        delta_enabled=bool(stage_data.get("delta_enabled", False)),
        delta_reference=stage_data.get("delta_reference"),
        plot_config=plot_config,
    )

def _build_pca_batch_configs(stage_data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]) -> List[PCAConfig]:
    base_dir = stage_data.get("base_dir") or defaults.get("data_root")
    if not base_dir:
        raise ValueError("pca_batch requires base_dir (or paths.data_root)")
    base_dir = _resolve_path(base_dir, base_path)

    output_dir = stage_data.get("output_dir") or _default_output_dir("pca", defaults.get("output_root"), base_dir)
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None

    protein_name = stage_data.get("protein_name") or defaults.get("protein_name") or "Protein"
    n_pcs = stage_data.get("n_pcs", 10)

    systems_data = stage_data.get("systems", [])
    if not systems_data:
        raise ValueError("pca_batch requires systems list")

    shared = {
        "protein_name": protein_name,
        "n_pcs": n_pcs,
        "plot_config": stage_data.get("plot_config", {}),
        "run_gromacs": stage_data.get("run_gromacs", {}),
    }

    output_pattern = stage_data.get("output_subdir_pattern", "{name}")
    configs: List[PCAConfig] = []

    for slot, sys_data in enumerate(systems_data):
        sys_name = sys_data.get("name")
        if not sys_name:
            raise ValueError("pca_batch systems entries require name")

        dir_name = sys_data.get("pca_dir") or sys_data.get("dir_name") or sys_data.get("dir_pattern")
        if not dir_name:
            raise ValueError(f"pca_batch system '{sys_name}' requires pca_dir, dir_name, or dir_pattern")

        sys_base_dir = (base_dir / dir_name).resolve()
        sys_output = output_dir / output_pattern.format(name=sys_name) if output_dir else None
        run_gromacs_data = dict(shared.get("run_gromacs") or {})
        run_gromacs_data.update(sys_data.get("run_gromacs", {}) or {})
        simulation_dir = sys_data.get("simulation_dir") or run_gromacs_data.get("input_dir")
        if simulation_dir:
            simulation_path = Path(str(simulation_dir))
            if not simulation_path.is_absolute():
                simulation_path = (base_dir / simulation_path).resolve()
            run_gromacs_data["input_dir"] = str(simulation_path)

        batch_data = {
            **shared,
            "base_dir": str(sys_base_dir),
            "output_dir": str(sys_output) if sys_output else None,
            "ligand_name": sys_data.get("ligand_name", sys_name),
            "apo_base_dir": sys_data.get("apo_base_dir"),
            "systems": sys_data.get("systems", []),
            "run_gromacs": run_gromacs_data,
            "color": sys_data.get("color") or categorical(slot),
        }

        configs.append(_build_pca_config(batch_data, base_path, defaults))

    return configs


def _build_interaction_config(data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]) -> InteractionConfig:
    raw_base_dir = data.get("base_dir") or defaults.get("data_root")
    if not raw_base_dir:
        raise ValueError("interactions config requires base_dir (or paths.data_root)")
    base_dir = _resolve_path(raw_base_dir, base_path)

    output_dir = data.get("output_dir") or _default_output_dir("interactions", defaults.get("output_root"), base_dir)
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None

    protein_name = data.get("protein_name") or defaults.get("protein_name") or "Protein"
    ligand_name = data.get("ligand_name") or "Ligand"

    topology = data.get("topology")
    if not topology:
        raise ValueError("interactions config requires topology")
    topology = _resolve_stage_input_path(topology, base_path, base_dir, raw_base_dir)

    trajectories = data.get("trajectories", [])
    if not trajectories:
        raise ValueError("interactions config requires trajectories")
    trajectories = [
        _resolve_stage_input_path(p, base_path, base_dir, raw_base_dir)
        for p in trajectories
    ]

    plot_data = data.get("plot_config", {})
    plot_config = PlotConfig(
        template=plot_data.get("template", "ggplot2"),
        style=plot_data.get("style", "simple"),
        width=plot_data.get("width", 1200),
        height=plot_data.get("height", 600),
        scale=plot_data.get("scale", 2),
        font_family=plot_data.get("font_family", "Times New Roman"),
        font_size=plot_data.get("font_size", 20),
        colors={k: tuple(v) for k, v in plot_data.get("colors", {}).items()},
        save_formats=plot_data.get("save_formats", ["html", "svg"]),
    )

    return InteractionConfig(
        base_dir=base_dir,
        output_dir=output_dir,
        protein_name=protein_name,
        ligand_name=ligand_name,
        topology=topology,
        trajectories=trajectories,
        ligand_selection=data.get("ligand_selection", "resname UNK"),
        protein_selection=data.get("protein_selection", "protein"),
        interaction_types=data.get("interaction_types", ["hbond", "salt_bridge", "hydrophobic", "pi_stacking"]),
        hbond_cutoff=data.get("hbond_cutoff", 3.5),
        salt_bridge_cutoff=data.get("salt_bridge_cutoff", 4.0),
        hydrophobic_cutoff=data.get("hydrophobic_cutoff", 4.5),
        pi_stack_cutoff=data.get("pi_stack_cutoff", 5.0),
        stride=data.get("stride", 1),
        ligand_pos_sel=data.get("ligand_pos_sel"),
        ligand_neg_sel=data.get("ligand_neg_sel"),
        ligand_ring_sel=data.get("ligand_ring_sel"),
        plot_config=plot_config,
    )


def _build_interaction_batch_configs(stage_data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]) -> List[InteractionConfig]:
    systems_data = stage_data.get("systems", [])
    if not systems_data:
        raise ValueError("interactions_batch requires systems list")

    raw_base_dir = stage_data.get("base_dir") or defaults.get("data_root")
    if not raw_base_dir:
        raise ValueError("interactions_batch requires base_dir (or paths.data_root)")
    base_dir = _resolve_path(raw_base_dir, base_path)

    output_dir = stage_data.get("output_dir") or _default_output_dir("interactions_batch", defaults.get("output_root"), base_dir)
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None
    output_pattern = stage_data.get("output_subdir_pattern", "output_{name}")

    configs: List[InteractionConfig] = []
    for sys_data in systems_data:
        ligand_name = sys_data.get("ligand_name") or sys_data.get("name") or "Ligand"
        batch_output = output_dir / output_pattern.format(name=ligand_name) if output_dir else None

        cfg_data = {
            **stage_data,
            **sys_data,
            "base_dir": str(sys_data.get("base_dir") or raw_base_dir),
            "output_dir": str(batch_output) if batch_output else None,
            "ligand_name": ligand_name,
        }
        configs.append(_build_interaction_config(cfg_data, base_path, defaults))

    return configs


def _build_interaction_compare_config(
    stage_data: Dict[str, Any],
    base_path: Path,
    defaults: Dict[str, Any],
    interaction_batch_configs: Optional[List[InteractionConfig]],
) -> InteractionCompareConfig:
    from_batch = bool(stage_data.get("from_interactions_batch", False))

    if from_batch:
        if not interaction_batch_configs:
            raise ValueError("interactions_compare from_interactions_batch requires interactions_batch stage")
        systems = interaction_batch_configs
    elif stage_data.get("systems"):
        systems = _build_interaction_batch_configs(stage_data, base_path, defaults)
    else:
        raise ValueError("interactions_compare requires from_interactions_batch or systems list")

    base_dir = stage_data.get("base_dir") or defaults.get("data_root")
    if base_dir:
        base_dir = _resolve_path(base_dir, base_path)
    else:
        base_dir = systems[0].base_dir

    output_dir = stage_data.get("output_dir") or _default_output_dir(
        "interactions_compare",
        defaults.get("output_root"),
        base_dir,
    )
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None
    if output_dir is None:
        raise ValueError("interactions_compare requires output_dir or output_root")

    protein_name = stage_data.get("protein_name") or defaults.get("protein_name") or systems[0].protein_name

    plot_data = stage_data.get("plot_config", {})
    base_plot = systems[0].plot_config
    plot_config = PlotConfig(
        template=plot_data.get("template", base_plot.template),
        style=plot_data.get("style", base_plot.style),
        width=plot_data.get("width", base_plot.width),
        height=plot_data.get("height", base_plot.height),
        scale=plot_data.get("scale", base_plot.scale),
        font_family=plot_data.get("font_family", base_plot.font_family),
        font_size=plot_data.get("font_size", base_plot.font_size),
        colors={k: tuple(v) for k, v in plot_data.get("colors", {}).items()},
        save_formats=plot_data.get("save_formats", base_plot.save_formats),
    )

    return InteractionCompareConfig(
        output_dir=output_dir,
        protein_name=protein_name,
        systems=systems,
        compare_mode=stage_data.get("compare_mode", "intersection"),
        plot_config=plot_config,
    )


def _build_ranking_config(
    stage_data: Dict[str, Any],
    base_path: Path,
    defaults: Dict[str, Any],
    reference_systems: Optional[List[str]] = None,
) -> RankingConfig:
    output_dir = stage_data.get("output_dir") or _default_output_dir(
        "ranking",
        defaults.get("output_root"),
        defaults.get("data_root") or base_path,
    )
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None
    if output_dir is None:
        raise ValueError("ranking requires output_dir or output_root")

    protein_name = stage_data.get("protein_name") or defaults.get("protein_name") or "Protein"

    md_stats_file = stage_data.get("md_stats_file")
    if md_stats_file:
        md_stats_file = _resolve_path(md_stats_file, base_path)

    mmpbsa_data_dir = stage_data.get("mmpbsa_data_dir")
    if mmpbsa_data_dir:
        mmpbsa_data_dir = _resolve_path(mmpbsa_data_dir, base_path)

    distance_summary_file = stage_data.get("distance_summary_file")
    if distance_summary_file:
        distance_summary_file = _resolve_path(distance_summary_file, base_path)

    metric_directions = stage_data.get("metric_directions")
    kwargs = {
        "output_dir": output_dir,
        "protein_name": protein_name,
        "md_stats_file": md_stats_file,
        "mmpbsa_data_dir": mmpbsa_data_dir,
        "distance_summary_file": distance_summary_file,
        "exclude_systems": list(stage_data.get("exclude_systems") or []) + list(reference_systems or []),
        "n_bootstrap": int(stage_data.get("n_bootstrap", 2000)),
        "seed": int(stage_data.get("seed", 0)),
        "n_permutations": int((stage_data.get("tests") or {}).get("n_permutations", 10000)),
    }
    if metric_directions is not None:
        kwargs["metric_directions"] = metric_directions

    return RankingConfig(**kwargs)

def _build_pca_compare_config(
    stage_data: Dict[str, Any],
    base_path: Path,
    defaults: Dict[str, Any],
    pca_batch_configs: Optional[List[PCAConfig]],
) -> PCACompareConfig:
    from_batch = bool(stage_data.get("from_pca_batch", False))

    if from_batch:
        if not pca_batch_configs:
            raise ValueError("pca_compare from_pca_batch requires pca_batch stage")
        systems = pca_batch_configs
    elif stage_data.get("systems"):
        systems = _build_pca_batch_configs(stage_data, base_path, defaults)
    else:
        raise ValueError("pca_compare requires from_pca_batch or systems list")

    base_dir = stage_data.get("base_dir") or defaults.get("data_root")
    if base_dir:
        base_dir = _resolve_path(base_dir, base_path)
    else:
        base_dir = systems[0].base_dir

    output_dir = stage_data.get("output_dir") or _default_output_dir(
        "pca_compare",
        defaults.get("output_root"),
        base_dir,
    )
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None
    if output_dir is None:
        raise ValueError("pca_compare requires output_dir or output_root")

    protein_name = stage_data.get("protein_name") or defaults.get("protein_name") or systems[0].protein_name

    plot_data = stage_data.get("plot_config", {})
    base_plot = systems[0].plot_config
    plot_config = PlotConfig(
        template=plot_data.get("template", base_plot.template),
        style=plot_data.get("style", base_plot.style),
        width=plot_data.get("width", base_plot.width),
        height=plot_data.get("height", base_plot.height),
        scale=plot_data.get("scale", base_plot.scale),
        font_family=plot_data.get("font_family", base_plot.font_family),
        font_size=plot_data.get("font_size", base_plot.font_size),
        colors={k: tuple(v) for k, v in plot_data.get("colors", {}).items()},
        save_formats=plot_data.get("save_formats", base_plot.save_formats),
    )

    return PCACompareConfig(
        output_dir=output_dir,
        protein_name=protein_name,
        systems=systems,
        compare_pairs=stage_data.get("compare_pairs", ["12"]),
        compare_scree=bool(stage_data.get("compare_scree", True)),
        compare_cosine=bool(stage_data.get("compare_cosine", False)),
        downsample=stage_data.get("downsample"),
        max_points=stage_data.get("max_points"),
        plot_config=plot_config,
    )


def _build_md_batch_configs(stage_data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]) -> List[MDConfig]:
    base_dir = stage_data.get("base_dir") or defaults.get("data_root")
    if not base_dir:
        raise ValueError("md_batch requires base_dir (or paths.data_root)")
    base_dir = _resolve_path(base_dir, base_path)

    output_dir = stage_data.get("output_dir") or _default_output_dir("md", defaults.get("output_root"), base_dir)
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None

    protein_name = stage_data.get("protein_name") or defaults.get("protein_name")
    if not protein_name:
        raise ValueError("md_batch requires protein_name (or project.protein_name)")

    common_reps = stage_data.get("replicates")
    apo_system = stage_data.get("apo_system") or stage_data.get("apo")
    holo_systems = stage_data.get("holo_systems", [])
    if not holo_systems:
        raise ValueError("md_batch requires holo_systems list")

    # Share metrics/plot settings across batches
    shared = {
        "base_dir": str(base_dir),
        "protein_name": protein_name,
        "metrics": stage_data.get("metrics", []),
        "plot_config": stage_data.get("plot_config", {}),
        "residue_range": stage_data.get("residue_range"),
        "amino_acids": stage_data.get("amino_acids"),
        "time_unit_input": stage_data.get("time_unit_input", "auto"),
        "time_unit_output": stage_data.get("time_unit_output", "ns"),
        "time_step_ps": stage_data.get("time_step_ps"),
        "time_scale": stage_data.get("time_scale"),
    }

    if shared.get("residue_range") is None:
        shared.pop("residue_range", None)
    if shared.get("amino_acids") is None:
        shared.pop("amino_acids", None)

    output_pattern = stage_data.get("output_subdir_pattern", "{name}")
    configs: List[MDConfig] = []

    for holo in holo_systems:
        holo_entry = dict(holo)
        holo_entry.setdefault("replicates", common_reps or 1)
        holo_entry.setdefault("is_apo", False)

        systems = [holo_entry]
        if apo_system:
            apo_entry = dict(apo_system)
            apo_entry.setdefault("replicates", common_reps or 1)
            apo_entry.setdefault("is_apo", True)
            systems.append(apo_entry)

        batch_output = output_dir / output_pattern.format(name=holo_entry["name"]) if output_dir else None

        batch_data = {
            **shared,
            "base_dir": str(base_dir),
            "output_dir": str(batch_output) if batch_output else None,
            "systems": systems,
        }

        configs.append(_build_md_config(batch_data, base_path, defaults))

    return configs


def _build_md_compare_config(
    stage_data: Dict[str, Any],
    base_path: Path,
    defaults: Dict[str, Any],
    md_batch_configs: Optional[List[MDConfig]],
) -> MDConfig:
    from_batch = bool(stage_data.get("from_md_batch", False))
    palette_data = stage_data.get("palette")
    palette = [tuple(c) for c in palette_data] if palette_data else None

    if from_batch:
        if not md_batch_configs:
            raise ValueError("md_compare from_md_batch requires md_batch stage")

        base_cfg = md_batch_configs[0]
        include_apo = bool(stage_data.get("include_apo", False))

        base_dir = stage_data.get("base_dir") or str(base_cfg.base_dir)
        base_dir = _resolve_path(base_dir, base_path)

        output_dir = stage_data.get("output_dir") or _default_output_dir("md_compare", defaults.get("output_root"), base_dir)
        output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None

        protein_name = stage_data.get("protein_name") or defaults.get("protein_name") or base_cfg.protein_name or "Protein"

        holo_systems: List[SystemConfig] = []
        apo_system: Optional[SystemConfig] = None
        for cfg in md_batch_configs:
            for sys in cfg.systems:
                if sys.is_apo:
                    if apo_system is None:
                        apo_system = sys
                else:
                    holo_systems.append(sys)

        if not holo_systems:
            raise ValueError("md_compare from_md_batch found no holo systems to compare")

        systems: List[SystemConfig] = []
        seen = set()
        if include_apo and apo_system:
            systems.append(_clone_system(apo_system))
            seen.add(apo_system.name)

        for sys in holo_systems:
            if sys.name in seen:
                continue
            systems.append(_clone_system(sys))
            seen.add(sys.name)

        metrics_data = stage_data.get("metrics", [])
        metrics: List[AnalysisMetric] = []
        if metrics_data:
            for m in metrics_data:
                metrics.append(
                    AnalysisMetric(
                        name=m["name"],
                        file_pattern=m["file_pattern"],
                        title=m["title"],
                        ylabel=m["ylabel"],
                        is_holo_only=m.get("is_holo_only", False),
                        data_format=m.get("data_format", "two_column"),
                    )
                )
        else:
            metrics = base_cfg.metrics

        amino_acids = base_cfg.amino_acids
        if "amino_acids" in stage_data:
            amino_acids = parse_amino_acids_config(stage_data.get("amino_acids"))

        residue_range = base_cfg.residue_range
        if stage_data.get("residue_range"):
            residue_range = tuple(stage_data.get("residue_range"))

        plot_data = stage_data.get("plot_config", {})
        base_plot = base_cfg.plot_config
        colors = dict(base_plot.colors)
        colors.update({k: tuple(v) for k, v in plot_data.get("colors", {}).items()})
        plot_config = PlotConfig(
            template=plot_data.get("template", base_plot.template),
            style=plot_data.get("style", base_plot.style),
            width=plot_data.get("width", base_plot.width),
            height=plot_data.get("height", base_plot.height),
            scale=plot_data.get("scale", base_plot.scale),
            font_family=plot_data.get("font_family", base_plot.font_family),
            font_size=plot_data.get("font_size", base_plot.font_size),
            colors=colors,
            save_formats=plot_data.get("save_formats", base_plot.save_formats),
        )

        _assign_palette_if_needed(systems, palette)

        return MDConfig(
            base_dir=base_dir,
            output_dir=output_dir,
            protein_name=protein_name,
            systems=systems,
            metrics=metrics,
            plot_config=plot_config,
            residue_range=residue_range,
            amino_acids=amino_acids,
            time_unit_input=stage_data.get("time_unit_input", base_cfg.time_unit_input),
            time_unit_output=stage_data.get("time_unit_output", base_cfg.time_unit_output),
            time_step_ps=stage_data.get("time_step_ps", base_cfg.time_step_ps),
            time_scale=stage_data.get("time_scale", base_cfg.time_scale),
            smoothing_window=stage_data.get("smoothing_window", base_cfg.smoothing_window),
            outlier_threshold=stage_data.get("outlier_threshold", base_cfg.outlier_threshold),
            ligands=stage_data.get("ligands", base_cfg.ligands),
            rmsf_terminals=base_cfg.rmsf_terminals,
            convert_units=stage_data.get("convert_units", base_cfg.convert_units),
        )

    config = _build_md_config(stage_data, base_path, defaults)
    _assign_palette_if_needed(config.systems, palette)
    return config


def _reference_system_names(config: PipelineConfig) -> List[str]:
    """APO systems declared by already-built MD stages (excluded from ranking)."""
    names: List[str] = []
    for stage_name in ("md", "md_compare", "md_batch"):
        stage = config.stages.get(stage_name)
        md_configs = stage.config if stage is not None else None
        if md_configs is None:
            continue
        for md_config in md_configs if isinstance(md_configs, list) else [md_configs]:
            names.extend(s.name for s in getattr(md_config, "systems", []) if s.is_apo)
    return sorted(set(names))


def load_pipeline_config(yaml_path: Path) -> PipelineConfig:
    """Load pipeline configuration from YAML."""
    base_path = yaml_path.parent
    with open(yaml_path, "r") as f:
        data = yaml.safe_load(f) or {}

    project = data.get("project", {})
    paths = data.get("paths", {})

    data_root = _resolve_path(paths.get("data_root"), base_path)
    output_root = _resolve_path(paths.get("output_root"), base_path)

    config = PipelineConfig(
        project_name=project.get("name", "GROMACS Pipeline"),
        protein_name=project.get("protein_name"),
        data_root=data_root,
        output_root=output_root,
        run=data.get("run", []),
        stages={},
        source_path=Path(yaml_path).resolve(),
    )

    defaults = {
        "protein_name": config.protein_name,
        "data_root": config.data_root,
        "output_root": config.output_root,
    }

    md_batch_configs: Optional[List[MDConfig]] = None
    pca_batch_configs: Optional[List[PCAConfig]] = None
    distance_batch_configs: Optional[List[DistanceConfig]] = None
    interaction_batch_configs: Optional[List[InteractionConfig]] = None
    prolif_batch_configs: Optional[List[ProlifConfig]] = None
    ttclust_batch_configs: Optional[List[TTClustConfig]] = None

    for stage_name in (
        "qc",
        "md",
        "md_batch",
        "md_compare",
        "mmpbsa",
        "pca",
        "pca_batch",
        "pca_compare",
        "pca_joint",
        "distance",
        "distance_batch",
        "distance_compare",
        "interactions",
        "interactions_batch",
        "interactions_compare",
        "prolif",
        "prolif_batch",
        "prolif_compare",
        "ttclust",
        "ttclust_batch",
        "ranking",
        "utils",
    ):
        stage_data = data.get(stage_name)
        if not stage_data:
            continue

        enabled = stage_data.get("enabled", True)
        stage = PipelineStage(name=stage_name, enabled=enabled)

        config_file = stage_data.get("config_file")
        config_file_supported = ("md", "md_compare", "mmpbsa", "pca", "distance")
        if config_file and stage_name not in config_file_supported:
            raise ValueError(
                f"config_file is not supported for stage '{stage_name}'. "
                "Use inline pipeline YAML for this stage."
            )
        if config_file and stage_name in config_file_supported:
            config_path = _resolve_path(config_file, base_path)
            if stage_name == "md":
                stage.config = load_md_yaml(config_path)
            elif stage_name == "md_compare":
                stage.config = load_md_yaml(config_path)
            elif stage_name == "mmpbsa":
                stage.config = load_mmpbsa_yaml(config_path)
            elif stage_name == "pca":
                stage.config = load_pca_yaml(config_path)
            elif stage_name == "distance":
                stage.config = load_distance_yaml(config_path)

            stage.config_file = config_path
            stage.inline = False
            _apply_relative_paths(stage.config, config_path.parent)
        else:
            if stage_name == "md":
                stage.config = _build_md_config(stage_data, base_path, defaults)
            elif stage_name == "qc":
                stage.config = _build_qc_config(stage_data, base_path, defaults)
            elif stage_name == "md_batch":
                stage.config = _build_md_batch_configs(stage_data, base_path, defaults)
                md_batch_configs = stage.config
            elif stage_name == "md_compare":
                stage.config = _build_md_compare_config(stage_data, base_path, defaults, md_batch_configs)
            elif stage_name == "mmpbsa":
                stage.config = _build_mmpbsa_config(stage_data, base_path, defaults)
            elif stage_name == "pca":
                stage.config = _build_pca_config(stage_data, base_path, defaults)
            elif stage_name == "pca_batch":
                stage.config = _build_pca_batch_configs(stage_data, base_path, defaults)
                pca_batch_configs = stage.config
            elif stage_name == "pca_compare":
                stage.config = _build_pca_compare_config(stage_data, base_path, defaults, pca_batch_configs)
            elif stage_name == "pca_joint":
                stage.config = _build_pca_joint_config(stage_data, base_path, defaults)
            elif stage_name == "distance":
                stage.config = _build_distance_config(stage_data, base_path, defaults)
            elif stage_name == "distance_batch":
                stage.config = _build_distance_batch_configs(stage_data, base_path, defaults)
                distance_batch_configs = stage.config
            elif stage_name == "distance_compare":
                stage.config = _build_distance_compare_config(stage_data, base_path, defaults, distance_batch_configs)
            elif stage_name == "interactions":
                stage.config = _build_interaction_config(stage_data, base_path, defaults)
            elif stage_name == "interactions_batch":
                stage.config = _build_interaction_batch_configs(stage_data, base_path, defaults)
                interaction_batch_configs = stage.config
            elif stage_name == "interactions_compare":
                stage.config = _build_interaction_compare_config(stage_data, base_path, defaults, interaction_batch_configs)
            elif stage_name == "prolif":
                stage.config = _build_prolif_config(stage_data, base_path, defaults)
            elif stage_name == "prolif_batch":
                stage.config = _build_prolif_batch_configs(stage_data, base_path, defaults)
                prolif_batch_configs = stage.config
            elif stage_name == "prolif_compare":
                stage.config = _build_prolif_compare_config(stage_data, base_path, defaults, prolif_batch_configs)
            elif stage_name == "ttclust":
                stage.config = _build_ttclust_config(stage_data, base_path, defaults)
            elif stage_name == "ttclust_batch":
                stage.config = _build_ttclust_batch_configs(stage_data, base_path, defaults)
                ttclust_batch_configs = stage.config
            elif stage_name == "ranking":
                stage.config = _build_ranking_config(
                    stage_data, base_path, defaults, _reference_system_names(config)
                )
            elif stage_name == "utils":
                stage.config = stage_data
            stage.inline = True

        stage.cache_enabled = bool(stage_data.get("cache", False))

        # Optional output_dir override
        if stage.config and "output_dir" in stage_data and stage_name not in (
            "utils",
            "md_batch",
            "pca_batch",
            "distance_batch",
            "interactions_batch",
            "prolif_batch",
            "ttclust_batch",
        ):
            stage.config.output_dir = _resolve_path(stage_data["output_dir"], base_path)

        config.stages[stage_name] = stage

    logger.info("Loaded pipeline configuration")
    return config


def _md_config_to_dict(config: MDConfig) -> Dict[str, Any]:
    return {
        "protein_name": config.protein_name,
        "base_dir": str(config.base_dir),
        "output_dir": str(config.output_dir),
        "systems": [
            {
                "name": s.name,
                "dir_pattern": s.dir_pattern,
                "dir_names": s.dir_names,
                "is_apo": s.is_apo,
                "replicates": s.replicates,
                "color": list(s.color) if s.color else None,
            }
            for s in config.systems
        ],
        "metrics": [
            {
                "name": m.name,
                "file_pattern": m.file_pattern,
                "title": m.title,
                "ylabel": m.ylabel,
                "is_holo_only": m.is_holo_only,
                "data_format": m.data_format,
            }
            for m in config.metrics
        ],
        "plot_config": {
            "template": config.plot_config.template,
            "style": config.plot_config.style,
            "width": config.plot_config.width,
            "height": config.plot_config.height,
            "scale": config.plot_config.scale,
            "font_family": config.plot_config.font_family,
            "font_size": config.plot_config.font_size,
            "colors": {k: list(v) for k, v in config.plot_config.colors.items()}
            if config.plot_config.colors
            else None,
            "save_formats": config.plot_config.save_formats,
        },
        "residue_range": list(config.residue_range) if config.residue_range else None,
        "amino_acids": {"custom": config.amino_acids} if config.amino_acids else None,
        "stability": {
            "enabled": config.stability.enabled,
            "metric": config.stability.metric,
            "window": config.stability.window,
            "std_threshold": config.stability.std_threshold,
            "slope_threshold": config.stability.slope_threshold,
            "min_points": config.stability.min_points,
            "apply_to_all_metrics": config.stability.apply_to_all_metrics,
            "method": config.stability.method,
            "n_candidates": config.stability.n_candidates,
        },
        "sequence_topology": str(config.sequence_topology) if config.sequence_topology else None,
        "sequence_selection": config.sequence_selection,
        "time_unit_input": config.time_unit_input,
        "time_unit_output": config.time_unit_output,
        "time_step_ps": config.time_step_ps,
        "time_scale": config.time_scale,
        "smoothing_window": config.smoothing_window,
        "outlier_threshold": config.outlier_threshold,
    }


def _mmpbsa_config_to_dict(config: MMPBSAConfig) -> Dict[str, Any]:
    return {
        "protein_name": config.protein_name,
        "ligand_name": config.ligand_name,
        "base_dir": str(config.base_dir),
        "output_dir": str(config.output_dir),
        "compare_systems": config.compare_systems,
        "compare_binding_key": config.compare_binding_key,
        "compare_components": config.compare_components,
        "compare_component_order": config.compare_component_order,
        "incomplete_replicates_policy": config.incomplete_replicates_policy,
        "require_decomp_replicates": config.require_decomp_replicates,
        "systems": [
            {
                "name": s.name,
                "dir_pattern": s.dir_pattern,
                "replicates": s.replicates,
                "results_file_pattern": s.results_file_pattern,
                "decomp_file_pattern": s.decomp_file_pattern,
            }
            for s in config.systems
        ],
        "file_format": config.file_format,
        "plot_config": {
            "template": config.plot_config.template,
            "style": config.plot_config.style,
            "width": config.plot_config.width,
            "height": config.plot_config.height,
            "scale": config.plot_config.scale,
            "font_family": config.plot_config.font_family,
            "font_size": config.plot_config.font_size,
            "save_formats": config.plot_config.save_formats,
        },
        "amino_acids": {"custom": config.amino_acids} if config.amino_acids else None,
    }


def _qc_config_to_dict(config: QCConfig) -> Dict[str, Any]:
    return {
        "base_dir": str(config.base_dir),
        "output_dir": str(config.output_dir),
        "topology": str(config.topology) if config.topology else None,
        "trajectories": [str(p) for p in config.trajectories] if config.trajectories else [],
        "index_file": str(config.index_file) if config.index_file else None,
        "expected_groups": config.expected_groups or None,
        "protein_selection": config.protein_selection,
        "ligand_selection": config.ligand_selection,
        "expected_replicates": config.expected_replicates,
        "min_distance_cutoff": config.min_distance_cutoff,
        "dissociation_cutoff": config.dissociation_cutoff,
        "check_frame_times": config.check_frame_times,
        "sample_stride": config.sample_stride,
        "max_frames": config.max_frames,
        "fail_on_error": config.fail_on_error,
    }


def _pca_config_to_dict(config: PCAConfig) -> Dict[str, Any]:
    return {
        "base_dir": str(config.base_dir),
        "output_dir": str(config.output_dir) if config.output_dir else None,
        "protein_name": config.protein_name,
        "ligand_name": config.ligand_name,
        "apo_base_dir": str(config.apo_base_dir) if config.apo_base_dir else None,
        "n_pcs": config.n_pcs,
        "fel": {
            "enabled": config.fel.enabled,
            "pairs": config.fel.pairs,
            "contours": config.fel.contours,
            "step": config.fel.step,
            "colorscale": config.fel.colorscale,
            "overlay_projection": config.fel.overlay_projection,
            "include_probability": config.fel.include_probability,
            "include_entropy": config.fel.include_entropy,
            "include_enthalpy": config.fel.include_enthalpy,
        },
        "clustering": {
            "enabled": config.clustering.enabled,
            "method": config.clustering.method,
            "pair": config.clustering.pair,
            "n_clusters": config.clustering.n_clusters,
            "eps": config.clustering.eps,
            "min_samples": config.clustering.min_samples,
            "downsample": config.clustering.downsample,
            "max_points": config.clustering.max_points,
        },
        "run_gromacs": {
            "enabled": config.run_gromacs.enabled,
            "input_dir": str(config.run_gromacs.input_dir) if config.run_gromacs.input_dir else None,
            "gmx": config.run_gromacs.gmx,
            "tpr": config.run_gromacs.tpr,
            "trajectory": config.run_gromacs.trajectory,
            "no_pbc_trajectory": config.run_gromacs.no_pbc_trajectory,
            "index": config.run_gromacs.index,
            "protein_group": config.run_gromacs.protein_group,
            "output_group": config.run_gromacs.output_group,
            "core_group_name": config.run_gromacs.core_group_name,
            "calpha_group_name": config.run_gromacs.calpha_group_name,
            "n_pcs": config.run_gromacs.n_pcs,
            "window_cosine": {
                "enabled": config.run_gromacs.window_cosine_enabled,
                "max_ns": config.run_gromacs.window_max_ns,
                "step_ns": config.run_gromacs.window_step_ns,
            },
            "extreme_nframes": config.run_gromacs.extreme_nframes,
            "overwrite": config.run_gromacs.overwrite,
            "terminal_detection": {
                "enabled": config.run_gromacs.terminal_detection.enabled,
                "rmsf_file": config.run_gromacs.terminal_detection.rmsf_file,
                "pdb_file": config.run_gromacs.terminal_detection.pdb_file,
                "gro_file": config.run_gromacs.terminal_detection.gro_file,
                "smoothing_window": config.run_gromacs.terminal_detection.smoothing_window,
                "mad_multiplier": config.run_gromacs.terminal_detection.mad_multiplier,
                "stable_run": config.run_gromacs.terminal_detection.stable_run,
                "max_terminal_fraction": config.run_gromacs.terminal_detection.max_terminal_fraction,
                "core_residue_range": list(config.run_gromacs.terminal_detection.core_residue_range)
                if config.run_gromacs.terminal_detection.core_residue_range
                else None,
            },
        },
        "plot_config": {
            "template": config.plot_config.template,
            "width": config.plot_config.width,
            "height": config.plot_config.height,
            "scale": config.plot_config.scale,
            "font_family": config.plot_config.font_family,
            "font_size": config.plot_config.font_size,
            "save_formats": config.plot_config.save_formats,
        },
        "systems": [
            {
                "name": s.name,
                "is_apo": s.is_apo,
                "ligand_name": s.ligand_name,
            }
            for s in config.systems
        ]
        if config.systems
        else None,
    }


def _pca_compare_config_to_dict(config: PCACompareConfig) -> Dict[str, Any]:
    return {
        "output_dir": str(config.output_dir),
        "protein_name": config.protein_name,
        "compare_pairs": config.compare_pairs,
        "compare_scree": config.compare_scree,
        "compare_cosine": config.compare_cosine,
        "downsample": config.downsample,
        "max_points": config.max_points,
        "plot_config": {
            "template": config.plot_config.template,
            "style": config.plot_config.style,
            "width": config.plot_config.width,
            "height": config.plot_config.height,
            "scale": config.plot_config.scale,
            "font_family": config.plot_config.font_family,
            "font_size": config.plot_config.font_size,
            "save_formats": config.plot_config.save_formats,
        },
        "systems": [
            _pca_config_to_dict(cfg)
            for cfg in config.systems
        ],
    }


def _distance_config_to_dict(config: DistanceConfig) -> Dict[str, Any]:
    return {
        "base_dir": str(config.base_dir),
        "output_dir": str(config.output_dir),
        "protein_name": config.protein_name,
        "ligand_name": config.ligand_name,
        "topology": str(config.topology),
        "trajectories": [str(p) for p in config.trajectories],
        "ligand_selection": config.ligand_selection,
        "residue_ids": config.residue_ids,
        "residue_selection": config.residue_selection,
        "auto_detect_residues": config.auto_detect_residues,
        "auto_detect_cutoff": config.auto_detect_cutoff,
        "auto_detect_stride": config.auto_detect_stride,
        "auto_detect_max_frames": config.auto_detect_max_frames,
        "auto_detect_selection": config.auto_detect_selection,
        "auto_detect_method": config.auto_detect_method,
        "method": config.method,
        "distance_scale": config.distance_scale,
        "time_unit": config.time_unit,
        "time_step_ns": config.time_step_ns,
        "per_residue_plots": config.per_residue_plots,
        "combined_plot": config.combined_plot,
        "combined_show_std": config.combined_show_std,
        "include_overall": config.include_overall,
        "plot_config": {
            "template": config.plot_config.template,
            "style": config.plot_config.style,
            "width": config.plot_config.width,
            "height": config.plot_config.height,
            "scale": config.plot_config.scale,
            "font_family": config.plot_config.font_family,
            "font_size": config.plot_config.font_size,
            "save_formats": config.plot_config.save_formats,
        },
    }


def _distance_compare_config_to_dict(config: DistanceCompareConfig) -> Dict[str, Any]:
    return {
        "output_dir": str(config.output_dir),
        "protein_name": config.protein_name,
        "compare_mode": config.compare_mode,
        "show_std": config.show_std,
        "include_overall": config.include_overall,
        "heatmap_enabled": config.heatmap_enabled,
        "heatmap_downsample": config.heatmap_downsample,
        "heatmap_max_frames": config.heatmap_max_frames,
        "delta_enabled": config.delta_enabled,
        "delta_reference": config.delta_reference,
        "plot_config": {
            "template": config.plot_config.template,
            "style": config.plot_config.style,
            "width": config.plot_config.width,
            "height": config.plot_config.height,
            "scale": config.plot_config.scale,
            "font_family": config.plot_config.font_family,
            "font_size": config.plot_config.font_size,
            "save_formats": config.plot_config.save_formats,
        },
        "systems": [
            _distance_config_to_dict(cfg)
            for cfg in config.systems
        ],
    }


def _interaction_config_to_dict(config: InteractionConfig) -> Dict[str, Any]:
    return {
        "base_dir": str(config.base_dir),
        "output_dir": str(config.output_dir),
        "protein_name": config.protein_name,
        "ligand_name": config.ligand_name,
        "topology": str(config.topology),
        "trajectories": [str(p) for p in config.trajectories],
        "ligand_selection": config.ligand_selection,
        "protein_selection": config.protein_selection,
        "interaction_types": config.interaction_types,
        "hbond_cutoff": config.hbond_cutoff,
        "salt_bridge_cutoff": config.salt_bridge_cutoff,
        "hydrophobic_cutoff": config.hydrophobic_cutoff,
        "pi_stack_cutoff": config.pi_stack_cutoff,
        "stride": config.stride,
        "ligand_pos_sel": config.ligand_pos_sel,
        "ligand_neg_sel": config.ligand_neg_sel,
        "ligand_ring_sel": config.ligand_ring_sel,
        "plot_config": {
            "template": config.plot_config.template,
            "style": config.plot_config.style,
            "width": config.plot_config.width,
            "height": config.plot_config.height,
            "scale": config.plot_config.scale,
            "font_family": config.plot_config.font_family,
            "font_size": config.plot_config.font_size,
            "save_formats": config.plot_config.save_formats,
        },
    }


def _interaction_compare_config_to_dict(config: InteractionCompareConfig) -> Dict[str, Any]:
    return {
        "output_dir": str(config.output_dir),
        "protein_name": config.protein_name,
        "compare_mode": config.compare_mode,
        "plot_config": {
            "template": config.plot_config.template,
            "style": config.plot_config.style,
            "width": config.plot_config.width,
            "height": config.plot_config.height,
            "scale": config.plot_config.scale,
            "font_family": config.plot_config.font_family,
            "font_size": config.plot_config.font_size,
            "save_formats": config.plot_config.save_formats,
        },
        "systems": [
            _interaction_config_to_dict(cfg)
            for cfg in config.systems
        ],
    }


def _ttclust_config_to_dict(config: TTClustConfig) -> Dict[str, Any]:
    return {
        "base_dir": str(config.base_dir),
        "output_dir": str(config.output_dir),
        "protein_name": config.protein_name,
        "system_name": config.system_name,
        "topology": str(config.topology) if config.topology else None,
        "trajectories": [str(p) for p in config.trajectories],
        "executable": config.executable,
        "stride": config.stride,
        "logfile": config.logfile,
        "select_traj": config.select_traj,
        "select_alignment": config.select_alignment,
        "select_rmsd": config.select_rmsd,
        "method": config.method,
        "cutoff": config.cutoff,
        "n_groups": config.n_groups,
        "autoclust": config.autoclust,
        "interactive_matrix": config.interactive_matrix,
        "axis": config.axis,
        "limit_matrix": config.limit_matrix,
        "extra_args": config.extra_args,
        "environment": config.environment,
    }


def save_pipeline_config(config: PipelineConfig, output_path: Path) -> None:
    """Save pipeline configuration to YAML."""
    data: Dict[str, Any] = {
        "project": {
            "name": config.project_name,
            "protein_name": config.protein_name,
        },
        "paths": {
            "data_root": str(config.data_root) if config.data_root else None,
            "output_root": str(config.output_root) if config.output_root else None,
        },
        "run": config.run or list(config.stages.keys()),
    }

    for stage_name, stage in config.stages.items():
        stage_block: Dict[str, Any] = {"enabled": stage.enabled}
        if stage.cache_enabled:
            stage_block["cache"] = True
        if stage.config_file and not stage.inline:
            stage_block["config_file"] = str(stage.config_file)
            if stage.config and getattr(stage.config, "output_dir", None):
                stage_block["output_dir"] = str(stage.config.output_dir)
        else:
            if stage_name in ("md", "md_compare") and isinstance(stage.config, MDConfig):
                stage_block.update(_md_config_to_dict(stage.config))
            elif stage_name == "qc" and isinstance(stage.config, QCConfig):
                stage_block.update(_qc_config_to_dict(stage.config))
            elif stage_name == "mmpbsa" and isinstance(stage.config, MMPBSAConfig):
                stage_block.update(_mmpbsa_config_to_dict(stage.config))
            elif stage_name == "pca" and isinstance(stage.config, PCAConfig):
                stage_block.update(_pca_config_to_dict(stage.config))
            elif stage_name == "pca_compare" and isinstance(stage.config, PCACompareConfig):
                stage_block.update(_pca_compare_config_to_dict(stage.config))
            elif stage_name == "distance" and isinstance(stage.config, DistanceConfig):
                stage_block.update(_distance_config_to_dict(stage.config))
            elif stage_name == "distance_compare" and isinstance(stage.config, DistanceCompareConfig):
                stage_block.update(_distance_compare_config_to_dict(stage.config))
            elif stage_name == "interactions" and isinstance(stage.config, InteractionConfig):
                stage_block.update(_interaction_config_to_dict(stage.config))
            elif stage_name == "interactions_compare" and isinstance(stage.config, InteractionCompareConfig):
                stage_block.update(_interaction_compare_config_to_dict(stage.config))
            elif stage_name == "ttclust" and isinstance(stage.config, TTClustConfig):
                stage_block.update(_ttclust_config_to_dict(stage.config))
            elif stage_name == "ttclust_batch" and isinstance(stage.config, list):
                if stage.config:
                    first_cfg: TTClustConfig = stage.config[0]
                    stage_block.update({
                        "base_dir": str(first_cfg.base_dir),
                        "output_dir": str(first_cfg.output_dir.parent),
                        "stride": first_cfg.stride,
                        "logfile": first_cfg.logfile,
                        "select_traj": first_cfg.select_traj,
                        "select_alignment": first_cfg.select_alignment,
                        "select_rmsd": first_cfg.select_rmsd,
                        "method": first_cfg.method,
                        "cutoff": first_cfg.cutoff,
                        "n_groups": first_cfg.n_groups,
                        "autoclust": first_cfg.autoclust,
                        "interactive_matrix": first_cfg.interactive_matrix,
                        "axis": first_cfg.axis,
                        "limit_matrix": first_cfg.limit_matrix,
                        "extra_args": first_cfg.extra_args,
                    })
                    stage_block["systems"] = [
                        {
                            "name": cfg.system_name,
                            "system_name": cfg.system_name,
                            "topology": str(cfg.topology) if cfg.topology else None,
                            "trajectories": [str(p) for p in cfg.trajectories],
                        }
                        for cfg in stage.config
                    ]
            elif stage_name == "utils" and isinstance(stage.config, dict):
                stage_block.update(stage.config)
        data[stage_name] = stage_block

    # Clean out None values for nicer YAML
    def _clean(obj: Any) -> Any:
        if isinstance(obj, dict):
            return {k: _clean(v) for k, v in obj.items() if v is not None}
        if isinstance(obj, list):
            return [ _clean(v) for v in obj if v is not None ]
        return obj

    cleaned = _clean(data)
    output_path.write_text(yaml.dump(cleaned, default_flow_style=False, sort_keys=False))
    logger.info(f"Pipeline configuration saved to: {output_path}")


def generate_pipeline_template(output_path: Path) -> None:
    """Generate a pipeline YAML template."""
    template = PIPELINE_TEMPLATE
    output_path.write_text(template)
    logger.info(f"Pipeline template created: {output_path}")
