"""Builders for the per-ligand stages (ProLIF, TTClust) and joint PCA.

Kept separate from ``yaml_loader`` so each stage's configuration code can be read on
its own; ``yaml_loader`` re-exports these names for backward compatibility.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..md.config import PlotConfig
from ..prolif.config import ProlifConfig, ProlifCompareConfig
from ..ttclust.config import TTClustConfig
from .discovery import discover_systems, normalise_ligands
from .loader_utils import _default_output_dir, _resolve_path, _resolve_stage_input_path

logger = logging.getLogger(__name__)


def _batch_systems(stage_data: Dict[str, Any], base_path: Path, stage: str) -> List[Dict[str, Any]]:
    """Explicit ``systems`` list, or systems found by ``auto_discover`` (one per ligand)."""
    systems = list(stage_data.get("systems") or [])
    auto = stage_data.get("auto_discover")
    if auto:
        root_value = auto.get("root") or stage_data.get("base_dir")
        if not root_value:
            raise ValueError(f"{stage}.auto_discover requires root (or stage base_dir)")
        for key in ("pattern", "topology", "trajectory"):
            if not auto.get(key):
                raise ValueError(f"{stage}.auto_discover requires {key}")
        known = {s.get("name") for s in systems}
        enabled = stage_data.get("enabled", True)
        try:
            discovered = discover_systems(
                _resolve_path(str(root_value), base_path),
                auto["pattern"],
                auto["topology"],
                auto["trajectory"],
                include=auto.get("include"),
                exclude=auto.get("exclude"),
            )
        except ValueError as exc:
            if enabled and not systems:
                raise
            if enabled:
                logger.warning(f"{stage}: {exc}; using the explicit systems list only.")
            discovered = []
        # Explicit entries override discovered ones with the same name.
        systems.extend(s for s in discovered if s["name"] not in known)
    if not systems and stage_data.get("enabled", True):
        raise ValueError(f"{stage} requires a systems list or auto_discover")
    return systems


def _build_prolif_config(data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]) -> ProlifConfig:
    raw_base_dir = data.get("base_dir") or defaults.get("data_root")
    if not raw_base_dir:
        raise ValueError("prolif config requires base_dir (or paths.data_root)")
    base_dir = _resolve_path(raw_base_dir, base_path)

    output_dir = data.get("output_dir") or _default_output_dir("prolif", defaults.get("output_root"), base_dir)
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None

    protein_name = data.get("protein_name") or defaults.get("protein_name") or "Protein"
    ligand_name = data.get("ligand_name") or "Ligand"

    topology = data.get("topology")
    if not topology:
        raise ValueError("prolif config requires topology")
    topology = _resolve_stage_input_path(topology, base_path, base_dir, raw_base_dir)

    trajectories = data.get("trajectories", [])
    if not trajectories:
        raise ValueError("prolif config requires trajectories")
    trajectories = [
        _resolve_stage_input_path(p, base_path, base_dir, raw_base_dir)
        for p in trajectories
    ]

    plot_data = data.get("plot_config", {})
    barcode_cfg = data.get("barcode_config", {}) or {}
    raw_figsize = barcode_cfg.get("figsize", data.get("barcode_figsize", [12.0, 4.0]))
    if not isinstance(raw_figsize, (list, tuple)) or len(raw_figsize) != 2:
        raise ValueError("prolif barcode figsize must be a 2-item list/tuple, e.g. [12, 4]")
    barcode_figsize = (float(raw_figsize[0]), float(raw_figsize[1]))
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

    return ProlifConfig(
        base_dir=base_dir,
        output_dir=output_dir,
        protein_name=protein_name,
        ligand_name=ligand_name,
        topology=topology,
        trajectories=trajectories,
        ligand_selection=data.get("ligand_selection", "resname UNK"),
        protein_selection=data.get("protein_selection", "protein"),
        interaction_types=data.get("interaction_types", None) or [],
        vicinity_cutoff=data.get("vicinity_cutoff", 5.0),
        stride=data.get("stride", 1),
        n_jobs=data.get("n_jobs"),
        save_pickle=bool(data.get("save_pickle", True)),
        outputs=data.get("outputs", ["barcode", "lignetwork", "occupancy", "tanimoto"]),
        lignetwork_display_all=bool(data.get("lignetwork_display_all", True)),
        lignetwork_count_aware=bool(data.get("lignetwork_count_aware", True)),
        lignetwork_kind=data.get("lignetwork_kind", "aggregate"),
        lignetwork_threshold=float(data.get("lignetwork_threshold", 0.3)),
        occupancy_top_n=int(data.get("occupancy_top_n", 25)),
        barcode_n_frame_ticks=int(barcode_cfg.get("n_frame_ticks", data.get("barcode_n_frame_ticks", 10))),
        barcode_residues_tick_location=barcode_cfg.get(
            "residues_tick_location",
            data.get("barcode_residues_tick_location", "top"),
        ),
        barcode_xlabel=barcode_cfg.get("xlabel", data.get("barcode_xlabel", "Frame")),
        barcode_figsize=barcode_figsize,
        barcode_dpi=int(barcode_cfg.get("dpi", data.get("barcode_dpi", 150))),
        barcode_only_interacting_residues=bool(
            barcode_cfg.get("only_interacting_residues", data.get("barcode_only_interacting_residues", True))
        ),
        barcode_min_interaction_frames=int(
            barcode_cfg.get("min_interaction_frames", data.get("barcode_min_interaction_frames", 1))
        ),
        barcode_max_residues=(
            int(barcode_cfg.get("max_residues", data.get("barcode_max_residues")))
            if barcode_cfg.get("max_residues", data.get("barcode_max_residues")) is not None
            else None
        ),
        plot_config=plot_config,
    )


def _build_prolif_batch_configs(stage_data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]) -> List[ProlifConfig]:
    raw_base_dir = stage_data.get("base_dir") or defaults.get("data_root")
    if not raw_base_dir:
        raise ValueError("prolif_batch requires base_dir (or paths.data_root)")
    base_dir = _resolve_path(raw_base_dir, base_path)
    systems_data = _batch_systems({**stage_data, "base_dir": raw_base_dir}, base_path, "prolif_batch")

    output_dir = stage_data.get("output_dir") or _default_output_dir("prolif_batch", defaults.get("output_root"), base_dir)
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None
    output_pattern = stage_data.get("output_subdir_pattern", "output_{name}")

    configs: List[ProlifConfig] = []
    for sys_data in systems_data:
        system_name = sys_data.get("ligand_name") or sys_data.get("name") or "Ligand"
        sys_raw_base = str(sys_data.get("base_dir") or raw_base_dir)
        ligands_spec = sys_data.get("ligands", stage_data.get("ligands"))
        topology = sys_data.get("topology") or stage_data.get("topology")
        topology_path = (
            _resolve_stage_input_path(topology, base_path, _resolve_path(sys_raw_base, base_path), sys_raw_base)
            if topology else None
        )
        ligands = normalise_ligands(
            ligands_spec,
            topology_path,
            exclude=sys_data.get("ligand_exclude", stage_data.get("ligand_exclude", [])),
        )
        # One fingerprint run per ligand; a single ligand keeps the system name.
        runs = ligands or [{"name": None, "selection": None}]
        for lig in runs:
            name = system_name if len(runs) == 1 or not lig["name"] else f"{system_name}_{lig['name']}"
            batch_output = output_dir / output_pattern.format(name=name) if output_dir else None
            cfg_data = {
                **{k: v for k, v in stage_data.items() if k not in ("systems", "auto_discover", "ligands")},
                **{k: v for k, v in sys_data.items() if k not in ("ligands", "directories")},
                "base_dir": sys_raw_base,
                "output_dir": str(batch_output) if batch_output else None,
                "ligand_name": name,
            }
            if lig["selection"]:
                cfg_data["ligand_selection"] = lig["selection"]
            configs.append(_build_prolif_config(cfg_data, base_path, defaults))

    return configs


def _build_prolif_compare_config(
    stage_data: Dict[str, Any],
    base_path: Path,
    defaults: Dict[str, Any],
    prolif_batch_configs: Optional[List[ProlifConfig]],
) -> ProlifCompareConfig:
    from_batch = bool(stage_data.get("from_prolif_batch", False))

    if from_batch:
        if not prolif_batch_configs:
            raise ValueError("prolif_compare from_prolif_batch requires prolif_batch stage")
        systems = prolif_batch_configs
    elif stage_data.get("systems") or stage_data.get("auto_discover"):
        systems = _build_prolif_batch_configs(stage_data, base_path, defaults)
    else:
        raise ValueError("prolif_compare requires from_prolif_batch or systems list")

    base_dir = stage_data.get("base_dir") or defaults.get("data_root")
    if base_dir:
        base_dir = _resolve_path(base_dir, base_path)
    else:
        base_dir = systems[0].base_dir

    output_dir = stage_data.get("output_dir") or _default_output_dir(
        "prolif_compare",
        defaults.get("output_root"),
        base_dir,
    )
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None
    if output_dir is None:
        raise ValueError("prolif_compare requires output_dir or output_root")

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

    return ProlifCompareConfig(
        output_dir=output_dir,
        protein_name=protein_name,
        systems=systems,
        compare_mode=stage_data.get("compare_mode", "intersection"),
        plot_config=plot_config,
    )


def _build_ttclust_config(data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]) -> TTClustConfig:
    raw_base_dir = data.get("base_dir") or defaults.get("data_root")
    if not raw_base_dir:
        raise ValueError("ttclust config requires base_dir (or paths.data_root)")
    base_dir = _resolve_path(raw_base_dir, base_path)

    output_dir = data.get("output_dir") or _default_output_dir("ttclust", defaults.get("output_root"), base_dir)
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None
    if output_dir is None:
        raise ValueError("ttclust config requires output_dir or output_root")

    protein_name = data.get("protein_name") or defaults.get("protein_name") or "Protein"
    system_name = data.get("system_name") or data.get("ligand_name") or data.get("name") or "System"

    topology = data.get("topology")
    if topology:
        topology = _resolve_stage_input_path(topology, base_path, base_dir, raw_base_dir)

    trajectories = data.get("trajectories", [])
    if not trajectories:
        raise ValueError("ttclust config requires trajectories")
    trajectories = [
        _resolve_stage_input_path(p, base_path, base_dir, raw_base_dir)
        for p in trajectories
    ]

    return TTClustConfig(
        base_dir=base_dir,
        output_dir=output_dir,
        protein_name=protein_name,
        system_name=system_name,
        topology=topology,
        trajectories=trajectories,
        executable=data.get("executable"),
        stride=int(data.get("stride", 1)),
        logfile=data.get("logfile", "clustering.log"),
        select_traj=data.get("select_traj", "all"),
        select_alignment=data.get("select_alignment", "backbone"),
        select_rmsd=data.get("select_rmsd", "backbone"),
        method=data.get("method", "ward"),
        cutoff=data.get("cutoff"),
        n_groups=str(data["n_groups"]) if data.get("n_groups") is not None else None,
        autoclust=bool(data.get("autoclust", True)),
        interactive_matrix=bool(data.get("interactive_matrix", False)),
        axis=data.get("axis"),
        limit_matrix=data.get("limit_matrix"),
        extra_args=[str(v) for v in data.get("extra_args", [])],
        environment={str(k): str(v) for k, v in (data.get("environment", {}) or {}).items()},
        plot_config=PlotConfig(
            template=(data.get("plot_config") or {}).get("template", "plotly_white"),
            width=(data.get("plot_config") or {}).get("width", 1200),
            height=(data.get("plot_config") or {}).get("height", 600),
            scale=(data.get("plot_config") or {}).get("scale", 2),
            font_family=(data.get("plot_config") or {}).get("font_family", "Arial"),
            font_size=(data.get("plot_config") or {}).get("font_size", 16),
            save_formats=(data.get("plot_config") or {}).get("save_formats", ["html", "svg"]),
        ),
    )


def _build_ttclust_batch_configs(stage_data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]) -> List[TTClustConfig]:
    raw_base_dir = stage_data.get("base_dir") or defaults.get("data_root")
    if not raw_base_dir:
        raise ValueError("ttclust_batch requires base_dir (or paths.data_root)")
    base_dir = _resolve_path(raw_base_dir, base_path)
    systems_data = _batch_systems({**stage_data, "base_dir": raw_base_dir}, base_path, "ttclust_batch")

    output_dir = stage_data.get("output_dir") or _default_output_dir("ttclust_batch", defaults.get("output_root"), base_dir)
    output_dir = _resolve_path(str(output_dir), base_path) if output_dir else None
    output_pattern = stage_data.get("output_subdir_pattern", "output_{name}")

    configs: List[TTClustConfig] = []
    for sys_data in systems_data:
        sys_name = sys_data.get("name") or sys_data.get("system_name") or sys_data.get("ligand_name")
        if not sys_name:
            raise ValueError("ttclust_batch systems entries require name")
        batch_output = output_dir / output_pattern.format(name=sys_name) if output_dir else None

        cfg_data = {
            **{k: v for k, v in stage_data.items() if k not in ("systems", "auto_discover")},
            **{k: v for k, v in sys_data.items() if k != "directories"},
            "base_dir": str(sys_data.get("base_dir") or raw_base_dir),
            "output_dir": str(batch_output) if batch_output else None,
            "system_name": sys_data.get("system_name") or sys_name,
            "name": sys_name,
        }
        configs.append(_build_ttclust_config(cfg_data, base_path, defaults))

    return configs


def _build_pca_joint_config(stage_data: Dict[str, Any], base_path: Path, defaults: Dict[str, Any]):
    """Joint PCA over several systems (explicit ``systems`` or ``auto_discover``)."""
    from ..pca.config import PCAJointConfig, PCAJointSystem  # local import: pca imports plotting deps

    raw_base_dir = stage_data.get("base_dir") or defaults.get("data_root")
    if not raw_base_dir:
        raise ValueError("pca_joint requires base_dir (or paths.data_root)")
    base_dir = _resolve_path(str(raw_base_dir), base_path)
    output_dir = stage_data.get("output_dir") or _default_output_dir("pca_joint", defaults.get("output_root"), base_dir)
    if not output_dir:
        raise ValueError("pca_joint requires output_dir or paths.output_root")
    output_dir = _resolve_path(str(output_dir), base_path)

    systems_data = _batch_systems({**stage_data, "base_dir": raw_base_dir}, base_path, "pca_joint")
    if not systems_data:
        return None  # disabled stage without systems: nothing to configure
    systems = []
    for sys_data in systems_data:
        sys_raw_base = str(sys_data.get("base_dir") or raw_base_dir)
        sys_base = _resolve_path(sys_raw_base, base_path)
        topology = sys_data.get("topology") or stage_data.get("topology")
        if not topology:
            raise ValueError(f"pca_joint system {sys_data.get('name')} requires topology")
        systems.append(PCAJointSystem(
            name=str(sys_data.get("name")),
            topology=_resolve_stage_input_path(topology, base_path, sys_base, sys_raw_base),
            trajectories=[
                _resolve_stage_input_path(t, base_path, sys_base, sys_raw_base)
                for t in sys_data.get("trajectories", [])
            ],
            is_apo=bool(sys_data.get("is_apo", str(sys_data.get("name", "")).upper() == "APO")),
            color=sys_data.get("color"),
        ))

    plot_data = stage_data.get("plot_config") or {}
    plot_config = PlotConfig(
        template=plot_data.get("template", "plotly_white"),
        width=plot_data.get("width", 1200),
        height=plot_data.get("height", 600),
        scale=plot_data.get("scale", 2),
        font_family=plot_data.get("font_family", "Arial"),
        font_size=plot_data.get("font_size", 16),
        save_formats=plot_data.get("save_formats", ["html", "svg"]),
    )
    return PCAJointConfig(
        output_dir=output_dir,
        protein_name=stage_data.get("protein_name") or defaults.get("protein_name") or "Protein",
        systems=systems,
        base_dir=base_dir,
        selection=stage_data.get("selection", "protein and name CA"),
        n_components=int(stage_data.get("n_components", 10)),
        n_overlap=int(stage_data.get("n_overlap", 10)),
        stride=int(stage_data.get("stride", 1)),
        max_frames_per_system=stage_data.get("max_frames_per_system"),
        fit_iterations=int(stage_data.get("fit_iterations", 2)),
        convergence={k: float(v) for k, v in (stage_data.get("convergence") or {}).items()},
        plot_config=plot_config,
    )
