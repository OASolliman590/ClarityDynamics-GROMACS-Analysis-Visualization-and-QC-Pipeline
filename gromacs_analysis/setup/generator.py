"""Turn a :class:`Project` into a pipeline YAML that the existing loader accepts.

Every system is listed explicitly (no globbing), so display names flow into every
stage while selections keep the original residue names. The result is validated
with ``load_pipeline_config`` before it is written.
"""

from __future__ import annotations

import posixpath
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml

from .project import Project, StageChoice, SystemEntry

# Stages the setup tool can configure fully from discovery, in pipeline order.
STAGES: List[Tuple[str, str]] = [
    ("qc", "Preflight QC (frames, timestep, ligand contact)"),
    ("md", "MD metrics: RMSD/RMSF/Rg/SASA/H-bonds/COM plots and statistics"),
    ("mmpbsa", "MM/PBSA binding energy, decomposition and entropy"),
    ("pca_joint", "Joint PCA on shared axes with subspace overlap"),
    ("distance_batch", "Ligand–residue distances (pocket residues auto-detected)"),
    ("prolif_batch", "ProLIF interaction fingerprints, one run per ligand"),
    ("prolif_compare", "Cross-ligand ProLIF comparison"),
    ("ttclust_batch", "TTClust conformational clustering (slow)"),
    ("ranking", "Ligand ranking with bootstrap uncertainty"),
]
STAGE_LABELS = dict(STAGES)
HEAVY_STAGES = {"ttclust_batch"}


def stage_availability(project: Project) -> Dict[str, StageChoice]:
    """Which stages can run on the discovered inputs, with a reason when not."""
    systems = project.active_systems
    holo = [s for s in systems if not s.apo]
    with_ligands = [s for s in holo if s.active_ligands]
    has_metrics = any(s.inputs.get("md_metrics") for s in systems)
    has_mmpbsa = any(s.inputs.get("mmpbsa") for s in holo)
    ligand_runs = sum(len(s.active_ligands) for s in with_ligands)

    def choice(ok: bool, reason: str) -> StageChoice:
        return StageChoice(enabled=ok, available=ok, reason="" if ok else reason)

    avail = {
        "qc": choice(bool(systems), "no systems"),
        "md": choice(has_metrics, "no RMSD/RMSF/... metric files found (generate them with the VMD helper)"),
        "mmpbsa": choice(has_mmpbsa, "no FINAL_RESULTS_MMPBSA files found"),
        "pca_joint": choice(len(systems) >= 1, "no systems"),
        "distance_batch": choice(bool(with_ligands), "no ligand residues detected"),
        "prolif_batch": choice(bool(with_ligands), "no ligand residues detected"),
        "prolif_compare": choice(ligand_runs >= 2, "needs at least two ligand runs"),
        "ttclust_batch": choice(bool(systems), "no systems"),
        "ranking": choice(has_metrics or has_mmpbsa, "needs MD metrics or MM/PBSA results"),
    }
    for name in HEAVY_STAGES:
        avail[name].enabled = False
    return avail


def refresh_stages(project: Project) -> None:
    """Update availability; keep the user's on/off choice for stages that stay available."""
    fresh = stage_availability(project)
    for name, choice in fresh.items():
        previous = project.stages.get(name)
        if previous is not None and choice.available:
            choice.enabled = previous.enabled
        project.stages[name] = choice


def enabled_stages(project: Project) -> List[str]:
    if not project.stages:
        refresh_stages(project)
    return [name for name, _ in STAGES
            if project.stages.get(name) and project.stages[name].enabled and project.stages[name].available]


# ---------------------------------------------------------------- stage blocks

def _abs(project: Project, rel: Optional[str]) -> Optional[str]:
    if rel is None:
        return None
    return posixpath.normpath(posixpath.join(project.root, rel))


def _output(project: Project, stage: str) -> str:
    return posixpath.join(project.root, project.output_root, stage)


def _dir_pattern(system: SystemEntry) -> Tuple[str, int]:
    """'HOLO_LigA_{}' when replicate dirs are <base><sep>1..N, else the first replicate dir."""
    dirs = [r.dir for r in system.replicates]
    numbers = [r.number for r in system.replicates]
    if len(dirs) > 1 and numbers == list(range(1, len(dirs) + 1)):
        first = dirs[0]
        suffix = str(numbers[0])
        if first.endswith(suffix):
            pattern = first[: -len(suffix)] + "{}"
            if all(pattern.format(n) == d for n, d in zip(numbers, dirs)):
                return pattern, len(dirs)
    return dirs[0], 1


def _qc_block(project: Project) -> dict:
    reference = next((s for s in project.active_systems if not s.apo), project.active_systems[0])
    ligand = reference.active_ligands[0].selection if reference.active_ligands else "resname UNK"
    return {
        "enabled": True,
        "output_dir": _output(project, "qc"),
        "topology": _abs(project, reference.topology),
        "trajectories": [_abs(project, t) for t in reference.trajectories],
        "ligand_selection": ligand,
        "expected_replicates": len(reference.replicates),
    }


def _md_block(project: Project) -> dict:
    systems = []
    multi = sorted({lig.resname for s in project.active_systems for lig in s.active_ligands
                    if len(s.active_ligands) > 1})
    for s in project.active_systems:
        dirs = [r.metric_dir for r in s.replicates if r.metric_dir]
        if not dirs:
            continue
        systems.append({"name": s.name, "dir_names": dirs, "is_apo": s.apo, "replicates": len(dirs)})
    block = {
        "enabled": True,
        "base_dir": project.root,
        "output_dir": _output(project, "md"),
        "protein_name": project.protein_name,
        "time_unit_input": "auto",
        "time_unit_output": "ns",
        "systems": systems,
    }
    if multi:
        block["ligands"] = multi  # metric files are named by residue (RMSD_<RES>.dat)
    return block


def _mmpbsa_block(project: Project) -> dict:
    systems = []
    for s in project.active_systems:
        if s.apo or not s.inputs.get("mmpbsa"):
            continue
        pattern, n = _dir_pattern(s)
        systems.append({"name": s.name, "dir_pattern": pattern, "replicates": n,
                        "results_file_pattern": "FINAL_RESULTS_MMPBSA*.dat",
                        "decomp_file_pattern": "FINAL_DECOMP_MMPBSA*.dat"})
    return {
        "enabled": True,
        "base_dir": project.root,
        "output_dir": _output(project, "mmpbsa"),
        "protein_name": project.protein_name,
        "ligand_name": "ligands",
        "compare_systems": len(systems) > 1,
        "systems": systems,
    }


def _trajectory_systems(project: Project, holo_only: bool = False) -> List[dict]:
    entries = []
    for s in project.active_systems:
        if holo_only and s.apo:
            continue
        entries.append({
            "name": s.name,
            "is_apo": s.apo,
            "topology": _abs(project, s.topology),
            "trajectories": [_abs(project, t) for t in s.trajectories],
        })
    return entries


def _pca_joint_block(project: Project) -> dict:
    return {
        "enabled": True,
        "base_dir": project.root,
        "output_dir": _output(project, "pca_joint"),
        "selection": "protein and name CA",
        "stride": 10,
        "systems": _trajectory_systems(project),
    }


def _ligand_systems(project: Project) -> List[dict]:
    entries = []
    for s in project.active_systems:
        if s.apo or not s.active_ligands:
            continue
        entries.append({
            "name": s.name,
            "topology": _abs(project, s.topology),
            "trajectories": [_abs(project, t) for t in s.trajectories],
            "ligands": [{"name": lig.name, "selection": lig.selection} for lig in s.active_ligands],
        })
    return entries


def _distance_block(project: Project) -> dict:
    systems = []
    for s in _ligand_systems(project):
        for lig in s["ligands"]:
            name = s["name"] if len(s["ligands"]) == 1 else f"{s['name']}_{lig['name']}"
            systems.append({"name": name, "topology": s["topology"], "trajectories": s["trajectories"],
                            "ligand_selection": lig["selection"], "ligand_name": lig["name"]})
    return {
        "enabled": True,
        "base_dir": project.root,
        "output_dir": _output(project, "distance_batch"),
        "residues": {"auto": True, "cutoff": 4.5, "stride": 10},
        "systems": systems,
    }


def _prolif_block(project: Project) -> dict:
    return {
        "enabled": True,
        "base_dir": project.root,
        "output_dir": _output(project, "prolif_batch"),
        "stride": 10,
        "systems": _ligand_systems(project),
    }


def _ttclust_block(project: Project) -> dict:
    return {
        "enabled": True,
        "base_dir": project.root,
        "output_dir": _output(project, "ttclust_batch"),
        "stride": 10,
        "systems": [{k: v for k, v in s.items() if k != "is_apo"} for s in _trajectory_systems(project)],
    }


def _ranking_block(project: Project, stages: List[str]) -> dict:
    block = {"enabled": True, "output_dir": _output(project, "ranking"), "n_bootstrap": 2000, "seed": 0,
             "exclude_systems": [s.name for s in project.active_systems if s.apo]}
    if "md" in stages:
        block["md_stats_file"] = posixpath.join(_output(project, "md"), "statistics_summary.csv")
    if "mmpbsa" in stages:
        block["mmpbsa_data_dir"] = posixpath.join(_output(project, "mmpbsa"), "data")
    return block


def build_pipeline_dict(project: Project) -> dict:
    problems = project.validate()
    if problems:
        raise ValueError("Cannot generate configuration:\n- " + "\n- ".join(problems))
    stages = enabled_stages(project)
    data: dict = {
        "project": {"name": project.project_name, "protein_name": project.protein_name},
        "paths": {"data_root": project.root, "output_root": posixpath.join(project.root, project.output_root)},
        "run": stages,
    }
    builders = {
        "qc": lambda: _qc_block(project),
        "md": lambda: _md_block(project),
        "mmpbsa": lambda: _mmpbsa_block(project),
        "pca_joint": lambda: _pca_joint_block(project),
        "distance_batch": lambda: _distance_block(project),
        "prolif_batch": lambda: _prolif_block(project),
        "prolif_compare": lambda: {"enabled": True, "from_prolif_batch": True,
                                   "output_dir": _output(project, "prolif_compare")},
        "ttclust_batch": lambda: _ttclust_block(project),
        "ranking": lambda: _ranking_block(project, stages),
    }
    for stage in stages:
        data[stage] = builders[stage]()
    if "prolif_compare" in stages and "prolif_batch" not in stages:
        raise ValueError("prolif_compare needs prolif_batch")
    return data


def render_yaml(project: Project, project_file: Optional[str] = None) -> str:
    header = [
        "# Generated by clarity-setup; regenerate instead of editing by hand.",
        f"# Project file: {project_file or '(not saved)'}",
        f"# Data: {project.connection.describe()}:{project.root}",
    ]
    body = yaml.safe_dump(build_pipeline_dict(project), sort_keys=False, allow_unicode=True)
    return "\n".join(header) + "\n" + body


def validate_yaml(text: str) -> None:
    """Load with the real pipeline loader; raises ValueError with the loader's message."""
    from ..pipeline.yaml_loader import load_pipeline_config  # noqa: WPS433 - heavy import

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "pipeline.yaml"
        path.write_text(text)
        try:
            config = load_pipeline_config(path)
        except Exception as exc:  # noqa: BLE001
            raise ValueError(f"Generated configuration failed validation: {exc}") from exc
        missing = [s for s in yaml.safe_load(text).get("run", []) if s not in config.stages
                   or config.stages[s].config is None]
        if missing:
            raise ValueError(f"Stages not configured by the loader: {missing}")


def write_pipeline(project: Project, path: Path, project_file: Optional[str] = None) -> Path:
    text = render_yaml(project, project_file)
    validate_yaml(text)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path
