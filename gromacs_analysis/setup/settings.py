"""Editable analysis settings per stage: schema, validation and application to the YAML.

Only overrides are stored in the project file; each setting's default is the value the
generator or the stage loader uses when nothing is set. Keys may be dotted to reach
nested blocks (``residues.cutoff``, ``convergence.neff_min``, ``tests.n_permutations``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .project import Project, ProjectError


@dataclass(frozen=True)
class Setting:
    stage: str
    key: str
    kind: str                      # int | float | str | bool | choice
    default: Any
    help: str
    minimum: Optional[float] = None
    maximum: Optional[float] = None
    choices: Tuple[str, ...] = ()
    optional: bool = False         # "none" is allowed (e.g. no frame limit)

    @property
    def name(self) -> str:
        return f"{self.stage}.{self.key}"


TTCLUST_METHODS = ("ward", "single", "complete", "average", "weighted", "centroid", "median")

SCHEMA: List[Setting] = [
    Setting("qc", "sample_stride", "int", 10, "Analyse every Nth frame in the preflight checks", 1),
    Setting("qc", "max_frames", "int", 500, "Maximum frames sampled per trajectory", 10, optional=True),
    Setting("qc", "dissociation_cutoff", "float", 5.0, "Warn when the ligand moves further than this (Å)", 0.1,
            optional=True),
    Setting("md", "smoothing_window", "int", 1, "Moving-average window for time-series plots (1 = off)", 1),
    Setting("md", "time_unit_output", "choice", "ns", "Time unit on plot axes", choices=("ns", "ps", "frame", "auto")),
    Setting("md", "convergence.neff_min", "float", 20.0, "Flag replicates with fewer effective samples", 1),
    Setting("md", "convergence.half_z_max", "float", 3.0, "Flag drift when the two halves differ by more (z)", 0.5),
    Setting("mmpbsa", "start_frame", "int", None, "First frame used from per-frame MM/PBSA results", 0,
            optional=True),
    Setting("mmpbsa", "incomplete_replicates_policy", "choice", "warn",
            "What to do when a replicate's results are missing", choices=("ignore", "warn", "error")),
    Setting("pca_joint", "selection", "str", "protein and name CA", "Atoms used for joint PCA (MDAnalysis syntax)"),
    Setting("pca_joint", "n_components", "int", 10, "Principal components reported", 1, 100),
    Setting("pca_joint", "n_overlap", "int", 10, "Components used for RMSIP overlap", 1, 100),
    Setting("pca_joint", "stride", "int", 10, "Read every Nth frame", 1),
    Setting("pca_joint", "max_frames_per_system", "int", None, "Cap frames per system (evenly spaced)", 10,
            optional=True),
    Setting("pca_joint", "convergence.cosine_max", "float", 0.5, "Flag PC1 cosine content above this", 0.0, 1.0),
    Setting("pca_joint", "convergence.half_rmsip_min", "float", 0.5, "Flag half-split RMSIP below this", 0.0, 1.0),
    Setting("distance_batch", "residues.cutoff", "float", 4.5, "Pocket residues: within this distance of the ligand (Å)",
            0.5, 20.0),
    Setting("distance_batch", "residues.stride", "int", 10, "Frames sampled when detecting pocket residues", 1),
    Setting("prolif_batch", "stride", "int", 10, "Fingerprint every Nth frame", 1),
    Setting("prolif_batch", "vicinity_cutoff", "float", 5.0, "Residues considered around the ligand (Å)", 1.0, 20.0),
    Setting("prolif_batch", "lignetwork_threshold", "float", 0.3, "Minimum occupancy shown in the ligand network",
            0.0, 1.0),
    Setting("prolif_batch", "n_jobs", "int", None, "Parallel workers for ProLIF (none = automatic)", 1, optional=True),
    Setting("ttclust_batch", "stride", "int", 10, "Cluster every Nth frame", 1),
    Setting("ttclust_batch", "method", "choice", "ward", "Linkage method", choices=TTCLUST_METHODS),
    Setting("ttclust_batch", "select_alignment", "str", "backbone", "Atoms used for alignment"),
    Setting("ttclust_batch", "select_rmsd", "str", "backbone", "Atoms used for the RMSD matrix"),
    Setting("ranking", "n_bootstrap", "int", 2000, "Bootstrap resamples for rank probabilities (0 = off)", 0, 100000),
    Setting("ranking", "seed", "int", 0, "Random seed for bootstrap and permutation tests", 0),
    Setting("ranking", "tests.n_permutations", "int", 10000, "Monte-Carlo permutations for pairwise tests", 100,
            1000000),
]
_BY_NAME = {s.name: s for s in SCHEMA}


def settings_for(stage: Optional[str] = None) -> List[Setting]:
    return [s for s in SCHEMA if stage is None or s.stage == stage]


def lookup(name: str) -> Setting:
    try:
        return _BY_NAME[name]
    except KeyError:
        stages = sorted({s.stage for s in SCHEMA})
        raise ProjectError(f"Unknown setting {name!r}. Stages with settings: {', '.join(stages)}; "
                           "list them with `clarity-setup settings`") from None


def coerce(setting: Setting, raw: Any) -> Any:
    """Validate and convert a raw value (text from the CLI/interface, or a typed value)."""
    if isinstance(raw, str):
        text = raw.strip()
        if setting.optional and text.lower() in ("", "none", "null", "default"):
            return None
    else:
        text = raw
        if raw is None:
            if setting.optional:
                return None
            raise ProjectError(f"{setting.name} needs a value")
    try:
        if setting.kind == "int":
            if isinstance(text, float) and not text.is_integer():
                raise ValueError
            value: Any = int(text)
        elif setting.kind == "float":
            value = float(text)
        elif setting.kind == "bool":
            if isinstance(text, bool):
                value = text
            elif str(text).lower() in ("1", "true", "yes", "on"):
                value = True
            elif str(text).lower() in ("0", "false", "no", "off"):
                value = False
            else:
                raise ValueError
        elif setting.kind == "choice":
            value = str(text)
            if value not in setting.choices:
                raise ProjectError(f"{setting.name} must be one of {', '.join(setting.choices)} (got {value!r})")
        else:
            value = str(text)
            if not value:
                raise ProjectError(f"{setting.name} must not be empty")
    except ProjectError:
        raise
    except ValueError:
        raise ProjectError(f"{setting.name} expects {setting.kind}, got {raw!r}") from None
    if setting.kind in ("int", "float"):
        if setting.minimum is not None and value < setting.minimum:
            raise ProjectError(f"{setting.name} must be >= {setting.minimum:g} (got {value})")
        if setting.maximum is not None and value > setting.maximum:
            raise ProjectError(f"{setting.name} must be <= {setting.maximum:g} (got {value})")
    return value


def set_value(project: Project, name: str, raw: Any) -> Any:
    """Store an override (validated). Setting the default value removes the override."""
    setting = lookup(name)
    value = coerce(setting, raw)
    stage = project.settings.setdefault(setting.stage, {})
    if value == setting.default:
        stage.pop(setting.key, None)
    else:
        stage[setting.key] = value
    if not stage:
        project.settings.pop(setting.stage, None)
    return value


def reset(project: Project, name: str) -> None:
    setting = lookup(name)
    project.settings.get(setting.stage, {}).pop(setting.key, None)
    if not project.settings.get(setting.stage):
        project.settings.pop(setting.stage, None)


def current(project: Project, setting: Setting) -> Any:
    stage = project.settings.get(setting.stage, {})
    return stage[setting.key] if setting.key in stage else setting.default


def parse_assignments(items: Sequence[str]) -> List[Tuple[str, str]]:
    pairs = []
    for item in items:
        name, sep, value = item.partition("=")
        if not sep or "." not in name:
            raise ProjectError(f"--set needs STAGE.KEY=VALUE, got {item!r}")
        pairs.append((name.strip(), value))
    return pairs


def apply_to_block(block: Dict[str, Any], overrides: Dict[str, Any]) -> Dict[str, Any]:
    """Write overrides (dotted keys) into a generated stage block."""
    for key, value in overrides.items():
        parts = key.split(".")
        target = block
        for part in parts[:-1]:
            target = target.setdefault(part, {})
        target[parts[-1]] = value  # None is written as null: "no limit" for optional settings
    return block


def describe(project: Project, stage: Optional[str] = None) -> List[Dict[str, Any]]:
    rows = []
    for setting in settings_for(stage):
        value = current(project, setting)
        rows.append({
            "name": setting.name, "value": value, "default": setting.default,
            "changed": setting.key in project.settings.get(setting.stage, {}),
            "type": setting.kind if setting.kind != "choice" else "|".join(setting.choices),
            "help": setting.help,
        })
    return rows
