"""Project state: what was discovered plus the user's names, exclusions and stage choices.

Saved as ``clarity_project.yaml`` so a rescan keeps every rename. Selections always
use the original residue names; display names are only labels and output names.
"""

from __future__ import annotations

import datetime as _dt
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

from .scanner import ScanResult

PROJECT_FILE = "clarity_project.yaml"
PROJECT_VERSION = 1


class ProjectError(ValueError):
    """Invalid user edits (duplicate or empty names, unknown systems...)."""


@dataclass
class Connection:
    kind: str = "local"              # local | ssh | condorctl
    host: Optional[str] = None       # ssh alias / hostname
    profile: Optional[str] = None    # condorctl config path

    @property
    def is_remote(self) -> bool:
        return self.kind != "local"

    def describe(self) -> str:
        if self.kind == "ssh":
            return f"ssh {self.host}"
        if self.kind == "condorctl":
            return f"condorctl {self.profile or '(default config)'}"
        return "local"


@dataclass
class Replicate:
    dir: str
    number: int
    topology: Optional[str]
    trajectories: List[str]
    metric_dir: Optional[str] = None
    inputs: Dict[str, bool] = field(default_factory=dict)
    files: Dict[str, Optional[str]] = field(default_factory=dict)

    @property
    def trajectory(self) -> Optional[str]:
        return self.trajectories[0] if self.trajectories else None


@dataclass
class LigandEntry:
    resname: str
    name: str
    included: bool = True

    @property
    def selection(self) -> str:
        return f"resname {self.resname}"


@dataclass
class SystemEntry:
    key: str
    name: str
    included: bool = True
    apo: bool = False
    replicates: List[Replicate] = field(default_factory=list)
    ligands: List[LigandEntry] = field(default_factory=list)
    missing: bool = False

    @property
    def inputs(self) -> Dict[str, bool]:
        merged: Dict[str, bool] = {}
        for rep in self.replicates:
            for key, value in rep.inputs.items():
                merged[key] = merged.get(key, False) or bool(value)
        return merged

    @property
    def active_ligands(self) -> List[LigandEntry]:
        return [lig for lig in self.ligands if lig.included]

    @property
    def topology(self) -> Optional[str]:
        return next((r.topology for r in self.replicates if r.topology), None)

    @property
    def trajectories(self) -> List[str]:
        return [r.trajectory for r in self.replicates if r.trajectory]


@dataclass
class StageChoice:
    enabled: bool = False
    available: bool = True
    reason: str = ""


@dataclass
class Project:
    root: str
    connection: Connection = field(default_factory=Connection)
    protein_name: str = "Protein"
    project_name: str = "ClarityDynamics project"
    output_root: str = "analysis_outputs"
    systems: List[SystemEntry] = field(default_factory=list)
    stages: Dict[str, StageChoice] = field(default_factory=dict)
    scanned_at: Optional[str] = None
    version: int = PROJECT_VERSION
    settings: Dict[str, Dict[str, Any]] = field(default_factory=dict)   # stage -> {key: override}
    runs: List[Dict[str, Any]] = field(default_factory=list)            # run records (setup.jobs)

    # ------------------------------------------------------------ construction

    @classmethod
    def from_scan(cls, scan: ScanResult, connection: Optional[Connection] = None, **kwargs) -> "Project":
        project = cls(root=scan.root, connection=connection or Connection(), **kwargs)
        project.merge_scan(scan)
        return project

    def merge_scan(self, scan: ScanResult) -> None:
        """Refresh discovered fields, keep user edits, add new systems, flag vanished ones."""
        self.root = scan.root
        self.scanned_at = _dt.datetime.now().isoformat(timespec="seconds")
        existing = {s.key: s for s in self.systems}
        seen = set()
        for found in scan.systems:
            seen.add(found.key)
            entry = existing.get(found.key)
            replicates = [_replicate(r) for r in found.replicates]
            if entry is None:
                name = self._unique_name(found.default_name)
                self.systems.append(SystemEntry(
                    key=found.key, name=name, apo=found.apo, replicates=replicates,
                    ligands=[LigandEntry(resname=res, name=res) for res in found.ligands],
                ))
                continue
            entry.replicates = replicates
            entry.missing = False
            known = {lig.resname: lig for lig in entry.ligands}
            entry.ligands = [known.get(res) or LigandEntry(resname=res, name=res) for res in found.ligands]
        for entry in self.systems:
            if entry.key not in seen:
                entry.missing = True

    def _unique_name(self, base: str) -> str:
        taken = {s.name for s in self.systems}
        name = sanitize_name(base) or "system"
        candidate, i = name, 2
        while candidate in taken:
            candidate, i = f"{name}_{i}", i + 1
        return candidate

    # ------------------------------------------------------------ edits

    def system(self, key_or_name: str) -> SystemEntry:
        for entry in self.systems:
            if entry.key == key_or_name or entry.name == key_or_name:
                return entry
        raise ProjectError(f"No system named {key_or_name!r}")

    def rename_system(self, key_or_name: str, new_name: str) -> None:
        entry = self.system(key_or_name)
        name = sanitize_name(new_name)
        if not name:
            raise ProjectError("System name must not be empty")
        clash = next((s for s in self.systems if s is not entry and s.included and s.name == name), None)
        if clash is not None:
            raise ProjectError(f"Name {name!r} is already used by system {clash.key!r}")
        entry.name = name

    def rename_ligand(self, resname: str, new_name: str, system: Optional[str] = None) -> int:
        """Rename a ligand in one system, or in every system that has the residue. Returns count."""
        name = sanitize_name(new_name)
        if not name:
            raise ProjectError("Ligand name must not be empty")
        targets = [self.system(system)] if system else self.systems
        changed = 0
        for entry in targets:
            for lig in entry.ligands:
                if lig.resname == resname:
                    clash = [o for o in entry.ligands if o is not lig and o.included and o.name == name]
                    if clash:
                        raise ProjectError(f"Ligand name {name!r} already used in {entry.name}")
                    lig.name = name
                    changed += 1
        if not changed:
            raise ProjectError(f"No ligand with residue name {resname!r}")
        return changed

    def set_apo(self, key_or_name: str, apo: bool = True) -> None:
        self.system(key_or_name).apo = apo

    def set_included(self, key_or_name: str, included: bool) -> None:
        entry = self.system(key_or_name)
        if included:
            clash = next((s for s in self.systems if s is not entry and s.included and s.name == entry.name), None)
            if clash is not None:
                raise ProjectError(f"Cannot include {entry.name!r}: name also used by {clash.key!r}")
        entry.included = included

    def set_ligand_included(self, resname: str, included: bool, system: Optional[str] = None) -> None:
        targets = [self.system(system)] if system else self.systems
        for entry in targets:
            for lig in entry.ligands:
                if lig.resname == resname:
                    lig.included = included

    @property
    def active_systems(self) -> List[SystemEntry]:
        return [s for s in self.systems if s.included and not s.missing]

    def validate(self) -> List[str]:
        """Problems that block generation (empty list = OK)."""
        problems = []
        names = [s.name for s in self.active_systems]
        for name in sorted({n for n in names if names.count(n) > 1}):
            problems.append(f"Duplicate system name: {name}")
        for entry in self.active_systems:
            if not entry.topology:
                problems.append(f"{entry.name}: no topology found")
            if not entry.trajectories:
                problems.append(f"{entry.name}: no trajectory found")
            lig_names = [lig.name for lig in entry.active_ligands]
            for name in sorted({n for n in lig_names if lig_names.count(n) > 1}):
                problems.append(f"{entry.name}: duplicate ligand name {name}")
        if not self.active_systems:
            problems.append("No systems selected")
        return problems

    # ------------------------------------------------------------ persistence

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "Project":
        systems = []
        for s in data.get("systems", []):
            systems.append(SystemEntry(
                key=s["key"], name=s["name"], included=s.get("included", True), apo=s.get("apo", False),
                missing=s.get("missing", False),
                replicates=[Replicate(**r) for r in s.get("replicates", [])],
                ligands=[LigandEntry(**lig) for lig in s.get("ligands", [])],
            ))
        return cls(
            root=data["root"],
            connection=Connection(**(data.get("connection") or {})),
            protein_name=data.get("protein_name", "Protein"),
            project_name=data.get("project_name", "ClarityDynamics project"),
            output_root=data.get("output_root", "analysis_outputs"),
            systems=systems,
            stages={k: StageChoice(**v) for k, v in (data.get("stages") or {}).items()},
            scanned_at=data.get("scanned_at"),
            version=int(data.get("version", PROJECT_VERSION)),
            settings={k: dict(v or {}) for k, v in (data.get("settings") or {}).items()},
            runs=list(data.get("runs") or []),
        )

    def save(self, path: Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        header = (
            "# ClarityDynamics project file, written by clarity-setup.\n"
            "# Holds discovery results and your names/choices; edit with `clarity-setup tui`.\n"
        )
        path.write_text(header + yaml.safe_dump(self.to_dict(), sort_keys=False, allow_unicode=True))
        return path

    @classmethod
    def load(cls, path: Path) -> "Project":
        data = yaml.safe_load(Path(path).read_text()) or {}
        if "root" not in data:
            raise ProjectError(f"{path} is not a clarity project file (no 'root')")
        return cls.from_dict(data)


def _replicate(rep) -> Replicate:
    return Replicate(
        dir=rep.dir, number=rep.number, topology=rep.topology, trajectories=list(rep.trajectories),
        metric_dir=rep.metric_dir, inputs=dict(rep.inputs), files=dict(getattr(rep, "files", {}) or {}),
    )


def sanitize_name(name: str) -> str:
    """Names become folder and file names: spaces -> '_', other unsafe characters removed."""
    cleaned = re.sub(r"\s+", "_", str(name).strip())
    cleaned = re.sub(r"[^A-Za-z0-9_.+\-]", "", cleaned)
    return cleaned
