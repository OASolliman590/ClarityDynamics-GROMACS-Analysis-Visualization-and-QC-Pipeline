"""Discover simulation systems under a project root (locally or over SSH).

One ``find`` lists every relevant file; grouping into systems and replicates is
done locally. Ligands are then read with one short command per unique topology.
"""

from __future__ import annotations

import posixpath
import re
import shlex
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .backends import Backend, BackendError
from .ligands import filter_ligands, parse_resnames, pick_text_topology, resname_command

TRAJECTORY_EXT = (".xtc", ".trr", ".dcd", ".nc")
TOPOLOGY_EXT = (".psf", ".tpr", ".pdb", ".gro", ".top")
METRIC_PREFIXES = ("rmsd", "rmsf", "rog", "sasa", "hbond", "comcom", "com")
APO_TOKENS = ("apo", "free", "unbound", "noligand", "no_ligand", "ligandfree")
HOLO_PREFIXES = ("holo", "complex", "cplx")

FIND_NAMES = [
    "*.xtc", "*.trr", "*.dcd", "*.nc", "*.psf", "*.tpr", "*.pdb", "*.gro", "*.top",
    "RMSD*", "RMSF*", "rog*", "RoG*", "SASA*", "sasa*", "hbond*", "Hbond*", "comcom*", "COM*",
    "FINAL_RESULTS_MMPBSA*", "eigenval*.xvg",
]

_REP_SUFFIX = re.compile(r"^(?P<base>.+?)(?:[_\-.](?:rep|replica|run|r|R)?|(?:rep|replica|run|R))(?P<num>\d+)$")
_REP_COMPONENT = re.compile(r"^(?:rep|replica|run|r|R)?[_\-]?(?P<num>\d+)$")


@dataclass
class FileInfo:
    path: str
    size: int = 0

    @property
    def name(self) -> str:
        return posixpath.basename(self.path)


@dataclass
class Replicate:
    dir: str                      # relative to root
    number: int
    topology: Optional[str]       # relative path used by the pipeline
    text_topology: Optional[str]  # relative path read for residue names
    trajectories: List[str]       # relative paths, preferred first
    metric_dir: Optional[str] = None
    inputs: Dict[str, bool] = field(default_factory=dict)
    # Files needed to prepare missing inputs: tpr, top, raw_trajectory, processed_trajectory.
    files: Dict[str, Optional[str]] = field(default_factory=dict)


@dataclass
class DiscoveredSystem:
    key: str
    default_name: str
    apo: bool
    replicates: List[Replicate]
    ligands: List[str] = field(default_factory=list)
    all_resnames: List[str] = field(default_factory=list)

    @property
    def inputs(self) -> Dict[str, bool]:
        merged: Dict[str, bool] = {}
        for rep in self.replicates:
            for key, value in rep.inputs.items():
                merged[key] = merged.get(key, False) or value
        return merged


@dataclass
class ScanResult:
    root: str
    systems: List[DiscoveredSystem]
    skipped_dirs: List[str]
    commands_used: int
    backend: str

    def system(self, key: str) -> Optional[DiscoveredSystem]:
        return next((s for s in self.systems if s.key == key), None)


# ---------------------------------------------------------------- listing

def _find_command(root: str, depth: int, printf: bool = True) -> str:
    names = " -o ".join(f"-name {shlex.quote(n)}" for n in FIND_NAMES)
    tail = "-printf '%p\\t%s\\n'" if printf else "-print"
    return f"find {shlex.quote(root)} -maxdepth {int(depth)} -type f \\( {names} \\) {tail}"


def list_files(backend: Backend, root: str, depth: int = 4) -> List[FileInfo]:
    result = backend.run(_find_command(root, depth))
    if not result.ok or ("\t" not in result.stdout and result.stdout.strip()):
        # BSD/BusyBox find without -printf: names only.
        result = backend.run(_find_command(root, depth, printf=False))
        if not result.ok:
            raise BackendError(f"Cannot list {root} on {backend.label}: {result.stderr.strip()}")
        return [FileInfo(line.strip()) for line in result.stdout.splitlines() if line.strip()]
    files = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        path, _, size = line.rpartition("\t")
        try:
            files.append(FileInfo(path, int(size)))
        except ValueError:
            files.append(FileInfo(line.strip()))
    return files


# ---------------------------------------------------------------- grouping helpers

def split_replicate(rel_dir: str) -> Tuple[str, int]:
    """('HOLO_LigA_2') -> ('HOLO_LigA', 2); ('LigA/rep1') -> ('LigA', 1); no token -> (dir, 1)."""
    parent, last = posixpath.split(rel_dir)
    whole = _REP_COMPONENT.match(last)
    if whole and parent:
        return parent, int(whole.group("num"))
    match = _REP_SUFFIX.match(last)
    if match:
        base = match.group("base").rstrip("_-.")
        return (posixpath.join(parent, base) if parent else base), int(match.group("num"))
    return rel_dir, 1


def default_display_name(key: str, apo: bool) -> str:
    last = posixpath.basename(key) or key
    if apo:
        return "APO"
    lowered = last.lower()
    for prefix in HOLO_PREFIXES:
        for sep in ("_", "-", "."):
            if lowered.startswith(prefix + sep) and len(last) > len(prefix) + 1:
                return last[len(prefix) + 1:]
    return last


def looks_apo(key: str) -> bool:
    tokens = re.split(r"[_\-./]", key.lower())
    return any(token in APO_TOKENS for token in tokens) or any(t in key.lower() for t in ("noligand", "ligandfree"))


def _trajectory_rank(info: FileInfo) -> Tuple[Tuple[int, int, int], int]:
    name = info.name.lower()
    calpha = "calpha" in name or "_ca." in name or name.startswith("ca_")
    processed = any(t in name for t in ("nopbc", "nojump", "center", "fit", "whole"))
    production = "prod" in name or "md" in name
    return (1 if calpha else 0, 0 if processed else 1, 0 if production else 1), -info.size


def _choose_trajectories(files: List[FileInfo]) -> List[str]:
    trajs = [f for f in files if f.name.lower().endswith(TRAJECTORY_EXT)]
    trajs.sort(key=lambda f: (_trajectory_rank(f)[0], _trajectory_rank(f)[1], f.name))
    return [f.path for f in trajs]


PROCESSED_TOKENS = ("nopbc", "nojump", "center", "fit", "whole")


def _is_processed(name: str) -> bool:
    return any(t in name.lower() for t in PROCESSED_TOKENS)


def _prep_files(items: List[FileInfo]) -> Dict[str, Optional[str]]:
    """The .tpr/.top and raw/processed trajectories that input preparation needs."""
    def pick(candidates: List[FileInfo]) -> Optional[str]:
        if not candidates:
            return None
        ranked = sorted(candidates, key=lambda f: (0 if ("prod" in f.name.lower() or "md" in f.name.lower()) else 1,
                                                   -f.size, f.name))
        return ranked[0].path

    trajs = [f for f in items if f.name.lower().endswith(TRAJECTORY_EXT) and _trajectory_rank(f)[0][0] == 0]
    return {
        "tpr": pick([f for f in items if f.name.lower().endswith(".tpr")]),
        "top": pick([f for f in items if f.name.lower().endswith(".top")]),
        "raw_trajectory": pick([f for f in trajs if not _is_processed(f.name)]),
        "processed_trajectory": pick([f for f in trajs if _is_processed(f.name)]),
    }


def _choose_topology(paths: List[str]) -> Optional[str]:
    for ext in (".psf", ".tpr", ".pdb", ".gro"):
        candidates = sorted(p for p in paths if p.lower().endswith(ext))
        if candidates:
            # CHARMM-GUI: step3_input.* is the full system; prefer it.
            preferred = [p for p in candidates if "step3_input" in p or "input" in posixpath.basename(p)]
            return (preferred or candidates)[0]
    return None


def _is_metric(name: str) -> bool:
    lower = name.lower()
    return lower.endswith((".dat", ".xvg")) and lower.startswith(METRIC_PREFIXES) and not lower.startswith("eigenval")


# ---------------------------------------------------------------- scan

def scan(backend: Backend, root: str, depth: int = 4) -> ScanResult:
    """Discover systems, replicates, topologies, trajectories, ligands and analysis inputs."""
    root = root.rstrip("/") or "/"
    start = backend.commands_run
    files = list_files(backend, root, depth)

    by_dir: "OrderedDict[str, List[FileInfo]]" = OrderedDict()
    for info in sorted(files, key=lambda f: f.path):
        by_dir.setdefault(posixpath.dirname(info.path), []).append(info)

    sim_dirs = [
        d for d, items in by_dir.items()
        if any(i.name.lower().endswith(TRAJECTORY_EXT) for i in items)
        and any(i.name.lower().endswith(TOPOLOGY_EXT) for i in items)
    ]

    def owner(directory: str) -> Optional[str]:
        best = None
        for sim in sim_dirs:
            if directory == sim or directory.startswith(sim + "/"):
                if best is None or len(sim) > len(best):
                    best = sim
        return best

    extra: Dict[str, List[FileInfo]] = {sim: [] for sim in sim_dirs}
    skipped = []
    for directory, items in by_dir.items():
        if directory in extra:
            continue
        parent = owner(directory)
        if parent is None:
            skipped.append(posixpath.relpath(directory, root))
        else:
            extra[parent].extend(items)

    grouped: "OrderedDict[str, List[Replicate]]" = OrderedDict()
    for sim in sim_dirs:
        rel = posixpath.relpath(sim, root)
        key, number = split_replicate(rel)
        items = by_dir[sim]
        nested = extra[sim]
        topo_paths = [i.path for i in items if i.name.lower().endswith(TOPOLOGY_EXT)]
        metric_items = [i for i in items if _is_metric(i.name)]
        metric_dir = sim if metric_items else None
        if metric_dir is None:
            nested_metrics = sorted({posixpath.dirname(i.path) for i in nested if _is_metric(i.name)})
            metric_dir = nested_metrics[0] if nested_metrics else None
        all_items = items + nested
        replicate = Replicate(
            dir=rel,
            number=number,
            topology=_rel(_choose_topology(topo_paths), root),
            text_topology=_rel(pick_text_topology(topo_paths), root),
            trajectories=[posixpath.relpath(p, root) for p in _choose_trajectories(items)],
            metric_dir=_rel(metric_dir, root),
            inputs={
                "md_metrics": metric_dir is not None,
                "mmpbsa": any(i.name.startswith("FINAL_RESULTS_MMPBSA") for i in all_items),
                "pca_dir": any(i.name.lower().startswith("eigenval") for i in all_items),
                "processed_trajectory": any(
                    t in i.name.lower() for i in items for t in ("nopbc", "nojump")
                ),
            },
            files={**{k: _rel(v, root) for k, v in _prep_files(items).items()},
                   "structure": _rel(pick_text_topology(topo_paths), root)},
        )
        grouped.setdefault(key, []).append(replicate)

    systems: List[DiscoveredSystem] = []
    for key, reps in grouped.items():
        reps.sort(key=lambda r: r.number)
        apo = looks_apo(key)
        systems.append(DiscoveredSystem(key=key, default_name=default_display_name(key, apo), apo=apo,
                                        replicates=reps))

    _detect_ligands(backend, root, systems)

    return ScanResult(
        root=root,
        systems=systems,
        skipped_dirs=sorted(set(skipped)),
        commands_used=backend.commands_run - start,
        backend=backend.label,
    )


def _rel(path: Optional[str], root: str) -> Optional[str]:
    if path is None:
        return None
    return posixpath.relpath(path, root)


def _detect_ligands(backend: Backend, root: str, systems: List[DiscoveredSystem]) -> None:
    cache: Dict[str, List[str]] = {}
    for system in systems:
        topology = next((r.text_topology for r in system.replicates if r.text_topology), None)
        if topology is None:
            continue
        if topology not in cache:
            command = resname_command(posixpath.join(root, topology))
            result = backend.run(command) if command else None
            cache[topology] = parse_resnames(result.stdout) if result is not None and result.ok else []
        system.all_resnames = cache[topology]

    apo_names = {n for s in systems if s.apo for n in s.all_resnames}
    for system in systems:
        if system.apo:
            system.ligands = []
            continue
        system.ligands = filter_ligands(system.all_resnames, shared=apo_names)
        if system.all_resnames and not system.ligands and not apo_names:
            # A system with a readable topology but no ligand residue is a reference system.
            system.apo = True
            system.default_name = default_display_name(system.key, True) if looks_apo(system.key) else system.default_name
