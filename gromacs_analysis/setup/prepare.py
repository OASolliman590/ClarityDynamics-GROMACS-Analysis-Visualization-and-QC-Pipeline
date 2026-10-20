"""Prepare analysis inputs the scan reports as missing, one reviewed script per replicate.

Tasks per replicate (only when the input is missing, unless forced):

- ``trjconv``: PBC-corrected, centred and fitted trajectory (``<raw>_noPBC.xtc``) from
  the run input (``.tpr``) and the raw trajectory. Groups are passed by name
  (``Protein``, ``System``, ``Backbone``), never by number.
- ``metrics``: the MD metric files (RMSD/RMSF/Rg/SASA/H-bonds/COM) with the existing
  VMD helper, written next to the trajectory where the scan looks for them.
- ``mmpbsa``: gmx_MMPBSA for holo systems. The index is built with ``gmx select`` so
  it holds exactly two groups (0 = receptor, 1 = ligand) and ``-cg 0 1`` is always right.
  The ligand is the first included ligand of the system.

Scripts are written under ``<root>/clarity/prep/`` and run locally, in the background
on the data host, or as one HTCondor job per script.
"""

from __future__ import annotations

import posixpath
import re
import shlex
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from ..pipeline.upstream import _build_md_vmd_script
from .backends import Backend, BackendError
from .jobs import new_record
from .project import Project, Replicate, SystemEntry
from .runner import CLUSTER_RE, RUN_DIR, _classad, _memory_mb, pool_policy

TASKS = ("trjconv", "metrics", "mmpbsa")
PREP_DIR = "prep"
TEXT_TOPOLOGY = (".psf", ".pdb", ".gro")

MMPBSA_IN = """&general
  sys_name="{name}", startframe={start}, endframe={end}, interval={interval}, verbose=2,
  interaction_entropy=1, ie_segment=25, temperature=298.15,
/
&gb
  igb=5, saltcon=0.150,
/
"""


@dataclass
class PrepTask:
    system: str          # display name
    system_key: str
    replicate: str       # replicate dir (relative to root)
    kind: str            # trjconv | metrics | mmpbsa
    outputs: List[str] = field(default_factory=list)
    skipped_reason: Optional[str] = None
    commands: List[str] = field(default_factory=list)
    files: Dict[str, str] = field(default_factory=dict)  # extra files to write: relative name -> text

    @property
    def ok(self) -> bool:
        return self.skipped_reason is None


@dataclass
class PrepScript:
    system: str
    replicate: str
    path: str
    text: str
    tasks: List[PrepTask]
    extra_files: Dict[str, str] = field(default_factory=dict)  # absolute path -> text


def _name(path: Optional[str]) -> Optional[str]:
    return posixpath.basename(path) if path else None


def _processed_name(raw: str) -> str:
    stem, ext = posixpath.splitext(posixpath.basename(raw))
    return f"{stem}_noPBC{ext or '.xtc'}"


def _trjconv(entry: SystemEntry, rep: Replicate, gmx: str) -> PrepTask:
    task = PrepTask(entry.name, entry.key, rep.dir, "trjconv")
    tpr, raw = _name(rep.files.get("tpr")), _name(rep.files.get("raw_trajectory"))
    if not tpr:
        task.skipped_reason = "no .tpr run input in the replicate folder"
        return task
    if not raw:
        task.skipped_reason = "no raw trajectory to process"
        return task
    out = _processed_name(raw)
    tmp = f"clarity_tmp_whole{posixpath.splitext(out)[1]}"
    task.outputs = [out]
    q = shlex.quote
    task.commands = [
        "# PBC: make molecules whole and centre the protein, then remove rotation/translation.",
        f"printf 'Protein\\nSystem\\n' | {gmx} trjconv -s {q(tpr)} -f {q(raw)} -o {q(tmp)} -pbc mol -center",
        f"printf 'Backbone\\nSystem\\n' | {gmx} trjconv -s {q(tpr)} -f {q(tmp)} -o {q(out)} -fit rot+trans",
        f"rm -f {q(tmp)}",
    ]
    return task


def _metrics(entry: SystemEntry, rep: Replicate, processed: Optional[str], vmd: str) -> PrepTask:
    task = PrepTask(entry.name, entry.key, rep.dir, "metrics")
    structure = _name(rep.files.get("structure"))
    if not structure and rep.topology and rep.topology.lower().endswith(TEXT_TOPOLOGY):
        structure = _name(rep.topology)
    if not structure:
        task.skipped_reason = "no psf/pdb/gro structure that VMD can read"
        return task
    trajectory = processed or _name(rep.trajectory)
    if not trajectory:
        task.skipped_reason = "no trajectory"
        return task
    ligands = [lig.resname for lig in entry.active_ligands]
    vmd_task = {
        "topology": structure,
        "trajectory": trajectory,
        "ligand_selection": f"resname {' '.join(ligands)}" if ligands else "resname CLARITY_NO_LIGAND",
        "metric_output_dir": ".",
        "ligands": ligands if len(ligands) > 1 else [],
    }
    task.files["clarity_md_metrics.tcl"] = _build_md_vmd_script(vmd_task)
    task.outputs = ["RMSD_*.dat", "RMSF_.dat", "rog_.dat", "SASA_.dat", "hbonds_*.dat", "comcom_*.dat"]
    task.commands = [f"{vmd} -dispdev text -e clarity_md_metrics.tcl"]
    return task


def _mmpbsa(entry: SystemEntry, rep: Replicate, processed: Optional[str], gmx: str, mpi_ranks: int,
            frames: Dict[str, int]) -> PrepTask:
    task = PrepTask(entry.name, entry.key, rep.dir, "mmpbsa")
    if entry.apo:
        task.skipped_reason = "APO reference (no ligand)"
        return task
    if not entry.active_ligands:
        task.skipped_reason = "no included ligand"
        return task
    tpr, top = _name(rep.files.get("tpr")), _name(rep.files.get("top"))
    if not tpr:
        task.skipped_reason = "no .tpr run input in the replicate folder"
        return task
    if not top:
        task.skipped_reason = "no .top topology for gmx_MMPBSA"
        return task
    trajectory = processed or _name(rep.files.get("processed_trajectory")) or _name(rep.trajectory)
    if not trajectory:
        task.skipped_reason = "no trajectory"
        return task
    ligand = entry.active_ligands[0]
    q = shlex.quote
    selection = f'"Receptor" group "Protein"; "Ligand" resname {ligand.resname}'
    runner = f"mpirun -np {int(mpi_ranks)} gmx_MMPBSA MPI" if mpi_ranks > 1 else "gmx_MMPBSA"
    task.files["clarity_mmpbsa.in"] = MMPBSA_IN.format(name=re.sub(r"\W", "_", entry.name), **frames)
    task.outputs = ["FINAL_RESULTS_MMPBSA.dat", "FINAL_RESULTS_MMPBSA.csv"]
    task.commands = [
        f"# Index with exactly two groups: 0 = receptor (Protein), 1 = ligand ({ligand.resname}).",
        f"{gmx} select -s {q(tpr)} -on clarity_mmpbsa.ndx -select {q(selection)}",
        f"{runner} -O -i clarity_mmpbsa.in -cs {q(tpr)} -ci clarity_mmpbsa.ndx -cg 0 1 -cp {q(top)} "
        f"-ct {q(trajectory)} -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv -nogui",
    ]
    return task


def plan_preparation(project: Project, what: Sequence[str] = TASKS, force: bool = False, gmx: str = "gmx",
                     vmd: str = "vmd", mpi_ranks: int = 1, frames: Optional[Dict[str, int]] = None) -> List[PrepTask]:
    """Tasks for every included system and replicate (skipped ones carry the reason)."""
    unknown = set(what) - set(TASKS)
    if unknown:
        raise ValueError(f"Unknown preparation task(s): {', '.join(sorted(unknown))} (choose from {', '.join(TASKS)})")
    frames = {"start": 1, "end": 9999999, "interval": 1, **(frames or {})}
    tasks: List[PrepTask] = []
    for entry in project.active_systems:
        for rep in entry.replicates:
            has_processed = bool(rep.files.get("processed_trajectory")) or rep.inputs.get("processed_trajectory")
            processed: Optional[str] = None
            if "trjconv" in what and (force or not has_processed):
                task = _trjconv(entry, rep, gmx)
                tasks.append(task)
                if task.ok:
                    processed = task.outputs[0]
            if "metrics" in what and (force or not rep.inputs.get("md_metrics")):
                tasks.append(_metrics(entry, rep, processed, vmd))
            if "mmpbsa" in what and not entry.apo and (force or not rep.inputs.get("mmpbsa")):
                tasks.append(_mmpbsa(entry, rep, processed, gmx, mpi_ranks, frames))
    return tasks


def _safe(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", text)


def prep_dir(project: Project) -> str:
    return posixpath.join(project.root, RUN_DIR, PREP_DIR)


def build_scripts(project: Project, tasks: Sequence[PrepTask]) -> List[PrepScript]:
    """One bash script per replicate holding its runnable tasks in order."""
    by_rep: Dict[str, List[PrepTask]] = {}
    for task in tasks:
        if task.ok:
            by_rep.setdefault(task.replicate, []).append(task)
    scripts = []
    for rep_dir, rep_tasks in by_rep.items():
        rep_abs = posixpath.join(project.root, rep_dir)
        lines = [
            "#!/bin/bash",
            f"# Generated by clarity-setup: prepare inputs for {rep_tasks[0].system} / {rep_dir}.",
            f"# Tasks: {', '.join(t.kind for t in rep_tasks)}. Review before running on new data.",
            "set -euo pipefail",
            f"cd {shlex.quote(rep_abs)}",
        ]
        extra: Dict[str, str] = {}
        for task in rep_tasks:
            lines.append(f'echo "== {task.kind}: {rep_dir}"')
            lines.extend(task.commands)
            for name, text in task.files.items():
                extra[posixpath.join(rep_abs, name)] = text
        lines.append(f'echo "== done: {rep_dir}"')
        path = posixpath.join(prep_dir(project), f"{_safe(rep_dir)}.sh")
        scripts.append(PrepScript(rep_tasks[0].system, rep_dir, path, "\n".join(lines) + "\n", rep_tasks, extra))
    return scripts


def write_scripts(backend: Backend, scripts: Sequence[PrepScript]) -> None:
    for script in scripts:
        for path, text in script.extra_files.items():
            backend.write_text(path, text)
        backend.write_text(script.path, script.text)
    if scripts:
        quoted = " ".join(shlex.quote(s.path) for s in scripts)
        result = backend.run(f"chmod +x {quoted}")
        if not result.ok:
            raise BackendError(f"Could not make preparation scripts executable: {result.stderr.strip()}")


def run_preparation(backend: Backend, project: Project, scripts: Sequence[PrepScript], where: str = "host",
                    cpus: int = 1, memory="4GB") -> Dict:
    """Write and start the scripts; ``where`` is ``local``/``host`` (background, sequential) or ``condor``."""
    if not scripts:
        raise BackendError("Nothing to prepare: every included replicate already has these inputs "
                           "(use --force to redo them)")
    write_scripts(backend, scripts)
    base = prep_dir(project)
    log = posixpath.join(base, "prepare.log")
    if where in ("local", "host"):
        runner = posixpath.join(base, "run_all.sh")
        exit_file = posixpath.join(base, "prepare.exit")
        body = ["#!/bin/bash", "# Generated by clarity-setup: run every preparation script in turn.", "status=0"]
        body += [f"{shlex.quote(s.path)} || status=1" for s in scripts]
        body += [f"echo $status > {shlex.quote(exit_file)}", "exit $status", ""]
        backend.write_text(runner, "\n".join(body))
        result = backend.run(
            f"chmod +x {shlex.quote(runner)} && rm -f {shlex.quote(exit_file)} && "
            f"nohup {shlex.quote(runner)} > {shlex.quote(log)} 2>&1 < /dev/null & echo $!"
        )
        if not result.ok:
            raise BackendError(f"Could not start preparation: {result.stderr.strip()}")
        pid = result.stdout.strip().splitlines()[-1] if result.stdout.strip() else None
        record = new_record("prepare", "host", pid=pid, log=log, exit_file=exit_file, run_dir=base,
                            last_state="running", note=f"{len(scripts)} script(s)")
        return {"where": "host", "pid": pid, "log": log, "scripts": [s.path for s in scripts], "record": record}
    if where != "condor":
        raise ValueError(f"Unknown place to run preparation: {where}")
    policy = pool_policy(backend)
    listing = posixpath.join(base, "scripts.txt")
    backend.write_text(listing, "\n".join(s.path for s in scripts) + "\n")
    lines = [
        "# Generated by clarity-setup: one job per preparation script (shared filesystem).",
        "universe = vanilla",
        "executable = $(script)",
        f"initialdir = {base}",
        "getenv = True",
        "should_transfer_files = NO",
        f"request_cpus = {int(cpus)}",
        f"request_memory = {_memory_mb(memory)}MB",
        "log = logs/prepare.log",
        "output = logs/$(Cluster).$(Process).out",
        "error = logs/$(Cluster).$(Process).err",
        "notification = never",
    ]
    if policy["arch"]:
        lines.append(f'requirements = (TARGET.Arch == "{policy["arch"]}")')
    for key, value in policy["extra"].items():
        lines.append(f"+{key} = {_classad(value)}")
    lines.append(f"queue script from {listing}")
    submit = posixpath.join(base, "prepare.sub")
    backend.write_text(submit, "\n".join(lines) + "\n")
    result = backend.run(f"mkdir -p {shlex.quote(posixpath.join(base, 'logs'))} && cd {shlex.quote(base)} "
                         f"&& condor_submit {shlex.quote(posixpath.basename(submit))}")
    if not result.ok:
        raise BackendError(f"condor_submit failed: {(result.stderr or result.stdout).strip()}")
    match = CLUSTER_RE.search(result.stdout)
    if not match:
        raise BackendError(f"Could not read the cluster ID from condor_submit output:\n{result.stdout}")
    cluster = int(match.group(1))
    record = new_record("prepare", "condor", clusters=[cluster], log=posixpath.join(base, "logs", "prepare.log"),
                        run_dir=base, last_state="idle", note=f"{len(scripts)} job(s)")
    return {"where": "condor", "cluster_id": cluster, "scripts": [s.path for s in scripts], "record": record}


def summarise(tasks: Sequence[PrepTask]) -> List[str]:
    """One line per task: replicate, task, and either its main output or why it is skipped."""
    width = max([len(t.replicate) for t in tasks] + [9])
    lines = []
    for task in tasks:
        if task.ok:
            state = f"ready    -> {task.outputs[0]}" + (" ..." if len(task.outputs) > 1 else "")
        else:
            state = f"skipped  ({task.skipped_reason})"
        lines.append(f"{task.replicate:<{width}}  {task.kind:<8} {state}")
    return lines
