"""Check what the data host offers before trusting it with a run.

One shell script collects everything in a single round trip (connection, write access,
Python, the pipeline, GROMACS/VMD/gmx_MMPBSA, HTCondor tools, free space and the
filesystem type of the project folder). Each failed check carries a remedy; checks
are "required" only for the modes that need them.
"""

from __future__ import annotations

import shlex
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

from .backends import Backend

SHARED_FS = {"nfs", "nfs4", "lustre", "gpfs", "ceph", "cifs", "smb2", "smbfs", "beegfs", "panfs", "glusterfs",
             "fuse.glusterfs", "fuse.ceph", "fuse.sshfs", "autofs", "wekafs", "fuse.lustre"}
MIN_FREE_GB = 5.0

REQUIRED_FOR = {
    "connection": {"*"},
    "root": {"*"},
    "pipeline": {"local", "host", "condor", "dag"},
    "python": {"host", "condor", "dag"},
    "gmx": {"trjconv", "mmpbsa"},
    "vmd": {"metrics"},
    "gmx_MMPBSA": {"mmpbsa"},
    "condor_submit": {"condor", "dag"},
    "condor_q": {"condor", "dag"},
    "condor_submit_dag": {"dag"},
}

REMEDY = {
    "connection": "check the host name/alias, your SSH key and that the host is reachable",
    "root": "the folder must exist and be writable (results and scripts are written under it)",
    "pipeline": 'install the toolkit there: pip install --user "git+https://github.com/OASolliman590/'
                'ClarityDynamics-GROMACS-Analysis-Visualization-and-QC-Pipeline"',
    "python": "install Python 3.10 or newer, or load it (e.g. module load python)",
    "gmx": "load GROMACS (module load gromacs / source GMXRC) so gmx or gmx_mpi is on PATH",
    "vmd": "install VMD or load its module so `vmd` is on PATH",
    "gmx_MMPBSA": "pip install gmx_MMPBSA (needs AmberTools) or load its module",
    "condor_submit": "this host cannot submit HTCondor jobs: use the pool's submit node (e.g. a condorctl profile)",
    "condor_q": "HTCondor tools missing on this host",
    "condor_submit_dag": "DAGMan missing: run --where condor instead of dag",
    "disk": f"less than {MIN_FREE_GB:g} GB free under the project folder: trajectories and outputs need space",
    "shared_fs": "the project folder looks local to this host; HTCondor execute nodes may not see it. "
                 "Put the data on the pool's shared filesystem (should_transfer_files = NO is used)",
}


@dataclass
class Check:
    name: str
    ok: bool
    value: str
    required: bool
    remedy: str = ""

    @property
    def status(self) -> str:
        if self.ok:
            return "ok"
        return "MISSING" if self.required else "warning"


def _script(root: str) -> str:
    q = shlex.quote(root)
    probe = shlex.quote(root.rstrip("/") + "/.clarity_write_test")
    tools = ["gromacs-pipeline", "python3", "vmd", "gmx_MMPBSA", "condor_submit", "condor_q", "condor_submit_dag"]
    lines = ["echo connection=ok"]
    lines.append(f"if test -d {q} && touch {probe} 2>/dev/null; then rm -f {probe}; echo root=writable; "
                 f"elif test -d {q}; then echo root=read-only; else echo root=missing; fi")
    for tool in tools:
        lines.append(f'echo "{tool}=$(command -v {tool} 2>/dev/null)"')
    lines.append('echo "gmx=$(command -v gmx 2>/dev/null || command -v gmx_mpi 2>/dev/null)"')
    lines.append('echo "python_version=$(python3 -c \'import sys;print(sys.version.split()[0])\' 2>/dev/null)"')
    lines.append(f'echo "disk_kb=$(df -Pk {q} 2>/dev/null | tail -1 | awk \'{{print $4}}\')"')
    lines.append(f'echo "fstype=$(stat -f -c %T {q} 2>/dev/null || df -PT {q} 2>/dev/null | tail -1 | awk \'{{print $2}}\')"')
    return "; ".join(lines)


def _parse(text: str) -> Dict[str, str]:
    values = {}
    for line in text.splitlines():
        key, sep, value = line.partition("=")
        if sep:
            values[key.strip()] = value.strip()
    return values


def run_checks(backend: Backend, root: str, modes: Sequence[str] = ("host",)) -> List[Check]:
    """All checks for ``root``; ``modes`` (local/host/condor/dag/trjconv/metrics/mmpbsa) decide what is required."""
    wanted = set(modes)

    def required(name: str) -> bool:
        needs = REQUIRED_FOR.get(name, set())
        return "*" in needs or bool(needs & wanted)

    try:
        result = backend.run(_script(root), timeout=60)
        values = _parse(result.stdout)
    except Exception as exc:  # noqa: BLE001 - report instead of raising
        return [Check("connection", False, str(exc), True, REMEDY["connection"])]
    checks = [Check("connection", values.get("connection") == "ok", backend.label, True, REMEDY["connection"])]
    root_state = values.get("root", "unknown")
    checks.append(Check("root", root_state == "writable", f"{root} ({root_state})", True, REMEDY["root"]))
    version = values.get("python_version", "")
    py_ok = bool(version) and tuple(int(p) for p in version.split(".")[:2] if p.isdigit()) >= (3, 10)
    checks.append(Check("python", py_ok, version or "not found", required("python"), REMEDY["python"]))
    for name, key in (("pipeline", "gromacs-pipeline"), ("gmx", "gmx"), ("vmd", "vmd"), ("gmx_MMPBSA", "gmx_MMPBSA"),
                      ("condor_submit", "condor_submit"), ("condor_q", "condor_q"),
                      ("condor_submit_dag", "condor_submit_dag")):
        path = values.get(key, "")
        checks.append(Check(name, bool(path), path or "not found", required(name), REMEDY[name]))
    try:
        free_gb = int(values.get("disk_kb", "")) / 1024 / 1024
        checks.append(Check("disk", free_gb >= MIN_FREE_GB, f"{free_gb:.1f} GB free", False, REMEDY["disk"]))
    except ValueError:
        checks.append(Check("disk", False, "unknown", False, REMEDY["disk"]))
    fstype = values.get("fstype", "").lower()
    if wanted & {"condor", "dag"} or values.get("condor_submit"):
        shared = fstype in SHARED_FS or fstype.startswith(("nfs", "fuse."))
        checks.append(Check("shared_fs", shared, fstype or "unknown", False, REMEDY["shared_fs"]))
    return checks


def exit_code(checks: Sequence[Check]) -> int:
    return 4 if any(c.required and not c.ok for c in checks) else 0


def format_table(checks: Sequence[Check]) -> List[str]:
    width = max(len(c.name) for c in checks)
    lines = []
    for check in checks:
        lines.append(f"{check.name:<{width}}  {check.status:<8} {check.value}")
        if not check.ok:
            lines.append(f"{'':<{width}}  {'':<8} -> {check.remedy}")
    return lines


def modes_for(where: Optional[str], prepare: Sequence[str] = ()) -> List[str]:
    return [m for m in [where or "host", *prepare] if m]
