"""Follow runs after they start: records in the project file, status, and fetching results.

- host runs: process liveness (``kill -0``), exit code written by the wrapper, log tail;
- HTCondor jobs and DAG workflows: ``condor_q``, then ``condor_history`` once jobs
  have left the queue (JobStatus 1..7);
- fetch: one ``tar`` of the output folder on the remote side, one download, unpacked
  locally (no per-file round trips).
"""

from __future__ import annotations

import posixpath
import shlex
import tarfile
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from .backends import Backend, BackendError
from .project import Project, ProjectError

JOB_STATES = {1: "idle", 2: "running", 3: "removed", 4: "completed", 5: "held", 6: "transferring", 7: "suspended"}
FINISHED = {"completed", "failed", "removed"}


@dataclass
class RunRecord:
    id: str
    kind: str                         # analysis | prepare
    where: str                        # local | host | condor | dag
    started: str
    log: Optional[str] = None
    output_dir: Optional[str] = None  # analysis output root (what fetch downloads)
    run_dir: Optional[str] = None     # folder with the generated scripts/configs
    report: Optional[str] = None
    pid: Optional[str] = None
    exit_file: Optional[str] = None
    clusters: List[int] = field(default_factory=list)
    dag_cluster: Optional[int] = None
    jobs: Dict[str, Any] = field(default_factory=dict)
    last_state: Optional[str] = None
    fetched_to: Optional[str] = None
    note: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RunRecord":
        known = {f for f in cls.__dataclass_fields__}  # noqa: WPS609
        return cls(**{k: v for k, v in data.items() if k in known})


def new_record(kind: str, where: str, **fields: Any) -> RunRecord:
    stamp = datetime.now()
    return RunRecord(id=f"{kind}-{stamp:%Y%m%d-%H%M%S}", kind=kind, where=where,
                     started=stamp.isoformat(timespec="seconds"), **fields)


def add_run(project: Project, record: RunRecord) -> RunRecord:
    ids = {r.get("id") for r in project.runs}
    base, n = record.id, 2
    while record.id in ids:
        record.id, n = f"{base}-{n}", n + 1
    project.runs.append(record.to_dict())
    return record


def get_run(project: Project, run_id: Optional[str] = None) -> RunRecord:
    if not project.runs:
        raise ProjectError("No runs recorded in this project yet (start one with run or prepare)")
    if run_id is None:
        return RunRecord.from_dict(project.runs[-1])
    for data in project.runs:
        if data.get("id") == run_id:
            return RunRecord.from_dict(data)
    raise ProjectError(f"No run {run_id!r}; known: {', '.join(r['id'] for r in project.runs)}")


def update_run(project: Project, record: RunRecord) -> None:
    for i, data in enumerate(project.runs):
        if data.get("id") == record.id:
            project.runs[i] = record.to_dict()
            return
    project.runs.append(record.to_dict())


# ---------------------------------------------------------------- status

def _tail(backend: Backend, path: Optional[str], lines: int) -> str:
    if not path or lines <= 0:
        return ""
    result = backend.run(f"tail -n {int(lines)} {shlex.quote(path)} 2>/dev/null")
    return result.stdout if result.ok else ""


def _host_status(backend: Backend, record: RunRecord) -> Dict[str, Any]:
    """The wrapper writes the exit code when the run ends; that file wins over the PID check
    (a finished process can linger as a zombie and still answer ``kill -0``)."""
    exit_part = f"cat {shlex.quote(record.exit_file)} 2>/dev/null; " if record.exit_file else ""
    script = (f"echo exit=$({exit_part}true); "
              f"if kill -0 {shlex.quote(str(record.pid))} 2>/dev/null; then echo alive; else echo gone; fi")
    out = backend.run(script).stdout.split()
    code_text = out[0].partition("=")[2] if out and out[0].startswith("exit=") else ""
    if code_text.lstrip("-").isdigit():
        code = int(code_text)
        return {"state": "completed" if code == 0 else "failed", "detail": f"exit code {code}"}
    if "alive" in out:
        return {"state": "running", "detail": f"PID {record.pid} is running"}
    return {"state": "unknown", "detail": f"PID {record.pid} is not running and no exit code was recorded"}


def _parse_jobs(text: str, with_exit: bool) -> List[Dict[str, Any]]:
    jobs = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) < 3 or not parts[0].isdigit():
            continue
        status = int(parts[2]) if parts[2].isdigit() else 0
        job = {"id": f"{parts[0]}.{parts[1]}", "cluster": int(parts[0]), "state": JOB_STATES.get(status, "unknown")}
        if with_exit and len(parts) > 3 and parts[3].lstrip("-").isdigit():
            job["exit_code"] = int(parts[3])
        jobs.append(job)
    return jobs


def _condor_constraint(record: RunRecord) -> str:
    if record.dag_cluster is not None:
        return f"DAGManJobId == {int(record.dag_cluster)} || ClusterId == {int(record.dag_cluster)}"
    clusters = [int(c) for c in record.clusters] or [-1]
    return " || ".join(f"ClusterId == {c}" for c in clusters)


def _condor_status(backend: Backend, record: RunRecord) -> Dict[str, Any]:
    constraint = shlex.quote(_condor_constraint(record))
    queued = backend.run(f"condor_q -constraint {constraint} -af ClusterId ProcId JobStatus")
    jobs = _parse_jobs(queued.stdout, with_exit=False) if queued.ok else []
    in_queue = bool(jobs)
    done = backend.run(f"condor_history -constraint {constraint} -af ClusterId ProcId JobStatus ExitCode -limit 1000")
    seen = {j["id"] for j in jobs}
    if done.ok:
        jobs += [j for j in _parse_jobs(done.stdout, with_exit=True) if j["id"] not in seen]
    if not jobs:
        reason = (queued.stderr or done.stderr or "").strip().splitlines()
        return {"state": "unknown", "jobs": [], "counts": {},
                "detail": reason[-1] if reason else "no matching jobs in condor_q or condor_history"}
    counts: Dict[str, int] = {}
    for job in jobs:
        key = job["state"]
        if key == "completed" and job.get("exit_code", 0) != 0:
            key = job["state"] = "failed"
        counts[key] = counts.get(key, 0) + 1
    if in_queue:
        state = "held" if counts.get("held") and not counts.get("running") else (
            "running" if counts.get("running") else "idle")
    elif counts.get("removed"):
        state = "removed"
    elif counts.get("failed"):
        state = "failed"
    else:
        state = "completed"
    if record.dag_cluster is not None:
        dag_job = next((j for j in jobs if j["cluster"] == int(record.dag_cluster)), None)
        if dag_job and dag_job["state"] in ("running", "idle", "held"):
            state = dag_job["state"] if dag_job["state"] != "idle" else "running"
    summary = ", ".join(f"{n} {k}" for k, n in sorted(counts.items()))
    return {"state": state, "jobs": jobs, "counts": counts, "detail": summary}


def status(backend: Backend, record: RunRecord, lines: int = 20) -> Dict[str, Any]:
    """State of a recorded run plus the last ``lines`` of its log."""
    if record.where == "host":
        info = _host_status(backend, record)
    elif record.where in ("condor", "dag"):
        info = _condor_status(backend, record)
    else:
        info = {"state": record.last_state or "completed", "detail": "ran in the foreground"}
    info.setdefault("jobs", [])
    info.setdefault("counts", {})
    info["log_tail"] = _tail(backend, record.log, lines)
    record.last_state = info["state"]
    return info


# ---------------------------------------------------------------- fetch

def _safe_extract(archive: Path, dest: Path) -> None:
    dest = dest.resolve()
    with tarfile.open(archive, "r:gz") as tar:
        members = []
        for member in tar.getmembers():
            target = (dest / member.name).resolve()
            if not str(target).startswith(str(dest)) or member.issym() or member.islnk() or member.isdev():
                continue  # never write outside dest, never follow links
            members.append(member)
        if hasattr(tarfile, "data_filter"):  # Python >= 3.12 (and patched 3.10/3.11): also strip unsafe metadata
            tar.extractall(dest, members=members, filter="data")
        else:  # pragma: no cover - older Pythons
            tar.extractall(dest, members=members)  # noqa: S202 - members validated above


def fetch(backend: Backend, record: RunRecord, dest: Path, force: bool = False) -> Path:
    """Download the run's output folder as one archive and unpack it under ``dest/<run id>``."""
    if not record.output_dir:
        raise BackendError(f"Run {record.id} has no output folder to fetch")
    if record.last_state not in FINISHED and not force:
        raise BackendError(f"Run {record.id} is {record.last_state or 'not checked yet'}; "
                           "check status first, or fetch anyway with --force")
    if not backend.exists(record.output_dir):
        raise BackendError(f"{record.output_dir} does not exist on {backend.label}")
    archive = posixpath.join(record.run_dir or posixpath.dirname(record.output_dir), f"{record.id}-results.tgz")
    made = backend.run(f"tar czf {shlex.quote(archive)} -C {shlex.quote(record.output_dir)} .", timeout=3600)
    if not made.ok:
        raise BackendError(f"Could not archive {record.output_dir}: {made.stderr.strip()}")
    target = Path(dest).expanduser().resolve() / record.id
    target.mkdir(parents=True, exist_ok=True)
    local_archive = target.parent / f"{record.id}-results.tgz"
    try:
        backend.download(archive, str(local_archive))
        _safe_extract(local_archive, target)
    finally:
        backend.run(f"rm -f {shlex.quote(archive)}")
        if local_archive.exists():
            local_archive.unlink()
    record.fetched_to = str(target)
    report = target / "report.html"
    return report if report.exists() else target
