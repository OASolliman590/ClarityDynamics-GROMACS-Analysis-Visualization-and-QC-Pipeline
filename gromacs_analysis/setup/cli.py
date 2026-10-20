"""``clarity-setup``: discover systems, name them, and generate/run the pipeline without YAML.

    clarity-setup scan ROOT [--ssh HOST | --condorctl CONFIG]
    clarity-setup init ROOT [--rename OLD=NEW] [--ligand RES=NAME] [--apo NAME] [--exclude NAME] ...
    clarity-setup run  [--project clarity_project.yaml] [--where local|host|condor]
    clarity-setup tui  [ROOT] [--ssh HOST | --condorctl CONFIG]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional

from .backends import BackendError, connect
from .generator import STAGES, refresh_stages, render_yaml, validate_yaml, write_pipeline
from .project import PROJECT_FILE, Connection, Project, ProjectError
from .scanner import scan


def _connection(args) -> Connection:
    if getattr(args, "ssh", None):
        return Connection(kind="ssh", host=args.ssh)
    if getattr(args, "condorctl", None):
        return Connection(kind="condorctl", profile=args.condorctl)
    return Connection()


def _add_connection_args(parser: argparse.ArgumentParser) -> None:
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--ssh", metavar="HOST", help="SSH host or alias from ~/.ssh/config (key auth)")
    group.add_argument("--condorctl", metavar="CONFIG", help="condorctl config/profile (e.g. nmrbox.yaml)")
    parser.add_argument("--depth", type=int, default=4, help="Directory depth to search (default 4)")


def _resolve_root(root: str, connection: Connection) -> str:
    return root if connection.is_remote else str(Path(root).expanduser().resolve())


def _scan_project(root: str, connection: Connection, depth: int, existing: Optional[Project] = None) -> Project:
    with connect(connection.kind, connection.host, connection.profile) as backend:
        result = scan(backend, _resolve_root(root, connection), depth=depth)
    if existing is not None:
        existing.connection = connection
        existing.merge_scan(result)
        project = existing
    else:
        project = Project.from_scan(result, connection=connection)
    project._last_scan = result  # noqa: SLF001 - transient, for reporting only
    refresh_stages(project)
    return project


def _print_table(project: Project) -> None:
    rows = [("System", "Name", "Reps", "APO", "Ligands", "Inputs", "Topology")]
    for s in project.systems:
        inputs = ",".join(k for k, v in s.inputs.items() if v) or "-"
        ligands = ",".join(f"{l.resname}" + (f"→{l.name}" if l.name != l.resname else "") for l in s.ligands) or "-"
        flag = "" if s.included else " (excluded)"
        rows.append((s.key, s.name + flag, str(len(s.replicates)), "yes" if s.apo else "",
                     ligands, inputs, s.topology or "-"))
    widths = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    for i, row in enumerate(rows):
        print("  ".join(cell.ljust(widths[j]) for j, cell in enumerate(row)).rstrip())
        if i == 0:
            print("  ".join("-" * w for w in widths))


def _apply_edits(project: Project, args) -> None:
    if getattr(args, "protein", None):
        project.protein_name = args.protein
    for item in args.exclude or []:
        project.set_included(item, False)
    for item in args.apo or []:
        project.set_apo(item, True)
    for item in args.rename or []:
        old, _, new = item.partition("=")
        if not new:
            raise ProjectError(f"--rename needs OLD=NEW, got {item!r}")
        project.rename_system(old, new)
    for item in args.ligand or []:
        key, _, new = item.partition("=")
        if not new:
            raise ProjectError(f"--ligand needs RES=NAME or SYSTEM:RES=NAME, got {item!r}")
        system, _, res = key.rpartition(":")
        project.rename_ligand(res, new, system=system or None)
    refresh_stages(project)
    if getattr(args, "stages", None):
        wanted = {s.strip() for s in args.stages.split(",") if s.strip()}
        unknown = wanted - {name for name, _ in STAGES}
        if unknown:
            raise ProjectError(f"Unknown stages: {', '.join(sorted(unknown))}")
        for name, choice in project.stages.items():
            choice.enabled = name in wanted and choice.available


# ---------------------------------------------------------------- commands

def cmd_scan(args) -> int:
    project = _scan_project(args.root, _connection(args), args.depth)
    result = project._last_scan  # noqa: SLF001
    if args.json:
        print(json.dumps({
            "root": project.root,
            "backend": result.backend,
            "remote_commands": result.commands_used,
            "skipped_dirs": result.skipped_dirs,
            "systems": [
                {"key": s.key, "name": s.name, "apo": s.apo, "replicates": [r.dir for r in s.replicates],
                 "topology": s.topology, "trajectories": s.trajectories,
                 "ligands": [l.resname for l in s.ligands], "inputs": s.inputs}
                for s in project.systems
            ],
        }, indent=2))
    else:
        print(f"Scanned {project.root} on {result.backend} ({result.commands_used} command(s))\n")
        _print_table(project)
        if result.skipped_dirs:
            print(f"\nSkipped {len(result.skipped_dirs)} folder(s) without trajectories: "
                  + ", ".join(result.skipped_dirs[:8]) + (" ..." if len(result.skipped_dirs) > 8 else ""))
    return 0 if project.systems else 2


def cmd_init(args) -> int:
    project_path = Path(args.project)
    existing = Project.load(project_path) if project_path.exists() else None
    root = args.root or (existing.root if existing else None)
    if not root:
        print("error: give ROOT or an existing --project file", file=sys.stderr)
        return 1
    connection = _connection(args) if (args.ssh or args.condorctl or not existing) else existing.connection
    project = _scan_project(root, connection, args.depth, existing)
    _apply_edits(project, args)
    project.save(project_path)
    out = write_pipeline(project, Path(args.out), project_file=str(project_path))
    _print_table(project)
    enabled = [n for n, c in project.stages.items() if c.enabled]
    print(f"\nStages: {', '.join(enabled) or 'none'}")
    for name, choice in project.stages.items():
        if not choice.available:
            print(f"  - {name}: unavailable ({choice.reason})")
    print(f"\nProject file:  {project_path}\nPipeline YAML: {out} (validated)")
    return 0


def cmd_run(args) -> int:
    from . import runner  # noqa: WPS433

    project = Project.load(Path(args.project))
    refresh_stages(project) if not project.stages else None
    text = render_yaml(project, project_file=args.project)
    validate_yaml(text)
    where = args.where or ("host" if project.connection.is_remote else "local")
    if where == "local":
        if project.connection.is_remote:
            print("error: data is remote; use --where host or --where condor", file=sys.stderr)
            return 1
        out = Path(args.out)
        out.write_text(text)
        result = runner.run_local(out)
        print(f"Run finished. Report: {result['report']}")
        if result["failed_stages"]:
            print(f"Failed stages: {', '.join(result['failed_stages'])}")
            return 1
        return 0
    with connect(project.connection.kind, project.connection.host, project.connection.profile) as backend:
        if where == "host":
            result = runner.run_on_host(backend, project, text)
            print(f"Started on {backend.label} (pid {result['pid']}). Log: {result['log']}")
        else:
            result = runner.submit_condor(backend, project, text, cpus=args.cpus, memory=args.memory)
            print(f"Submitted cluster {result['cluster_id']}. Output: {result['log']}")
    print(f"Report when finished: {result['report']}")
    return 0


def cmd_tui(args) -> int:
    try:
        from .tui import run_app  # noqa: WPS433
    except ImportError as exc:
        print(f'The interface needs Textual: pip install "gromacs-analysis-toolkit[tui]" ({exc})',
              file=sys.stderr)
        return 3
    project_path = Path(args.project)
    project = Project.load(project_path) if project_path.exists() and not args.root else None
    return run_app(project=project, root=args.root, connection=_connection(args), depth=args.depth,
                   project_path=project_path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="clarity-setup",
        description="Discover simulations, name systems and ligands, and generate/run the pipeline.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("scan", help="List the systems found under ROOT")
    p.add_argument("root")
    _add_connection_args(p)
    p.add_argument("--json", action="store_true", help="Machine-readable output")
    p.set_defaults(func=cmd_scan)

    p = sub.add_parser("init", help="Scan, apply names/choices, write project file + pipeline YAML")
    p.add_argument("root", nargs="?")
    _add_connection_args(p)
    p.add_argument("--protein", help="Protein name used in titles")
    p.add_argument("--rename", action="append", metavar="OLD=NEW", help="Rename a system (key or name)")
    p.add_argument("--ligand", action="append", metavar="RES=NAME", help="Rename a ligand (RES or SYSTEM:RES)")
    p.add_argument("--apo", action="append", metavar="NAME", help="Mark a system as the APO reference")
    p.add_argument("--exclude", action="append", metavar="NAME", help="Leave a system out")
    p.add_argument("--stages", help="Comma-separated stages to enable (default: all available)")
    p.add_argument("--project", default=PROJECT_FILE, help=f"Project file (default {PROJECT_FILE})")
    p.add_argument("--out", default="pipeline.yaml", help="Pipeline YAML to write (default pipeline.yaml)")
    p.set_defaults(func=cmd_init)

    p = sub.add_parser("run", help="Run the project's pipeline locally, on the host, or via HTCondor")
    p.add_argument("--project", default=PROJECT_FILE)
    p.add_argument("--out", default="pipeline.yaml", help="Local copy of the generated YAML")
    p.add_argument("--where", choices=["local", "host", "condor"])
    p.add_argument("--cpus", type=int, default=4)
    p.add_argument("--memory", default="8GB")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("tui", help="Open the interactive terminal interface")
    p.add_argument("root", nargs="?")
    _add_connection_args(p)
    p.add_argument("--project", default=PROJECT_FILE)
    p.set_defaults(func=cmd_tui)
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return int(args.func(args) or 0)
    except (BackendError, ProjectError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
