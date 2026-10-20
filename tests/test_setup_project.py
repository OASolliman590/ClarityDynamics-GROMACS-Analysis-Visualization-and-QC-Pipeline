"""Project edits, configuration generation, runners and CLI (specs/002-pipeline-setup-tui, US2/US3/US5)."""

import json
from pathlib import Path

import pytest
import yaml
from setup_fixtures import build_project

from gromacs_analysis.setup import cli, runner
from gromacs_analysis.setup.backends import BackendError, CommandResult, LocalBackend, ShellRunner
from gromacs_analysis.setup.generator import (
    build_pipeline_dict,
    enabled_stages,
    refresh_stages,
    render_yaml,
    stage_availability,
    validate_yaml,
    write_pipeline,
)
from gromacs_analysis.setup.project import Project, ProjectError, sanitize_name
from gromacs_analysis.setup.scanner import scan


@pytest.fixture
def project(tmp_path: Path) -> Project:
    root = build_project(tmp_path / "proj")
    proj = Project.from_scan(scan(LocalBackend(), str(root)))
    refresh_stages(proj)
    return proj


# ---------------------------------------------------------------- US2 naming

def test_rename_system_and_uniqueness(project: Project):
    project.rename_system("HOLO_LigA", "Isoniazid complex")
    assert project.system("HOLO_LigA").name == "Isoniazid_complex"
    assert project.system("Isoniazid_complex").key == "HOLO_LigA"
    with pytest.raises(ProjectError, match="already used"):
        project.rename_system("HOLO_LigB", "Isoniazid_complex")
    with pytest.raises(ProjectError):
        project.rename_system("HOLO_LigB", "  ")
    with pytest.raises(ProjectError, match="No system"):
        project.rename_system("nope", "x")


def test_excluded_system_frees_its_name(project: Project):
    project.set_included("HOLO_LigA", False)
    name = project.system("HOLO_LigA").name
    project.rename_system("HOLO_LigB", name)  # allowed: the clash is excluded
    with pytest.raises(ProjectError, match="Cannot include"):
        project.set_included("HOLO_LigA", True)


def test_rename_ligand_everywhere_or_per_system(project: Project):
    assert project.rename_ligand("INH", "Isoniazid") == 1
    lig = next(lg for lg in project.system("HOLO_LigA").ligands if lg.resname == "INH")
    assert lig.name == "Isoniazid" and lig.selection == "resname INH"
    project.rename_ligand("LGC", "Compound C", system="LigC")
    assert project.system("LigC").ligands[0].name == "Compound_C"
    with pytest.raises(ProjectError, match="already used"):
        project.rename_ligand("NAD", "Isoniazid", system="HOLO_LigA")
    with pytest.raises(ProjectError, match="No ligand"):
        project.rename_ligand("XXX", "y")


def test_sanitize_name():
    assert sanitize_name(" My  ligand/(2) ") == "My_ligand2"


def test_save_load_round_trip(project: Project, tmp_path: Path):
    project.rename_system("APO", "Apo_enzyme")
    project.rename_ligand("NAD", "Cofactor")
    project.set_included("LigC", False)
    path = project.save(tmp_path / "clarity_project.yaml")
    assert path.read_text().startswith("# ClarityDynamics project file")
    loaded = Project.load(path)
    assert loaded.to_dict() == project.to_dict()


def test_load_rejects_foreign_yaml(tmp_path: Path):
    p = tmp_path / "x.yaml"
    p.write_text("run: [md]\n")
    with pytest.raises(ProjectError):
        Project.load(p)


def test_rescan_merge_keeps_edits_adds_and_flags_missing(project: Project, tmp_path: Path):
    import shutil

    root = Path(project.root)
    project.rename_system("HOLO_LigA", "Isoniazid")
    project.rename_ligand("LGB", "Ligand_B")
    shutil.rmtree(root / "LigC")
    shutil.copytree(root / "HOLO_LigB_rep1", root / "NEWSYS_1")
    project.merge_scan(scan(LocalBackend(), str(root)))
    assert project.system("HOLO_LigA").name == "Isoniazid"
    assert project.system("HOLO_LigB").ligands[0].name == "Ligand_B"
    assert project.system("LigC").missing
    assert "LigC" not in [s.key for s in project.active_systems]
    assert project.system("NEWSYS").name == "NEWSYS"


# ---------------------------------------------------------------- US3 configuration

def test_stage_availability_reasons(project: Project):
    avail = stage_availability(project)
    assert all(avail[s].available for s in ("qc", "md", "mmpbsa", "pca_joint", "prolif_batch", "ranking"))
    assert avail["ttclust_batch"].available and not avail["ttclust_batch"].enabled  # slow: opt-in
    for s in project.systems:
        for rep in s.replicates:
            rep.inputs["mmpbsa"] = False
    avail = stage_availability(project)
    assert not avail["mmpbsa"].available and "FINAL_RESULTS_MMPBSA" in avail["mmpbsa"].reason


def test_refresh_keeps_user_choice(project: Project):
    project.stages["pca_joint"].enabled = False
    project.stages["ttclust_batch"].enabled = True
    refresh_stages(project)
    assert not project.stages["pca_joint"].enabled and project.stages["ttclust_batch"].enabled


def test_generated_yaml_validates_with_names_and_selections(project: Project):
    project.rename_system("HOLO_LigA", "Isoniazid_cplx")
    project.rename_ligand("INH", "Isoniazid")
    project.stages["ttclust_batch"].enabled = True
    text = render_yaml(project, project_file="clarity_project.yaml")
    validate_yaml(text)
    data = yaml.safe_load(text)
    assert data["run"] == enabled_stages(project)

    md_names = [s["name"] for s in data["md"]["systems"]]
    assert "Isoniazid_cplx" in md_names
    assert next(s for s in data["md"]["systems"] if s["name"] == "APO")["dir_names"] == [
        "APO_1/analysis", "APO_2/analysis", "APO_3/analysis"]
    assert "Isoniazid_cplx" in [s["name"] for s in data["pca_joint"]["systems"]]
    assert "Isoniazid_cplx" in [s["name"] for s in data["ttclust_batch"]["systems"]]

    # One ProLIF entry per ligand run (independent work units), named <system>_<ligand>.
    prolif = {s["name"]: s for s in data["prolif_batch"]["systems"]}
    assert "APO" not in prolif
    assert prolif["Isoniazid_cplx_Isoniazid"]["ligands"] == [{"name": "Isoniazid", "selection": "resname INH"}]
    assert prolif["Isoniazid_cplx_NAD"]["ligands"] == [{"name": "NAD", "selection": "resname NAD"}]
    assert prolif["LigC"]["ligands"] == [{"name": "LGC", "selection": "resname LGC"}]
    assert all(t.startswith(project.root) for t in prolif["Isoniazid_cplx_NAD"]["trajectories"])
    assert all(data[stage]["cache"] is True for stage in data["run"])

    dist_names = [s["name"] for s in data["distance_batch"]["systems"]]
    assert "Isoniazid_cplx_Isoniazid" in dist_names and "Isoniazid_cplx_NAD" in dist_names

    assert data["ranking"]["exclude_systems"] == ["APO"]
    mm = data["mmpbsa"]["systems"]
    assert [s["name"] for s in mm] == ["Isoniazid_cplx"]
    assert mm[0]["dir_pattern"] == "HOLO_LigA_{}" and mm[0]["replicates"] == 3


def test_excluded_systems_and_ligands_are_absent(project: Project):
    project.set_included("LigC", False)
    project.set_ligand_included("NAD", False)
    data = build_pipeline_dict(project)
    names = [s["name"] for s in data["prolif_batch"]["systems"]]
    assert "LigC" not in names
    liga = next(s for s in data["prolif_batch"]["systems"] if s["name"] == "LigA")
    assert [lg["selection"] for lg in liga["ligands"]] == ["resname INH"]


def test_generation_refuses_invalid_project(project: Project, tmp_path: Path):
    project.system("HOLO_LigB").name = project.system("HOLO_LigA").name  # bypass edit checks
    with pytest.raises(ValueError, match="Duplicate system name"):
        write_pipeline(project, tmp_path / "pipeline.yaml")
    assert not (tmp_path / "pipeline.yaml").exists()


# ---------------------------------------------------------------- US5 run

class FakeHost(ShellRunner):
    """ShellRunner that fakes the pipeline executable and condor_submit."""

    def __init__(self, installed=True):
        super().__init__()
        self.installed = installed
        self.log = []

    def run(self, command, timeout=None):
        if command.startswith("command -v"):
            self.log.append(command)
            return CommandResult(0 if self.installed else 1, "/usr/bin/gromacs-pipeline\n" if self.installed else "", "")
        if "condor_submit" in command:
            self.log.append(command)
            return CommandResult(0, "1 job(s) submitted to cluster 4242.\n", "")
        if "nohup" in command:
            self.log.append(command)
            return CommandResult(0, "31337\n", "")
        return super().run(command, timeout)


def test_run_on_host_writes_config_and_starts_once(project: Project):
    host = FakeHost()
    result = runner.run_on_host(host, project, render_yaml(project))
    paths = runner.remote_paths(project)
    assert Path(paths["config"]).read_text().startswith("# Generated by clarity-setup")
    assert result["pid"] == "31337" and result["log"] == paths["log"]
    assert sum("nohup" in c for c in host.log) == 1


def test_run_on_host_reports_missing_install(project: Project):
    with pytest.raises(BackendError, match="not installed"):
        runner.run_on_host(FakeHost(installed=False), project, "run: []\n")


def test_submit_condor_builds_files_and_parses_cluster(project: Project):
    host = FakeHost()
    result = runner.submit_condor(host, project, render_yaml(project), cpus=8, memory="16GB")
    paths = runner.remote_paths(project)
    assert result["cluster_id"] == 4242
    sub = Path(paths["submit"]).read_text()
    assert "request_cpus = 8" in sub and "request_memory = 16384MB" in sub
    assert "should_transfer_files = NO" in sub and sub.rstrip().endswith("queue 1")
    script = Path(paths["script"]).read_text()
    assert script.startswith("#!/bin/bash") and "--config" in script


def test_submit_files_carry_pool_policy(project: Project):
    files = runner.build_submit_files(project, extra_attributes={"Production": True, "Group": "lab"},
                                      required_arch="X86_64")
    assert "+Production = True" in files["submit"] and '+Group = "lab"' in files["submit"]
    assert 'TARGET.Arch == "X86_64"' in files["submit"]


@pytest.mark.parametrize("text,mb", [("8GB", 8192), ("4096MB", 4096), ("2g", 2048), ("512", 512)])
def test_memory_parsing(text, mb):
    assert runner._memory_mb(text) == mb


def test_run_local_invokes_pipeline(project: Project, tmp_path: Path, monkeypatch):
    seen = {}

    def fake_run(config, continue_on_error=True):
        seen["stages"] = list(config.run)
        return {"md": {"success": True}, "qc": {"success": False}, "_report": "report.html"}

    import gromacs_analysis.pipeline as pipeline_mod

    monkeypatch.setattr(pipeline_mod, "run_pipeline", fake_run)
    out = write_pipeline(project, tmp_path / "pipeline.yaml")
    result = runner.run_local(out)
    assert seen["stages"] and result["failed_stages"] == ["qc"] and result["report"] == "report.html"


# ---------------------------------------------------------------- CLI

def test_cli_scan_json(project: Project, capsys):
    assert cli.main(["scan", project.root, "--json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert {s["key"] for s in data["systems"]} == {"HOLO_LigA", "HOLO_LigB", "APO", "LigC"}


def test_cli_scan_empty_exits_2(tmp_path: Path):
    assert cli.main(["scan", str(tmp_path)]) == 2


def test_cli_init_applies_edits_and_writes_valid_files(project: Project, tmp_path: Path, capsys):
    proj_file, out = tmp_path / "p.yaml", tmp_path / "pipeline.yaml"
    code = cli.main(["init", project.root, "--project", str(proj_file), "--out", str(out),
                     "--rename", "HOLO_LigA=Isoniazid_cplx", "--ligand", "INH=Isoniazid",
                     "--ligand", "LigC:LGC=CompoundC", "--exclude", "HOLO_LigB",
                     "--protein", "InhA", "--stages", "qc,md,prolif_batch"])
    assert code == 0, capsys.readouterr().err
    data = yaml.safe_load(out.read_text())
    assert data["run"] == ["qc", "md", "prolif_batch"]
    assert data["md"]["protein_name"] == "InhA"
    names = [s["name"] for s in data["prolif_batch"]["systems"]]
    assert names == ["Isoniazid_cplx_Isoniazid", "Isoniazid_cplx_NAD", "LigC"]
    saved = Project.load(proj_file)
    assert saved.system("LigC").ligands[0].name == "CompoundC"

    # Re-running init without ROOT reuses the project file and keeps the edits.
    assert cli.main(["init", "--project", str(proj_file), "--out", str(out)]) == 0
    assert Project.load(proj_file).system("HOLO_LigA").name == "Isoniazid_cplx"


def test_cli_init_reports_errors(project: Project, tmp_path: Path, capsys):
    code = cli.main(["init", project.root, "--project", str(tmp_path / "p.yaml"),
                     "--out", str(tmp_path / "o.yaml"), "--rename", "HOLO_LigA=LigB"])
    assert code == 1 and "already used" in capsys.readouterr().err
    code = cli.main(["init", project.root, "--project", str(tmp_path / "p.yaml"),
                     "--out", str(tmp_path / "o.yaml"), "--stages", "md,bogus"])
    assert code == 1 and "Unknown stages" in capsys.readouterr().err


def test_cli_tui_without_textual_exits_3(monkeypatch, capsys):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "tui" or name.endswith(".tui") or (level and fromlist and "run_app" in fromlist):
            raise ImportError("No module named 'textual'")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    assert cli.main(["tui"]) == 3
    assert "pip install" in capsys.readouterr().err


def test_prolif_split_keeps_loader_output_folders(project: Project, tmp_path: Path):
    """Per-ligand entries produce the same output folders the loader made for multi-ligand entries."""
    from gromacs_analysis.pipeline.yaml_loader import load_pipeline_config

    out = write_pipeline(project, tmp_path / "pipeline.yaml")
    cfgs = load_pipeline_config(out).stages["prolif_batch"].config
    folders = sorted(c.output_dir.name for c in cfgs)
    assert folders == ["output_LigA_INH", "output_LigA_NAD", "output_LigB", "output_LigC"]
    assert {c.ligand_selection for c in cfgs} == {"resname INH", "resname NAD", "resname LGB", "resname LGC"}
