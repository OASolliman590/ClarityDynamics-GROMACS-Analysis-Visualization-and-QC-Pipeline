"""Settings, parallel workflow, run tracking, input preparation and doctor (specs/003, US1-US6).

Remote flows run through ``ShellRunner`` (every file access is a shell command, as over SSH)
with fake ``condor_*``/``gmx``/``vmd``/``gromacs-pipeline`` executables first on PATH, so
the generated commands and scripts really execute.
"""

import posixpath
import time
from pathlib import Path

import pytest
import yaml
from setup_fixtures import add_prep_files, build_project, fake_bin

from gromacs_analysis.pipeline.cache import config_fingerprint
from gromacs_analysis.pipeline.yaml_loader import load_pipeline_config
from gromacs_analysis.setup import doctor, jobs, prepare, runner, settings, workflow
from gromacs_analysis.setup.backends import BackendError, LocalBackend, ShellRunner
from gromacs_analysis.setup.generator import refresh_stages, render_yaml, write_pipeline
from gromacs_analysis.setup.project import Project, ProjectError
from gromacs_analysis.setup.scanner import scan

SUBMIT_OK = 'echo "Submitting job(s)."; echo "1 job(s) submitted to cluster 777."'


@pytest.fixture
def project(tmp_path: Path) -> Project:
    root = add_prep_files(build_project(tmp_path / "proj"))
    proj = Project.from_scan(scan(LocalBackend(), str(root)))
    refresh_stages(proj)
    return proj


def wait_for(path: Path, timeout: float = 15.0) -> None:
    end = time.time() + timeout
    while time.time() < end:
        if path.exists() and path.read_text().strip():
            return
        time.sleep(0.05)
    raise AssertionError(f"{path} did not appear")


# ---------------------------------------------------------------- settings (US4)

def test_settings_set_validate_and_apply(project: Project, tmp_path: Path):
    settings.set_value(project, "prolif_batch.stride", "5")
    settings.set_value(project, "ranking.n_bootstrap", 500)
    settings.set_value(project, "distance_batch.residues.cutoff", "6")
    settings.set_value(project, "qc.max_frames", "none")
    settings.set_value(project, "md.convergence.neff_min", "40")
    assert project.settings["prolif_batch"] == {"stride": 5}
    out = write_pipeline(project, tmp_path / "p.yaml")
    data = yaml.safe_load(out.read_text())
    assert data["prolif_batch"]["stride"] == 5 and data["ranking"]["n_bootstrap"] == 500
    assert data["distance_batch"]["residues"] == {"auto": True, "cutoff": 6.0, "stride": 10}
    assert data["qc"]["max_frames"] is None
    cfg = load_pipeline_config(out)
    assert {c.stride for c in cfg.stages["prolif_batch"].config} == {5}
    assert cfg.stages["qc"].config.max_frames is None
    assert cfg.stages["md"].config.convergence == {"neff_min": 40.0}

    settings.set_value(project, "prolif_batch.stride", "10")  # back to default -> override removed
    assert "prolif_batch" not in project.settings


@pytest.mark.parametrize("name,value,message", [
    ("prolif_batch.stride", "0", ">= 1"),
    ("prolif_batch.stride", "five", "expects int"),
    ("ttclust_batch.method", "kmeans", "one of"),
    ("pca_joint.convergence.cosine_max", "1.5", "<= 1"),
    ("nope.key", "1", "Unknown setting"),
])
def test_settings_reject_invalid(project: Project, name, value, message):
    with pytest.raises(ProjectError, match=message):
        settings.set_value(project, name, value)


def test_settings_listing_and_assignment_parsing(project: Project):
    rows = settings.describe(project, "ranking")
    assert {r["name"] for r in rows} == {"ranking.n_bootstrap", "ranking.seed", "ranking.tests.n_permutations"}
    assert settings.parse_assignments(["md.smoothing_window=5"]) == [("md.smoothing_window", "5")]
    with pytest.raises(ProjectError):
        settings.parse_assignments(["stride=5"])


def test_settings_survive_save_and_load(project: Project, tmp_path: Path):
    settings.set_value(project, "pca_joint.stride", 20)
    project.save(tmp_path / "p.yaml")
    assert Project.load(tmp_path / "p.yaml").settings == {"pca_joint": {"stride": 20}}


# ---------------------------------------------------------------- workflow (US3)

def test_units_cover_every_entry_and_match_full_run(project: Project, tmp_path: Path):
    pipeline = yaml.safe_load(render_yaml(project))
    units = workflow.build_units(pipeline)
    names = [u.name for u in units]
    assert names[:4] == ["qc", "md", "mmpbsa", "pca_joint"]
    assert sum(u.stage == "prolif_batch" for u in units) == len(pipeline["prolif_batch"]["systems"])
    assert sum(u.stage == "distance_batch" for u in units) == len(pipeline["distance_batch"]["systems"])
    assert "prolif_compare" not in {u.stage for u in units} and "ranking" not in {u.stage for u in units}
    workflow.validate_units(units)

    full = tmp_path / "full.yaml"
    full.write_text(render_yaml(project))
    full_cfg = load_pipeline_config(full)
    for stage in ("prolif_batch", "distance_batch"):
        expected = {str(c.output_dir): config_fingerprint(c) for c in full_cfg.stages[stage].config}
        got = {}
        for unit in (u for u in units if u.stage == stage):
            path = tmp_path / f"{unit.name}.yaml"
            path.write_text(workflow.unit_yaml(unit))
            for c in load_pipeline_config(path).stages[stage].config:
                got[str(c.output_dir)] = config_fingerprint(c)
        assert got == expected  # same folders, same fingerprints -> the final job hits the cache


def test_dag_file_structure(project: Project):
    units = workflow.build_units(yaml.safe_load(render_yaml(project)))
    dag = workflow.render_dag(units, final_cpus=8)
    assert dag.count("JOB ") == len(units) + 1
    parent_line = next(line for line in dag.splitlines() if line.startswith("PARENT"))
    assert parent_line.endswith("CHILD final") and all(u.name in parent_line for u in units)
    assert 'VARS final config="pipeline.yaml" job="final" cpus="8"' in dag


def test_submit_dag_through_fake_pool(project: Project, tmp_path: Path, monkeypatch):
    fake_bin(tmp_path, {"gromacs-pipeline": "exit 0", "condor_submit_dag": SUBMIT_OK}, monkeypatch)
    result = workflow.submit_dag(ShellRunner(), project, render_yaml(project), cpus=6, memory="12GB")
    paths = workflow.dag_paths(project)
    assert result["cluster_id"] == 777 and result["record"].dag_cluster == 777
    assert Path(paths["dag"]).read_text().count("JOB ") == len(result["units"]) + 1
    sub = Path(paths["submit"]).read_text()
    assert "request_memory = 12288MB" in sub and "arguments = $(config)" in sub
    assert len(list(Path(paths["units"]).glob("*.yaml"))) == len(result["units"])
    script = Path(paths["script"]).read_text()
    assert "--continue-on-error" in script and 'config="$dag_dir/$1"' in script


# ---------------------------------------------------------------- runs (US1)

def test_host_run_status_and_exit_code(project: Project, tmp_path: Path, monkeypatch):
    fake_bin(tmp_path, {"gromacs-pipeline": 'if [ "$1" = "--version" ]; then exit 0; fi; '
                                             'echo "stage md"; sleep 1; echo "all done"; exit 3'}, monkeypatch)
    host = ShellRunner()
    started = runner.run_on_host(host, project, render_yaml(project))
    record = started["record"]
    jobs.add_run(project, record)
    assert jobs.status(host, record)["state"] == "running"
    wait_for(Path(record.exit_file))
    info = jobs.status(host, record, lines=5)
    assert info["state"] == "failed" and "exit code 3" in info["detail"]
    assert "all done" in info["log_tail"]
    project.save(tmp_path / "p.yaml")
    assert jobs.get_run(Project.load(tmp_path / "p.yaml")).id == record.id


def test_condor_and_dag_status_from_queue_and_history(tmp_path: Path, monkeypatch):
    state = tmp_path / "state"
    state.mkdir()
    fake_bin(tmp_path, {"condor_q": f'cat {state}/q.txt 2>/dev/null', "condor_history": f'cat {state}/h.txt 2>/dev/null'},
             monkeypatch)
    record = jobs.new_record("analysis", "dag", dag_cluster=777, clusters=[777])
    (state / "q.txt").write_text("777 0 2\n778 0 1\n779 0 2\n")
    info = jobs.status(ShellRunner(), record)
    assert info["state"] == "running" and info["counts"] == {"running": 2, "idle": 1}

    (state / "q.txt").write_text("")
    (state / "h.txt").write_text("777 0 4 0\n778 0 4 0\n779 0 4 1\n")
    info = jobs.status(ShellRunner(), record)
    assert info["state"] == "failed" and info["counts"] == {"completed": 2, "failed": 1}

    (state / "h.txt").write_text("777 0 4 0\n778 0 4 0\n")
    assert jobs.status(ShellRunner(), record)["state"] == "completed"
    (state / "h.txt").write_text("")
    assert jobs.status(ShellRunner(), record)["state"] == "unknown"


def test_fetch_downloads_and_unpacks_results(tmp_path: Path):
    out = tmp_path / "remote" / "analysis_outputs"
    (out / "md" / "plots").mkdir(parents=True)
    (out / "report.html").write_text("<html>report</html>")
    (out / "md" / "plots" / "rmsd.png").write_bytes(b"\x89PNG\x00\x01binary")
    record = jobs.new_record("analysis", "host", output_dir=str(out), run_dir=str(tmp_path / "remote" / "clarity"))
    host = ShellRunner()
    with pytest.raises(BackendError, match="not checked"):
        jobs.fetch(host, record, tmp_path / "local")
    record.last_state = "completed"
    (tmp_path / "remote" / "clarity").mkdir()
    report = jobs.fetch(host, record, tmp_path / "local")
    assert report.read_text() == "<html>report</html>"
    assert (report.parent / "md" / "plots" / "rmsd.png").read_bytes() == b"\x89PNG\x00\x01binary"
    assert not list((tmp_path / "remote" / "clarity").glob("*.tgz"))  # remote archive cleaned up


# ---------------------------------------------------------------- prepare (US2)

def test_plan_reports_ready_and_skipped_tasks(project: Project):
    tasks = prepare.plan_preparation(project)
    by = {(t.replicate, t.kind): t for t in tasks}
    assert by[("HOLO_LigB_rep1", "trjconv")].ok and by[("HOLO_LigB_rep1", "trjconv")].outputs == ["md_noPBC.xtc"]
    assert by[("HOLO_LigB_rep1", "metrics")].ok and by[("HOLO_LigB_rep1", "mmpbsa")].ok
    assert ("HOLO_LigA_1", "trjconv") not in by           # already processed
    assert ("HOLO_LigA_1", "mmpbsa") not in by            # results exist
    assert "no .tpr" in by[("HOLO_LigA_2", "mmpbsa")].skipped_reason
    assert by[("APO_1", "trjconv")].ok and "no .tpr" in by[("APO_2", "trjconv")].skipped_reason
    assert ("APO_1", "mmpbsa") not in by                  # APO has nothing to bind
    assert by[("LigC/rep1", "metrics")].ok                # structure + dcd are enough for VMD
    assert prepare.summarise(tasks)


def test_scripts_use_named_groups_and_generated_inputs(project: Project):
    scripts = prepare.build_scripts(project, prepare.plan_preparation(project, mpi_ranks=4))
    ligb = next(s for s in scripts if s.replicate == "HOLO_LigB_rep1")
    text = ligb.text
    assert "printf 'Protein\\nSystem\\n'" in text and "-pbc mol -center" in text and "-fit rot+trans" in text
    assert "resname LGB" in text and "-cg 0 1" in text and "-cp topol.top" in text and "mpirun -np 4" in text
    assert text.index("trjconv") < text.index("vmd") < text.index("gmx_MMPBSA")
    files = {posixpath.basename(p): t for p, t in ligb.extra_files.items()}
    assert 'set trajectory "md_noPBC.xtc"' in files["clarity_md_metrics.tcl"]
    assert "interaction_entropy=1" in files["clarity_mmpbsa.in"]


def test_preparation_runs_and_rescan_sees_new_inputs(project: Project, tmp_path: Path, monkeypatch):
    touch_o = 'out=""; while [ $# -gt 0 ]; do if [ "$1" = "-o" ] || [ "$1" = "-on" ]; then out="$2"; fi; shift; done; ' \
              'cat > /dev/null; [ -n "$out" ] && echo traj > "$out"; exit 0'
    fake_bin(tmp_path, {
        "gmx": touch_o,
        "vmd": "for f in RMSD_apo_.dat RMSF_.dat rog_.dat; do echo '0 1' > $f; done",
        "gmx_MMPBSA": "echo results > FINAL_RESULTS_MMPBSA.dat",
    }, monkeypatch)
    host = ShellRunner()
    scripts = prepare.build_scripts(project, prepare.plan_preparation(project, what=["trjconv", "metrics", "mmpbsa"]))
    result = prepare.run_preparation(host, project, scripts, where="host")
    wait_for(Path(result["record"].exit_file))
    assert Path(result["record"].exit_file).read_text().strip() == "0"
    project.merge_scan(scan(LocalBackend(), project.root))
    ligb = project.system("HOLO_LigB")
    assert ligb.inputs["processed_trajectory"] and ligb.inputs["md_metrics"] and ligb.inputs["mmpbsa"]
    assert ligb.replicates[0].trajectory.endswith("md_noPBC.xtc")


def test_preparation_as_condor_jobs(project: Project, tmp_path: Path, monkeypatch):
    fake_bin(tmp_path, {"condor_submit": SUBMIT_OK}, monkeypatch)
    scripts = prepare.build_scripts(project, prepare.plan_preparation(project, what=["trjconv"]))
    result = prepare.run_preparation(ShellRunner(), project, scripts, where="condor", cpus=2)
    sub = Path(prepare.prep_dir(project), "prepare.sub").read_text()
    assert result["cluster_id"] == 777 and "queue script from" in sub and "request_cpus = 2" in sub
    listed = Path(prepare.prep_dir(project), "scripts.txt").read_text().split()
    assert listed == [s.path for s in scripts]
    with pytest.raises(BackendError, match="Nothing to prepare"):
        prepare.run_preparation(ShellRunner(), project, [], where="condor")


# ---------------------------------------------------------------- doctor (US6)

def test_doctor_reports_missing_tools_with_remedies(project: Project, tmp_path: Path, monkeypatch):
    import os
    import sys

    bin_dir = fake_bin(tmp_path, {"gromacs-pipeline": "exit 0", "gmx": "exit 0", "condor_submit": "exit 0",
                                  "condor_q": "exit 0"}, monkeypatch)
    monkeypatch.setenv("PATH", os.pathsep.join([str(bin_dir), os.path.dirname(sys.executable), "/usr/bin", "/bin"]))
    checks = {c.name: c for c in doctor.run_checks(ShellRunner(), project.root, modes=["host"])}
    assert checks["connection"].ok and checks["root"].ok and checks["pipeline"].ok and checks["gmx"].ok
    assert not checks["vmd"].ok and not checks["vmd"].required and "VMD" in checks["vmd"].remedy
    assert doctor.exit_code(list(checks.values())) == 0
    dag = doctor.run_checks(ShellRunner(), project.root, modes=["dag", "metrics"])
    missing = {c.name for c in dag if c.required and not c.ok}
    assert missing == {"condor_submit_dag", "vmd"} and doctor.exit_code(dag) == 4
    assert any("shared" in line or "fs" in line for line in doctor.format_table(dag))
    bad = doctor.run_checks(ShellRunner(), str(tmp_path / "nope"))
    assert not next(c for c in bad if c.name == "root").ok


# ---------------------------------------------------------------- CLI (contracts/cli.md)

def _cli_project(tmp_path: Path, kind: str = "shell") -> Path:
    from gromacs_analysis.setup import cli

    root = add_prep_files(build_project(tmp_path / "proj"))
    proj_file = tmp_path / "clarity_project.yaml"
    assert cli.main(["init", str(root), "--project", str(proj_file), "--out", str(tmp_path / "pipeline.yaml")]) == 0
    project = Project.load(proj_file)
    project.connection.kind = kind  # commands go through the shell, as over SSH
    project.save(proj_file)
    return proj_file


def test_cli_settings_and_init_set(tmp_path: Path, capsys):
    from gromacs_analysis.setup import cli

    proj_file = _cli_project(tmp_path)
    assert cli.main(["settings", "--project", str(proj_file), "--stage", "prolif_batch",
                     "--set", "prolif_batch.stride=4"]) == 0
    out = capsys.readouterr().out
    assert "* prolif_batch.stride" in out and "4" in out
    assert cli.main(["settings", "--project", str(proj_file), "--set", "prolif_batch.stride=0"]) == 1
    assert ">= 1" in capsys.readouterr().err
    assert cli.main(["init", "--project", str(proj_file), "--out", str(tmp_path / "pipeline.yaml"),
                     "--set", "ranking.n_bootstrap=100"]) == 0
    data = yaml.safe_load((tmp_path / "pipeline.yaml").read_text())
    assert data["prolif_batch"]["stride"] == 4 and data["ranking"]["n_bootstrap"] == 100


def test_cli_prepare_dry_run_writes_nothing(tmp_path: Path, capsys):
    from gromacs_analysis.setup import cli

    proj_file = _cli_project(tmp_path)
    assert cli.main(["prepare", "--project", str(proj_file), "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "HOLO_LigB_rep1" in out and "skipped  (no .tpr" in out and "Dry run" in out
    assert not (tmp_path / "proj" / "clarity" / "prep").exists()


def test_cli_dag_status_and_fetch(tmp_path: Path, capsys, monkeypatch):
    from gromacs_analysis.setup import cli

    state = tmp_path / "state"
    state.mkdir()
    fake_bin(tmp_path, {
        "gromacs-pipeline": "exit 0", "condor_submit_dag": SUBMIT_OK,
        "condor_q": f"cat {state}/q.txt 2>/dev/null", "condor_history": f"cat {state}/h.txt 2>/dev/null",
    }, monkeypatch)
    proj_file = _cli_project(tmp_path)
    assert cli.main(["run", "--project", str(proj_file), "--where", "dag", "--unit-cpus", "2"]) == 0
    out = capsys.readouterr().out
    assert "DAG cluster 777" in out and "Run id: analysis-" in out
    record = jobs.get_run(Project.load(proj_file))
    assert record.where == "dag" and record.dag_cluster == 777

    (state / "q.txt").write_text("777 0 2\n778 0 2\n")
    assert cli.main(["status", "--project", str(proj_file)]) == 0
    assert "RUNNING" in capsys.readouterr().out
    assert cli.main(["fetch", "--project", str(proj_file), "--to", str(tmp_path / "local")]) == 1
    assert "running" in capsys.readouterr().err

    (state / "q.txt").write_text("")
    (state / "h.txt").write_text("777 0 4 0\n778 0 4 0\n")
    out_root = Path(record.output_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "report.html").write_text("<html>ok</html>")
    assert cli.main(["fetch", "--project", str(proj_file), "--to", str(tmp_path / "local")]) == 0
    assert (tmp_path / "local" / record.id / "report.html").read_text() == "<html>ok</html>"
    assert Project.load(proj_file).runs[-1]["fetched_to"]


def test_cli_doctor_exit_codes(tmp_path: Path, capsys, monkeypatch):
    import os
    import sys

    from gromacs_analysis.setup import cli

    bin_dir = fake_bin(tmp_path, {"gromacs-pipeline": "exit 0"}, monkeypatch)
    monkeypatch.setenv("PATH", os.pathsep.join([str(bin_dir), os.path.dirname(sys.executable), "/usr/bin", "/bin"]))
    proj_file = _cli_project(tmp_path)
    assert cli.main(["doctor", "--project", str(proj_file)]) == 0
    assert "All required checks passed" in capsys.readouterr().out
    assert cli.main(["doctor", "--project", str(proj_file), "--where", "dag"]) == 4
    assert "condor_submit_dag" in capsys.readouterr().out


def test_cli_status_without_runs(tmp_path: Path, capsys):
    from gromacs_analysis.setup import cli

    proj_file = _cli_project(tmp_path)
    assert cli.main(["status", "--project", str(proj_file)]) == 0
    assert "No runs recorded" in capsys.readouterr().out
