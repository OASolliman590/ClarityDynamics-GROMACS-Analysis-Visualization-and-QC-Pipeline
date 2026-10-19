"""Round-3 tests (specs/001-pipeline-hardening-round3)."""

import logging
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------- US2 CI

def test_ci_workflow_runs_pytest_on_pull_requests():
    workflow = yaml.safe_load((ROOT / ".github" / "workflows" / "tests.yml").read_text())
    triggers = workflow.get(True) or workflow.get("on")  # YAML 1.1 parses "on" as True
    assert "pull_request" in triggers
    assert triggers["push"]["branches"] == ["main"]
    job = workflow["jobs"]["pytest"]
    assert {"3.10", "3.12"} <= set(job["strategy"]["matrix"]["python-version"])
    assert any("pytest" in str(step.get("run", "")) for step in job["steps"])


# ---------------------------------------------------------------- US9 static export availability

@pytest.fixture
def broken_exporter(monkeypatch):
    import plotly.graph_objects as go

    from gromacs_analysis.utils import figure_export

    calls = {"n": 0}

    def fail(self, *args, **kwargs):
        calls["n"] += 1
        raise RuntimeError("Kaleido requires Google Chrome to be installed. Run plotly_get_chrome.")

    figure_export.reset_static_export_state()
    monkeypatch.setattr(go.Figure, "write_image", fail)
    yield calls
    figure_export.reset_static_export_state()


def test_unavailable_export_warns_once_and_falls_back(tmp_path: Path, broken_exporter, caplog):
    import numpy as np
    import plotly.graph_objects as go

    from gromacs_analysis.md import AnalysisMetric, MDConfig, MDPlotter, SystemConfig
    from gromacs_analysis.utils.figure_export import REMEDY

    cfg = MDConfig(base_dir=tmp_path, output_dir=tmp_path / "out", protein_name="P",
                   systems=[SystemConfig(name="LigA")])
    cfg.plot_config.save_formats = ["html", "svg", "png"]
    plotter = MDPlotter(cfg)
    metric = AnalysisMetric(name="rmsd_prot", file_pattern="", title="RMSD", ylabel="Å")
    data = {"LigA": {"mode": "single", "values": np.ones(5), "time": np.arange(5)}}
    fallback = {"kind": "metric", "metric": metric, "metric_data": data, "x_values": None}

    with caplog.at_level(logging.WARNING):
        for i in range(10):
            plotter.save_plot(go.Figure(go.Scatter(y=[1, 2])), tmp_path / "plots", f"fig{i}", fallback=fallback)

    export_warnings = [r for r in caplog.records if REMEDY[:40] in r.getMessage()]
    errors = [r for r in caplog.records if r.levelno >= logging.ERROR]
    assert len(export_warnings) == 1
    assert errors == []
    assert broken_exporter["n"] == 1  # later figures skip the exporter entirely
    assert len(list((tmp_path / "plots").glob("*.html"))) == 10
    assert len(list((tmp_path / "plots").glob("*.svg"))) == 10  # matplotlib fallback


def test_non_dependency_export_errors_still_raise(tmp_path: Path, monkeypatch):
    import plotly.graph_objects as go

    from gromacs_analysis.utils import figure_export

    figure_export.reset_static_export_state()
    monkeypatch.setattr(go.Figure, "write_image", lambda self, *a, **k: (_ for _ in ()).throw(ValueError("bad layout")))
    with pytest.raises(ValueError):
        figure_export.write_static(go.Figure(), tmp_path / "x.svg")
    assert figure_export._EXPORT_STATE["available"] is None
    figure_export.reset_static_export_state()


# ---------------------------------------------------------------- US11 run report

def _pipeline_config(tmp_path: Path, stages):
    from gromacs_analysis.pipeline.config import PipelineConfig, PipelineStage

    cfg = PipelineConfig(project_name="Demo", protein_name="ProteinX", output_root=tmp_path / "out",
                         run=[name for name, _ in stages])
    for name, stage_cfg in stages:
        cfg.stages[name] = PipelineStage(name=name, enabled=True, config=stage_cfg)
    return cfg


class _Cfg:
    def __init__(self, output_dir):
        self.output_dir = output_dir


def test_run_report_lists_stages_figures_tables_and_warnings(tmp_path: Path):
    from gromacs_analysis.utils.report import write_run_report

    md_out = tmp_path / "out" / "md"
    (md_out / "plots").mkdir(parents=True)
    (md_out / "plots" / "rmsd.svg").write_text("<svg xmlns='http://www.w3.org/2000/svg'/>")
    (md_out / "plots" / "rmsd.html").write_text("<html></html>")
    (md_out / "plots" / "overview.html").write_text("<html></html>")
    (md_out / "statistics_summary.csv").write_text("System,Metric,Mean\nLigA,rmsd_prot,1.5\n")
    qc_out = tmp_path / "out" / "qc"
    qc_out.mkdir(parents=True)
    cfg = _pipeline_config(tmp_path, [("qc", _Cfg(qc_out)), ("md", _Cfg(md_out)), ("ranking", _Cfg(tmp_path / "out" / "r"))])
    results = {
        "qc": {"success": True, "warnings": ["traj.xtc: ligand loses contact"], "errors": []},
        "md": {"success": True},
        "ranking": {"success": False, "error": "No metrics found for ranking"},
    }

    report = write_run_report(cfg, results)
    text = report.read_text()

    assert report == tmp_path / "out" / "report.html"
    assert text.index("ligand loses contact") < text.index("<section class='card' id='qc'>")
    assert "src='md/plots/rmsd.svg'" in text and "href='md/plots/rmsd.html'" in text
    assert "href='md/plots/overview.html'" in text
    assert "statistics_summary.csv" in text and "rmsd_prot" in text
    assert "No metrics found for ranking" in text and "failed" in text
    assert "No figures were produced" in text  # ranking stage
    assert str(tmp_path) not in text  # links are relative


def test_run_pipeline_writes_report(tmp_path: Path):
    from gromacs_analysis.pipeline import run_pipeline
    from gromacs_analysis.qc import QCConfig

    traj = tmp_path / "rmsd.xvg"
    traj.write_text("0 0.1\n10 0.1\n")
    cfg = _pipeline_config(tmp_path, [("qc", QCConfig(base_dir=tmp_path, output_dir=tmp_path / "out" / "qc",
                                                       trajectories=[traj]))])
    results = run_pipeline(cfg)
    assert Path(results["_report"]).exists()


# ---------------------------------------------------------------- US12 cleanup

def test_split_loader_keeps_public_and_private_imports(tmp_path: Path):
    from gromacs_analysis.pipeline import stage_builders, templates
    from gromacs_analysis.pipeline.yaml_loader import (  # noqa: F401
        _batch_systems,
        _build_pca_joint_config,
        _build_prolif_batch_configs,
        _build_prolif_compare_config,
        _build_prolif_config,
        _build_ttclust_batch_configs,
        _build_ttclust_config,
        _default_output_dir,
        _resolve_path,
        _resolve_stage_input_path,
        generate_pipeline_template,
        load_pipeline_config,
    )

    assert _build_prolif_config is stage_builders._build_prolif_config
    out = tmp_path / "template.yaml"
    generate_pipeline_template(out)
    assert out.read_text() == templates.PIPELINE_TEMPLATE
    config = load_pipeline_config(out)
    assert "pca_joint" in config.stages and "prolif_batch" in config.stages


def test_repo_template_loads_with_new_stages():
    from gromacs_analysis.pipeline.yaml_loader import load_pipeline_config

    config = load_pipeline_config(ROOT / "gromacs_pipeline_template.yaml")
    assert "pca_joint" in config.stages
    assert config.stages["ranking"].config.exclude_systems == ["APO"]


@pytest.mark.parametrize("script", ["md_analysis_base.py", "holo_md_analysis.py", "apo_md_analysis.py"])
def test_legacy_scripts_warn_deprecation(script):
    import subprocess
    import sys

    proc = subprocess.run([sys.executable, "-W", "error::DeprecationWarning", "-c",
                           f"import runpy; runpy.run_path({str(ROOT / script)!r}, run_name='legacy')"],
                          capture_output=True, text=True, cwd=ROOT)
    assert proc.returncode != 0 and "DeprecationWarning" in proc.stderr
    assert "gromacs-md" in proc.stderr or "gromacs_analysis.md" in proc.stderr
