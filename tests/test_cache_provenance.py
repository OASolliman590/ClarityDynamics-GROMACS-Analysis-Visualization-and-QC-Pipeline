"""Per-unit caching with config fingerprints, and report provenance (specs/003, US5/US10)."""

import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

import pandas as pd
import pytest

from gromacs_analysis.pipeline import runner
from gromacs_analysis.pipeline.cache import config_fingerprint
from gromacs_analysis.pipeline.config import PipelineConfig, PipelineStage


@dataclass
class FakeUnit:
    output_dir: Path
    topology: Path
    trajectories: List[Path]
    stride: int = 10
    extras: dict = field(default_factory=dict)


class CountingAnalyzer:
    calls: List[str] = []

    def __init__(self, cfg):
        self.cfg = cfg

    def run_analysis(self):
        CountingAnalyzer.calls.append(self.cfg.output_dir.name)
        return {"success": True}


class FailingAnalyzer(CountingAnalyzer):
    def run_analysis(self):
        CountingAnalyzer.calls.append(self.cfg.output_dir.name)
        return {"success": False, "error": "boom"}


@pytest.fixture
def batch(tmp_path: Path, monkeypatch):
    CountingAnalyzer.calls = []
    monkeypatch.setitem(runner.BATCH_UNITS, "distance_batch", ("distance", CountingAnalyzer))
    monkeypatch.setattr(runner, "write_run_report", lambda *a, **k: tmp_path / "report.html")
    units = []
    for name in ("A", "B"):
        top, traj = tmp_path / f"{name}.pdb", tmp_path / f"{name}.xtc"
        top.write_text("x")
        traj.write_text("y")
        units.append(FakeUnit(tmp_path / "out" / name, top, [traj]))
    stage = PipelineStage(name="distance_batch", config=units, cache_enabled=True)
    return PipelineConfig(run=["distance_batch"], stages={"distance_batch": stage}), units


def test_fingerprint_is_stable_and_sensitive(tmp_path: Path):
    a = FakeUnit(tmp_path, tmp_path / "t", [tmp_path / "x"], extras={"b": 1, "a": {2, 1}})
    b = FakeUnit(tmp_path, tmp_path / "t", [tmp_path / "x"], extras={"a": {1, 2}, "b": 1})
    assert config_fingerprint(a) == config_fingerprint(b)
    b.stride = 5
    assert config_fingerprint(a) != config_fingerprint(b)


def test_batch_units_cache_individually(batch):
    config, units = batch
    runner.run_pipeline(config)
    assert sorted(CountingAnalyzer.calls) == ["A", "B"]
    result = runner.run_pipeline(config)
    assert sorted(CountingAnalyzer.calls) == ["A", "B"]  # both cached
    assert all(r.get("skipped") for r in result["distance_batch"])

    time.sleep(0.01)
    traj = units[1].trajectories[0]
    traj.write_text("changed input")
    os.utime(traj, (time.time() + 5, time.time() + 5))
    runner.run_pipeline(config)
    assert CountingAnalyzer.calls[-1] == "B" and len(CountingAnalyzer.calls) == 3

    units[0].stride = 2  # a setting change invalidates that unit
    runner.run_pipeline(config)
    assert CountingAnalyzer.calls[-1] == "A" and len(CountingAnalyzer.calls) == 4


def test_no_cache_and_failures_are_not_cached(batch, monkeypatch):
    config, _ = batch
    runner.run_pipeline(config)
    runner.run_pipeline(config, use_cache=False)
    assert len(CountingAnalyzer.calls) == 4
    monkeypatch.setitem(runner.BATCH_UNITS, "distance_batch", ("distance", FailingAnalyzer))
    for unit in config.stages["distance_batch"].config:
        (unit.output_dir / ".cache.json").unlink()
    runner.run_pipeline(config, continue_on_error=True)
    runner.run_pipeline(config, continue_on_error=True)
    assert len(CountingAnalyzer.calls) == 8  # failed units re-run every time


def test_stage_cache_includes_settings(tmp_path: Path, monkeypatch):
    calls = []

    class FakeQC:
        def __init__(self, cfg):
            pass

        def run_analysis(self):
            calls.append(1)
            return {"success": True}

    @dataclass
    class QCCfg:
        output_dir: Path
        topology: Path
        trajectories: list
        threshold: float = 5.0
        index_file: object = None

    monkeypatch.setattr(runner, "QCAnalyzer", FakeQC)
    monkeypatch.setattr(runner, "write_run_report", lambda *a, **k: None)
    top = tmp_path / "t.pdb"
    top.write_text("x")
    cfg = QCCfg(tmp_path / "qc", top, [])
    config = PipelineConfig(run=["qc"], stages={"qc": PipelineStage(name="qc", config=cfg, cache_enabled=True)})
    runner.run_pipeline(config)
    runner.run_pipeline(config)
    assert len(calls) == 1
    cfg.threshold = 3.0
    runner.run_pipeline(config)
    assert len(calls) == 2


def test_cli_no_cache_flag(monkeypatch, tmp_path: Path):
    from gromacs_analysis.cli import pipeline_cli

    seen = {}
    yaml_path = tmp_path / "p.yaml"
    yaml_path.write_text("run: []\n")
    monkeypatch.setattr(pipeline_cli, "run_pipeline", lambda cfg, **kw: seen.update(kw))
    monkeypatch.setattr("sys.argv", ["gromacs-pipeline", "--config", str(yaml_path), "--no-cache"])
    pipeline_cli.main()
    assert seen["use_cache"] is False


def test_report_shows_sampling_flags_and_provenance(tmp_path: Path):
    from gromacs_analysis.utils.provenance import collect
    from gromacs_analysis.utils.report import write_run_report

    out = tmp_path / "out" / "md"
    (out / "data").mkdir(parents=True)
    pd.DataFrame([
        {"System": "LigA", "Metric": "rmsd_lig", "Replicate": 2, "N_eff": 7.0, "Half_Z": 5.2, "Flags": "low_neff;drift"},
        {"System": "LigA", "Metric": "rmsd_lig", "Replicate": 1, "N_eff": 90.0, "Half_Z": 0.4, "Flags": ""},
    ]).to_csv(out / "data" / "convergence.csv", index=False)
    cfg_file = tmp_path / "pipeline.yaml"
    cfg_file.write_text("project:\n  name: Demo\nrun: [md]\n")

    @dataclass
    class MDCfg:
        output_dir: Path

    config = PipelineConfig(project_name="Demo", output_root=tmp_path / "out", run=["md"],
                            stages={"md": PipelineStage(name="md", config=MDCfg(out))}, source_path=cfg_file)
    html = write_run_report(config, {"md": {"success": True}}).read_text()
    assert "Sampling checks (1 flagged)" in html and "low_neff;drift" in html and "still drifting" in html
    assert "Provenance" in html and "Configuration that ran" in html and "run: [md]" in html

    info = collect(cfg_file)
    assert info["python"] and info["packages"].get("numpy") and info["config_text"].startswith("project")


def test_loaded_config_records_source_path(tmp_path: Path):
    from gromacs_analysis.pipeline.yaml_loader import load_pipeline_config

    path = tmp_path / "p.yaml"
    path.write_text("run: []\n")
    assert load_pipeline_config(path).source_path == path.resolve()
