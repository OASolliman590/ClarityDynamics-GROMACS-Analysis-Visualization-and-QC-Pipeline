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
