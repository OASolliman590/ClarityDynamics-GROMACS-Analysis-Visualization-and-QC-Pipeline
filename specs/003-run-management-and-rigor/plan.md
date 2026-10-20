# Implementation Plan: Run Management, Input Preparation and Statistical Rigour

**Branch**: `ccr-deb837be-iavlqx` | **Date**: 2026-10-20 | **Spec**: [spec.md](./spec.md)

## Summary

The work extends `gromacs_analysis.setup` and three analysis modules (no rewrites):

- **Runs**: `setup/jobs.py` holds run records and status/fetch for host, condor and DAG runs.
  The project file gains `runs`.
- **Inputs**: `setup/prepare.py` writes per-replicate scripts for trjconv, the VMD metrics
  helper (reusing `pipeline.upstream`) and gmx_MMPBSA, and runs them. The scanner records
  `tpr`/`top`/raw trajectory per replicate.
- **Parallel**: `setup/workflow.py` splits the pipeline dict into unit configs plus a DAGMan
  file. The generator emits one ProLIF entry per ligand.
- **Settings**: `setup/settings.py` holds the per-stage schema, validation and application in
  the generator. The interface gets a settings modal; the CLI gets `--set`.
- **Caching**:
  - `pipeline/cache.py` gains `config_fingerprint`.
  - `pipeline/runner.py` caches per unit inside batch stages.
  - The generator writes `cache: true`.
  - `gromacs-pipeline` gains `--no-cache`.
- **Doctor**: `setup/doctor.py` holds checks as data with a remedy each.
- **Convergence**: `utils/convergence.py` provides block SEM, statistical inefficiency, half
  agreement, cosine content and half-split RMSIP. These are hooked into the MD analyzer and
  joint PCA, and the report shows the flags.
- **Statistics**: `ranking/stats.py` provides Welch, exact or Monte-Carlo permutation, Hedges'
  g, Holm and BH. It is hooked into ranking, with MM/PBSA replicate means read from
  `*_frame_results.csv`.
- **Packaging/CI/provenance**:
  - `pyproject.toml` with `setup.py` reduced to a shim;
  - `environment.yml`;
  - ruff and mypy configuration;
  - CI coverage;
  - `utils/provenance.py` and a report section.

## Technical Context

**Language/Version**: Python 3.10+ for setup/TUI (3.9 core unchanged)
**Primary Dependencies**: numpy, scipy (t-distribution; already a dependency), pandas, PyYAML; optional textual, paramiko, condorctl
**Testing**: pytest. Remote paths use `ShellRunner` plus fake executables (`condor_q`, `condor_submit_dag`, `gmx`) placed on `PATH` in a temp dir, so the generated commands really execute.
**Constraints**: No new required dependencies. Convergence and statistics are vectorised numpy. The DAG relies on a shared filesystem.

## Constitution Check

| Principle | Compliance |
|-----------|-----------|
| I. Scientific correctness | Standard estimators are used: Flyvbjerg-Petersen blocking, Chodera statistical inefficiency, Hess cosine content, Welch, Holm and BH. Each is tested against analytic cases. |
| II. Replicates | Convergence is computed per replicate. Tests use replicate means as the unit. |
| III. Flag, don't discard | Convergence flags, prep skip reasons and doctor failures are all reported. Nothing is dropped. |
| IV. Tests | Every new module has unit tests, and the remote flows run against fake executables. |
| V. Provenance | Run records, the report provenance section, and generated scripts carry headers. |
| VI. Visualisation | Report additions are text and tables; there are no colour-only signals. |
| VII. Explicit changes | Caching is on in generated configs and the change is documented. Hand-written YAML keeps opt-in caching. The ProLIF entry layout produces the same output folders as before. |

**Gate**: PASS.

## Project Structure

```text
gromacs_analysis/setup/{jobs,prepare,workflow,settings,doctor}.py   (new)
gromacs_analysis/setup/{scanner,project,generator,runner,cli,tui}.py (extended)
gromacs_analysis/utils/{convergence,provenance}.py                 (new)
gromacs_analysis/ranking/stats.py                                   (new)
gromacs_analysis/pipeline/{cache,runner}.py, cli/pipeline_cli.py    (extended)
gromacs_analysis/md/analyzer.py, pca/{joint,joint_analyzer}.py, utils/report.py (extended)
pyproject.toml, environment.yml, setup.py (shim), .github/workflows/tests.yml
tests/test_setup_runs.py, tests/test_convergence_stats.py, tests/test_cache_provenance.py
```
