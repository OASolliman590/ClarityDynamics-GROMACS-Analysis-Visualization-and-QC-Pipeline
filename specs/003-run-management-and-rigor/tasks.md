# Tasks: Run Management, Input Preparation and Statistical Rigour

**Input**: `specs/003-run-management-and-rigor/` (spec, plan, research, data-model, contracts)
**Tests**: Required (Constitution IV).

## Phase 1: Setup
- [ ] T001 Create test files `tests/test_convergence_stats.py`, `tests/test_cache_provenance.py` and `tests/test_setup_runs.py`, plus a fake-executable helper in `tests/setup_fixtures.py` (`fake_bin(tmp, {name: script})` prepended to `PATH`).

## Phase 2: Foundational (analysis-side, independent of setup)
- [X] T002 [P] [US7] Write tests for `utils/convergence.py`:
  - block SEM on white noise ≈ σ/√N;
  - a random walk is flagged by half-z and cosine content;
  - "too short" handling;
  - half-split RMSIP of identical halves = 1.
- [X] T003 [US7] Implement `gromacs_analysis/utils/convergence.py`:
  - `block_sem`, `n_effective`, `half_agreement`, `series_convergence(values, thresholds)`;
  - `cosine_content(proj, i)`, `half_split_rmsip(flat, k)`.
- [X] T004 [US7] Hook convergence into the MD analyzer:
  - `md/analyzer.py` writes `data/convergence.csv` from `replicates_raw`/`values`;
  - thresholds come from `MDConfig.convergence` (dict, defaults), parsed in the loader.
- [X] T005 [US7] Hook convergence into joint PCA:
  - `pca/joint.py` adds `half_rmsip` to the result;
  - `pca/joint_analyzer.py` records replicate segments and writes `data/convergence.csv`;
  - the config gains `convergence` thresholds.
- [X] T006 [P] [US8] Write tests for `ranking/stats.py`:
  - Welch against scipy;
  - exact permutation p for 3 vs 3 equals 0.1 at full separation;
  - Hedges' g sign and size;
  - Holm/BH against hand-computed values;
  - insufficient replicates.
- [X] T007 [US8] Implement `gromacs_analysis/ranking/stats.py` (`pairwise_tests(groups, directions)`, `holm`, `bh`), hook it into `RankingAnalyzer` (writing `pairwise_tests.csv`), and read MM/PBSA replicate totals from `*_frame_results.csv`.
- [X] T008 [P] [US5] Write tests:
  - `config_fingerprint` is stable across reloads and changes when a setting changes;
  - a batch stage caches per unit, and touching one input re-runs one unit;
  - `--no-cache` disables caching.
- [X] T009 [US5] Implement caching:
  - `pipeline/cache.py` `config_fingerprint`;
  - `pipeline/runner.py` uses the fingerprint and caches per unit for batch stages, plus `run_pipeline(use_cache=...)`;
  - `_collect_stage_files` covers `pca_joint`;
  - `--no-cache` in `cli/pipeline_cli.py`.
- [X] T010 [P] [US10] Write tests: provenance has version, python and the config text; the report shows Provenance and Convergence sections.
- [X] T011 [US10] Implement `gromacs_analysis/utils/provenance.py` and extend `utils/report.py` with a convergence-flags card and a provenance card; the runner passes the config path.

## Phase 3: Setup-side
- [ ] T012 [US2] Extend the scanner: `Replicate.files` (tpr, top, raw/processed trajectory). Carry it through `project.Replicate` with defaults for old project files.
- [ ] T013 [P] [US4] Write tests for settings: schema listing, `--set` parsing/coercion/bounds, application in the YAML, and invalid values rejected.
- [ ] T014 [US4] Implement `gromacs_analysis/setup/settings.py`; the generator applies settings, emits `cache: true`, and emits one ProLIF entry per ligand.
- [ ] T015 [P] [US3] Write workflow tests:
  - unit count and names;
  - every unit validates;
  - unit output dirs equal the full run's;
  - the DAG file has all units as parents of the final job;
  - fake `condor_submit_dag` records the cluster.
- [ ] T016 [US3] Implement `gromacs_analysis/setup/workflow.py` (`build_units`, `write_dag`, `submit_dag`).
- [ ] T017 [P] [US1] Write tests for jobs:
  - host status alive/dead with the log tail;
  - condor status from fake `condor_q`/`condor_history`;
  - fetch through `ShellRunner` unpacks the report;
  - records persist in the project file.
- [ ] T018 [US1] Implement `gromacs_analysis/setup/jobs.py` (`RunRecord`, `status`, `fetch`); runner functions return records, and `Project.runs` is added.
- [ ] T019 [P] [US2] Write tests for prepare: tasks per missing input, skip reasons, script content (named groups, `gmx select` index, mmpbsa.in), the dry run, and submission through a fake host.
- [ ] T020 [US2] Implement `gromacs_analysis/setup/prepare.py` (`plan_preparation`, `write_scripts`, `run_preparation`).
- [ ] T021 [P] [US6] Write tests for doctor: all checks against a fake host; a missing vmd is reported with a remedy; the exit code depends on what is required.
- [ ] T022 [US6] Implement `gromacs_analysis/setup/doctor.py`.
- [ ] T023 [P] Write CLI tests for `settings`, `--set`, `prepare --dry-run`, `run --where dag`, `status`, `fetch` and `doctor`.
- [ ] T024 Extend `gromacs_analysis/setup/cli.py` per `contracts/cli.md`.
- [ ] T025 [P] [US1] Write TUI pilot tests: Jobs tab lists runs; settings modal edits a value; prepare dry-run notice; DAG key without a remote connection gives an error.
- [ ] T026 [US1] Extend `gromacs_analysis/setup/tui.py`: Jobs tab, settings modal, keys `p`/`g`/`d`/`j`/`f`/`D`.

## Phase 4: Packaging & quality
- [ ] T027 [US9] Add `pyproject.toml` (PEP 621, extras, scripts, ruff/mypy/pytest/coverage config), reduce `setup.py` to a shim, and add `environment.yml`.
- [ ] T028 [US10] Update CI: `pip install -e .[tui]` plus test deps, `ruff check`, `mypy gromacs_analysis/setup`, and `pytest --cov` with a threshold.

## Phase 5: Polish
- [ ] T029 Docs: a README section; USAGE sections for run management, preparation, parallel runs, settings, caching, doctor, convergence and statistics; the cache behaviour change.
- [ ] T030 Run the full suite and render a screenshot of the Jobs/settings interface.
- [ ] T031 Mark tasks done, commit and push.

## Dependencies
- T003 precedes T004/T005.
- T009 precedes T016, because the final job relies on cache hits.
- T012 precedes T020.
- T014 precedes T016.
- T018 precedes T024/T026.
- Every implementation task follows its test task.
