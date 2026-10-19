# Tasks: Pipeline Hardening Round 3

**Input**: Design documents from `specs/001-pipeline-hardening-round3/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, contracts/

**Tests**: Required. Constitution Principle IV makes tests mandatory for every scientific
behaviour, so each story starts with a failing test that has a known answer.

**Organization**: Tasks are grouped by user story (US1–US12 from spec.md).

## Format: `[ID] [P?] [Story] Description`

## Phase 1: Setup

- [X] T001 Add `pytest` and `MDAnalysis` to `requirements-dev.txt` and confirm `pip install -e .` works for the CI install path in `requirements-dev.txt`
- [X] T002 [P] Create empty test modules `tests/test_round3_scientific.py`, `tests/test_round3_engineering.py`, `tests/test_end_to_end_trajectory.py` with module docstrings

---

## Phase 2: Foundational (blocking)

- [X] T003 Implement `StaticExportUnavailable` and a cached availability check in `gromacs_analysis/utils/figure_export.py`. On the first failing `write_image`, decide whether the error means a missing browser/Kaleido. If so, log one warning naming the remedies (`plotly_get_chrome`, `BROWSER_PATH`, `pip install "kaleido<1"`) and raise `StaticExportUnavailable` without retrying on later calls. Every plotter depends on this (US9 and all figure-producing stories).
- [X] T004 Create `gromacs_analysis/utils/units.py` with `parse_xvg_units(content) -> (x_unit, y_unit)`, `display_unit_from_label(label)` and `conversion_factor(declared, display) -> (factor, status)`. Statuses: `converted`, `matched`, `assumed`, `unrecognised`. Recognise `nm`, `nm^2`, `nm\S2\N`, `nm²`, `Å`, `A`, `Angstrom`, `Å²`, `A^2` (US3, used by US6).

---

## Phase 3: User Story 1 - Joint PCA on shared axes (P1) 🎯 MVP

**Goal**: One PCA basis for all systems, every system projected on it, RMSIP overlap.
**Independent Test**: Synthetic systems with a known displacement direction recover it as PC1 (|cos| ≥ 0.95); self RMSIP = 1.

- [X] T005 [P] [US1] Write tests in `tests/test_round3_scientific.py` for:
  - Kabsch fit removing a pure rotation;
  - PC1 recovery on synthetic data;
  - shared projections having identical eigenvectors;
  - RMSIP equal to 1 on the diagonal, symmetric, and in [0, 1];
  - an atom-count mismatch raising an error that names the systems.
- [X] T006 [US1] Implement `gromacs_analysis/pca/joint.py`:
  - `kabsch_fit(coords, ref)`;
  - `iterative_fit(systems, iterations=2)`;
  - `joint_pca(coords_by_system, n_components)`, returning a `JointPCAResult`. Eigenvalues are sorted descending and eigenvector signs normalised so the largest-|component| is positive.
  - `rmsip(a, b, k)`;
  - `load_coordinates(topology, trajectories, selection, stride, max_frames)` using MDAnalysis, imported lazily.
  - Validation: "S ≥ 1; every system has ≥ 2 frames; n_overlap ≤ min(3N, frames−1)".
- [X] T007 [US1] Add `PCAJointConfig` (fields per contracts/config-schema.md) and `PCAJointSystem` to `gromacs_analysis/pca/config.py`, and export them from `gromacs_analysis/pca/__init__.py`.
- [X] T008 [US1] Implement `PCAJointAnalyzer` in `gromacs_analysis/pca/joint_analyzer.py`. It writes the CSV/JSON listed in contracts/outputs.md and draws `joint_projection_12`, `joint_density_12`, `rmsip_heatmap` and `joint_scree` with the shared palette (APO grey) via `write_static`.
- [X] T009 [US1] Add the `pca_joint` stage. The builder goes in `gromacs_analysis/pipeline/stage_builders.py`, supporting `systems` or `auto_discover`. Register it in `gromacs_analysis/pipeline/yaml_loader.py` (stage list and dispatch) and in `gromacs_analysis/pipeline/runner.py`.
- [X] T010 [US1] Add a test in `tests/test_round3_scientific.py` that runs `PCAJointAnalyzer` on synthetic coordinates with HTML-only output and asserts the files listed in contracts/outputs.md exist.

---

## Phase 4: User Story 2 - CI on every PR (P1)

**Independent Test**: The workflow file is valid YAML with a `pull_request` trigger and a Python 3.10/3.12 matrix running pytest.

- [X] T011 [US2] Create `.github/workflows/tests.yml`: triggers on push to main and on pull_request; matrix Python 3.10 and 3.12; installs `requirements.txt` plus `requirements-dev.txt`; runs `pytest -q`.
- [X] T012 [P] [US2] Add a test in `tests/test_round3_engineering.py` that parses the workflow and checks the triggers, matrix and pytest step.

---

## Phase 5: User Story 3 - Units (P2)

**Independent Test**: The same RMSD given once as an nm `.xvg` and once in Å is reported within 1%.

- [X] T013 [P] [US3] Write tests in `tests/test_round3_scientific.py`:
  - nm→Å factor 10;
  - nm²→Å² factor 100;
  - GROMACS `nm\S2\N` parsing;
  - a file without a header leaves values unchanged with status `assumed`;
  - an unknown unit gives `unrecognised` and no change;
  - an MD processor end-to-end case for two systems.
- [X] T014 [US3] In `gromacs_analysis/md/data_processor.py`, record the header y-unit per file, apply the conversion to values (not to the RMSF residue axis), and store a `unit_metadata` dict on each metric result. Add `convert_units: bool = True` to `MDConfig` in `gromacs_analysis/md/config.py` and an optional `AnalysisMetric.unit`.
- [X] T015 [US3] Write the `Value_Unit_Declared`, `Value_Unit`, `Value_Unit_Factor` and `Value_Unit_Status` columns in `gromacs_analysis/md/analyzer.py::_save_statistics_csv`, and pass `convert_units` through `gromacs_analysis/pipeline/yaml_loader.py`.

---

## Phase 6: User Story 4 - Equilibration (P2)

**Independent Test**: On relaxation-plus-noise series, the detected start lies within ±15% of the true end in at least 9 of 10 seeds.

- [X] T016 [P] [US4] Write tests in `tests/test_round3_scientific.py`:
  - statistical inefficiency ≈ 1 for white noise and ≈ (1+φ)/(1−φ) for AR(1) within 25%;
  - detection accuracy over 10 seeds;
  - a short series reports `insufficient_data`;
  - `equilibration.json` is written;
  - MM/PBSA `start_frame` filtering;
  - the VMD `rmsf_first_frame`.
- [X] T017 [US4] Implement `gromacs_analysis/utils/equilibration.py` with `statistical_inefficiency(x, mintime=3)` and `detect_equilibration(x, n_candidates=100, max_fraction=0.5)`, returning `(t0, g, n_eff)`. Use pymbar when `use_pymbar=True` and it is importable.
- [X] T018 [US4] Add `method` ("chodera" | "window", default "chodera") and `n_candidates` to `StabilityConfig` in `gromacs_analysis/md/config.py`. In `gromacs_analysis/md/data_processor.py::_apply_stability_window`, dispatch on the method and store records. Write `data/equilibration.json` from `gromacs_analysis/md/analyzer.py`. Wire the new keys in `gromacs_analysis/pipeline/yaml_loader.py` and `gromacs_analysis/config/yaml_loader.py`.
- [X] T019 [US4] Add `start_frame: Optional[int] = None` to `MMPBSAConfig` in `gromacs_analysis/mmpbsa/config.py`. In `gromacs_analysis/mmpbsa/processor.py`, recompute delta components and binding energy from frame-level data with `Frame #` ≥ start_frame and record `frames_used` in the QC metadata. Wire it in `gromacs_analysis/pipeline/yaml_loader.py`.
- [X] T020 [US4] Add `rmsf_first_frame` to the VMD helper in `gromacs_analysis/pipeline/upstream.py` (`measure rmsf $sel first N`).

---

## Phase 7: User Story 5 - Ranking with uncertainty (P2)

**Independent Test**: A 10-SEM separation gives P ≥ 0.99; identical ligands give P in [0.4, 0.6].

- [X] T021 [P] [US5] Write tests in `tests/test_round3_scientific.py`:
  - the separation probabilities;
  - determinism with a fixed seed;
  - APO exclusion (exact name and configured list) while "Apomorphine" is kept;
  - all systems excluded raises an error;
  - `Replicate_Means` is used when present.
- [X] T022 [US5] Write a `Replicate_Means` column (semicolon-separated) in `gromacs_analysis/md/analyzer.py::_save_statistics_csv`.
- [X] T023 [US5] Add `exclude_systems`, `n_bootstrap` (default 2000) and `seed` (default 0) to `RankingConfig` in `gromacs_analysis/ranking/config.py`. Implement exclusion and the bootstrap in `gromacs_analysis/ranking/analyzer.py`, writing `ligand_ranking.csv` (new columns), `rank_pairwise.csv` and `excluded_systems.csv`.
- [X] T024 [US5] In `gromacs_analysis/pipeline/yaml_loader.py::_build_ranking_config`, merge the APO system names from the md/md_compare stage configs into `exclude_systems`.

---

## Phase 8: User Story 6 - Protein–ligand H-bonds (P3)

- [X] T025 [P] [US6] Write tests in `tests/test_round3_scientific.py`: the VMD script contains `hbonds_lig_.dat` and a per-ligand `hbonds_lig_<name>.dat`; the `hbonds_lig` default metric exists and is split per ligand.
- [X] T026 [US6] Add a `calc_hbonds_ligand` proc and its calls in `gromacs_analysis/pipeline/upstream.py`. Add the `hbonds_lig` default metric and a `LIGAND_METRICS` entry in `gromacs_analysis/md/config.py`, and give `hbonds_lig` the H-bond axis limits in `gromacs_analysis/md/plotter.py`.

---

## Phase 9: User Story 7 - Flat-profile terminal detection (P3)

- [X] T027 [P] [US7] Write tests in `tests/test_round3_scientific.py`: an exactly flat core with 6/10 elevated termini is trimmed at least 6/10; a non-degenerate profile gives an unchanged threshold.
- [X] T028 [US7] Update `robust_threshold` in `gromacs_analysis/pca/prep.py`: when MAD = 0, scale = 0.7413·IQR, and if that is also 0, max(0.05·|median|, 1e-6).

---

## Phase 10: User Story 8 - Multi-chain distances (P3)

- [X] T029 [P] [US8] Write tests in `tests/test_round3_scientific.py` on a two-chain MDAnalysis universe built in memory:
  - "A:45" selects only chain A;
  - a bare 45 returns two entries with a warning;
  - labels are chain-prefixed;
  - an unknown chain raises an error.
- [X] T030 [US8] Re-key residues by `(chain, resid)` in `gromacs_analysis/distance/processor.py` (resolution, distance loop, labels) and adapt consumers in `gromacs_analysis/distance/analyzer.py`, `plotter.py` and `compare.py` that look residues up by `resid`.

---

## Phase 11: User Story 9 - Export availability (P2)

- [ ] T031 [P] [US9] Write tests in `tests/test_round3_engineering.py`: when the simulated exporter raises a Chrome error, ten saves log exactly one warning and later calls raise `StaticExportUnavailable` immediately; the MD plotter still writes HTML and the matplotlib SVG fallback.
- [ ] T032 [US9] Make the plotter save loops in `gromacs_analysis/md/plotter.py`, `gromacs_analysis/pca/plotter.py`, `gromacs_analysis/pca/compare.py`, `gromacs_analysis/mmpbsa/plotter.py`, `gromacs_analysis/distance/plotter.py`, `gromacs_analysis/prolif/plotter.py` and `gromacs_analysis/ttclust/analyzer.py` treat `StaticExportUnavailable` as a quiet fallback (debug log) instead of an error.

---

## Phase 12: User Story 10 - End-to-end trajectory test (P2)

- [ ] T033 [US10] In `tests/test_end_to_end_trajectory.py`, build a small protein (two chains of glycine-like residues with CA) plus a ligand residue and write PDB+DCD (or XTC). Run `QCAnalyzer`, then `DistanceAnalyzer` on a drifting ligand (distance increasing), then `PCAJointAnalyzer` on two systems. Use `importorskip("MDAnalysis")`; run ProLIF only if importable.

---

## Phase 13: User Story 11 - Run report (P3)

- [ ] T034 [P] [US11] Write tests in `tests/test_round3_engineering.py`: the report lists stages with status, QC warnings first, figure links relative, and a failed stage's error shown; it works when there are no figures.
- [ ] T035 [US11] Implement `write_run_report(config, results, output_root)` in `gromacs_analysis/utils/report.py` and call it at the end of `run_pipeline` in `gromacs_analysis/pipeline/runner.py` (and on `continue_on_error`). Return the path in the results under `_report`.

---

## Phase 14: User Story 12 - Cleanup (P3)

- [ ] T036 [P] [US12] Write tests in `tests/test_round3_engineering.py`:
  - importing the moved names from `gromacs_analysis.pipeline.yaml_loader` still works;
  - the template YAML loads identically;
  - each legacy script module emits a `DeprecationWarning` when run as `__main__` (check the source contains the warning call).
- [ ] T037 [US12] Move the template text to `gromacs_analysis/pipeline/templates.py`. Move `_batch_systems`, `_build_prolif_config`, `_build_prolif_batch_configs`, `_build_prolif_compare_config`, `_build_ttclust_config`, `_build_ttclust_batch_configs` and the `pca_joint` builder to `gromacs_analysis/pipeline/stage_builders.py`. Re-export them from `gromacs_analysis/pipeline/yaml_loader.py`.
- [ ] T038 [P] [US12] Add a `DeprecationWarning` and a stderr notice naming `gromacs-md` / `gromacs-pipeline` to `md_analysis_base.py`, `holo_md_analysis.py` and `apo_md_analysis.py`.

---

## Phase 15: Polish

- [ ] T039 Update `SCIENTIFIC_VALIDITY.md` (move the resolved caveats into the fixed table), `USAGE.md`, `README.md` and both YAML templates (`gromacs_pipeline_template.yaml` and the in-code template) for the new keys and stage.
- [ ] T040 Run `pytest -q` and render the joint PCA figures from synthetic data for a visual check.
- [ ] T041 Mark the tasks complete in this file, then commit and push to `ccr-deb837be-iavlqx`.

---

## Dependencies & Execution Order

- Phase 1 → Phase 2 → all stories.
- T003 (export) precedes T008, T032. T004 (units) precedes T014.
- US1 (T005–T010) precedes US10 (T033, which uses the joint PCA) and US12 (T037, which moves the `pca_joint` builder).
- US5 T022 depends on the MD statistics writer that US3 T015 also edits, so do T015 before T022 (same file).
- US11 T035 and US1 T009 both edit `runner.py`; do them sequentially.
- All other stories are independent.

## Parallel Example

```text
T005 [US1] tests  |  T013 [US3] tests  |  T016 [US4] tests  |  T021 [US5] tests  |  T027 [US7]  |  T029 [US8]
```

## Implementation Strategy

- MVP: Phases 1–3 (joint PCA) plus Phase 4 (CI), so that every later step is tested automatically.
- Then the P2 stories in this order: US3, US4, US5, US9, US10.
- Then the P3 stories: US6, US7, US8, US11, US12.
- Polish last.
