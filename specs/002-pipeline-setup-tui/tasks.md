# Tasks: Pipeline Setup TUI

**Input**: `specs/002-pipeline-setup-tui/` (plan, spec, research, data-model, contracts)
**Tests**: Required (Constitution IV).

## Phase 1: Setup
- [ ] T001 Add the `tui` extra (`textual>=0.60`, `paramiko>=3.4`) and the console script `clarity-setup=gromacs_analysis.setup.cli:main` to `setup.py`.
- [ ] T002 [P] Create the `gromacs_analysis/setup/__init__.py` package and the test files `tests/test_setup_scanner.py`, `tests/test_setup_project.py` and `tests/test_setup_tui.py`.

## Phase 2: Foundational
- [ ] T003 Implement `gromacs_analysis/setup/backends.py`:
  - a `Backend` base with `run`, `read_text`, `write_text`, `exists`, `is_remote` and `label`;
  - `LocalBackend`;
  - `ShellRunner(LocalBackend)`, which marks `is_remote=True` to exercise remote code paths in tests;
  - `SSHBackend`: paramiko, `~/.ssh/config` aware, key-only, TOFU host keys;
  - `CondorctlBackend`: wraps `condorctl.config.Config` and `RemoteConnection`;
  - `connect(kind, host, profile)`.

## Phase 3: US1 Discover (P1) 🎯 MVP
**Independent test**: a scan of a fixture tree yields LigA(3), LigB(2), APO(3) with ligands, locally and via `ShellRunner`.
- [ ] T004 [P] [US1] In `tests/test_setup_scanner.py`, write a fixture builder for CHARMM-GUI-like trees (psf/pdb/gro topologies, `step5_production.xtc` + `*_noPBC.xtc`, metric files, MMPBSA results) and tests for:
  - grouping (all naming variants and nested `LigA/rep1`);
  - APO detection;
  - multi-ligand detection, with shared ions filtered;
  - trajectory preference;
  - inputs detection;
  - skipped-directory count;
  - the remote command budget (count `run` calls on `ShellRunner`).
- [ ] T005 [US1] Implement `gromacs_analysis/setup/ligands.py`: `resname_command(path)` per format, `parse_resnames(output)`, `filter_ligands(resnames, apo_resnames)` reusing `pipeline.discovery` sets.
- [ ] T006 [US1] Implement `gromacs_analysis/setup/scanner.py`:
  - `scan(backend, root, depth=4) -> ScanResult` (systems, skipped, commands_used);
  - the `find` listing with a non-`-printf` fallback;
  - replicate-token stripping, file classification, trajectory choice, APO tokens and input detection.

## Phase 4: US2 Naming (P1)
**Independent test**: renames propagate to every generated stage; duplicate names are rejected; excluded systems are absent.
- [ ] T007 [P] [US2] In `tests/test_setup_project.py`, write tests for the rename/exclude/APO edits, the uniqueness validation, project save/load round-trip, and rescan merge (kept edits, new systems added, vanished systems marked missing).
- [ ] T008 [US2] Implement `gromacs_analysis/setup/project.py`: `Project`, `SystemEntry`, `Replicate`, `LigandEntry`, `StageChoice`, `from_scan`, `merge_scan`, `rename_system`, `rename_ligand`, `set_apo`, `set_included`, `validate()`, `save()` and `load()`.

## Phase 5: US3 Configuration (P1)
**Independent test**: the generated YAML loads with `load_pipeline_config` and a dry run; renames appear in md, prolif_batch, ttclust_batch and pca_joint; selections use residue names.
- [ ] T009 [P] [US3] In `tests/test_setup_project.py`, write generator tests:
  - stage availability and pre-selection reasons;
  - the YAML validates;
  - names and selections;
  - multi-ligand ProLIF runs;
  - APO excluded from ranking;
  - mmpbsa included only when results exist.
- [ ] T010 [US3] Implement `gromacs_analysis/setup/generator.py`: `stage_availability(project)`, `build_pipeline_dict(project)`, `write_pipeline(project, path)`, and validation through `load_pipeline_config` in a temp file, with a header comment.

## Phase 6: US4 Interface (P2)
**Independent test**: headless pilot: load project, rename via cell edit, toggle APO and a stage, save; files reflect the actions.
- [ ] T011 [P] [US4] In `tests/test_setup_tui.py`, write pilot tests (skip if textual missing) for: tabs present, system rename, APO toggle, exclude toggle, ligand rename, stage toggle, preview content, save writes files, an error notice on an invalid duplicate name, and a connect-screen scan of a local folder.
- [ ] T012 [US4] Implement `gromacs_analysis/setup/tui.py`:
  - `ClaritySetupApp(project | root/backend)` with a header and footer;
  - `TabbedContent`: Systems (`DataTable`), Ligands (`DataTable`), Stages (checkboxes with reasons), Preview (YAML);
  - an edit modal for names;
  - key bindings `e` edit, `a` APO, `x` include, `s` save, `r` run, `c` condor submit, `ctrl+r` rescan, `q` quit;
  - a connection screen when no project is given;
  - errors shown via notify, never crashing.

## Phase 7: US5 Run (P2)
**Independent test**:
- a local run executes the pipeline and returns the report path;
- a remote run writes `clarity/pipeline.yaml` through the backend and issues one `nohup` start command;
- condor submission builds `run_clarity.sh` and calls condorctl, mocked.
- [ ] T013 [P] [US5] Write runner tests in `tests/test_setup_project.py` using `ShellRunner` and a fake condorctl submit.
- [ ] T014 [US5] Implement `gromacs_analysis/setup/runner.py`: `run_local`, `run_on_host` (checks `command -v gromacs-pipeline` and reports a clear error) and `submit_condor` (condorctl API when importable; otherwise a clear error).

## Phase 8: CLI
- [ ] T015 [P] Write CLI tests in `tests/test_setup_project.py`:
  - `scan --json`;
  - `init` with renames, which writes valid files;
  - `tui` without textual exits 3 with a hint (simulated).
- [ ] T016 Implement `gromacs_analysis/setup/cli.py` per `contracts/cli.md`.

## Phase 9: Polish
- [ ] T017 Document the workflow in `README.md` and `USAGE.md` (the new "Setup without YAML" section), mention `clarity-setup tui` in `gromacs-pipeline --help`, and install textual/paramiko in CI so the TUI tests run.
- [ ] T018 Run the full suite and render a TUI screenshot (SVG via Textual) of the Systems tab for a visual check.
- [ ] T019 Mark tasks done, then commit and push.

## Dependencies
- T003 precedes everything that touches a backend.
- T005 and T006 (scanner) precede T008, which precedes T010, which precedes T012, T014 and T016.
- The TUI (T012) uses the generator (T010) and the runner (T014).

## MVP
Phases 1–5 make the CLI `init` usable without the TUI. Phase 6 adds the interface, and Phase 7 adds remote runs.
