# Implementation Plan: Pipeline Setup TUI

**Branch**: `ccr-deb837be-iavlqx` | **Date**: 2026-10-19 | **Spec**: [spec.md](./spec.md)

## Summary

A new `gromacs_analysis.setup` package turns discovery, naming and stage choice into a
validated pipeline YAML. It has four layers:

1. **Filesystem backends**: `LocalBackend`, `SSHBackend` over paramiko, and `CondorctlBackend`
   reusing condorctl's `RemoteConnection`. All three expose the same small command interface.
   Discovery runs as one remote `find` plus one residue-name `awk` per unique topology.
2. **Scanner**: groups replicate directories into systems, picks topology and trajectories,
   detects APO systems, ligands and existing analysis inputs.
3. **Project model and generator**: a `Project` holds the scan plus the user's names, exclusions
   and stage choices. It is saved as `clarity_project.yaml`. The generator emits a pipeline YAML
   dict and validates it with the existing `load_pipeline_config`.
4. **Front ends**:
   - a non-interactive CLI, `clarity-setup scan|init|run`;
   - a Textual TUI, `clarity-setup tui`, with Connect → Systems → Ligands → Stages → Preview →
     Save/Run/Submit;
   - a runner for local runs, remote background runs and condorctl submission.

## Technical Context

**Language/Version**: Python ≥ 3.9 for core; TUI and condorctl paths need 3.10+ (their own requirement)
**Primary Dependencies**: PyYAML (existing); optional `textual` (TUI), `paramiko` (SSH), `condorctl`
**Storage**: `clarity_project.yaml` and generated `pipeline.yaml` (local; uploaded for remote runs)
**Testing**:
- pytest. The SSH backend is tested through a `ShellRunner` that runs the same command strings
  locally, so remote scan logic is exercised byte for byte.
- The TUI is tested headlessly with Textual's `App.run_test()` pilot.

**Target Platform**: Linux/macOS terminals, including over SSH
**Project Type**: library + CLI + TUI
**Performance Goals**: Scan of 50 replicate dirs ≤ 1 + #topologies remote commands; TUI responsive (<100 ms per action)
**Constraints**: Key-based SSH only; no passwords stored; core package import unaffected when optional deps absent
**Scale/Scope**: ~8 new modules, 1 console script, ~40 tests

## Constitution Check

| Principle | Compliance |
|-----------|-----------|
| I. Scientific correctness | Selections always use real residue names. The generator validates against the real loader. APO detection is explicit and editable. |
| II. Replicates | Replicates are grouped and passed as such to every stage. |
| III. Flag, don't discard | Skipped directories and unavailable stages are reported with reasons. Exclusions are explicit user choices, saved in the project file. |
| IV. Tested behaviour | Scanner, generator, backends (through the shell runner) and every TUI action have tests. |
| V. Provenance | The project file records the root, connection, scan time and choices. The generated YAML carries a header naming the project file. |
| VI. Visualisation | The TUI uses the shared palette for system colour chips and has no colour-only indicators: status always has text. |
| VII. Explicit changes | Additive: new commands and an optional extra. Existing YAML workflows are unchanged. |

**Gate**: PASS. Post-design: PASS.

## Project Structure

```text
gromacs_analysis/setup/
├── __init__.py
├── backends.py      # Backend protocol, LocalBackend, SSHBackend, CondorctlBackend, ShellRunner
├── scanner.py       # discovery: listing → systems, replicates, files, APO, inputs
├── ligands.py       # residue-name extraction commands + ligand filtering
├── project.py       # Project / SystemEntry / LigandEntry / StageChoice; save/load/merge
├── generator.py     # Project → pipeline YAML dict; validation via load_pipeline_config
├── runner.py        # run local / remote background / condorctl submit
├── cli.py           # clarity-setup scan|init|run|tui
└── tui.py           # Textual app (optional dependency)
tests/
├── test_setup_scanner.py
├── test_setup_project.py
└── test_setup_tui.py
```

**Structure Decision**: This is a separate subpackage, so the analysis packages are untouched. It
reuses `pipeline.discovery` filters (protein, water and ion residue names) and
`pipeline.yaml_loader.load_pipeline_config` for validation.
