# Implementation Plan: Pipeline Hardening Round 3

**Branch**: `ccr-deb837be-iavlqx` | **Date**: 2026-10-19 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `specs/001-pipeline-hardening-round3/spec.md`

## Summary

Twelve improvements, grouped into seven scientific items and five engineering items. The
scientific items are:

- joint PCA in a shared subspace, with RMSIP overlap;
- unit detection and conversion from `.xvg` headers;
- equilibration detection by maximum effective sample size, with the start propagated to RMSF
  and MM/PBSA;
- bootstrap rank probabilities, with reference systems excluded;
- protein–ligand H-bonds;
- a terminal-threshold fix for flat profiles;
- chain-aware distances.

The engineering items are CI, a single static-export availability check with fallback, a real
trajectory end-to-end test, a per-run HTML report, and cleanup.

The design adds two numeric core modules, `pca/joint.py` and `utils/equilibration.py`. They work
on plain numpy arrays, so the statistics can be tested without trajectory files. Thin MDAnalysis
loaders sit on top of them, and the rest is extensions of existing modules.

## Technical Context

**Language/Version**: Python ≥ 3.9 (CI on 3.10 and 3.12)
**Primary Dependencies**: numpy, scipy, pandas, plotly, matplotlib, PyYAML; optional MDAnalysis
(joint PCA loader, distances, QC, e2e test), ProLIF/RDKit (skipped when absent)
**Storage**: Files (CSV/JSON/HTML/SVG) under each stage's `output_dir`; run report in `output_root`
**Testing**: pytest; CI via GitHub Actions
**Target Platform**: Linux/macOS workstations and HPC login nodes
**Project Type**: Python library + CLI (`gromacs_analysis`)
**Performance Goals**: Joint PCA over 3–6 systems × ≤ 10k frames × ≤ 1,000 Cα in under a minute
on a laptop (covariance 3N × 3N ≤ 3000²); bootstrap 2,000 resamples in under 5 s for ≤ 20 ligands
**Constraints**: No new required dependency; disabled stages load without file-system access;
existing configs load unchanged
**Scale/Scope**: ~12 modules touched, 2 new core modules, 1 new pipeline stage (`pca_joint`)

## Constitution Check

| Principle | How this plan complies |
|-----------|------------------------|
| I. Scientific correctness | Units are read from headers and recorded (FR-007/008). Joint PCA uses a shared basis (FR-001–004). Distances are chain-aware (FR-020). |
| II. Replicates as unit of inference | Bootstrap resamples replicate means (FR-014/016). Equilibration maximises the effective sample size (FR-009). |
| III. Flag, do not discard | Equilibration trimming and the MM/PBSA start frame are recorded in output files. Distance chain ambiguity raises a warning instead of a silent merge. |
| IV. Tested behaviour | Every FR has a unit test with a known answer, CI enforces them (FR-006), and an e2e test covers real readers (FR-022). |
| V. Provenance | `equilibration.json`, the unit records, the joint PCA parameters and the run report (FR-023). |
| VI. Visualisation | Joint PCA uses the shared palette (APO grey) and a sequential overlap scale. Export failure is reported once (FR-021). |
| VII. Explicit changes | The stability default method changes to `chodera` when enabled. The ranking table gains columns. `SCIENTIFIC_VALIDITY.md`, `USAGE.md` and the templates are updated. |

**Gate result**: PASS (no violations; Complexity Tracking not needed).
**Post-design re-check**: PASS. The `pca_joint` stage is additive, and the loader split keeps
imports, so nothing breaks.

## Project Structure

### Documentation (this feature)

```text
specs/001-pipeline-hardening-round3/
├── plan.md
├── research.md
├── data-model.md
├── quickstart.md
├── contracts/
│   ├── config-schema.md     # new/changed YAML keys
│   └── outputs.md           # new output files and columns
├── checklists/requirements.md
└── tasks.md
```

### Source Code (repository root)

```text
gromacs_analysis/
├── pca/
│   ├── joint.py             # NEW: fit, covariance, projection, RMSIP (numpy core + MDAnalysis loader)
│   ├── joint_analyzer.py    # NEW: pca_joint stage (figures + CSV)
│   ├── config.py            # + PCAJointConfig
│   └── prep.py              # robust_threshold fallback (FR-019)
├── utils/
│   ├── equilibration.py     # NEW: statistical inefficiency, detect_equilibration
│   ├── units.py             # NEW: header unit parsing + conversion factors
│   ├── figure_export.py     # availability check + StaticExportUnavailable (FR-021)
│   └── report.py            # NEW: run report (FR-023)
├── md/                      # units, chodera stability, equilibration.json, replicate means, hbonds_lig metric
├── mmpbsa/                  # start_frame filtering of frame-level results
├── ranking/                 # bootstrap, reference exclusion
├── distance/                # chain-aware residues
├── pipeline/
│   ├── yaml_loader.py       # keeps public API; imports from split modules
│   ├── templates.py         # NEW: template YAML text
│   ├── stage_builders.py    # NEW: prolif/ttclust/pca_joint builders (+ _batch_systems)
│   ├── runner.py            # pca_joint stage + report generation
│   └── upstream.py          # protein–ligand H-bonds, rmsf first frame
.github/workflows/tests.yml  # NEW: CI
md_analysis_base.py, holo_md_analysis.py, apo_md_analysis.py  # deprecation warnings
tests/
├── test_round3_scientific.py    # NEW
├── test_round3_engineering.py   # NEW
└── test_end_to_end_trajectory.py# NEW (MDAnalysis required, ProLIF optional)
```

**Structure Decision**: Single Python package. New logic goes into focused modules beside the
code that consumes it. The YAML loader is split by moving its template text and the newer stage
builders out, and the moved names are re-exported, so old imports keep working.

## Complexity Tracking

No constitution violations.
