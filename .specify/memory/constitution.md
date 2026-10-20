# ClarityDynamics Constitution

## Core Principles

### I. Scientific Correctness Before Convenience
Every number the pipeline reports must have a defined physical meaning, an explicit unit and
a documented method. Units are read from the data where the data declares them (e.g. `.xvg`
headers) and are never assumed silently; when they must be assumed, a warning is logged and the
assumption is recorded in the output metadata. Quantities are labelled for what they are
(e.g. MM/PBSA ΔTOTAL is not "binding affinity"). Periodic boundary conditions are respected in
every geometric computation.

### II. Replicates Are the Unit of Inference
Uncertainty of a mean is reported as a standard error that reflects independent samples:
between-replicate SD/√n when replicates exist, an autocorrelation-aware estimate (block
averaging) otherwise. Correlated frames, residues or other non-independent units are never
treated as independent samples in statistical tests. Comparisons and rankings must expose
their uncertainty, not only point estimates.

### III. Flag, Do Not Silently Discard
Data are never removed without the removal being visible in plots and recorded in output files.
Filtering that can remove real signal (outliers, terminal residues, equilibration) is opt-in
or, when on by default, is shown (faded/shaded regions, QC warnings) and logged with the rule
and parameters used.

### IV. Tested Scientific Behaviour (NON-NEGOTIABLE)
Every scientific behaviour (statistics, parsing of tool outputs, unit handling, geometry,
selection logic) has a regression test with a known answer. Bug fixes add a test that fails on
the old behaviour. The full test suite runs in CI on every pull request and must pass before
merge.

### V. Reproducible, Provenance-Tracked Outputs
Each stage writes the data behind its figures (CSV/JSON) next to the figures, together with the
parameters that produced them. Runs are deterministic for fixed inputs (seeded randomness).
Outputs of one run can be browsed from a single report.

### VI. Consistent, Accessible Visualisation
Systems keep one colour across all stages (APO neutral grey, ligands in a fixed
colour-blind-checked order); hues are never cycled without a secondary encoding. No dual-axis
charts. Perceptually uniform scales for continuous data. Static exports preserve figure size
and use scale only as resolution. Missing export dependencies are reported once and clearly.

### VII. Explicit Behaviour Changes
Changing a default that alters results requires a changelog/PR note and documentation update
(`SCIENTIFIC_VALIDITY.md`, `USAGE.md`, templates). Legacy entry points that bypass these
principles are deprecated with a warning before removal.

## Technical Constraints

- Python ≥ 3.9 package `gromacs_analysis`; core dependencies numpy, scipy, pandas, plotly,
  matplotlib, PyYAML. MDAnalysis, ProLIF, RDKit, TTClust and pymbar are optional and every
  feature that needs them must degrade with a clear message when they are absent.
- GROMACS/VMD/gmx_MMPBSA are external: the toolkit generates or parses their inputs/outputs and
  only executes them when a stage is explicitly configured to do so.
- Configuration is YAML-driven; disabled stages must load without touching the file system.

## Development Workflow

- Features are specified with Spec Kit (`specs/<NNN-feature>/spec.md`, `plan.md`, `tasks.md`)
  before implementation.
- Every PR states which principles it touches, which defaults change, and how it was tested.
- Figures added or changed in a PR are rendered once and checked visually.

## Governance

This constitution supersedes ad-hoc conventions. Amendments are made by PR that updates this
file, increments the version (MAJOR: principle removed/redefined; MINOR: principle or section
added/expanded; PATCH: wording) and lists affected templates and docs. Reviewers verify
compliance with Principles I–VII for every PR that touches analysis code.

**Version**: 1.0.0 | **Ratified**: 2026-10-19 | **Last Amended**: 2026-10-19
