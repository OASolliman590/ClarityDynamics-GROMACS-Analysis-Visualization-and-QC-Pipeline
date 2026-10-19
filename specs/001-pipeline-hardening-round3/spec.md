# Feature Specification: Pipeline Hardening Round 3

**Feature Branch**: `ccr-deb837be-iavlqx` (spec directory `001-pipeline-hardening-round3`)

**Created**: 2026-10-19

**Status**: Draft

**Input**: User description: "Pipeline improvements round 3, covering 12 items. Scientific: (1) joint PCA in a common essential subspace across systems with RMSIP subspace overlap; (2) automatic nm/Å (and nm²/Å²) unit detection and conversion from .xvg headers, recorded in provenance; (3) statistically principled equilibration detection (Chodera's maximum effective-sample-size method), with the detected start written out so RMSF and MM/PBSA can use the same window; (4) ranking with uncertainty via bootstrap rank probabilities and pairwise P(A better than B), excluding APO/reference systems; (5) protein–ligand hydrogen bonds in the VMD helper and as an MD metric (per ligand for multi-ligand complexes); (6) robust terminal-detection threshold when the RMSF core is flat; (7) chain-aware residue identity in the distance module. Engineering: (8) CI running the tests on every pull request; (9) static figure export dependency check reported once with fallback; (10) real end-to-end test on a small generated trajectory; (11) a single HTML report per pipeline run; (12) deprecate legacy scripts and split the pipeline YAML loader while keeping the public API."

## User Scenarios & Testing *(mandatory)*

The users are computational chemists running replicate MD simulations of one protein with
several ligands (and an APO reference), who use the pipeline to compare ligands and report
results in papers.

### User Story 1 - Compare conformational landscapes on shared axes (Priority: P1)

A researcher wants to know whether ligands shift the protein's large-scale motions. Today each
system is analysed on its own principal components, so overlays are only qualitative. The
researcher runs a joint principal component analysis: all systems' Cα trajectories are fitted to
one common reference, analysed together, and every system is projected on the same components.
The output also states how similar each pair of systems' dominant motions are (subspace overlap).

**Why this priority**: It turns the most-requested cross-ligand figure from a caveated sketch
into a quantitative comparison; without it, PCA comparisons cannot be published.

**Independent Test**: Run the joint PCA on two synthetic systems with known, different dominant
motions and check that both are projected on identical axes, that their separation along PC1 is
recovered, and that the overlap between a system and itself is 1.

**Acceptance Scenarios**:

1. **Given** Cα trajectories for an APO system and two ligand systems with identical atom
   selections, **When** the joint PCA stage runs, **Then** one set of eigenvalues/eigenvectors is
   produced and every system's frames are projected on the same PC1–PC3 axes in one figure.
2. **Given** the joint PCA output, **When** the user opens the overlap table, **Then** each pair
   of systems has a subspace-overlap value between 0 and 1 for the configured number of PCs, with
   1.0 on the diagonal.
3. **Given** systems whose atom selections differ in size, **When** the stage runs, **Then** it
   stops with a message naming the systems and atom counts, and produces no figures.

---

### User Story 2 - Every pull request is tested automatically (Priority: P1)

A maintainer opens a pull request and the full test suite runs automatically, reporting pass or
fail on the pull request before anyone reviews or merges it.

**Why this priority**: The project's scientific guarantees live in its regression tests; they
only protect results if they run on every change.

**Independent Test**: Open a pull request and observe a test check that runs the suite and
reports its status.

**Acceptance Scenarios**:

1. **Given** a pull request against the default branch, **When** it is opened or updated,
   **Then** the test suite runs on at least two supported Python versions and its status is shown
   on the pull request.
2. **Given** a change that breaks a test, **When** CI runs, **Then** the check fails and names the
   failing test.

---

### User Story 3 - Correct units without manual bookkeeping (Priority: P2)

A researcher mixes metric files produced by GROMACS tools (lengths in nm, areas in nm²) and by
the VMD helper (Å, Å²). The pipeline reads the unit each file declares and converts values to
the unit shown on the plot axis, recording what it did.

**Why this priority**: A silent factor-of-10 (or 100) error in RMSD, RMSF, Rg, distances or SASA
invalidates every comparison built on it.

**Independent Test**: Feed one `.xvg` RMSD file declared in nm and one plain file in Å for two
systems and check both are reported in Å with values on the same scale, and that the conversion
is recorded.

**Acceptance Scenarios**:

1. **Given** an RMSD `.xvg` whose header declares nm and a metric displayed in Å, **When** it is
   processed, **Then** values are multiplied by 10 and the conversion (from-unit, to-unit, factor)
   is written to the statistics output.
2. **Given** a SASA `.xvg` declared in nm², **When** displayed in Å², **Then** values are
   multiplied by 100.
3. **Given** a file without a unit declaration, **When** processed, **Then** values are left
   unchanged and the output records the unit as assumed rather than detected.

---

### User Story 4 - Defensible equilibration cut-off applied consistently (Priority: P2)

A researcher enables automatic equilibration detection. The pipeline determines, per system, the
start of the production window by maximising the number of effectively independent samples,
trims all time-series metrics to that window, and writes the start time so that RMSF and MM/PBSA
can be computed over the same window.

**Why this priority**: Including equilibration biases means and uncertainties; ad-hoc thresholds
are unit-dependent and hard to defend in review.

**Independent Test**: Run detection on a synthetic series with a known relaxation followed by a
stationary segment and check the detected start falls at the end of the relaxation.

**Acceptance Scenarios**:

1. **Given** a series that relaxes over its first 20% then fluctuates stationarily, **When**
   detection runs, **Then** the detected start lies within the last part of the relaxation
   (not in the first 10%, not beyond 35%) and the effective sample size is reported.
2. **Given** detected start times, **When** the MD stage finishes, **Then** a machine-readable
   file lists, per system, the start time, start frame and effective sample size.
3. **Given** frame-level MM/PBSA results and a configured start frame, **When** MM/PBSA is
   processed, **Then** averages exclude earlier frames and the number of frames used is reported.
4. **Given** a VMD helper task with a first frame, **When** the script is generated, **Then** RMSF
   is computed only from that frame onward.

---

### User Story 5 - Rankings that show their uncertainty (Priority: P2)

A researcher ranks ligands. Besides the point ranking, the output gives each ligand's probability
of being ranked first, its median rank with an interval, and, for every pair, the probability
that one outperforms the other. The APO/reference system never appears in the ranking.

**Why this priority**: With 3–5 ligands and replicate noise, point rankings overstate certainty;
reviewers ask for this.

**Independent Test**: Rank two ligands whose ΔTOTAL differ by far more than their uncertainty and
two that overlap; the first pair yields a probability near 1, the second near 0.5.

**Acceptance Scenarios**:

1. **Given** ligands with replicate-level values or a mean and SEM, **When** ranking runs,
   **Then** a table reports P(rank 1), median rank and a 95% rank interval for each ligand, and
   a second table reports P(A better than B) for every pair.
2. **Given** MD statistics that include an APO system, **When** ranking runs, **Then** the APO
   system is excluded and this is logged.
3. **Given** a fixed random seed, **When** ranking runs twice, **Then** results are identical.

---

### User Story 6 - Protein–ligand hydrogen bonds (Priority: P3)

A researcher wants the number of hydrogen bonds between the protein and each ligand over time,
not only hydrogen bonds within the protein.

**Why this priority**: Protein–ligand hydrogen bonds are the relevant binding descriptor; the
intra-protein count is mostly a stability descriptor.

**Independent Test**: Generate the VMD helper for a two-ligand complex and check one protein–ligand
hydrogen-bond output per ligand is requested; process matching files and check one plot per ligand.

**Acceptance Scenarios**:

1. **Given** a single-ligand task, **When** the VMD helper is generated, **Then** it writes a
   protein–ligand hydrogen-bond file in addition to the intra-protein one.
2. **Given** a complex with ligands INH and NAD, **When** the MD stage runs on the generated files,
   **Then** a protein–ligand hydrogen-bond plot exists for each ligand.

---

### User Story 7 - Robust terminal detection on flat profiles (Priority: P3)

When the RMSF core is nearly constant, terminal overhangs must still be fully detected instead of
inflating the threshold and being under-trimmed.

**Why this priority**: Under-trimming leaks flexible termini into the PCA core and RMSF summaries.

**Independent Test**: A flat core with 6 elevated N-terminal and 10 elevated C-terminal residues
is trimmed by at least 6 and 10 residues respectively.

**Acceptance Scenarios**:

1. **Given** an exactly flat core with elevated termini, **When** terminals are detected,
   **Then** all elevated residues are flagged and no core residue is flagged.
2. **Given** a profile whose spread is not degenerate, **When** terminals are detected, **Then**
   results are unchanged from the current behaviour.

---

### User Story 8 - Multi-chain proteins in distance analysis (Priority: P3)

For a dimer where both chains contain residue 45, distances to "A:45" and "B:45" are reported
separately and labelled with their chain.

**Why this priority**: Merging residues of different chains silently produces wrong distances.

**Independent Test**: A two-chain system with duplicated residue numbers yields separate,
chain-labelled residue entries.

**Acceptance Scenarios**:

1. **Given** residues specified as chain-qualified identifiers, **When** distances are computed,
   **Then** only the residue in the named chain is used and its label includes the chain.
2. **Given** an unqualified residue number present in two chains, **When** distances are
   computed, **Then** both chains' residues are reported as separate entries and a warning
   explains the ambiguity.

---

### User Story 9 - Clear behaviour when image export is unavailable (Priority: P2)

On a machine without the browser needed for static image export, the user sees one clear message
explaining how to enable it, interactive figures are still written, and static images fall back
to an alternative renderer where one exists, instead of one error per figure.

**Why this priority**: Silent or repeated failures waste time and lead to missing figures in
reports.

**Independent Test**: Simulate an unavailable exporter and check that exactly one warning is
logged across many figures and that fallback images are still produced where available.

**Acceptance Scenarios**:

1. **Given** static export is unavailable, **When** a stage saves ten figures, **Then** one
   warning with the remedy is logged and interactive files are still written for all ten.
2. **Given** a figure that has an alternative renderer, **When** static export is unavailable,
   **Then** its SVG is produced by the alternative renderer.

---

### User Story 10 - End-to-end confidence on real trajectory files (Priority: P2)

Maintainers can run a test that writes a small, physically plausible protein–ligand trajectory
to disk and runs QC, distance analysis, joint PCA and (when installed) ProLIF on it.

**Why this priority**: Unit tests with stand-ins cannot catch integration failures with the real
trajectory readers.

**Independent Test**: The test passes on a clean environment with the trajectory library
installed and skips ProLIF parts when ProLIF is not installed.

**Acceptance Scenarios**:

1. **Given** the generated trajectory, **When** QC runs, **Then** frame count, timestep and
   ligand–protein distance are reported correctly.
2. **Given** the generated trajectory with the ligand drifting away, **When** the distance stage
   runs, **Then** the reported distance increases over time as constructed.

---

### User Story 11 - One report per run (Priority: P3)

After a pipeline run, the researcher opens a single HTML page listing every executed stage with
its figures, key tables and QC warnings.

**Why this priority**: Outputs are currently spread across many folders, making review and sharing
slow.

**Independent Test**: Run a pipeline with two stages and check the report links every produced
figure and shows the QC warnings.

**Acceptance Scenarios**:

1. **Given** a completed run, **When** the user opens the report, **Then** each stage has a
   section with its figures (interactive or image), its summary tables, and its status.
2. **Given** QC warnings, **When** the report is generated, **Then** they appear at the top.
3. **Given** a stage that failed, **When** the report is generated, **Then** the failure and its
   message are shown.

---

### User Story 12 - Maintainable code base (Priority: P3)

Maintainers can find stage configuration code by stage, and users of the old standalone scripts
are told to move to the package.

**Why this priority**: Reduces the risk of future changes; no direct scientific effect.

**Independent Test**: Existing configurations load identically before and after the
reorganisation; running a legacy script prints a deprecation notice.

**Acceptance Scenarios**:

1. **Given** existing pipeline YAML files and imports from the configuration loader, **When**
   loaded after the reorganisation, **Then** they produce identical configurations.
2. **Given** a user runs a legacy standalone script, **When** it starts, **Then** it warns that
   it is deprecated and names the replacement.

### Edge Cases

- Joint PCA with only one system: runs, overlap table has a single entry, a warning notes that
  comparison needs two or more systems.
- Joint PCA with very different trajectory lengths: every frame contributes equally by default; a
  per-system frame cap can balance systems and is recorded.
- Unit header present but unrecognised (e.g. "a.u."): value left unchanged, unit recorded as
  declared but not converted, warning logged.
- Equilibration detection on a series shorter than the minimum length: no trimming, recorded as
  "insufficient data".
- Bootstrap with a ligand missing one metric: that metric is ignored for that ligand in every
  bootstrap sample and the number of scored metrics is reported.
- Ranking where every system is a reference: stops with a clear message.
- Distance residue identifier naming a chain that does not exist: error naming the identifier.
- HTML report when no figures were produced: report still lists stages and states that no
  figures exist.

## Requirements *(mandatory)*

### Functional Requirements

**Joint PCA (US1)**

- **FR-001**: The system MUST provide a joint PCA analysis that takes multiple systems' Cα
  coordinate trajectories, fits all frames to one common reference, and computes one set of
  principal components from all systems together.
- **FR-002**: The system MUST project every system's frames onto the shared components and write
  per-system projections, eigenvalues and explained-variance fractions.
- **FR-003**: The system MUST compute the subspace overlap (root mean square inner product) over
  a configurable number of components between every pair of systems' individual subspaces and
  between each system and the joint subspace.
- **FR-004**: The system MUST produce a shared-axes projection figure (systems coloured with the
  project palette, APO grey), a per-system density figure on PC1–PC2, and an overlap heatmap.
- **FR-005**: The system MUST refuse to run when systems have different numbers of selected atoms,
  naming the systems and counts.

**CI (US2)**

- **FR-006**: Every pull request and push to the default branch MUST run the full test suite on
  at least two supported Python versions and report status on the pull request.

**Units (US3)**

- **FR-007**: The MD stage MUST detect the y-axis unit declared in `.xvg` headers (nm, nm², Å, Å²,
  including GROMACS formatting codes) and convert to the metric's display unit.
- **FR-008**: The statistics output MUST record, per system and metric, the detected unit, the
  display unit, the factor applied and whether the unit was detected or assumed.

**Equilibration (US4)**

- **FR-009**: The stability option MUST offer a method that selects the production start
  maximising the effective number of uncorrelated samples, and this MUST be the default method
  when stability detection is enabled; the previous threshold method stays available.
- **FR-010**: The MD stage MUST write per-system start time, start index and effective sample size
  to a machine-readable file.
- **FR-011**: The MM/PBSA stage MUST accept a start frame and, when frame-level results are
  available, exclude earlier frames from averages and report frames used.
- **FR-012**: The VMD helper MUST accept a first frame from which RMSF is computed.

**Ranking (US5)**

- **FR-013**: Ranking MUST compute, by bootstrap with a configurable number of resamples and seed,
  each ligand's probability of rank 1, median rank and 95% rank interval, and pairwise
  probabilities of one ligand outranking another.
- **FR-014**: Bootstrap samples MUST use replicate-level values when available, otherwise a normal
  distribution with the reported mean and standard error.
- **FR-015**: Reference systems (configured, the MD stage's APO systems, or named exactly "APO",
  case-insensitive) MUST be excluded from ranking.
- **FR-016**: The MD statistics output MUST include per-replicate means so ranking can resample
  them.

**Hydrogen bonds (US6)**

- **FR-017**: The VMD helper MUST write a protein–ligand hydrogen-bond time series for the ligand
  selection and for each listed ligand.
- **FR-018**: The MD stage MUST provide a protein–ligand hydrogen-bond metric, split per ligand for
  multi-ligand complexes like the other ligand metrics.

**Terminal detection (US7)**

- **FR-019**: When the median absolute deviation of the RMSF profile is zero, the robust threshold
  MUST use a scale that terminal outliers cannot inflate (interquartile-based, then a small
  fraction of the median) instead of the standard deviation.

**Multi-chain distances (US8)**

- **FR-020**: The distance module MUST identify residues by chain/segment and number, accept
  chain-qualified identifiers ("A:45"), label residues with their chain when more than one chain
  is present, and warn when an unqualified number matches several chains.

**Export (US9)**

- **FR-021**: The system MUST determine once per run whether static export works; when it does
  not, it MUST log a single warning with remedies, skip further attempts, and use the alternative
  renderer where one exists.

**End-to-end test (US10)**

- **FR-022**: The test suite MUST include a test that writes a small topology and multi-frame
  trajectory and runs QC, distance and joint PCA on them, and ProLIF when available.

**Report (US11)**

- **FR-023**: The pipeline MUST write one HTML report per run in the output root listing stages
  with status, figures, summary tables (first rows) and QC warnings, using relative links.

**Cleanup (US12)**

- **FR-024**: Legacy standalone scripts MUST emit a deprecation warning naming the replacement.
- **FR-025**: The pipeline configuration loader MUST be split so the template text and the builders
  for the ProLIF, TTClust and joint-PCA stages live in separate modules, while all existing imports
  from the loader keep working.

### Key Entities

- **Joint PCA result**: shared eigenvalues/eigenvectors, per-system projections, explained
  variance, overlap matrix, parameters (selection, reference, number of components, frame caps).
- **Unit conversion record**: system, metric, declared unit, display unit, factor, detected/assumed.
- **Equilibration record**: system, start time, start index, effective sample size, statistical
  inefficiency, method.
- **Rank distribution**: ligand, P(rank 1), median rank, rank interval, number of scored metrics;
  pairwise table of P(A better than B).
- **Run report**: stage list with status, figure paths, table excerpts, QC warnings.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Joint PCA on synthetic systems with a known displacement along one direction
  recovers that direction as PC1 (absolute cosine ≥ 0.95) and gives self-overlap of 1.0 ± 1e-6.
- **SC-002**: 100% of pull requests show a test status before merge.
- **SC-003**: Lengths declared in nm and in Å for the same structure are reported within 1% of
  each other after processing.
- **SC-004**: On synthetic relaxation-plus-noise series, the detected production start lies within
  ±15% of the series length of the true end of relaxation in at least 9 of 10 seeds.
- **SC-005**: For two ligands separated by 10 standard errors, P(better) ≥ 0.99; for identical
  ligands, P(better) is between 0.4 and 0.6.
- **SC-006**: When image export is unavailable, a run producing 10 figures logs exactly one export
  warning.
- **SC-007**: A reviewer can open one file per run and reach every figure in at most two clicks.
- **SC-008**: All previously passing tests still pass and existing example configurations load to
  identical configurations.

## Assumptions

- Joint PCA uses Cα atoms with equal weights, fits to the average structure after an initial fit to
  the first frame of the first system (two iterations), and does not depend on GROMACS being
  installed; GROMACS-based per-system PCA remains available.
- Overlap is computed on the first 10 components by default.
- Display units are taken from the metric's axis label (Å, Å²); metrics without a length or area
  unit are never converted.
- Equilibration detection uses an in-house implementation of the published method so no new
  required dependency is added.
- The bootstrap default is 2,000 resamples with seed 0.
- CI installs the optional trajectory library so the end-to-end test runs; ProLIF/TTClust are not
  installed in CI and their tests skip or use stand-ins.
- Pushing CI workflow files requires repository permissions that the automation account may lack;
  if so, the workflow file is delivered for a maintainer to add.
- The configuration loader split is internal; no configuration keys are renamed.
