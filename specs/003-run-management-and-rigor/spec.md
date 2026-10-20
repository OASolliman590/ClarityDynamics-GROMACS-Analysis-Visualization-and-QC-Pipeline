# Feature Specification: Run Management, Input Preparation and Statistical Rigour

**Feature Branch**: `ccr-deb837be-iavlqx` (spec directory `003-run-management-and-rigor`)

**Created**: 2026-10-20

**Status**: Draft

**Input**: User description: "Implement all of them accordingly". This covers the ten follow-ups
proposed after feature 002:

1. Follow a run after it starts.
2. Produce missing inputs from the interface.
3. Run per system or ligand in parallel on HTCondor.
4. Edit stage settings without YAML.
5. Skip unchanged work.
6. Check the real cluster environment.
7. Report sampling sufficiency.
8. Compare ligands statistically.
9. Modern packaging.
10. Code-quality gates in CI, plus run provenance in the report.

## User Scenarios & Testing *(mandatory)*

Users run `clarity-setup` against simulation folders on a workstation, an SSH host or an
HTCondor pool reached with condorctl (for example NMRbox). Feature 002 lets them discover,
name and start a run. This feature covers what happens after that: preparing missing inputs,
running fast, following runs, and trusting the numbers.

### User Story 1 - Follow my runs and get the report back (Priority: P1)

After starting or submitting a run, the user sees its state (running, idle/queued, held,
completed, failed) and the last lines of its log. When it finishes, the user downloads the
report folder to their computer, from the interface or the command line.

**Independent Test**: Start a run through a fake host. `status` shows running and the log tail.
After the fake job finishes, `status` shows completed, and `fetch` downloads `report.html` and
the figures it links to into a local folder.

**Acceptance Scenarios**:
1. **Given** a background run on a host, **When** the user asks for status, **Then** they see
   whether the process is alive and the last N log lines.
2. **Given** an HTCondor submission (single job or workflow), **When** the user asks for status,
   **Then** they see each job's state, counted by state, read from `condor_q`. Once jobs leave
   the queue, state comes from `condor_history`.
3. **Given** a finished run, **When** the user fetches results, **Then** the report and every
   file under the output folder are copied locally, and the local path is shown.
4. Every run is recorded in the project file, so status survives closing the interface.

### User Story 2 - Produce inputs the scan says are missing (Priority: P1)

For each replicate the scan reports what is missing:

- a processed (PBC-corrected) trajectory;
- MD metric files;
- MM/PBSA results.

The user asks the tool to prepare them. It writes one reviewed script per replicate and task,
using the tools already in the environment: GROMACS `trjconv`, the existing VMD metrics helper,
and gmx_MMPBSA. It then runs the scripts locally, on the host, or as HTCondor jobs.

**Independent Test**: For a fixture project with tpr files and raw trajectories, `prepare`
writes:
- a trjconv script per replicate without a processed trajectory;
- a VMD script per replicate without metrics;
- an MM/PBSA script plus index/input files per holo replicate without results.

Running them through a fake host issues one start or submit per batch, and the dry run lists the
scripts without running anything.

**Acceptance Scenarios**:
1. The scripts use named groups (`Protein`, `System`, and the ligand residue selection) rather
   than hard-coded index numbers. The MM/PBSA index file is built with `gmx select`, so group
   numbers are known.
2. Replicates that already have an input are not touched unless the user forces it.
3. A replicate missing a prerequisite (e.g. no `.tpr` for trjconv) is listed with the reason,
   not silently skipped.
4. After preparation, a rescan shows the new inputs and the related stages become available.

### User Story 3 - Run analysis in parallel on HTCondor (Priority: P1)

The user submits the analysis as an HTCondor workflow:

- **Parallel jobs**:
  - one job per system and ligand for ProLIF, distances and TTClust;
  - one job each for QC, MD, MM/PBSA and joint PCA.
- **Final assembly job**: runs after all of them. It builds the cross-ligand comparison, the
  ranking and the report, and reuses everything the parallel jobs produced.

**Independent Test**: For the fixture project, building the workflow gives:
- one unit job per ProLIF ligand run, distance entry and TTClust system, each with its own
  validated configuration;
- one job per single stage;
- a final job that depends on all of them.

Submitting through a fake host calls `condor_submit_dag` once and records the cluster.

**Acceptance Scenarios**:
1. Every unit configuration loads with the real pipeline loader.
2. Each unit writes to exactly the output folder the full run would use, so the final job finds
   the results.
3. The final job re-runs nothing that a unit already completed (cache hit). It still runs
   anything that failed or was not produced.

### User Story 4 - Change analysis settings without YAML (Priority: P2)

From the interface or the command line, the user changes the documented settings of each stage:

- stride;
- contact cutoff;
- PCA selection and component count;
- bootstrap count;
- TTClust clustering method;
- MM/PBSA start frame;
- convergence thresholds.

Values are type-checked and stored in the project file, and the generated YAML uses them.

**Independent Test**: Set `prolif_batch.stride=5` and `ranking.n_bootstrap=500` via `init --set`.
The generated YAML carries both values and validates. An invalid value (e.g. `stride=0` or text
where a number is expected) is rejected with a message.

### User Story 5 - Do not redo unchanged work (Priority: P2)

Re-running after changing one system re-runs only what depends on it. Caching is on by default
for generated configurations and works per system or ligand within batch stages. A change to a
stage's settings invalidates that stage. `--no-cache` forces a full run.

**Independent Test**: Run a two-system batch stage twice. The second run skips both units.
Touching one system's input re-runs only that unit. Changing a stage setting re-runs the stage.

### User Story 6 - Check the cluster before I trust it (Priority: P2)

`clarity-setup doctor` checks the connection and reports what is available, with a remedy for
each missing item:

- write access to the project root;
- Python and the pipeline executable;
- `gmx`, `vmd`, `gmx_MMPBSA`;
- `condor_q`, `condor_submit`, `condor_submit_dag`;
- free disk space;
- whether the HTCondor execute nodes can see the project folder (shared filesystem).

**Independent Test**: Against a fake host missing `vmd`, doctor reports every tool's state, marks
VMD as missing with an install hint, and exits non-zero only when a tool required by the chosen
stages or run mode is missing.

### User Story 7 - Know whether sampling is sufficient (Priority: P1)

For every system, metric and replicate the pipeline reports:

- the block-averaged standard error and effective sample size;
- whether the first and second halves of the production window agree;
- for joint PCA, the cosine content of the first principal components per replicate, and the
  overlap (RMSIP) between the principal-component spaces of the two halves of each system's
  frames.

Replicates that fail a check are flagged in a `convergence.csv` table and at the top of the run
report.

**Independent Test**:
- A white-noise series passes the checks. A random walk is flagged by half-agreement and by
  cosine content near 1.
- A series that is too short has fewer effective samples than the threshold and is flagged.
- Joint PCA on a converged fixture gives cosine content below 0.5 for PC1.

### User Story 8 - Say which ligand is better, with statistics (Priority: P2)

The ranking stage writes `pairwise_tests.csv`. For each ranked metric, including MM/PBSA ΔTOTAL
when per-replicate values exist, it compares every pair of ligands with:

- a Welch t-test on replicate means;
- an exact permutation test when there are few replicates;
- Hedges' g effect size;
- Holm-adjusted and Benjamini-Hochberg-adjusted p-values, across all comparisons of each metric.

**Independent Test**: Two ligands with well-separated replicate means give an adjusted p-value
below 0.05 and a large |g|. Identical distributions give p ≈ 1. With a single replicate per
ligand, the rows are reported as "insufficient replicates" rather than tested.

### Engineering Stories (Priority: P3)

- **US9 Packaging**:
  - `pyproject.toml` (PEP 621) replaces the metadata in `setup.py`, and editable installs work
    with modern pip.
  - A `environment.yml` describes a conda environment with GROMACS-side tools, kept optional.
- **US10 Quality and provenance**:
  - CI runs `ruff` and `mypy` on the setup package, and coverage with a threshold.
  - The run report shows provenance: package version, git commit (when available), Python and
    key library versions, host, start/end time, and the exact configuration that ran.

### Edge Cases

- A host where `condor_q` reports no job and `condor_history` is unavailable: state is
  "unknown", with the reason.
- A run whose log does not exist yet (still queued): status says so instead of failing.
- Fetching a large output folder: one archive is made remotely and downloaded once, not one file
  at a time.
- A unit job fails: the final job still runs, because the workflow uses retry-free `PARENT/CHILD`
  with the final job marked to always run. It reports the failed units, which have no cache
  entry.
- Settings for a stage that is not enabled are kept but have no effect.
- Statistical tests on metrics with no defined direction are still reported, with "two-sided".
- Convergence on a series shorter than 8 points: "too short", not an error.

## Requirements *(mandatory)*

### Functional Requirements

**Runs (US1)**
- **FR-001**: The system MUST record every started run in the project file. Each record holds:
  where it runs, the identifiers (PID or cluster), log, output folder, start time and the set of
  submitted jobs.
- **FR-002**: The system MUST report each run's state:
  - a host run: from the process table and its log;
  - an HTCondor run: from `condor_q`, falling back to `condor_history`, mapping status codes
    1–7 to names.
- **FR-003**: The system MUST fetch a finished run's output folder as one archive and unpack it
  locally.

**Inputs (US2)**
- **FR-004**: The system MUST detect, per replicate, the files needed for preparation: `.tpr`,
  `.top`, the raw trajectory, and whether a processed trajectory exists.
- **FR-005**: The system MUST generate preparation scripts per replicate for:
  - processed trajectory (`gmx trjconv`, named groups);
  - MD metrics (the existing VMD helper);
  - MM/PBSA (index via `gmx select`, a default `mmpbsa.in`, `gmx_MMPBSA`).

  It MUST list replicates skipped for missing prerequisites, with the reason.
- **FR-006**: The system MUST run preparation locally, on the host, or as one HTCondor job per
  script, and offer a dry run.

**Parallel workflow (US3)**
- **FR-007**: The system MUST split the generated configuration into unit configurations:
  - per entry for the batch stages (ProLIF, distances, TTClust);
  - per stage for the single stages.

  Output folders MUST be identical to the full run's.
- **FR-008**: The system MUST write a DAGMan file with every unit as a parent of the final job,
  and submit it with `condor_submit_dag`.
- **FR-009**: The generator MUST list one ProLIF entry per ligand, with stable names, so units
  map one-to-one to runs.

**Settings (US4)**
- **FR-010**: The system MUST provide a schema of editable settings per stage, giving name, type,
  bounds and description. Values are stored per project and applied by the generator.
- **FR-011**: The interface MUST let the user edit settings per stage; the CLI MUST accept
  `--set stage.key=value`.

**Caching (US5)**
- **FR-012**: The cache signature MUST include a fingerprint of the stage configuration, not
  only the input files.
- **FR-013**: Batch stages MUST cache per unit, using each unit's own output folder.
- **FR-014**: Generated configurations MUST enable caching. `gromacs-pipeline --no-cache` MUST
  disable it for a run.

**Doctor (US6)**
- **FR-015**: `doctor` MUST check the connection, write access, executables, HTCondor tools,
  free disk space and shared-filesystem visibility. It MUST print a table with a remedy per
  failure and use a meaningful exit code.

**Convergence (US7)**
- **FR-016**: The MD stage MUST write `data/convergence.csv` with these columns:
  `System, Metric, Replicate, N, SEM_block, N_eff, Half_Diff, Half_Z, Flags`.
- **FR-017**: Joint PCA MUST write `data/convergence.csv` with cosine content of PC1–3 per
  system and replicate, and a half-split RMSIP per system.
- **FR-018**: The run report MUST list flagged items at the top, next to the QC warnings.

**Statistics (US8)**
- **FR-019**: Ranking MUST write `data/pairwise_tests.csv` with these columns:
  `Metric, LigandA, LigandB, n_A, n_B, Mean_A, Mean_B, Diff, Welch_t, Welch_p, Perm_p,
  Hedges_g, p_holm, p_bh, Direction, Note`.

**Engineering (US9, US10)**
- **FR-020**: Packaging MUST use `pyproject.toml` with the existing extras and console scripts.
- **FR-021**: CI MUST run ruff, mypy (setup package) and pytest with coverage ≥ the threshold.
- **FR-022**: The run report MUST include a provenance section and the configuration that ran.

### Key Entities

- **RunRecord**: `id`, `where` (local|host|condor|dag), `ids` (PID or clusters), `log`,
  `output_dir`, `started`, `jobs` (name → cluster), `kind` (analysis|prepare).
- **PrepTask**: `system`, `replicate`, `kind` (trjconv|metrics|mmpbsa), `script`, `inputs`,
  `outputs`, `skipped_reason`.
- **WorkUnit**: `name`, `stage`, `config` (dict), `parents`.
- **StageSetting**: `stage`, `key`, `type`, `default`, `min`, `max`, `choices`, `help`.

## Success Criteria *(mandatory)*

- **SC-001**: A user goes from raw simulations (tpr + raw trajectories) to a finished report
  without editing a file by hand: `prepare` → rescan → `run --where dag` → `status` → `fetch`.
- **SC-002**: For N ligand runs, the parallel workflow's wall time is bounded by the slowest
  unit plus the final job, not the sum of all units. The structure is verified by the workflow
  having N independent parents.
- **SC-003**: A second identical run executes zero analysis units, all cache hits. Changing one
  system re-executes only the units for that system.
- **SC-004**: Every flagged convergence problem appears in the report without opening a CSV.
- **SC-005**: All new behaviour is covered by tests. The suite stays green on Python 3.10 and 3.12.

## Assumptions

- HTCondor pools share the project folder with the submit host (as on NMRbox). There is no
  file transfer.
- GROMACS, VMD and gmx_MMPBSA are installed by the user. The tool generates and runs scripts
  but does not install them.
- MM/PBSA preparation uses the first included ligand of each holo system as the ligand group.
  Users who need another choice exclude the other ligands or edit the generated script.
- Live cluster validation (US6) is done by the user with `doctor`. The tests use a fake host.
