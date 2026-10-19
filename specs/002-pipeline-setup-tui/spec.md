# Feature Specification: Pipeline Setup TUI (Discover → Name → Run)

**Feature Branch**: `ccr-deb837be-iavlqx` (spec directory `002-pipeline-setup-tui`)

**Created**: 2026-10-19

**Status**: Draft

**Input**: User description: "Work on pipeline usability. The YAML is a headache. Use condorctl or any other SSH connection to detect systems, rename them, rename ligands, with an intuitive TUI tool that eases use of the pipeline."

## User Scenarios & Testing *(mandatory)*

Users are computational chemists whose simulation folders live on a workstation or on a
remote cluster (often an HTCondor pool such as NMRbox, reached with condorctl). Today they
hand-write a 450-line YAML file: they list every replicate directory, topology and trajectory,
copy ligand residue names out of topology files, and keep system names consistent across stages.

### User Story 1 - Discover my simulations automatically (Priority: P1)

The user points the tool at a project folder, either local or on a remote host over SSH. The tool
lists the simulation systems it found, with:

- replicates grouped together;
- topology and trajectory files identified;
- APO (ligand-free) systems recognised;
- ligand residues detected from each topology;
- existing analysis inputs noted (MD metric files, MM/PBSA results, PCA folders).

**Why this priority**: Discovery removes most of the YAML work and is useful even without the
interactive interface (`scan` prints the result).

**Independent Test**: Scan a folder containing `HOLO_LigA_1..3`, `HOLO_LigB_1..2` and `APO_1..3`
directories with topologies and trajectories. Expect three systems with 3, 2 and 3 replicates,
APO flagged, and the ligand residue names listed for each holo system.

**Acceptance Scenarios**:

1. **Given** a local project root, **When** the user runs a scan, **Then** each system is listed
   with its replicate count, topology, trajectories, APO flag and detected ligand residue names.
2. **Given** an SSH host (or a condorctl config) and a remote root, **When** the user scans, **Then**
   the same result is produced. Discovery uses a small number of remote commands rather than one
   per file, and nothing is downloaded beyond the topology residue names.
3. **Given** a complex with two ligand residues (inhibitor plus cofactor), **When** scanned, **Then**
   both are listed for that system.
4. **Given** directories that do not look like simulations (no trajectory), **When** scanned,
   **Then** they are ignored and the scan reports how many directories were skipped.

---

### User Story 2 - Name systems and ligands the way I want them in figures (Priority: P1)

The user renames systems (e.g. `HOLO_LigA` → `Isoniazid`) and ligands (residue `INH` →
`Isoniazid`, `NAD` → `NADH`). They also mark which system is the APO reference and exclude
systems or ligands they don't want analysed. The names are used consistently in every stage's
outputs.

**Why this priority**: Correct, consistent names are what appear in papers. Renaming by hand
across YAML sections is error-prone.

**Independent Test**: Rename one system and one ligand, then generate the configuration. The new
names appear as system names in MD, ProLIF, TTClust and joint PCA, while residue selections still
use the original residue name.

**Acceptance Scenarios**:

1. **Given** discovered systems, **When** the user renames a system, **Then** every generated stage
   uses the new name, and folder paths are unchanged.
2. **Given** a ligand residue `INH` renamed to `Isoniazid`, **When** the configuration is
   generated, **Then** ligand labels use `Isoniazid` and selections still use `resname INH`.
3. **Given** two systems given the same name, **When** the user tries to save, **Then** the tool
   refuses and names the conflict.
4. **Given** a system marked as excluded, **When** generating, **Then** it appears in no stage.

---

### User Story 3 - Choose analyses and get a valid configuration without writing YAML (Priority: P1)

The user ticks the analyses they want (QC, MD, ProLIF, TTClust, joint PCA, ranking and so on)
and sees a preview. They save a pipeline configuration that the existing pipeline accepts
unchanged, plus a small project file that remembers the naming choices for next time.

**Why this priority**: This closes the loop. Discovery plus naming becomes a runnable pipeline.

**Independent Test**: Generate a configuration from a scanned project. It loads and validates
with the existing pipeline loader (dry run) with no edits, and reopening the project file
restores every rename.

**Acceptance Scenarios**:

1. **Given** a scanned and named project, **When** the user saves, **Then** a pipeline YAML and a
   project file are written. The YAML passes the pipeline's own validation (dry run).
2. **Given** only stages whose inputs exist are suggested, **When** a stage's inputs are missing
   (e.g. no MM/PBSA results), **Then** that stage is unticked by default and shows why.
3. **Given** a saved project file, **When** the user reopens it and rescans, **Then** previous
   renames and choices are kept, and new systems are added as new entries.

---

### User Story 4 - Interactive terminal interface (Priority: P2)

The user runs one command and works through a keyboard-driven terminal interface:

1. Connect (local folder, SSH host, or condorctl profile).
2. Review the systems table: rename, toggle APO, exclude.
3. Review the ligands table: rename, exclude.
4. Pick analyses.
5. Preview the YAML.
6. Save, run, or submit.

It works over SSH sessions and in plain terminals.

**Why this priority**: This is the most intuitive front end. Stories 1–3 already deliver value
through non-interactive commands.

**Independent Test**: Drive the interface headlessly. Load a project, rename a system in the
table, toggle a stage and save. Expect the saved YAML to reflect each action.

**Acceptance Scenarios**:

1. **Given** the interface is open on a scanned project, **When** the user edits a name cell and
   presses Enter, **Then** the name updates, and the preview refreshes on the next view.
2. **Given** the user presses the save key, **Then** the files are written and their paths are
   shown.
3. **Given** a connection or scan error, **Then** it is shown inside the interface with the cause,
   and the interface stays usable.

---

### User Story 5 - Run where the data lives (Priority: P2)

From the interface or the command line, the user runs the pipeline on the configured project. Run
locally for local data. For remote data, upload the configuration and either start it on the
remote host or submit it as a job through condorctl.

**Why this priority**: Remote data should not have to be copied back. condorctl already provides
the job machinery.

**Independent Test**: For a remote project, the run action uploads the configuration next to the
data and issues exactly one remote start command (or one condorctl submission). The returned job
identifier or log path is shown.

**Acceptance Scenarios**:

1. **Given** a local project, **When** the user chooses Run, **Then** the pipeline runs and the path
   of the run report is shown.
2. **Given** a remote project, **When** the user chooses Run on host, **Then** the configuration is
   uploaded and the pipeline starts in the background on the host, with a log file path shown.
3. **Given** a remote project and condorctl installed, **When** the user chooses Submit, **Then** a
   job is submitted and its cluster ID is shown.

### Edge Cases

- Replicate naming variants: `_1`, `-2`, `rep3`, `r4`, `run5`, `replica_6`, `R7`, and nested
  layouts `LigA/rep1/`. All group under one system.
- One directory holds several trajectories (e.g. `step5_1.xtc … step5_10.xtc` parts): treated as
  one replicate with multiple parts. The tool proposes the processed (`*noPBC*`) trajectory when
  present.
- Topology formats: `.psf`, `.pdb`, `.gro`; `.tpr` is listed but not parsed for residue names.
- Very large topologies on remote hosts: only unique residue names are transferred.
- A residue name appearing in both APO and holo systems (e.g. a crystallographic ion) is not
  reported as a ligand. Ions, water and lipids are filtered.
- SSH host key unknown: trust on first use is reported. Authentication failure gives a clear
  message naming the key.
- No systems found: the tool reports which patterns were searched and how many directories were
  inspected.

## Requirements *(mandatory)*

### Functional Requirements

**Discovery (US1)**

- **FR-001**: The system MUST scan a root directory, locally or over SSH, to a configurable depth
  and find simulation directories, i.e. directories containing at least one trajectory file and
  one structure/topology file.
- **FR-002**: The system MUST group replicate directories into systems by stripping replicate
  tokens from directory names (`_N`, `-N`, `repN`, `rN`, `runN`, `replicaN`, `RN`, and a numeric
  last path component).
- **FR-003**: The system MUST flag a system as APO when its name contains an APO token (`apo`,
  `free`, `unbound`, `noligand`) or its topology contains no ligand residue.
- **FR-004**: The system MUST detect ligand residue names from each system's topology. It excludes
  protein, water, ion and lipid residues, and residue names shared with the APO system.
- **FR-005**: The system MUST record which analysis inputs exist per system: MD metric files,
  MM/PBSA result files, PCA folders and processed trajectories.
- **FR-006**: Remote discovery MUST use a bounded number of remote commands: one listing command,
  plus one residue-name command per unique topology.

**Naming and selection (US2)**

- **FR-007**: Users MUST be able to set a display name for each system and each ligand, mark
  APO/reference systems, and exclude systems or ligands.
- **FR-008**: Display names MUST be unique among included systems. Ligand display names MUST be
  unique within a system.
- **FR-009**: Generated selections MUST use original residue names. Labels and output names MUST
  use display names.

**Configuration output (US3)**

- **FR-010**: The system MUST generate a pipeline YAML for the selected stages that the existing
  loader accepts without edits, and validate it before writing.
- **FR-011**: The system MUST pre-select only stages whose inputs were found, and show the reason
  for each unselected stage.
- **FR-012**: The system MUST save and reload a project file containing the connection, root, scan
  results, names, exclusions and stage choices. Rescans MUST merge with saved choices.

**Interface (US4)**

- **FR-013**: The system MUST provide a keyboard-driven terminal interface with screens for:
  - connection;
  - systems table (editable);
  - ligands table (editable);
  - stage selection;
  - YAML preview;
  - save, run and submit actions.
- **FR-014**: The system MUST provide non-interactive commands with the same capabilities:
  - scan;
  - init (generate, with rename options);
  - run.

**Execution (US5)**

- **FR-015**: For local projects, the system MUST run the pipeline and report the run-report path.
- **FR-016**: For remote projects, the system MUST upload the configuration and start the pipeline
  on the host in the background, writing a log, or submit it through condorctl when available.
  It MUST report the log path or the job ID.

**Connection**

- **FR-017**: The system MUST connect over SSH using key-based authentication, either from a
  condorctl config/profile or from an SSH host alias. It MUST never prompt for or store passwords.

### Key Entities

- **Project**: connection (local / ssh host / condorctl profile), remote or local root, protein
  name, output root, systems, stage choices.
- **System**: discovered key, display name, included flag, APO flag, replicate directories,
  topology, trajectories per replicate, ligands, detected inputs.
- **Ligand**: residue name, display name, included flag, selection.
- **Stage choice**: stage name, enabled, reason when unavailable, key options.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: For a standard multi-ligand project (3 systems × 3 replicates), a user produces a
  valid pipeline configuration without editing any text file, in under 2 minutes.
- **SC-002**: The generated configuration passes the pipeline's dry-run validation in 100% of the
  supported layouts covered by the tests.
- **SC-003**: Remote discovery of a project with 50 replicate directories completes with no more
  than 1 + (number of unique topologies) remote commands.
- **SC-004**: Renaming a system or ligand once changes its name in every generated stage (0
  manual edits).
- **SC-005**: Every interface action (rename, toggle, stage selection, save) is covered by a
  headless interface test.

## Assumptions

- The interactive interface uses an optional dependency, installed with the `tui` extra, so the
  core package stays light. The non-interactive commands need no new dependency.
- SSH uses condorctl's connection when condorctl is installed and a profile is given. Otherwise a
  direct SSH connection reads host aliases from `~/.ssh/config`. Both use key authentication only.
- Remote hosts provide POSIX tools (`find`, `awk`, `sort`). The pipeline package must be installed
  on the remote host to run there. The tool checks this and reports it if missing.
- Generated configurations target the stages that can be fully configured from discovery: qc, md,
  prolif_batch, prolif_compare, ttclust_batch, pca_joint, distance_batch and ranking. MM/PBSA
  is included when result files are found.
- The project file is plain YAML and is human-readable but never needs hand-editing.
