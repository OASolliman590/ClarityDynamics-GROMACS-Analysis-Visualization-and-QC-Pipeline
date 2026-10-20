# ClarityDynamics Usage

This guide explains practical use cases for the maintained package pipeline, from quick checks to full multi-ligand reporting.

## 0) Setup Without YAML

`clarity-setup` builds the pipeline configuration from what is on disk, so you only
choose names and stages. Install the interface extra first:

```bash
pip install -e ".[tui]"     # textual (interface) + paramiko (SSH)
```

### Where the data can be

| Data location | Option | Notes |
|---|---|---|
| This computer | (default) | `clarity-setup tui /path/to/simulations` |
| Any SSH host | `--ssh HOST` | Host or `~/.ssh/config` alias; key or agent authentication only. Unknown host keys are trusted on first use and remembered. |
| HTCondor pool via condorctl | `--condorctl CONFIG` | Reuses the condorctl config (host, user, key, pool attributes such as `+Production`). |

Discovery runs one `find` plus one residue listing per topology, so a remote scan is
quick even on a slow link. Nothing is copied back.

### What discovery recognises

- **Systems and replicates**: folders such as `HOLO_LigA_1..3`, `LigB_rep1`, `APO-r2`,
  or `LigC/rep1` are grouped into one system with numbered replicates.
- **APO reference**: names containing `apo`, `free`, `unbound` or `noligand`. If no
  folder is named that way, a system without ligands is used. Toggle it with `a`.
- **Ligands**: residues in the topology (psf, pdb or gro) that are not protein,
  solvent, ions or lipids. Residues that also occur in the APO system (cofactors
  such as HEM) are not counted as ligands.
- **Files**: the topology (psf > tpr > pdb > gro); a processed trajectory
  (`*noPBC*`, `*nojump*`, `*center*`, `*fit*`) over the raw one; metric files
  (`RMSD*`, `RMSF*`...); `FINAL_RESULTS_MMPBSA*`; and PCA `eigenval*`.
- **Stages**: only the stages the inputs support are offered. For example, `mmpbsa`
  needs `FINAL_RESULTS_MMPBSA` files and `md` needs metric files. Each unavailable
  stage says why. `ttclust_batch` is off by default because it is slow.

### Interactive interface

```bash
clarity-setup tui /path/to/simulations
clarity-setup tui                         # reopen clarity_project.yaml, or ask where to connect
```

| Key | Action |
|---|---|
| `e` | Rename the selected system or ligand. Ligand names are labels only; selections keep `resname XXX`. |
| `a` | Mark or unmark the selected system as the APO reference (excluded from ranking). |
| `x` | Include or exclude the selected system or ligand. |
| `s` | Save `clarity_project.yaml` and a validated `pipeline.yaml`. |
| `r` | Run: locally, or in the background on the remote host (`nohup`, log in `<root>/clarity/pipeline.log`). |
| `c` | Submit to HTCondor as one job (remote connections only). |
| `d` | Submit as a parallel HTCondor workflow (see [Parallel runs](#parallel-runs-on-htcondor)). |
| `p` | Prepare missing inputs (see [Preparing inputs](#preparing-missing-inputs)). |
| `g` | Settings of the focused stage on the Stages tab (or pick a stage). |
| `j` | Jobs tab: refresh the selected run's state and log. |
| `f` | Fetch the selected run's results into `clarity_results/<run id>/`. |
| `D` | Doctor: check the tools available where the data lives. |
| `ctrl+r` | Rescan. Your names and choices are kept, new folders are added, and vanished ones are marked missing. |
| `o` | Connect to a different location. |
| `q` | Quit. |

The **Stages** tab lists every stage with a checkbox. The **Preview** tab shows the
exact YAML that will run. The **Jobs** tab lists every run and preparation you started,
with its state and the end of its log.

### Command line (scripts, notebooks, CI)

```bash
clarity-setup scan /data/inhA --json
clarity-setup init /data/inhA \
    --protein InhA \
    --rename HOLO_LigA=Isoniazid_NAD --ligand INH=Isoniazid --ligand LigC:LGC=Compound_C \
    --apo APO --exclude HOLO_LigB --stages qc,md,mmpbsa,prolif_batch,ranking
clarity-setup run                                   # local data: runs here
clarity-setup run --where host                      # remote data: start on that host
clarity-setup run --where condor --cpus 8 --memory 16GB
```

`init` without ROOT re-reads `clarity_project.yaml` and keeps earlier edits. Exit
codes:

- 0: success;
- 1: an invalid edit or an invalid configuration (the message says which);
- 2: `scan` found no systems;
- 3: Textual is not installed.

HTCondor submission writes three files to `<root>/clarity/`: `pipeline.yaml`,
`run_clarity.sh` and `run_clarity.sub`. It relies on the pool's shared filesystem
(`should_transfer_files = NO`) and needs `gromacs-pipeline` installed on the
execute nodes.

### Check the environment first

```bash
clarity-setup doctor /data/inhA --condorctl nmrbox.yaml --where dag --prepare trjconv,metrics,mmpbsa
```

`doctor` checks, in one round trip:

- the connection and write access to the project folder;
- Python (3.10 or newer) and `gromacs-pipeline`;
- `gmx` (or `gmx_mpi`), `vmd` and `gmx_MMPBSA`;
- `condor_submit`, `condor_q` and `condor_submit_dag`;
- free disk space;
- whether the project folder is on a shared filesystem (NFS, Lustre, GPFS...), which
  HTCondor execute nodes need.

Each failed check prints a remedy. A check is *required* only for what you plan to do
(`--where`, `--prepare`). The exit code is 4 when a required check fails, otherwise 0.

### Preparing missing inputs

The scan records, for each replicate, the `.tpr`, `.top`, raw trajectory and processed
trajectory. `prepare` writes one reviewed script per replicate under
`<root>/clarity/prep/`. It only covers inputs that are missing, unless you pass `--force`:

| Task | What it does | Needs |
|---|---|---|
| `trjconv` | `gmx trjconv -pbc mol -center` and then `-fit rot+trans` produce `<raw>_noPBC.xtc`. Groups are passed by name (`Protein`, `System`, `Backbone`). | `.tpr` and a raw trajectory |
| `metrics` | The VMD metrics helper writes `RMSD_*.dat`, `RMSF_.dat`, `rog_.dat`, `SASA_.dat`, `hbonds_*.dat` and `comcom_*.dat` next to the trajectory. Multi-ligand complexes get per-ligand files. | a psf/pdb/gro structure and a trajectory |
| `mmpbsa` | `gmx select` builds an index with exactly two groups (0 = Protein, 1 = ligand), so `-cg 0 1` is always right. A default `clarity_mmpbsa.in` (GB, igb=5, 0.15 M salt, interaction entropy) is written, then `gmx_MMPBSA` runs. The ligand is the system's first included ligand. | `.tpr`, `.top` and a trajectory; holo systems only |

```bash
clarity-setup prepare --dry-run                    # list tasks and show the scripts; nothing is written
clarity-setup prepare --where condor --cpus 4      # one HTCondor job per replicate script
clarity-setup prepare --what metrics --where host  # only the VMD metrics, in the background on the host
clarity-setup prepare --gmx gmx_mpi --mpi-ranks 8  # site-specific executables
```

A replicate that lacks a prerequisite is listed with the reason (for example "no .tpr
run input in the replicate folder"); it is never skipped silently. When the jobs
finish, rescan with `init` or `ctrl+r`. The new inputs appear, and the stages that need
them become available.

### Analysis settings

Every stage has a few documented settings: strides, cutoffs, the PCA selection and
component count, the TTClust method, bootstrap and permutation counts, the MM/PBSA start
frame and the convergence thresholds. Defaults match the values the pipeline uses when
nothing is set. Only your changes are stored, in `clarity_project.yaml`.

```bash
clarity-setup settings                                   # list all settings, * = changed
clarity-setup settings --stage pca_joint
clarity-setup settings --set prolif_batch.stride=5 --set pca_joint.selection="name CA and resid 1-250"
clarity-setup init --set ranking.n_bootstrap=5000        # change while regenerating
clarity-setup settings --set qc.max_frames=none          # optional settings accept "none" (no limit)
```

Values are type-checked and range-checked, and invalid values are rejected with the
reason. In the interface, press `g` on a stage.

### Parallel runs on HTCondor

```bash
clarity-setup run --where dag --cpus 8 --unit-cpus 2 --memory 16GB
```

This submits a DAGMan workflow from `<root>/clarity/dag/`:

- **Parallel jobs**:
  - one per ProLIF ligand run, distance entry and TTClust system;
  - one each for QC, MD metrics, MM/PBSA and joint PCA.
- **A final job**: runs the whole configuration. Generated stages use `cache: true` and
  batch stages are cached per unit, so the final job skips everything the parallel jobs
  produced. It then builds what needs all of them (ProLIF comparison, ranking, report).
- **Failed units**: a unit that failed has no cache entry, so the final job re-runs it and
  reports it.

The wall time is roughly that of the slowest unit plus the final job, instead of the sum of
all units. Each generated ProLIF entry is one ligand (`<system>_<ligand>`). The output folders
are the same as before.

### Following runs and getting results

Every run and preparation you start is recorded in `clarity_project.yaml`:

```bash
clarity-setup status                 # latest run: state, job counts, last log lines
clarity-setup status --all --lines 5
clarity-setup fetch --to results/    # one archive of the output folder, downloaded and unpacked
```

How `status` determines the state:

- **Host runs**: the exit code written by the wrapper, or the process if it is still
  running.
- **HTCondor jobs and DAGs**: `condor_q` while jobs are queued, then `condor_history`.
  States are idle, running, held, completed, failed (non-zero exit) and removed.

`fetch` refuses unfinished runs unless you pass `--force`. The report opens from
`results/<run id>/report.html`.

### Caching

Generated configurations set `cache: true` on every stage. A stage is skipped only when both
its inputs (path, size, modification time) and its settings (a fingerprint of the full stage
configuration) are unchanged. Batch stages are cached per system or ligand, so changing one
system re-runs only that system. Failed results are never cached. Force a full run with
`gromacs-pipeline --config pipeline.yaml --no-cache`. Hand-written YAML keeps caching off
unless it sets `cache: true`.

## 1) Setup

```bash
cd streamlined_pipeline/gromacs_analysis_toolkit
pip install -e .
pip install -e ".[all]"
```

Create your config from either template:

```bash
cp gromacs_pipeline_template.yaml my_pipeline.yaml
# or
cp pipeline_guide_empty.yaml my_pipeline.yaml
```

Validate and run:

```bash
CONFIG=my_pipeline.yaml ./run_pipeline.sh --dry-run
CONFIG=my_pipeline.yaml ./run_pipeline.sh
```

## 2) Core Use Cases

### Use Case A: Preflight Scientific QC Before Any Analysis

Enable `qc` first in `run`.
Recommended when ingesting new trajectories or when replicate quality is uncertain.

What it checks:
- frame counts and atom counts
- timestep/time-axis detectability
- box/PBC metadata
- protein/ligand selection sanity
- minimum ligand-protein distance
- optional index/group consistency

### Use Case B: Single-System MD Stability and Dynamics

Enable `md`.
Typical outputs: RMSD, RMSF, Rg, SASA, H-bonds, COM distance, summary stats, and plots.

Key controls:
- `time_unit_input`, `time_unit_output`, optional `time_step_ps`
- `smoothing_window`
- `outlier_threshold`
- RMSF labeling via `amino_acids.sequence_start` or explicit labels

### Use Case C: Binding Energy Analysis with Replicate Integrity

Enable `mmpbsa`.

Key controls:
- `incomplete_replicates_policy`: `ignore|warn|error`
- `require_decomp_replicates`
- `compare_binding_key`

Recommended pattern:
- run `qc` and `md` first
- inspect MMPBSA QC metadata before ranking

### Use Case D: PCA and FEL for Conformational Landscape

Enable `pca`.
Optional: `run_gromacs.enabled: true` for automated `covar/anaeig` preparation.

Use for:
- scree and PC projection
- cosine content over windows
- FEL maps (`fel.enabled: true`)
- clustering in PCA space (`clustering.enabled: true`)

### Use Case E: Residue Distance Mapping

Enable `distance` for one ligand or `distance_batch` for many ligands.
Use `distance_compare` for cross-ligand heatmaps and deltas.

### Use Case F: Classical Interaction Fingerprints

Enable `interactions` or `interactions_batch`.
Use for basic interaction occupancy profiling independent of ProLIF chemistry typing.

### Use Case G: ProLIF Contact Chemistry and Visualization

Enable `prolif` or `prolif_batch`, then `prolif_compare` for cross-ligand analysis.

Current visualization defaults support:
- ligand network with full occurrence display (`lignetwork_display_all: true`)
- count-aware rendering (`lignetwork_count_aware: true`)
- barcode controls (`n_frame_ticks`, `residues_tick_location`, `xlabel`, `figsize`, `dpi`)
- residue-level barcode QC filtering (`only_interacting_residues`, `min_interaction_frames`, optional `max_residues`)

Expected barcode outputs:
- static SVG barcode (height grows with the number of interaction rows)
- interactive HTML barcode
- `barcode_qc_repN.csv` with per-residue interaction-frame counts

Per-ligand figures (`outputs` includes `occupancy`):
- `occupancy_heatmap` — residue × interaction occupancy (mean ± SD over replicates in the hover), sized to the residue count
- `occupancy_top_residues` — the `occupancy_top_n` most contacted residues with SD error bars
- `lignetwork` — aggregated over the trajectory (`lignetwork_kind: aggregate`, interactions present in ≥ `lignetwork_threshold` of frames); use `frame` for the old first-frame view

Automated per-ligand runs (`prolif_batch`):

```yaml
prolif_batch:
  base_dir: "./1-MD/analysis"
  auto_discover:                     # HOLO_<ligand>_<rep>/ -> one system per ligand, reps as replicates
    pattern: "HOLO_{ligand}_{rep}"
    topology: "step3_input.psf"
    trajectory: "step5_production_noPBC.xtc"
  ligands: auto                      # several ligands in one complex -> one run per ligand residue
  # ligand_exclude: ["NAD"]
```

`ligands: auto` keeps every residue that is not protein, water, ion or lipid. Use an explicit list
(`["INH", {name: NAD, selection: "resname NAD and segid B"}]`) when you need custom selections.
Runs are named `<system>_<ligand>` and feed `prolif_compare` (`from_prolif_batch: true`) unchanged.

### Use Case H: TTClust Conformational Clustering

Enable `ttclust` or `ttclust_batch`.

Typical outputs:
- dendrogram
- distance matrix
- linear and histogram views
- representative cluster PDBs

Use `stride` and selection fields to control runtime and clustering scope.

`ttclust_batch` accepts the same `auto_discover` block as `prolif_batch`, so every ligand found
on disk is clustered without listing systems. Besides TTClust's own images, each run writes
`data/clusters_summary.csv`, `data/cluster_assignment.csv` and two figures drawn from
`clustering.log`: `ttclust_populations` (cluster share of frames) and `ttclust_timeline`
(cluster membership along the concatenated, strided frames).

### Figures, colours and multi-ligand systems

- Every system keeps one colour across all stages: APO/reference in neutral grey, ligands in a
  fixed, colour-blind-checked order (blue, orange, aqua, yellow, magenta, green, violet, red).
  Beyond eight ligands a dash pattern is added rather than reusing a colour alone.
- Static exports (SVG/PNG) keep each figure's own size and use `plot_config.scale` as resolution,
  so fonts stay readable and content-sized figures (heatmaps, barcodes) are not squeezed.
- Several ligand-bound systems: each gets its own line, legend entry and band in every MD plot;
  `rmsd_combined` overlays the APO backbone with every complex.
- Several ligands in one complex: set `md.ligands: ["INH", "NAD"]`. Ligand RMSD and COM-COM become
  `rmsd_lig_<name>` / `comcom_<name>` metrics; the VMD helper writes matching
  `RMSD_<name>.dat` / `comcom_<name>.dat` when its task lists the same `ligands`.
- RMSF: flexible termini (smoothed RMSF above median + k·MAD, the rule used for the PCA core) are
  drawn faded inside grey bands with two dashed boundary lines and an explanatory note; the
  detected region is saved to `data/rmsf_terminal_detection.json`. Configure via `md.rmsf_terminals`.
- PCA: ligand colour vs grey APO, light→dark within a system for trajectory progress, diamond
  markers for means, Viridis free-energy surfaces, PC time series on the real time axis.
- MM/PBSA entropy: interaction entropy and C2 entropy are read from `FINAL_RESULTS*.dat`
  (also when results come from CSV). `entropy_<system>` shows ΔTOTAL, −TΔS and
  ΔG = ΔTOTAL + (−TΔS); `entropy_compare` compares systems; `data/entropy_summary.csv` holds the
  numbers. Interaction entropy with σ(ΔE_int) > 3.6 kcal/mol is hatched as unreliable.

### Use Case I: Cross-Ligand Scoring and Prioritization

Enable `ranking` after upstream stages complete.
Recommended inputs:
- MD compare stats
- MMPBSA QC-aware summaries
- distance compare summaries

### Use Case J: Joint PCA, Units, Equilibration and Ranking Uncertainty

- **Joint PCA** (`pca_joint`): every system's Cα trajectory is fitted to one reference and
  analysed together, so all systems share the same PC axes. Use the same `auto_discover`
  block as `prolif_batch`, or list `systems` (mark APO with `is_apo: true`). The stage writes:
  - `joint_projection_12` (all systems on shared PC1–PC2);
  - `joint_density_12`;
  - `rmsip_heatmap` (subspace overlap, 1 = identical motions);
  - `joint_scree`;
  - per-system projection CSVs.

  It needs MDAnalysis but not GROMACS.
- **Units**: `.xvg` files that declare nm or nm² are converted to the axis unit (Å/Å²)
  automatically. The `Value_Unit_*` columns of `statistics_summary.csv` show what was done.
  Set `md.convert_units: false` to disable this.
- **Equilibration**: `md.stability.enabled: true` (method `chodera`) detects the production
  start per system and writes `data/equilibration.json`. Pass its start frame to
  `mmpbsa.start_frame` and to the VMD helper's `rmsf_first_frame`, so all metrics use the
  same window.
- **Ranking uncertainty**: `ranking` adds these outputs:
  - `P_rank1`, `median_rank` and the 95% rank interval for each ligand;
  - `rank_pairwise.csv` with P(A better than B), from `n_bootstrap` resamples of the
    replicate means (or mean ± SEM).

  APO and other reference systems are excluded and listed in `excluded_systems.csv`.
- **Protein–ligand H-bonds**: the VMD helper writes `hbonds_lig_.dat`, plus
  `hbonds_lig_<name>.dat` per listed ligand. The MD stage plots them as `hbonds_lig`.
- **Multi-chain distances**: give residues as `"A:45"`. A bare number found in several chains
  is reported once per chain, with a warning.
- **Run report**: every pipeline run writes `<output_root>/report.html`. It shows these
  sections:
  - QC warnings;
  - sampling-check flags;
  - each stage with its status, figures and key tables;
  - a **Provenance** section: package versions, Python, host, the toolkit's git commit and
    whether it has uncommitted changes, and the exact configuration that ran.
- **Sampling checks**: the MD stage writes `data/convergence.csv` per system, metric and
  replicate over the analysed window. Flags:
  - `low_neff`: the block-averaged effective sample count is below `convergence.neff_min`
    (default 20);
  - `drift`: the first and second halves differ by more than `convergence.half_z_max`
    (default 3) combined standard errors.

  Joint PCA writes `data/convergence.csv` with these per-replicate checks:
  - the Hess cosine content of PC1–3 (flag `pc1_diffusive` above 0.5, meaning motion that
    looks like random diffusion);
  - the RMSIP between the first and second half of each system's frames (flag
    `unstable_subspace` below 0.5).

  Flagged rows are listed at the top of the report with an explanation of each flag.
- **Pairwise ligand tests**: `ranking` writes `data/pairwise_tests.csv`. It covers every
  scored metric, including MM/PBSA ΔTOTAL when per-replicate values exist. Every pair of
  ligands is compared on replicate means, with:
  - Welch's t-test;
  - an exact permutation test (Monte-Carlo above 20 000 relabellings);
  - Hedges' g;
  - Holm and Benjamini–Hochberg adjusted p-values within each metric.

  With 3 vs 3 replicates the smallest attainable permutation p is 0.1. The `Note` column
  says so, so a non-significant result is not misread as "no difference". Ligands with fewer
  than two replicates are listed as "insufficient replicates".
- **Static images**: if Kaleido cannot find a browser, the pipeline logs one warning with the
  fix (`plotly_get_chrome` or `BROWSER_PATH`). It keeps writing interactive HTML and falls back
  to matplotlib for SVG where it can.

## 3) Recommended Stage Recipes

### Conservative full scientific workflow

```yaml
run: ["qc", "md_batch", "md_compare", "mmpbsa", "pca_batch", "pca_compare", "prolif_batch", "prolif_compare", "ttclust_batch", "ranking"]
```

### Fast screening workflow

```yaml
run: ["qc", "md_batch", "mmpbsa", "prolif_batch", "ranking"]
```

### Visualization-focused workflow

```yaml
run: ["prolif_batch", "prolif_compare", "ttclust_batch", "pca_batch", "pca_compare"]
```

## 4) Utility Workflows (`utils` stage)

Available helper task patterns include:
- `examine_csv`
- `convert_dat`
- `renumber_gro`
- `modevectors`
- `generate_upstream_scripts`
- `duivy_inspect_guide`

`duivy_inspect_guide` is intended for raw `.xvg/.xpm` inspection guidance and not a replacement for maintained toolkit outputs.

## 5) Validation Gates

Run these before sharing results:

```bash
python -m compileall -q gromacs_analysis
pytest -q
python -m gromacs_analysis.cli.pipeline_cli --generate-template
CONFIG=my_pipeline.yaml ./run_pipeline.sh --dry-run
```

## 6) Data and Provenance Practices

- Keep raw trajectories and bulky generated outputs outside source control.
- Set `paths.output_root` to project storage.
- Keep project YAMLs free of personal absolute paths when sharing publicly.
- Store system-specific private paths in local, non-committed configs.
