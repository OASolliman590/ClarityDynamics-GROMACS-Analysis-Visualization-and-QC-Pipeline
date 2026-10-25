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
| `c` | Submit to HTCondor (remote connections only). |
| `ctrl+r` | Rescan. Your names and choices are kept, new folders are added, and vanished ones are marked missing. |
| `o` | Connect to a different location. |
| `q` | Quit. |

The **Stages** tab lists every stage with a checkbox. The **Preview** tab shows the
exact YAML that will run.

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
- **Run report**: every pipeline run writes `<output_root>/report.html`. It shows QC warnings
  first, then each stage with its status, figures and key tables.
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

## Publication PCA and MM/GBSA reports

The PCA, joint PCA, and MM/GBSA analyzers add a `publication/` directory by
default. Existing figures and filenames remain available. Publication figures
use a white theme and a consistent palette, with interactive HTML, SVG, PNG at
2× resolution, and PDF. Static exports require working Kaleido/Chrome (see
`plotly_get_chrome` or set `BROWSER_PATH`). Each format's actual success or failure
is recorded in JSON; a missing exporter does not prevent HTML or source tables.

Use the same `publication` block in a direct PCA/MMGBSA YAML file or beneath the
`pca`, `pca_batch` defaults, `pca_joint`, or `mmpbsa` stage in pipeline YAML:

```yaml
publication:
  enabled: true
  subdirectory: publication
  formats: [html, svg, png, pdf]
  scree_components: 10
  moving_window: 11
  top_n: 15
```

`moving_window` is a positive odd sample count. `top_n: 8` selects the eight most
negative residue TOTAL contributions; `top_n: all` retains every parsed residue.
The default selects 15, with stable sorting from most negative upward. The
all-residues figure grows with its input; HTML supports inspecting residue labels.
`enabled: false` suppresses only the additional publication outputs.

| Stage | Figures | Numerical sources and metadata |
|---|---|---|
| PCA / joint PCA | `pca_composite` with PC1–PC2, PC1–PC2–PC3, scree bars and cumulative variance | `pca_projections.csv`, `pca_variance.csv`, `pca_diagnostics.csv`, `pca_metadata.json` |
| PCA with incompatible bases | Independent `pca_state_N` composites; shared overlay recorded unavailable | Same source/diagnostic tables and reason |
| MM/GBSA | `components`, per-system `saved_total_*`, `residues` | Component, saved-sample, selected residue CSV sources and publication JSON metadata |

PCA fractions use the full eigenvalue denominator, even when only ten leading
components are shown. Native joint PCA records a fingerprint of the shared
basis and centering reference. Imported PCA files do **not** acquire shared-basis
status from matching filenames or equal numbers of atoms. Provide an explicit
`publication_inputs.json` in `base_dir` (`publication.pca_manifest` changes its
name) to attest externally prepared states:

```json
{
  "states": [
    {"name": "Apo", "projections": "apo_shared.csv", "eigenvalues": "eigenvals_shared.xvg",
     "basis_id": "joint-fit-001", "basis_verified": true, "basis_source": "joint-fit-log.json", "unit": "nm"},
    {"name": "Bound", "projections": "bound_shared.csv", "eigenvalues": "eigenvals_shared.xvg",
     "basis_id": "joint-fit-001", "basis_verified": true, "basis_source": "joint-fit-log.json", "unit": "nm"}
  ],
  "diagnostics": "publication_diagnostics.csv"
}
```

Each projections CSV has `PC1,PC2,PC3` columns in saved-sample order. PC3 may be
absent; the 3D panel then states its unavailability. Eigenvalues are two-column
XVG component/eigenvalue data containing the **full** spectrum, not a truncated
spectrum renormalized to 100%. The manifest paths resolve relative to the manifest.
Shared overlays require all states to attest the same verified basis identifier,
include basis provenance, and agree on variance fractions and units. The attestation
must describe the same eigenvectors, atom order, fitting reference and centering;
ClarityDynamics cannot independently reconstruct those facts from projection CSVs.
Single-state composites do not require a multi-state attestation. Legacy PCA inputs
are also supported, but their filenames alone never authorize a shared overlay.

Optional dedicated diagnostic CSV rows have columns
`metric,value,compatible,state_pair,window,source`. Supported metric keys are
`rmsip`, `window_sensitivity`, `split_half_rmsip`, and `cosine_content`. Values
must be finite fractions in [0,1], with `compatible: true` and a nonempty source.
Window sensitivity also requires a window label. For example:

```csv
metric,value,compatible,state_pair,window,source
rmsip,0.8,true,Apo|Bound,first 10 PCs,overlap-log.csv
window_sensitivity,0.7,true,Bound,0-10 ns,window-overlap.csv
split_half_rmsip,0.9,true,Bound,first half|second half,split-overlap.csv
```

Native joint RMSIP comes from compatible aligned subspaces. Parsed `cc_all.xvg`
cosine values are included when available. Every missing diagnostic remains in
the long summary with `status=unavailable`, a reason and no numeric value. The
schema includes metric, value, unit, status, reason, state pair, window and source.
PCA projections describe sampled coordinate covariance. They do not establish
equilibrium populations, a free-energy landscape, or convergence.

The MM/GBSA component composite distinguishes detailed VDWAALS/EEL/EGB/ESURF
terms, grouped GGAS/GSOLV sums, and the TOTAL endpoint estimate. They are
algebraically related quantities, not independent evidence. Panels retain independent
scales. A single results file uses saved-sample standard deviation (ddof=1);
independently aggregated replicate files use dispersion of replicate means.
Neither is a claim of experimental binding free energy.

Saved TOTAL traces show every estimate, a centered moving mean (partial windows
at the ends), and the overall mean. Replicates are kept separate when smoothing.
The default x-axis is **Saved sample index**. Frame numbers alone never supply time.
To use physical time, give a validated one-to-one frame/time map in the same
publication block; frame identities must match the retained results, and the mapped
times must be finite and strictly increasing in saved-sample order:

```yaml
publication:
  moving_window: 5
  top_n: all
  sampling_metadata:
    Bound:
      validated: true
      frames: [10, 20, 30]
      time: [0.1, 0.2, 0.3]
      unit: ns
      source: saved-frame-map.csv
```

Invalid or incomplete metadata leaves the sample-index axis in place and records
why physical time is unavailable. Negative residue values are favorable within
the endpoint decomposition. Optional missing decomposition, saved samples, or PCA
diagnostics are recorded without aborting independent outputs. Source CSV and JSON
files sit beside the figures and record the actual numerical values, selection,
smoothing window, availability, units, state labels, and input provenance. No MD,
docking, or energy calculations are performed by this reporting layer.
