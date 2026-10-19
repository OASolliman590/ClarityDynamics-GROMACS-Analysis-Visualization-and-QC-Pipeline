# ClarityDynamics Usage

This guide explains practical use cases for the maintained package pipeline, from quick checks to full multi-ligand reporting.

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
