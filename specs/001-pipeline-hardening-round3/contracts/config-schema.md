# Contract: Configuration keys (pipeline YAML)

Only new or changed keys are listed. Existing keys keep their meaning.

## `pca_joint` (new stage)

```yaml
pca_joint:
  enabled: true
  base_dir: "./1-MD/analysis"
  output_dir: "./analysis_outputs/pca_joint"
  selection: "protein and name CA"     # same atom count required in every system
  n_components: 10                     # PCs written/projected
  n_overlap: 10                        # PCs used for RMSIP
  stride: 1
  max_frames_per_system: null          # optional cap to balance systems
  fit_iterations: 2
  systems:                             # or auto_discover (same block as prolif_batch)
    - name: "APO"
      is_apo: true
      topology: "APO_1/step3_input.psf"
      trajectories: ["APO_1/step5_production_noPBC.xtc", "APO_2/step5_production_noPBC.xtc"]
    - name: "LigandA"
      topology: "HOLO_LigandA_1/step3_input.psf"
      trajectories: ["HOLO_LigandA_1/step5_production_noPBC.xtc"]
```

## `md`

| Key | Change |
|-----|--------|
| `stability.method` | NEW: `chodera` (default) or `window` (previous behaviour) |
| `stability.n_candidates` | NEW: candidate start points for `chodera` (default 100) |
| `convert_units` | NEW: bool, default `true` (FR-007) |
| metric `unit` | NEW optional display unit; inferred from `ylabel` when omitted |

New default metric: `hbonds_lig` (pattern `hbonds_lig*.dat`, holo only, split per ligand).

## `mmpbsa`

| Key | Change |
|-----|--------|
| `start_frame` | NEW int, default null. With frame-level CSV results, frames with `Frame #` < start_frame are excluded |

## `ranking`

| Key | Change |
|-----|--------|
| `exclude_systems` | NEW list; merged with md APO systems and exact name "APO" |
| `n_bootstrap` | NEW int, default 2000 (0 disables) |
| `seed` | NEW int, default 0 |

## `distance`

`residues.ids` accepts ints or chain-qualified strings (`"A:45"`).

## Upstream VMD task (`utils` / orchestrator)

| Key | Change |
|-----|--------|
| `rmsf_first_frame` | NEW int, default 0 |
| `ligands` | existing; now also drives per-ligand `hbonds_lig_<name>.dat` |
