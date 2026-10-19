# Contract: Output files

## pca_joint (`<output_dir>/`)

| Path | Content |
|------|---------|
| `data/eigenvalues.csv` | `PC, Eigenvalue, ExplainedVariance, Cumulative` |
| `data/projections_<system>.csv` | `Frame, PC1..PCk` |
| `data/rmsip_pairwise.csv` | square matrix, systems × systems |
| `data/rmsip_to_joint.csv` | `System, RMSIP` |
| `data/parameters.json` | selection, n_atoms, frames per system, fit iterations, k, n_overlap |
| `plots/joint_projection_12.(html,svg)` | all systems on shared PC1–PC2 |
| `plots/joint_density_12.(html,svg)` | per-system density contours on shared axes |
| `plots/rmsip_heatmap.(html,svg)` | pairwise overlap |
| `plots/joint_scree.(html,svg)` | eigenvalues + cumulative variance (stacked) |

## md

- `statistics_summary.csv`: new columns `Replicate_Means` (semicolon-separated), `Value_Unit_Declared`,
  `Value_Unit`, `Value_Unit_Factor`, `Value_Unit_Status`.
- `data/equilibration.json`: list of EquilibrationRecord (see data-model.md), written when stability is enabled.

## mmpbsa

- `<system>_results.csv` / statistics: unchanged columns; `QC_` columns gain `Frames_Used` when
  `start_frame` is set.

## ranking

- `data/ligand_ranking.csv`: adds `P_rank1, median_rank, rank_ci_low, rank_ci_high` when bootstrap is enabled.
- `data/rank_pairwise.csv`: `LigandA, LigandB, P_A_better`.
- `data/excluded_systems.csv`: `System, Reason`.

## pipeline

- `<output_root>/report.html`: run report (see data-model.md RunReport).
