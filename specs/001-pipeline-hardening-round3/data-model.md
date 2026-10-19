# Data Model: Pipeline Hardening Round 3

## JointPCAResult

| Field | Type | Notes |
|-------|------|-------|
| system_names | list[str] | Order of systems as configured |
| n_atoms | int | Must be identical for all systems (FR-005) |
| eigenvalues | float[k] | Descending; unit nm² or Å² (follows coordinates, Å for MDAnalysis) |
| explained_variance | float[k] | eigenvalue / trace |
| eigenvectors | float[3N, k] | Sign-normalised (largest-|component| positive) |
| projections | dict[str, float[n_frames_i, k]] | Per system, on shared eigenvectors |
| rmsip_pairwise | float[S, S] | Between per-system subspaces (first `n_overlap` PCs); diagonal = 1 |
| rmsip_to_joint | float[S] | Each system's subspace vs joint subspace |
| parameters | dict | selection, fit iterations, n_components, n_overlap, frame caps, stride |

Validation: S ≥ 1; every system has ≥ 2 frames; n_overlap ≤ min(3N, frames−1).

## UnitConversionRecord (per system × metric)

| Field | Type | Notes |
|-------|------|-------|
| declared_unit | str or null | From header; null if absent |
| display_unit | str or null | From metric ylabel |
| factor | float | 1.0 when no conversion |
| status | enum | `converted`, `matched`, `assumed`, `unrecognised` |

## EquilibrationRecord (per system)

| Field | Type | Notes |
|-------|------|-------|
| system | str | |
| metric | str | Metric used for detection |
| method | enum | `chodera` (default) or `window` |
| start_index | int | Index into the metric's series |
| start_time | float | In the metric's output time unit |
| time_unit | str | |
| n_effective | float | (n − t0)/g |
| statistical_inefficiency | float | g ≥ 1 |
| status | enum | `detected`, `insufficient_data` |

## RankDistribution

Per ligand: `P_rank1`, `median_rank`, `rank_ci_low`, `rank_ci_high`, `n_scored_metrics`,
`n_bootstrap`, `seed`.
Pairwise: `LigandA`, `LigandB`, `P_A_better` (lower composite score is better).

## ResidueKey (distance)

`(chain: str, resid: int)`; label `RESNAME<resid>` or `<chain>:RESNAME<resid>` when several chains exist.

## RunReport

Stages in run order: `name`, `status` (`ok`, `skipped`, `failed`), `error`, `output_dirs`,
`figures` (relative paths), `tables` (relative path + first rows), `warnings` (QC).
