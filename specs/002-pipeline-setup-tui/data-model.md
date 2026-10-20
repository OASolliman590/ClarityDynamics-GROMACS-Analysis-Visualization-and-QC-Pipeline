# Data Model: Pipeline Setup TUI

## Project
| Field | Type | Notes |
|-------|------|-------|
| version | int | 1 |
| connection | Connection | how to reach the data |
| root | str | project root on that host |
| protein_name | str | default "Protein" |
| output_root | str | relative to root, default `analysis_outputs` |
| scanned_at | str | ISO time of last scan |
| systems | list[SystemEntry] | discovered + user edits |
| stages | dict[str, StageChoice] | |

## Connection
`kind`: `local` \| `ssh` \| `condorctl`; `host` (ssh alias or hostname); `profile` (condorctl config path).

## SystemEntry
| Field | Type | Notes |
|-------|------|-------|
| key | str | relative path with replicate tokens removed (stable id) |
| name | str | display name; unique among included systems |
| included | bool | default true |
| apo | bool | detected, user-editable |
| replicates | list[Replicate] | sorted by replicate number |
| ligands | list[LigandEntry] | |
| inputs | dict[str, bool] | md_metrics, mmpbsa, pca_dir, processed_trajectory |
| missing | bool | true when absent from the latest scan |

## Replicate
`dir` (relative to root), `number` (int), `topology`, `structure` (gro/pdb if different), `trajectories` (list, chosen first).

## LigandEntry
`resname` (str, selection source, immutable), `name` (display; unique within system), `included` (bool), `selection` (default `resname <resname>`).

## StageChoice
`enabled` (bool), `available` (bool), `reason` (str when unavailable), `options` (dict).

## Validation rules
- Included systems need ≥1 replicate with topology and trajectory.
- Display names: non-empty, `[A-Za-z0-9_.+-]` plus spaces converted to `_` for paths, unique (FR-008).
- At most one APO system is used as the reference in md stages; more are allowed and all are excluded from ranking.
