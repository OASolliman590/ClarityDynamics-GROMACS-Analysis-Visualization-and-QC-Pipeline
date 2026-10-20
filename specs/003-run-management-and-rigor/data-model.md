# Data model: feature 003

These are added to `clarity_project.yaml`; all are backwards compatible.

- `settings: {stage: {key: value}}` holds overrides only. Defaults live in the schema.
- `runs: [RunRecord]` contains:
  - `id`: short id, `<kind>-<YYYYmmdd-HHMMSS>`;
  - `kind`: `analysis` or `prepare`;
  - `where`: `local`, `host`, `condor` or `dag`;
  - `pid`: optional;
  - `clusters: [int]`;
  - `dag_cluster`: optional;
  - `log`: path;
  - `output_dir`: path;
  - `run_dir`: remote folder holding the scripts;
  - `started`: ISO time;
  - `jobs: {name: cluster|null}`;
  - `last_state`: as of the last status check;
  - `fetched_to`: optional.
- `Replicate.files: {tpr, top, raw_trajectory, processed_trajectory}` holds relative paths or
  null.

The following are added to the pipeline YAML:

- `cache: true` per stage (generated).
- `convergence: {neff_min: 20, half_z_max: 3.0}` under `md`.
- `convergence: {cosine_max: 0.5, half_rmsip_min: 0.5}` under `pca_joint`.
- `ranking.tests: {alpha: 0.05, n_permutations: 10000}`.

Output files:

| File | Columns |
|---|---|
| `md/data/convergence.csv` | System, Metric, Replicate, N, Mean, SEM_block, N_eff, Half_Diff, Half_Z, Flags |
| `pca_joint/data/convergence.csv` | System, Replicate, Frames, Cos_PC1, Cos_PC2, Cos_PC3, Half_RMSIP, Flags |
| `ranking/data/pairwise_tests.csv` | as FR-019 |
