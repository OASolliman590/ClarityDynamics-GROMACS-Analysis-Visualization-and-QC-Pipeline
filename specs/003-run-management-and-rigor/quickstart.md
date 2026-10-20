# Quickstart validation (feature 003)

```bash
clarity-setup doctor /data/inhA --condorctl nmrbox.yaml
clarity-setup prepare --dry-run
clarity-setup prepare --where condor
clarity-setup init --set prolif_batch.stride=5
clarity-setup run --where dag
clarity-setup status
clarity-setup fetch --to results/
```

Expected:
- `doctor` lists the available tools.
- `prepare` writes `clarity/prep/<system>/<rep>/*.sh` and submits one job per script.
- After the prep jobs finish, a rescan shows `md_metrics` and `mmpbsa` as present.
- `run --where dag` prints the DAG cluster.
- `status` shows the per-state job counts.
- `fetch` downloads `report.html`. The report shows Convergence flags, Provenance and pairwise
  tests in the ranking section.
