# CLI contract additions (feature 003)

```
clarity-setup init ... [--set STAGE.KEY=VALUE ...]
clarity-setup settings [--project P] [--stage S]          # list settings with current values
clarity-setup prepare  [--project P] [--what trjconv,metrics,mmpbsa] [--where local|host|condor]
                       [--force] [--dry-run]
clarity-setup run      ... [--where local|host|condor|dag]
clarity-setup status   [--project P] [--run ID] [--lines 20]
clarity-setup fetch    [--project P] [--run ID] [--to DIR]
clarity-setup doctor   [ROOT] [--ssh HOST | --condorctl CONFIG] [--project P]
gromacs-pipeline --config X [--no-cache]
```

Exit codes:

| Command | Codes |
|---|---|
| `status` | 0 |
| `fetch` | 0 on success; 1 when the run is not finished or the output is missing |
| `doctor` | 0 when everything needed is OK; 4 when a required item is missing |
| `prepare --dry-run` | 0, and lists the planned scripts |

TUI key additions:

| Key | Action |
|---|---|
| `p` | Prepare missing inputs |
| `g` | Settings for the selected stage (Stages tab) |
| `d` | Submit as a DAG |
| `j` | Jobs tab |
| `f` | Fetch the selected run |
| `D` | Doctor |
