# Contract: `clarity-setup` command line

```
clarity-setup scan  ROOT [--ssh HOST | --condorctl CONFIG] [--depth N] [--json]
clarity-setup init  ROOT [--ssh HOST | --condorctl CONFIG] [--protein NAME]
                         [--rename OLD=NEW ...] [--ligand RES=NAME ...] [--apo NAME ...]
                         [--exclude NAME ...] [--stages a,b,c] [--project clarity_project.yaml]
                         [--out pipeline.yaml]
clarity-setup run   [--project clarity_project.yaml] [--out pipeline.yaml]
                         [--where local|host|condor] [--cpus N] [--memory 8GB]
clarity-setup tui   [ROOT] [--ssh HOST | --condorctl CONFIG] [--project clarity_project.yaml]
```

- `scan`: prints a table (system, replicates, APO, ligands, inputs, topology), or JSON with `--json`. Exit 0 if ≥1 system, 2 if none.
- `init`: scans (or loads `--project` and rescans), applies the edits, and writes the project file and the pipeline YAML. It validates with the pipeline loader and exits 1 with a message on validation errors.
- `run`: loads the project, regenerates and validates, then runs according to `--where` (default: `local` for local projects, `host` for remote). Prints the report path, the log path, or the condor cluster ID.
- `tui`: opens the interface. It requires the `tui` extra; otherwise it prints the install hint and exits 3.

Rename keys accept a system key or a current display name; `--ligand` keys are residue names (applied to every system) or `SYSTEM:RES`.
