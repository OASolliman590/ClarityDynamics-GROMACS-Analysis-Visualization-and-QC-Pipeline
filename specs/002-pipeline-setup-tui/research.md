# Research: Pipeline Setup TUI

## R1. TUI framework
- **Decision**: Textual, as an optional extra (`pip install -e ".[tui]"`).
- **Rationale**:
  - It works in plain terminals and over SSH.
  - It has editable `DataTable`s, tabs, and checkboxes.
  - Its headless test pilot (`App.run_test`) lets every interaction be tested in CI.
  - It is actively maintained and pure Python.
- **Alternatives**:
  - curses: too low-level, with no testing harness.
  - prompt_toolkit / questionary wizards: linear prompts, with no table editing, and they repeat
    the existing `--interactive` wizard problem.
  - A web UI: needs a browser and port forwarding on clusters.

## R2. Remote access
- **Decision**: A `Backend` protocol with `run(cmd) -> (status, stdout, stderr)`, `read_text`,
  `write_text` and `is_remote`. There are three implementations:
  - `LocalBackend`: subprocess calls plus the local filesystem.
  - `SSHBackend`: paramiko, honouring `~/.ssh/config` (HostName, User, Port, IdentityFile).
    Key-based only, with trust on first use, recording host keys in `~/.ssh/known_hosts`
    (the same policy as condorctl).
  - `CondorctlBackend`: wraps `condorctl.connection.RemoteConnection` built from a condorctl
    config or profile, so NMRbox users reuse their working setup.
- **Rationale**: The scanner only needs to run commands, so one interface serves all three. Tests
  use a `ShellRunner` backend that executes commands locally. This exercises the exact remote
  command strings.

## R3. Discovery strategy
- **Decision**: One command per scan:

  ```
  find <root> -maxdepth D -type f \( -name '*.xtc' -o -name '*.trr' -o -name '*.dcd' -o -name '*.nc'
  -o -name '*.tpr' -o -name '*.gro' -o -name '*.pdb' -o -name '*.psf' -o -name '*.top'
  -o -name 'FINAL_RESULTS_MMPBSA*' -o -name 'RMSD*' -o -name 'RMSF*' -o -name 'eigenval*.xvg' \)
  -printf '%p\t%s\n'
  ```

  It falls back to `-exec stat` / `ls`-free parsing when `-printf` is unavailable (BusyBox/macOS),
  detected by the exit status. Grouping happens locally.
- **Directory classification**: a directory is a simulation directory if it has at least one
  trajectory and one structure/topology. The trajectory choice prefers names containing `noPBC`
  or `nojump`, then `production`, then the largest file. All `*.part*`/numbered parts are listed.
- **Replicate grouping**: the system key is the directory path relative to the root, with
  replicate tokens removed:
  - in the last component: `(_|-|\.)?(rep|replica|run|r|R)?\d+$`;
  - or a whole last component matching `^(rep|replica|run|r)?\d+$`, in which case the key is the
    parent.

  Within a key, replicates are sorted by their numeric suffix.

## R4. Ligand detection
- **Decision**: Remote-efficient unique residue-name extraction per topology format:
  - `.gro`: `awk 'NR>2 && length($0)>=10 {print substr($0,6,5)}' | sort -u` (the last line is the
    box, which has fewer than 10 characters before the coordinates; it is filtered by an
    atom-line check);
  - `.pdb`: `awk '/^(ATOM|HETATM)/{print substr($0,18,4)}' | sort -u`;
  - `.psf`: the `!NATOM` section, field 4.

  The output is filtered with the `pipeline.discovery` sets (protein/solvent/ion/lipid). Residue
  names that also occur in an APO system's topology are dropped (e.g. crystallographic ions).
- **Topology preference for residue names**: `.psf` > `.pdb` > `.gro` (CHARMM-GUI layouts). `.tpr`
  is binary and is skipped.

## R5. Stage availability
- `md`: metric files (`RMSD*`, `RMSF*`, `rog*`, `SASA*`, `hbonds*`, `comcom*`) found.
- `qc`, `distance_batch`, `prolif_batch`, `ttclust_batch`, `pca_joint`: trajectory + topology
  (always available after discovery). `prolif_batch` and `distance_batch` need at least one
  ligand. `ttclust_batch` needs the `ttclust` executable, a local check only.
- `mmpbsa`: `FINAL_RESULTS_MMPBSA*` found.
- `ranking`: md or mmpbsa enabled.
- Pre-selection is the set of available stages, minus heavy stages (`ttclust_batch`), which start
  unticked.

## R6. Remote run
- **Decision**:
  - **Run on host**: upload `pipeline.yaml` to `<root>/clarity/pipeline.yaml`, then
    `cd <root> && nohup gromacs-pipeline --config clarity/pipeline.yaml > clarity/pipeline.log 2>&1 &`.
    First check `command -v gromacs-pipeline`.
  - **Submit (condorctl)**: build a small `run_clarity.sh` that runs the same command, and submit
    it with condorctl's `submit_job` or CLI (`condorctl submit run_clarity.sh --cpus N --memory M`).
    Report the cluster ID.
- **Rationale**: Data never moves, and condorctl already handles pool policy.

## R7. Project file
- **Decision**: `clarity_project.yaml`. It records the version, connection
  `{kind, host, profile, root}`, protein, output_root, scanned_at, systems (key, name, included,
  apo, replicates [dir, topology, trajectories], ligands [resname, name, included], inputs) and
  stages {name: enabled}.
- **Merge on rescan**: match by system key and ligand resname; keep the user fields (name,
  included, apo, ligands' name/included); refresh the discovered fields; add new entries; mark
  vanished ones `missing: true` (shown, not generated).
