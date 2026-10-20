# Quickstart: Pipeline Setup TUI

```bash
pip install -e ".[tui]"          # textual + paramiko; add condorctl for HTCondor pools
```

## Local project
```bash
clarity-setup scan ./1-MD/analysis
clarity-setup init ./1-MD/analysis --protein PBP2 --rename HOLO_LigA=Isoniazid --ligand INH=Isoniazid
CONFIG=pipeline.yaml ./run_pipeline.sh --dry-run      # validates unchanged
```

## Remote (SSH alias or condorctl profile)
```bash
clarity-setup tui /data/project --ssh mycluster
clarity-setup tui ~/projects/pbp2 --condorctl ~/condorctl/config.yaml
```
In the TUI:
1. Review the **Systems** tab: edit names, toggle APO (`a`), exclude (`x`).
2. Review the **Ligands** tab: rename or exclude.
3. Pick analyses in **Stages**.
4. Check the YAML in **Preview**.
5. Press `s` to save, then `r` to run on the host or `c` to submit with condorctl.

## Validation
- `pytest tests/test_setup_*.py` covers scanning (local, plus the remote command path through the shell runner), naming, generation and validation, merging, and headless TUI interactions.
