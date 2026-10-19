# ClarityDynamics: GROMACS Analysis, Visualization, and QC Pipeline

Maintained package for molecular-dynamics analysis and visualization from GROMACS-derived outputs.

## Install

```bash
pip install -e .
pip install -e ".[all]"   # optional: trajectory, ProLIF, clustering extras
```

## Main Pipeline

Create a project config:

```bash
cp gromacs_pipeline_template.yaml my_pipeline.yaml
# or start from an empty anonymized guide:
cp pipeline_guide_empty.yaml my_pipeline.yaml
```

Validate:

```bash
CONFIG=my_pipeline.yaml ./run_pipeline.sh --dry-run
```

Run:

```bash
CONFIG=my_pipeline.yaml ./run_pipeline.sh
```

## Console Scripts

```bash
gromacs-md --help
gromacs-mmpbsa --help
gromacs-pca --help
gromacs-pipeline --help
gromacs-orchestrate --help
```

## Usage Guide

See [USAGE.md](./USAGE.md) for end-to-end use cases, stage combinations, expected outputs, and recommended run patterns.

## Scientific Validity

See [SCIENTIFIC_VALIDITY.md](./SCIENTIFIC_VALIDITY.md) for the statistical conventions (replicate SEM, block averaging), corrected defects, and methodological caveats to read before interpreting results.

## Orchestrated Upstream Workflow

Template:

```bash
cp orchestrated_workflow_template.yaml my_orchestrated_workflow.yaml
```

Dry run:

```bash
python -m gromacs_analysis.cli.orchestrate_cli --config my_orchestrated_workflow.yaml --dry-run
```

The orchestrator can build a PCA-specific index with RMSF-guided terminal trimming (`core_no_terminals_Calpha`) before PCA extraction.

## Package Layout

```text
gromacs_analysis/
├── cli/           # command-line interfaces
├── config/        # shared YAML helpers
├── md/            # RMSD/RMSF/Rg/SASA/H-bond/COM-COM analysis
├── mmpbsa/        # binding-energy parsing and plotting
├── pca/           # PCA parsing, plotting, comparison
├── distance/      # ligand distance analysis
├── interactions/  # basic interaction fingerprints
├── prolif/        # ProLIF fingerprints
├── ranking/       # cross-metric ranking
├── pipeline/      # multi-stage runner and YAML loader
└── utils/         # shared helpers
```

## Supported Stages

- `qc`
- `md`, `md_batch`, `md_compare`
- `mmpbsa`
- `pca`, `pca_batch`, `pca_compare`
- `distance`, `distance_batch`, `distance_compare`
- `interactions`, `interactions_batch`, `interactions_compare`
- `prolif`, `prolif_batch`, `prolif_compare`
- `ttclust`, `ttclust_batch`
- `ranking`
- `utils`

## Output Policy

Keep generated analysis outputs outside the source repository. Set `paths.output_root` in your pipeline YAML to a project-specific external folder.

## Path Rules

For stage input files, absolute paths stay absolute. Relative paths that already include the stage `base_dir` resolve from the YAML file location; shorter paths resolve under `base_dir`.

## Upstream Helpers

The `utils` stage can generate VMD, MM-PBSA, and PCA helper scripts with `type: generate_upstream_scripts`. These scripts are written for manual review and execution. The PCA stage can also run an opt-in GROMACS PCA generator when `run_gromacs.enabled: true`; it remains disabled by default.

## Legacy Scripts

Older standalone scripts are archived in `../../_legacy_scripts/`. They are retained for reference only, except for a few utility helpers used by the optional `utils` stage.

## Contributors

- Omar Ahmed Soliman (project author and maintainer)
- Codex (AI engineering collaborator by OpenAI)
