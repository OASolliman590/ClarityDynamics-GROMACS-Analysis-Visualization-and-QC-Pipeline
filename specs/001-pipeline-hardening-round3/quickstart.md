# Quickstart / Validation: Pipeline Hardening Round 3

## Prerequisites

```bash
pip install -e . pytest MDAnalysis
```

## 1. Test suite (all user stories)

```bash
pytest -q
```

Expected: all tests pass. The ProLIF end-to-end checks are skipped unless ProLIF is installed.

## 2. Joint PCA (US1)

Add a `pca_joint` block (see [contracts/config-schema.md](./contracts/config-schema.md)), then:

```bash
CONFIG=my_pipeline.yaml ./run_pipeline.sh
```

Expected:
- `analysis_outputs/pca_joint/plots/joint_projection_12.html` shows every system on one set of
  axes, with APO in grey.
- `data/rmsip_pairwise.csv` has 1.0 on the diagonal.

## 3. Units (US3)

Point an MD metric at a GROMACS `.xvg` file in nm. In `statistics_summary.csv`, expect
`Value_Unit_Status = converted` and `Value_Unit_Factor = 10`.

## 4. Equilibration (US4)

Set `md.stability.enabled: true`. Expect `analysis_outputs/md/data/equilibration.json`.
Pass the reported `start_frame` to `mmpbsa.start_frame` and to the VMD task's `rmsf_first_frame`.

## 5. Ranking (US5)

Enable `ranking`. Expect these files:
- `data/ligand_ranking.csv` with `P_rank1`;
- `data/rank_pairwise.csv`;
- `data/excluded_systems.csv`, which lists APO.

## 6. Report (US11)

After any run, open `<output_root>/report.html`.

## 7. CI (US2)

Open a pull request. Expect a "tests" check running on Python 3.10 and 3.12.
