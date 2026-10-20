# Research: Pipeline Hardening Round 3

## R1. Joint PCA basis and fitting

- **Decision**: Fit all frames of all systems to a common reference with Kabsch superposition
  on the Cα. The first pass fits to frame 0 of the first system; the second pass fits to the mean
  structure of that first pass. Then compute the 3N × 3N covariance of the concatenated
  coordinates around their joint mean, eigendecompose it (`numpy.linalg.eigh`), sort the
  eigenvalues in descending order, and project each system's frames onto the top-k eigenvectors.
  Eigenvector signs are fixed so that the largest-magnitude component is positive, which makes
  the output deterministic.
- **Rationale**: This is the standard "essential dynamics in a combined subspace" approach
  (Amadei 1993; van Aalten 1997). The iterative fit removes any bias toward the arbitrary first
  frame. Equal atom weights suit Cα-only selections.
- **Alternatives considered**:
  - Generate a `gmx trjcat` + `gmx covar` script. Rejected as the default because it needs
    GROMACS and is hard to test; the existing GROMACS runner remains available for per-system PCA.
  - Mass weighting. Not useful for Cα, where all masses are equal.

## R2. Subspace overlap

- **Decision**: RMSIP over the first k components (default 10):
  RMSIP = sqrt((1/k) Σᵢ Σⱼ (vᵢ·wⱼ)²).
  It is computed between the per-system PCA subspaces (each obtained from the same commonly
  fitted coordinates) and between each system and the joint subspace.
- **Rationale**: This is the standard measure (Amadei 1999). It lies between 0 and 1 and is
  invariant to eigenvector sign.
- **Alternatives considered**: Covariance overlap (Hess 2002). It also works, but it requires the
  eigenvalues and is less familiar to readers; it can be added later.

## R3. Unit detection

- **Decision**: Parse `@ yaxis label "…"` and any `@ s0 legend` units. Recognise
  `nm`, `nm^2`, `nm\S2\N`, `nm²`, `Å`/`A`/`Angstrom`, `Å²`/`A^2`. Infer the display unit from the
  metric's `ylabel` (an `Å²` axis label means area, an `Å` axis label means length). The factors
  are nm→Å ×10 and nm²→Å² ×100. Anything else is not converted.
- **Rationale**: GROMACS writes these headers consistently. Inferring the target from the axis
  label keeps the plot axis and the data in agreement.
- **Alternatives considered**: A global "input unit" switch. Rejected because mixed sources are
  common.

## R4. Equilibration detection

- **Decision**: Implement Chodera's method (J. Chem. Theory Comput. 2016, 12, 1799):
  1. For candidate start indices t0 (every n/`n_candidates` points, at most 50% of the series),
     compute the statistical inefficiency g(t0) of the series from t0 onward. Use the
     autocorrelation-function integral with the pymbar-style cutoff: sum until the autocorrelation
     first becomes ≤ 0, and use mintime 3.
  2. Compute Neff(t0) = (n − t0) / g(t0).
  3. Choose t0 with the maximum Neff.
- **Rationale**: This is the published, widely used method; it needs no unit-dependent
  thresholds and reports Neff. An in-house implementation avoids adding pymbar as a dependency.
- **Alternatives considered**: Calling `pymbar.timeseries.detect_equilibration` directly. That
  adds a heavy optional dependency. The implementation will use pymbar when it is installed and
  configured, but the default is the in-house version.

## R5. Bootstrap ranking

- **Decision**: For B resamples (default 2000, seed 0):
  1. For each ligand and metric, draw a value. If replicate means are available, resample the
     mean of n replicate means with replacement. Otherwise draw from Normal(mean, SEM), using
     SEM = 0 when it is unknown.
  2. Recompute the z-scores and the composite score, then rank the ligands.
  3. Report P(rank = 1), the median rank, the 2.5–97.5% rank interval, and pairwise
     P(score_A < score_B), where a lower score is better.
- **Rationale**: This propagates replicate uncertainty through the same scoring rule used for
  the point ranking.
- **Alternatives considered**: Pairwise t-tests per metric. Rejected because they don't map onto
  the composite score.

## R6. Reference exclusion

- **Decision**: The set of excluded systems is the union of:
  - the `ranking.exclude_systems` list;
  - the names of APO systems in the `md` / `md_compare` stage configs, when those stages exist;
  - any system whose name is exactly "APO" (case-insensitive).
- **Rationale**: This avoids substring heuristics, so a ligand named "Apomorphine" is never
  excluded.

## R7. Protein–ligand H-bonds in VMD

- **Decision**: Run the VMD `hbonds` plugin with `-sel1 protein -sel2 ligand
  -writefile yes -outfile hbonds_lig_.dat`, and per ligand `hbonds_lig_<name>.dat`.
  Add a new MD metric `hbonds_lig` (holo only, pattern `hbonds_lig*.dat`) and include it in
  `LIGAND_METRICS` so it is split per ligand.
- **Rationale**: This is a direct counterpart to the existing intra-protein call.

## R8. Flat-profile threshold

- **Decision**: In `robust_threshold`, when MAD = 0, set scale = 0.7413 × IQR. If that is also 0,
  set scale = max(0.05 × |median|, 1e-6). Keep 1.4826 × MAD otherwise, so non-degenerate profiles
  are unchanged.
- **Rationale**: The IQR is robust to up to 25% outliers on each side; terminal overhangs are
  capped at 25% per end by `max_terminal_fraction`.

## R9. Chain-aware residues

- **Decision**: Key residues by `(chain, resid)`, where the chain is the `segid`, falling back to
  `chainID` and then "". Accept residue identifiers that are ints or "A:45" strings. An
  unqualified int matching several chains gives all matches plus one warning. Labels are
  `RESNAME45` when there is a single chain and `A:RESNAME45` when there are several.

## R10. Static export availability

- **Decision**: On the first static write, attempt it. If Kaleido or Chrome fails with an error
  indicating the browser or Kaleido is missing, cache "unavailable", log one warning naming the
  remedies (`plotly_get_chrome`, `BROWSER_PATH`, or `pip install "kaleido<1"`), and raise
  `StaticExportUnavailable` on later calls without retrying. The plotters already route SVG
  failures to matplotlib fallbacks; they log this exception at debug level instead of error.

## R11. Run report

- **Decision**: After the runner finishes (success or `continue_on_error`), write
  `<output_root>/report.html`. It has a QC warnings banner; for each stage, its status, error
  and output dir; for each plots dir found under the stage outputs, links/embeds (`<img>` for
  png/svg, links for html); and the first 10 rows of key CSVs (statistics_summary, ranking,
  results). It is self-contained HTML with relative paths, a light/dark-safe inline style, and
  no external assets.

## R12. CI

- **Decision**: `.github/workflows/tests.yml` runs on push to main and on pull requests, with a
  matrix of Python 3.10 and 3.12. It installs `requirements.txt`, `pytest` and `MDAnalysis`, then
  runs `pytest -q`. ProLIF is not installed, so those tests skip.
- **Risk**: The pushing token may lack the `workflow` scope. Mitigation: if the push is
  rejected, keep the file under `docs/ci/tests.yml` and tell the maintainer.

## R13. Loader split

- **Decision**: Move the template YAML text constant(s) into `pipeline/templates.py`. Move the
  ProLIF, TTClust and joint-PCA builders and `_batch_systems` into `pipeline/stage_builders.py`.
  `yaml_loader.py` imports and re-exports them. Older builders (md, mmpbsa, pca, distance)
  stay where they are, which limits the risk of the move.
