# Scientific Validity Review

This document records a review of the analysis methods in ClarityDynamics. It lists
what was corrected in the code and the caveats that remain for users. Severity
describes the effect on conclusions: **Critical** findings could reverse a result,
**Major** ones bias values or uncertainties, and **Minor** ones affect labelling or
robustness.

## Fixed in code

| # | Severity | Module | Problem | Fix |
|---|----------|--------|---------|-----|
| 1 | Critical | `ranking` | The MM/PBSA term used for ranking matched any component containing `"G"`, so it picked `Binding_GGAS` (gas-phase MM energy) before `Binding_TOTAL`. Ligands were ranked on the wrong quantity. | Ranking now selects `Binding_TOTAL`/`ΔTOTAL` by exact name and skips the term with a warning when it is missing. |
| 2 | Critical | `md` | The z-score outlier filter was on by default (`|z| ≥ 3`) and also ran on RMSF. In an RMSF profile it deleted the most flexible residues (termini, loops) and shifted the residue axis. | Outlier filtering is off by default (`outlier_threshold: null`) and is never applied to per-residue profiles. |
| 3 | Major | `md` | In replicate mode, outliers were *deleted* from each replicate separately, then the replicates were aligned by index. Frames from different times were averaged together. | Outliers are masked as `NaN` so frames stay aligned. Averages and bands are NaN-aware. |
| 4 | Major | `md` | Replicate "Std" was the time-fluctuation of the averaged trace, so no between-replicate uncertainty was reported. Band SDs used `ddof=0`. | Added per-replicate means, `between_replicate_sd` (ddof=1) and `SEM = SD/√n`. Bands use the sample SD. |
| 5 | Major | `md` | Single-run statistics treated correlated MD frames as independent and gave no usable uncertainty. | Added a block-averaging SEM (Flyvbjerg–Petersen) with effective sample size and statistical inefficiency. |
| 6 | Major | `md` | Statistics were computed after smoothing, which reduced the SD and changed min/max. | Statistics use unsmoothed data. Smoothing affects plots only. |
| 7 | Major | `md` | Enabling the stability (equilibration) window stored integers in the per-system metric dict, which crashed the statistics and report steps. It also trimmed the RMSF *residue* axis by a *frame* index. | Window metadata is stored on each metric. RMSF is not trimmed (recompute it upstream over the stable window). Other metrics are trimmed by time on their own axes. |
| 8 | Major | `md`, `qc` | The time unit was guessed from magnitude (`max ≤ 1000 → ns`), so a 1 ns run saved in ps was read as 1000 ns. | The unit declared in the `.xvg` header (`xaxis label "Time (ps)"`) takes precedence. A warning is logged when only the heuristic is available. The RMSF residue axis is never rescaled as time. |
| 9 | Major | `mmpbsa` | The DAT parser's first regex captured only the average, so all SD/SEM values were silently 0. The replicate SD used `ddof=0`. A single run reported ΔTOTAL ± 0. | The parser reads the gmx_MMPBSA ≥1.5 columns (Average, SD(Prop.), SD, SEM(Prop.), SEM) and the older 3-column layout, plus PB terms. Replicate statistics use ddof=1 and report SEM. Every value is tagged `frame_sd` or `between_replicates`. |
| 10 | Major | `mmpbsa` | The report labelled ΔTOTAL as "Binding Affinity (ΔG_bind)". | It is now labelled as the effective binding energy without entropy. If IE/C2 entropy is available, ΔTOTAL + (−TΔS) is shown separately, with a caveat on interpreting end-point methods. |
| 11 | Major | `ranking` | The composite score included protein-wide descriptors with arbitrary "better" directions (lower SASA, Rg, RMSF, protein RMSD; more intra-protein H-bonds). Unknown metrics defaulted to "lower is better". | Only metrics with a defensible direction are scored by default (`rmsd_lig`, `rmsd_complex`, `comcom`, `distance_mean`, `mmpbsa_total`). Others are reported but not scored unless configured. Added `n_scored_metrics`. |
| 12 | Major | `ranking` | ANOVA and t-tests pooled per-residue distances as if residues were independent replicates of a ligand (pseudo-replication), with no multiple-comparison correction. | Replaced with a repeated-measures design (residues matched across ligands): a Friedman test plus pairwise Wilcoxon signed-rank tests with Holm correction. |
| 13 | Major | `distance` | Distances ignored periodic boundary conditions. | Minimum-atom distances pass the box to `distance_array`. COM distances use the minimum-image convention. |
| 14 | Minor | `distance` | The docs suggested `distance_scale: 10` for nm→Å, but MDAnalysis coordinates are already in Å (that setting gives a 10× error). | Docstring corrected. |
| 15 | Major | `qc` | QC did not check whether the ligand stays bound, and the clash check only logged (it never reached the report). | Records the min and max of the per-frame ligand–protein minimum distance. Warns on clashes (`min_distance_cutoff`) and on dissociation (`dissociation_cutoff`, default 5 Å), reporting the fraction of detached frames. |
| 16 | Minor | `qc` | No check for irregular frame spacing (restarts, duplicate or overlapping frames). | Added a sampling-uniformity check that writes a warning to the report. |
| 17 | Minor | `pca` | Comparison plots overlay projections, each computed on its own eigenvectors, as if they shared axes. | Added a log warning and an on-plot annotation. See the caveats below. |
| 18 | Minor | upstream VMD | Rg used a mass-weighted centre with an unweighted sum. | Rg is now fully mass-weighted, consistent with `gmx gyrate`. |
| 19 | Minor | `md` | The H-bond y-axis was labelled "Change in the number of bonds". | Relabelled "Number of hydrogen bonds". |

## Remaining methodological caveats

These are matters of experimental design. Code alone cannot resolve them.

1. **Replicates are the unit of inference.** Compare systems using the replicate SEM
   (`SEM` with `SEM_Method = between_replicates`). With one trajectory, the
   block-averaged SEM only reflects sampling within that run. It cannot detect a
   trajectory trapped in one basin. Use ≥3 independent replicates (different
   initial velocities) for any comparative claim.
2. **Equilibration.** Statistics include every frame unless the `stability` window is
   enabled or the inputs are already trimmed. Prefer to discard equilibration
   explicitly upstream (`gmx ... -b <t>`), and recompute RMSF and MM/PBSA over the
   same window as the other metrics.
3. **Units.** The VMD helper writes Å. Native GROMACS `.xvg` outputs use nm. Metric
   y-labels default to Å, so check the source of each file or set `ylabel`
   explicitly.
4. **PCA across systems.** Each system's PCs are separate bases with arbitrary
   sign. For quantitative comparison, fit all trajectories to one reference,
   concatenate them (`gmx trjcat -cat`), run a single `gmx covar`, and project every
   system onto those eigenvectors. Report subspace overlap (RMSIP,
   `gmx anaeig -over`) instead of comparing PC1 vs PC1. A cosine content near 1 for
   the first PCs indicates random-diffusion-like, unconverged sampling. Gibbs free
   energy landscapes from `gmx sham` are only meaningful for converged sampling.
5. **MM/PBSA.** ΔTOTAL is an end-point estimate. It is useful for ranking congeneric
   ligands but not as an absolute affinity. The SEM printed by gmx_MMPBSA assumes
   uncorrelated frames. Interaction-entropy results depend strongly on the
   trajectory length and on the σ of the interaction energy (unreliable when
   σ > ~3.6 kcal/mol), so present them separately.
6. **Composite ranking.** Averaging z-scores across heterogeneous metrics with only a
   handful of ligands is a heuristic for triage, not a statistical test. Inspect the
   individual metrics and their uncertainties before drawing conclusions.
7. **PBC preprocessing.** Centre-of-mass metrics require whole molecules. Process
   trajectories with `gmx trjconv -pbc mol -center` (with the ligand kept with the
   protein) before computing COM, distance or interaction metrics.
8. **Residue numbering.** The distance module keys residues by `resid`. In
   multi-chain systems with repeated residue numbers, use a `residue_selection`
   that includes the segment or chain.
9. **Legacy scripts.** `md_analysis_base.py`, `holo_md_analysis.py` and
   `apo_md_analysis.py` keep their original z-score outlier removal. Prefer the
   `gromacs_analysis` package.
