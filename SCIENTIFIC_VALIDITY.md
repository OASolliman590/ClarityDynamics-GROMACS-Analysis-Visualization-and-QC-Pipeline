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
| 20 | Major | `pca` | Cross-system PCA overlays compared axes from separate per-system analyses. | New `pca_joint` stage: all systems are fitted to one reference, analysed together and projected on shared eigenvectors, with RMSIP overlap between systems (`specs/001-pipeline-hardening-round3`). |
| 21 | Major | `md` | Values in nm (native GROMACS `.xvg`) were labelled Å. | Units declared in `.xvg` headers are converted to the axis unit (nm→Å ×10, nm²→Å² ×100) and recorded in `statistics_summary.csv` (`Value_Unit_*`). |
| 22 | Major | `md`, `mmpbsa`, VMD | Equilibration used unit-dependent std/slope thresholds, and RMSF and MM/PBSA could not use the same window. | `stability.method: chodera` (default when enabled) picks the start with maximum effective sample size and writes `data/equilibration.json`. `mmpbsa.start_frame` and the VMD `rmsf_first_frame` apply the same window. |
| 23 | Major | `ranking` | Rankings were point estimates, and APO could be ranked as a ligand. | Bootstrap rank probabilities, rank intervals and pairwise P(A better than B). Reference systems are excluded (`excluded_systems.csv`). |
| 24 | Major | `distance` | Residues with the same number in different chains were merged. | Residues are keyed by (chain, number). `"A:45"` identifiers are accepted, and an ambiguous bare number warns and reports each chain. |
| 25 | Minor | `pca` | On a flat RMSF core (MAD = 0), the threshold fell back to the SD, which the termini themselves inflated, so termini were under-trimmed. | IQR-based fallback scale. |
| 26 | Minor | VMD, `md` | Only intra-protein H-bonds were counted. | Protein–ligand H-bonds (`hbonds_lig`, per ligand in multi-ligand complexes). |

## Remaining methodological caveats

These are matters of experimental design. Code alone cannot resolve them.

1. **Replicates are the unit of inference.** Compare systems using the replicate SEM
   (`SEM` with `SEM_Method = between_replicates`). With one trajectory, the
   block-averaged SEM only reflects sampling within that run. It cannot detect a
   trajectory trapped in one basin. Use ≥3 independent replicates (different
   initial velocities) for any comparative claim. Bootstrap rank probabilities are
   only as good as the replicates behind them.
2. **Equilibration.** Statistics include every frame unless `md.stability` is enabled
   or the inputs are already trimmed. When it is enabled, reuse the start written to
   `data/equilibration.json` for RMSF (`rmsf_first_frame`) and MM/PBSA
   (`start_frame`) so every quantity covers the same window.
3. **Units of files without headers.** Plain `.dat` files declare no unit and are
   taken as given (`Value_Unit_Status = assumed`). Make sure their source matches the
   axis label.
4. **PCA.** Use `pca_joint` for quantitative cross-system comparison. Per-system PCA
   (`pca`, `pca_compare`) remains useful for each system's own motions only. A cosine
   content near 1 for the first PCs indicates random-diffusion-like, unconverged
   sampling. Gibbs free energy landscapes from `gmx sham` are only meaningful for
   converged sampling.
5. **MM/PBSA.** ΔTOTAL is an end-point estimate. It is useful for ranking congeneric
   ligands but not as an absolute affinity. The SEM printed by gmx_MMPBSA assumes
   uncorrelated frames. Interaction-entropy results depend strongly on the
   trajectory length and on the σ of the interaction energy (unreliable when
   σ > ~3.6 kcal/mol), so present them separately.
6. **Composite ranking.** Averaging z-scores across heterogeneous metrics with only a
   handful of ligands is a heuristic for triage, not a statistical test. Inspect the
   individual metrics, their uncertainties and the pairwise probabilities before
   drawing conclusions.
7. **PBC preprocessing.** Centre-of-mass metrics require whole molecules. Process
   trajectories with `gmx trjconv -pbc mol -center` (with the ligand kept with the
   protein) before computing COM, distance or interaction metrics.
8. **Legacy scripts.** `md_analysis_base.py`, `holo_md_analysis.py` and
   `apo_md_analysis.py` are deprecated (they warn when run) and keep their original
   z-score outlier removal. Use the `gromacs_analysis` package.

## Publication reporting

PCA and MM/GBSA now write additional publication composites and CSV/JSON sources
by default under `publication/`; legacy plots remain available. HTML, SVG, PNG
and PDF are attempted separately and unavailable exports include a reason.
The scree bars and cumulative curve share a percentage scale, without dual axes.

PCA variance fractions use the full covariance/eigenvalue denominator. Imported
multi-state projections require an explicit attestation of one verified basis,
including provenance; filenames alone cannot verify common eigenvectors, atom
order, reference fitting or centering. Native joint PCA records a basis fingerprint.
Projections represent sampled coordinate covariance, not equilibrium populations
or a free-energy landscape. RMSIP requires compatible subspaces. Dedicated
window and split-half results and parsed cosine values are included only when
supported; unavailable diagnostics have no numeric value and carry a reason.

MM/GBSA uses the canonical VDWAALS/EEL/EGB/ESURF/GGAS/GSOLV/TOTAL vocabulary.
The first four are detailed terms, the next two grouped sums, and TOTAL the
endpoint estimate. These are not independent energetic observations. The
report uses the gmx_MMPBSA energy convention (kcal/mol) and records it in metadata.
Within-file error bars are saved-sample SD (ddof=1), not replicate uncertainty;
true replicate summaries use dispersion of replicate means. Saved TOTAL traces
use sample index unless explicit validated metadata maps retained frames to
finite, strictly increasing times. Smoothing is a centered moving mean with
partial edge windows and configurable positive odd width (default 11); replicate
boundaries are preserved. Missing samples retain their original positions and
visible gaps, with identities recorded in metadata. Sample SD is unavailable for
one observation; a zero SD is never invented for that case. Residue panels select the 15 most negative parsed
TOTAL values by default, support custom positive top-N and all-residues modes,
and identify negative values as favorable within endpoint decomposition.

`USAGE.md` documents the imported-input and validated-time contracts. Synthetic
parser, plotter, analyzer and configuration tests verify these behaviors without
production trajectories or molecular/energy calculations.
