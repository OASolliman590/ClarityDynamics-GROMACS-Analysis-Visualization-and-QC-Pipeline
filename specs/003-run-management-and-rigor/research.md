# Research: feature 003

**R1. Job state.** `condor_q -af ClusterId ProcId JobStatus` gives the state; `condor_history`
gives the state of jobs that have left the queue. The JobStatus codes are:

| Code | State |
|---|---|
| 1 | Idle |
| 2 | Running |
| 3 | Removed |
| 4 | Completed |
| 5 | Held |
| 6 | Transferring output |
| 7 | Suspended |

For a DAG, the jobs are the DAGMan cluster plus node clusters with `DAGManJobId == <dag
cluster>`. We query `condor_q -constraint 'DAGManJobId == N || ClusterId == N'`. Host runs use
`kill -0 PID` and `tail -n`.

**R2. Fetching results.** `tar czf` of the output root is written to a temporary file on the
remote, downloaded once (SFTP `get` or a local copy), then extracted. This was chosen over
per-file SFTP, which costs one round trip per file.

**R3. Index groups without numbers.** `gmx select -s tpr -on idx.ndx -select '"Receptor" group
"Protein"; "Ligand" resname X'` writes exactly two groups, numbered 0 and 1. That makes
`gmx_MMPBSA -cg 0 1` safe. trjconv accepts group names on stdin (`echo Protein System`).

**R4. Unit-to-full equivalence.**
- Batch stage outputs are `output_dir/output_{name}`, where `name` is the system name, or
  `system_ligand` for multi-ligand ProLIF.
- Emitting one ProLIF entry per ligand, named `system_ligand` only when the system has more than
  one ligand, reproduces exactly the folders the current loader creates.
- Unit configurations keep the full stage block minus the other entries, so every per-entry
  dataclass is identical. Therefore the cache fingerprints are identical, and the final job's
  cache hits.

**R5. Cache fingerprint.** We serialise the dataclass with `dataclasses.asdict`. Paths and other
non-JSON values are converted with `str`, sets are sorted, and keys are sorted. The resulting
JSON is hashed together with the input signatures (path, size and mtime). Plot settings are
included because figures are outputs.

**R6. Convergence estimators.**
- **Block SEM**: Flyvbjerg-Petersen, with the plateau taken as the maximum over levels that
  keep at least 4 blocks. This is the same as the existing MD processor, which we refactor to
  share it.
- **N_eff**: N/g, with g taken from `utils.equilibration.statistical_inefficiency`.
- **Half agreement**: z = |m1 − m2| / √(s1² + s2²) with block SEMs. Flagged when z > 3.
- **Cosine content** (Hess 2002): c_i = (2/T)(Σ cos(iπ(t+0.5)/T) p_i(t))² / Σ p_i(t)².
  Flagged when PC1 > 0.5.
- **Half-split RMSIP**: flagged below 0.5 (k = min(10, …)).

The thresholds are configurable.

**R7. Tests.**
- **Welch**: scipy.stats.ttest_ind(equal_var=False).
- **Permutation**: exact enumeration when C(nA+nB, nA) ≤ 20 000, otherwise 10 000 Monte-Carlo
  draws with a seed. The statistic is the difference of means, two-sided. With 3 vs 3
  replicates the smallest possible p is 0.1, which is stated in the `Note` column so readers
  are not misled.
- **Hedges' g**: pooled SD with the small-sample correction J = 1 − 3/(4(nA+nB)−9).
- **Adjustment**: Holm and BH within each metric, across all pairs.

**R8. Packaging.** PEP 621 `pyproject.toml` with the setuptools backend. `setup.py` is kept as
`setup()` for old pips. Debian's patched setuptools broke `setup.py develop`; PEP 660 editable
installs avoid that path.

**R9. Provenance.** `importlib.metadata` gives versions. `git rev-parse` runs only if the
package is in a git checkout, and so does `git status --porcelain` for the dirty flag. The
configuration is shown as YAML text from the loaded file path, when it is known.
