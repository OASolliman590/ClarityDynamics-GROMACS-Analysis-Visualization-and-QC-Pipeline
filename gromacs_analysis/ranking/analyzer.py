"""Ranking analyzer across ligands."""

from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

from .config import RankingConfig

logger = logging.getLogger(__name__)


class RankingAnalyzer:
    def __init__(self, config: RankingConfig):
        self.config = config

    def run_analysis(self) -> Dict:
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        data_dir = self.config.output_dir / "data"
        data_dir.mkdir(parents=True, exist_ok=True)

        metrics: Dict[str, Dict[str, float]] = {}
        uncertainty: Dict[str, Dict[str, Dict[str, Any]]] = {}
        qc_metadata: Dict[str, Dict[str, Any]] = {}

        # MD metrics
        if self.config.md_stats_file and self.config.md_stats_file.exists():
            md_df = pd.read_csv(self.config.md_stats_file)
            for _, row in md_df.iterrows():
                lig = row.get("System")
                metric = row.get("Metric")
                mean = row.get("Mean")
                if lig is None or metric is None:
                    continue
                metrics.setdefault(lig, {})[metric] = float(mean)
            for lig, values in self._uncertainty_from_md(md_df).items():
                uncertainty.setdefault(lig, {}).update(values)

        # MMPBSA total binding
        if self.config.mmpbsa_data_dir and self.config.mmpbsa_data_dir.exists():
            for file_path in self.config.mmpbsa_data_dir.glob("*_results.csv"):
                lig = file_path.stem.replace("_results", "")
                df = pd.read_csv(file_path)
                qc_from_results = self._extract_qc_from_results_df(df)
                total_val = self._select_total_binding(df)
                if total_val is None:
                    logger.warning(f"No ΔTOTAL / Binding_TOTAL row found in {file_path}; skipping MM/PBSA term.")
                else:
                    metrics.setdefault(lig, {})["mmpbsa_total"] = total_val
                    sem = self._select_total_sem(df)
                    uncertainty.setdefault(lig, {})["mmpbsa_total"] = {"mean": total_val, "sem": sem}
                    reps = self._select_total_replicates(df)
                    if len(reps) > 1:
                        uncertainty[lig]["mmpbsa_total"]["replicates"] = reps

                qc_sidecar = self._extract_qc_from_sidecar(
                    self.config.mmpbsa_data_dir / f"{lig}_qc.csv"
                )
                merged_qc = qc_sidecar or qc_from_results
                if qc_sidecar and qc_from_results:
                    merged_qc = {**qc_from_results, **qc_sidecar}
                if merged_qc:
                    qc_metadata[lig] = merged_qc

        # Distance summary (mean across residues)
        if self.config.distance_summary_file and self.config.distance_summary_file.exists():
            dist_df = pd.read_csv(self.config.distance_summary_file)
            if {"Ligand", "MeanDistance"}.issubset(dist_df.columns):
                for lig, group in dist_df.groupby("Ligand"):
                    metrics.setdefault(lig, {})["distance_mean"] = float(group["MeanDistance"].mean())

        if not metrics:
            raise ValueError("No metrics found for ranking")

        excluded = self._excluded_systems(metrics)
        pd.DataFrame(
            [{"System": name, "Reason": reason} for name, reason in excluded.items()],
            columns=["System", "Reason"],
        ).to_csv(data_dir / "excluded_systems.csv", index=False)
        if excluded:
            logger.info(f"Excluded from ranking (reference systems): {', '.join(excluded)}")
        metrics = {lig: vals for lig, vals in metrics.items() if lig not in excluded}
        if not metrics:
            raise ValueError("All systems are reference systems; nothing left to rank")

        ranking_df = self._build_ranking_table(metrics, qc_metadata)
        pairwise_path = None
        if self.config.n_bootstrap > 0 and len(metrics) > 1:
            rank_dist, pairwise = self._bootstrap_ranks(metrics, uncertainty)
            ranking_df = ranking_df.merge(rank_dist, on="Ligand", how="left")
            pairwise_path = data_dir / "rank_pairwise.csv"
            pairwise.to_csv(pairwise_path, index=False)
        ranking_path = data_dir / "ligand_ranking.csv"
        ranking_df.to_csv(ranking_path, index=False)

        stats_path = self._compute_stats(data_dir)
        tests_path = self._write_pairwise_tests(metrics, uncertainty, data_dir)

        return {
            "success": True,
            "output_dir": str(self.config.output_dir),
            "ranking_file": str(ranking_path),
            "stats_file": str(stats_path) if stats_path else None,
            "pairwise_file": str(pairwise_path) if pairwise_path else None,
            "tests_file": str(tests_path) if tests_path else None,
            "excluded_systems": sorted(excluded),
        }

    def _build_ranking_table(
        self,
        metrics: Dict[str, Dict[str, float]],
        qc_metadata: Optional[Dict[str, Dict[str, Any]]] = None,
    ) -> pd.DataFrame:
        ligand_names = sorted(metrics.keys())
        all_metrics = sorted({m for vals in metrics.values() for m in vals.keys()})

        data = {"Ligand": ligand_names}
        for metric in all_metrics:
            data[metric] = [metrics[lig].get(metric, np.nan) for lig in ligand_names]

        df = pd.DataFrame(data)

        df = self._add_scores(df, all_metrics)

        z_cols = [c for c in df.columns if c.endswith("_z")]
        if z_cols:
            df["n_scored_metrics"] = df[z_cols].notna().sum(axis=1)
        else:
            df["n_scored_metrics"] = 0

        return self._finish_ranking_table(df, qc_metadata)

    def _add_scores(self, df: pd.DataFrame, all_metrics: List[str], log: bool = True) -> pd.DataFrame:
        """Direction-aware z-scores per metric and their mean (lower rank_score = better)."""
        # Z-scores. Only metrics with an explicit, physically meaningful direction
        # contribute to the score; descriptive metrics (Rg, SASA, RMSF, protein RMSD,
        # intra-protein H-bonds) are reported but have no "better" direction.
        unscored = [m for m in all_metrics if m not in self.config.metric_directions]
        if unscored and log:
            logger.info(f"Metrics reported but not scored (no direction configured): {', '.join(unscored)}")
        for metric in all_metrics:
            if metric not in self.config.metric_directions:
                continue
            vals = df[metric].astype(float)
            if vals.isna().all():
                continue
            z = (vals - vals.mean()) / (vals.std(ddof=0) if vals.std(ddof=0) else 1.0)
            direction = self.config.metric_directions[metric]
            if direction == "higher":
                z = -z
            df[f"{metric}_z"] = z

        z_cols = [c for c in df.columns if c.endswith("_z")]
        df["rank_score"] = df[z_cols].mean(axis=1) if z_cols else np.nan
        return df

    def _finish_ranking_table(self, df: pd.DataFrame, qc_metadata) -> pd.DataFrame:
        qc_metadata = qc_metadata or {}
        if qc_metadata:
            df["mmpbsa_expected_replicates"] = df["Ligand"].map(
                lambda lig: self._as_float(qc_metadata.get(lig, {}).get("expected_replicates"))
            )
            df["mmpbsa_found_results_replicates"] = df["Ligand"].map(
                lambda lig: self._as_float(qc_metadata.get(lig, {}).get("found_results_replicates"))
            )
            df["mmpbsa_missing_results_replicates"] = df["Ligand"].map(
                lambda lig: self._as_float(qc_metadata.get(lig, {}).get("missing_results_replicates"))
            )
            df["mmpbsa_replicate_complete"] = df["Ligand"].map(
                lambda lig: self._to_bool(qc_metadata.get(lig, {}).get("replicate_complete"))
            )
            df["mmpbsa_group_identity_ok"] = df["Ligand"].map(
                lambda lig: self._to_bool(qc_metadata.get(lig, {}).get("group_identity_ok"))
            )
            df["qc_penalty"] = df.apply(self._qc_penalty, axis=1)
            df["rank_score_qc"] = df["rank_score"].fillna(0.0) + df["qc_penalty"]
            df["qc_status"] = df["qc_penalty"].apply(lambda penalty: "warn" if penalty > 0 else "pass")
            df = df.sort_values("rank_score_qc")
        else:
            df = df.sort_values("rank_score", na_position="last")

        return df

    # ------------------------------------------------------------ uncertainty

    @staticmethod
    def _uncertainty_from_md(md_df: pd.DataFrame) -> Dict[str, Dict[str, Dict[str, Any]]]:
        """Per ligand/metric: replicate means when exported, else mean and SEM."""
        out: Dict[str, Dict[str, Dict[str, Any]]] = {}
        for _, row in md_df.iterrows():
            lig, metric = row.get("System"), row.get("Metric")
            if lig is None or metric is None or pd.isna(row.get("Mean")):
                continue
            entry: Dict[str, Any] = {"mean": float(row["Mean"])}
            reps = row.get("Replicate_Means")
            if isinstance(reps, str) and reps.strip():
                values = [float(v) for v in reps.split(";") if v.strip()]
                if len(values) > 1:
                    entry["replicates"] = values
            sem = row.get("SEM")
            entry["sem"] = float(sem) if sem is not None and not pd.isna(sem) else 0.0
            out.setdefault(lig, {})[metric] = entry
        return out

    def _write_pairwise_tests(self, metrics, uncertainty, data_dir: Path) -> Optional[Path]:
        """Replicate-level pairwise tests for every metric shared by two or more ligands."""
        from .stats import pairwise_tests  # noqa: WPS433

        values: Dict[str, Dict[str, List[float]]] = {}
        for lig, lig_metrics in metrics.items():
            for metric in lig_metrics:
                reps = uncertainty.get(lig, {}).get(metric, {}).get("replicates") or []
                values.setdefault(metric, {})[lig] = list(reps)
        values = {m: groups for m, groups in values.items() if len(groups) > 1}
        if not values:
            return None
        table = pairwise_tests(values, directions=self.config.metric_directions,
                               n_permutations=self.config.n_permutations, seed=self.config.seed)
        path = data_dir / "pairwise_tests.csv"
        table.to_csv(path, index=False)
        return path

    @staticmethod
    def _select_total_replicates(df: pd.DataFrame) -> List[float]:
        if "Replicate_Values" not in df.columns or "Component" not in df.columns:
            return []
        components = df["Component"].astype(str).str.strip()
        for name in ("Binding_TOTAL", "ΔTOTAL", "ΔGTOTAL", "DELTA TOTAL", "TOTAL"):
            match = df[components == name]
            if not match.empty:
                raw = match.iloc[0]["Replicate_Values"]
                if isinstance(raw, str) and raw.strip():
                    return [float(v) for v in raw.split(";") if v.strip()]
                return []
        return []

    @staticmethod
    def _select_total_sem(df: pd.DataFrame) -> float:
        if "SEM" not in df.columns or "Component" not in df.columns:
            return 0.0
        components = df["Component"].astype(str).str.strip()
        for name in ("Binding_TOTAL", "ΔTOTAL", "ΔGTOTAL", "DELTA TOTAL", "TOTAL"):
            match = df[components == name]
            if not match.empty:
                value = match.iloc[0]["SEM"]
                return float(value) if not pd.isna(value) else 0.0
        return 0.0

    def _excluded_systems(self, metrics: Dict[str, Dict[str, float]]) -> Dict[str, str]:
        configured = {str(name) for name in self.config.exclude_systems}
        excluded: Dict[str, str] = {}
        for name in metrics:
            if name in configured:
                excluded[name] = "configured reference system"
            elif str(name).strip().upper() == "APO":
                excluded[name] = "APO reference"
        return excluded

    def _bootstrap_ranks(self, metrics, uncertainty):
        """Rank distribution by resampling each ligand's metrics from their uncertainty.

        Uses the same direction-aware z-score composite as the point ranking, vectorised
        over resamples: samples has shape (bootstrap, ligand, scored metric).
        """
        rng = np.random.default_rng(self.config.seed)
        ligands = sorted(metrics)
        scored_metrics = [m for m in sorted({m for v in metrics.values() for m in v})
                          if m in self.config.metric_directions]
        n_boot = int(self.config.n_bootstrap)
        samples = np.full((n_boot, len(ligands), max(1, len(scored_metrics))), np.nan)
        for j, lig in enumerate(ligands):
            for k, metric in enumerate(scored_metrics):
                if metric not in metrics[lig]:
                    continue
                info = uncertainty.get(lig, {}).get(metric, {})
                if info.get("replicates"):
                    reps = np.asarray(info["replicates"], dtype=float)
                    draws = rng.choice(reps, size=(n_boot, len(reps)), replace=True).mean(axis=1)
                else:
                    sem = float(info.get("sem") or 0.0)
                    draws = rng.normal(metrics[lig][metric], sem, n_boot) if sem > 0 else metrics[lig][metric]
                samples[:, j, k] = draws

        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.nanmean(samples, axis=1, keepdims=True)
            std = np.nanstd(samples, axis=1, keepdims=True)
            std = np.where((std == 0) | ~np.isfinite(std), 1.0, std)
            z = (samples - mean) / std
            for k, metric in enumerate(scored_metrics):
                if self.config.metric_directions[metric] == "higher":
                    z[:, :, k] = -z[:, :, k]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)  # ligands with no scored metric
                scores = np.nanmean(z, axis=2)
        ranks = np.vstack([
            pd.Series(row).rank(method="average", na_option="bottom").to_numpy() for row in scores
        ])

        # Ties for first place share the probability equally.
        is_best = ranks == ranks.min(axis=1, keepdims=True)
        p_rank1 = (is_best / is_best.sum(axis=1, keepdims=True)).mean(axis=0)
        dist = pd.DataFrame({
            "Ligand": ligands,
            "P_rank1": p_rank1,
            "median_rank": np.median(ranks, axis=0),
            "rank_ci_low": np.percentile(ranks, 2.5, axis=0),
            "rank_ci_high": np.percentile(ranks, 97.5, axis=0),
            "n_bootstrap": n_boot,
            "seed": self.config.seed,
        })
        pairs = []
        for i, a in enumerate(ligands):
            for j, bname in enumerate(ligands):
                if i == j:
                    continue
                better = np.mean(scores[:, i] < scores[:, j]) + 0.5 * np.mean(scores[:, i] == scores[:, j])
                pairs.append({"LigandA": a, "LigandB": bname, "P_A_better": float(better)})
        return dist, pd.DataFrame(pairs)

    @staticmethod
    def _select_total_binding(df: pd.DataFrame) -> Optional[float]:
        """Return the total effective binding energy (ΔTOTAL), never a partial term like ΔGGAS."""
        if "Component" not in df.columns or "Mean" not in df.columns:
            return None
        components = df["Component"].astype(str).str.strip()
        for name in ("Binding_TOTAL", "ΔTOTAL", "ΔGTOTAL", "DELTA TOTAL", "TOTAL"):
            match = df[components == name]
            if not match.empty:
                return float(match.iloc[0]["Mean"])
        return None

    @staticmethod
    def _as_float(value: Any) -> Optional[float]:
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _to_bool(value: Any) -> Optional[bool]:
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return None
        if isinstance(value, bool):
            return value
        if isinstance(value, (int, np.integer)):
            return bool(value)
        text = str(value).strip().lower()
        if text in {"true", "1", "yes", "y"}:
            return True
        if text in {"false", "0", "no", "n"}:
            return False
        return None

    def _qc_penalty(self, row: pd.Series) -> float:
        penalty = 0.0

        missing = row.get("mmpbsa_missing_results_replicates")
        expected = row.get("mmpbsa_expected_replicates")
        if pd.notna(missing) and pd.notna(expected) and float(expected) > 0:
            penalty += float(missing) / float(expected)

        replicate_complete = row.get("mmpbsa_replicate_complete")
        if replicate_complete is False:
            penalty += 0.5

        group_identity_ok = row.get("mmpbsa_group_identity_ok")
        if group_identity_ok is False:
            penalty += 0.5

        return float(penalty)

    def _extract_qc_from_results_df(self, df: pd.DataFrame) -> Optional[Dict[str, Any]]:
        if df.empty:
            return None
        required = {
            "QC_ExpectedReplicates": "expected_replicates",
            "QC_FoundResultsReplicates": "found_results_replicates",
            "QC_MissingResultsReplicates": "missing_results_replicates",
            "QC_ReplicateComplete": "replicate_complete",
            "QC_GroupIdentityOK": "group_identity_ok",
        }
        if not any(col in df.columns for col in required):
            return None

        first = df.iloc[0]
        return {
            target: first.get(source)
            for source, target in required.items()
            if source in df.columns
        }

    def _extract_qc_from_sidecar(self, path: Path) -> Optional[Dict[str, Any]]:
        if not path.exists():
            return None
        try:
            df = pd.read_csv(path)
        except Exception as exc:
            logger.warning(f"Failed to parse MMPBSA QC sidecar {path}: {exc}")
            return None
        if df.empty:
            return None
        row = df.iloc[0]
        return {
            "expected_replicates": row.get("expected_replicates"),
            "found_results_replicates": row.get("found_results_replicates"),
            "missing_results_replicates": row.get("missing_results_replicates"),
            "replicate_complete": row.get("replicate_complete"),
            "group_identity_ok": row.get("group_identity_ok"),
        }

    def _compute_stats(self, data_dir: Path) -> Optional[Path]:
        """Compare per-residue distance profiles between ligands.

        Residues are matched across ligands, so this is a repeated-measures design:
        a Friedman test across all ligands and pairwise Wilcoxon signed-rank tests
        with Holm correction. The unit of observation is a residue's time- and
        replicate-averaged distance; the tests ask whether one ligand sits
        systematically closer to the selected residues than another. They are not
        a substitute for replicate-level uncertainty of a binding metric.
        """
        if not self.config.distance_summary_file or not self.config.distance_summary_file.exists():
            return None

        dist_df = pd.read_csv(self.config.distance_summary_file)
        if not {"Ligand", "MeanDistance"}.issubset(dist_df.columns):
            return None
        residue_col = "Residue" if "Residue" in dist_df.columns else None
        if residue_col is None:
            logger.warning("Distance summary has no Residue column; skipping paired statistics.")
            return None

        wide = dist_df.pivot_table(index=residue_col, columns="Ligand", values="MeanDistance", aggfunc="mean")
        wide = wide.dropna(axis=0, how="any")
        ligands = list(wide.columns)
        if len(ligands) < 2 or len(wide) < 3:
            return None

        stats_rows = []
        if len(ligands) >= 3:
            try:
                chi2, p_val = scipy_stats.friedmanchisquare(*[wide[lig].values for lig in ligands])
                stats_rows.append({
                    "Test": "Friedman",
                    "Statistic": chi2,
                    "p": p_val,
                    "p_holm": p_val,
                    "N_Residues": len(wide),
                })
            except Exception as exc:
                logger.warning(f"Friedman test failed: {exc}")

        pairwise = []
        for i in range(len(ligands)):
            for j in range(i + 1, len(ligands)):
                a = wide[ligands[i]].values
                b = wide[ligands[j]].values
                try:
                    stat, p_val = scipy_stats.wilcoxon(a, b)
                except Exception:
                    continue
                pairwise.append({
                    "Test": "Wilcoxon signed-rank (paired by residue)",
                    "GroupA": ligands[i],
                    "GroupB": ligands[j],
                    "Statistic": stat,
                    "MedianDiff": float(np.median(a - b)),
                    "p": p_val,
                    "N_Residues": len(wide),
                })

        # Holm-Bonferroni correction across pairwise comparisons.
        m = len(pairwise)
        running = 0.0
        for rank, idx in enumerate(sorted(range(m), key=lambda k: pairwise[k]["p"])):
            running = max(running, min(1.0, (m - rank) * pairwise[idx]["p"]))
            pairwise[idx]["p_holm"] = running
        stats_rows.extend(pairwise)

        if not stats_rows:
            return None

        stats_path = data_dir / "stats_summary.csv"
        pd.DataFrame(stats_rows).to_csv(stats_path, index=False)
        return stats_path
