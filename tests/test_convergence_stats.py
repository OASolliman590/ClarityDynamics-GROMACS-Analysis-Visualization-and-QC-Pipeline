"""Sampling checks and pairwise ligand tests (specs/003-run-management-and-rigor, US7/US8)."""

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gromacs_analysis.ranking.stats import (
    benjamini_hochberg, hedges_g, holm, pairwise_tests, permutation_p, welch,
)
from gromacs_analysis.utils.convergence import (
    block_average, cosine_content, half_agreement, half_split_rmsip, projection_convergence,
    series_convergence,
)


# ---------------------------------------------------------------- convergence

def test_block_sem_white_noise_matches_naive():
    rng = np.random.default_rng(1)
    x = rng.normal(0, 2.0, 4096)
    res = block_average(x)
    assert res["sem"] == pytest.approx(2.0 / math.sqrt(4096), rel=0.35)
    assert res["n_effective"] > 1000


def test_block_sem_correlated_series_has_fewer_effective_samples():
    rng = np.random.default_rng(2)
    noise = rng.normal(size=4096)
    ar = np.zeros_like(noise)
    for i in range(1, len(ar)):
        ar[i] = 0.95 * ar[i - 1] + noise[i]
    assert block_average(ar)["n_effective"] < block_average(noise)["n_effective"] / 5


def test_series_convergence_flags_drift_and_low_neff():
    rng = np.random.default_rng(3)
    flat = rng.normal(2.0, 0.1, 2000)
    assert series_convergence(flat)["Flags"] == ""
    drifting = np.linspace(1.0, 3.0, 2000) + rng.normal(0, 0.05, 2000)
    assert "drift" in series_convergence(drifting)["Flags"].split(";")
    assert half_agreement(drifting)["half_diff"] > 0.9
    short = rng.normal(size=12)
    assert "low_neff" in series_convergence(short)["Flags"]
    assert series_convergence([1.0, 2.0, 3.0])["Flags"] == "too_short"


def test_cosine_content_random_walk_vs_noise():
    rng = np.random.default_rng(4)
    walk = np.cumsum(rng.normal(size=(1000, 60)), axis=0)
    centred = walk - walk.mean(axis=0)
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    proj = centred @ vt[:3].T
    cos = projection_convergence(proj)
    assert cos["Cos_PC1"] > 0.7
    noise = rng.normal(size=1000)
    assert cosine_content(noise, 1) < 0.05
    assert cosine_content(np.cos(np.pi * (np.arange(500) + 0.5) / 500), 1) == pytest.approx(1.0, abs=1e-6)


def test_half_split_rmsip():
    rng = np.random.default_rng(5)
    basis = np.linalg.qr(rng.normal(size=(30, 3)))[0]
    flat = rng.normal(size=(400, 3)) * [5, 3, 2] @ basis.T + rng.normal(0, 0.05, (400, 30))
    assert half_split_rmsip(flat, k=3) > 0.95
    first = rng.normal(size=(200, 3)) * [5, 3, 2] @ basis.T
    other = np.linalg.qr(rng.normal(size=(30, 3)))[0]
    second = rng.normal(size=(200, 3)) * [5, 3, 2] @ other.T
    assert half_split_rmsip(np.vstack([first, second]), k=3) < 0.6


def test_md_stage_writes_convergence_table(tmp_path: Path):
    from gromacs_analysis.md import MDAnalyzer, MDConfig, SystemConfig

    cfg = MDConfig(base_dir=tmp_path, output_dir=tmp_path / "out", protein_name="P",
                   systems=[SystemConfig(name="A")], convergence={"neff_min": 10})
    analyzer = MDAnalyzer(cfg)
    rng = np.random.default_rng(6)
    data = {
        "LigA": {
            "rmsd_prot": {"mode": "replicate", "replicates_raw": [rng.normal(2, 0.1, 500), np.linspace(1, 3, 500)]},
            "rmsf": {"mode": "single", "values": rng.normal(size=300)},
        },
        "APO": {"rmsd_prot": {"mode": "single", "values": rng.normal(2, 0.1, 500)}},
    }
    path = analyzer._write_convergence(data)
    table = pd.read_csv(path, keep_default_na=False)
    assert set(table["Metric"]) == {"rmsd_prot"}  # per-residue profiles are not time series
    liga = table[table["System"] == "LigA"].set_index("Replicate")
    assert liga.loc[1, "Flags"] == "" and "drift" in liga.loc[2, "Flags"]
    assert list(table.columns[:4]) == ["System", "Metric", "Replicate", "N"]


def test_joint_pca_writes_per_replicate_convergence(tmp_path: Path):
    from gromacs_analysis.pca import PCAJointAnalyzer
    from gromacs_analysis.pca.config import PCAJointConfig, PCAJointSystem

    rng = np.random.default_rng(7)
    ref = rng.normal(size=(20, 3)) * 5
    stable = ref + rng.normal(0, 0.3, (300, 20, 3))
    walk = ref + np.cumsum(rng.normal(0, 0.3, (300, 20, 3)), axis=0)
    cfg = PCAJointConfig(
        output_dir=tmp_path / "pj", protein_name="P",
        systems=[PCAJointSystem("Stable", "t.pdb", ["a.xtc", "b.xtc"]), PCAJointSystem("Walk", "t.pdb", ["c.xtc"])],
    )
    cfg.plot_config.save_formats = ["html"]
    analyzer = PCAJointAnalyzer(cfg)
    analyzer.segments = {"Stable": [150, 150], "Walk": [300]}
    analyzer.run_analysis({"Stable": stable, "Walk": walk})
    table = pd.read_csv(tmp_path / "pj" / "data" / "convergence.csv", keep_default_na=False)
    assert list(table["System"]) == ["Stable", "Stable", "Walk"]
    walk_row = table[table["System"] == "Walk"].iloc[0]
    assert "pc1_diffusive" in walk_row["Flags"] and walk_row["Cos_PC1"] > 0.5
    assert (table[table["System"] == "Stable"]["Cos_PC1"] < 0.5).all()


# ---------------------------------------------------------------- statistics

def test_welch_matches_scipy_and_handles_constant_groups():
    from scipy import stats

    a, b = [1.0, 1.4, 0.9, 1.2], [2.0, 2.5, 1.8]
    expected = stats.ttest_ind(a, b, equal_var=False)
    assert welch(a, b)["p"] == pytest.approx(expected.pvalue)
    assert welch([1.0, 1.0], [1.0, 1.0])["p"] == 1.0
    assert welch([1.0, 1.0], [2.0, 2.0])["p"] == 0.0


def test_exact_permutation_floor_for_three_vs_three():
    assert permutation_p([1, 2, 3], [10, 11, 12]) == pytest.approx(0.1)
    assert permutation_p([1, 2, 3], [1, 2, 3]) == pytest.approx(1.0)
    big_a, big_b = np.arange(15), np.arange(15) + 100
    assert permutation_p(big_a, big_b, n_permutations=2000) < 0.01  # Monte-Carlo branch


def test_hedges_g_sign_and_size():
    assert hedges_g([1, 2, 3], [4, 5, 6]) < -2
    assert hedges_g([4, 5, 6], [1, 2, 3]) > 2
    assert hedges_g([1, 2, 3], [1, 2, 3]) == 0


def test_holm_and_bh_adjustments():
    p = [0.01, 0.04, 0.03, 0.005]
    assert holm(p) == pytest.approx([0.03, 0.06, 0.06, 0.02])
    assert benjamini_hochberg(p) == pytest.approx([0.02, 0.04, 0.04, 0.02])
    assert np.isnan(holm([np.nan, 0.01])[0]) and holm([np.nan, 0.01])[1] == pytest.approx(0.01)


def test_pairwise_tests_table():
    table = pairwise_tests(
        {"rmsd_lig": {"A": [1.0, 1.1, 0.9, 1.05], "B": [3.0, 3.2, 2.9, 3.1], "C": [1.0]}},
        directions={"rmsd_lig": "lower"},
    )
    ab = table[(table.LigandA == "A") & (table.LigandB == "B")].iloc[0]
    assert ab.p_holm < 0.05 and ab.Hedges_g < -5 and ab.Direction == "lower"
    ac = table[(table.LigandA == "A") & (table.LigandB == "C")].iloc[0]
    assert "insufficient" in ac.Note and np.isnan(ac.Welch_p)
    small = pairwise_tests({"m": {"A": [1, 2, 3], "B": [7, 8, 9]}})
    assert "0.1" in small.iloc[0].Note


def test_ranking_writes_pairwise_tests(tmp_path: Path):
    from gromacs_analysis.ranking import RankingAnalyzer
    from gromacs_analysis.ranking.config import RankingConfig

    md = tmp_path / "statistics_summary.csv"
    pd.DataFrame([
        {"System": "LigA", "Metric": "rmsd_lig", "Mean": 1.0, "SEM": 0.05, "Replicate_Means": "0.9;1.0;1.1"},
        {"System": "LigB", "Metric": "rmsd_lig", "Mean": 3.0, "SEM": 0.05, "Replicate_Means": "2.9;3.0;3.1"},
        {"System": "APO", "Metric": "rmsd_lig", "Mean": 9.0, "SEM": 0.05, "Replicate_Means": "9;9;9"},
    ]).to_csv(md, index=False)
    mm = tmp_path / "mm"
    mm.mkdir()
    for lig, vals in {"LigA": "-30;-31;-29", "LigB": "-20;-21;-19"}.items():
        pd.DataFrame([{"Component": "Binding_TOTAL", "Mean": float(vals.split(";")[0]), "SEM": 0.5,
                       "Replicate_Values": vals}]).to_csv(mm / f"{lig}_results.csv", index=False)
    result = RankingAnalyzer(RankingConfig(output_dir=tmp_path / "rank", protein_name="P", md_stats_file=md,
                                           mmpbsa_data_dir=mm, n_bootstrap=200)).run_analysis()
    tests = pd.read_csv(result["tests_file"])
    assert set(tests.Metric) == {"rmsd_lig", "mmpbsa_total"}
    assert set(tests.LigandA) | set(tests.LigandB) == {"LigA", "LigB"}  # APO excluded
    mmrow = tests[tests.Metric == "mmpbsa_total"].iloc[0]
    assert mmrow.Mean_A == pytest.approx(-30) and mmrow.Diff < 0
