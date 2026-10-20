"""Regression tests for the statistical and physical correctness fixes."""

from pathlib import Path

import numpy as np
import pandas as pd

from gromacs_analysis.md import AnalysisMetric, MDAnalyzer, MDConfig, MDDataProcessor, SystemConfig
from gromacs_analysis.md.config import StabilityConfig
from gromacs_analysis.mmpbsa import MMPBSAParser
from gromacs_analysis.qc import QCAnalyzer, QCConfig
from gromacs_analysis.ranking import RankingAnalyzer, RankingConfig


def _metric(name: str, pattern: str, data_format: str = "two_column") -> AnalysisMetric:
    return AnalysisMetric(name=name, file_pattern=pattern, title=name, ylabel="value", data_format=data_format)


def test_ranking_uses_total_not_gas_phase_term(tmp_path: Path):
    mmpbsa_dir = tmp_path / "mmpbsa"
    mmpbsa_dir.mkdir()
    # LigA has the better ΔTOTAL but the worse ΔGGAS; ranking must follow ΔTOTAL.
    for lig, ggas, total in (("LigA", -10.0, -40.0), ("LigB", -90.0, -20.0)):
        pd.DataFrame(
            [
                {"Component": "Binding_GGAS", "Mean": ggas},
                {"Component": "Binding_GSOLV", "Mean": total - ggas},
                {"Component": "Binding_TOTAL", "Mean": total},
            ]
        ).to_csv(mmpbsa_dir / f"{lig}_results.csv", index=False)

    result = RankingAnalyzer(
        RankingConfig(output_dir=tmp_path / "rank", protein_name="P", mmpbsa_data_dir=mmpbsa_dir)
    ).run_analysis()
    df = pd.read_csv(result["ranking_file"]).set_index("Ligand")

    assert df.loc["LigA", "mmpbsa_total"] == -40.0
    assert df.loc["LigB", "mmpbsa_total"] == -20.0
    assert df.index[0] == "LigA"


def test_ranking_does_not_score_directionless_metrics(tmp_path: Path):
    md_stats = tmp_path / "md.csv"
    pd.DataFrame(
        [
            {"System": "LigA", "Metric": "sasa", "Mean": 100.0},
            {"System": "LigB", "Metric": "sasa", "Mean": 200.0},
            {"System": "LigA", "Metric": "rmsd_lig", "Mean": 3.0},
            {"System": "LigB", "Metric": "rmsd_lig", "Mean": 1.0},
        ]
    ).to_csv(md_stats, index=False)

    result = RankingAnalyzer(
        RankingConfig(output_dir=tmp_path / "rank", protein_name="P", md_stats_file=md_stats)
    ).run_analysis()
    df = pd.read_csv(result["ranking_file"])

    assert "sasa_z" not in df.columns
    assert df.iloc[0]["Ligand"] == "LigB"


def test_rmsf_profile_is_not_outlier_filtered(tmp_path: Path):
    sys_dir = tmp_path / "data/sys"
    sys_dir.mkdir(parents=True)
    rmsf = np.full(60, 1.0)
    rmsf[-1] = 12.0  # flexible terminus: real signal, not an outlier
    (sys_dir / "RMSF_.dat").write_text("\n".join(f"{i + 1} {v}" for i, v in enumerate(rmsf)))

    config = MDConfig(
        base_dir=tmp_path / "data",
        output_dir=tmp_path / "out",
        protein_name="P",
        systems=[SystemConfig(name="LigA", dir_pattern="sys", replicates=1)],
        metrics=[_metric("rmsf", "RMSF_.dat")],
        outlier_threshold=3.0,
    )
    data = MDDataProcessor(config).load_all_data()["LigA"]["rmsf"]

    assert len(data["values"]) == 60
    assert data["values"][-1] == 12.0
    assert np.allclose(data["time"], np.arange(1, 61))


def test_replicate_outliers_are_masked_without_shifting_frames(tmp_path: Path):
    series = {1: np.ones(40), 2: np.ones(40)}
    series[1][10] = 50.0
    for rep, values in series.items():
        rep_dir = tmp_path / "data" / f"sys_{rep}"
        rep_dir.mkdir(parents=True)
        (rep_dir / "m.dat").write_text("\n".join(f"{t} {v}" for t, v in enumerate(values)))

    config = MDConfig(
        base_dir=tmp_path / "data",
        output_dir=tmp_path / "out",
        protein_name="P",
        systems=[SystemConfig(name="LigA", dir_pattern="sys_{}", replicates=2)],
        metrics=[_metric("m", "m.dat")],
        time_unit_input="frame",
        time_unit_output="frame",
        outlier_threshold=3.0,
    )
    data = MDDataProcessor(config).load_all_data()["LigA"]["m"]

    assert data["n_points"] == 40
    assert np.isnan(data["replicates_raw"][0][10])
    assert np.allclose(data["mean"], 1.0)


def test_replicate_statistics_report_between_replicate_sem(tmp_path: Path):
    for rep, level in ((1, 1.0), (2, 2.0), (3, 3.0)):
        rep_dir = tmp_path / "data" / f"sys_{rep}"
        rep_dir.mkdir(parents=True)
        (rep_dir / "m.dat").write_text("\n".join(f"{t} {level}" for t in range(20)))

    config = MDConfig(
        base_dir=tmp_path / "data",
        output_dir=tmp_path / "out",
        protein_name="P",
        systems=[SystemConfig(name="LigA", dir_pattern="sys_{}", replicates=3)],
        metrics=[_metric("m", "m.dat")],
    )
    stats = MDDataProcessor(config).load_all_data()["LigA"]["m"]["stats"]

    assert stats["mean"] == 2.0
    assert np.isclose(stats["between_replicate_sd"], 1.0)
    assert np.isclose(stats["sem"], 1.0 / np.sqrt(3))
    assert stats["sem_method"] == "between_replicates"


def test_block_sem_exceeds_naive_sem_for_correlated_series(tmp_path: Path):
    rng = np.random.default_rng(0)
    x = np.zeros(4096)
    for i in range(1, len(x)):
        x[i] = 0.98 * x[i - 1] + rng.normal()

    processor = MDDataProcessor(MDConfig(base_dir=tmp_path, output_dir=tmp_path, protein_name="P", systems=[]))
    result = processor._block_average_sem(x)
    naive = np.std(x, ddof=1) / np.sqrt(len(x))

    assert result["sem"] > 3 * naive
    assert result["statistical_inefficiency"] > 10


def test_stability_window_keeps_report_working_and_skips_rmsf(tmp_path: Path):
    sys_dir = tmp_path / "data/sys"
    sys_dir.mkdir(parents=True)
    rmsd = np.concatenate([np.linspace(0.0, 3.0, 100), np.full(200, 3.0)])
    (sys_dir / "RMSD_apo_.dat").write_text("\n".join(f"{t} {v}" for t, v in enumerate(rmsd)))
    (sys_dir / "RMSF_.dat").write_text("\n".join(f"{i + 1} 1.0" for i in range(50)))

    config = MDConfig(
        base_dir=tmp_path / "data",
        output_dir=tmp_path / "out",
        protein_name="P",
        systems=[SystemConfig(name="LigA", dir_pattern="sys", replicates=1)],
        metrics=[_metric("rmsd_prot", "RMSD_apo_.dat"), _metric("rmsf", "RMSF_.dat")],
        time_unit_input="frame",
        time_unit_output="frame",
        stability=StabilityConfig(enabled=True, metric="rmsd_prot", window=20, min_points=50),
    )
    analyzer = MDAnalyzer(config)
    data = analyzer.processor.load_all_data()["LigA"]

    assert set(data) == {"rmsd_prot", "rmsf"}
    assert data["rmsd_prot"]["stable_start_time"] >= 90
    assert len(data["rmsf"]["values"]) == 50

    analyzer._create_output_dirs()
    stats_file = analyzer._save_statistics_csv({"LigA": data})
    stats = pd.read_csv(stats_file)
    assert "SEM" in stats.columns


def test_xvg_header_time_unit_overrides_magnitude_heuristic(tmp_path: Path):
    sys_dir = tmp_path / "data/sys"
    sys_dir.mkdir(parents=True)
    # 1 ns run written in ps: magnitude alone (max <= 1000) would be read as ns.
    lines = ['@    xaxis  label "Time (ps)"'] + [f"{t} 0.1" for t in range(0, 1001, 100)]
    (sys_dir / "m.xvg").write_text("\n".join(lines))

    config = MDConfig(
        base_dir=tmp_path / "data",
        output_dir=tmp_path / "out",
        protein_name="P",
        systems=[SystemConfig(name="LigA", dir_pattern="sys", replicates=1)],
        metrics=[_metric("m", "m.xvg")],
        time_unit_output="ns",
    )
    data = MDDataProcessor(config).load_all_data()["LigA"]["m"]

    assert np.isclose(data["time"][-1], 1.0)
    assert data["time_metadata"]["inferred_input_unit"] == "ps"


def test_mmpbsa_dat_parser_reads_sd_and_sem_columns(tmp_path: Path):
    content = """Delta (Complex - Receptor - Ligand):
Energy Component       Average     SD(Prop.)         SD   SEM(Prop.)        SEM
-------------------------------------------------------------------------------------
ΔVDWAALS                -45.12          3.10       2.90         0.31       0.29
ΔGGAS                   -57.62          6.00       5.50         0.60       0.55
ΔGSOLV                   19.50          4.10       4.00         0.41       0.40
ΔTOTAL                  -38.12          4.75       3.41         0.48       0.34
-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------
"""
    components = MMPBSAParser()._extract_delta_components(content)

    assert components["ΔTOTAL"] == {"value": -38.12, "std_dev": 3.41, "std_err": 0.34}
    assert components["ΔVDWAALS"]["std_dev"] == 2.90


def test_qc_flags_irregular_frame_spacing(tmp_path: Path):
    traj = tmp_path / "rmsd.xvg"
    traj.write_text('@ xaxis label "Time (ps)"\n0 0.1\n10 0.1\n20 0.1\n20 0.1\n40 0.1\n')

    result = QCAnalyzer(
        QCConfig(base_dir=tmp_path, output_dir=tmp_path / "qc", trajectories=[traj])
    ).run_analysis()

    assert any("irregular frame interval" in w for w in result["warnings"])
