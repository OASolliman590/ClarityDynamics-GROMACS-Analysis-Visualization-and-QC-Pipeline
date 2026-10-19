"""Round-3 tests (specs/001-pipeline-hardening-round3)."""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


# ---------------------------------------------------------------- US1 joint PCA

def _rotation(angle: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _synthetic_systems(seed: int = 0, n_frames: int = 300, n_atoms: int = 12):
    """Two systems that differ by a displacement of atom 0 along x, plus noise and random rotations."""
    rng = np.random.default_rng(seed)
    base = rng.normal(0, 5, (n_atoms, 3))
    systems = {}
    for name, shift in (("APO", 0.0), ("LigA", 3.0)):
        frames = np.repeat(base[None], n_frames, axis=0) + rng.normal(0, 0.15, (n_frames, n_atoms, 3))
        frames[:, 0, 0] += shift + rng.normal(0, 0.3, n_frames)
        rotated = np.einsum("fij,fkj->fki", np.stack([_rotation(a) for a in rng.uniform(0, 6, n_frames)]), frames)
        systems[name] = rotated + rng.normal(0, 4, (n_frames, 1, 3))  # random translations too
    return systems


def test_kabsch_fit_removes_rotation_and_translation():
    from gromacs_analysis.pca.joint import kabsch_fit

    rng = np.random.default_rng(1)
    ref = rng.normal(0, 3, (10, 3))
    moved = ref @ _rotation(1.1).T + np.array([5.0, -2.0, 1.0])
    fitted = kabsch_fit(moved, ref)
    assert np.allclose(fitted - fitted.mean(0), ref - ref.mean(0), atol=1e-8)


def test_joint_pca_recovers_displacement_on_shared_axes():
    from gromacs_analysis.pca.joint import joint_pca

    result = joint_pca(_synthetic_systems(), n_components=5, n_overlap=3)
    pc1 = result.eigenvectors[:, 0].reshape(-1, 3)
    # PC1 is dominated by atom 0; after fitting the direction is shared with small rigid corrections.
    assert np.linalg.norm(pc1[0]) ** 2 > 0.8
    proj = result.projections
    assert set(proj) == {"APO", "LigA"}
    separation = abs(proj["APO"][:, 0].mean() - proj["LigA"][:, 0].mean())
    assert separation > 2.5
    assert np.all(np.diff(result.eigenvalues) <= 1e-12)
    assert abs(result.explained_variance.sum() - result.eigenvalues.sum() / result.total_variance) < 1e-9


def test_rmsip_properties():
    from gromacs_analysis.pca.joint import joint_pca, rmsip

    result = joint_pca(_synthetic_systems(seed=3), n_components=5, n_overlap=4)
    m = result.rmsip_pairwise
    assert np.allclose(np.diag(m), 1.0, atol=1e-6)
    assert np.allclose(m, m.T)
    assert np.all((m >= -1e-9) & (m <= 1 + 1e-9))
    v = np.linalg.qr(np.random.default_rng(0).normal(size=(30, 4)))[0]
    assert abs(rmsip(v, -v, 4) - 1.0) < 1e-9  # sign invariant


def test_joint_pca_rejects_mismatched_atom_counts():
    from gromacs_analysis.pca.joint import joint_pca

    systems = {"APO": np.zeros((5, 10, 3)), "LigA": np.zeros((5, 11, 3))}
    with pytest.raises(ValueError, match=r"APO.*10.*LigA.*11|LigA.*11.*APO.*10"):
        joint_pca(systems, n_components=3)


def test_joint_pca_analyzer_writes_contract_outputs(tmp_path: Path):
    from gromacs_analysis.pca import PCAJointAnalyzer, PCAJointConfig, PCAJointSystem

    systems = _synthetic_systems(seed=5, n_frames=60)
    cfg = PCAJointConfig(
        output_dir=tmp_path / "joint",
        protein_name="P",
        systems=[
            PCAJointSystem(name="APO", topology="t.pdb", trajectories=["a.xtc"], is_apo=True),
            PCAJointSystem(name="LigA", topology="t.pdb", trajectories=["b.xtc"]),
        ],
        n_components=4,
        n_overlap=3,
    )
    cfg.plot_config.save_formats = ["html"]
    result = PCAJointAnalyzer(cfg).run_analysis(systems)

    data = tmp_path / "joint" / "data"
    for name in ("eigenvalues.csv", "projections_APO.csv", "projections_LigA.csv",
                 "rmsip_pairwise.csv", "rmsip_to_joint.csv", "parameters.json"):
        assert (data / name).exists(), name
    for name in ("joint_projection_12", "joint_density_12", "rmsip_heatmap", "joint_scree"):
        assert (tmp_path / "joint" / "plots" / f"{name}.html").exists(), name
    eig = pd.read_csv(data / "eigenvalues.csv")
    assert list(eig.columns) == ["PC", "Eigenvalue", "ExplainedVariance", "Cumulative"]
    assert result["success"]


# ---------------------------------------------------------------- US3 units

def test_unit_parsing_and_factors():
    from gromacs_analysis.utils.units import conversion_factor, display_unit_from_label, parse_xvg_units

    header = '@    xaxis  label "Time (ps)"\n@    yaxis  label "Area (nm\\S2\\N)"\n'
    assert parse_xvg_units(header)[1] == "nm^2"
    assert parse_xvg_units('@ yaxis label "RMSD (nm)"\n')[1] == "nm"
    assert display_unit_from_label("Å") == "Å" and display_unit_from_label("Å²") == "Å^2"
    assert display_unit_from_label("Number of hydrogen bonds") is None
    assert conversion_factor("nm", "Å") == (10.0, "converted")
    assert conversion_factor("nm^2", "Å^2") == (100.0, "converted")
    assert conversion_factor(None, "Å") == (1.0, "assumed")
    assert conversion_factor("kJ/mol", "Å") == (1.0, "unrecognised")
    assert conversion_factor("Å", "Å") == (1.0, "matched")


def _metric(name, pattern, ylabel="Å", data_format="two_column"):
    from gromacs_analysis.md import AnalysisMetric

    return AnalysisMetric(name=name, file_pattern=pattern, title=name, ylabel=ylabel, data_format=data_format)


def test_md_converts_nm_xvg_to_angstrom_and_records_it(tmp_path: Path):
    from gromacs_analysis.md import MDAnalyzer, MDConfig, SystemConfig

    gmx = tmp_path / "data" / "gmx"
    vmd = tmp_path / "data" / "vmd"
    gmx.mkdir(parents=True)
    vmd.mkdir(parents=True)
    values_a = np.linspace(1.0, 3.0, 50)
    (gmx / "rmsd.xvg").write_text(
        '@    xaxis  label "Time (ps)"\n@    yaxis  label "RMSD (nm)"\n'
        + "\n".join(f"{i * 10} {v / 10:.5f}" for i, v in enumerate(values_a))
    )
    (vmd / "rmsd.xvg").write_text("\n".join(f"{i * 10} {v:.5f}" for i, v in enumerate(values_a)))
    (gmx / "sasa.xvg").write_text('@ yaxis label "Area (nm\\S2\\N)"\n' + "\n".join(f"{i} 150.0" for i in range(50)))

    cfg = MDConfig(
        base_dir=tmp_path / "data", output_dir=tmp_path / "out", protein_name="P",
        systems=[SystemConfig(name="GMX", dir_pattern="gmx"), SystemConfig(name="VMD", dir_pattern="vmd")],
        metrics=[_metric("rmsd_prot", "rmsd.xvg"), _metric("sasa", "sasa.xvg", ylabel="Å²")],
        time_unit_input="ps",
    )
    analyzer = MDAnalyzer(cfg)
    data = analyzer.processor.load_all_data()
    gmx_mean = data["GMX"]["rmsd_prot"]["stats"]["mean"]
    vmd_mean = data["VMD"]["rmsd_prot"]["stats"]["mean"]
    assert abs(gmx_mean - vmd_mean) / vmd_mean < 0.01
    assert data["GMX"]["sasa"]["stats"]["mean"] == pytest.approx(15000.0)

    analyzer._create_output_dirs()
    stats = pd.read_csv(analyzer._save_statistics_csv(data)).set_index(["System", "Metric"])
    assert stats.loc[("GMX", "rmsd_prot"), "Value_Unit_Status"] == "converted"
    assert stats.loc[("GMX", "rmsd_prot"), "Value_Unit_Factor"] == 10
    assert stats.loc[("VMD", "rmsd_prot"), "Value_Unit_Status"] == "assumed"
    assert stats.loc[("GMX", "sasa"), "Value_Unit_Factor"] == 100


def test_md_unit_conversion_can_be_disabled(tmp_path: Path):
    from gromacs_analysis.md import MDConfig, MDDataProcessor, SystemConfig

    d = tmp_path / "data" / "s"
    d.mkdir(parents=True)
    (d / "r.xvg").write_text('@ yaxis label "RMSD (nm)"\n' + "\n".join(f"{i} 0.2" for i in range(20)))
    cfg = MDConfig(base_dir=tmp_path / "data", output_dir=tmp_path / "o", protein_name="P",
                   systems=[SystemConfig(name="S", dir_pattern="s")], metrics=[_metric("rmsd_prot", "r.xvg")],
                   convert_units=False)
    assert MDDataProcessor(cfg).load_all_data()["S"]["rmsd_prot"]["stats"]["mean"] == pytest.approx(0.2)


# ---------------------------------------------------------------- US4 equilibration

def _ar1(n, phi, rng):
    x = np.zeros(n)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + rng.normal()
    return x


def test_statistical_inefficiency_white_noise_and_ar1():
    from gromacs_analysis.utils.equilibration import statistical_inefficiency

    rng = np.random.default_rng(0)
    assert statistical_inefficiency(rng.normal(size=20000)) == pytest.approx(1.0, abs=0.15)
    phi = 0.8
    expected = (1 + phi) / (1 - phi)
    assert statistical_inefficiency(_ar1(50000, phi, rng)) == pytest.approx(expected, rel=0.25)


def test_detect_equilibration_finds_end_of_relaxation():
    from gromacs_analysis.utils.equilibration import detect_equilibration

    n, relax = 2000, 400
    hits = 0
    for seed in range(10):
        rng = np.random.default_rng(seed)
        t = np.arange(n)
        series = 3.0 * np.exp(-t / (relax / 4)) + 0.3 * _ar1(n, 0.5, rng)
        t0, g, n_eff = detect_equilibration(series)
        assert g >= 1.0 and 0 < n_eff <= n
        if abs(t0 - relax) <= 0.15 * n:
            hits += 1
    assert hits >= 9


def test_detect_equilibration_short_series_is_reported():
    from gromacs_analysis.utils.equilibration import detect_equilibration

    t0, g, n_eff = detect_equilibration(np.arange(5.0))
    assert t0 == 0 and np.isnan(n_eff)


def test_md_stability_chodera_writes_equilibration_json(tmp_path: Path):
    import json

    from gromacs_analysis.md import MDAnalyzer, MDConfig, SystemConfig
    from gromacs_analysis.md.config import StabilityConfig

    d = tmp_path / "data" / "s"
    d.mkdir(parents=True)
    rng = np.random.default_rng(1)
    t = np.arange(1500)
    series = 2.0 + 2.0 * np.exp(-t / 75.0) + 0.1 * _ar1(1500, 0.5, rng)
    (d / "r.dat").write_text("\n".join(f"{i} {v:.5f}" for i, v in zip(t, series)))
    cfg = MDConfig(
        base_dir=tmp_path / "data", output_dir=tmp_path / "out", protein_name="P",
        systems=[SystemConfig(name="S", dir_pattern="s")], metrics=[_metric("rmsd_prot", "r.dat")],
        time_unit_input="frame", time_unit_output="frame",
        stability=StabilityConfig(enabled=True, metric="rmsd_prot", min_points=50),
    )
    result = MDAnalyzer(cfg).run_analysis()
    records = json.loads((tmp_path / "out" / "data" / "equilibration.json").read_text())
    assert result["success"] and records[0]["system"] == "S"
    assert records[0]["method"] == "chodera" and records[0]["status"] == "detected"
    assert 150 <= records[0]["start_index"] <= 600
    assert records[0]["n_effective"] > 0


def test_mmpbsa_start_frame_excludes_early_frames(tmp_path: Path):
    from gromacs_analysis.mmpbsa import MMPBSAConfig, MMPBSAProcessor, MMPBSASystemConfig

    d = tmp_path / "LigA"
    d.mkdir()
    frames = np.arange(1, 11)
    total = np.where(frames <= 5, -10.0, -30.0)
    pd.DataFrame({"Frame #": frames, "VDWAALS": total - 5, "GGAS": total - 2, "GSOLV": np.full(10, 2.0),
                  "TOTAL": total}).to_csv(d / "FINAL_RESULTS_MMPBSA.csv", index=False)
    cfg = MMPBSAConfig(
        base_dir=tmp_path, output_dir=tmp_path / "out", protein_name="P", ligand_name="L",
        systems=[MMPBSASystemConfig(name="LigA", dir_pattern="LigA", results_file_pattern="FINAL_RESULTS_MMPBSA.csv")],
        start_frame=6,
    )
    data = MMPBSAProcessor(cfg).load_all_data()["LigA"]
    assert data["results"]["binding_energy"]["TOTAL"]["mean"] == pytest.approx(-30.0)
    assert data["qc"]["frames_used"] == [5]


def test_vmd_helper_rmsf_first_frame():
    from gromacs_analysis.pipeline.upstream import _build_md_vmd_script

    script = _build_md_vmd_script({"rmsf_first_frame": 250})
    assert "measure rmsf $sel first $rmsf_first_frame" in script
    assert "set rmsf_first_frame 250" in script


# ---------------------------------------------------------------- US5 ranking with uncertainty

def _mmpbsa_rows(tmp_path: Path, values: dict) -> Path:
    d = tmp_path / "mm"
    d.mkdir(exist_ok=True)
    for lig, (mean, sem) in values.items():
        pd.DataFrame([{"Component": "Binding_TOTAL", "Mean": mean, "Std": sem * 3 ** 0.5, "SEM": sem, "N": 3,
                       "Uncertainty": "between_replicates"}]).to_csv(d / f"{lig}_results.csv", index=False)
    return d


def test_bootstrap_rank_probabilities(tmp_path: Path):
    from gromacs_analysis.ranking import RankingAnalyzer, RankingConfig

    mm = _mmpbsa_rows(tmp_path, {"LigA": (-40.0, 1.0), "LigB": (-30.0, 1.0), "LigC": (-30.0, 1.0)})
    result = RankingAnalyzer(RankingConfig(output_dir=tmp_path / "r", protein_name="P", mmpbsa_data_dir=mm,
                                           n_bootstrap=2000, seed=0)).run_analysis()
    ranking = pd.read_csv(result["ranking_file"]).set_index("Ligand")
    pairs = pd.read_csv(tmp_path / "r" / "data" / "rank_pairwise.csv")
    p = {(a, b): v for a, b, v in pairs[["LigandA", "LigandB", "P_A_better"]].itertuples(index=False)}

    assert ranking.loc["LigA", "P_rank1"] >= 0.99
    assert p[("LigA", "LigB")] >= 0.99
    assert 0.4 <= p[("LigB", "LigC")] <= 0.6
    assert ranking.loc["LigA", "median_rank"] == 1
    assert ranking.loc["LigB", "rank_ci_high"] <= 3


def test_bootstrap_is_deterministic_with_seed(tmp_path: Path):
    from gromacs_analysis.ranking import RankingAnalyzer, RankingConfig

    mm = _mmpbsa_rows(tmp_path, {"LigA": (-31.0, 2.0), "LigB": (-30.0, 2.0)})
    outputs = []
    for run in range(2):
        cfg = RankingConfig(output_dir=tmp_path / f"r{run}", protein_name="P", mmpbsa_data_dir=mm, seed=7)
        outputs.append(pd.read_csv(RankingAnalyzer(cfg).run_analysis()["ranking_file"]))
    pd.testing.assert_frame_equal(outputs[0], outputs[1])


def test_reference_systems_are_excluded_but_similar_names_kept(tmp_path: Path):
    from gromacs_analysis.ranking import RankingAnalyzer, RankingConfig

    md = tmp_path / "md.csv"
    pd.DataFrame([
        {"System": name, "Metric": "rmsd_lig", "Mean": v, "SEM": 0.1}
        for name, v in (("APO", 0.5), ("Apomorphine", 2.0), ("Ref", 1.0), ("LigA", 1.5))
    ]).to_csv(md, index=False)
    cfg = RankingConfig(output_dir=tmp_path / "r", protein_name="P", md_stats_file=md, exclude_systems=["Ref"])
    ranking = pd.read_csv(RankingAnalyzer(cfg).run_analysis()["ranking_file"])
    excluded = pd.read_csv(tmp_path / "r" / "data" / "excluded_systems.csv")

    assert set(ranking["Ligand"]) == {"Apomorphine", "LigA"}
    assert set(excluded["System"]) == {"APO", "Ref"}


def test_ranking_with_only_reference_systems_fails(tmp_path: Path):
    from gromacs_analysis.ranking import RankingAnalyzer, RankingConfig

    md = tmp_path / "md.csv"
    pd.DataFrame([{"System": "APO", "Metric": "rmsd_lig", "Mean": 1.0}]).to_csv(md, index=False)
    with pytest.raises(ValueError, match="reference"):
        RankingAnalyzer(RankingConfig(output_dir=tmp_path / "r", protein_name="P", md_stats_file=md)).run_analysis()


def test_bootstrap_uses_replicate_means(tmp_path: Path):
    from gromacs_analysis.ranking import RankingAnalyzer

    md = tmp_path / "md.csv"
    pd.DataFrame([
        {"System": "LigA", "Metric": "rmsd_lig", "Mean": 1.0, "SEM": 0.0, "Replicate_Means": "0.5;1.0;1.5"},
        {"System": "LigB", "Metric": "rmsd_lig", "Mean": 1.0, "SEM": 0.0, "Replicate_Means": "1.0;1.0;1.0"},
    ]).to_csv(md, index=False)
    analyzer = RankingAnalyzer.__new__(RankingAnalyzer)
    samples = analyzer._uncertainty_from_md(pd.read_csv(md))
    assert samples["LigA"]["rmsd_lig"]["replicates"] == [0.5, 1.0, 1.5]


def test_md_statistics_include_replicate_means(tmp_path: Path):
    from gromacs_analysis.md import MDAnalyzer, MDConfig, SystemConfig

    for rep, level in ((1, 1.0), (2, 3.0)):
        d = tmp_path / "data" / f"s_{rep}"
        d.mkdir(parents=True)
        (d / "m.dat").write_text("\n".join(f"{t} {level}" for t in range(20)))
    cfg = MDConfig(base_dir=tmp_path / "data", output_dir=tmp_path / "o", protein_name="P",
                   systems=[SystemConfig(name="S", dir_pattern="s_{}", replicates=2)],
                   metrics=[_metric("rmsd_prot", "m.dat")])
    analyzer = MDAnalyzer(cfg)
    analyzer._create_output_dirs()
    stats = pd.read_csv(analyzer._save_statistics_csv(analyzer.processor.load_all_data()))
    assert stats.loc[0, "Replicate_Means"] == "1;3"


# ---------------------------------------------------------------- US6 protein-ligand H-bonds

def test_vmd_helper_writes_protein_ligand_hbonds():
    from gromacs_analysis.pipeline.upstream import _build_md_vmd_script

    script = _build_md_vmd_script({"ligands": ["INH", {"name": "NAD", "selection": "resname NAD"}]})
    assert 'calc_hbonds_ligand $ligand_selection "hbonds_lig_.dat"' in script
    assert '"hbonds_lig_INH.dat"' in script and '"hbonds_lig_NAD.dat"' in script
    assert "-sel2 $lig" in script


def test_hbonds_lig_metric_is_default_and_split_per_ligand(tmp_path: Path):
    from gromacs_analysis.md import MDConfig, SystemConfig

    single = MDConfig(base_dir=tmp_path, output_dir=tmp_path, protein_name="P", systems=[SystemConfig(name="L")])
    metric = next(m for m in single.metrics if m.name == "hbonds_lig")
    assert metric.is_holo_only and "hbonds_lig" in metric.file_pattern
    multi = MDConfig(base_dir=tmp_path, output_dir=tmp_path, protein_name="P", systems=[SystemConfig(name="L")],
                     ligands=["INH", "NAD"])
    names = {m.name for m in multi.metrics}
    assert {"hbonds_lig_INH", "hbonds_lig_NAD"} <= names and "hbonds_lig" not in names
    assert "hbonds" in names  # intra-protein count is kept


# ---------------------------------------------------------------- US7 flat-profile terminal threshold

def test_flat_core_terminals_fully_trimmed():
    from gromacs_analysis.pca.config import PCATerminalDetectionConfig
    from gromacs_analysis.pca.prep import ResidueMapEntry, detect_terminals

    values = np.full(120, 1.0)
    values[:6] = np.linspace(8, 3, 6)
    values[-10:] = np.linspace(3, 9, 10)
    entries = [ResidueMapEntry(rmsf_row=i + 1, rmsf_resid=i + 1, rmsf_value=float(v), gro_resid=i + 1,
                               native_resid=i + 1, resname="ALA") for i, v in enumerate(values)]
    det = detect_terminals(entries, PCATerminalDetectionConfig(smoothing_window=1))
    assert det.n_terminal_trimmed >= 6 and det.c_terminal_trimmed >= 10
    assert det.n_terminal_trimmed <= 7 and det.c_terminal_trimmed <= 11


def test_nondegenerate_threshold_unchanged():
    from gromacs_analysis.pca.prep import robust_threshold

    values = list(np.random.default_rng(0).normal(1.0, 0.1, 200))
    threshold, median, mad = robust_threshold(values, 3.0)
    raw_mad = float(np.median(np.abs(np.array(values) - np.median(values))))
    assert mad == pytest.approx(1.4826 * raw_mad)
    assert threshold == pytest.approx(median + 3.0 * mad)


def test_intra_protein_hbonds_pattern_ignores_ligand_files(tmp_path: Path):
    from gromacs_analysis.md import MDConfig, MDDataProcessor, SystemConfig

    d = tmp_path / "s"
    d.mkdir()
    (d / "hbonds_lig_.dat").write_text("0 2\n1 2\n")
    (d / "hbonds_.dat").write_text("0 50\n1 50\n")
    cfg = MDConfig(base_dir=tmp_path, output_dir=tmp_path / "o", protein_name="P",
                   systems=[SystemConfig(name="L", dir_pattern="s")])
    data = MDDataProcessor(cfg).load_all_data()["L"]
    assert data["hbonds"]["stats"]["mean"] == 50
    assert data["hbonds_lig"]["stats"]["mean"] == 2


# ---------------------------------------------------------------- US8 chain-aware distances

def _distance_processor(tmp_path: Path, residue_ids):
    pytest.importorskip("MDAnalysis")
    from trajectory_factory import build_system

    from gromacs_analysis.distance.config import DistanceConfig
    from gromacs_analysis.distance.processor import DistanceProcessor

    topology, trajectory = build_system(tmp_path / "traj", n_frames=5)
    cfg = DistanceConfig(base_dir=tmp_path, output_dir=tmp_path / "out", protein_name="P", ligand_name="LIG",
                         topology=topology, trajectories=[trajectory], ligand_selection="resname LIG",
                         residue_ids=residue_ids, method="min")
    return DistanceProcessor(cfg)


def test_chain_qualified_residue_selects_one_chain(tmp_path: Path):
    data = _distance_processor(tmp_path, ["A:4"]).load_data()
    assert [r["label"] for r in data["residues"]] == ["A:GLY4"]
    assert data["residues"][0]["resid"] == "A:4" and data["residues"][0]["chain"] == "A"


def test_unqualified_residue_in_two_chains_reported_separately(tmp_path: Path, caplog):
    with caplog.at_level(logging.WARNING):
        data = _distance_processor(tmp_path, [4]).load_data()
    labels = [r["label"] for r in data["residues"]]
    assert labels == ["A:GLY4", "B:GLY4"]
    assert not np.allclose(data["residues"][0]["mean"], data["residues"][1]["mean"])
    assert any("exists in chains" in rec.message for rec in caplog.records)


def test_unknown_chain_raises(tmp_path: Path):
    with pytest.raises(ValueError, match="chain 'Z'"):
        _distance_processor(tmp_path, ["Z:4"]).load_data()
