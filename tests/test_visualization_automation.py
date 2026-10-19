"""Tests for per-ligand automation, multi-ligand visualization and entropy plots."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import yaml

from gromacs_analysis.md import AnalysisMetric, MDConfig, MDPlotter, SystemConfig
from gromacs_analysis.mmpbsa import MMPBSAConfig, MMPBSAParser, MMPBSASystemConfig
from gromacs_analysis.mmpbsa.plotter import MMPBSAPlotter
from gromacs_analysis.mmpbsa.processor import MMPBSAProcessor
from gromacs_analysis.pca import PCAConfig, PCAPlotter
from gromacs_analysis.pipeline.discovery import detect_ligands, discover_systems
from gromacs_analysis.pipeline.yaml_loader import load_pipeline_config
from gromacs_analysis.prolif import ProlifAnalyzer, ProlifConfig
from gromacs_analysis.ttclust import TTClustAnalyzer, TTClustConfig
from gromacs_analysis.utils.figure_export import figure_size
from gromacs_analysis.utils.palette import APO_HEX, CATEGORICAL_HEX, hex_to_rgb

PDB = """ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  HSD A   2       3.800   0.000   0.000  1.00  0.00           C
HETATM    3  C1  INH A 101       5.000   1.000   0.000  1.00  0.00           C
HETATM    4  C1  NAD A 102       6.000   1.000   0.000  1.00  0.00           C
HETATM    5  OH2 TIP3W   1       9.000   1.000   0.000  1.00  0.00           O
HETATM    6 SOD  SOD I   1      12.000   1.000   0.000  1.00  0.00          NA
END
"""


def _make_layout(root: Path, ligands=("LigA", "LigB"), reps=(1, 2)):
    for lig in ligands:
        for rep in reps:
            d = root / f"HOLO_{lig}_{rep}"
            d.mkdir(parents=True)
            (d / "top.pdb").write_text(PDB)
            (d / "traj.xtc").write_text("")
    (root / "APO_1").mkdir()


def test_discover_systems_groups_replicates_per_ligand(tmp_path: Path):
    _make_layout(tmp_path, reps=(2, 1, 10))
    systems = discover_systems(tmp_path, "HOLO_{ligand}_{rep}", "top.pdb", "traj.xtc")

    assert [s["name"] for s in systems] == ["LigA", "LigB"]
    assert [Path(t).parent.name for t in systems[0]["trajectories"]] == ["HOLO_LigA_1", "HOLO_LigA_2", "HOLO_LigA_10"]


def test_detect_ligands_skips_protein_solvent_and_ions(tmp_path: Path):
    pdb = tmp_path / "top.pdb"
    pdb.write_text(PDB)
    assert [lig["name"] for lig in detect_ligands(pdb)] == ["INH", "NAD"]


def test_prolif_batch_auto_discovers_and_splits_ligands(tmp_path: Path):
    data = tmp_path / "data"
    _make_layout(data)
    cfg_path = tmp_path / "pipeline.yaml"
    cfg_path.write_text(yaml.safe_dump({
        "project": {"name": "demo", "protein_name": "ProteinX"},
        "paths": {"data_root": "./data", "output_root": "./out"},
        "prolif_batch": {
            "enabled": True,
            "base_dir": "./data",
            "auto_discover": {"pattern": "HOLO_{ligand}_{rep}", "topology": "top.pdb", "trajectory": "traj.xtc"},
            "ligands": "auto",
        },
        "prolif_compare": {"enabled": True, "from_prolif_batch": True},
    }))

    config = load_pipeline_config(cfg_path)
    stage = config.stages["prolif_batch"]
    names = [c.ligand_name for c in stage.config]
    selections = {c.ligand_name: c.ligand_selection for c in stage.config}

    assert names == ["LigA_INH", "LigA_NAD", "LigB_INH", "LigB_NAD"]
    assert selections["LigB_NAD"] == "resname NAD"
    assert all(len(c.trajectories) == 2 for c in stage.config)
    assert stage.config[0].output_dir.name == "output_LigA_INH"


def test_disabled_auto_discover_stage_does_not_break_loading(tmp_path: Path):
    cfg_path = tmp_path / "pipeline.yaml"
    cfg_path.write_text(yaml.safe_dump({
        "project": {"name": "demo", "protein_name": "ProteinX"},
        "paths": {"data_root": "./missing", "output_root": "./out"},
        "ttclust_batch": {
            "enabled": False,
            "auto_discover": {"pattern": "HOLO_{ligand}_{rep}", "topology": "t.pdb", "trajectory": "x.xtc"},
        },
    }))
    load_pipeline_config(cfg_path)


def test_prolif_occupancy_uses_protein_residue_not_ligand(tmp_path: Path):
    cfg = ProlifConfig(
        base_dir=tmp_path, output_dir=tmp_path / "o", protein_name="P", ligand_name="L",
        topology=tmp_path / "t.pdb", trajectories=[tmp_path / "x.xtc"],
    )
    analyzer = ProlifAnalyzer.__new__(ProlifAnalyzer)
    analyzer.config = cfg
    columns = pd.MultiIndex.from_tuples(
        [("UNK1", "ASP129.A", "HBDonor"), ("UNK1", "PHE45.A", "Hydrophobic")],
        names=["ligand", "protein", "interaction"],
    )
    df = pd.DataFrame([[1, 0], [1, 1], [0, 1], [1, 1]], columns=columns).astype(bool)

    occ = analyzer._compute_occupancy(df, cfg.interaction_types).set_index("Residue")
    assert occ.loc["ASP129.A", "Occupancy"] == 75.0
    assert occ.loc["PHE45.A", "Interaction"] == "Hydrophobic"
    assert analyzer._extract_protein_residue_from_column(("UNK1", "ASP129.A", "HBDonor")) == "ASP129"


def test_ttclust_log_parsing_and_summary_figures(tmp_path: Path):
    out = tmp_path / "tt"
    out.mkdir()
    (out / "clustering.log").write_text(
        "------- CLUSTERS -------\n"
        "cluster 1\n    size = 3\n    representative frame= 2\n    Members : 1 2 3 \n    spread  : 1.20 \n"
        "cluster 2\n    size = 1\n    representative frame=4\n    Members : 4 \n    spread  : 0.00 \n"
    )
    cfg = TTClustConfig(base_dir=tmp_path, output_dir=out, protein_name="P", system_name="LigA",
                        trajectories=[tmp_path / "x.xtc"])
    cfg.plot_config.save_formats = ["html"]
    analyzer = TTClustAnalyzer(cfg)
    (out / "data").mkdir()
    (out / "plots").mkdir()

    figures = analyzer._summarize_clusters(out / "data", out / "plots")
    summary = pd.read_csv(out / "data" / "clusters_summary.csv")

    assert list(summary["Size"]) == [3, 1]
    assert summary["Population_percent"].tolist() == [75.0, 25.0]
    assert {p.name for p in figures} == {"ttclust_populations.html", "ttclust_timeline.html"}


def _md_config(tmp_path: Path, **kwargs) -> MDConfig:
    systems = [SystemConfig(name="APO", is_apo=True)] + [SystemConfig(name=f"Lig{i}") for i in range(3)]
    return MDConfig(base_dir=tmp_path, output_dir=tmp_path / "out", protein_name="P", systems=systems, **kwargs)


def test_md_ligand_systems_get_distinct_palette_colors(tmp_path: Path):
    plotter = MDPlotter(_md_config(tmp_path))
    colors = [plotter._get_system_color(n) for n in ("APO", "Lig0", "Lig1", "Lig2")]

    assert colors[0] == hex_to_rgb(APO_HEX)
    assert colors[1:] == [hex_to_rgb(h) for h in CATEGORICAL_HEX[:3]]
    assert len({plotter._format_legend_name(n, "rmsd_prot") for n in ("Lig0", "Lig1")}) == 2


def test_md_multi_ligand_complex_splits_ligand_metrics(tmp_path: Path):
    config = _md_config(tmp_path, ligands=["INH", {"name": "NAD", "selection": "resname NAD"}])
    names = [m.name for m in config.metrics]

    assert "rmsd_lig" not in names and "comcom" not in names
    assert {"rmsd_lig_INH", "rmsd_lig_NAD", "comcom_INH", "comcom_NAD"} <= set(names)
    metric = next(m for m in config.metrics if m.name == "comcom_NAD")
    assert metric.base_name == "comcom" and "comcom_NAD.dat" in metric.file_pattern


def _series(values):
    return {"mode": "single", "values": np.asarray(values, float), "time": np.arange(1, len(values) + 1)}


def test_combined_rmsd_includes_every_complex(tmp_path: Path):
    plotter = MDPlotter(_md_config(tmp_path))
    data = {"APO": {"rmsd_prot": _series([1, 2])}}
    for i in range(3):
        data[f"Lig{i}"] = {"rmsd_complex": _series([1, 2])}

    fig = plotter.create_combined_rmsd_plot(data)
    assert [t.name for t in fig.data] == ["P · Apo", "P–Lig0 complex", "P–Lig1 complex", "P–Lig2 complex"]


def test_rmsf_terminals_are_faded_with_two_boundaries(tmp_path: Path):
    config = _md_config(tmp_path)
    plotter = MDPlotter(config)
    profile = 1.0 + 0.05 * np.random.default_rng(0).standard_normal(120)
    profile[:6] = np.linspace(8, 2, 6)
    profile[-10:] = np.linspace(2, 9, 10)
    metric = AnalysisMetric(name="rmsf", file_pattern="", title="RMSF", ylabel="Å")

    fig = plotter.create_plot({"APO": _series(profile)}, metric)
    boundaries = [s for s in fig.layout.shapes if s.type == "line" and s.x0 == s.x1]
    faded = [t for t in fig.data if t.line.color and t.line.color.endswith(",0.3)")]

    assert len(boundaries) == 2
    assert len(faded) == 2
    assert "flexible terminal" in fig.layout.annotations[0].text

    plotter._write_terminal_detection({"APO": _series(profile)}, None)
    detection = json.loads((config.output_dir / "data" / "rmsf_terminal_detection.json").read_text())
    assert detection["n_trimmed"] >= 6 and detection["c_trimmed"] >= 10


ENTROPY_DAT = """ENTROPY RESULTS (INTERACTION ENTROPY):
Model           σ(Int. Energy)      Average       SD        SEM
-------------------------------------------------------------------------------
GB                    4.20            7.50        0.30       0.05

"""


def test_mmpbsa_entropy_from_sibling_dat_and_dg(tmp_path: Path):
    sys_dir = tmp_path / "LigA"
    sys_dir.mkdir()
    pd.DataFrame({"Frame #": [1, 2], "VDWAALS": [-30, -32], "GGAS": [-40, -42], "GSOLV": [10, 10],
                  "TOTAL": [-30.0, -32.0]}).to_csv(sys_dir / "FINAL_RESULTS_MMPBSA.csv", index=False)
    (sys_dir / "FINAL_RESULTS_MMPBSA.dat").write_text(ENTROPY_DAT)
    config = MMPBSAConfig(
        base_dir=tmp_path, output_dir=tmp_path / "out", protein_name="P", ligand_name="L",
        systems=[MMPBSASystemConfig(name="LigA", dir_pattern="LigA", results_file_pattern="FINAL_RESULTS_MMPBSA.csv")],
    )

    results = MMPBSAProcessor(config).load_all_data()["LigA"]["results"]
    rows = MMPBSAPlotter(config).entropy_terms(results)

    assert results["entropy"]["interaction_entropy"]["mean"] == 7.5
    assert rows[0]["dg"] == -31.0 + 7.5
    assert rows[0]["reliable"] is False
    fig = MMPBSAPlotter(config).create_entropy_plot(results, "LigA")
    assert list(fig.data[0].x) == ["ΔTOTAL", "−TΔS (IE)", "ΔG (IE)"]


def test_mmpbsa_c2_entropy_has_no_sigma_verdict():
    text = """ENTROPY RESULTS (C2 ENTROPY):
Model           σ(Int. Energy)      C2 Value      SD        C2 Std. Dev. (95%)
-------------------------------------------------------------------------------
GB                    4.13            5.29        0.12        2.05-2.53

"""
    c2 = MMPBSAParser()._extract_entropy(text)["c2_entropy"]
    assert c2["value"] == 5.29 and c2["confidence_interval"] == "2.05-2.53" and c2["reliable"] is None


def test_pca_time_series_converts_ps_to_ns_and_uses_ligand_color(tmp_path: Path):
    plotter = PCAPlotter(PCAConfig(base_dir=tmp_path, output_dir=tmp_path, color="#eb6834", n_pcs=2))
    proj = {"xs": [[0, 1000, 2000]] * 2, "ys": [[0.1, 0.2, 0.3]] * 2, "titles": ["a", "b"],
            "x_titles": ["Time (ps)"] * 2}

    fig = plotter.create_pc_timeseries(proj)
    assert list(fig.data[0].x) == [0.0, 1.0, 2.0]
    assert fig.data[0].line.color == "#eb6834"


def test_static_export_respects_figure_size():
    fig = go.Figure()
    assert figure_size(fig, 1200, 600) == (1200, 600)
    fig.update_layout(width=900, height=1500)
    assert figure_size(fig, 1200, 600) == (900, 1500)
