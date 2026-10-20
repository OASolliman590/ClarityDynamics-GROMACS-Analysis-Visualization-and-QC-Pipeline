from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import plotly.io as pio
import pytest

from gromacs_analysis.md import AnalysisMetric, MDAnalyzer, MDConfig, MDDataProcessor, MDPlotter, SystemConfig
from gromacs_analysis.mmpbsa import MMPBSAConfig, MMPBSASystemConfig, MMPBSAParser, MMPBSAProcessor
from gromacs_analysis.pca import PCAConfig, PCAPlotter
from gromacs_analysis.pca.runner import PCAGromacsRunner
from gromacs_analysis.pipeline.utils import run_utils_tasks
from gromacs_analysis.prolif import ProlifAnalyzer, ProlifComparisonAnalyzer, ProlifCompareConfig, ProlifConfig
from gromacs_analysis.prolif.plotter import ProlifPlotter
from gromacs_analysis.ttclust import TTClustAnalyzer, TTClustConfig
from gromacs_analysis.pipeline import run_pipeline
from gromacs_analysis.pipeline.upstream import generate_upstream_scripts
from gromacs_analysis.pipeline.yaml_loader import load_pipeline_config
from gromacs_analysis.ranking import RankingAnalyzer, RankingConfig
from gromacs_analysis.utils.sequence import build_residue_labels, load_sequence_from_pdb


def test_pipeline_input_paths_do_not_duplicate_base_dir(tmp_path: Path):
    config_path = tmp_path / "pipeline.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "project": {"name": "demo", "protein_name": "ProteinX"},
                "paths": {"data_root": "./data", "output_root": "./out"},
                "run": ["distance", "interactions", "prolif"],
                "distance": {
                    "enabled": True,
                    "base_dir": "./data",
                    "output_dir": "./out/distance",
                    "topology": "./data/sys/top.psf",
                    "trajectories": ["./data/sys/full.xtc", "sys/base_relative.xtc"],
                    "residues": {"range": [1, 2]},
                },
                "interactions": {
                    "enabled": True,
                    "base_dir": "./data",
                    "output_dir": "./out/interactions",
                    "topology": "./data/sys/top.psf",
                    "trajectories": ["./data/sys/full.xtc", "sys/base_relative.xtc"],
                },
                "prolif": {
                    "enabled": True,
                    "base_dir": "./data",
                    "output_dir": "./out/prolif",
                    "topology": "./data/sys/top.psf",
                    "trajectories": ["./data/sys/full.xtc", "sys/base_relative.xtc"],
                },
            },
            sort_keys=False,
        )
    )

    cfg = load_pipeline_config(config_path)

    for stage_name in ["distance", "interactions", "prolif"]:
        stage_cfg = cfg.stages[stage_name].config
        assert stage_cfg.topology == tmp_path / "data/sys/top.psf"
        assert stage_cfg.trajectories[0] == tmp_path / "data/sys/full.xtc"
        assert stage_cfg.trajectories[1] == tmp_path / "data/sys/base_relative.xtc"


def test_mmpbsa_replicate_wildcards_are_preserved(tmp_path: Path):
    base = tmp_path / "mmpbsa"
    for idx in [1, 2, 3]:
        rep_dir = base / f"rep{idx}"
        rep_dir.mkdir(parents=True)
        (rep_dir / f"FINAL_RESULTS_MMPBSA_T{idx}.dat").write_text("result")
        (rep_dir / f"FINAL_DECOMP_MMPBSA_T{idx}.dat").write_text("decomp")

    config = MMPBSAConfig(
        base_dir=base,
        output_dir=tmp_path / "out",
        protein_name="ProteinX",
        ligand_name="LigandA",
        systems=[
            MMPBSASystemConfig(
                name="LigandA",
                dir_pattern="rep{}",
                replicates=3,
                results_file_pattern="FINAL_RESULTS_MMPBSA*.dat",
                decomp_file_pattern="FINAL_DECOMP_MMPBSA*.dat",
            )
        ],
    )

    results_files, decomp_files = MMPBSAProcessor(config)._find_system_files(config.systems[0])

    assert len(results_files) == 3
    assert len(decomp_files) == 3
    assert results_files[0].name == "FINAL_RESULTS_MMPBSA_T1.dat"


def test_mmpbsa_csv_results_preserve_frame_data(tmp_path: Path):
    csv_path = tmp_path / "FINAL_RESULTS_MMPBSA.csv"
    pd.DataFrame(
        {
            "Frame #": [1, 2],
            "VDWAALS": [-10.0, -12.0],
            "GGAS": [-15.0, -17.0],
            "GSOLV": [4.0, 5.0],
            "TOTAL": [-11.0, -12.0],
        }
    ).to_csv(csv_path, index=False)

    parsed = MMPBSAParser().parse_results(csv_path, file_format="csv")
    formatted = MMPBSAProcessor(
        MMPBSAConfig(
            base_dir=tmp_path,
            output_dir=tmp_path / "out",
            protein_name="ProteinX",
            ligand_name="LigandA",
            systems=[MMPBSASystemConfig(name="LigandA", dir_pattern=".")],
        )
    )._format_single_results(parsed)

    assert parsed["metadata"]["num_frames"] == 2
    assert list(parsed["frame_data"]["TOTAL"]) == [-11.0, -12.0]
    assert formatted["frame_data"] is not None


def test_mmpbsa_qc_metadata_captures_missing_replicates(tmp_path: Path):
    base = tmp_path / "mmpbsa"
    for idx in (1, 2):
        rep_dir = base / f"rep{idx}"
        rep_dir.mkdir(parents=True)
        pd.DataFrame(
            {
                "Frame #": [1, 2],
                "VDWAALS": [-10.0, -11.0],
                "GGAS": [-15.0, -16.0],
                "GSOLV": [4.0, 5.0],
                "TOTAL": [-11.0, -11.5],
            }
        ).to_csv(rep_dir / f"FINAL_RESULTS_MMPBSA_{idx}.csv", index=False)

    config = MMPBSAConfig(
        base_dir=base,
        output_dir=tmp_path / "out",
        protein_name="ProteinX",
        ligand_name="LigandA",
        file_format="csv",
        systems=[
            MMPBSASystemConfig(
                name="LigandA",
                dir_pattern="rep{}",
                replicates=3,
                results_file_pattern="FINAL_RESULTS_MMPBSA*.csv",
                decomp_file_pattern="FINAL_DECOMP_MMPBSA*.csv",
            )
        ],
    )

    data = MMPBSAProcessor(config).load_all_data()
    qc = data["LigandA"]["qc"]
    assert qc["expected_replicates"] == 3
    assert qc["found_results_replicates"] == 2
    assert qc["missing_results_replicates"] == 1
    assert qc["replicate_complete"] is False


def test_mmpbsa_policy_error_on_incomplete_results(tmp_path: Path):
    base = tmp_path / "mmpbsa"
    (base / "rep1").mkdir(parents=True)
    pd.DataFrame(
        {
            "Frame #": [1],
            "VDWAALS": [-10.0],
            "GGAS": [-15.0],
            "GSOLV": [4.0],
            "TOTAL": [-11.0],
        }
    ).to_csv(base / "rep1" / "FINAL_RESULTS_MMPBSA_1.csv", index=False)

    config = MMPBSAConfig(
        base_dir=base,
        output_dir=tmp_path / "out",
        protein_name="ProteinX",
        ligand_name="LigandA",
        file_format="csv",
        incomplete_replicates_policy="error",
        systems=[
            MMPBSASystemConfig(
                name="LigandA",
                dir_pattern="rep{}",
                replicates=2,
                results_file_pattern="FINAL_RESULTS_MMPBSA*.csv",
            )
        ],
    )

    with pytest.raises(ValueError):
        MMPBSAProcessor(config).load_all_data()


def test_md_smoothing_window_is_applied(tmp_path: Path):
    data_dir = tmp_path / "data/sys"
    data_dir.mkdir(parents=True)
    (data_dir / "values.dat").write_text("1\n2\n3\n4\n")

    config = MDConfig(
        base_dir=tmp_path / "data",
        output_dir=tmp_path / "out",
        protein_name="ProteinX",
        systems=[SystemConfig(name="LigandA", dir_pattern="sys", replicates=1)],
        metrics=[
            AnalysisMetric(
                name="toy",
                file_pattern="values.dat",
                title="Toy",
                ylabel="value",
                data_format="single_column",
            )
        ],
        smoothing_window=3,
        outlier_threshold=999.0,
    )

    data = MDDataProcessor(config).load_all_data()

    assert np.allclose(data["LigandA"]["toy"]["values"], [1.5, 2.0, 3.0, 3.5])


def test_md_two_column_time_axis_is_preserved(tmp_path: Path):
    data_dir = tmp_path / "data/sys"
    data_dir.mkdir(parents=True)
    (data_dir / "metric.dat").write_text("@ title \"metric\"\n0 1.0\n10 2.0\n20 3.0\n")

    config = MDConfig(
        base_dir=tmp_path / "data",
        output_dir=tmp_path / "out",
        protein_name="ProteinX",
        systems=[SystemConfig(name="LigandA", dir_pattern="sys", replicates=1)],
        metrics=[
            AnalysisMetric(
                name="metric",
                file_pattern="metric.dat",
                title="Metric",
                ylabel="value",
                data_format="two_column",
            )
        ],
        outlier_threshold=999.0,
    )

    data = MDDataProcessor(config).load_all_data()
    assert np.allclose(data["LigandA"]["metric"]["time"], [0.0, 10.0, 20.0])


def test_md_replicate_time_axis_is_preserved(tmp_path: Path):
    for rep in (1, 2):
        rep_dir = tmp_path / "data" / f"sys_{rep}"
        rep_dir.mkdir(parents=True)
        (rep_dir / "metric.dat").write_text("0 1.0\n5 2.0\n10 3.0\n")

    config = MDConfig(
        base_dir=tmp_path / "data",
        output_dir=tmp_path / "out",
        protein_name="ProteinX",
        systems=[SystemConfig(name="LigandA", dir_pattern="sys_{}", replicates=2)],
        metrics=[
            AnalysisMetric(
                name="metric",
                file_pattern="metric.dat",
                title="Metric",
                ylabel="value",
                data_format="two_column",
            )
        ],
        outlier_threshold=999.0,
    )

    data = MDDataProcessor(config).load_all_data()
    assert np.allclose(data["LigandA"]["metric"]["time"], [0.0, 5.0, 10.0])


def test_md_time_unit_conversion_and_metadata(tmp_path: Path):
    data_dir = tmp_path / "data/sys"
    data_dir.mkdir(parents=True)
    (data_dir / "metric.dat").write_text("0 1.0\n1000 2.0\n2000 3.0\n")

    config = MDConfig(
        base_dir=tmp_path / "data",
        output_dir=tmp_path / "out",
        protein_name="ProteinX",
        systems=[SystemConfig(name="LigandA", dir_pattern="sys", replicates=1)],
        metrics=[
            AnalysisMetric(
                name="metric",
                file_pattern="metric.dat",
                title="Metric",
                ylabel="value",
                data_format="two_column",
            )
        ],
        time_unit_input="ps",
        time_unit_output="ns",
        outlier_threshold=999.0,
    )

    data = MDDataProcessor(config).load_all_data()
    metric = data["LigandA"]["metric"]
    assert np.allclose(metric["time"], [0.0, 1.0, 2.0])
    assert metric["time_metadata"]["output_unit"] == "ns"
    assert metric["time_metadata"]["inferred_input_unit"] == "ps"
    assert metric["time_metadata"]["applied_scale"] == 0.001


def test_qc_stage_generates_report(tmp_path: Path):
    top = tmp_path / "sys.pdb"
    top.write_text(
        "ATOM      1  N   ALA A   1      10.0  10.0  10.0  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1      10.5  10.5  10.5  1.00  0.00           C\n"
    )
    ndx = tmp_path / "index.ndx"
    ndx.write_text("[ Protein ]\n1 2\n")
    traj1 = tmp_path / "rep1.xvg"
    traj2 = tmp_path / "rep2.xvg"
    traj1.write_text("@ xaxis label \"Time\"\n0 1.0\n100 2.0\n200 3.0\n")
    traj2.write_text("0 1.5\n100 2.5\n200 3.5\n")

    pipeline_yaml = tmp_path / "pipeline.yaml"
    pipeline_yaml.write_text(
        yaml.safe_dump(
            {
                "project": {"name": "qc_demo", "protein_name": "ProteinX"},
                "run": ["qc"],
                "qc": {
                    "enabled": True,
                    "base_dir": str(tmp_path),
                    "output_dir": str(tmp_path / "qc_out"),
                    "topology": str(top),
                    "trajectories": [str(traj1), str(traj2)],
                    "index_file": str(ndx),
                    "expected_groups": ["Protein", "Ligand"],
                    "expected_replicates": 3,
                },
            },
            sort_keys=False,
        )
    )

    cfg = load_pipeline_config(pipeline_yaml)
    results = run_pipeline(cfg)

    assert "qc" in results
    assert results["qc"]["success"] is True
    report_file = Path(results["qc"]["report_file"])
    assert report_file.exists()
    report = yaml.safe_load(report_file.read_text())
    assert report["checks"]["topology"]["atom_count"] == 2
    assert report["checks"]["replicates"]["missing_replicates"] == 1
    assert report["checks"]["trajectories"][0]["frame_count"] == 3
    assert report["checks"]["trajectories"][0]["timestep"] == 100.0
    assert report["checks"]["index"]["missing_expected_groups"] == ["Ligand"]


def test_plotly_export_defaults_use_pio_defaults(tmp_path: Path):
    md_cfg = MDConfig(
        base_dir=tmp_path,
        output_dir=tmp_path / "out",
        protein_name="ProteinX",
        systems=[SystemConfig(name="LigandA", dir_pattern="sys", replicates=1)],
    )
    MDPlotter(md_cfg)
    assert pio.defaults.default_format == "svg"
    assert pio.defaults.default_width == md_cfg.plot_config.width
    assert pio.defaults.default_height == md_cfg.plot_config.height
    assert pio.defaults.default_scale == md_cfg.plot_config.scale

    pca_cfg = PCAConfig(base_dir=tmp_path / "pca", output_dir=tmp_path / "pca_out")
    PCAPlotter(pca_cfg)
    assert pio.defaults.default_format == "svg"


def test_prolif_replicate_occupancy_summary_stats(tmp_path: Path):
    analyzer = ProlifAnalyzer(
        ProlifConfig(
            base_dir=tmp_path,
            output_dir=tmp_path / "prolif",
            protein_name="ProteinX",
            ligand_name="LigandA",
            topology=tmp_path / "top.pdb",
            trajectories=[tmp_path / "traj1.xtc", tmp_path / "traj2.xtc"],
        )
    )

    replicate_df = pd.DataFrame(
        [
            {"Interaction": "Hydrophobic", "Residue": "ALA123", "Occupancy": 20.0, "Replicate": 1, "Trajectory": "a"},
            {"Interaction": "Hydrophobic", "Residue": "ALA123", "Occupancy": 40.0, "Replicate": 2, "Trajectory": "b"},
            {"Interaction": "HBAcceptor", "Residue": "SER88", "Occupancy": 10.0, "Replicate": 1, "Trajectory": "a"},
        ]
    )
    summary = analyzer._summarize_replicate_occupancy(replicate_df)

    hydro = summary[(summary["Interaction"] == "Hydrophobic") & (summary["Residue"] == "ALA123")].iloc[0]
    assert hydro["OccupancyMean"] == 30.0
    assert round(float(hydro["OccupancySD"]), 6) == round(np.std([20.0, 40.0], ddof=1), 6)
    assert hydro["ReplicateCount"] == 2


def test_prolif_compare_table_accepts_summary_format(tmp_path: Path):
    system_a = ProlifConfig(
        base_dir=tmp_path,
        output_dir=tmp_path / "A",
        protein_name="ProteinX",
        ligand_name="LigA",
        topology=tmp_path / "top.pdb",
        trajectories=[tmp_path / "traj1.xtc"],
    )
    system_b = ProlifConfig(
        base_dir=tmp_path,
        output_dir=tmp_path / "B",
        protein_name="ProteinX",
        ligand_name="LigB",
        topology=tmp_path / "top.pdb",
        trajectories=[tmp_path / "traj2.xtc"],
    )
    cmp = ProlifComparisonAnalyzer(
        ProlifCompareConfig(
            output_dir=tmp_path / "cmp",
            protein_name="ProteinX",
            systems=[system_a, system_b],
            compare_mode="union",
        )
    )

    systems_data = {
        "LigA": pd.DataFrame(
            [{"Interaction": "Hydrophobic", "Residue": "ALA123", "OccupancyMean": 45.0, "OccupancySD": 2.0, "ReplicateCount": 3}]
        ),
        "LigB": pd.DataFrame(
            [{"Interaction": "Hydrophobic", "Residue": "ALA123", "OccupancyMean": 10.0, "OccupancySD": 1.0, "ReplicateCount": 3}]
        ),
    }
    compare_table = cmp._build_compare_table(systems_data)
    assert set(compare_table.columns) >= {"Occupancy", "OccupancySD", "ReplicateCount"}
    assert len(compare_table) == 2


def test_prolif_visualization_defaults_loaded_from_yaml(tmp_path: Path):
    config_path = tmp_path / "pipeline.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "project": {"name": "demo", "protein_name": "ProteinX"},
                "run": ["prolif"],
                "prolif": {
                    "enabled": True,
                    "base_dir": "./data",
                    "output_dir": "./out/prolif",
                    "topology": "./data/sys/top.psf",
                    "trajectories": ["./data/sys/traj.xtc"],
                },
            },
            sort_keys=False,
        )
    )

    cfg = load_pipeline_config(config_path)
    prolif_cfg = cfg.stages["prolif"].config

    assert prolif_cfg.lignetwork_display_all is True
    assert prolif_cfg.lignetwork_count_aware is True
    assert prolif_cfg.barcode_n_frame_ticks == 10
    assert prolif_cfg.barcode_residues_tick_location == "top"
    assert prolif_cfg.barcode_xlabel == "Frame"
    assert prolif_cfg.barcode_figsize == (12.0, 4.0)
    assert prolif_cfg.barcode_dpi == 150
    assert prolif_cfg.barcode_only_interacting_residues is True
    assert prolif_cfg.barcode_min_interaction_frames == 1
    assert prolif_cfg.barcode_max_residues is None


def test_prolif_uses_count_aware_lignetwork_and_barcode_controls(tmp_path: Path):
    config = ProlifConfig(
        base_dir=tmp_path,
        output_dir=tmp_path / "prolif",
        protein_name="ProteinX",
        ligand_name="LigandA",
        topology=tmp_path / "top.pdb",
        trajectories=[tmp_path / "traj1.xtc"],
        outputs=["barcode", "lignetwork", "occupancy"],
        lignetwork_display_all=True,
        lignetwork_count_aware=True,
        barcode_n_frame_ticks=7,
        barcode_residues_tick_location="bottom",
        barcode_xlabel="Time (frame)",
        barcode_figsize=(9.0, 3.0),
        barcode_dpi=220,
    )
    analyzer = ProlifAnalyzer(config)

    class DummyFigure:
        def to_html(self):
            return "<html></html>"

        def savefig(self, path, format=None, bbox_inches=None):
            Path(path).write_text("dummy")

    class DummyFingerprint:
        def __init__(self, label: str, calls: dict):
            self.label = label
            self.calls = calls

        def to_dataframe(self):
            return pd.DataFrame(
                {
                    ("LIG1", "ALA123", "Hydrophobic"): [1, 0, 1],
                }
            )

        def to_pickle(self, path: str):
            Path(path).write_text("pickle")

        def plot_barcode(self, **kwargs):
            self.calls.setdefault("barcode_kwargs", []).append(kwargs)
            return DummyFigure()

        def plot_lignetwork(self, ligand_mol, **kwargs):
            self.calls["lignetwork_fp"] = self.label
            self.calls["lignetwork_kwargs"] = kwargs
            self.calls["lignetwork_ligand"] = ligand_mol
            return DummyFigure()

    calls = {}
    fp = DummyFingerprint("default_fp", calls)
    fp_count = DummyFingerprint("count_fp", calls)

    analyzer.processor.run = lambda: {
        "replicates": [
            {
                "replicate": 1,
                "trajectory": "traj1.xtc",
                "fingerprint": fp,
                "fingerprint_count": fp_count,
                "ligand_mol": "ligand",
                "protein_mol": "protein",
                "chemistry": {},
            }
        ],
        "chemistry_requirements": {},
        "chemistry_summary": {},
    }

    result = analyzer.run_analysis()
    assert result["success"] is True
    assert calls["lignetwork_fp"] == "count_fp"
    assert calls["lignetwork_kwargs"]["kind"] == "aggregate"
    assert calls["lignetwork_kwargs"]["threshold"] == 0.3
    assert calls["lignetwork_kwargs"]["display_all"] is True
    assert calls["barcode_kwargs"][0] == {
        "n_frame_ticks": 7,
        "residues_tick_location": "bottom",
        "xlabel": "Time (frame)",
        "figsize": (9.0, 3.0),
        "dpi": 220,
    }
    plots_dir = config.output_dir / "plots"
    data_dir = config.output_dir / "data"
    assert (plots_dir / "barcode_rep1.svg").exists()
    assert (plots_dir / "barcode_rep1.html").exists()
    assert (plots_dir / "barcode_rep1_interactive.html").exists()
    qc = pd.read_csv(data_dir / "barcode_qc_rep1.csv")
    assert list(qc["Residue"]) == ["ALA123"]
    assert int(qc.loc[0, "InteractingFrames"]) == 2


def test_ttclust_batch_config_paths_and_defaults(tmp_path: Path):
    config_path = tmp_path / "pipeline.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "project": {"name": "tt_demo", "protein_name": "KEAP1"},
                "paths": {"data_root": "./data", "output_root": "./out"},
                "run": ["ttclust_batch"],
                "ttclust_batch": {
                    "enabled": True,
                    "base_dir": "./data",
                    "output_dir": "./out/ttclust_batch",
                    "systems": [
                        {
                            "name": "Benzyl",
                            "topology": "./data/sys/top.psf",
                            "trajectories": ["./data/sys/r1.xtc", "sys/r2.xtc"],
                        }
                    ],
                },
            },
            sort_keys=False,
        )
    )

    cfg = load_pipeline_config(config_path)
    tt_cfg = cfg.stages["ttclust_batch"].config[0]
    assert tt_cfg.topology == tmp_path / "data/sys/top.psf"
    assert tt_cfg.trajectories[0] == tmp_path / "data/sys/r1.xtc"
    assert tt_cfg.trajectories[1] == tmp_path / "data/sys/r2.xtc"
    assert tt_cfg.method == "ward"
    assert tt_cfg.autoclust is True
    assert tt_cfg.interactive_matrix is False


def test_ttclust_analyzer_runs_and_collects_figures(tmp_path: Path):
    data_dir = tmp_path / "data" / "sys"
    data_dir.mkdir(parents=True)
    (data_dir / "top.psf").write_text("PSF")
    (data_dir / "r1.xtc").write_text("XTC")

    fake_ttclust = tmp_path / "fake_ttclust.py"
    fake_ttclust.write_text(
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import pathlib",
                "import sys",
                "args = sys.argv[1:]",
                "logfile = 'clustering.log'",
                "if '-l' in args and args.index('-l') + 1 < len(args):",
                "    logfile = args[args.index('-l') + 1]",
                "pathlib.Path(logfile).write_text('ok\\n')",
                "folder = pathlib.Path(pathlib.Path(logfile).stem)",
                "folder.mkdir(exist_ok=True)",
                "for name in ['dendrogram.png', 'linear_projection.png', 'barplot.png', 'distance_projection.png', 'distance_matrix.png']:",
                "    (folder / name).write_text('png')",
                "cluster_dir = pathlib.Path('Cluster_PDB')",
                "cluster_dir.mkdir(exist_ok=True)",
                "(cluster_dir / 'C1-f1-s1.pdb').write_text('MODEL\\nEND\\n')",
                "pathlib.Path('argv.txt').write_text('\\n'.join(args) + '\\n')",
                "sys.exit(0)",
            ]
        )
        + "\n"
    )
    fake_ttclust.chmod(0o755)

    cfg = TTClustConfig(
        base_dir=tmp_path / "data",
        output_dir=tmp_path / "out" / "ttclust",
        protein_name="KEAP1",
        system_name="Benzyl",
        topology=tmp_path / "data" / "sys" / "top.psf",
        trajectories=[tmp_path / "data" / "sys" / "r1.xtc"],
        executable=str(fake_ttclust),
        stride=10,
        logfile="benzyl.log",
        select_traj="all",
        select_alignment="backbone",
        select_rmsd="backbone",
        method="ward",
        autoclust=True,
        interactive_matrix=False,
    )

    result = TTClustAnalyzer(cfg).run_analysis()

    assert result["success"] is True
    assert result["figure_count"] == 5
    assert Path(result["log_file"]).exists()
    assert Path(result["cluster_pdb_dir"]).exists()
    assert Path(result["figures_index_file"]).exists()
    command_line = Path(result["command_file"]).read_text()
    assert "-f" in command_line
    assert "-t" in command_line
    assert "-i n" in command_line
    assert "-aa Y" in command_line


def test_prolif_plotter_saves_barcode_axes_and_lignetwork_like_objects(tmp_path: Path):
    cfg = ProlifConfig(
        base_dir=tmp_path,
        output_dir=tmp_path / "prolif",
        protein_name="ProteinX",
        ligand_name="LigandA",
        topology=tmp_path / "top.pdb",
        trajectories=[tmp_path / "traj.xtc"],
    )
    plotter = ProlifPlotter(cfg.plot_config, cfg.protein_name)
    out_dir = tmp_path / "plots"

    class DummyFigure:
        def savefig(self, path: str, format=None, bbox_inches=None):
            Path(path).write_text("svg")

    class DummyAxes:
        def __init__(self):
            self.figure = DummyFigure()

    plotter.save_figure(DummyAxes(), out_dir, "barcode_axes", formats=["svg"])
    assert (out_dir / "barcode_axes.svg").exists()

    class DummyLigNetwork:
        def __init__(self):
            self.saved = []

        def save(self, path: str):
            Path(path).write_text("<html></html>")
            self.saved.append(path)

    ln = DummyLigNetwork()
    plotter.save_figure(ln, out_dir, "lignetwork", formats=["html"])
    assert (out_dir / "lignetwork.html").exists()
    assert ln.saved


def test_pca_preflight_fails_on_missing_prepared_files(tmp_path: Path):
    cfg = PCAConfig(base_dir=tmp_path / "pca", output_dir=tmp_path / "out")
    runner = PCAGromacsRunner(cfg)
    runner.pca_dir.mkdir(parents=True, exist_ok=True)
    (runner.pca_dir / "ref.pdb").write_text(
        "ATOM      1  CA  ALA A   1      10.0  10.0  10.0  1.00  0.00           C\n"
    )

    with pytest.raises(RuntimeError):
        runner._preflight_prepared_trajectory()

    preflight = runner.pca_dir / "pca_preflight.json"
    assert preflight.exists()


def test_utils_duivy_inspection_guide_generation(tmp_path: Path):
    out_dir = tmp_path / "utils"
    result = run_utils_tasks(
        {
            "output_dir": str(out_dir),
            "tasks": [
                {
                    "type": "duivy_inspect_guide",
                    "xvg_files": ["/tmp/a.xvg", "/tmp/b.xvg"],
                    "xpm_files": ["/tmp/a.xpm"],
                }
            ],
        }
    )
    assert result["success"] is True
    task_result = result["tasks"][0]
    assert task_result["success"] is True
    guide = Path(task_result["guide_file"])
    assert guide.exists()
    text = guide.read_text()
    assert "dit xvg_show -f /tmp/a.xvg" in text
    assert "does not replace the maintained Python toolkit outputs" in text


def test_ranking_consumes_mmpbsa_qc_metadata(tmp_path: Path):
    md_stats = pd.DataFrame(
        [
            {"System": "LigA", "Metric": "rmsd_prot", "Mean": 1.0},
            {"System": "LigB", "Metric": "rmsd_prot", "Mean": 1.0},
        ]
    )
    md_stats_file = tmp_path / "md_stats.csv"
    md_stats.to_csv(md_stats_file, index=False)

    mmpbsa_dir = tmp_path / "mmpbsa_data"
    mmpbsa_dir.mkdir(parents=True)
    pd.DataFrame(
        [
            {
                "Component": "Binding_TOTAL",
                "Mean": -20.0,
                "QC_ExpectedReplicates": 3,
                "QC_FoundResultsReplicates": 3,
                "QC_MissingResultsReplicates": 0,
                "QC_ReplicateComplete": True,
                "QC_GroupIdentityOK": True,
            }
        ]
    ).to_csv(mmpbsa_dir / "LigA_results.csv", index=False)
    pd.DataFrame(
        [
            {"Component": "Binding_TOTAL", "Mean": -20.0},
        ]
    ).to_csv(mmpbsa_dir / "LigB_results.csv", index=False)
    pd.DataFrame(
        [
            {
                "expected_replicates": 3,
                "found_results_replicates": 2,
                "missing_results_replicates": 1,
                "replicate_complete": False,
                "group_identity_ok": False,
            }
        ]
    ).to_csv(mmpbsa_dir / "LigB_qc.csv", index=False)

    analyzer = RankingAnalyzer(
        RankingConfig(
            output_dir=tmp_path / "ranking",
            protein_name="ProteinX",
            md_stats_file=md_stats_file,
            mmpbsa_data_dir=mmpbsa_dir,
        )
    )
    result = analyzer.run_analysis()
    df = pd.read_csv(result["ranking_file"])

    assert "qc_penalty" in df.columns
    row_a = df[df["Ligand"] == "LigA"].iloc[0]
    row_b = df[df["Ligand"] == "LigB"].iloc[0]
    assert float(row_a["qc_penalty"]) == 0.0
    assert float(row_b["qc_penalty"]) > 0.0
    assert row_a["qc_status"] == "pass"
    assert row_b["qc_status"] == "warn"


def test_upstream_script_generation_is_write_only(tmp_path: Path):
    result = generate_upstream_scripts(
        {
            "outputs": ["md_vmd", "mmpbsa", "pca"],
            "topology": "top.psf",
            "trajectory": "traj.xtc",
            "groups": [1, 2],
        },
        tmp_path,
    )

    script_paths = [Path(path) for path in result["scripts"]]
    assert {path.name for path in script_paths} == {
        "generate_md_metrics.tcl",
        "run_mmpbsa.sh",
        "prepare_pca.sh",
    }
    md_script = (tmp_path / "generate_md_metrics.tcl").read_text()
    assert "RMSD_complex_.dat" in md_script
    assert "RMSF_.dat" in md_script
    assert "SASA_.dat" in md_script
    assert "comcom_.dat" in md_script
    assert "gmx_MMPBSA" in (tmp_path / "run_mmpbsa.sh").read_text()
    assert "gmx covar" in (tmp_path / "prepare_pca.sh").read_text()


def test_sequence_labels_from_pdb_are_extracted_as_strings(tmp_path: Path):
    pdb_path = tmp_path / "step3_input.pdb"
    pdb_path.write_text(
        "\n".join(
            [
                "ATOM      1  N   ALA A 321      10.000  10.000  10.000  1.00  0.00           N",
                "ATOM      2  CA  ALA A 321      11.000  10.000  10.000  1.00  0.00           C",
                "ATOM      3  N   VAL A 322      12.000  10.000  10.000  1.00  0.00           N",
                "ATOM      4  CA  VAL A 322      13.000  10.000  10.000  1.00  0.00           C",
            ]
        )
        + "\n"
    )

    residues = load_sequence_from_pdb(pdb_path)
    labels = build_residue_labels(residues, include_names=True)
    assert labels == ["ALA321", "VAL322"]


def test_md_rmsf_axis_auto_detects_step3_pdb_labels(tmp_path: Path):
    system_dir = tmp_path / "data" / "sys_1"
    system_dir.mkdir(parents=True)
    (system_dir / "step3_input.pdb").write_text(
        "\n".join(
            [
                "ATOM      1  N   GLY A 100      10.000  10.000  10.000  1.00  0.00           N",
                "ATOM      2  CA  GLY A 100      11.000  10.000  10.000  1.00  0.00           C",
                "ATOM      3  N   LEU A 101      12.000  10.000  10.000  1.00  0.00           N",
                "ATOM      4  CA  LEU A 101      13.000  10.000  10.000  1.00  0.00           C",
            ]
        )
        + "\n"
    )

    config = MDConfig(
        base_dir=tmp_path / "data",
        output_dir=tmp_path / "out",
        protein_name="ProteinX",
        systems=[SystemConfig(name="LigandA", dir_names=["sys_1"])],
        metrics=[
            AnalysisMetric(
                name="rmsf",
                file_pattern="RMSF_.dat",
                title="RMSF",
                ylabel="A",
                data_format="two_column",
            )
        ],
    )

    labels = MDDataProcessor(config).prepare_amino_acid_axis()
    assert labels == ["GLY100", "LEU101"]


def test_md_run_exports_metric_csv_to_data_folder(tmp_path: Path, monkeypatch):
    config = MDConfig(
        base_dir=tmp_path / "data",
        output_dir=tmp_path / "out",
        protein_name="ProteinX",
        systems=[SystemConfig(name="LigandA", dir_pattern="sys", replicates=1)],
        metrics=[
            AnalysisMetric(
                name="rmsd_prot",
                file_pattern="RMSD_apo_.dat",
                title="Protein RMSD",
                ylabel="A",
                data_format="single_column",
            )
        ],
    )

    analyzer = MDAnalyzer(config)
    fake_all_data = {
        "LigandA": {
            "rmsd_prot": {
                "time": np.array([0.0, 1.0, 2.0]),
                "values": np.array([1.0, 1.5, 2.0]),
                "mode": "single",
                "stats": {"mean": 1.5, "std": 0.5, "min": 1.0, "max": 2.0},
                "n_points": 3,
            }
        }
    }

    monkeypatch.setattr(analyzer.processor, "load_all_data", lambda: fake_all_data)
    monkeypatch.setattr(analyzer.plotter, "create_all_plots", lambda data: [])

    result = analyzer.run_analysis()
    assert result["success"] is True
    exported = config.output_dir / "data" / "LigandA_rmsd_prot.csv"
    assert exported.exists()
    df = pd.read_csv(exported)
    assert list(df.columns) == ["x", "value"]
    assert len(df) == 3
