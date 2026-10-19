"""End-to-end tests on real trajectory files written to disk (US10).

QC, distance analysis and joint PCA run on a generated two-chain protein-ligand
system through the real MDAnalysis readers. ProLIF runs on ProLIF's bundled
protein-ligand test trajectory (it needs explicit hydrogens and chemistry that a
synthetic system cannot provide) and is skipped when ProLIF is not installed.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("MDAnalysis")

from trajectory_factory import build_system  # noqa: E402


@pytest.fixture
def drifting_system(tmp_path: Path):
    return build_system(tmp_path / "sys", n_frames=20, drift_per_frame=0.5, noise=0.05, seed=1)


def test_qc_reads_trajectory_metadata(tmp_path: Path, drifting_system):
    from gromacs_analysis.qc import QCAnalyzer, QCConfig

    topology, trajectory = drifting_system
    result = QCAnalyzer(QCConfig(
        base_dir=tmp_path, output_dir=tmp_path / "qc", topology=topology, trajectories=[trajectory],
        ligand_selection="resname LIG", sample_stride=1, dissociation_cutoff=8.0,
    )).run_analysis()

    import json
    report = json.loads(Path(result["report_file"]).read_text())
    entry = report["checks"]["trajectories"][0]
    assert entry["frame_count"] == 20
    assert entry["timestep"] == pytest.approx(10.0)
    assert entry["uniform_sampling"] is True
    assert entry["ligand_selection_count"] == 3
    # Ligand starts ~3 Å above the protein plane and drifts 9.5 Å away.
    assert entry["min_ligand_protein_distance"] < entry["max_ligand_protein_min_distance"]
    assert any("loses contact" in w for w in result["warnings"])


def test_distance_increases_as_ligand_drifts(tmp_path: Path, drifting_system):
    from gromacs_analysis.distance import DistanceAnalyzer
    from gromacs_analysis.distance.config import DistanceConfig

    topology, trajectory = drifting_system
    cfg = DistanceConfig(
        base_dir=tmp_path, output_dir=tmp_path / "dist", protein_name="P", ligand_name="LIG",
        topology=topology, trajectories=[trajectory], ligand_selection="resname LIG",
        residue_ids=["A:4", "B:4"], method="min", time_unit="ps",
    )
    cfg.plot_config.save_formats = ["html"]
    result = DistanceAnalyzer(cfg).run_analysis()
    assert result["success"]

    csvs = sorted((tmp_path / "dist" / "data").glob("distance_*.csv"))
    assert {p.stem for p in csvs} >= {"distance_A_GLY4", "distance_B_GLY4"}
    series = pd.read_csv(tmp_path / "dist" / "data" / "distance_A_GLY4.csv")["mean_distance"].to_numpy()
    assert np.all(np.diff(series) > 0)  # monotonic drift as constructed
    assert series[-1] - series[0] > 5.0


def test_joint_pca_on_real_trajectories(tmp_path: Path):
    from gromacs_analysis.pca import PCAJointAnalyzer, PCAJointConfig, PCAJointSystem

    apo_top, apo_traj = build_system(tmp_path / "apo", n_frames=15, noise=0.2, seed=2, name="apo")
    lig_top, lig_traj = build_system(tmp_path / "lig", n_frames=15, noise=0.2, seed=3, name="lig")
    cfg = PCAJointConfig(
        output_dir=tmp_path / "joint", protein_name="P",
        systems=[
            PCAJointSystem(name="APO", topology=apo_top, trajectories=[apo_traj], is_apo=True),
            PCAJointSystem(name="LigA", topology=lig_top, trajectories=[lig_traj]),
        ],
        n_components=5, n_overlap=3,
    )
    cfg.plot_config.save_formats = ["html"]
    result = PCAJointAnalyzer(cfg).run_analysis()

    params = pd.read_json(tmp_path / "joint" / "data" / "parameters.json", typ="series")
    assert result["success"]
    assert params["n_atoms"] == 16  # two chains x 8 CA
    assert (tmp_path / "joint" / "plots" / "joint_projection_12.html").exists()


def test_prolif_on_bundled_trajectory(tmp_path: Path):
    pytest.importorskip("prolif")
    pytest.importorskip("rdkit")
    from prolif.datafiles import TOP, TRAJ

    from gromacs_analysis.prolif import ProlifAnalyzer, ProlifConfig

    cfg = ProlifConfig(
        base_dir=tmp_path, output_dir=tmp_path / "prolif", protein_name="P", ligand_name="LIG",
        topology=Path(TOP), trajectories=[Path(TRAJ)], ligand_selection="resname LIG",
        protein_selection="protein and byres around 6.0 resname LIG", stride=10, n_jobs=1,
        outputs=["occupancy", "barcode"], save_pickle=False,
    )
    cfg.plot_config.save_formats = ["html"]
    result = ProlifAnalyzer(cfg).run_analysis()

    summary = pd.read_csv(result["occupancy_summary_file"])
    assert "index" not in summary.columns
    assert not summary.empty
    assert summary["Residue"].str.match(r"^[A-Z]{3}\d+").all()  # protein residues, not the ligand
    assert (tmp_path / "prolif" / "plots" / "occupancy_heatmap.html").exists()
