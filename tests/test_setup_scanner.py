"""Discovery of systems, replicates and ligands (specs/002-pipeline-setup-tui, US1)."""

from pathlib import Path

import pytest

from gromacs_analysis.setup.backends import LocalBackend, ShellRunner
from gromacs_analysis.setup.scanner import scan, split_replicate
from setup_fixtures import build_project


@pytest.mark.parametrize("backend_cls", [LocalBackend, ShellRunner])
def test_scan_groups_systems_and_replicates(tmp_path: Path, backend_cls):
    root = build_project(tmp_path / "proj")
    result = scan(backend_cls(), str(root))
    systems = {s.key: s for s in result.systems}

    assert set(systems) == {"HOLO_LigA", "HOLO_LigB", "APO", "LigC"}
    assert [r.number for r in systems["HOLO_LigA"].replicates] == [1, 2, 3]
    assert [r.dir for r in systems["HOLO_LigB"].replicates] == ["HOLO_LigB_rep1", "HOLO_LigB_rep2"]
    assert [r.dir for r in systems["LigC"].replicates] == ["LigC/rep1", "LigC/rep2"]
    assert systems["APO"].apo and not systems["HOLO_LigA"].apo
    # Folders with simulation files but no trajectory are reported; unrelated folders are not listed.
    assert result.skipped_dirs == ["topology_only"]


def test_scan_detects_ligands_and_filters_shared_and_solvent(tmp_path: Path):
    systems = {s.key: s for s in scan(LocalBackend(), str(build_project(tmp_path / "p"))).systems}
    assert systems["HOLO_LigA"].ligands == ["INH", "NAD"]
    assert systems["HOLO_LigB"].ligands == ["LGB"]  # HEM is also in APO -> cofactor, not ligand
    assert systems["LigC"].ligands == ["LGC"]
    assert systems["APO"].ligands == []


def test_scan_picks_topology_trajectory_and_inputs(tmp_path: Path):
    systems = {s.key: s for s in scan(LocalBackend(), str(build_project(tmp_path / "p"))).systems}
    rep = systems["HOLO_LigA"].replicates[0]
    assert rep.topology == "HOLO_LigA_1/step3_input.psf"
    assert rep.trajectories[0] == "HOLO_LigA_1/step5_production_noPBC.xtc"
    assert rep.metric_dir == "HOLO_LigA_1"
    assert systems["HOLO_LigA"].inputs == {
        "md_metrics": True, "mmpbsa": True, "pca_dir": False, "processed_trajectory": True,
    }
    apo_rep = systems["APO"].replicates[0]
    assert apo_rep.metric_dir == "APO_1/analysis"  # metric files in a subfolder are found
    assert systems["LigC"].inputs["md_metrics"] is False
    assert systems["HOLO_LigA"].default_name == "LigA" and systems["APO"].default_name == "APO"


def test_remote_scan_uses_bounded_commands(tmp_path: Path):
    root = build_project(tmp_path / "p")
    backend = ShellRunner()
    result = scan(backend, str(root))
    unique_topologies = 4  # one text topology per system
    assert result.commands_used <= 1 + unique_topologies


@pytest.mark.parametrize("rel,expected", [
    ("HOLO_LigA_1", ("HOLO_LigA", 1)),
    ("sysX-rep2", ("sysX", 2)),
    ("SysY_r3", ("SysY", 3)),
    ("SysY_run4", ("SysY", 4)),
    ("Lig_replica5", ("Lig", 5)),
    ("Lig_R6", ("Lig", 6)),
    ("LigC/rep1", ("LigC", 1)),
    ("LigC/7", ("LigC", 7)),
    ("Cmpd12", ("Cmpd12", 1)),       # no separator: a compound name, not a replicate
    ("single_system", ("single_system", 1)),
])
def test_split_replicate_variants(rel, expected):
    assert split_replicate(rel) == expected


def test_scan_without_systems_reports_skipped(tmp_path: Path):
    (tmp_path / "empty" / "a").mkdir(parents=True)
    (tmp_path / "empty" / "a" / "x.pdb").write_text("END\n")
    result = scan(LocalBackend(), str(tmp_path / "empty"))
    assert result.systems == [] and result.skipped_dirs == ["a"]
