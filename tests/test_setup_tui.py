"""Headless pilot tests of the setup interface (specs/002-pipeline-setup-tui, US4)."""

import asyncio
from pathlib import Path

import pytest
import yaml

pytest.importorskip("textual")

from textual.widgets import Checkbox, DataTable, Input, Select, Static, TabbedContent  # noqa: E402

from gromacs_analysis.setup.backends import LocalBackend  # noqa: E402
from gromacs_analysis.setup.generator import refresh_stages  # noqa: E402
from gromacs_analysis.setup.project import Project  # noqa: E402
from gromacs_analysis.setup.scanner import scan  # noqa: E402
from gromacs_analysis.setup.tui import ClaritySetupApp, ConnectScreen, EditNameScreen  # noqa: E402
from setup_fixtures import build_project  # noqa: E402

SIZE = (160, 45)


def make_app(tmp_path: Path) -> ClaritySetupApp:
    root = build_project(tmp_path / "proj")
    project = Project.from_scan(scan(LocalBackend(), str(root)))
    refresh_stages(project)
    return ClaritySetupApp(project=project, project_path=tmp_path / "clarity_project.yaml")


def run(coro):
    return asyncio.run(coro)


def cursor_to(table: DataTable, key: str) -> None:
    table.move_cursor(row=table.get_row_index(key))


async def rename_selected(app, pilot, new_name: str) -> None:
    await pilot.press("e")
    await pilot.pause()
    assert isinstance(app.screen, EditNameScreen)
    field = app.screen.query_one(Input)
    field.value = new_name
    await pilot.press("enter")
    await pilot.pause()


def test_tabs_and_tables_are_filled(tmp_path: Path):
    async def go():
        app = make_app(tmp_path)
        async with app.run_test(size=SIZE) as pilot:
            await pilot.pause()
            assert {p.id for p in app.query("TabPane")} == {"tab-systems", "tab-ligands", "tab-stages",
                                                           "tab-preview"}
            assert app.query_one("#systems", DataTable).row_count == 4
            assert app.query_one("#ligands", DataTable).row_count == 4  # INH, NAD, LGB, LGC
            assert len(app.query(Checkbox)) == 9
            assert "prolif_batch" in str(app.query_one("#preview", Static).render())
    run(go())


def test_rename_system_apo_and_exclude(tmp_path: Path):
    async def go():
        app = make_app(tmp_path)
        async with app.run_test(size=SIZE) as pilot:
            await pilot.pause()
            table = app.query_one("#systems", DataTable)
            cursor_to(table, "HOLO_LigA")
            await rename_selected(app, pilot, "Isoniazid complex")
            assert app.project.system("HOLO_LigA").name == "Isoniazid_complex"

            cursor_to(table, "LigC")
            await pilot.press("a")
            await pilot.press("x")
            await pilot.pause()
            assert app.project.system("LigC").apo and not app.project.system("LigC").included
            assert table.get_row("HOLO_LigA")[0] == "Isoniazid_complex"
            assert table.get_row("LigC")[0].endswith("(excluded)") and table.get_row("LigC")[3] == "APO ref"
    run(go())


def test_duplicate_name_shows_error_and_keeps_old_name(tmp_path: Path):
    async def go():
        app = make_app(tmp_path)
        notes = []
        app.notify = lambda message, **kw: notes.append((message, kw.get("severity")))
        async with app.run_test(size=SIZE) as pilot:
            await pilot.pause()
            old = app.project.system("HOLO_LigB").name
            cursor_to(app.query_one("#systems", DataTable), "HOLO_LigB")
            await rename_selected(app, pilot, app.project.system("HOLO_LigA").name)
            assert app.project.system("HOLO_LigB").name == old
            assert any(sev == "error" and "already used" in msg for msg, sev in notes)
    run(go())


def test_ligand_rename_and_exclude(tmp_path: Path):
    async def go():
        app = make_app(tmp_path)
        async with app.run_test(size=SIZE) as pilot:
            await pilot.pause()
            app.query_one("#tabs", TabbedContent).active = "tab-ligands"
            await pilot.pause()
            table = app.query_one("#ligands", DataTable)
            table.focus()
            cursor_to(table, "HOLO_LigA\tINH")
            await rename_selected(app, pilot, "Isoniazid")
            lig = next(l for l in app.project.system("HOLO_LigA").ligands if l.resname == "INH")
            assert lig.name == "Isoniazid" and lig.selection == "resname INH"
            cursor_to(table, "HOLO_LigA\tNAD")
            await pilot.press("x")
            await pilot.pause()
            nad = next(l for l in app.project.system("HOLO_LigA").ligands if l.resname == "NAD")
            assert not nad.included
            assert "Isoniazid" in str(app.query_one("#preview", Static).render())
    run(go())


def test_stage_toggle_updates_preview_and_save_writes_files(tmp_path: Path):
    async def go():
        app = make_app(tmp_path)
        async with app.run_test(size=SIZE) as pilot:
            await pilot.pause()
            app.query_one("#tabs", TabbedContent).active = "tab-stages"
            await pilot.pause()
            box = app.query_one("#stage-pca_joint", Checkbox)
            box.value = False
            await pilot.pause()
            assert not app.project.stages["pca_joint"].enabled
            assert "pca_joint:" not in str(app.query_one("#preview", Static).render())
            await pilot.press("s")
            await pilot.pause()
        data = yaml.safe_load((tmp_path / "pipeline.yaml").read_text())
        assert "pca_joint" not in data["run"] and "md" in data["run"]
        assert Project.load(tmp_path / "clarity_project.yaml").stages["pca_joint"].enabled is False
    run(go())


def test_connect_screen_scans_local_folder(tmp_path: Path):
    root = build_project(tmp_path / "proj")

    async def go():
        app = ClaritySetupApp(project_path=tmp_path / "clarity_project.yaml")
        async with app.run_test(size=SIZE) as pilot:
            await pilot.pause()
            assert isinstance(app.screen, ConnectScreen)
            assert app.screen.query_one("#kind", Select).value == "local"
            app.screen.query_one("#root", Input).value = str(root)
            await pilot.click("#scan")
            await app.workers.wait_for_complete()
            await pilot.pause()
            assert app.project is not None and len(app.project.systems) == 4
            assert app.query_one("#systems", DataTable).row_count == 4
    run(go())


def test_submit_without_remote_connection_is_an_error_not_a_crash(tmp_path: Path):
    async def go():
        app = make_app(tmp_path)
        notes = []
        app.notify = lambda message, **kw: notes.append((message, kw.get("severity")))
        async with app.run_test(size=SIZE) as pilot:
            await pilot.pause()
            await pilot.press("c")
            await pilot.pause()
            assert any(sev == "error" and "HTCondor" in msg for msg, sev in notes)
            assert app.is_running
    run(go())
