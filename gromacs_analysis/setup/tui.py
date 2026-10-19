"""Keyboard-driven terminal interface for clarity-setup (Textual).

Connect -> review systems and ligands (rename, APO, include) -> pick stages ->
preview the generated YAML -> save / run / submit. Every action goes through the
same :class:`Project`, generator and runner used by the CLI; errors become
notifications instead of crashing the interface.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Optional

from textual import on
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.screen import ModalScreen
from textual.widgets import (
    Button, Checkbox, DataTable, Footer, Header, Input, Label, Select, Static, TabbedContent, TabPane,
)

from .backends import connect
from .generator import STAGE_LABELS, STAGES, refresh_stages, render_yaml, validate_yaml
from .project import PROJECT_FILE, Connection, Project, ProjectError
from .scanner import scan

PIPELINE_FILE = "pipeline.yaml"


# ---------------------------------------------------------------- modal screens

class EditNameScreen(ModalScreen[Optional[str]]):
    """Ask for a new display name; dismisses with the text or None on cancel."""

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def __init__(self, title: str, value: str) -> None:
        super().__init__()
        self._title = title
        self._value = value

    def compose(self) -> ComposeResult:
        with Vertical(id="dialog"):
            yield Label(self._title, id="dialog-title")
            yield Input(value=self._value, id="name-input")
            with Horizontal(id="dialog-buttons"):
                yield Button("Save", variant="primary", id="ok")
                yield Button("Cancel", id="cancel")

    def on_mount(self) -> None:
        self.query_one("#name-input", Input).focus()

    @on(Input.Submitted)
    @on(Button.Pressed, "#ok")
    def _accept(self) -> None:
        self.dismiss(self.query_one("#name-input", Input).value)

    @on(Button.Pressed, "#cancel")
    def action_cancel(self) -> None:
        self.dismiss(None)


class ConnectScreen(ModalScreen[Optional[tuple]]):
    """Where the simulations are: dismisses with (root, Connection) or None."""

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def __init__(self, root: str = "", connection: Optional[Connection] = None) -> None:
        super().__init__()
        self._root = root
        self._connection = connection or Connection()

    def compose(self) -> ComposeResult:
        c = self._connection
        with Vertical(id="dialog"):
            yield Label("Where are your simulations?", id="dialog-title")
            yield Select([("This computer", "local"), ("SSH host", "ssh"), ("condorctl profile", "condorctl")],
                         value=c.kind, allow_blank=False, id="kind")
            yield Input(value=self._root, placeholder="Folder that contains the system folders", id="root")
            yield Input(value=c.host or "", placeholder="SSH host or ~/.ssh/config alias", id="host")
            yield Input(value=c.profile or "", placeholder="condorctl config, e.g. ~/condorctl.yaml",
                        id="profile")
            with Horizontal(id="dialog-buttons"):
                yield Button("Scan", variant="primary", id="scan")
                yield Button("Cancel", id="cancel")

    def on_mount(self) -> None:
        self._sync_fields()
        self.query_one("#root", Input).focus()

    @on(Select.Changed, "#kind")
    def _sync_fields(self) -> None:
        kind = self.query_one("#kind", Select).value
        self.query_one("#host", Input).display = kind == "ssh"
        self.query_one("#profile", Input).display = kind == "condorctl"

    @on(Input.Submitted)
    @on(Button.Pressed, "#scan")
    def _accept(self) -> None:
        root = self.query_one("#root", Input).value.strip()
        if not root:
            self.notify("Enter the folder to scan", severity="error")
            return
        kind = self.query_one("#kind", Select).value
        connection = Connection(
            kind=kind,
            host=self.query_one("#host", Input).value.strip() or None if kind == "ssh" else None,
            profile=self.query_one("#profile", Input).value.strip() or None if kind == "condorctl" else None,
        )
        if kind == "ssh" and not connection.host:
            self.notify("Enter the SSH host", severity="error")
            return
        self.dismiss((root, connection))

    @on(Button.Pressed, "#cancel")
    def action_cancel(self) -> None:
        self.dismiss(None)


# ---------------------------------------------------------------- main app

class ClaritySetupApp(App):
    TITLE = "ClarityDynamics setup"
    CSS = """
    #dialog { width: 70; height: auto; padding: 1 2; border: thick $accent; background: $surface; }
    EditNameScreen, ConnectScreen { align: center middle; }
    #dialog-title { text-style: bold; margin-bottom: 1; }
    #dialog-buttons { height: auto; margin-top: 1; }
    #dialog-buttons Button { margin-right: 2; }
    #tabs { height: 1fr; }
    #status { height: 1; padding: 0 1; background: $panel; color: $text-muted; }
    #stages { padding: 1 1; }
    #stages Checkbox { border: none; padding: 0; height: 1; width: 100%; background: transparent; }
    #stages Checkbox:focus { background: $boost; }
    .hint { color: $text-muted; margin: 0 0 1 0; padding: 0 1; }
    .stage-reason { color: $warning; margin: 0 0 0 4; }
    #preview { padding: 0 1; }
    """
    BINDINGS = [
        Binding("e", "edit", "Rename"),
        Binding("a", "toggle_apo", "APO"),
        Binding("x", "toggle_include", "Include"),
        Binding("s", "save", "Save"),
        Binding("r", "run", "Run"),
        Binding("c", "submit", "Condor"),
        Binding("ctrl+r", "rescan", "Rescan"),
        Binding("o", "connect", "Connect"),
        Binding("q", "quit", "Quit"),
    ]

    def __init__(self, project: Optional[Project] = None, root: Optional[str] = None,
                 connection: Optional[Connection] = None, depth: int = 4,
                 project_path: Path = Path(PROJECT_FILE), pipeline_path: Optional[Path] = None,
                 scan_fn: Callable = scan) -> None:
        super().__init__()
        self.project = project
        self._root = root
        self._connection = connection or (project.connection if project else Connection())
        self.depth = depth
        self.project_path = Path(project_path)
        self.pipeline_path = Path(pipeline_path) if pipeline_path else self.project_path.with_name(PIPELINE_FILE)
        self._scan_fn = scan_fn
        self.last_result: Optional[dict] = None

    # ------------------------------------------------------------ layout

    def compose(self) -> ComposeResult:
        yield Header()
        with TabbedContent(id="tabs"):
            with TabPane("Systems", id="tab-systems"):
                yield Static("e rename · a mark APO reference · x include/exclude · ctrl+r rescan",
                             classes="hint")
                yield DataTable(id="systems", cursor_type="row", zebra_stripes=True)
            with TabPane("Ligands", id="tab-ligands"):
                yield Static("e rename (labels only; selections keep the residue name) · x include/exclude",
                             classes="hint")
                yield DataTable(id="ligands", cursor_type="row", zebra_stripes=True)
            with TabPane("Stages", id="tab-stages"):
                yield VerticalScroll(id="stages")
            with TabPane("Preview", id="tab-preview"):
                with VerticalScroll():
                    yield Static("", id="preview")
        yield Static("", id="status")
        yield Footer()

    def on_mount(self) -> None:
        systems = self.query_one("#systems", DataTable)
        systems.add_columns("Name", "Folder", "Reps", "Role", "Ligands", "Inputs", "Topology")
        ligands = self.query_one("#ligands", DataTable)
        ligands.add_columns("Label", "Residue", "System", "Included", "Selection")
        if self.project is not None:
            refresh_stages(self.project) if not self.project.stages else None
            self.refresh_all()
        elif self._root:
            self.start_scan(self._root, self._connection)
        else:
            self.action_connect()

    # ------------------------------------------------------------ scanning

    def action_connect(self) -> None:
        root = self.project.root if self.project else (self._root or "")
        self.push_screen(ConnectScreen(root, self._connection), self._connected)

    def _connected(self, value: Optional[tuple]) -> None:
        if value:
            root, connection = value
            self.start_scan(root, connection)

    def action_rescan(self) -> None:
        if self.project is None:
            self.action_connect()
        else:
            self.start_scan(self.project.root, self.project.connection)

    def start_scan(self, root: str, connection: Connection) -> None:
        self._connection = connection
        self.set_status(f"Scanning {root} ({connection.describe()})...")
        self.run_worker(lambda: self._scan(root, connection), thread=True, exclusive=True, group="scan")

    def _scan(self, root: str, connection: Connection) -> None:
        try:
            if not connection.is_remote:
                root = str(Path(root).expanduser().resolve())
            with connect(connection.kind, connection.host, connection.profile) as backend:
                result = self._scan_fn(backend, root, depth=self.depth)
        except Exception as exc:  # noqa: BLE001 - shown to the user
            self.call_from_thread(self._error, f"Scan failed: {exc}")
            return
        self.call_from_thread(self._apply_scan, result, connection)

    def _apply_scan(self, result, connection: Connection) -> None:
        if self.project is None or self.project.root != result.root:
            self.project = Project.from_scan(result, connection=connection)
        else:
            self.project.connection = connection
            self.project.merge_scan(result)
        refresh_stages(self.project)
        self.refresh_all()
        msg = f"Found {len(result.systems)} system(s) in {result.root}"
        if result.skipped_dirs:
            msg += f"; skipped {len(result.skipped_dirs)} folder(s) without trajectories"
        self.notify(msg, severity="information" if result.systems else "warning")

    # ------------------------------------------------------------ rendering

    def refresh_all(self) -> None:
        self._fill_systems()
        self._fill_ligands()
        self._fill_stages()
        self._fill_preview()
        if self.project is not None:
            n = len(self.project.active_systems)
            self.sub_title = f"{self.project.connection.describe()}:{self.project.root}"
            self.set_status(f"{n} system(s) selected · project file {self.project_path}")

    def _keep_cursor(self, table: DataTable, fill: Callable[[], None]) -> None:
        row = table.cursor_row
        table.clear()
        fill()
        if table.row_count:
            table.move_cursor(row=min(max(row, 0), table.row_count - 1))

    def _fill_systems(self) -> None:
        table = self.query_one("#systems", DataTable)

        def fill() -> None:
            for s in (self.project.systems if self.project else []):
                role = "APO ref" if s.apo else "complex"
                state = "missing" if s.missing else ("" if s.included else "excluded")
                name = s.name + (f" ({state})" if state else "")
                ligs = ", ".join(l.name for l in s.ligands) or "-"
                inputs = ", ".join(k.replace("_", " ") for k, v in s.inputs.items() if v) or "-"
                table.add_row(name, s.key, str(len(s.replicates)), role, ligs, inputs,
                              s.topology or "-", key=s.key)

        self._keep_cursor(table, fill)

    def _fill_ligands(self) -> None:
        table = self.query_one("#ligands", DataTable)

        def fill() -> None:
            for s in (self.project.systems if self.project else []):
                for lig in s.ligands:
                    table.add_row(lig.name, lig.resname, s.name, "yes" if lig.included else "no",
                                  lig.selection, key=f"{s.key}\t{lig.resname}")

        self._keep_cursor(table, fill)

    def _fill_stages(self) -> None:
        box = self.query_one("#stages", VerticalScroll)
        if not box.children:
            widgets = []
            for name, label in STAGES:
                widgets.append(Checkbox(f"{name}: {label}", id=f"stage-{name}"))
                widgets.append(Static("", id=f"reason-{name}", classes="stage-reason"))
            box.mount_all(widgets)
        for name, _ in STAGES:
            choice = self.project.stages.get(name) if self.project else None
            checkbox = box.query_one(f"#stage-{name}", Checkbox)
            reason = box.query_one(f"#reason-{name}", Static)
            available = bool(choice and choice.available)
            with checkbox.prevent(Checkbox.Changed):
                checkbox.value = bool(choice and choice.enabled and available)
            checkbox.disabled = not available
            reason.update(_escape(f"unavailable: {choice.reason}") if choice and not available else "")
            reason.display = bool(choice and not available)

    def _fill_preview(self) -> None:
        preview = self.query_one("#preview", Static)
        if self.project is None:
            preview.update("Nothing scanned yet (press o to connect).")
            return
        try:
            text = render_yaml(self.project, project_file=str(self.project_path))
        except ValueError as exc:
            preview.update(f"[b red]Cannot generate yet[/]\n{exc}")
            return
        preview.update(_escape(text))

    def set_status(self, text: str) -> None:
        self.query_one("#status", Static).update(_escape(text))

    def _error(self, message: str) -> None:
        self.set_status(message.splitlines()[0])
        self.notify(message, severity="error", timeout=8)

    # ------------------------------------------------------------ selection helpers

    def _active_table(self) -> Optional[str]:
        active = self.query_one("#tabs", TabbedContent).active
        return {"tab-systems": "systems", "tab-ligands": "ligands"}.get(active)

    def _selected_key(self, table_id: str) -> Optional[str]:
        table = self.query_one(f"#{table_id}", DataTable)
        if not table.row_count:
            return None
        return table.coordinate_to_cell_key((table.cursor_row, 0)).row_key.value

    def _selected_system(self):
        key = self._selected_key("systems")
        return self.project.system(key) if key and self.project else None

    def _selected_ligand(self):
        key = self._selected_key("ligands")
        if not key or self.project is None:
            return None, None
        system_key, resname = key.split("\t")
        system = self.project.system(system_key)
        return system, next(l for l in system.ligands if l.resname == resname)

    def _edit(self, fn: Callable[[], None]) -> bool:
        try:
            fn()
        except ProjectError as exc:
            self._error(str(exc))
            return False
        refresh_stages(self.project)
        self.refresh_all()
        return True

    # ------------------------------------------------------------ actions

    def action_edit(self) -> None:
        table = self._active_table()
        if table == "systems":
            system = self._selected_system()
            if system is None:
                return

            def done(value: Optional[str]) -> None:
                if value is not None and self._edit(lambda: self.project.rename_system(system.key, value)):
                    self.notify(f"{system.key} is now {system.name}")

            self.push_screen(EditNameScreen(f"Name for {system.key}", system.name), done)
        elif table == "ligands":
            system, lig = self._selected_ligand()
            if lig is None:
                return

            def done_lig(value: Optional[str]) -> None:
                if value is not None and self._edit(
                        lambda: self.project.rename_ligand(lig.resname, value, system=system.key)):
                    self.notify(f"{lig.resname} in {system.name} is labelled {lig.name}")

            self.push_screen(EditNameScreen(f"Label for {lig.resname} in {system.name}", lig.name), done_lig)

    def action_toggle_apo(self) -> None:
        if self._active_table() != "systems":
            return
        system = self._selected_system()
        if system is not None:
            self._edit(lambda: self.project.set_apo(system.key, not system.apo))

    def action_toggle_include(self) -> None:
        table = self._active_table()
        if table == "systems":
            system = self._selected_system()
            if system is not None:
                self._edit(lambda: self.project.set_included(system.key, not system.included))
        elif table == "ligands":
            system, lig = self._selected_ligand()
            if lig is not None:
                self._edit(lambda: self.project.set_ligand_included(lig.resname, not lig.included,
                                                                    system=system.key))

    @on(Checkbox.Changed)
    def _stage_toggled(self, event: Checkbox.Changed) -> None:
        name = (event.checkbox.id or "").removeprefix("stage-")
        if self.project is None or name not in self.project.stages:
            return
        self.project.stages[name].enabled = bool(event.value)
        self._fill_preview()

    def action_save(self) -> bool:
        if self.project is None:
            self._error("Nothing to save yet; press o to connect and scan")
            return False
        try:
            text = render_yaml(self.project, project_file=str(self.project_path))
            validate_yaml(text)
        except ValueError as exc:
            self._error(str(exc))
            return False
        self.project.save(self.project_path)
        self.pipeline_path.parent.mkdir(parents=True, exist_ok=True)
        self.pipeline_path.write_text(text)
        self.notify(f"Saved {self.project_path} and {self.pipeline_path}")
        self.set_status(f"Saved {self.project_path} · {self.pipeline_path}")
        return True

    def action_run(self) -> None:
        self._launch("host" if self.project and self.project.connection.is_remote else "local")

    def action_submit(self) -> None:
        if self.project is None or not self.project.connection.is_remote:
            self._error("HTCondor submission needs an SSH or condorctl connection (press o)")
            return
        self._launch("condor")

    def _launch(self, where: str) -> None:
        if not self.action_save():
            return
        project, text = self.project, self.pipeline_path.read_text()
        self.set_status(f"Starting ({where})...")
        self.run_worker(lambda: self._run(where, project, text), thread=True, exclusive=True, group="run")

    def _run(self, where: str, project: Project, text: str) -> None:
        from . import runner  # noqa: WPS433

        try:
            if where == "local":
                result = runner.run_local(self.pipeline_path)
                failed = result["failed_stages"]
                message = f"Run finished. Report: {result['report']}"
                if failed:
                    message += f" · failed: {', '.join(failed)}"
            else:
                c = project.connection
                with connect(c.kind, c.host, c.profile) as backend:
                    if where == "host":
                        result = runner.run_on_host(backend, project, text)
                        message = f"Started on {backend.label} (pid {result['pid']}); log {result['log']}"
                    else:
                        result = runner.submit_condor(backend, project, text)
                        message = f"Submitted cluster {result['cluster_id']}; output {result['log']}"
        except Exception as exc:  # noqa: BLE001 - shown to the user
            self.call_from_thread(self._error, f"{where} run failed: {exc}")
            return
        self.last_result = result
        self.call_from_thread(self.notify, message, timeout=10)
        self.call_from_thread(self.set_status, message)


def _escape(text: str) -> str:
    return text.replace("[", r"\[")


def run_app(project: Optional[Project] = None, root: Optional[str] = None,
            connection: Optional[Connection] = None, depth: int = 4,
            project_path: Path = Path(PROJECT_FILE)) -> int:
    ClaritySetupApp(project=project, root=root, connection=connection, depth=depth,
                    project_path=project_path).run()
    return 0


__all__ = ["ClaritySetupApp", "ConnectScreen", "EditNameScreen", "run_app", "STAGE_LABELS"]
