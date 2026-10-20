"""Keyboard-driven terminal interface for clarity-setup (Textual).

Connect -> review systems and ligands (rename, APO, include) -> prepare missing
inputs -> pick stages and their settings -> preview the generated YAML -> save / run /
submit (single job or parallel DAG) -> follow runs and fetch results. Every action goes
through the same :class:`Project`, generator, runner and job helpers used by the CLI;
errors become notifications instead of crashing the interface.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, List, Optional

from textual import on
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.coordinate import Coordinate
from textual.screen import ModalScreen
from textual.widget import Widget
from textual.widgets import (
    Button,
    Checkbox,
    DataTable,
    Footer,
    Header,
    Input,
    Label,
    Select,
    Static,
    TabbedContent,
    TabPane,
)

from . import jobs
from . import settings as stage_settings
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
        kind = str(self.query_one("#kind", Select).value)
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


class SettingsScreen(ModalScreen[Optional[dict]]):
    """Edit one stage's settings; dismisses with {setting name: text} or None."""

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def __init__(self, project: Project, stage: str) -> None:
        super().__init__()
        self._project = project
        self._stage = stage

    def compose(self) -> ComposeResult:
        rows = stage_settings.settings_for(self._stage)
        with Vertical(id="dialog", classes="wide"):
            yield Label(f"Settings: {self._stage} ({STAGE_LABELS.get(self._stage, '')})", id="dialog-title")
            with VerticalScroll(id="settings-list"):
                for setting in rows:
                    value = stage_settings.current(self._project, setting)
                    kind = setting.kind if setting.kind != "choice" else " | ".join(setting.choices)
                    default = "none" if setting.default is None else setting.default
                    yield Label(f"{setting.key}  [{kind}; default {default}]", classes="setting-label")
                    yield Static(setting.help, classes="hint")
                    yield Input(value="none" if value is None else str(value), id=f"set-{_widget_id(setting.key)}",
                                name=setting.name)
            with Horizontal(id="dialog-buttons"):
                yield Button("Save", variant="primary", id="ok")
                yield Button("Defaults", id="defaults")
                yield Button("Cancel", id="cancel")

    def on_mount(self) -> None:
        inputs = list(self.query(Input))
        if inputs:
            inputs[0].focus()

    @on(Input.Submitted)
    @on(Button.Pressed, "#ok")
    def _accept(self) -> None:
        self.dismiss({field.name: field.value for field in self.query(Input)})

    @on(Button.Pressed, "#defaults")
    def _defaults(self) -> None:
        for setting in stage_settings.settings_for(self._stage):
            field = self.query_one(f"#set-{_widget_id(setting.key)}", Input)
            field.value = "none" if setting.default is None else str(setting.default)

    @on(Button.Pressed, "#cancel")
    def action_cancel(self) -> None:
        self.dismiss(None)


class ChoiceScreen(ModalScreen[Optional[str]]):
    """Text plus a row of choices (button ids); dismisses with the chosen id or None."""

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def __init__(self, title: str, body: str, choices: List[tuple]) -> None:
        super().__init__()
        self._title, self._body, self._choices = title, body, choices

    def compose(self) -> ComposeResult:
        with Vertical(id="dialog", classes="wide"):
            yield Label(self._title, id="dialog-title")
            with VerticalScroll(id="choice-body"):
                yield Static(_escape(self._body), id="choice-text")
            with Horizontal(id="dialog-buttons"):
                for index, (choice_id, label) in enumerate(self._choices):
                    yield Button(label, id=f"choice-{choice_id}", variant="primary" if index == 0 else "default")

    def on_mount(self) -> None:
        buttons = list(self.query(Button))
        if buttons:
            buttons[0].focus()

    @on(Button.Pressed)
    def _pressed(self, event: Button.Pressed) -> None:
        choice = (event.button.id or "").removeprefix("choice-")
        self.dismiss(None if choice == "cancel" else choice)

    def action_cancel(self) -> None:
        self.dismiss(None)


def _widget_id(key: str) -> str:
    return key.replace(".", "__")


# ---------------------------------------------------------------- main app

class ClaritySetupApp(App):
    TITLE = "ClarityDynamics setup"
    CSS = """
    #dialog { width: 70; height: auto; padding: 1 2; border: thick $accent; background: $surface; }
    EditNameScreen, ConnectScreen, SettingsScreen, ChoiceScreen { align: center middle; }
    #dialog.wide { width: 120; max-width: 95%; height: 85%; }
    #settings-list, #choice-body { height: 1fr; }
    #settings-list Input { margin-bottom: 0; }
    .setting-label { text-style: bold; margin-top: 1; }
    #jobs { height: 1fr; }
    #job-log { height: auto; max-height: 14; padding: 0 1; border-top: solid $panel; }
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
        Binding("d", "submit_dag", "Parallel"),
        Binding("p", "prepare", "Prepare"),
        Binding("g", "settings", "Settings"),
        Binding("j", "jobs", "Jobs"),
        Binding("f", "fetch", "Fetch"),
        Binding("D", "doctor", "Doctor"),
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
                yield Static("space toggle · g settings of the focused stage", classes="hint")
                yield VerticalScroll(id="stages")
            with TabPane("Preview", id="tab-preview"):
                with VerticalScroll():
                    yield Static("", id="preview")
            with TabPane("Jobs", id="tab-jobs"):
                yield Static("j refresh status · f fetch results of the selected run · "
                             "r run · c submit · d parallel DAG · p prepare inputs", classes="hint")
                yield DataTable(id="jobs", cursor_type="row", zebra_stripes=True)
                yield Static("", id="job-log")
        yield Static("", id="status")
        yield Footer()

    def on_mount(self) -> None:
        systems = self.query_one("#systems", DataTable)
        systems.add_columns("Name", "Folder", "Reps", "Role", "Ligands", "Inputs", "Topology")
        ligands = self.query_one("#ligands", DataTable)
        ligands.add_columns("Label", "Residue", "System", "Included", "Selection")
        self.query_one("#jobs", DataTable).add_columns("Run", "Kind", "Where", "State", "Started", "Detail")
        self.job_info: dict = {}
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
        self._fill_jobs()
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
                ligs = ", ".join(lig.name for lig in s.ligands) or "-"
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
            widgets: List[Widget] = []
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

    def _fill_jobs(self) -> None:
        table = self.query_one("#jobs", DataTable)

        def fill() -> None:
            for data in reversed(self.project.runs if self.project else []):
                record = jobs.RunRecord.from_dict(data)
                info = self.job_info.get(record.id, {})
                detail = info.get("detail") or record.note or ""
                table.add_row(record.id, record.kind, record.where, (record.last_state or "?").upper(),
                              record.started.replace("T", " "), detail, key=record.id)

        self._keep_cursor(table, fill)
        self._show_job_log()

    def _show_job_log(self) -> None:
        key = self._selected_key("jobs")
        tail = self.job_info.get(key, {}).get("log_tail", "") if key else ""
        self.query_one("#job-log", Static).update(_escape(tail.rstrip()) if tail.strip() else
                                                  "Select a run and press j to refresh its status and log.")

    @on(DataTable.RowHighlighted, "#jobs")
    def _job_highlighted(self) -> None:
        self._show_job_log()

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
        return table.coordinate_to_cell_key(Coordinate(table.cursor_row, 0)).row_key.value

    def _selected_system(self):
        key = self._selected_key("systems")
        return self.project.system(key) if key and self.project else None

    def _selected_ligand(self):
        key = self._selected_key("ligands")
        if not key or self.project is None:
            return None, None
        system_key, resname = key.split("\t")
        system = self.project.system(system_key)
        return system, next(lig for lig in system.ligands if lig.resname == resname)

    def _edit(self, fn: Callable[[], object]) -> bool:
        try:
            fn()
        except ProjectError as exc:
            self._error(str(exc))
            return False
        refresh_stages(self._proj())
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
                if value is not None and self._edit(lambda: self._proj().rename_system(system.key, value)):
                    self.notify(f"{system.key} is now {system.name}")

            self.push_screen(EditNameScreen(f"Name for {system.key}", system.name), done)
        elif table == "ligands":
            system, lig = self._selected_ligand()
            if lig is None:
                return

            def done_lig(value: Optional[str]) -> None:
                if value is not None and self._edit(
                        lambda: self._proj().rename_ligand(lig.resname, value, system=system.key)):
                    self.notify(f"{lig.resname} in {system.name} is labelled {lig.name}")

            self.push_screen(EditNameScreen(f"Label for {lig.resname} in {system.name}", lig.name), done_lig)

    def action_toggle_apo(self) -> None:
        if self._active_table() != "systems":
            return
        system = self._selected_system()
        if system is not None:
            self._edit(lambda: self._proj().set_apo(system.key, not system.apo))

    def action_toggle_include(self) -> None:
        table = self._active_table()
        if table == "systems":
            system = self._selected_system()
            if system is not None:
                self._edit(lambda: self._proj().set_included(system.key, not system.included))
        elif table == "ligands":
            system, lig = self._selected_ligand()
            if lig is not None:
                self._edit(lambda: self._proj().set_ligand_included(lig.resname, not lig.included,
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

    def action_submit_dag(self) -> None:
        if self.project is None or not self.project.connection.is_remote:
            self._error("The parallel HTCondor workflow needs an SSH or condorctl connection (press o)")
            return
        self._launch("dag")

    def _launch(self, where: str) -> None:
        if not self.action_save():
            return
        project, text = self._proj(), self.pipeline_path.read_text()
        self.set_status(f"Starting ({where})...")
        self.run_worker(lambda: self._run(where, project, text), thread=True, exclusive=True, group="run")

    def _run(self, where: str, project: Project, text: str) -> None:
        from . import runner, workflow  # noqa: WPS433

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
                    elif where == "dag":
                        result = workflow.submit_dag(backend, project, text)
                        message = (f"Submitted DAG cluster {result['cluster_id']}: {len(result['units'])} "
                                   "parallel unit(s) + final job")
                    else:
                        result = runner.submit_condor(backend, project, text)
                        message = f"Submitted cluster {result['cluster_id']}; output {result['log']}"
        except Exception as exc:  # noqa: BLE001 - shown to the user
            self.call_from_thread(self._error, f"{where} run failed: {exc}")
            return
        self.last_result = result
        self.call_from_thread(self._record_run, result.get("record"), message)

    def _record_run(self, record, message: str) -> None:
        if record is not None and self.project is not None:
            jobs.add_run(self.project, record)
            try:
                self.project.save(self.project_path)
            except OSError as exc:
                self._error(f"Could not save the run record: {exc}")
            self._fill_jobs()
        self.notify(message, timeout=10)
        self.set_status(message)

    # ------------------------------------------------------------ settings

    def _focused_stage(self) -> Optional[str]:
        focused = self.focused
        name = (getattr(focused, "id", None) or "").removeprefix("stage-")
        return name if name in STAGE_LABELS else None

    def action_settings(self, stage: Optional[str] = None) -> None:
        if self.project is None:
            self._error("Nothing to configure yet; press o to connect and scan")
            return
        stage = stage or self._focused_stage()
        if stage is None:
            with_settings = [name for name, _ in STAGES if stage_settings.settings_for(name)]
            def chosen_stage(chosen: Optional[str]) -> None:
                if chosen:
                    self.action_settings(chosen)

            self.push_screen(ChoiceScreen("Settings for which stage?", "Pick a stage to edit its settings.",
                                          [(name, name) for name in with_settings] + [("cancel", "Cancel")]),
                             chosen_stage)
            return
        if not stage_settings.settings_for(stage):
            self.notify(f"{stage} has no adjustable settings")
            return
        self.push_screen(SettingsScreen(self.project, stage), lambda values: self._apply_settings(stage, values))

    def _apply_settings(self, stage: str, values: Optional[dict]) -> None:
        if not values:
            return
        try:
            for name, raw in values.items():
                stage_settings.set_value(self._proj(), name, raw)
        except ProjectError as exc:
            self._error(str(exc))
            return
        self._fill_preview()
        changed = self._proj().settings.get(stage, {})
        self.notify(f"{stage}: " + (", ".join(f"{k}={v}" for k, v in changed.items()) or "defaults"))

    # ------------------------------------------------------------ jobs

    def action_jobs(self) -> None:
        self.query_one("#tabs", TabbedContent).active = "tab-jobs"
        if self.project is None or not self.project.runs:
            self.set_status("No runs recorded yet")
            return
        key = self._selected_key("jobs")
        record = jobs.get_run(self.project, key)
        self.set_status(f"Checking {record.id}...")
        self.run_worker(lambda: self._check(record), thread=True, exclusive=True, group="jobs")

    def _check(self, record) -> None:
        try:
            with self._project_backend() as backend:
                info = jobs.status(backend, record, lines=30)
        except Exception as exc:  # noqa: BLE001
            self.call_from_thread(self._error, f"Status failed: {exc}")
            return
        self.call_from_thread(self._apply_status, record, info)

    def _apply_status(self, record, info: dict) -> None:
        self.job_info[record.id] = info
        project = self._proj()
        jobs.update_run(project, record)
        try:
            project.save(self.project_path)
        except OSError:
            pass
        self._fill_jobs()
        self.set_status(f"{record.id}: {info['state']} {info.get('detail', '')}")

    def action_fetch(self) -> None:
        if self.project is None or not self.project.runs:
            self._error("No runs to fetch")
            return
        record = jobs.get_run(self.project, self._selected_key("jobs"))
        dest = self.project_path.parent / "clarity_results"
        self.set_status(f"Fetching {record.id}...")
        self.run_worker(lambda: self._fetch(record, dest), thread=True, exclusive=True, group="jobs")

    def _fetch(self, record, dest: Path) -> None:
        try:
            with self._project_backend() as backend:
                if record.last_state not in jobs.FINISHED:
                    jobs.status(backend, record, lines=0)
                report = jobs.fetch(backend, record, dest)
        except Exception as exc:  # noqa: BLE001
            self.call_from_thread(self._error, f"Fetch failed: {exc}")
            return
        self.call_from_thread(self._apply_status, record, {"state": record.last_state, "detail": f"fetched to {report}"})
        self.call_from_thread(self.notify, f"Results: {report}", timeout=12)

    def _proj(self) -> Project:
        """The loaded project (actions that need one check ``self.project`` first)."""
        if self.project is None:
            raise ProjectError("No project loaded")
        return self.project

    def _project_backend(self):
        c = self._proj().connection
        return connect(c.kind, c.host, c.profile)

    # ------------------------------------------------------------ prepare

    def action_prepare(self) -> None:
        from . import prepare  # noqa: WPS433

        if self.project is None:
            self._error("Nothing scanned yet; press o to connect")
            return
        tasks = prepare.plan_preparation(self.project)
        if not tasks:
            self.notify("Nothing to prepare: every included replicate has its inputs")
            return
        runnable = [t for t in tasks if t.ok]
        body = "\n".join(prepare.summarise(tasks))
        choices = []
        if runnable:
            choices.append(("host", "Run on host" if self.project.connection.is_remote else "Run here"))
            if self.project.connection.is_remote:
                choices.append(("condor", "Submit to HTCondor"))
        choices.append(("cancel", "Cancel"))
        title = f"Prepare inputs: {len(runnable)} task(s) ready, {len(tasks) - len(runnable)} skipped"
        def chosen(where: Optional[str]) -> None:
            if where:
                self._run_prepare(where, tasks)

        self.push_screen(ChoiceScreen(title, body, choices), chosen)

    def _run_prepare(self, where: str, tasks) -> None:
        from . import prepare  # noqa: WPS433

        project = self._proj()
        scripts = prepare.build_scripts(project, tasks)

        def work() -> None:
            try:
                with self._project_backend() as backend:
                    result = prepare.run_preparation(backend, project, scripts, where=where)
            except Exception as exc:  # noqa: BLE001
                self.call_from_thread(self._error, f"Preparation failed: {exc}")
                return
            message = (f"Submitted {len(scripts)} preparation job(s) (cluster {result['cluster_id']})"
                       if result["where"] == "condor" else f"Started {len(scripts)} preparation script(s)")
            self.call_from_thread(self._record_run, result.get("record"), message + "; rescan (ctrl+r) when done")

        self.run_worker(work, thread=True, exclusive=True, group="run")

    # ------------------------------------------------------------ doctor

    def action_doctor(self) -> None:
        from . import doctor  # noqa: WPS433

        if self.project is None:
            self._error("Connect first (press o)")
            return
        project = self.project
        modes = ["host" if project.connection.is_remote else "local"]
        if project.connection.is_remote:
            modes.append("dag")

        def work() -> None:
            try:
                with self._project_backend() as backend:
                    checks = doctor.run_checks(backend, project.root, modes=modes)
            except Exception as exc:  # noqa: BLE001
                self.call_from_thread(self._error, f"Doctor failed: {exc}")
                return
            body = "\n".join(doctor.format_table(checks))
            title = "Environment check: " + ("all required checks passed" if doctor.exit_code(checks) == 0
                                             else "some required checks failed")
            self.call_from_thread(self.push_screen, ChoiceScreen(title, body, [("cancel", "Close")]))

        self.set_status("Checking the environment...")
        self.run_worker(work, thread=True, exclusive=True, group="doctor")


def _escape(text: str) -> str:
    return text.replace("[", r"\[")


def run_app(project: Optional[Project] = None, root: Optional[str] = None,
            connection: Optional[Connection] = None, depth: int = 4,
            project_path: Path = Path(PROJECT_FILE)) -> int:
    ClaritySetupApp(project=project, root=root, connection=connection, depth=depth,
                    project_path=project_path).run()
    return 0


__all__ = ["ChoiceScreen", "ClaritySetupApp", "ConnectScreen", "EditNameScreen", "SettingsScreen", "run_app",
           "STAGE_LABELS"]
