"""Where the simulation data lives: local disk, an SSH host, or a condorctl pool.

All backends expose the same small interface so discovery and runs are written
once. Discovery only needs to run shell commands and read/write small text files,
which keeps remote scans to a handful of round trips.

SSH is key-based only (no passwords are prompted for or stored) and records host
keys on first use, the same policy as condorctl.
"""

from __future__ import annotations

import os
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


class BackendError(RuntimeError):
    """Raised when a connection cannot be made or a remote operation fails."""


@dataclass
class CommandResult:
    status: int
    stdout: str
    stderr: str

    @property
    def ok(self) -> bool:
        return self.status == 0


class Backend:
    """Run commands and move small text files where the project lives."""

    kind = "base"
    is_remote = False

    def __init__(self) -> None:
        self.commands_run = 0

    @property
    def label(self) -> str:
        return self.kind

    def run(self, command: str, timeout: float = 120.0) -> CommandResult:
        raise NotImplementedError

    def read_text(self, path: str) -> str:
        result = self.run(f"cat {shlex.quote(path)}")
        if not result.ok:
            raise BackendError(f"Cannot read {path} on {self.label}: {result.stderr.strip()}")
        return result.stdout

    def write_text(self, path: str, text: str) -> None:
        raise NotImplementedError

    def exists(self, path: str) -> bool:
        return self.run(f"test -e {shlex.quote(path)}").ok

    def download(self, remote_path: str, local_path: str) -> None:
        """Copy one (possibly binary) file to this machine; base64 over ``run`` by default."""
        import base64  # noqa: WPS433

        result = self.run(f"base64 < {shlex.quote(remote_path)}", timeout=3600)
        if not result.ok:
            raise BackendError(f"Cannot download {remote_path} from {self.label}: {result.stderr.strip()}")
        target = Path(local_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(base64.b64decode("".join(result.stdout.split())))

    def close(self) -> None:
        """Release the connection (no-op for local)."""

    def __enter__(self) -> "Backend":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


class LocalBackend(Backend):
    """Data on this machine."""

    kind = "local"

    def run(self, command: str, timeout: float = 120.0) -> CommandResult:
        self.commands_run += 1
        proc = subprocess.run(command, shell=True, capture_output=True, text=True, timeout=timeout)
        return CommandResult(proc.returncode, proc.stdout, proc.stderr)

    def read_text(self, path: str) -> str:
        return Path(path).read_text(errors="replace")

    def write_text(self, path: str, text: str) -> None:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text)

    def exists(self, path: str) -> bool:
        return Path(path).exists()

    def download(self, remote_path: str, local_path: str) -> None:
        import shutil  # noqa: WPS433

        Path(local_path).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(remote_path, local_path)


class ShellRunner(LocalBackend):
    """Local commands, but treated as remote: every file access goes through ``run``.

    Used in tests so the exact command strings sent to SSH hosts are exercised.
    """

    kind = "shell"
    is_remote = True

    def read_text(self, path: str) -> str:
        return Backend.read_text(self, path)

    def write_text(self, path: str, text: str) -> None:
        directory = os.path.dirname(path) or "."
        quoted = shlex.quote(text)
        result = self.run(f"mkdir -p {shlex.quote(directory)} && printf %s {quoted} > {shlex.quote(path)}")
        if not result.ok:
            raise BackendError(f"Cannot write {path}: {result.stderr.strip()}")

    def exists(self, path: str) -> bool:
        return Backend.exists(self, path)

    def download(self, remote_path: str, local_path: str) -> None:
        Backend.download(self, remote_path, local_path)


class SSHBackend(Backend):
    """Data on a host reachable with ``ssh <alias>`` (reads ``~/.ssh/config``)."""

    kind = "ssh"
    is_remote = True

    def __init__(self, host: str, user: Optional[str] = None, port: Optional[int] = None,
                 key_file: Optional[str] = None) -> None:
        super().__init__()
        try:
            import paramiko  # noqa: WPS433
        except ImportError as exc:
            raise BackendError('SSH needs paramiko: pip install "gromacs-analysis-toolkit[tui]"') from exc
        self._paramiko = paramiko
        self.alias = host
        settings = self._ssh_config(host)
        self.host = settings.get("hostname", host)
        self.user = user or settings.get("user") or os.environ.get("USER")
        self.port = int(port or settings.get("port", 22))
        identity = key_file or (settings.get("identityfile") or [None])[0]
        self.key_file = os.path.expanduser(identity) if identity else None
        self._client = None
        self._sftp = None

    @property
    def label(self) -> str:
        return f"{self.user}@{self.host}" if self.user else self.host

    @staticmethod
    def _ssh_config(host: str) -> dict:
        path = Path.home() / ".ssh" / "config"
        if not path.exists():
            return {}
        import paramiko  # noqa: WPS433

        config = paramiko.SSHConfig()
        with open(path) as handle:
            config.parse(handle)
        return config.lookup(host)

    def _connect(self):
        if self._client is not None:
            return self._client
        paramiko = self._paramiko
        known_hosts = Path.home() / ".ssh" / "known_hosts"
        known_hosts.parent.mkdir(parents=True, exist_ok=True)
        known_hosts.touch(exist_ok=True)
        client = paramiko.SSHClient()
        client.load_host_keys(str(known_hosts))
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())  # trust on first use
        try:
            client.connect(
                hostname=self.host, port=self.port, username=self.user,
                key_filename=self.key_file, look_for_keys=self.key_file is None,
                allow_agent=True, timeout=20,
            )
        except paramiko.AuthenticationException as exc:
            key = self.key_file or "your default SSH keys/agent"
            raise BackendError(
                f"SSH key authentication failed for {self.label}. Check that the public half of "
                f"{key} is in ~/.ssh/authorized_keys on the host."
            ) from exc
        except Exception as exc:  # noqa: BLE001 - socket/SSH errors
            raise BackendError(f"Could not connect to {self.host}:{self.port}: {exc}") from exc
        client.save_host_keys(str(known_hosts))
        self._client = client
        return client

    def run(self, command: str, timeout: float = 120.0) -> CommandResult:
        client = self._connect()
        self.commands_run += 1
        try:
            _, stdout, stderr = client.exec_command(command, timeout=timeout)
            out = stdout.read().decode(errors="replace")
            err = stderr.read().decode(errors="replace")
            status = stdout.channel.recv_exit_status()
        except Exception as exc:  # noqa: BLE001
            raise BackendError(f"Remote command failed on {self.label}: {exc}") from exc
        return CommandResult(status, out, err)

    def _sftp_client(self):
        if self._sftp is None:
            self._sftp = self._connect().open_sftp()
        return self._sftp

    def write_text(self, path: str, text: str) -> None:
        directory = os.path.dirname(path)
        if directory:
            self.run(f"mkdir -p {shlex.quote(directory)}")
        with self._sftp_client().open(path, "w") as handle:
            handle.write(text)

    def download(self, remote_path: str, local_path: str) -> None:
        Path(local_path).parent.mkdir(parents=True, exist_ok=True)
        self._sftp_client().get(remote_path, local_path)

    def close(self) -> None:
        if self._sftp is not None:
            self._sftp.close()
            self._sftp = None
        if self._client is not None:
            self._client.close()
            self._client = None


class CondorctlBackend(Backend):
    """Data on an HTCondor pool configured for condorctl (e.g. profiles/nmrbox.yaml)."""

    kind = "condorctl"
    is_remote = True

    def __init__(self, profile: Optional[str] = None) -> None:
        super().__init__()
        try:
            from condorctl.config import Config  # noqa: WPS433
            from condorctl.connection import RemoteConnection  # noqa: WPS433
        except ImportError as exc:
            raise BackendError(
                "condorctl is not installed: pip install git+https://github.com/OASolliman590/condorctl"
            ) from exc
        try:
            self.config = Config.load(profile)
        except Exception as exc:  # noqa: BLE001 - ConfigError carries a helpful message
            raise BackendError(str(exc)) from exc
        self.profile = profile
        self._conn = RemoteConnection(self.config)
        self._connected = False

    @property
    def label(self) -> str:
        return f"{self.config.username}@{self.config.host} (condorctl)"

    def _ensure(self):
        if not self._connected:
            try:
                self._conn.connect()
            except Exception as exc:  # noqa: BLE001
                raise BackendError(str(exc)) from exc
            self._connected = True
        return self._conn

    def run(self, command: str, timeout: float = 120.0) -> CommandResult:
        conn = self._ensure()
        self.commands_run += 1
        result = conn.run(command, timeout=timeout)
        return CommandResult(result.exit_status, result.stdout, result.stderr)

    def write_text(self, path: str, text: str) -> None:
        conn = self._ensure()
        directory = os.path.dirname(path)
        if directory:
            conn.run(f"mkdir -p {shlex.quote(directory)}")
        with conn.sftp.open(path, "w") as handle:
            handle.write(text)

    def download(self, remote_path: str, local_path: str) -> None:
        Path(local_path).parent.mkdir(parents=True, exist_ok=True)
        self._ensure().sftp.get(remote_path, local_path)

    def close(self) -> None:
        if self._connected:
            self._conn.close()
            self._connected = False


def connect(kind: str = "local", host: Optional[str] = None, profile: Optional[str] = None) -> Backend:
    """Backend for a project connection (``local``, ``ssh`` with ``host``, ``condorctl`` with ``profile``)."""
    if kind == "local":
        return LocalBackend()
    if kind == "ssh":
        if not host:
            raise BackendError("SSH connection needs a host (alias from ~/.ssh/config or hostname)")
        return SSHBackend(host)
    if kind == "condorctl":
        return CondorctlBackend(profile)
    if kind == "shell":
        return ShellRunner()
    raise BackendError(f"Unknown connection kind: {kind}")
