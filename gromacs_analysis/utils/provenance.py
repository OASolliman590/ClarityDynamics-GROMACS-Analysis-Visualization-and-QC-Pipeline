"""What ran, where and with which software: recorded in every run report."""

from __future__ import annotations

import getpass
import platform
import socket
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

PACKAGES = ("gromacs-analysis-toolkit", "numpy", "scipy", "pandas", "plotly", "MDAnalysis", "prolif", "rdkit",
            "kaleido", "PyYAML")


def _version(name: str) -> Optional[str]:
    try:
        from importlib.metadata import PackageNotFoundError, version  # noqa: WPS433
    except ImportError:  # pragma: no cover
        return None
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def git_state(path: Path) -> Dict[str, Optional[str]]:
    """Commit and dirty flag of the git checkout containing ``path`` (empty when not a checkout)."""
    def _git(*args: str) -> Optional[str]:
        try:
            out = subprocess.run(["git", "-C", str(path), *args], capture_output=True, text=True, timeout=5)
        except (OSError, subprocess.SubprocessError):
            return None
        return out.stdout.strip() if out.returncode == 0 else None

    commit = _git("rev-parse", "HEAD")
    if not commit:
        return {}
    dirty = _git("status", "--porcelain", "--untracked-files=no")
    return {"commit": commit, "dirty": "yes" if dirty else "no"}


def collect(config_path: Optional[Path] = None, started: Optional[datetime] = None) -> Dict:
    package_dir = Path(__file__).resolve().parents[2]
    info: Dict = {
        "started": started.isoformat(timespec="seconds") if started else None,
        "finished": datetime.now().isoformat(timespec="seconds"),
        "host": socket.gethostname(),
        "user": _safe_user(),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "command": " ".join(sys.argv) if sys.argv else None,
        "packages": {name: v for name in PACKAGES if (v := _version(name))},
        "toolkit_git": git_state(package_dir),
        "config_path": str(config_path) if config_path else None,
        "config_text": None,
    }
    if config_path and Path(config_path).exists():
        try:
            info["config_text"] = Path(config_path).read_text()
        except OSError:
            pass
    return info


def _safe_user() -> Optional[str]:
    try:
        return getpass.getuser()
    except Exception:  # noqa: BLE001 - containers without a passwd entry
        return None
