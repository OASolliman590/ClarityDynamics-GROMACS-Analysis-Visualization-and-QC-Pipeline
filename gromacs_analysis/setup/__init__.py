"""Pipeline setup without hand-written YAML.

Discover simulation systems (locally or over SSH / condorctl), name systems and
ligands, choose analyses, and generate a validated pipeline configuration, from
the ``clarity-setup`` command line or its terminal interface.
"""

from .backends import Backend, LocalBackend, ShellRunner, connect
from .scanner import ScanResult, scan

__all__ = ["Backend", "LocalBackend", "ShellRunner", "connect", "ScanResult", "scan"]
