"""Simple file-based caching for pipeline stages."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Iterable, List, Optional


def _file_signature(path: Path) -> str:
    try:
        stat = path.stat()
        return f"{path}:{stat.st_size}:{int(stat.st_mtime)}"
    except FileNotFoundError:
        return f"{path}:missing"


def compute_signature(files: Iterable[Path], extra: Optional[dict] = None) -> str:
    items = [_file_signature(Path(p)) for p in files]
    if extra:
        items.append(json.dumps(extra, sort_keys=True))
    payload = "|".join(items)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def cache_path(output_dir: Path) -> Path:
    return output_dir / ".cache.json"


def load_cache(output_dir: Path) -> Optional[dict]:
    path = cache_path(output_dir)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def save_cache(output_dir: Path, signature: str, files: List[Path]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "signature": signature,
        "files": [str(p) for p in files],
    }
    cache_path(output_dir).write_text(json.dumps(payload, indent=2))


def should_skip(output_dir: Path, signature: str) -> bool:
    cache = load_cache(output_dir)
    if not cache:
        return False
    return cache.get("signature") == signature


def _normalise(value):
    """JSON-safe, order-stable form of a stage configuration."""
    import dataclasses
    from enum import Enum

    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {f.name: _normalise(getattr(value, f.name)) for f in dataclasses.fields(value)}
    if isinstance(value, dict):
        return {str(k): _normalise(v) for k, v in sorted(value.items(), key=lambda kv: str(kv[0]))}
    if isinstance(value, (list, tuple)):
        return [_normalise(v) for v in value]
    if isinstance(value, (set, frozenset)):
        return sorted(_normalise(v) for v in value)
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def config_fingerprint(config) -> str:
    """Hash of a stage (or unit) configuration: any setting change invalidates the cache."""
    payload = json.dumps(_normalise(config), sort_keys=True, default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()
