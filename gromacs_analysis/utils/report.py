"""Single self-contained HTML report for a pipeline run.

Lists every executed stage with its status, figures (images inline, interactive
HTML linked) and the first rows of its key tables, with QC warnings at the top.
All links are relative to the report so the output folder can be moved or shared.
"""

from __future__ import annotations

import html
import json
import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import pandas as pd

logger = logging.getLogger(__name__)

KEY_TABLES = (
    "statistics_summary.csv",
    "ligand_ranking.csv",
    "rank_pairwise.csv",
    "excluded_systems.csv",
    "entropy_summary.csv",
    "occupancy_summary.csv",
    "clusters_summary.csv",
    "eigenvalues.csv",
    "rmsip_pairwise.csv",
    "distance_compare_summary.csv",
    "stats_summary.csv",
)
MAX_TABLES_PER_STAGE = 4
MAX_ROWS = 10
REPORT_NAME = "report.html"

_STYLE = """
:root { --bg:#fcfcfb; --fg:#0b0b0b; --muted:#52514e; --card:#ffffff; --border:#e2e1dc;
        --ok:#0a7d32; --warn:#9a6700; --fail:#c62828; --link:#1c5cab; }
@media (prefers-color-scheme: dark) {
  :root { --bg:#1a1a19; --fg:#f0efec; --muted:#c3c2b7; --card:#232322; --border:#3a3a37;
          --ok:#5cc47a; --warn:#e0b44c; --fail:#ef6b6b; --link:#86b6ef; }
}
body { background:var(--bg); color:var(--fg); font:15px/1.5 system-ui,-apple-system,Segoe UI,sans-serif;
       margin:0; padding:24px 16px; }
main { max-width:1100px; margin:0 auto; }
h1 { font-size:24px; margin:0 0 4px; } h2 { font-size:19px; margin:0; } h3 { font-size:15px; margin:16px 0 6px; }
.meta { color:var(--muted); margin-bottom:20px; }
.card { background:var(--card); border:1px solid var(--border); border-radius:10px; padding:16px; margin:16px 0; }
.status { font-weight:600; margin-left:8px; } .ok { color:var(--ok); } .failed { color:var(--fail); }
.skipped { color:var(--muted); }
.warnings { border-left:4px solid var(--warn); } .warnings li { margin:4px 0; }
.figs { display:grid; grid-template-columns:repeat(auto-fill,minmax(300px,1fr)); gap:12px; }
.figs figure { margin:0; } .figs img { width:100%; border:1px solid var(--border); border-radius:6px; background:#fff; }
figcaption { font-size:13px; color:var(--muted); word-break:break-all; }
a { color:var(--link); }
.table-wrap { overflow-x:auto; } table { border-collapse:collapse; font-size:13px; }
th, td { border-bottom:1px solid var(--border); padding:4px 8px; text-align:left; white-space:nowrap; }
.links { columns:2 280px; font-size:13px; } code { font-size:13px; }
"""


def _rel(path: Path, root: Path) -> str:
    return Path(os.path.relpath(path, root)).as_posix()


def _stage_output_dirs(config: Any) -> List[Path]:
    configs = config if isinstance(config, list) else [config]
    dirs = []
    for cfg in configs:
        out = getattr(cfg, "output_dir", None)
        if out:
            dirs.append(Path(out))
    return dirs


def _status(result: Any) -> (str, Optional[str]):
    if result is None:
        return "not run", None
    items = result if isinstance(result, list) else [result]
    errors = [r.get("error") for r in items if isinstance(r, dict) and r.get("success") is False]
    if errors:
        return "failed", "; ".join(str(e) for e in errors if e)
    if all(isinstance(r, dict) and r.get("skipped") for r in items):
        return "skipped", "cache hit"
    return "ok", None


def _collect_figures(dirs: Iterable[Path]) -> Dict[str, Dict[str, Path]]:
    """stem -> {ext: path} for figures below the stage output directories."""
    figures: Dict[str, Dict[str, Path]] = {}
    for base in dirs:
        if not base.exists():
            continue
        for path in sorted(base.rglob("*")):
            if path.suffix.lower() not in {".png", ".svg", ".html"} or path.name == REPORT_NAME:
                continue
            key = str(path.with_suffix(""))
            figures.setdefault(key, {})[path.suffix.lower()] = path
    return figures


def _collect_tables(dirs: Iterable[Path]) -> List[Path]:
    found: List[Path] = []
    for base in dirs:
        if not base.exists():
            continue
        for name in KEY_TABLES:
            found.extend(sorted(base.rglob(name)))
    unique = list(dict.fromkeys(found))
    return unique[:MAX_TABLES_PER_STAGE]


def _qc_warnings(results: Dict[str, Any], dirs_by_stage: Dict[str, List[Path]]) -> List[str]:
    warnings: List[str] = []
    qc = results.get("qc")
    if isinstance(qc, dict):
        warnings.extend(str(w) for w in qc.get("warnings", []) or [])
        warnings.extend(f"ERROR: {e}" for e in qc.get("errors", []) or [])
    if not warnings:
        for path in (p / "qc_report.json" for p in dirs_by_stage.get("qc", [])):
            if path.exists():
                try:
                    data = json.loads(path.read_text())
                    warnings.extend(data.get("warnings", []))
                    warnings.extend(f"ERROR: {e}" for e in data.get("errors", []))
                except (OSError, ValueError):
                    pass
    return warnings


def _table_html(path: Path, root: Path) -> str:
    try:
        frame = pd.read_csv(path, nrows=MAX_ROWS)
    except Exception as exc:  # noqa: BLE001
        return f"<p>Could not read <code>{html.escape(_rel(path, root))}</code>: {html.escape(str(exc))}</p>"
    table = frame.to_html(index=False, float_format=lambda v: f"{v:.4g}", border=0, na_rep="")
    return (f'<h3><a href="{html.escape(_rel(path, root))}">{html.escape(path.name)}</a> '
            f'<small>(first {min(len(frame), MAX_ROWS)} rows)</small></h3><div class="table-wrap">{table}</div>')


def write_run_report(config: Any, results: Dict[str, Any], output_root: Optional[Path] = None) -> Path:
    """Write ``report.html`` into ``output_root`` (or the common parent of stage outputs)."""
    stages = [s for s in config.enabled_stages()] if hasattr(config, "enabled_stages") else []
    dirs_by_stage = {s.name: _stage_output_dirs(s.config) for s in stages}
    root = Path(output_root or getattr(config, "output_root", None) or "")
    if not str(root):
        all_dirs = [d.resolve() for dirs in dirs_by_stage.values() for d in dirs]
        root = Path(os.path.commonpath(all_dirs)) if all_dirs else Path.cwd()
    root.mkdir(parents=True, exist_ok=True)

    parts = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        f"<title>Pipeline Report</title><style>{_STYLE}</style></head><body><main>",
        f"<h1>{html.escape(str(getattr(config, 'project_name', 'Pipeline run')))}</h1>",
        f"<div class='meta'>Protein: {html.escape(str(getattr(config, 'protein_name', '') or '—'))} · "
        f"generated {datetime.now().strftime('%Y-%m-%d %H:%M')} · {len(stages)} stage(s)</div>",
    ]

    warnings = _qc_warnings(results, dirs_by_stage)
    if warnings:
        items = "".join(f"<li>{html.escape(w)}</li>" for w in warnings)
        parts.append(f"<section class='card warnings'><h2>QC warnings ({len(warnings)})</h2><ul>{items}</ul></section>")

    if not stages:
        parts.append("<p>No stages were enabled.</p>")
    for stage in stages:
        status, detail = _status(results.get(stage.name))
        dirs = dirs_by_stage.get(stage.name, [])
        parts.append(f"<section class='card' id='{html.escape(stage.name)}'>")
        parts.append(f"<h2>{html.escape(stage.name)}<span class='status {status.replace(' ', '-')}'>{status}</span></h2>")
        if detail:
            parts.append(f"<p class='{status}'>{html.escape(detail)}</p>")
        for d in dirs:
            parts.append(f"<div class='meta'>Output: <a href='{html.escape(_rel(d, root))}'>"
                         f"{html.escape(_rel(d, root))}</a></div>")

        figures = _collect_figures(dirs)
        images, links = [], []
        for stem, variants in figures.items():
            image = variants.get(".png") or variants.get(".svg")
            interactive = variants.get(".html")
            label = html.escape(Path(stem).name)
            if image is not None:
                target = interactive or image
                images.append(
                    f"<figure><a href='{html.escape(_rel(target, root))}'><img loading='lazy' "
                    f"src='{html.escape(_rel(image, root))}' alt='{label}'></a>"
                    f"<figcaption>{label}</figcaption></figure>"
                )
            elif interactive is not None:
                links.append(f"<li><a href='{html.escape(_rel(interactive, root))}'>{label}</a></li>")
        if images:
            parts.append(f"<h3>Figures</h3><div class='figs'>{''.join(images)}</div>")
        if links:
            parts.append(f"<h3>Interactive figures</h3><ul class='links'>{''.join(links)}</ul>")
        if not images and not links:
            parts.append("<p class='skipped'>No figures were produced for this stage.</p>")

        for table in _collect_tables(dirs):
            parts.append(_table_html(table, root))
        parts.append("</section>")

    parts.append("</main></body></html>")
    path = root / REPORT_NAME
    path.write_text("\n".join(parts), encoding="utf-8")
    logger.info(f"Run report: {path}")
    return path
