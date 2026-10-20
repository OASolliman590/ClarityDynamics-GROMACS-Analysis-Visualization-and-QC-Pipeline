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
    "pairwise_tests.csv",
    "excluded_systems.csv",
    "entropy_summary.csv",
    "occupancy_summary.csv",
    "clusters_summary.csv",
    "eigenvalues.csv",
    "rmsip_pairwise.csv",
    "distance_compare_summary.csv",
    "stats_summary.csv",
)
MAX_TABLES_PER_STAGE = 5
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
pre { white-space:pre-wrap; font-size:12px; } details { margin-top:12px; }
#provenance th { color:var(--muted); font-weight:500; }
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


MAX_FLAG_ROWS = 40


def _convergence_flags(dirs_by_stage: Dict[str, List[Path]]) -> List[Dict[str, str]]:
    """Flagged rows from every stage's data/convergence.csv (MD metrics, joint PCA)."""
    flagged: List[Dict[str, str]] = []
    for stage, dirs in dirs_by_stage.items():
        for base in dirs:
            path = base / "data" / "convergence.csv"
            if not path.exists():
                continue
            try:
                frame = pd.read_csv(path, dtype={"Flags": str}).fillna("")
            except Exception:  # noqa: BLE001
                continue
            for _, row in frame[frame["Flags"].astype(str).str.len() > 0].iterrows():
                what = " / ".join(str(row[c]) for c in ("System", "Metric") if c in frame.columns)
                detail = ", ".join(
                    f"{c}={row[c]:.3g}" for c in ("N_eff", "Half_Z", "Cos_PC1", "Half_RMSIP")
                    if c in frame.columns and isinstance(row[c], (int, float)) and row[c] == row[c]
                )
                flagged.append({"stage": stage, "what": what, "replicate": str(row.get("Replicate", "")),
                                "flags": str(row["Flags"]), "detail": detail, "file": str(path)})
    return flagged


_FLAG_HELP = {
    "low_neff": "fewer effective independent samples than the threshold: run longer or report wider uncertainty",
    "drift": "first and second half of the production window disagree: the property is still drifting",
    "too_short": "series too short to assess",
    "pc1_diffusive": "PC1 projection looks like random diffusion (cosine content high): conformational sampling not converged",
    "unstable_subspace": "dominant motions differ between first and second half of the frames",
}


def _convergence_html(flagged: List[Dict[str, str]], root: Path) -> str:
    if not flagged:
        return ""
    kinds = sorted({f for row in flagged for f in row["flags"].split(";") if f})
    legend = "".join(f"<li><code>{html.escape(k)}</code>: {html.escape(_FLAG_HELP.get(k, ''))}</li>" for k in kinds)
    rows = "".join(
        f"<tr><td>{html.escape(r['stage'])}</td><td>{html.escape(r['what'])}</td><td>{html.escape(r['replicate'])}</td>"
        f"<td>{html.escape(r['flags'])}</td><td>{html.escape(r['detail'])}</td></tr>"
        for r in flagged[:MAX_FLAG_ROWS]
    )
    more = f"<p>…and {len(flagged) - MAX_FLAG_ROWS} more.</p>" if len(flagged) > MAX_FLAG_ROWS else ""
    files = sorted({r["file"] for r in flagged})
    links = ", ".join(f"<a href='{html.escape(_rel(Path(f), root))}'>{html.escape(_rel(Path(f), root))}</a>"
                      for f in files)
    return (f"<section class='card warnings' id='sampling'><h2>Sampling checks ({len(flagged)} flagged)</h2>"
            f"<ul>{legend}</ul><div class='table-wrap'><table><tr><th>Stage</th><th>System / metric</th>"
            f"<th>Replicate</th><th>Flags</th><th>Values</th></tr>{rows}</table></div>{more}"
            f"<p class='meta'>Full tables: {links}</p></section>")


def _provenance_html(info: Dict[str, Any]) -> str:
    items = [
        ("Started", info.get("started")), ("Finished", info.get("finished")), ("Host", info.get("host")),
        ("User", info.get("user")), ("Platform", info.get("platform")), ("Python", info.get("python")),
        ("Command", info.get("command")), ("Configuration", info.get("config_path")),
    ]
    git = info.get("toolkit_git") or {}
    if git:
        items.append(("Toolkit commit", f"{git.get('commit')} (uncommitted changes: {git.get('dirty')})"))
    rows = "".join(f"<tr><th>{html.escape(k)}</th><td>{html.escape(str(v))}</td></tr>" for k, v in items if v)
    packages = "".join(f"<tr><th>{html.escape(k)}</th><td>{html.escape(v)}</td></tr>"
                       for k, v in (info.get("packages") or {}).items())
    config_text = info.get("config_text")
    config_block = (f"<details><summary>Configuration that ran</summary><pre>{html.escape(config_text)}</pre></details>"
                    if config_text else "")
    return (f"<section class='card' id='provenance'><h2>Provenance</h2><div class='table-wrap'><table>{rows}</table>"
            f"</div><h3>Software</h3><div class='table-wrap'><table>{packages}</table></div>{config_block}</section>")


def _table_html(path: Path, root: Path) -> str:
    try:
        frame = pd.read_csv(path, nrows=MAX_ROWS)
    except Exception as exc:  # noqa: BLE001
        return f"<p>Could not read <code>{html.escape(_rel(path, root))}</code>: {html.escape(str(exc))}</p>"
    table = frame.to_html(index=False, float_format=lambda v: f"{v:.4g}", border=0, na_rep="")
    return (f'<h3><a href="{html.escape(_rel(path, root))}">{html.escape(path.name)}</a> '
            f'<small>(first {min(len(frame), MAX_ROWS)} rows)</small></h3><div class="table-wrap">{table}</div>')


def write_run_report(config: Any, results: Dict[str, Any], output_root: Optional[Path] = None,
                     started: Optional[datetime] = None) -> Path:
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

    parts.append(_convergence_html(_convergence_flags(dirs_by_stage), root))

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

    try:
        from .provenance import collect  # noqa: WPS433

        parts.append(_provenance_html(collect(getattr(config, "source_path", None), started)))
    except Exception as exc:  # noqa: BLE001 - provenance must never break the report
        logger.debug(f"Provenance unavailable: {exc}")
    parts.append("</main></body></html>")
    path = root / REPORT_NAME
    path.write_text("\n".join(parts), encoding="utf-8")
    logger.info(f"Run report: {path}")
    return path
