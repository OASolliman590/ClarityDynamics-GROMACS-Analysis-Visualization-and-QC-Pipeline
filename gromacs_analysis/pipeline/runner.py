"""
Pipeline runner for multi-stage analyses.
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from .config import PipelineConfig
from ..md import MDAnalyzer
from ..mmpbsa import MMPBSAAnalyzer
from ..pca import PCAAnalyzer, PCAComparisonAnalyzer, PCAJointAnalyzer
from ..distance import DistanceAnalyzer, DistanceComparisonAnalyzer
from ..interactions import InteractionAnalyzer, InteractionComparisonAnalyzer
from ..ranking import RankingAnalyzer
from ..prolif import ProlifAnalyzer, ProlifComparisonAnalyzer
from ..ttclust import TTClustAnalyzer
from ..qc import QCAnalyzer
from .cache import compute_signature, config_fingerprint, should_skip, save_cache
from .utils import run_utils_tasks
from ..utils.report import write_run_report

logger = logging.getLogger(__name__)


BATCH_UNITS = {
    "md_batch": ("md", MDAnalyzer),
    "pca_batch": ("pca", PCAAnalyzer),
    "distance_batch": ("distance", DistanceAnalyzer),
    "interactions_batch": ("interactions", InteractionAnalyzer),
    "prolif_batch": ("prolif", ProlifAnalyzer),
    "ttclust_batch": ("ttclust", TTClustAnalyzer),
}


def _signature(stage_name: str, config: Any) -> (str, List[Path]):
    """Input files (path, size, mtime) plus a fingerprint of the configuration itself."""
    files = _collect_stage_files(stage_name, config)
    return compute_signature(files, extra={"stage": stage_name, "config": config_fingerprint(config)}), files


def _succeeded(result: Any) -> bool:
    items = result if isinstance(result, list) else [result]
    return not any(isinstance(r, dict) and r.get("success") is False for r in items)


def _run_batch(name: str, configs: List[Any], use_cache: bool) -> List[Dict[str, Any]]:
    """Run each unit of a batch stage; with caching, unchanged units are skipped individually."""
    unit_stage, analyzer_cls = BATCH_UNITS[name]
    out: List[Dict[str, Any]] = []
    for idx, cfg in enumerate(configs):
        output_dir = getattr(cfg, "output_dir", None)
        signature = files = None
        if use_cache and output_dir:
            signature, files = _signature(unit_stage, cfg)
            if should_skip(Path(output_dir), signature):
                logger.info(f"  - {name} unit {idx + 1}/{len(configs)}: cache hit, skipped")
                out.append({"success": True, "skipped": True, "output_dir": str(output_dir)})
                continue
        if output_dir:
            Path(output_dir).mkdir(parents=True, exist_ok=True)
        result = analyzer_cls(cfg).run_analysis()
        out.append(result)
        if signature is not None and _succeeded(result):
            save_cache(Path(output_dir), signature, files)
        logger.info(f"  ✓ Completed {name} unit {idx + 1}/{len(configs)}")
    return out


def run_pipeline(config: PipelineConfig, continue_on_error: bool = False,
                 use_cache: bool = True) -> Dict[str, Any]:
    """Run enabled pipeline stages in order and return results.

    Stages with ``cache: true`` are skipped when neither their inputs nor their settings
    changed; batch stages are cached per unit. ``use_cache=False`` forces a full run.
    """
    results: Dict[str, Any] = {}
    started = datetime.now()

    for stage in config.enabled_stages():
        name = stage.name
        logger.info(f"Starting stage: {name}")

        try:
            caching = use_cache and stage.cache_enabled
            signature = files = None
            if caching and getattr(stage.config, "output_dir", None):
                output_dir = stage.config.output_dir
                signature, files = _signature(name, stage.config)
                if should_skip(output_dir, signature):
                    logger.info(f"Skipping stage '{name}' (cache hit)")
                    results[name] = {"success": True, "skipped": True}
                    continue
            if name in BATCH_UNITS:
                results[name] = _run_batch(name, stage.config, caching)
            elif name == "md":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = MDAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "qc":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = QCAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "md_compare":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = MDAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "mmpbsa":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = MMPBSAAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "pca":
                if stage.config.output_dir:
                    stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = PCAAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "pca_compare":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = PCAComparisonAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "pca_joint":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                results[name] = PCAJointAnalyzer(stage.config).run_analysis()
            elif name == "distance":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = DistanceAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "distance_compare":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = DistanceComparisonAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "interactions":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = InteractionAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "interactions_compare":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = InteractionComparisonAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "prolif":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = ProlifAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "prolif_compare":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = ProlifComparisonAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "ttclust":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = TTClustAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "ranking":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = RankingAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "utils":
                results[name] = run_utils_tasks(stage.config)
            else:
                logger.warning(f"Unknown stage '{name}' - skipping")
                continue

            if signature is not None and _succeeded(results.get(name)):
                save_cache(stage.config.output_dir, signature, files)

            logger.info(f"Finished stage: {name}")

        except Exception as exc:
            logger.error(f"Stage '{name}' failed: {exc}")
            results[name] = {"success": False, "error": str(exc)}
            if not continue_on_error:
                _write_report(config, results, started)
                raise

    report = _write_report(config, results, started)
    if report is not None:
        results["_report"] = str(report)
    return results


def _write_report(config: PipelineConfig, results: Dict[str, Any], started: Optional[datetime] = None):
    """Write the run report; a report failure never fails the pipeline."""
    try:
        return write_run_report(config, results, started=started)
    except Exception as exc:  # noqa: BLE001
        logger.warning(f"Could not write run report: {exc}")
        return None


def _collect_stage_files(stage_name: str, config: Any) -> List[Path]:
    """Collect input files for caching purposes."""
    from pathlib import Path

    files: List[Path] = []

    if stage_name in ("md", "md_compare"):
        from ..md.data_processor import MDDataProcessor  # noqa: WPS433

        processor = MDDataProcessor(config)
        found = processor.find_data_files()
        for system_data in found.values():
            for metric_info in system_data.values():
                files.extend(metric_info.get("files", []))

    elif stage_name == "md_batch":
        for cfg in config:
            files.extend(_collect_stage_files("md", cfg))

    elif stage_name == "qc":
        if getattr(config, "topology", None):
            files.append(Path(config.topology))
        if getattr(config, "index_file", None):
            files.append(Path(config.index_file))
        files.extend([Path(p) for p in getattr(config, "trajectories", [])])

    elif stage_name == "mmpbsa":
        from ..mmpbsa.processor import MMPBSAProcessor  # noqa: WPS433

        processor = MMPBSAProcessor(config)
        for system in config.systems:
            results_files, decomp_files = processor._find_system_files(system)  # noqa: WPS437
            files.extend(results_files)
            files.extend(decomp_files)

    elif stage_name == "pca":
        if config.run_gromacs.enabled:
            input_dir = Path(config.run_gromacs.input_dir) if config.run_gromacs.input_dir else Path(config.base_dir)
            terminal = config.run_gromacs.terminal_detection
            for name in (
                terminal.rmsf_file,
                terminal.gro_file,
                terminal.pdb_file,
                config.run_gromacs.tpr,
                config.run_gromacs.trajectory,
            ):
                path = Path(name)
                files.append(path if path.is_absolute() else input_dir / path)
        from ..pca.processor import PCAProcessor  # noqa: WPS433

        processor = PCAProcessor(config)
        detected = processor.auto_detect_files()
        for val in detected.values():
            if isinstance(val, dict):
                files.extend([Path(p) for p in val.values() if p])
            elif val:
                files.append(Path(val))

    elif stage_name == "pca_batch":
        for cfg in config:
            files.extend(_collect_stage_files("pca", cfg))

    elif stage_name == "distance":
        files.append(Path(config.topology))
        files.extend([Path(p) for p in config.trajectories])

    elif stage_name == "distance_batch":
        for cfg in config:
            files.extend(_collect_stage_files("distance", cfg))

    elif stage_name == "interactions":
        files.append(Path(config.topology))
        files.extend([Path(p) for p in config.trajectories])

    elif stage_name == "interactions_batch":
        for cfg in config:
            files.extend(_collect_stage_files("interactions", cfg))

    elif stage_name == "pca_joint":
        for system in config.systems:
            base = Path(config.base_dir) if config.base_dir else None
            for path in [system.topology, *system.trajectories]:
                path = Path(path)
                files.append(path if path.is_absolute() or base is None else base / path)
    elif stage_name == "ranking":
        if config.md_stats_file:
            files.append(Path(config.md_stats_file))
        if config.mmpbsa_data_dir and Path(config.mmpbsa_data_dir).exists():
            files.extend(list(Path(config.mmpbsa_data_dir).glob("*_results.csv")))
        if config.distance_summary_file:
            files.append(Path(config.distance_summary_file))
    elif stage_name == "prolif":
        files.append(Path(config.topology))
        files.extend([Path(p) for p in config.trajectories])
    elif stage_name == "prolif_batch":
        for cfg in config:
            files.extend(_collect_stage_files("prolif", cfg))
    elif stage_name == "prolif_compare":
        for cfg in config.systems:
            occ_path = Path(cfg.output_dir) / "data" / "occupancy.csv"
            files.append(occ_path)
    elif stage_name == "ttclust":
        if getattr(config, "topology", None):
            files.append(Path(config.topology))
        files.extend([Path(p) for p in getattr(config, "trajectories", [])])
    elif stage_name == "ttclust_batch":
        for cfg in config:
            files.extend(_collect_stage_files("ttclust", cfg))

    return [Path(p) for p in files]
