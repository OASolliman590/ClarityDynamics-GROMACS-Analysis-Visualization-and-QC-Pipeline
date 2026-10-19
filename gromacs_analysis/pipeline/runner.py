"""
Pipeline runner for multi-stage analyses.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List

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
from .cache import compute_signature, should_skip, save_cache
from .utils import run_utils_tasks

logger = logging.getLogger(__name__)


def run_pipeline(config: PipelineConfig, continue_on_error: bool = False) -> Dict[str, Any]:
    """Run enabled pipeline stages in order and return results."""
    results: Dict[str, Any] = {}

    for stage in config.enabled_stages():
        name = stage.name
        logger.info(f"Starting stage: {name}")

        try:
            if stage.cache_enabled and getattr(stage.config, "output_dir", None):
                output_dir = stage.config.output_dir
                files = _collect_stage_files(name, stage.config)
                signature = compute_signature(files, extra={"stage": name})
                if should_skip(output_dir, signature):
                    logger.info(f"Skipping stage '{name}' (cache hit)")
                    results[name] = {"success": True, "skipped": True}
                    continue
            if name == "md":
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
            elif name == "md_batch":
                batch_results = []
                for idx, md_config in enumerate(stage.config):
                    md_config.output_dir.mkdir(parents=True, exist_ok=True)
                    analyzer = MDAnalyzer(md_config)
                    batch_results.append(analyzer.run_analysis())
                    logger.info(f"  ✓ Completed batch {idx + 1}/{len(stage.config)}")
                results[name] = batch_results
            elif name == "mmpbsa":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = MMPBSAAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "pca":
                if stage.config.output_dir:
                    stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = PCAAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "pca_batch":
                batch_results = []
                for idx, pca_config in enumerate(stage.config):
                    if pca_config.output_dir:
                        pca_config.output_dir.mkdir(parents=True, exist_ok=True)
                    analyzer = PCAAnalyzer(pca_config)
                    batch_results.append(analyzer.run_analysis())
                    logger.info(f"  ✓ Completed PCA batch {idx + 1}/{len(stage.config)}")
                results[name] = batch_results
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
            elif name == "distance_batch":
                batch_results = []
                for idx, dist_config in enumerate(stage.config):
                    dist_config.output_dir.mkdir(parents=True, exist_ok=True)
                    analyzer = DistanceAnalyzer(dist_config)
                    batch_results.append(analyzer.run_analysis())
                    logger.info(f"  ✓ Completed distance batch {idx + 1}/{len(stage.config)}")
                results[name] = batch_results
            elif name == "distance_compare":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = DistanceComparisonAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "interactions":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = InteractionAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "interactions_batch":
                batch_results = []
                for idx, inter_cfg in enumerate(stage.config):
                    inter_cfg.output_dir.mkdir(parents=True, exist_ok=True)
                    analyzer = InteractionAnalyzer(inter_cfg)
                    batch_results.append(analyzer.run_analysis())
                    logger.info(f"  ✓ Completed interactions batch {idx + 1}/{len(stage.config)}")
                results[name] = batch_results
            elif name == "interactions_compare":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = InteractionComparisonAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "prolif":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = ProlifAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "prolif_batch":
                batch_results = []
                for idx, cfg in enumerate(stage.config):
                    cfg.output_dir.mkdir(parents=True, exist_ok=True)
                    analyzer = ProlifAnalyzer(cfg)
                    batch_results.append(analyzer.run_analysis())
                    logger.info(f"  ✓ Completed prolif batch {idx + 1}/{len(stage.config)}")
                results[name] = batch_results
            elif name == "prolif_compare":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = ProlifComparisonAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "ttclust":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = TTClustAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "ttclust_batch":
                batch_results = []
                for idx, cfg in enumerate(stage.config):
                    cfg.output_dir.mkdir(parents=True, exist_ok=True)
                    analyzer = TTClustAnalyzer(cfg)
                    batch_results.append(analyzer.run_analysis())
                    logger.info(f"  ✓ Completed ttclust batch {idx + 1}/{len(stage.config)}")
                results[name] = batch_results
            elif name == "ranking":
                stage.config.output_dir.mkdir(parents=True, exist_ok=True)
                analyzer = RankingAnalyzer(stage.config)
                results[name] = analyzer.run_analysis()
            elif name == "utils":
                results[name] = run_utils_tasks(stage.config)
            else:
                logger.warning(f"Unknown stage '{name}' - skipping")
                continue

            if stage.cache_enabled and getattr(stage.config, "output_dir", None):
                output_dir = stage.config.output_dir
                files = _collect_stage_files(name, stage.config)
                signature = compute_signature(files, extra={"stage": name})
                save_cache(output_dir, signature, files)

            logger.info(f"Finished stage: {name}")

        except Exception as exc:
            logger.error(f"Stage '{name}' failed: {exc}")
            results[name] = {"success": False, "error": str(exc)}
            if not continue_on_error:
                raise

    return results


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
