"""
Pipeline CLI Entry Point
========================

Run a multi-stage analysis pipeline (MD -> MMPBSA -> PCA) from one YAML config.
"""

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import List, Optional

from ..pipeline import (
    PipelineConfig,
    PipelineStage,
    load_pipeline_config,
    save_pipeline_config,
    generate_pipeline_template,
    run_pipeline,
)
from ..cli.interactive import interactive_prompt as md_interactive_prompt
from ..cli.mmpbsa_cli import interactive_prompt as mmpbsa_interactive_prompt
from ..cli.pca_cli import interactive_prompt as pca_interactive_prompt
from ..config.yaml_loader import load_yaml_config as load_md_yaml
from ..mmpbsa.yaml_loader import load_yaml_config as load_mmpbsa_yaml
from ..pca.yaml_loader import load_yaml_config as load_pca_yaml

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

VALID_STAGES = [
    "qc",
    "md",
    "md_batch",
    "md_compare",
    "mmpbsa",
    "pca",
    "pca_batch",
    "pca_compare",
    "distance",
    "distance_batch",
    "distance_compare",
    "interactions",
    "interactions_batch",
    "interactions_compare",
    "prolif",
    "prolif_batch",
    "prolif_compare",
    "ttclust",
    "ttclust_batch",
    "ranking",
    "utils",
]


def _parse_stage_list(raw: str) -> List[str]:
    if not raw:
        return []
    tokens = [t.strip().lower() for t in re.split(r"[,\s]+", raw) if t.strip()]
    if "all" in tokens:
        return VALID_STAGES[:]
    stages = [t for t in tokens if t in VALID_STAGES]
    return stages


def _resolve_config_paths(config, base_path: Path) -> None:
    for attr in ("base_dir", "output_dir", "apo_base_dir"):
        if hasattr(config, attr):
            current = getattr(config, attr)
            if isinstance(current, Path) and not current.is_absolute():
                setattr(config, attr, (base_path / current).resolve())


def _interactive_pipeline() -> PipelineConfig:
    print("\n" + "=" * 70)
    print(" " * 18 + "GROMACS ANALYSIS PIPELINE")
    print(" " * 18 + "Interactive Configuration")
    print("=" * 70 + "\n")

    project_name = input("Project name [GROMACS Pipeline]: ").strip() or "GROMACS Pipeline"
    protein_name = input("Protein name (shared default, optional): ").strip() or None

    data_root_input = input("Data root (optional): ").strip()
    data_root = Path(data_root_input).expanduser().resolve() if data_root_input else None

    output_root_input = input("Output root [./analysis_outputs]: ").strip() or "./analysis_outputs"
    output_root = Path(output_root_input).expanduser().resolve() if output_root_input else None

    print("\nAvailable stages:")
    print("  1. qc       - Scientific preflight checks")
    print("  2. md       - MD trajectory analysis")
    print("  3. md_batch - Holo vs Apo batches (one output per holo)")
    print("  4. md_compare - Cross-ligand comparison (overlay all holos)")
    print("  5. mmpbsa   - Binding energy analysis")
    print("  6. pca      - Principal component analysis")
    print("  7. pca_batch- PCA per system (one output per PCA dir)")
    print("  8. pca_compare - Cross-ligand PCA comparison")
    print("  9. distance - Ligand-residue distance analysis")
    print("  10. distance_batch - Distance per system (one output per system)")
    print("  11. distance_compare - Cross-ligand distance comparison")
    print("  12. interactions - Interaction fingerprint analysis")
    print("  13. interactions_batch - Interactions per system (one output per system)")
    print("  14. interactions_compare - Cross-ligand interaction comparison")
    print("  15. prolif - ProLIF fingerprint analysis")
    print("  16. prolif_batch - ProLIF per system (one output per system)")
    print("  17. prolif_compare - Cross-ligand ProLIF comparison")
    print("  18. ttclust - TTClust trajectory clustering")
    print("  19. ttclust_batch - TTClust per system (one output per system)")
    print("  20. ranking - Cross-ligand ranking table")
    print("  21. utils    - Utility tasks (CSV inspect, DAT convert, renumber, modevectors)")
    stage_input = input("\nStages to include (comma-separated) [qc,md,mmpbsa,pca]: ").strip()
    stage_names = _parse_stage_list(stage_input) or VALID_STAGES[:]

    pipeline = PipelineConfig(
        project_name=project_name,
        protein_name=protein_name,
        data_root=data_root,
        output_root=output_root,
        run=stage_names,
        stages={},
    )

    for stage_name in stage_names:
        print("\n" + "-" * 70)
        print(f"Configure stage: {stage_name}")
        print("-" * 70)

        use_existing = input(f"Use existing {stage_name} YAML config? (y/n) [n]: ").strip().lower()
        if use_existing in ("y", "yes"):
            path_input = input("Path to YAML config: ").strip()
            if not path_input:
                print("  ✗ Config path is required.")
                sys.exit(1)
            config_path = Path(path_input).expanduser().resolve()
            if not config_path.exists():
                print(f"  ✗ Config file not found: {config_path}")
                sys.exit(1)

            if stage_name == "md":
                config = load_md_yaml(config_path)
            elif stage_name == "mmpbsa":
                config = load_mmpbsa_yaml(config_path)
            elif stage_name == "pca":
                config = load_pca_yaml(config_path)
            else:
                print(f"  ⚠ YAML loading for '{stage_name}' is not supported here.")
                print("    Add it directly in the pipeline YAML if needed.")
                continue

            _resolve_config_paths(config, config_path.parent)

            if output_root:
                use_root = input(f"Use output_root for {stage_name}? (y/n) [n]: ").strip().lower()
                if use_root in ("y", "yes"):
                    config.output_dir = output_root / stage_name

            pipeline.stages[stage_name] = PipelineStage(
                name=stage_name,
                enabled=True,
                config=config,
                config_file=config_path,
                inline=False,
            )
            continue

        default_output = output_root / stage_name if output_root else None
        if stage_name == "md":
            config = md_interactive_prompt(
                default_base_dir=data_root,
                default_output_dir=default_output,
                default_protein_name=protein_name,
                save_yaml=False,
            )
        elif stage_name == "mmpbsa":
            config = mmpbsa_interactive_prompt(
                default_base_dir=data_root,
                default_output_dir=default_output,
                default_protein_name=protein_name,
                save_yaml=False,
            )
        elif stage_name == "pca":
            config = pca_interactive_prompt(
                default_base_dir=data_root,
                default_output_dir=default_output,
                default_protein_name=protein_name,
            )
        else:
            print(f"  ⚠ Interactive setup for '{stage_name}' is not supported yet.")
            print("    Add it directly in the pipeline YAML if needed.")
            continue

        pipeline.stages[stage_name] = PipelineStage(
            name=stage_name,
            enabled=True,
            config=config,
            inline=True,
        )

    print("\n" + "=" * 70)
    print("Pipeline configuration complete.")
    print("=" * 70 + "\n")

    save = input("Save this pipeline configuration to YAML? (y/n) [y]: ").strip().lower()
    if save in ("", "y", "yes"):
        output_path = Path(f"{project_name.lower().replace(' ', '_')}_pipeline.yaml")
        save_pipeline_config(pipeline, output_path)
        print(f"  ✓ Saved: {output_path}")
    else:
        print("  ℹ Not saved. You can still run the pipeline now.")

    return pipeline


def _print_summary(config: PipelineConfig) -> None:
    print("\nPipeline Summary")
    print("-" * 70)
    print(f"Project: {config.project_name}")
    if config.protein_name:
        print(f"Protein: {config.protein_name}")
    if config.data_root:
        print(f"Data root: {config.data_root}")
    if config.output_root:
        print(f"Output root: {config.output_root}")
    print(f"Run order: {', '.join(config.ordered_stage_names())}")
    for stage in config.ordered_stage_names():
        stage_cfg = config.stages.get(stage)
        if not stage_cfg:
            continue
        status = "enabled" if stage_cfg.enabled else "disabled"
        if isinstance(stage_cfg.config, list):
            out_dir = [str(getattr(cfg, "output_dir", "")) for cfg in stage_cfg.config]
        elif isinstance(stage_cfg.config, dict):
            out_dir = stage_cfg.config.get("output_dir")
        else:
            out_dir = getattr(stage_cfg.config, "output_dir", None) if stage_cfg.config else None
        print(f"  - {stage}: {status} | output_dir={out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run GROMACS analyses from a single pipeline config",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate a pipeline template
  gromacs-pipeline --generate-template

  # Run with pipeline config
  gromacs-pipeline --config project_pipeline.yaml

  # Interactive wizard
  gromacs-pipeline --interactive

  # No YAML by hand: discover systems (local, SSH or condorctl), name them, run
  clarity-setup tui /path/to/simulations
        """,
    )

    config_group = parser.add_mutually_exclusive_group()
    config_group.add_argument("--config", "-c", type=Path, help="Path to pipeline YAML config")
    config_group.add_argument("--interactive", "-i", action="store_true", help="Interactive pipeline wizard")
    config_group.add_argument("--generate-template", "-g", action="store_true", help="Generate pipeline YAML template")

    parser.add_argument(
        "--run",
        help=(
            "Comma-separated stages to run (qc,md,md_batch,md_compare,mmpbsa,pca,"
            "pca_batch,pca_compare,distance,distance_batch,distance_compare,"
            "interactions,interactions_batch,interactions_compare,"
            "prolif,prolif_batch,prolif_compare,ttclust,ttclust_batch,ranking,utils)"
        ),
    )
    parser.add_argument("--skip", help="Comma-separated stages to skip")
    parser.add_argument("--dry-run", action="store_true", help="Validate and print summary only")
    parser.add_argument("--continue-on-error", action="store_true", help="Continue even if a stage fails")
    parser.add_argument("--no-cache", action="store_true",
                        help="Ignore stage caches (cache: true) and recompute everything")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.generate_template:
        output_path = Path("gromacs_pipeline_template.yaml")
        generate_pipeline_template(output_path)
        print(f"✓ Template created: {output_path}")
        return

    if args.interactive:
        pipeline_config = _interactive_pipeline()
    elif args.config:
        if not args.config.exists():
            print(f"Error: Config file not found: {args.config}")
            sys.exit(1)
        pipeline_config = load_pipeline_config(args.config)
    else:
        parser.print_help()
        sys.exit(1)

    # Apply run/skip overrides
    if args.run:
        pipeline_config.run = _parse_stage_list(args.run)
    if args.skip:
        for stage in _parse_stage_list(args.skip):
            if stage in pipeline_config.stages:
                pipeline_config.stages[stage].enabled = False

    _print_summary(pipeline_config)

    if args.dry_run:
        print("\nDry run complete. No analysis executed.")
        return

    try:
        run_pipeline(pipeline_config, continue_on_error=args.continue_on_error, use_cache=not args.no_cache)
    except Exception:
        sys.exit(1)


if __name__ == "__main__":
    main()
