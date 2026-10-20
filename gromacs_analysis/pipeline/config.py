"""
Pipeline configuration dataclasses.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclass
class PipelineStage:
    """Configuration wrapper for a single pipeline stage."""
    name: str
    enabled: bool = True
    config: Optional[Any] = None
    config_file: Optional[Path] = None
    inline: bool = True
    cache_enabled: bool = False


@dataclass
class PipelineConfig:
    """Top-level configuration for multi-stage analysis."""
    project_name: str = "GROMACS Pipeline"
    protein_name: Optional[str] = None
    data_root: Optional[Path] = None
    output_root: Optional[Path] = None
    run: List[str] = field(default_factory=list)
    stages: Dict[str, PipelineStage] = field(default_factory=dict)
    # YAML file the configuration was loaded from (provenance in the run report).
    source_path: Optional[Path] = None

    def ordered_stage_names(self) -> List[str]:
        """Return stage names in execution order."""
        if self.run:
            return self.run
        return list(self.stages.keys())

    def enabled_stages(self) -> List[PipelineStage]:
        """Return enabled stages with configs in execution order."""
        stages: List[PipelineStage] = []
        for name in self.ordered_stage_names():
            stage = self.stages.get(name)
            if not stage:
                continue
            if stage.enabled and stage.config is not None:
                stages.append(stage)
        return stages
