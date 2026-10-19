"""
PCA Analysis Module
==================

Pattern-based PCA parsing, processing, and visualization.
"""

from .config import (
    PCAConfig,
    PCASystemConfig,
    PCACompareConfig,
    PCAGromacsRunConfig,
    PCATerminalDetectionConfig,
    PCAJointConfig,
    PCAJointSystem,
)
from .parser import PCAParser
from .plotter import PCAPlotter
from .analyzer import PCAAnalyzer
from .compare import PCAComparisonAnalyzer
from .runner import PCAGromacsRunner
from .joint import JointPCAResult, joint_pca, rmsip
from .joint_analyzer import PCAJointAnalyzer

__all__ = [
    'PCAConfig',
    'PCASystemConfig',
    'PCACompareConfig',
    'PCAGromacsRunConfig',
    'PCATerminalDetectionConfig',
    'PCAParser',
    'PCAPlotter',
    'PCAAnalyzer',
    'PCAComparisonAnalyzer',
    'PCAGromacsRunner',
    'PCAJointConfig',
    'PCAJointSystem',
    'PCAJointAnalyzer',
    'JointPCAResult',
    'joint_pca',
    'rmsip',
]
