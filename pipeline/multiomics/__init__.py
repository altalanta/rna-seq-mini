"""
Multi-omics integration framework for RNASEQ-MINI.
Provides unified analysis across multiple omics data types.
"""

from .omics_types import OmicsType, OmicsData
from .integration import MultiOmicsIntegrator
from .normalization import CrossOmicsNormalizer
from .visualization import MultiOmicsVisualizer
from .statistics import JointStatisticalTester

__all__ = [
    'OmicsType',
    'OmicsData',
    'MultiOmicsIntegrator',
    'CrossOmicsNormalizer',
    'MultiOmicsVisualizer',
    'JointStatisticalTester'
]

