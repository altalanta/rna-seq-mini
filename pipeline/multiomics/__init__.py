"""
Multi-omics integration framework for RNASEQ-MINI.
Provides unified analysis across multiple omics data types.
"""

from .omics_types import OmicsType, OmicsData
from .normalization import CrossOmicsNormalizer

# Note: The following modules are planned but not yet implemented:
# - integration (MultiOmicsIntegrator)
# - visualization (MultiOmicsVisualizer)  
# - statistics (JointStatisticalTester)

__all__ = [
    'OmicsType',
    'OmicsData',
    'CrossOmicsNormalizer',
]
