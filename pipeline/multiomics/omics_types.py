"""
Defines omics data types and their characteristics for multi-omics integration.
"""

from enum import Enum
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, field
import pandas as pd
import numpy as np


class OmicsType(Enum):
    """Enumeration of supported omics data types."""
    RNASEQ = "rnaseq"
    ATACSEQ = "atacseq"
    CHIPSEQ = "chipseq"
    METHYLATION = "methylation"
    PROTEOMICS = "proteomics"
    METABOLOMICS = "metabolomics"
    MICROBIOME = "microbiome"


@dataclass
class OmicsData:
    """Container for omics data with metadata and processing information."""

    omics_type: OmicsType
    data: pd.DataFrame
    sample_metadata: pd.DataFrame
    feature_metadata: pd.DataFrame
    normalization_method: Optional[str] = None
    batch_correction_method: Optional[str] = None
    quality_metrics: Dict[str, float] = field(default_factory=dict)

    def __post_init__(self):
        """Validate data structure and compute basic metrics."""
        self._validate_data_structure()
        self._compute_basic_metrics()

    def _validate_data_structure(self):
        """Validate that the data structure is consistent."""
        # Check that samples match between data and metadata
        data_samples = set(self.data.columns)
        metadata_samples = set(self.sample_metadata.index)

        if not data_samples.issubset(metadata_samples):
            missing_samples = data_samples - metadata_samples
            raise ValueError(f"Samples in data but not in metadata: {missing_samples}")

        # Check for required metadata columns
        required_columns = ['sample_id', 'condition']
        missing_columns = [col for col in required_columns if col not in self.sample_metadata.columns]
        if missing_columns:
            raise ValueError(f"Missing required metadata columns: {missing_columns}")

    def _compute_basic_metrics(self):
        """Compute basic quality metrics."""
        self.quality_metrics.update({
            'n_samples': len(self.data.columns),
            'n_features': len(self.data.index),
            'missing_rate': self.data.isnull().sum().sum() / (self.data.shape[0] * self.data.shape[1]),
            'mean_expression': self.data.mean().mean(),
            'median_expression': self.data.median().median()
        })


@dataclass
class MultiOmicsDataset:
    """Container for multiple omics datasets with integrated metadata."""

    datasets: Dict[OmicsType, OmicsData]
    integrated_metadata: pd.DataFrame
    integration_method: str = "inner"
    harmonized_features: Optional[pd.DataFrame] = None

    def __post_init__(self):
        """Validate and harmonize the multi-omics dataset."""
        self._validate_compatibility()
        self._harmonize_metadata()

    def _validate_compatibility(self):
        """Validate that datasets are compatible for integration."""
        sample_sets = []
        for omics_type, dataset in self.datasets.items():
            sample_sets.append(set(dataset.sample_metadata.index))

        # Check sample overlap
        common_samples = set.intersection(*sample_sets) if sample_sets else set()

        if not common_samples:
            raise ValueError("No common samples found across omics datasets")

        # Check for consistent metadata structure
        required_columns = ['sample_id', 'condition']
        for omics_type, dataset in self.datasets.items():
            missing_cols = [col for col in required_columns if col not in dataset.sample_metadata.columns]
            if missing_cols:
                raise ValueError(f"Dataset {omics_type.value} missing columns: {missing_cols}")

    def _harmonize_metadata(self):
        """Harmonize sample metadata across all datasets."""
        # Find common samples
        sample_sets = [set(dataset.sample_metadata.index) for dataset in self.datasets.values()]
        common_samples = set.intersection(*sample_sets)

        # Create integrated metadata from common samples
        metadata_frames = []
        for omics_type, dataset in self.datasets.items():
            # Add omics type identifier to avoid column conflicts
            metadata = dataset.sample_metadata.loc[common_samples].copy()
            metadata['omics_type'] = omics_type.value
            metadata_frames.append(metadata)

        self.integrated_metadata = pd.concat(metadata_frames, axis=0)

        # Create harmonized features if data dimensions allow
        self._harmonize_features()

    def _harmonize_features(self):
        """Harmonize feature names across omics types where possible."""
        # For now, just store original feature metadata
        # In a more sophisticated implementation, this could map genes to common identifiers
        feature_frames = []
        for omics_type, dataset in self.datasets.items():
            features = dataset.feature_metadata.copy()
            features['omics_type'] = omics_type.value
            feature_frames.append(features)

        self.harmonized_features = pd.concat(feature_frames, axis=0) if feature_frames else None

    def get_common_samples(self) -> List[str]:
        """Get list of samples common across all omics datasets."""
        sample_sets = [set(dataset.sample_metadata.index) for dataset in self.datasets.values()]
        return list(set.intersection(*sample_sets))

    def get_omics_types(self) -> List[OmicsType]:
        """Get list of omics types in this dataset."""
        return list(self.datasets.keys())

    def subset_to_common_samples(self) -> 'MultiOmicsDataset':
        """Return new dataset containing only common samples."""
        common_samples = self.get_common_samples()

        subset_datasets = {}
        for omics_type, dataset in self.datasets.items():
            # Filter data to common samples
            data_subset = dataset.data.loc[:, common_samples]

            # Filter metadata to common samples
            metadata_subset = dataset.sample_metadata.loc[common_samples]

            subset_datasets[omics_type] = OmicsData(
                omics_type=omics_type,
                data=data_subset,
                sample_metadata=metadata_subset,
                feature_metadata=dataset.feature_metadata,
                normalization_method=dataset.normalization_method,
                batch_correction_method=dataset.batch_correction_method,
                quality_metrics=dataset.quality_metrics
            )

        return MultiOmicsDataset(
            datasets=subset_datasets,
            integrated_metadata=self.integrated_metadata.loc[common_samples],
            integration_method=self.integration_method
        )




