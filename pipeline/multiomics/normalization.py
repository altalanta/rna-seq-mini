"""
Cross-omics normalization and batch correction methods for multi-omics integration.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Union
from sklearn.preprocessing import StandardScaler, QuantileTransformer
from scipy import stats
import logging

logger = logging.getLogger(__name__)


class CrossOmicsNormalizer:
    """Handles normalization across different omics data types."""

    def __init__(self):
        self.normalization_params = {}
        self.batch_effects = {}

    def normalize_omics_data(self, omics_data: Dict, method: str = 'quantile') -> Dict:
        """
        Normalize multiple omics datasets using specified method.

        Args:
            omics_data: Dictionary of OmicsData objects
            method: Normalization method ('quantile', 'zscore', 'combat', 'limma')
        """
        normalized_data = {}

        for omics_type, data in omics_data.items():
            logger.info(f"Normalizing {omics_type.value} data using {method} method")

            if method == 'quantile':
                normalized = self._quantile_normalize(data.data)
            elif method == 'zscore':
                normalized = self._zscore_normalize(data.data)
            elif method == 'combat':
                normalized = self._combat_normalize(data.data, data.sample_metadata)
            elif method == 'limma':
                normalized = self._limma_normalize(data.data, data.sample_metadata)
            else:
                raise ValueError(f"Unknown normalization method: {method}")

            # Update data with normalized values
            data.data = normalized
            data.normalization_method = method
            normalized_data[omics_type] = data

        return normalized_data

    def _quantile_normalize(self, data: pd.DataFrame) -> pd.DataFrame:
        """Apply quantile normalization to make distributions comparable."""
        # Sort each column
        sorted_data = np.sort(data.values, axis=0)

        # Compute row means of sorted data
        row_means = np.mean(sorted_data, axis=1)

        # For each original value, find its rank and replace with corresponding row mean
        normalized = data.copy()
        for col in data.columns:
            sorted_indices = np.argsort(data[col].values)
            normalized.loc[data.index[sorted_indices], col] = row_means

        return normalized

    def _zscore_normalize(self, data: pd.DataFrame) -> pd.DataFrame:
        """Apply z-score normalization (standardization)."""
        scaler = StandardScaler()
        normalized_values = scaler.fit_transform(data.values)
        return pd.DataFrame(normalized_values, index=data.index, columns=data.columns)

    def _combat_normalize(self, data: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
        """
        Apply ComBat batch correction for cross-platform normalization.
        This is a simplified implementation - in practice, you'd use the ComBat R package.
        """
        # For now, implement a simple batch correction using linear regression
        normalized = data.copy()

        # Identify batch columns in metadata
        batch_columns = [col for col in metadata.columns if 'batch' in col.lower()]

        if not batch_columns:
            logger.warning("No batch columns found for ComBat normalization")
            return normalized

        for batch_col in batch_columns:
            if batch_col in metadata.columns:
                # Simple batch correction using group means
                for batch_value in metadata[batch_col].unique():
                    if pd.isna(batch_value):
                        continue

                    batch_mask = metadata[batch_col] == batch_value
                    if batch_mask.sum() > 1:
                        # Subtract batch mean for each feature
                        batch_samples = metadata[batch_mask].index
                        batch_mean = data[batch_samples].mean(axis=1)
                        normalized[batch_samples] = data[batch_samples].subtract(batch_mean, axis=0)

        return normalized

    def _limma_normalize(self, data: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
        """
        Apply LIMMA-style normalization for microarray-like data.
        This removes systematic biases between arrays/chips.
        """
        # Calculate normalization factors for each sample
        norm_factors = []

        for sample in data.columns:
            # Use median absolute deviation as a measure of scale
            sample_data = data[sample].dropna()
            if len(sample_data) > 0:
                # LIMMA uses median absolute deviation from median
                median = sample_data.median()
                mad = np.median(np.abs(sample_data - median))
                # Avoid division by zero
                norm_factor = mad if mad > 0 else 1.0
            else:
                norm_factor = 1.0

            norm_factors.append(norm_factor)

        # Apply normalization factors
        norm_factors = np.array(norm_factors)
        normalized = data.div(norm_factors, axis=1)

        return normalized

    def harmonize_scales(self, omics_data: Dict) -> Dict:
        """
        Harmonize scales across different omics types for better comparability.

        This ensures that different omics types contribute equally to downstream analyses.
        """
        harmonized_data = {}

        # Calculate scale factors for each omics type
        scale_factors = {}
        for omics_type, data in omics_data.items():
            # Use median absolute deviation as scale measure
            mad = np.median(np.abs(data.data.values - np.median(data.data.values)))
            scale_factors[omics_type] = 1.0 / mad if mad > 0 else 1.0

        # Apply scale harmonization
        for omics_type, data in omics_data.items():
            scaled_data = data.data * scale_factors[omics_type]

            # Update data object
            data.data = scaled_data
            harmonized_data[omics_type] = data

        return harmonized_data


class BatchCorrector:
    """Handles batch effect correction across omics datasets."""

    def __init__(self):
        self.batch_models = {}

    def detect_batch_effects(self, data: pd.DataFrame, metadata: pd.DataFrame) -> Dict[str, float]:
        """
        Detect potential batch effects in the data.

        Returns dictionary with batch variables and their effect sizes.
        """
        batch_effects = {}

        # Look for batch-related columns in metadata
        potential_batch_cols = [col for col in metadata.columns
                              if any(term in col.lower() for term in ['batch', 'plate', 'run', 'date'])]

        for batch_col in potential_batch_cols:
            if metadata[batch_col].nunique() > 1:
                # Simple ANOVA-like test for batch effects
                try:
                    # For each feature, test if batch explains variance
                    batch_effect_sizes = []

                    for feature in data.index[:1000]:  # Sample first 1000 features for speed
                        if feature in data.index:
                            feature_data = data.loc[feature].dropna()

                            if len(feature_data) > 10:  # Need minimum samples
                                # Simple F-statistic calculation
                                groups = [feature_data[metadata[batch_col] == val] for val in metadata[batch_col].unique()]
                                groups = [g for g in groups if len(g) > 1]  # Remove single-sample groups

                                if len(groups) > 1:
                                    f_stat, p_val = stats.f_oneway(*groups)
                                    if p_val < 0.05:  # Significant batch effect
                                        batch_effect_sizes.append(abs(f_stat))

                    if batch_effect_sizes:
                        batch_effects[batch_col] = np.mean(batch_effect_sizes)

                except Exception as e:
                    logger.warning(f"Error detecting batch effects for {batch_col}: {e}")

        return batch_effects

    def correct_batch_effects(self, data: pd.DataFrame, metadata: pd.DataFrame,
                            batch_columns: List[str] = None) -> pd.DataFrame:
        """
        Correct for batch effects using linear modeling.

        Args:
            data: Expression data (features x samples)
            metadata: Sample metadata
            batch_columns: List of batch variables to correct for
        """
        if batch_columns is None:
            batch_columns = [col for col in metadata.columns
                           if any(term in col.lower() for term in ['batch', 'plate', 'run'])]

        corrected_data = data.copy()

        for feature in data.index:
            feature_data = data.loc[feature].dropna()

            if len(feature_data) < 10:  # Skip features with insufficient data
                continue

            # Build linear model: expression ~ batch + condition
            try:
                # Create design matrix
                design_df = metadata.loc[feature_data.index].copy()

                # Add condition as covariate
                if 'condition' in design_df.columns:
                    design_df['condition_encoded'] = pd.Categorical(design_df['condition']).codes

                # Fit linear model and extract residuals (batch-corrected values)
                batch_corrected = self._fit_batch_correction_model(feature_data, design_df, batch_columns)

                if batch_corrected is not None:
                    corrected_data.loc[feature, feature_data.index] = batch_corrected

            except Exception as e:
                logger.warning(f"Error correcting batch effects for {feature}: {e}")

        return corrected_data

    def _fit_batch_correction_model(self, feature_data: pd.Series, design_df: pd.DataFrame,
                                  batch_columns: List[str]) -> Optional[np.ndarray]:
        """Fit linear model for batch correction."""
        try:
            import statsmodels.api as sm

            # Create design matrix
            X = pd.DataFrame(index=feature_data.index)

            # Add batch variables (one-hot encoded)
            for batch_col in batch_columns:
                if batch_col in design_df.columns:
                    batch_dummies = pd.get_dummies(design_df[batch_col], prefix=batch_col)
                    X = pd.concat([X, batch_dummies], axis=1)

            # Add condition as covariate
            if 'condition_encoded' in design_df.columns:
                X['condition'] = design_df['condition_encoded']

            # Add intercept
            X = sm.add_constant(X)

            # Fit model and get residuals
            model = sm.OLS(feature_data.values, X.values)
            results = model.fit()

            # Return residuals (batch-corrected values)
            return results.resid

        except ImportError:
            logger.warning("statsmodels not available for batch correction")
            return None
        except Exception as e:
            logger.warning(f"Error in batch correction model: {e}")
            return None


class MultiOmicsAligner:
    """Aligns multiple omics datasets for integrated analysis."""

    def __init__(self):
        self.alignment_matrices = {}

    def align_datasets(self, datasets: Dict, alignment_method: str = 'correlation') -> Dict:
        """
        Align multiple omics datasets for integrated analysis.

        Args:
            datasets: Dictionary of normalized omics data
            alignment_method: Method for aligning datasets ('correlation', 'mutual_information')
        """
        aligned_datasets = {}

        if alignment_method == 'correlation':
            aligned_datasets = self._align_by_correlation(datasets)
        elif alignment_method == 'mutual_information':
            aligned_datasets = self._align_by_mutual_information(datasets)
        else:
            raise ValueError(f"Unknown alignment method: {alignment_method}")

        return aligned_datasets

    def _align_by_correlation(self, datasets: Dict) -> Dict:
        """Align datasets by maximizing correlations between shared features."""
        # For simplicity, return datasets as-is
        # In a more sophisticated implementation, this would:
        # 1. Identify overlapping features across omics types
        # 2. Compute correlation matrices
        # 3. Apply rotation/transformation to align datasets
        return datasets

    def _align_by_mutual_information(self, datasets: Dict) -> Dict:
        """Align datasets by maximizing mutual information."""
        # Placeholder for mutual information alignment
        return datasets
