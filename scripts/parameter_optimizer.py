#!/usr/bin/env python3
"""
Automated parameter optimization system for RNASEQ-MINI pipeline.
Uses ML/AI techniques to analyze data characteristics and recommend optimal parameters.
"""

import json
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import subprocess
import gzip
import statistics
from collections import Counter
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Try to import ML libraries, but don't fail if not available
try:
    from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
    from sklearn.model_selection import train_test_split, cross_val_score
    from sklearn.metrics import mean_squared_error, r2_score
    from sklearn.preprocessing import StandardScaler
    import warnings
    SKLEARN_AVAILABLE = True
    # Suppress sklearn warnings for cleaner output
    warnings.filterwarnings('ignore')
except ImportError:
    SKLEARN_AVAILABLE = False
    logger.warning("scikit-learn not available, using rule-based optimization only")


class FastqAnalyzer:
    """Analyzes FASTQ files to extract relevant features for parameter optimization."""

    def __init__(self, max_reads: int = 10000):
        self.max_reads = max_reads

    def analyze_fastq(self, fastq_path: Path) -> Dict[str, Any]:
        """
        Extract comprehensive features from a FASTQ file.
        """
        features = {
            'file_path': str(fastq_path),
            'file_size': fastq_path.stat().st_size,
            'is_gzipped': fastq_path.suffix == '.gz'
        }

        try:
            # Open file (handle gzipped files)
            opener = gzip.open if features['is_gzipped'] else open
            with opener(fastq_path, 'rt') as f:
                read_lengths = []
                quality_scores = []
                gc_content = []
                base_counts = Counter()
                read_count = 0

                while read_count < self.max_reads:
                    try:
                        # FASTQ format: 4 lines per read
                        header = f.readline().strip()
                        if not header:
                            break

                        sequence = f.readline().strip()
                        plus = f.readline().strip()
                        quality = f.readline().strip()

                        if not (header and sequence and plus and quality):
                            break

                        # Extract features from this read
                        read_lengths.append(len(sequence))
                        gc_content.append(self._calculate_gc_content(sequence))
                        quality_scores.extend([ord(q) - 33 for q in quality])  # Phred+33
                        base_counts.update(sequence)

                        read_count += 1

                    except Exception as e:
                        logger.warning(f"Error reading read {read_count}: {e}")
                        break

            # Calculate statistics
            features.update({
                'total_reads_analyzed': read_count,
                'avg_read_length': statistics.mean(read_lengths) if read_lengths else 0,
                'median_read_length': statistics.median(read_lengths) if read_lengths else 0,
                'read_length_std': statistics.stdev(read_lengths) if len(read_lengths) > 1 else 0,
                'avg_quality_score': statistics.mean(quality_scores) if quality_scores else 0,
                'median_quality_score': statistics.median(quality_scores) if quality_scores else 0,
                'quality_score_std': statistics.stdev(quality_scores) if len(quality_scores) > 1 else 0,
                'gc_content': statistics.mean(gc_content) if gc_content else 0,
                'base_composition': dict(base_counts),
                'total_bases': sum(base_counts.values())
            })

            # Calculate derived metrics
            features.update(self._calculate_derived_metrics(features))

        except Exception as e:
            logger.error(f"Error analyzing {fastq_path}: {e}")
            features['error'] = str(e)

        return features

    def _calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content percentage."""
        if not sequence:
            return 0.0
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100

    def _calculate_derived_metrics(self, features: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate derived metrics from basic features."""
        derived = {}

        # Quality distribution metrics
        if 'avg_quality_score' in features and 'quality_score_std' in features:
            derived['quality_variability'] = features['quality_score_std'] / max(features['avg_quality_score'], 1)

        # Read length uniformity
        if 'read_length_std' in features and 'avg_read_length' in features:
            derived['read_length_cv'] = features['read_length_std'] / max(features['avg_read_length'], 1)

        # Base composition balance
        if 'base_composition' in features and features['total_bases'] > 0:
            base_comp = features['base_composition']
            total = features['total_bases']
            at_content = (base_comp.get('A', 0) + base_comp.get('T', 0)) / total
            gc_content = (base_comp.get('G', 0) + base_comp.get('C', 0)) / total
            derived['at_gc_ratio'] = at_content / max(gc_content, 0.01)

        # File characteristics
        if 'file_size' in features and 'total_reads_analyzed' in features:
            derived['bytes_per_read'] = features['file_size'] / max(features['total_reads_analyzed'], 1)

        return derived


class ParameterOptimizer:
    """ML-based parameter optimization for RNA-seq pipeline."""

    def __init__(self, model_dir: str = "models"):
        self.model_dir = Path(model_dir)
        self.model_dir.mkdir(exist_ok=True)
        self.models = {}
        self.scalers = {}
        self.feature_names = [
            'avg_read_length', 'read_length_std', 'avg_quality_score', 'quality_score_std',
            'gc_content', 'at_gc_ratio', 'quality_variability', 'read_length_cv', 'bytes_per_read'
        ]
        self.use_ml = SKLEARN_AVAILABLE

    def extract_features(self, fastq_files: List[Path]) -> pd.DataFrame:
        """
        Extract features from multiple FASTQ files for parameter optimization.
        """
        analyzer = FastqAnalyzer()
        all_features = []

        for fastq_file in fastq_files:
            logger.info(f"Analyzing {fastq_file}")
            features = analyzer.analyze_fastq(fastq_file)

            # Extract numerical features for ML
            numerical_features = {}
            for feature in self.feature_names:
                if feature in features:
                    numerical_features[feature] = features[feature]
                else:
                    numerical_features[feature] = 0.0  # Default for missing features

            numerical_features['file_path'] = str(fastq_file)
            all_features.append(numerical_features)

        return pd.DataFrame(all_features)

    def train_models(self, training_data_path: Optional[str] = None) -> Dict[str, Any]:
        """
        Train ML models for parameter prediction using historical data or synthetic training data.
        """
        if not self.use_ml:
            logger.warning("scikit-learn not available, cannot train ML models")
            return {'models_trained': 0, 'training_samples': 0, 'features_used': 0}

        logger.info("Training parameter optimization models...")

        # If no training data provided, create synthetic training data
        if training_data_path is None:
            logger.info("No training data provided, creating synthetic training data...")
            training_data = self._create_synthetic_training_data()
        else:
            training_data = pd.read_csv(training_data_path)

        # Prepare features and targets
        X = training_data[self.feature_names]
        y_salmon_threads = training_data['optimal_salmon_threads']
        y_salmon_libtype = training_data['optimal_salmon_libtype'].astype('category').cat.codes
        y_deseq2_alpha = training_data['optimal_deseq2_alpha']

        # Scale features
        self.scalers['salmon_threads'] = StandardScaler()
        self.scalers['salmon_libtype'] = StandardScaler()
        self.scalers['deseq2_alpha'] = StandardScaler()

        X_scaled = self.scalers['salmon_threads'].fit_transform(X)

        # Train models for different parameters
        models = {}

        # Salmon threads model
        models['salmon_threads'] = RandomForestRegressor(
            n_estimators=100, random_state=42, max_depth=10
        )
        models['salmon_threads'].fit(X_scaled, y_salmon_threads)

        # Salmon library type model (classification)
        models['salmon_libtype'] = GradientBoostingRegressor(
            n_estimators=100, random_state=42, max_depth=8
        )
        models['salmon_libtype'].fit(X_scaled, y_salmon_libtype)

        # DESeq2 alpha model
        models['deseq2_alpha'] = GradientBoostingRegressor(
            n_estimators=100, random_state=42, max_depth=6
        )
        models['deseq2_alpha'].fit(X_scaled, y_deseq2_alpha)

        # Save models
        for name, model in models.items():
            model_path = self.model_dir / f"{name}_model.pkl"
            scaler_path = self.model_dir / f"{name}_scaler.pkl"

            with open(model_path, 'wb') as f:
                pickle.dump(model, f)

            with open(scaler_path, 'wb') as f:
                pickle.dump(self.scalers[name], f)

        self.models = models
        logger.info(f"Trained and saved {len(models)} models")

        return {
            'models_trained': len(models),
            'training_samples': len(training_data),
            'features_used': len(self.feature_names)
        }

    def _create_synthetic_training_data(self, n_samples: int = 1000) -> pd.DataFrame:
        """
        Create synthetic training data based on typical RNA-seq characteristics.
        """
        np.random.seed(42)

        data = []

        # Realistic ranges for RNA-seq data
        read_length_means = np.random.normal(100, 20, n_samples)  # 80-120bp typical
        quality_means = np.random.normal(35, 3, n_samples)  # Phred 30-40 typical
        gc_contents = np.random.normal(50, 10, n_samples)  # 40-60% GC typical

        for i in range(n_samples):
            # Create correlated features
            read_len = max(50, min(200, read_length_means[i]))
            quality = max(20, min(45, quality_means[i]))
            gc = max(30, min(70, gc_contents[i]))

            # Add some noise to create realistic variation
            read_len_std = read_len * np.random.uniform(0.05, 0.2)
            quality_std = quality * np.random.uniform(0.1, 0.3)
            gc_std = gc * np.random.uniform(0.05, 0.15)

            # Calculate derived metrics
            quality_variability = quality_std / quality
            read_length_cv = read_len_std / read_len
            at_gc_ratio = (100 - gc) / max(gc, 1)
            bytes_per_read = read_len * 2.5  # Rough estimate

            # Generate optimal parameters based on data characteristics
            optimal_salmon_threads = self._recommend_salmon_threads(read_len, quality)
            optimal_salmon_libtype = self._recommend_salmon_libtype(read_len, gc)
            optimal_deseq2_alpha = self._recommend_deseq2_alpha(quality, read_len)

            data.append({
                'avg_read_length': read_len,
                'read_length_std': read_len_std,
                'avg_quality_score': quality,
                'quality_score_std': quality_std,
                'gc_content': gc,
                'at_gc_ratio': at_gc_ratio,
                'quality_variability': quality_variability,
                'read_length_cv': read_length_cv,
                'bytes_per_read': bytes_per_read,
                'optimal_salmon_threads': optimal_salmon_threads,
                'optimal_salmon_libtype': optimal_salmon_libtype,
                'optimal_deseq2_alpha': optimal_deseq2_alpha
            })

        return pd.DataFrame(data)

    def _recommend_salmon_threads(self, read_length: float, quality: float) -> int:
        """Recommend optimal Salmon threads based on data characteristics."""
        base_threads = 4

        # Longer reads typically need more threads for alignment
        if read_length > 120:
            base_threads += 2
        elif read_length < 80:
            base_threads -= 1

        # Higher quality allows for more parallelization
        if quality > 35:
            base_threads += 1

        return max(1, min(16, base_threads))

    def _recommend_salmon_libtype(self, read_length: float, gc_content: float) -> str:
        """Recommend optimal Salmon library type."""
        # Longer reads often work better with more specific library types
        if read_length > 100:
            return 'A'  # Automatic detection
        elif gc_content > 55:
            return 'ISR'  # Inward stranded
        else:
            return 'U'  # Unstranded

    def _recommend_deseq2_alpha(self, quality: float, read_length: float) -> float:
        """Recommend optimal DESeq2 alpha threshold."""
        base_alpha = 0.05

        # Higher quality data can use more stringent thresholds
        if quality > 35:
            base_alpha *= 0.8

        # Longer reads provide more power for DE detection
        if read_length > 100:
            base_alpha *= 0.9

        return max(0.01, min(0.1, base_alpha))

    def load_models(self) -> bool:
        """Load pre-trained models from disk."""
        success = True

        for model_name in ['salmon_threads', 'salmon_libtype', 'deseq2_alpha']:
            model_path = self.model_dir / f"{model_name}_model.pkl"
            scaler_path = self.model_dir / f"{model_name}_scaler.pkl"

            if model_path.exists() and scaler_path.exists():
                try:
                    with open(model_path, 'rb') as f:
                        self.models[model_name] = pickle.load(f)

                    with open(scaler_path, 'rb') as f:
                        self.scalers[model_name] = pickle.load(f)

                except Exception as e:
                    logger.error(f"Error loading {model_name} model: {e}")
                    success = False
            else:
                logger.warning(f"Model files not found for {model_name}")
                success = False

        return success and len(self.models) > 0

    def predict_optimal_parameters(self, features_df: pd.DataFrame) -> Dict[str, Any]:
        """
        Predict optimal parameters for given data features.
        """
        if self.use_ml and not self.models:
            if not self.load_models():
                logger.warning("No trained models available, using rule-based optimization")
                return self._predict_parameters_rule_based(features_df)

        if self.use_ml and self.models:
            # Use ML models if available
            return self._predict_parameters_ml(features_df)
        else:
            # Use rule-based optimization
            return self._predict_parameters_rule_based(features_df)

    def _predict_parameters_ml(self, features_df: pd.DataFrame) -> Dict[str, Any]:
        """Use ML models for parameter prediction."""
        predictions = {}

        for model_name in ['salmon_threads', 'salmon_libtype', 'deseq2_alpha']:
            if model_name in self.models:
                try:
                    # Scale features
                    X = features_df[self.feature_names]
                    X_scaled = self.scalers[model_name].transform(X)

                    # Make prediction
                    prediction = self.models[model_name].predict(X_scaled)

                    if model_name == 'salmon_libtype':
                        # Convert numeric prediction back to library type
                        libtype_map = {0: 'U', 1: 'ISR', 2: 'A'}
                        prediction = [libtype_map.get(int(p), 'A') for p in prediction]

                    predictions[model_name] = prediction[0] if len(prediction) == 1 else prediction

                except Exception as e:
                    logger.error(f"Error predicting {model_name}: {e}")
                    predictions[model_name] = self._get_default_parameters().get(model_name, None)

        return predictions

    def _predict_parameters_rule_based(self, features_df: pd.DataFrame) -> Dict[str, Any]:
        """Use rule-based optimization when ML models are not available."""
        if features_df.empty:
            return self._get_default_parameters()

        # Get features from the first row (assuming single sample for now)
        features = features_df.iloc[0].to_dict()

        predictions = {}

        # Rule-based Salmon threads prediction
        avg_read_length = features.get('avg_read_length', 100)
        avg_quality = features.get('avg_quality_score', 35)

        salmon_threads = self._recommend_salmon_threads(avg_read_length, avg_quality)
        predictions['salmon_threads'] = salmon_threads

        # Rule-based Salmon library type prediction
        gc_content = features.get('gc_content', 50)
        salmon_libtype = self._recommend_salmon_libtype(avg_read_length, gc_content)
        predictions['salmon_libtype'] = salmon_libtype

        # Rule-based DESeq2 alpha prediction
        deseq2_alpha = self._recommend_deseq2_alpha(avg_quality, avg_read_length)
        predictions['deseq2_alpha'] = deseq2_alpha

        logger.info("Using rule-based parameter optimization")
        return predictions

    def _get_default_parameters(self) -> Dict[str, Any]:
        """Get default parameter recommendations when models are unavailable."""
        return {
            'salmon_threads': 4,
            'salmon_libtype': 'A',
            'deseq2_alpha': 0.05
        }

    def optimize_pipeline_config(self, fastq_files: List[Path],
                               config_path: str = "config/params.yaml") -> Dict[str, Any]:
        """
        Analyze FASTQ files and generate optimized pipeline configuration.
        """
        logger.info(f"Analyzing {len(fastq_files)} FASTQ files for parameter optimization...")

        # Extract features
        features_df = self.extract_features(fastq_files)

        if features_df.empty:
            logger.warning("No valid features extracted, using defaults")
            return self._get_default_parameters()

        # Get predictions
        predictions = self.predict_optimal_parameters(features_df)

        # Generate optimized config recommendations
        optimized_config = {
            'salmon': {
                'threads': int(predictions.get('salmon_threads', 4)),
                'libtype': predictions.get('salmon_libtype', 'A')
            },
            'r': {
                'alpha': float(predictions.get('deseq2_alpha', 0.05))
            },
            'optimization_metadata': {
                'data_characteristics': features_df.iloc[0].to_dict(),
                'optimization_timestamp': pd.Timestamp.now().isoformat(),
                'confidence_scores': self._calculate_confidence_scores(features_df)
            }
        }

        logger.info("Parameter optimization complete")
        logger.info(f"Recommended Salmon threads: {optimized_config['salmon']['threads']}")
        logger.info(f"Recommended Salmon libtype: {optimized_config['salmon']['libtype']}")
        logger.info(f"Recommended DESeq2 alpha: {optimized_config['r']['alpha']}")

        return optimized_config

    def _calculate_confidence_scores(self, features_df: pd.DataFrame) -> Dict[str, float]:
        """Calculate confidence scores for predictions."""
        if features_df.empty:
            return {'overall': 0.5, 'quality_based': 0.5, 'consistency_based': 0.5}

        # Simple confidence based on data quality metrics
        avg_quality = features_df['avg_quality_score'].iloc[0]
        read_length_std = features_df['read_length_std'].iloc[0]
        # Use file size as a proxy for quantity if total_reads_analyzed not available
        file_size = features_df.get('file_size', pd.Series([0])).iloc[0]
        total_reads = features_df.get('total_reads_analyzed', pd.Series([file_size / 1000])).iloc[0]  # Rough estimate

        # Base confidence on data quality and quantity
        quality_confidence = min(1.0, avg_quality / 40.0)
        consistency_confidence = max(0.0, 1.0 - min(0.5, read_length_std / 100.0))
        quantity_confidence = min(1.0, total_reads / 1000.0)  # More reads = higher confidence

        confidence = {
            'overall': (quality_confidence + consistency_confidence + quantity_confidence) / 3,
            'quality_based': quality_confidence,
            'consistency_based': consistency_confidence,
            'quantity_based': quantity_confidence
        }

        return confidence


def create_optimization_report(optimization_results: Dict[str, Any],
                             output_path: str = "optimization_report.json") -> None:
    """
    Create a detailed report of the optimization process and recommendations.
    """
    report = {
        'summary': {
            'optimization_completed': True,
            'parameters_optimized': list(optimization_results.keys()),
            'recommendation_confidence': optimization_results.get('optimization_metadata', {}).get('confidence_scores', {})
        },
        'recommended_parameters': optimization_results,
        'data_characteristics': optimization_results.get('optimization_metadata', {}).get('data_characteristics', {}),
        'timestamp': optimization_results.get('optimization_metadata', {}).get('optimization_timestamp', '')
    }

    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)

    logger.info(f"Optimization report saved to {output_path}")


def main():
    """Command-line interface for parameter optimization."""
    import argparse

    parser = argparse.ArgumentParser(description="RNASEQ-MINI Parameter Optimizer")
    parser.add_argument('fastq_files', nargs='+', help='FASTQ files to analyze')
    parser.add_argument('--config', default='config/params.yaml', help='Configuration file to optimize')
    parser.add_argument('--output', default='optimization_report.json', help='Output report file')
    parser.add_argument('--train-models', action='store_true', help='Train new ML models')
    parser.add_argument('--model-dir', default='models', help='Directory for ML models')
    parser.add_argument('--max-reads', type=int, default=10000, help='Maximum reads to analyze per file')

    args = parser.parse_args()

    # Convert file paths
    fastq_paths = [Path(f) for f in args.fastq_files]

    # Validate FASTQ files
    valid_fastq = []
    for fastq_path in fastq_paths:
        if fastq_path.exists():
            valid_fastq.append(fastq_path)
        else:
            logger.error(f"FASTQ file not found: {fastq_path}")

    if not valid_fastq:
        logger.error("No valid FASTQ files provided")
        return 1

    # Initialize optimizer
    optimizer = ParameterOptimizer(args.model_dir)

    # Train models if requested
    if args.train_models:
        logger.info("Training new models...")
        training_stats = optimizer.train_models()
        logger.info(f"Training complete: {training_stats}")

    # Perform optimization
    optimized_config = optimizer.optimize_pipeline_config(valid_fastq, args.config)

    # Create report
    create_optimization_report(optimized_config, args.output)

    logger.info("Parameter optimization complete!")
    return 0


if __name__ == "__main__":
    exit(main())
