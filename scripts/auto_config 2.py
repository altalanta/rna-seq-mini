#!/usr/bin/env python3
"""
Automated configuration generator for RNASEQ-MINI pipeline.
Integrates parameter optimization with configuration management.
"""

import yaml
import json
import sys
from pathlib import Path
from typing import Dict, Any, List
import argparse
import logging

# Import our parameter optimizer
try:
    from parameter_optimizer import ParameterOptimizer
except ImportError:
    print("Error: parameter_optimizer.py not found. Please ensure it's in the scripts directory.")
    sys.exit(1)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class AutoConfigGenerator:
    """Generates optimized pipeline configurations based on data analysis."""

    def __init__(self, model_dir: str = "models"):
        self.model_dir = Path(model_dir)
        self.optimizer = ParameterOptimizer(model_dir)

    def load_base_config(self, config_path: str = "config/params.yaml") -> Dict[str, Any]:
        """Load the base configuration file."""
        config_path = Path(config_path)

        if not config_path.exists():
            logger.error(f"Configuration file not found: {config_path}")
            return {}

        try:
            with open(config_path, 'r') as f:
                return yaml.safe_load(f)
        except Exception as e:
            logger.error(f"Error loading configuration: {e}")
            return {}

    def extract_fastq_files_from_samples(self, samples_path: str = "config/samples.tsv") -> List[Path]:
        """Extract FASTQ file paths from samples configuration."""
        samples_path = Path(samples_path)

        if not samples_path.exists():
            logger.error(f"Samples file not found: {samples_path}")
            return []

        try:
            import pandas as pd
            samples_df = pd.read_csv(samples_path, sep='\t')

            fastq_files = []
            for _, row in samples_df.iterrows():
                fastq1 = row.get('fastq_1')
                fastq2 = row.get('fastq_2')

                if fastq1:
                    fastq_files.append(Path(fastq1))
                if fastq2:
                    fastq_files.append(Path(fastq2))

            return list(set(fastq_files))  # Remove duplicates

        except Exception as e:
            logger.error(f"Error reading samples file: {e}")
            return []

    def generate_optimized_config(self, fastq_files: List[Path],
                                base_config_path: str = "config/params.yaml",
                                output_path: str = "config/params_optimized.yaml") -> bool:
        """
        Generate an optimized configuration based on FASTQ file analysis.
        """
        logger.info("Starting automated configuration optimization...")

        # Load base configuration
        base_config = self.load_base_config(base_config_path)
        if not base_config:
            logger.error("Failed to load base configuration")
            return False

        # Extract FASTQ files if not provided
        if not fastq_files:
            fastq_files = self.extract_fastq_files_from_samples()
            if not fastq_files:
                logger.error("No FASTQ files found in samples configuration")
                return False

        # Validate FASTQ files exist
        valid_fastq = []
        for fastq_file in fastq_files:
            if fastq_file.exists():
                valid_fastq.append(fastq_file)
            else:
                logger.warning(f"FASTQ file not found: {fastq_file}")

        if not valid_fastq:
            logger.error("No valid FASTQ files found")
            return False

        # Perform parameter optimization
        logger.info(f"Analyzing {len(valid_fastq)} FASTQ files...")
        optimized_params = self.optimizer.optimize_pipeline_config(valid_fastq, base_config_path)

        if not optimized_params:
            logger.error("Parameter optimization failed")
            return False

        # Merge optimized parameters with base configuration
        optimized_config = self._merge_configurations(base_config, optimized_params)

        # Add metadata about optimization
        optimized_config['optimization_metadata'] = {
            'auto_generated': True,
            'source_config': base_config_path,
            'optimization_timestamp': optimized_params.get('optimization_metadata', {}).get('optimization_timestamp', ''),
            'data_characteristics': optimized_params.get('optimization_metadata', {}).get('data_characteristics', {}),
            'confidence_scores': optimized_params.get('optimization_metadata', {}).get('confidence_scores', {})
        }

        # Save optimized configuration
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            with open(output_path, 'w') as f:
                yaml.dump(optimized_config, f, default_flow_style=False, sort_keys=False)

            logger.info(f"Optimized configuration saved to: {output_path}")

            # Display summary of changes
            self._display_optimization_summary(base_config, optimized_config)

            return True

        except Exception as e:
            logger.error(f"Error saving optimized configuration: {e}")
            return False

    def _merge_configurations(self, base_config: Dict[str, Any],
                            optimized_params: Dict[str, Any]) -> Dict[str, Any]:
        """Merge optimized parameters with base configuration."""
        merged_config = base_config.copy()

        # Update Salmon parameters
        if 'salmon' in optimized_params:
            if 'salmon' not in merged_config:
                merged_config['salmon'] = {}

            merged_config['salmon'].update(optimized_params['salmon'])

        # Update R/DESeq2 parameters
        if 'r' in optimized_params:
            if 'r' not in merged_config:
                merged_config['r'] = {}

            merged_config['r'].update(optimized_params['r'])

        return merged_config

    def _display_optimization_summary(self, base_config: Dict[str, Any],
                                    optimized_config: Dict[str, Any]) -> None:
        """Display a summary of configuration changes."""
        logger.info("\n" + "="*50)
        logger.info("PARAMETER OPTIMIZATION SUMMARY")
        logger.info("="*50)

        changes = []

        # Compare Salmon parameters
        base_salmon = base_config.get('salmon', {})
        opt_salmon = optimized_config.get('salmon', {})

        if base_salmon.get('threads') != opt_salmon.get('threads'):
            changes.append(f"Salmon threads: {base_salmon.get('threads', 4)} → {opt_salmon.get('threads', 4)}")

        if base_salmon.get('libtype') != opt_salmon.get('libtype'):
            changes.append(f"Salmon libtype: {base_salmon.get('libtype', 'A')} → {opt_salmon.get('libtype', 'A')}")

        # Compare DESeq2 parameters
        base_r = base_config.get('r', {})
        opt_r = optimized_config.get('r', {})

        if base_r.get('alpha') != opt_r.get('alpha'):
            changes.append(f"DESeq2 alpha: {base_r.get('alpha', 0.05)} → {opt_r.get('alpha', 0.05)}")

        if changes:
            logger.info("Configuration changes:")
            for change in changes:
                logger.info(f"  • {change}")
        else:
            logger.info("No parameter changes recommended (using existing values)")

        # Show confidence scores
        confidence = optimized_config.get('optimization_metadata', {}).get('confidence_scores', {})
        if confidence:
            logger.info("\nConfidence scores:")
            for metric, score in confidence.items():
                logger.info(f"  • {metric}: {score:.2%}")

        logger.info("="*50)

    def validate_optimized_config(self, config_path: str) -> bool:
        """Validate the generated configuration file."""
        try:
            from validate_config import run_comprehensive_validation
            success, issues = run_comprehensive_validation(config_path)

            if success:
                logger.info("✅ Optimized configuration is valid")
                return True
            else:
                logger.error("❌ Optimized configuration has validation errors:")
                for issue in issues:
                    logger.error(f"  • {issue}")
                return False

        except ImportError:
            logger.warning("Configuration validation script not found, skipping validation")
            return True


def main():
    """Command-line interface for automated configuration generation."""
    parser = argparse.ArgumentParser(description="RNASEQ-MINI Automated Configuration Generator")
    parser.add_argument('fastq_files', nargs='*', help='FASTQ files to analyze (if not provided, reads from samples.tsv)')
    parser.add_argument('--base-config', default='config/params.yaml', help='Base configuration file')
    parser.add_argument('--output', default='config/params_optimized.yaml', help='Output configuration file')
    parser.add_argument('--samples-file', default='config/samples.tsv', help='Samples configuration file')
    parser.add_argument('--model-dir', default='models', help='Directory containing ML models')
    parser.add_argument('--validate', action='store_true', help='Validate the generated configuration')

    args = parser.parse_args()

    # Convert file paths
    fastq_paths = [Path(f) for f in args.fastq_files] if args.fastq_files else []

    # Initialize configuration generator
    generator = AutoConfigGenerator(args.model_dir)

    # Generate optimized configuration
    success = generator.generate_optimized_config(
        fastq_paths,
        args.base_config,
        args.output
    )

    if not success:
        logger.error("Configuration generation failed")
        return 1

    # Validate if requested
    if args.validate:
        if not generator.validate_optimized_config(args.output):
            logger.error("Generated configuration failed validation")
            return 1

    logger.info("✅ Automated configuration generation complete!")
    logger.info(f"Optimized configuration: {args.output}")

    # Show usage instructions
    logger.info("\nTo use the optimized configuration:")
    logger.info(f"  cp {args.output} config/params.yaml")
    logger.info("  make run")

    return 0


if __name__ == "__main__":
    exit(main())
