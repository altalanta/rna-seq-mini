#!/usr/bin/env python3
"""
Dynamic resource allocation integration for Snakemake.
Integrates resource estimation with Snakemake rule execution.
"""

import os
import sys
import json
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import pandas as pd
import numpy as np

# Import the resource estimator
from resource_estimator import ResourceEstimator, create_optimized_config

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class DynamicResourceManager:
    """Manages dynamic resource allocation for Snakemake workflows."""

    def __init__(self, params_file: str, samples_file: str):
        self.params_file = Path(params_file)
        self.samples_file = Path(samples_file)
        self.resource_cache_file = Path("results/.resource_cache.json")

        # Load current configuration
        self.config = self._load_config()

        # Check if dynamic optimization is enabled
        self.optimization_enabled = self.config.get("optimization", {}).get("enabled", False)

    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from params.yaml."""
        try:
            import yaml
            with open(self.params_file, 'r') as f:
                return yaml.safe_load(f)
        except Exception as e:
            logger.error(f"Error loading config: {e}")
            return {}

    def should_estimate_resources(self) -> bool:
        """Check if resource estimation should be performed."""
        if not self.optimization_enabled:
            return False

        # Check if samples file exists and is newer than cache
        if not self.samples_file.exists():
            return False

        # Check cache validity
        if self.resource_cache_file.exists():
            cache_mtime = self.resource_cache_file.stat().st_mtime
            samples_mtime = self.samples_file.stat().st_mtime

            # If cache is newer than samples file, use cache
            if cache_mtime > samples_mtime:
                logger.info("Using cached resource estimates")
                return False

        return True

    def estimate_and_cache_resources(self) -> Dict[str, Any]:
        """Estimate resources and cache the results."""
        logger.info("Estimating resources for dynamic allocation...")

        try:
            # Load samples
            samples_df = pd.read_csv(self.samples_file, sep='\t')

            # Get optimization config
            opt_config = self.config.get("optimization", {}).get("estimation", {})

            # Create estimator with custom config
            estimator = ResourceEstimator()
            if opt_config:
                # Update estimator config with custom parameters
                estimator.config.base_cores = opt_config.get("base_cores", 4)
                estimator.config.base_memory_gb = opt_config.get("base_memory_gb", 16)
                estimator.config.complexity_multiplier = opt_config.get("complexity_multiplier", 1.5)
                estimator.config.size_multiplier = opt_config.get("size_multiplier", 1.2)
                estimator.config.min_cores = opt_config.get("min_cores", 1)
                estimator.config.max_cores = opt_config.get("max_cores", 32)
                estimator.config.min_memory_gb = opt_config.get("min_memory_gb", 4)
                estimator.config.max_memory_gb = opt_config.get("max_memory_gb", 128)

            # Estimate resources
            estimates = estimator.estimate_project_resources(samples_df)

            # Create optimized configuration
            optimized_config = create_optimized_config(estimates)

            # Cache the results
            cache_data = {
                'timestamp': pd.Timestamp.now().isoformat(),
                'samples_hash': str(samples_df.to_dict()),
                'estimates': estimates,
                'optimized_config': optimized_config
            }

            with open(self.resource_cache_file, 'w') as f:
                json.dump(cache_data, f, indent=2)

            logger.info(f"Resource estimates cached to {self.resource_cache_file}")
            return optimized_config

        except Exception as e:
            logger.error(f"Error estimating resources: {e}")
            return {}

    def get_cached_resources(self) -> Dict[str, Any]:
        """Get cached resource estimates."""
        try:
            if not self.resource_cache_file.exists():
                return {}

            with open(self.resource_cache_file, 'r') as f:
                cache_data = json.load(f)

            logger.info("Loaded cached resource estimates")
            return cache_data.get('optimized_config', {})

        except Exception as e:
            logger.error(f"Error loading cached resources: {e}")
            return {}

    def update_snakemake_config(self, optimized_config: Dict[str, Any]):
        """Update Snakemake configuration with optimized resources."""
        if not optimized_config:
            return

        try:
            # Update main config with optimized settings
            snakemake_config = optimized_config.get('snakemake', {})

            # Update cores and resources in the main config
            if 'cores' in snakemake_config:
                self.config['threads'] = snakemake_config['cores']

            if 'resources' in snakemake_config:
                resources = snakemake_config['resources']
                if 'mem_gb' in resources:
                    self.config['memory_gb'] = resources['mem_gb']

            # Save updated config
            import yaml
            with open(self.params_file, 'w') as f:
                yaml.safe_dump(self.config, f, default_flow_style=False)

            logger.info("Updated Snakemake configuration with optimized resources")

        except Exception as e:
            logger.error(f"Error updating Snakemake config: {e}")

    def get_rule_resources(self, rule_name: str) -> Dict[str, Any]:
        """Get dynamic resources for a specific Snakemake rule."""
        if not self.optimization_enabled:
            return {}

        # Load cached estimates
        cached_data = self.get_cached_resources()
        if not cached_data:
            return {}

        # Map rule names to resource requirements
        rule_resource_map = {
            'salmon_quant': {
                'cores': min(8, self.config.get('threads', 4)),
                'memory': f"{min(32, self.config.get('memory_gb', 16))}GB",
                'time': '4h'
            },
            'fastqc': {
                'cores': 2,
                'memory': '2GB',
                'time': '1h'
            },
            'multiqc': {
                'cores': 1,
                'memory': '4GB',
                'time': '30m'
            },
            'tximport_deseq2': {
                'cores': min(4, self.config.get('threads', 4)),
                'memory': f"{min(16, self.config.get('memory_gb', 16) // 2)}GB",
                'time': '2h'
            },
            'singlecell_cellranger_count': {
                'cores': min(16, self.config.get('threads', 8)),
                'memory': f"{min(64, self.config.get('memory_gb', 32))}GB",
                'time': '8h'
            },
            'singlecell_kallisto_bus': {
                'cores': min(8, self.config.get('threads', 4)),
                'memory': f"{min(32, self.config.get('memory_gb', 16))}GB",
                'time': '2h'
            }
        }

        return rule_resource_map.get(rule_name, {
            'cores': min(4, self.config.get('threads', 4)),
            'memory': f"{min(8, self.config.get('memory_gb', 16) // 4)}GB",
            'time': '2h'
        })

def main():
    """Main function for Snakemake integration."""
    parser = argparse.ArgumentParser(description='Dynamic resource allocation for Snakemake')
    parser.add_argument('--params', required=True, help='Path to params.yaml file')
    parser.add_argument('--samples', required=True, help='Path to samples.tsv file')
    parser.add_argument('--rule', help='Snakemake rule name for resource query')
    parser.add_argument('--estimate-only', action='store_true', help='Only estimate resources, don\'t update config')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose logging')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Initialize resource manager
    manager = DynamicResourceManager(args.params, args.samples)

    # Estimate resources if needed
    if manager.should_estimate_resources():
        logger.info("Performing resource estimation...")
        optimized_config = manager.estimate_and_cache_resources()

        if not args.estimate_only:
            manager.update_snakemake_config(optimized_config)
    else:
        logger.info("Using existing resource estimates")
        optimized_config = manager.get_cached_resources()

    # Handle rule-specific resource query
    if args.rule:
        resources = manager.get_rule_resources(args.rule)
        print(json.dumps(resources))
        return

    # Print summary
    if optimized_config:
        snakemake_config = optimized_config.get('snakemake', {})
        print("
📊 Dynamic Resource Allocation Summary:"        print(f"   Optimized cores: {snakemake_config.get('cores', 'N/A')}")
        print(f"   Optimized memory: {snakemake_config.get('resources', {}).get('mem_gb', 'N/A')} GB")
        print(f"   Estimated runtime: {snakemake_config.get('resources', {}).get('time_min', 'N/A')} minutes")

if __name__ == "__main__":
    main()




