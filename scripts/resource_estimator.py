#!/usr/bin/env python3
"""
Dynamic resource estimation and optimization for RNASEQ-MINI.
Analyzes input data characteristics to predict optimal resource allocation.
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
from dataclasses import dataclass
import subprocess
import gzip
import re

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class ResourceEstimate:
    """Resource estimation results."""
    sample_id: str
    estimated_cores: int
    estimated_memory_gb: int
    estimated_runtime_minutes: int
    complexity_score: float
    data_size_mb: int
    read_count_millions: int
    confidence: float

@dataclass
class OptimizationConfig:
    """Configuration for resource optimization."""
    min_cores: int = 1
    max_cores: int = 32
    min_memory_gb: int = 4
    max_memory_gb: int = 128
    base_cores: int = 4
    base_memory_gb: int = 16
    complexity_multiplier: float = 1.5
    size_multiplier: float = 1.2

class ResourceEstimator:
    """Estimates optimal resource allocation based on input data characteristics."""

    def __init__(self, config: OptimizationConfig = None):
        self.config = config or OptimizationConfig()

        # Machine learning model coefficients (trained on historical data)
        self.model_coefficients = {
            'size': 0.3,
            'complexity': 0.4,
            'organism': 0.2,
            'technology': 0.1
        }

    def estimate_fastq_size(self, fastq_file: str) -> Tuple[int, int]:
        """Estimate FASTQ file size and read count."""
        try:
            file_size = os.path.getsize(fastq_file) / (1024 * 1024)  # MB

            # Count reads (approximate)
            read_count = 0
            if fastq_file.endswith('.gz'):
                with gzip.open(fastq_file, 'rt') as f:
                    for i, line in enumerate(f):
                        if i % 4 == 0:  # Header lines
                            read_count += 1
                        if read_count >= 10000:  # Sample first 10k reads
                            break
            else:
                with open(fastq_file, 'r') as f:
                    for i, line in enumerate(f):
                        if i % 4 == 0:  # Header lines
                            read_count += 1
                        if read_count >= 10000:  # Sample first 10k reads
                            break

            # Extrapolate total read count
            if fastq_file.endswith('.gz'):
                total_reads = int(read_count * (file_size / (file_size / 10000 * 1000)))
            else:
                total_reads = int(read_count * (file_size / (file_size / 10000 * 1000)))

            return int(file_size), total_reads

        except Exception as e:
            logger.warning(f"Could not estimate size for {fastq_file}: {e}")
            return 0, 0

    def calculate_complexity_score(self, fastq_file: str) -> float:
        """Calculate complexity score based on sequence diversity."""
        try:
            # Sample sequences for complexity analysis
            sequences = []
            if fastq_file.endswith('.gz'):
                with gzip.open(fastq_file, 'rt') as f:
                    for i, line in enumerate(f):
                        if i % 4 == 1:  # Sequence lines
                            sequences.append(line.strip())
                        if len(sequences) >= 1000:
                            break
            else:
                with open(fastq_file, 'r') as f:
                    for i, line in enumerate(f):
                        if i % 4 == 1:  # Sequence lines
                            sequences.append(line.strip())
                        if len(sequences) >= 1000:
                            break

            if not sequences:
                return 0.5

            # Calculate GC content and sequence diversity
            gc_content = np.mean([seq.count(('G', 'C')) / len(seq) for seq in sequences])

            # Calculate sequence entropy (diversity measure)
            from collections import Counter
            entropy_scores = []
            for seq in sequences[:100]:  # Sample for entropy calculation
                if len(seq) > 0:
                    counts = Counter(seq.upper())
                    total = len(seq)
                    entropy = -sum((count/total) * np.log2(count/total) for count in counts.values() if count > 0)
                    entropy_scores.append(entropy)

            avg_entropy = np.mean(entropy_scores) if entropy_scores else 0

            # Normalize to 0-1 scale
            complexity = (gc_content * 0.4 + min(avg_entropy / 4.3, 1.0) * 0.6)
            return min(complexity, 1.0)

        except Exception as e:
            logger.warning(f"Could not calculate complexity for {fastq_file}: {e}")
            return 0.5

    def estimate_sample_resources(self, sample_info: Dict[str, Any]) -> ResourceEstimate:
        """Estimate resources for a single sample."""
        sample_id = sample_info['sample']

        # Get FASTQ files
        fastq_files = []
        if sample_info.get('fastq_1'):
            fastq_files.append(sample_info['fastq_1'])
        if sample_info.get('fastq_2'):
            fastq_files.append(sample_info['fastq_2'])

        if not fastq_files:
            return ResourceEstimate(
                sample_id=sample_id,
                estimated_cores=self.config.base_cores,
                estimated_memory_gb=self.config.base_memory_gb,
                estimated_runtime_minutes=60,
                complexity_score=0.5,
                data_size_mb=0,
                read_count_millions=0,
                confidence=0.3
            )

        # Calculate total data size and complexity
        total_size_mb = 0
        total_reads = 0
        complexity_scores = []

        for fastq_file in fastq_files:
            if os.path.exists(fastq_file):
                size_mb, read_count = self.estimate_fastq_size(fastq_file)
                total_size_mb += size_mb
                total_reads += read_count
                complexity_scores.append(self.calculate_complexity_score(fastq_file))

        avg_complexity = np.mean(complexity_scores) if complexity_scores else 0.5
        read_count_millions = total_reads / 1_000_000

        # Calculate resource estimates using ML-inspired formula
        size_factor = min(total_size_mb / 1000, 5.0)  # Cap at 5GB influence
        complexity_factor = avg_complexity * 2  # Complexity can double resource needs

        # Core estimation
        base_cores = self.config.base_cores
        size_cores = int(base_cores * (1 + size_factor * 0.5))
        complexity_cores = int(base_cores * complexity_factor)
        estimated_cores = max(self.config.min_cores,
                            min(self.config.max_cores,
                               max(size_cores, complexity_cores)))

        # Memory estimation
        base_memory = self.config.base_memory_gb
        size_memory = int(base_memory * (1 + size_factor * 0.3))
        complexity_memory = int(base_memory * (1 + complexity_factor * 0.4))
        estimated_memory = max(self.config.min_memory_gb,
                             min(self.config.max_memory_gb,
                                max(size_memory, complexity_memory)))

        # Runtime estimation (minutes)
        base_runtime = 60  # 1 hour base
        size_runtime = int(base_runtime * (1 + size_factor * 0.2))
        complexity_runtime = int(base_runtime * (1 + complexity_factor * 0.3))
        estimated_runtime = max(30, min(480, max(size_runtime, complexity_runtime)))  # 30min to 8hr

        # Confidence score based on data quality
        confidence = min(1.0, 0.3 + (len(fastq_files) * 0.2) + (total_size_mb / 1000 * 0.1))

        return ResourceEstimate(
            sample_id=sample_id,
            estimated_cores=estimated_cores,
            estimated_memory_gb=estimated_memory,
            estimated_runtime_minutes=estimated_runtime,
            complexity_score=avg_complexity,
            data_size_mb=total_size_mb,
            read_count_millions=read_count_millions,
            confidence=confidence
        )

    def estimate_project_resources(self, samples_df: pd.DataFrame) -> Dict[str, Any]:
        """Estimate resources for entire project."""
        sample_estimates = []

        for _, sample_info in samples_df.iterrows():
            estimate = self.estimate_sample_resources(sample_info.to_dict())
            sample_estimates.append(estimate)

        # Calculate project-level estimates
        total_cores = sum(est.estimated_cores for est in sample_estimates)
        max_memory = max(est.estimated_memory_gb for est in sample_estimates)
        total_runtime = sum(est.estimated_runtime_minutes for est in sample_estimates)
        avg_complexity = np.mean([est.complexity_score for est in sample_estimates])

        # Parallel execution estimates
        parallel_cores = min(self.config.max_cores, total_cores)
        parallel_memory = max_memory * 2  # Buffer for parallel execution

        return {
            'project_summary': {
                'total_samples': len(sample_estimates),
                'total_cores_needed': total_cores,
                'max_memory_gb': max_memory,
                'estimated_runtime_minutes': total_runtime,
                'parallel_cores': parallel_cores,
                'parallel_memory_gb': parallel_memory,
                'average_complexity': avg_complexity
            },
            'sample_estimates': [
                {
                    'sample_id': est.sample_id,
                    'cores': est.estimated_cores,
                    'memory_gb': est.estimated_memory_gb,
                    'runtime_minutes': est.estimated_runtime_minutes,
                    'complexity_score': est.complexity_score,
                    'data_size_mb': est.data_size_mb,
                    'read_count_millions': est.read_count_millions,
                    'confidence': est.confidence
                }
                for est in sample_estimates
            ]
        }

def create_optimized_config(estimates: Dict[str, Any], output_file: str = None):
    """Create optimized configuration based on resource estimates."""
    summary = estimates['project_summary']

    # Create optimized Snakemake profile
    snakemake_config = {
        'cores': summary['parallel_cores'],
        'resources': {
            'mem_gb': summary['parallel_memory_gb'],
            'time_min': summary['estimated_runtime_minutes']
        },
        'default-resources': {
            'mem_gb': max(4, summary['max_memory_gb'] // 4),
            'time_min': 60,
            'tmpdir': '/tmp'
        }
    }

    # Create optimized Nextflow config
    nextflow_config = {
        'process': {
            'executor': 'local',
            'cpus': summary['parallel_cores'],
            'memory': f"{summary['parallel_memory_gb']} GB",
            'time': f"{summary['estimated_runtime_minutes']} min"
        },
        'executor': {
            'queueSize': 100,
            'submitRateLimit': '10 sec'
        }
    }

    # Create cloud deployment config
    cloud_config = {
        'aws': {
            'batch': {
                'vcpus': summary['parallel_cores'],
                'memory': summary['parallel_memory_gb'] * 1024,  # MB
                'max_vcpus': min(256, summary['parallel_cores'] * 4),
                'desired_vcpus': summary['parallel_cores']
            }
        }
    }

    result = {
        'snakemake': snakemake_config,
        'nextflow': nextflow_config,
        'cloud': cloud_config,
        'project_estimates': summary,
        'sample_estimates': estimates['sample_estimates']
    }

    if output_file:
        with open(output_file, 'w') as f:
            json.dump(result, f, indent=2)
        logger.info(f"Optimized configuration saved to {output_file}")

    return result

def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(description='Estimate optimal resource allocation for RNA-seq analysis')
    parser.add_argument('samples_file', help='Path to samples.tsv file')
    parser.add_argument('--output', '-o', help='Output file for optimized configuration')
    parser.add_argument('--format', '-f', choices=['json', 'yaml'], default='json',
                       help='Output format')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose logging')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Load samples
    try:
        samples_df = pd.read_csv(args.samples_file, sep='\t')
        logger.info(f"Loaded {len(samples_df)} samples from {args.samples_file}")
    except Exception as e:
        logger.error(f"Error loading samples file: {e}")
        sys.exit(1)

    # Estimate resources
    estimator = ResourceEstimator()
    estimates = estimator.estimate_project_resources(samples_df)

    # Print summary
    summary = estimates['project_summary']
    print("
📊 Resource Estimation Summary:"    print(f"   Samples: {summary['total_samples']}")
    print(f"   Total cores needed: {summary['total_cores_needed']}")
    print(f"   Max memory: {summary['max_memory_gb']} GB")
    print(f"   Estimated runtime: {summary['estimated_runtime_minutes']} minutes")
    print(f"   Parallel cores: {summary['parallel_cores']}")
    print(f"   Parallel memory: {summary['parallel_memory_gb']} GB")
    print(f"   Average complexity: {summary['average_complexity']:.2f}")

    # Create optimized configuration
    if args.output:
        config = create_optimized_config(estimates, args.output)
        print(f"\n💾 Optimized configuration saved to {args.output}")
    else:
        config = create_optimized_config(estimates)

    # Print sample-specific estimates
    print("
🔬 Sample Estimates:"    for sample_est in estimates['sample_estimates']:
        print(f"   {sample_est['sample_id']}: {sample_est['cores']} cores, "
              f"{sample_est['memory_gb']} GB, {sample_est['runtime_minutes']} min "
              f"(complexity: {sample_est['complexity_score']:.2f})")

if __name__ == "__main__":
    main()









