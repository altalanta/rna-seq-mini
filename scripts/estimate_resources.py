#!/usr/bin/env python3
"""
Resource estimation and auto-scaling script for RNASEQ-MINI.
Estimates optimal resource allocation based on dataset size and available hardware.

Usage: python scripts/estimate_resources.py [config_file] [samples_file]
"""

import argparse
import sys
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
import psutil
import os
import json


class ResourceEstimator:
    """Estimates optimal resource allocation for RNA-seq pipeline."""

    def __init__(self):
        # Resource scaling factors based on empirical data
        self.scaling_factors = {
            'fastqc': {
                'cpu_per_sample': 0.5,  # CPU cores per sample
                'memory_per_sample_gb': 0.5,  # GB per sample
                'time_per_sample_min': 5  # Minutes per sample
            },
            'salmon': {
                'cpu_per_sample': 2,  # CPU cores per sample
                'memory_per_sample_gb': 4,  # GB per sample
                'time_per_sample_min': 30  # Minutes per sample
            },
            'deseq2': {
                'cpu_cores': 4,  # Fixed CPU requirement
                'memory_gb': 8,  # Fixed memory requirement
                'time_min': 15  # Fixed time estimate
            },
            'fgsea': {
                'cpu_cores': 2,  # Fixed CPU requirement
                'memory_gb': 4,  # Fixed memory requirement
                'time_min': 10  # Fixed time estimate
            }
        }

        # Overhead factors
        self.overhead_factors = {
            'memory_buffer': 1.2,  # 20% buffer
            'cpu_efficiency': 0.8,  # 80% efficiency
            'parallelization_limit': 8  # Max parallel samples
        }

    def get_system_resources(self):
        """Get available system resources."""
        try:
            cpu_count = os.cpu_count() or 4
            memory_gb = psutil.virtual_memory().total / (1024**3)

            # Estimate disk space (rough approximation)
            disk_usage = psutil.disk_usage('/')
            disk_gb = disk_usage.free / (1024**3)

            return {
                'cpu_cores': cpu_count,
                'memory_gb': memory_gb,
                'disk_gb': disk_gb,
                'cpu_efficiency': 0.85  # Conservative estimate
            }
        except Exception as e:
            print(f"Warning: Could not detect system resources: {e}")
            return {
                'cpu_cores': 8,
                'memory_gb': 32,
                'disk_gb': 100,
                'cpu_efficiency': 0.8
            }

    def estimate_fastq_sizes(self, samples_file):
        """Estimate FASTQ file sizes from samples table."""
        try:
            samples_df = pd.read_csv(samples_file, sep='\t')

            total_samples = len(samples_df)
            estimated_sizes = []

            # Check if FASTQ files exist and get their sizes
            for _, row in samples_df.iterrows():
                fastq1 = Path(row['fastq_1'])
                size_mb = 0

                if fastq1.exists():
                    size_mb += fastq1.stat().st_size / (1024 * 1024)

                    # Add R2 if paired-end
                    if 'fastq_2' in row and row['fastq_2'] and Path(row['fastq_2']).exists():
                        fastq2 = Path(row['fastq_2'])
                        size_mb += fastq2.stat().st_size / (1024 * 1024)

                # If files don't exist, estimate based on typical RNA-seq file sizes
                if size_mb == 0:
                    # Typical RNA-seq sample: ~50M reads, ~100bp = ~5GB per sample for paired-end
                    size_mb = 5000 if 'fastq_2' in row and row['fastq_2'] else 2500

                estimated_sizes.append(size_mb)

            return {
                'total_samples': total_samples,
                'total_size_gb': sum(estimated_sizes) / 1024,
                'avg_size_gb': np.mean(estimated_sizes) / 1024,
                'max_size_gb': max(estimated_sizes) / 1024
            }

        except Exception as e:
            print(f"Warning: Could not estimate FASTQ sizes: {e}")
            # Return conservative estimates
            return {
                'total_samples': 6,  # Default from test data
                'total_size_gb': 3,  # Conservative estimate
                'avg_size_gb': 0.5,
                'max_size_gb': 1
            }

    def calculate_optimal_resources(self, samples_file, custom_config=None):
        """Calculate optimal resource allocation."""
        system_resources = self.get_system_resources()
        fastq_info = self.estimate_fastq_sizes(samples_file)

        print("📊 System Resources Detected:")
        print(f"   CPU cores: {system_resources['cpu_cores']}")
        print(f"   Memory: {system_resources['memory_gb']:.1f} GB")
        print(f"   Disk space: {system_resources['disk_gb']:.1f} GB")

        print("\n🧬 Dataset Information:")
        print(f"   Samples: {fastq_info['total_samples']}")
        print(f"   Total size: {fastq_info['total_size_gb']:.1f} GB")
        print(f"   Average per sample: {fastq_info['avg_size_gb']:.1f} GB")

        # Calculate stage-specific requirements
        stage_requirements = {}

        # FastQC requirements
        fastqc_cpu = min(
            fastq_info['total_samples'] * self.scaling_factors['fastqc']['cpu_per_sample'],
            system_resources['cpu_cores'] * self.overhead_factors['cpu_efficiency']
        )
        fastqc_memory = fastq_info['total_samples'] * self.scaling_factors['fastqc']['memory_per_sample_gb']
        stage_requirements['fastqc'] = {
            'cpu_cores': max(1, int(fastqc_cpu)),
            'memory_gb': max(1, int(fastqc_memory * self.overhead_factors['memory_buffer'])),
            'time_min': fastq_info['total_samples'] * self.scaling_factors['fastqc']['time_per_sample_min']
        }

        # Salmon requirements (most demanding)
        max_parallel_samples = min(
            int(system_resources['cpu_cores'] * self.overhead_factors['cpu_efficiency'] / self.scaling_factors['salmon']['cpu_per_sample']),
            self.overhead_factors['parallelization_limit'],
            fastq_info['total_samples']
        )

        salmon_cpu = max_parallel_samples * self.scaling_factors['salmon']['cpu_per_sample']
        salmon_memory = max_parallel_samples * self.scaling_factors['salmon']['memory_per_sample_gb']
        stage_requirements['salmon'] = {
            'cpu_cores': max(1, int(salmon_cpu)),
            'memory_gb': max(4, int(salmon_memory * self.overhead_factors['memory_buffer'])),
            'time_min': (fastq_info['total_samples'] / max_parallel_samples) * self.scaling_factors['salmon']['time_per_sample_min'],
            'parallel_samples': max_parallel_samples
        }

        # DESeq2 requirements (fixed)
        stage_requirements['deseq2'] = {
            'cpu_cores': min(self.scaling_factors['deseq2']['cpu_cores'], system_resources['cpu_cores']),
            'memory_gb': max(4, int(self.scaling_factors['deseq2']['memory_gb'] * self.overhead_factors['memory_buffer'])),
            'time_min': self.scaling_factors['deseq2']['time_min']
        }

        # FGSEA requirements (fixed)
        stage_requirements['fgsea'] = {
            'cpu_cores': min(self.scaling_factors['fgsea']['cpu_cores'], system_resources['cpu_cores']),
            'memory_gb': max(2, int(self.scaling_factors['fgsea']['memory_gb'] * self.overhead_factors['memory_buffer'])),
            'time_min': self.scaling_factors['fgsea']['time_min']
        }

        # Overall recommendations
        max_cpu = max(req['cpu_cores'] for req in stage_requirements.values())
        max_memory = max(req['memory_gb'] for req in stage_requirements.values())
        total_time = sum(req['time_min'] for req in stage_requirements.values())

        # Ensure we don't exceed system limits
        recommended_cpu = min(max_cpu, int(system_resources['cpu_cores'] * 0.9))
        recommended_memory = min(max_memory, int(system_resources['memory_gb'] * 0.85))

        return {
            'system': system_resources,
            'dataset': fastq_info,
            'stages': stage_requirements,
            'recommendations': {
                'threads': recommended_cpu,
                'memory_gb': recommended_memory,
                'estimated_time_min': total_time,
                'estimated_time_hours': total_time / 60,
                'disk_space_needed_gb': fastq_info['total_size_gb'] * 3,  # Conservative multiplier
                'warnings': self._generate_warnings(system_resources, fastq_info, recommended_cpu, recommended_memory)
            }
        }

    def _generate_warnings(self, system, dataset, recommended_cpu, recommended_memory):
        """Generate warnings about resource constraints."""
        warnings = []

        if recommended_cpu < system['cpu_cores'] * 0.5:
            warnings.append("Very low CPU utilization - consider increasing thread count for faster execution")

        if recommended_memory > system['memory_gb'] * 0.8:
            warnings.append(f"High memory usage ({recommended_memory}GB vs {system['memory_gb']:.1f}GB available) - monitor memory usage")

        if dataset['total_size_gb'] > system['disk_gb'] * 0.1:
            warnings.append(f"Large dataset ({dataset['total_size_gb']:.1f}GB) - ensure sufficient disk space")

        if dataset['total_samples'] > 50:
            warnings.append("Large number of samples - consider batch processing or HPC cluster")

        return warnings

    def generate_optimized_config(self, samples_file, output_config=None):
        """Generate an optimized configuration file."""
        estimates = self.calculate_optimal_resources(samples_file)

        # Load existing config if available
        config_file = Path("config/params.yaml")
        if config_file.exists():
            with open(config_file, 'r') as f:
                config = yaml.safe_load(f)
        else:
            config = {}

        # Update with optimized values
        config['threads'] = estimates['recommendations']['threads']
        config['memory_gb'] = estimates['recommendations']['memory_gb']

        # Optimize Salmon threads (typically 1/2 to 1/3 of total threads)
        config['salmon'] = config.get('salmon', {})
        config['salmon']['threads'] = max(1, estimates['recommendations']['threads'] // 3)

        # Update output file if specified
        if output_config:
            output_path = Path(output_config)
        else:
            output_path = Path("config/params_optimized.yaml")

        # Write optimized config
        with open(output_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)

        return estimates, output_path

    def print_resource_report(self, estimates):
        """Print a detailed resource report."""
        print("\n🎯 Resource Recommendations:")
        print(f"   Recommended threads: {estimates['recommendations']['threads']}")
        print(f"   Recommended memory: {estimates['recommendations']['memory_gb']} GB")
        print(f"   Estimated runtime: {estimates['recommendations']['estimated_time_hours']:.1f} hours")
        print(f"   Disk space needed: {estimates['recommendations']['disk_space_needed_gb']:.1f} GB")

        print("\n📈 Stage-by-stage breakdown:")
        for stage, req in estimates['stages'].items():
            print(f"   {stage.upper()}:")
            print(f"     CPU cores: {req['cpu_cores']}")
            print(f"     Memory: {req['memory_gb']} GB")
            print(f"     Time: {req['time_min']:.0f} min")

        if estimates['recommendations']['warnings']:
            print("\n⚠️  Warnings:")
            for warning in estimates['recommendations']['warnings']:
                print(f"   {warning}")

        print("\n💡 Tips:")
        print("   - Use 'make run' to execute with these settings")
        print("   - Monitor resource usage during execution")
        print("   - Consider HPC cluster for large datasets (>50 samples)")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Estimate optimal resource allocation for RNA-seq pipeline")
    parser.add_argument('samples_file', nargs='?', default='config/samples.tsv', help='Samples file to analyze')
    parser.add_argument('--config', '-c', help='Configuration file to optimize')
    parser.add_argument('--output', '-o', help='Output file for optimized config')
    parser.add_argument('--report-only', action='store_true', help='Only show report, don\'t generate config')

    args = parser.parse_args()

    estimator = ResourceEstimator()

    try:
        if args.report_only:
            # Just show the report
            estimates = estimator.calculate_optimal_resources(args.samples_file)
            estimator.print_resource_report(estimates)
        else:
            # Generate optimized config
            estimates, output_path = estimator.generate_optimized_config(
                args.samples_file,
                args.output
            )

            print(f"✅ Optimized configuration written to: {output_path}")
            estimator.print_resource_report(estimates)

    except Exception as e:
        print(f"❌ Error during resource estimation: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
