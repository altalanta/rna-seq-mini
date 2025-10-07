#!/usr/bin/env python3
"""
Interactive progress tracking and quality metrics dashboard for RNASEQ-MINI.
Provides real-time monitoring of pipeline execution and quality metrics.

Usage: python scripts/monitor_progress.py [--live] [--web] [log_dir]
"""

import argparse
import sys
import time
import json
import yaml
import psutil
import threading
from pathlib import Path
from datetime import datetime, timedelta
import subprocess
import os
from collections import defaultdict, deque
import re


class PipelineMonitor:
    """Monitors RNA-seq pipeline progress and quality metrics."""

    def __init__(self, log_dir="logs", update_interval=2):
        self.log_dir = Path(log_dir)
        self.update_interval = update_interval
        self.start_time = time.time()

        # Progress tracking
        self.stage_progress = {
            'fastqc': {'status': 'pending', 'progress': 0, 'samples': 0, 'total_samples': 0},
            'salmon': {'status': 'pending', 'progress': 0, 'samples': 0, 'total_samples': 0},
            'counts': {'status': 'pending', 'progress': 0, 'completed': False},
            'deseq2': {'status': 'pending', 'progress': 0, 'contrasts': 0, 'total_contrasts': 0},
            'fgsea': {'status': 'pending', 'progress': 0, 'contrasts': 0, 'total_contrasts': 0},
            'report': {'status': 'pending', 'progress': 0, 'completed': False}
        }

        # Quality metrics
        self.quality_metrics = {
            'fastqc': {'total_samples': 0, 'passed_samples': 0, 'failed_samples': 0},
            'multiqc': {'generated': False, 'report_exists': False},
            'salmon': {'total_samples': 0, 'successful_samples': 0, 'failed_samples': 0},
            'deseq2': {'total_contrasts': 0, 'successful_contrasts': 0, 'failed_contrasts': 0},
            'fgsea': {'total_contrasts': 0, 'successful_contrasts': 0, 'failed_contrasts': 0}
        }

        # Resource monitoring
        self.resource_history = {
            'cpu_percent': deque(maxlen=100),
            'memory_percent': deque(maxlen=100),
            'disk_usage': deque(maxlen=100)
        }

        # Pipeline status
        self.pipeline_status = 'unknown'
        self.current_stage = None

    def detect_pipeline_type(self):
        """Detect which pipeline engine is being used."""
        # Check for Snakemake files
        if (self.log_dir / 'snakemake.log').exists():
            return 'snakemake'

        # Check for Nextflow files
        if (self.log_dir / '.nextflow.log').exists() or (self.log_dir / 'nextflow.log').exists():
            return 'nextflow'

        return 'unknown'

    def parse_snakemake_logs(self):
        """Parse Snakemake log files for progress information."""
        log_files = [
            self.log_dir / 'snakemake.log',
            self.log_dir / 'slurm' / 'snakemake.log'
        ]

        for log_file in log_files:
            if log_file.exists():
                try:
                    with open(log_file, 'r') as f:
                        lines = f.readlines()[-100:]  # Last 100 lines

                    for line in reversed(lines):
                        # Look for stage completion messages
                        if 'Finished job' in line:
                            # Extract stage information
                            if 'fastqc' in line.lower():
                                self.update_stage_progress('fastqc', 'running')
                            elif 'salmon' in line.lower():
                                self.update_stage_progress('salmon', 'running')
                            elif 'tximport' in line.lower():
                                self.update_stage_progress('counts', 'running')
                            elif 'deseq2' in line.lower():
                                self.update_stage_progress('deseq2', 'running')
                            elif 'fgsea' in line.lower():
                                self.update_stage_progress('fgsea', 'running')
                            elif 'report' in line.lower():
                                self.update_stage_progress('report', 'running')

                        # Look for error messages
                        if 'error' in line.lower() or 'failed' in line.lower():
                            # Extract which stage failed
                            if 'fastqc' in line.lower():
                                self.stage_progress['fastqc']['status'] = 'failed'
                            elif 'salmon' in line.lower():
                                self.stage_progress['salmon']['status'] = 'failed'
                            elif 'tximport' in line.lower():
                                self.stage_progress['counts']['status'] = 'failed'
                            elif 'deseq2' in line.lower():
                                self.stage_progress['deseq2']['status'] = 'failed'
                            elif 'fgsea' in line.lower():
                                self.stage_progress['fgsea']['status'] = 'failed'
                            elif 'report' in line.lower():
                                self.stage_progress['report']['status'] = 'failed'

                except Exception as e:
                    print(f"Warning: Could not parse Snakemake logs: {e}")

    def parse_nextflow_logs(self):
        """Parse Nextflow log files for progress information."""
        log_files = [
            self.log_dir / '.nextflow.log',
            self.log_dir / 'nextflow.log'
        ]

        for log_file in log_files:
            if log_file.exists():
                try:
                    with open(log_file, 'r') as f:
                        lines = f.readlines()[-50:]  # Last 50 lines

                    for line in reversed(lines):
                        # Look for task completion messages
                        if 'task completed' in line.lower() or 'process completed' in line.lower():
                            # Extract stage information
                            if 'fastqc' in line.lower():
                                self.update_stage_progress('fastqc', 'running')
                            elif 'salmon' in line.lower():
                                self.update_stage_progress('salmon', 'running')
                            elif 'tximport' in line.lower():
                                self.update_stage_progress('counts', 'running')
                            elif 'deseq2' in line.lower():
                                self.update_stage_progress('deseq2', 'running')
                            elif 'fgsea' in line.lower():
                                self.update_stage_progress('fgsea', 'running')
                            elif 'report' in line.lower():
                                self.update_stage_progress('report', 'running')

                except Exception as e:
                    print(f"Warning: Could not parse Nextflow logs: {e}")

    def check_output_files(self):
        """Check for output files to determine stage completion."""
        results_dir = Path("results")

        # Check FastQC completion
        fastqc_dir = results_dir / "qc" / "fastqc"
        if fastqc_dir.exists():
            html_files = list(fastqc_dir.glob("*_fastqc.html"))
            self.stage_progress['fastqc']['total_samples'] = self._get_total_samples()
            self.stage_progress['fastqc']['samples'] = len(html_files)

            if self.stage_progress['fastqc']['samples'] > 0:
                progress = min(100, (self.stage_progress['fastqc']['samples'] / self.stage_progress['fastqc']['total_samples']) * 100)
                self.stage_progress['fastqc']['progress'] = progress

                if progress >= 100:
                    self.stage_progress['fastqc']['status'] = 'completed'

        # Check MultiQC
        multiqc_report = results_dir / "qc" / "multiqc" / "multiqc_report.html"
        if multiqc_report.exists():
            self.quality_metrics['multiqc']['report_exists'] = True

        # Check Salmon completion
        salmon_dir = results_dir / "salmon"
        if salmon_dir.exists():
            quant_dirs = [d for d in salmon_dir.iterdir() if d.is_dir()]
            self.stage_progress['salmon']['total_samples'] = self._get_total_samples()
            self.stage_progress['salmon']['samples'] = len(quant_dirs)

            if self.stage_progress['salmon']['samples'] > 0:
                progress = min(100, (self.stage_progress['salmon']['samples'] / self.stage_progress['salmon']['total_samples']) * 100)
                self.stage_progress['salmon']['progress'] = progress

                if progress >= 100:
                    self.stage_progress['salmon']['status'] = 'completed'

        # Check counts completion
        counts_file = results_dir / "counts" / "counts.tsv"
        if counts_file.exists():
            self.stage_progress['counts']['status'] = 'completed'
            self.stage_progress['counts']['progress'] = 100

        # Check DE results
        de_dir = results_dir / "de"
        if de_dir.exists():
            de_files = list(de_dir.glob("DE_*.tsv"))
            total_contrasts = self._get_total_contrasts()

            self.stage_progress['deseq2']['total_contrasts'] = total_contrasts
            self.stage_progress['deseq2']['contrasts'] = len(de_files)

            if len(de_files) > 0:
                progress = min(100, (len(de_files) / total_contrasts) * 100)
                self.stage_progress['deseq2']['progress'] = progress

                if progress >= 100:
                    self.stage_progress['deseq2']['status'] = 'completed'

        # Check FGSEA results
        fgsea_dir = results_dir / "fgsea"
        if fgsea_dir.exists():
            fgsea_files = list(fgsea_dir.glob("fgsea_*.tsv"))
            total_contrasts = self._get_total_contrasts()

            self.stage_progress['fgsea']['total_contrasts'] = total_contrasts
            self.stage_progress['fgsea']['contrasts'] = len(fgsea_files)

            if len(fgsea_files) > 0:
                progress = min(100, (len(fgsea_files) / total_contrasts) * 100)
                self.stage_progress['fgsea']['progress'] = progress

                if progress >= 100:
                    self.stage_progress['fgsea']['status'] = 'completed'

        # Check final report
        report_file = results_dir / "report.html"
        if report_file.exists():
            self.stage_progress['report']['status'] = 'completed'
            self.stage_progress['report']['progress'] = 100

    def _get_total_samples(self):
        """Get total number of samples from config."""
        try:
            samples_file = Path("config/samples.tsv")
            if samples_file.exists():
                import pandas as pd
                df = pd.read_csv(samples_file, sep='\t')
                return len(df)
        except:
            pass
        return 6  # Default for test data

    def _get_total_contrasts(self):
        """Get total number of contrasts from config."""
        try:
            contrasts_file = Path("config/contrasts.tsv")
            if contrasts_file.exists():
                import pandas as pd
                df = pd.read_csv(contrasts_file, sep='\t')
                return len(df)
        except:
            pass
        return 2  # Default for test data

    def update_stage_progress(self, stage, status):
        """Update progress for a specific stage."""
        if stage in self.stage_progress:
            self.stage_progress[stage]['status'] = status
            self.current_stage = stage

    def collect_resource_metrics(self):
        """Collect current system resource usage."""
        try:
            cpu_percent = psutil.cpu_percent(interval=1)
            memory_percent = psutil.virtual_memory().percent
            disk_percent = psutil.disk_usage('/').percent

            self.resource_history['cpu_percent'].append(cpu_percent)
            self.resource_history['memory_percent'].append(memory_percent)
            self.resource_history['disk_usage'].append(disk_percent)

        except Exception as e:
            print(f"Warning: Could not collect resource metrics: {e}")

    def print_progress_dashboard(self):
        """Print the progress dashboard."""
        os.system('clear' if os.name == 'posix' else 'cls')  # Clear screen

        elapsed_time = time.time() - self.start_time
        elapsed_str = str(timedelta(seconds=int(elapsed_time)))

        print("🚀 RNA-seq Pipeline Progress Dashboard")
        print("=" * 60)
        print(f"⏱️  Elapsed Time: {elapsed_str}")
        print(f"🔧 Pipeline Engine: {self.detect_pipeline_type().upper()}")
        print(f"📊 Current Stage: {self.current_stage or 'Unknown'}")
        print()

        # Stage progress
        print("📈 Stage Progress:")
        for stage, info in self.stage_progress.items():
            status_icon = {
                'pending': '⏳',
                'running': '🔄',
                'completed': '✅',
                'failed': '❌'
            }.get(info['status'], '❓')

            if stage in ['fastqc', 'salmon'] and info['total_samples'] > 0:
                progress_bar = self._create_progress_bar(info['progress'])
                print(f"   {status_icon} {stage.upper():<8}: {progress_bar} ({info['samples']}/{info['total_samples']} samples)")
            elif stage in ['deseq2', 'fgsea'] and info['total_contrasts'] > 0:
                progress_bar = self._create_progress_bar(info['progress'])
                print(f"   {status_icon} {stage.upper():<8}: {progress_bar} ({info['contrasts']}/{info['total_contrasts']} contrasts)")
            else:
                progress_bar = self._create_progress_bar(info['progress'])
                print(f"   {status_icon} {stage.upper():<8}: {progress_bar}")

        print()

        # Resource usage
        print("💻 System Resources:")
        if self.resource_history['cpu_percent']:
            avg_cpu = sum(self.resource_history['cpu_percent']) / len(self.resource_history['cpu_percent'])
            avg_memory = sum(self.resource_history['memory_percent']) / len(self.resource_history['memory_percent'])
            avg_disk = sum(self.resource_history['disk_usage']) / len(self.resource_history['disk_usage'])

            print(f"   CPU Usage:     {avg_cpu:.1f}%")
            print(f"   Memory Usage:  {avg_memory:.1f}%")
            print(f"   Disk Usage:    {avg_disk:.1f}%")

        print()

        # Quality metrics
        print("🔍 Quality Metrics:")
        print(f"   FastQC: {self.quality_metrics['fastqc']['passed_samples']}/{self.quality_metrics['fastqc']['total_samples']} samples passed")
        print(f"   MultiQC Report: {'✅ Generated' if self.quality_metrics['multiqc']['report_exists'] else '⏳ Pending'}")
        print(f"   Salmon: {self.quality_metrics['salmon']['successful_samples']}/{self.quality_metrics['salmon']['total_samples']} samples successful")
        print(f"   DESeq2: {self.quality_metrics['deseq2']['successful_contrasts']}/{self.quality_metrics['deseq2']['total_contrasts']} contrasts completed")
        print(f"   FGSEA: {self.quality_metrics['fgsea']['successful_contrasts']}/{self.quality_metrics['fgsea']['total_contrasts']} contrasts completed")

        print()
        print("💡 Tips:")
        print("   - Press Ctrl+C to exit monitoring")
        print("   - Check 'results/' directory for output files")
        print("   - Use 'make validate' to verify cross-engine reproducibility")

    def _create_progress_bar(self, progress):
        """Create a visual progress bar."""
        width = 20
        filled = int(width * progress / 100)
        bar = '█' * filled + '░' * (width - filled)
        return f"[{bar}] {progress:.1f}%"

    def run_monitoring_loop(self):
        """Run the main monitoring loop."""
        print("Starting pipeline monitoring... (Press Ctrl+C to stop)")

        try:
            while True:
                # Update information
                self.parse_snakemake_logs()
                self.parse_nextflow_logs()
                self.check_output_files()
                self.collect_resource_metrics()

                # Print dashboard
                self.print_progress_dashboard()

                # Check if pipeline is complete
                all_complete = all(
                    info['status'] in ['completed', 'failed']
                    for info in self.stage_progress.values()
                )

                if all_complete:
                    print("\n🎉 Pipeline execution finished!")
                    break

                time.sleep(self.update_interval)

        except KeyboardInterrupt:
            print("\n\n👋 Monitoring stopped by user")
            return True


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Monitor RNA-seq pipeline progress")
    parser.add_argument('--log-dir', '-l', default='logs', help='Log directory to monitor')
    parser.add_argument('--interval', '-i', type=int, default=2, help='Update interval in seconds')
    parser.add_argument('--once', action='store_true', help='Run once and exit')

    args = parser.parse_args()

    monitor = PipelineMonitor(args.log_dir, args.interval)

    if args.once:
        # Single update
        monitor.parse_snakemake_logs()
        monitor.parse_nextflow_logs()
        monitor.check_output_files()
        monitor.collect_resource_metrics()
        monitor.print_progress_dashboard()
    else:
        # Continuous monitoring
        monitor.run_monitoring_loop()


if __name__ == "__main__":
    main()
