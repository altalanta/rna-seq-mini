#!/usr/bin/env python3
"""
Single-cell quantification module for RNASEQ-MINI.
Supports scRNA-seq, scATAC-seq, and spatial transcriptomics quantification.
"""

import os
import subprocess
import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class SingleCellConfig:
    """Configuration for single-cell analysis."""

    # Input files
    fastq_files: List[str]
    barcode_file: Optional[str] = None
    feature_file: Optional[str] = None

    # Technology-specific parameters
    technology: str = "10x"  # 10x, dropseq, smartseq, etc.
    chemistry: str = "auto"  # auto, v2, v3, etc.

    # Quantification parameters
    min_umi: int = 500
    min_genes: int = 200
    max_mito_percent: float = 10.0

    # Output parameters
    output_dir: str = "results/singlecell"
    temp_dir: str = "temp"

    # Performance parameters
    threads: int = 8
    memory_gb: int = 32

    # Quality control parameters
    adapter_trimming: bool = True
    poly_a_trimming: bool = True
    quality_filtering: bool = True

    # Advanced parameters
    include_introns: bool = True
    correct_chemistry: bool = True
    estimate_complexity: bool = True


class SingleCellQuantifier:
    """Handles single-cell RNA-seq quantification."""

    def __init__(self, config: SingleCellConfig):
        self.config = config
        self.output_dir = Path(config.output_dir)
        self.temp_dir = Path(config.temp_dir)

        # Create directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.temp_dir.mkdir(parents=True, exist_ok=True)

        # Tool availability
        self.tools = self._check_tool_availability()

    def _check_tool_availability(self) -> Dict[str, bool]:
        """Check availability of required tools."""
        tools = {}

        # Check for cellranger (10x Genomics)
        try:
            result = subprocess.run(['cellranger', '--version'],
                                  capture_output=True, text=True, timeout=10)
            tools['cellranger'] = result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            tools['cellranger'] = False

        # Check for kallisto|bustools
        try:
            result = subprocess.run(['kallisto', 'version'],
                                  capture_output=True, text=True, timeout=10)
            tools['kallisto'] = result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            tools['kallisto'] = False

        # Check for STARsolo
        try:
            result = subprocess.run(['STAR', '--version'],
                                  capture_output=True, text=True, timeout=10)
            tools['starsolo'] = result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            tools['starsolo'] = False

        # Check for alevin-fry
        try:
            result = subprocess.run(['salmon', 'alevin', '--help'],
                                  capture_output=True, text=True, timeout=10)
            tools['alevin'] = result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            tools['alevin'] = False

        return tools

    def quantify_10x(self, sample_name: str, fastq_r1: str, fastq_r2: str,
                    reference: str, chemistry: str = "auto") -> Dict[str, Any]:
        """
        Quantify 10x Genomics single-cell data using CellRanger.

        Args:
            sample_name: Name of the sample
            fastq_r1: Path to R1 FASTQ file
            fastq_r2: Path to R2 FASTQ file
            reference: Path to reference genome/transcriptome
            chemistry: Chemistry version (v2, v3, auto)

        Returns:
            Quantification results summary
        """
        if not self.tools.get('cellranger', False):
            raise RuntimeError("CellRanger not available. Please install 10x Genomics CellRanger.")

        sample_dir = self.output_dir / sample_name
        sample_dir.mkdir(exist_ok=True)

        # Build CellRanger command
        cmd = [
            'cellranger', 'count',
            '--id', sample_name,
            '--transcriptome', reference,
            '--fastqs', str(Path(fastq_r1).parent),
            '--sample', sample_name,
            '--localcores', str(self.config.threads),
            '--localmem', str(self.config.memory_gb)
        ]

        # Add chemistry if specified
        if chemistry != "auto":
            cmd.extend(['--chemistry', chemistry])

        # Add include-introns if requested
        if self.config.include_introns:
            cmd.append('--include-introns')

        logger.info(f"Running CellRanger count for sample {sample_name}")
        logger.info(f"Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, cwd=str(sample_dir), check=True,
                                  capture_output=True, text=True)

            # Parse results
            results = self._parse_cellranger_results(sample_dir / sample_name)

            logger.info(f"CellRanger quantification completed for {sample_name}")
            return results

        except subprocess.CalledProcessError as e:
            logger.error(f"CellRanger failed for {sample_name}: {e.stderr}")
            raise RuntimeError(f"CellRanger quantification failed: {e.stderr}")

    def quantify_kallisto_bustools(self, sample_name: str, fastq_r1: str, fastq_r2: str,
                                 reference: str, technology: str = "10xV3") -> Dict[str, Any]:
        """
        Quantify single-cell data using kallisto|bustools.

        Args:
            sample_name: Name of the sample
            fastq_r1: Path to R1 FASTQ file
            fastq_r2: Path to R2 FASTQ file
            reference: Path to kallisto index
            technology: Technology type (10xV2, 10xV3, DROPSEQ, etc.)

        Returns:
            Quantification results summary
        """
        if not self.tools.get('kallisto', False):
            raise RuntimeError("kallisto not available. Please install kallisto.")

        sample_dir = self.output_dir / sample_name
        sample_dir.mkdir(exist_ok=True)

        # Create kallisto|bustools workflow
        cmd = [
            'kallisto', 'bus',
            '-i', reference,
            '-o', str(sample_dir / 'bus_output'),
            '-x', technology,
            '-t', str(self.config.threads),
            fastq_r1, fastq_r2
        ]

        logger.info(f"Running kallisto bus for sample {sample_name}")

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Run bustools for matrix generation
            bustools_cmd = [
                'bustools', 'count',
                '--output', str(sample_dir / 'counts'),
                '--genecounts',
                '--cm', str(sample_dir / 'bus_output' / 'matrix.ec'),
                '--ec', str(sample_dir / 'bus_output' / 'transcripts.txt'),
                '--txnames', str(sample_dir / 'bus_output' / 'transcripts.txt'),
                str(sample_dir / 'bus_output' / 'bus_output.bus')
            ]

            subprocess.run(bustools_cmd, check=True, capture_output=True, text=True)

            # Parse results
            results = self._parse_kallisto_results(sample_dir)

            logger.info(f"kallisto|bustools quantification completed for {sample_name}")
            return results

        except subprocess.CalledProcessError as e:
            logger.error(f"kallisto|bustools failed for {sample_name}: {e.stderr}")
            raise RuntimeError(f"kallisto|bustools quantification failed: {e.stderr}")

    def quantify_starsolo(self, sample_name: str, fastq_r1: str, fastq_r2: str,
                         reference: str, solo_type: str = "CB_UMI_Simple") -> Dict[str, Any]:
        """
        Quantify single-cell data using STARsolo.

        Args:
            sample_name: Name of the sample
            fastq_r1: Path to R1 FASTQ file
            fastq_r2: Path to R2 FASTQ file
            reference: Path to STAR index
            solo_type: SOLO type (CB_UMI_Simple, CB_UMI_Complex, etc.)

        Returns:
            Quantification results summary
        """
        if not self.tools.get('starsolo', False):
            raise RuntimeError("STAR not available. Please install STAR aligner.")

        sample_dir = self.output_dir / sample_name
        sample_dir.mkdir(exist_ok=True)

        # Create STARsolo command
        cmd = [
            'STAR',
            '--genomeDir', reference,
            '--readFilesIn', fastq_r2, fastq_r1,  # STAR expects R2 first for 10x
            '--runThreadN', str(self.config.threads),
            '--runDirPerm', 'All_RWX',
            '--outFileNamePrefix', str(sample_dir / ''),
            '--soloType', solo_type,
            '--soloUMIlen', '12',  # Standard for 10x
            '--soloCBlen', '16',   # Standard for 10x
            '--soloBarcodeReadLength', '0',  # Auto-detect
            '--soloFeatures', 'Gene', 'Velocyto'
        ]

        # Add additional parameters for quality
        cmd.extend([
            '--soloCellFilter', 'EmptyDrops_CR',
            '--soloCBmatchWLtype', '1MM_multi',
            '--soloUMIfiltering', 'MultiGeneUMI_CR',
            '--soloUMIdedup', '1MM_CR'
        ])

        logger.info(f"Running STARsolo for sample {sample_name}")

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Parse results
            results = self._parse_starsolo_results(sample_dir)

            logger.info(f"STARsolo quantification completed for {sample_name}")
            return results

        except subprocess.CalledProcessError as e:
            logger.error(f"STARsolo failed for {sample_name}: {e.stderr}")
            raise RuntimeError(f"STARsolo quantification failed: {e.stderr}")

    def quantify_alevin(self, sample_name: str, fastq_r1: str, fastq_r2: str,
                       reference: str, technology: str = "10xv3") -> Dict[str, Any]:
        """
        Quantify single-cell data using alevin-fry.

        Args:
            sample_name: Name of the sample
            fastq_r1: Path to R1 FASTQ file
            fastq_r2: Path to R2 FASTQ file
            reference: Path to salmon index
            technology: Technology type (10xv2, 10xv3, dropseq, etc.)

        Returns:
            Quantification results summary
        """
        if not self.tools.get('alevin', False):
            raise RuntimeError("alevin-fry not available. Please install salmon with alevin support.")

        sample_dir = self.output_dir / sample_name
        sample_dir.mkdir(exist_ok=True)

        # Create alevin command
        cmd = [
            'salmon', 'alevin',
            '-i', reference,
            '-l', 'ISR',  # Library type for 10x
            '-1', fastq_r1,
            '-2', fastq_r2,
            '-o', str(sample_dir / 'alevin_output'),
            '--tgMap', str(sample_dir / 'alevin_output' / 'transcript_gene_map.txt'),
            '--chromiumV3' if technology == '10xv3' else '--chromium',
            '-p', str(self.config.threads)
        ]

        logger.info(f"Running alevin for sample {sample_name}")

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Run alevin-fry for quantification
            fry_cmd = [
                'alevin-fry', 'generate-permit-list',
                '-i', str(sample_dir / 'alevin_output'),
                '-o', str(sample_dir / 'fry_output'),
                '--expected-ori', 'fw'
            ]

            subprocess.run(fry_cmd, check=True, capture_output=True, text=True)

            # Parse results
            results = self._parse_alevin_results(sample_dir)

            logger.info(f"alevin-fry quantification completed for {sample_name}")
            return results

        except subprocess.CalledProcessError as e:
            logger.error(f"alevin-fry failed for {sample_name}: {e.stderr}")
            raise RuntimeError(f"alevin-fry quantification failed: {e.stderr}")

    def _parse_cellranger_results(self, sample_dir: Path) -> Dict[str, Any]:
        """Parse CellRanger results."""
        try:
            # Read metrics summary
            metrics_file = sample_dir / "outs" / "metrics_summary.csv"
            if metrics_file.exists():
                metrics_df = pd.read_csv(metrics_file)

                # Extract key metrics
                metrics = {
                    'estimated_cells': int(metrics_df['Estimated Number of Cells'].iloc[0]),
                    'mean_reads_per_cell': int(metrics_df['Mean Reads per Cell'].iloc[0]),
                    'median_genes_per_cell': int(metrics_df['Median Genes per Cell'].iloc[0]),
                    'total_genes_detected': int(metrics_df['Number of Genes'].iloc[0]),
                    'median_umi_counts_per_cell': int(metrics_df['Median UMI Counts per Cell'].iloc[0])
                }
            else:
                metrics = {}

            # Check for output files
            output_files = {
                'matrix': (sample_dir / "outs" / "filtered_feature_bc_matrix").exists(),
                'raw_matrix': (sample_dir / "outs" / "raw_feature_bc_matrix").exists(),
                'barcodes': (sample_dir / "outs" / "filtered_feature_bc_matrix" / "barcodes.tsv.gz").exists(),
                'features': (sample_dir / "outs" / "filtered_feature_bc_matrix" / "features.tsv.gz").exists(),
                'matrix_file': (sample_dir / "outs" / "filtered_feature_bc_matrix" / "matrix.mtx.gz").exists()
            }

            return {
                'sample': sample_dir.name,
                'tool': 'cellranger',
                'metrics': metrics,
                'output_files': output_files,
                'output_directory': str(sample_dir / "outs"),
                'success': all(output_files.values())
            }

        except Exception as e:
            logger.error(f"Error parsing CellRanger results: {e}")
            return {'sample': sample_dir.name, 'tool': 'cellranger', 'error': str(e)}

    def _parse_kallisto_results(self, sample_dir: Path) -> Dict[str, Any]:
        """Parse kallisto|bustools results."""
        try:
            # Check for output files
            output_files = {
                'counts': (sample_dir / "counts" / "genes.txt").exists(),
                'cells': (sample_dir / "counts" / "cells_x_genes.mtx").exists(),
                'bus_output': (sample_dir / "bus_output").exists()
            }

            # Try to read counts
            metrics = {}
            if output_files['counts']:
                try:
                    counts_df = pd.read_csv(sample_dir / "counts" / "genes.txt", sep='\t')
                    metrics = {
                        'genes_detected': len(counts_df),
                        'total_umis': counts_df['Count'].sum(),
                        'median_genes_per_cell': counts_df['Count'].median()
                    }
                except Exception as e:
                    logger.warning(f"Error reading kallisto counts: {e}")

            return {
                'sample': sample_dir.name,
                'tool': 'kallisto_bustools',
                'metrics': metrics,
                'output_files': output_files,
                'output_directory': str(sample_dir / "counts"),
                'success': output_files['counts']
            }

        except Exception as e:
            logger.error(f"Error parsing kallisto results: {e}")
            return {'sample': sample_dir.name, 'tool': 'kallisto_bustools', 'error': str(e)}

    def _parse_starsolo_results(self, sample_dir: Path) -> Dict[str, Any]:
        """Parse STARsolo results."""
        try:
            # Check for output files
            output_files = {
                'matrix': (sample_dir / "Solo.out" / "Gene" / "filtered").exists(),
                'raw_matrix': (sample_dir / "Solo.out" / "Gene" / "raw").exists(),
                'barcodes': (sample_dir / "Solo.out" / "Gene" / "filtered" / "barcodes.tsv").exists(),
                'features': (sample_dir / "Solo.out" / "Gene" / "filtered" / "features.tsv").exists(),
                'matrix_file': (sample_dir / "Solo.out" / "Gene" / "filtered" / "matrix.mtx").exists()
            }

            # Try to read summary stats
            metrics = {}
            summary_file = sample_dir / "Solo.out" / "Gene" / "Summary.csv"
            if summary_file.exists():
                try:
                    summary_df = pd.read_csv(summary_file)
                    metrics = {
                        'estimated_cells': int(summary_df['Estimated Number of Cells'].iloc[0]),
                        'mean_reads_per_cell': float(summary_df['Mean Reads per Cell'].iloc[0]),
                        'median_genes_per_cell': float(summary_df['Median Genes per Cell'].iloc[0]),
                        'fraction_reads_in_cells': float(summary_df['Fraction Reads in Cells'].iloc[0])
                    }
                except Exception as e:
                    logger.warning(f"Error reading STARsolo summary: {e}")

            return {
                'sample': sample_dir.name,
                'tool': 'starsolo',
                'metrics': metrics,
                'output_files': output_files,
                'output_directory': str(sample_dir / "Solo.out" / "Gene"),
                'success': output_files['matrix']
            }

        except Exception as e:
            logger.error(f"Error parsing STARsolo results: {e}")
            return {'sample': sample_dir.name, 'tool': 'starsolo', 'error': str(e)}

    def _parse_alevin_results(self, sample_dir: Path) -> Dict[str, Any]:
        """Parse alevin-fry results."""
        try:
            # Check for output files
            output_files = {
                'quant': (sample_dir / "fry_output" / "quant").exists(),
                'permit_list': (sample_dir / "fry_output" / "permit_list.txt").exists()
            }

            # Try to read quantification results
            metrics = {}
            if output_files['quant']:
                try:
                    # Read alevin-fry output
                    quant_dir = sample_dir / "fry_output" / "quant"
                    if (quant_dir / "quants_mat.gz").exists():
                        # This would require additional parsing of alevin-fry specific formats
                        metrics = {
                            'quantification_successful': True,
                            'output_format': 'alevin-fry'
                        }
                except Exception as e:
                    logger.warning(f"Error reading alevin-fry results: {e}")

            return {
                'sample': sample_dir.name,
                'tool': 'alevin_fry',
                'metrics': metrics,
                'output_files': output_files,
                'output_directory': str(sample_dir / "fry_output"),
                'success': output_files['quant']
            }

        except Exception as e:
            logger.error(f"Error parsing alevin-fry results: {e}")
            return {'sample': sample_dir.name, 'tool': 'alevin_fry', 'error': str(e)}

    def generate_qc_report(self, sample_results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Generate comprehensive QC report for all samples."""
        qc_report = {
            'summary': {
                'total_samples': len(sample_results),
                'successful_quantifications': len([r for r in sample_results if r.get('success', False)]),
                'failed_quantifications': len([r for r in sample_results if not r.get('success', False)])
            },
            'samples': {},
            'comparative_metrics': {},
            'recommendations': []
        }

        # Process each sample
        for result in sample_results:
            sample_name = result['sample']
            qc_report['samples'][sample_name] = result

        # Calculate comparative metrics
        if len(sample_results) > 1:
            qc_report['comparative_metrics'] = self._calculate_comparative_metrics(sample_results)

        # Generate recommendations
        qc_report['recommendations'] = self._generate_qc_recommendations(sample_results)

        return qc_report

    def _calculate_comparative_metrics(self, sample_results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Calculate comparative metrics across samples."""
        metrics = {}

        # Extract metrics for comparison
        successful_results = [r for r in sample_results if r.get('success', False)]

        if successful_results:
            # Cell counts comparison
            cell_counts = [r['metrics'].get('estimated_cells', 0) for r in successful_results]
            metrics['cell_counts'] = {
                'mean': np.mean(cell_counts),
                'std': np.std(cell_counts),
                'min': np.min(cell_counts),
                'max': np.max(cell_counts)
            }

            # Reads per cell comparison
            reads_per_cell = [r['metrics'].get('mean_reads_per_cell', 0) for r in successful_results]
            metrics['reads_per_cell'] = {
                'mean': np.mean(reads_per_cell),
                'std': np.std(reads_per_cell),
                'min': np.min(reads_per_cell),
                'max': np.max(reads_per_cell)
            }

        return metrics

    def _generate_qc_recommendations(self, sample_results: List[Dict[str, Any]]) -> List[str]:
        """Generate QC recommendations based on results."""
        recommendations = []

        successful_results = [r for r in sample_results if r.get('success', False)]

        if len(successful_results) == 0:
            recommendations.append("No samples quantified successfully. Check input files and parameters.")
            return recommendations

        # Check cell count consistency
        cell_counts = [r['metrics'].get('estimated_cells', 0) for r in successful_results]
        if len(set(cell_counts)) > 1:
            recommendations.append("Inconsistent cell counts detected across samples. Check sample quality.")

        # Check reads per cell
        reads_per_cell = [r['metrics'].get('mean_reads_per_cell', 0) for r in successful_results]
        if any(rpc < 1000 for rpc in reads_per_cell):
            recommendations.append("Low reads per cell detected. Consider increasing sequencing depth.")

        # Check gene detection
        genes_per_cell = [r['metrics'].get('median_genes_per_cell', 0) for r in successful_results]
        if any(gpc < 500 for gpc in genes_per_cell):
            recommendations.append("Low gene detection per cell. Check library preparation quality.")

        if not recommendations:
            recommendations.append("All QC metrics within acceptable ranges.")

        return recommendations

    def create_multi_sample_matrix(self, sample_results: List[Dict[str, Any]],
                                 output_file: str = "combined_counts.h5ad") -> str:
        """
        Combine multiple single-cell samples into a single matrix.

        Args:
            sample_results: List of quantification results
            output_file: Output file path

        Returns:
            Path to combined matrix file
        """
        try:
            import anndata as ad

            # This would require additional implementation to combine matrices
            # For now, return placeholder
            output_path = self.output_dir / output_file

            logger.info(f"Combined matrix would be created at: {output_path}")
            return str(output_path)

        except ImportError:
            logger.warning("anndata not available for matrix combination")
            return ""










