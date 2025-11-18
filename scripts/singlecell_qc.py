#!/usr/bin/env python3
"""
Single-cell RNA-seq quality control script for Snakemake integration.
"""

import os
import sys
import logging
from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_singlecell_data(matrix_file, barcodes_file):
    """Load single-cell count matrix."""
    try:
        # Load count matrix
        from scipy.io import mmread
        matrix = mmread(matrix_file).tocsr()

        # Load barcodes
        barcodes = pd.read_csv(barcodes_file, header=None, names=['barcode'])

        # Load genes (assuming genes.txt exists in same directory)
        genes_file = matrix_file.replace('counts.mtx', 'genes.txt')
        genes = pd.read_csv(genes_file, header=None, names=['gene_id', 'gene_symbol'])

        return matrix, barcodes, genes

    except Exception as e:
        logger.error(f"Error loading single-cell data: {e}")
        sys.exit(1)

def run_qc(matrix, barcodes, genes, output_dir, min_genes=200, max_genes=6000, mito_threshold=10.0):
    """Run quality control on single-cell data."""
    try:
        # Create AnnData object
        adata = sc.AnnData(X=matrix)
        adata.obs_names = barcodes['barcode'].values
        adata.var_names = genes['gene_symbol'].values

        # Calculate QC metrics
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

        # Filter cells based on gene counts
        logger.info(f"Filtering cells: min_genes={min_genes}, max_genes={max_genes}")

        # Mitochondrial gene filtering
        mito_genes = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
        adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1) * 100

        # Apply filters
        adata = adata[adata.obs.n_genes_by_counts.between(min_genes, max_genes), :]
        adata = adata[adata.obs.percent_mito < mito_threshold, :]

        logger.info(f"After QC: {adata.n_obs} cells, {adata.n_vars} genes")

        # Save filtered data
        os.makedirs(output_dir, exist_ok=True)

        # Save filtered matrix in various formats
        from scipy.io import mmwrite
        mmwrite(f"{output_dir}/filtered_counts.mtx", adata.X.T)

        # Save barcodes and genes
        adata.obs_names.to_series().to_csv(f"{output_dir}/filtered_barcodes.txt", index=False, header=False)
        adata.var_names.to_series().to_csv(f"{output_dir}/filtered_genes.txt", index=False, header=False)

        # Save QC statistics
        qc_stats = adata.obs.describe()
        qc_stats.to_csv(f"{output_dir}/qc_statistics.csv")

        # Create QC plots
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # Gene count distribution
        sns.histplot(adata.obs['n_genes_by_counts'], ax=axes[0, 0], bins=50)
        axes[0, 0].set_xlabel('Genes per cell')
        axes[0, 0].set_ylabel('Number of cells')
        axes[0, 0].set_title('Genes per Cell Distribution')

        # UMI count distribution
        sns.histplot(adata.obs['total_counts'], ax=axes[0, 1], bins=50)
        axes[0, 1].set_xlabel('UMIs per cell')
        axes[0, 1].set_ylabel('Number of cells')
        axes[0, 1].set_title('UMIs per Cell Distribution')

        # Mitochondrial percentage
        sns.histplot(adata.obs['percent_mito'], ax=axes[1, 0], bins=50)
        axes[1, 0].set_xlabel('Mitochondrial percentage')
        axes[1, 0].set_ylabel('Number of cells')
        axes[1, 0].set_title('Mitochondrial Percentage Distribution')

        # Violin plots
        qc_df = adata.obs[['n_genes_by_counts', 'total_counts', 'percent_mito']]
        qc_df.columns = ['Genes', 'UMIs', 'Mitochondrial %']
        sns.violinplot(data=qc_df, ax=axes[1, 1])
        axes[1, 1].set_title('QC Metrics Distribution')

        plt.tight_layout()
        plt.savefig(f"{output_dir}/qc_violin_plots.png", dpi=300, bbox_inches='tight')

        # Create HTML report
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Single-Cell QC Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .stats {{ background-color: #f5f5f5; padding: 15px; margin: 10px 0; border-radius: 5px; }}
                img {{ max-width: 100%; height: auto; margin: 10px 0; }}
            </style>
        </head>
        <body>
            <h1>Single-Cell RNA-seq Quality Control Report</h1>

            <div class="stats">
                <h2>Dataset Statistics</h2>
                <p><strong>Cells after filtering:</strong> {adata.n_obs:,}</p>
                <p><strong>Genes:</strong> {adata.n_vars:,}</p>
                <p><strong>Median genes per cell:</strong> {adata.obs['n_genes_by_counts'].median():.1f}</p>
                <p><strong>Median UMIs per cell:</strong> {adata.obs['total_counts'].median():.1f}</p>
            </div>

            <h2>Quality Control Plots</h2>
            <img src="qc_violin_plots.png" alt="QC Violin Plots">

            <div class="stats">
                <h2>Filtering Summary</h2>
                <p><strong>Min genes per cell:</strong> {min_genes}</p>
                <p><strong>Max genes per cell:</strong> {max_genes}</p>
                <p><strong>Mitochondrial threshold:</strong> {mito_threshold}%</p>
            </div>
        </body>
        </html>
        """

        with open(f"{output_dir}/qc_report.html", 'w') as f:
            f.write(html_content)

        logger.info(f"QC completed. Results saved to {output_dir}")

        return adata

    except Exception as e:
        logger.error(f"Error in QC: {e}")
        sys.exit(1)

def main():
    """Main function for Snakemake integration."""
    # Get input files from Snakemake
    matrix_file = snakemake.input.matrix
    barcodes_file = snakemake.input.barcodes
    output_dir = snakemake.params.output_dir
    min_genes = snakemake.params.min_genes
    max_genes = snakemake.params.max_genes
    mito_threshold = snakemake.params.mito_threshold

    logger.info(f"Starting single-cell QC for {matrix_file}")

    # Load data
    matrix, barcodes, genes = load_singlecell_data(matrix_file, barcodes_file)

    # Run QC
    adata = run_qc(matrix, barcodes, genes, output_dir, min_genes, max_genes, mito_threshold)

    logger.info("Single-cell QC completed successfully")

if __name__ == "__main__":
    main()











