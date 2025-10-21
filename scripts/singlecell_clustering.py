#!/usr/bin/env python3
"""
Single-cell RNA-seq clustering script for Snakemake integration.
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

def load_filtered_data(matrix_file, barcodes_file, genes_file):
    """Load filtered single-cell data."""
    try:
        # Load count matrix
        from scipy.io import mmread
        matrix = mmread(matrix_file).tocsr()

        # Load barcodes and genes
        barcodes = pd.read_csv(barcodes_file, header=None, names=['barcode'])
        genes = pd.read_csv(genes_file, header=None, names=['gene_id', 'gene_symbol'])

        # Create AnnData object
        adata = sc.AnnData(X=matrix)
        adata.obs_names = barcodes['barcode'].values
        adata.var_names = genes['gene_symbol'].values

        logger.info(f"Loaded data: {adata.n_obs} cells, {adata.n_vars} genes")
        return adata

    except Exception as e:
        logger.error(f"Error loading filtered data: {e}")
        sys.exit(1)

def run_clustering(adata, output_dir, normalization="log_normalize", n_pcs=50,
                  n_neighbors=10, resolution=1.0, method="leiden"):
    """Run clustering analysis on single-cell data."""
    try:
        # Normalization
        logger.info(f"Normalizing data using {normalization} method")
        if normalization == "log_normalize":
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        elif normalization == "scran":
            # Would need scran R package for this
            logger.warning("SCRAN normalization not implemented, using log normalization")
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        elif normalization == "sct":
            # Would need sctransform R package for this
            logger.warning("SCT normalization not implemented, using log normalization")
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)

        # Find highly variable genes
        logger.info("Finding highly variable genes")
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata.var.highly_variable]

        # Scale data
        logger.info("Scaling data")
        sc.pp.scale(adata, max_value=10)

        # PCA
        logger.info(f"Running PCA with {n_pcs} components")
        sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)

        # Compute neighborhood graph
        logger.info(f"Computing neighborhood graph with {n_neighbors} neighbors")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

        # UMAP
        logger.info("Computing UMAP")
        sc.tl.umap(adata)

        # t-SNE
        logger.info("Computing t-SNE")
        sc.tl.tsne(adata, n_pcs=n_pcs)

        # Clustering
        logger.info(f"Running {method} clustering with resolution {resolution}")
        if method == "leiden":
            sc.tl.leiden(adata, resolution=resolution)
        elif method == "louvain":
            sc.tl.louvain(adata, resolution=resolution)
        else:
            logger.warning(f"Unknown clustering method {method}, using leiden")
            sc.tl.leiden(adata, resolution=resolution)

        # Find marker genes
        logger.info("Finding marker genes")
        sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

        # Save results
        os.makedirs(output_dir, exist_ok=True)

        # Save AnnData object
        adata.write(f"{output_dir}/clustering_results.h5ad")

        # Save clustering results in RDS-like format (for compatibility)
        import pickle
        with open(f"{output_dir}/clustering_results.rds", 'wb') as f:
            pickle.dump({
                'adata': adata,
                'parameters': {
                    'normalization': normalization,
                    'n_pcs': n_pcs,
                    'n_neighbors': n_neighbors,
                    'resolution': resolution,
                    'method': method
                }
            }, f)

        # Create visualization plots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))

        # UMAP colored by clusters
        sc.pl.umap(adata, color='leiden', ax=axes[0, 0], show=False)
        axes[0, 0].set_title('UMAP - Leiden Clusters')

        # t-SNE colored by clusters
        sc.pl.tsne(adata, color='leiden', ax=axes[0, 1], show=False)
        axes[0, 1].set_title('t-SNE - Leiden Clusters')

        # PCA colored by clusters
        sc.pl.pca(adata, color='leiden', ax=axes[1, 0], show=False)
        axes[1, 0].set_title('PCA - Leiden Clusters')

        # QC metrics on UMAP
        sc.pl.umap(adata, color='total_counts', ax=axes[1, 1], show=False)
        axes[1, 1].set_title('UMAP - Total Counts')

        plt.tight_layout()
        plt.savefig(f"{output_dir}/clustering_plots.png", dpi=300, bbox_inches='tight')

        # Create UMAP plot
        plt.figure(figsize=(8, 6))
        sc.pl.umap(adata, color='leiden', show=False)
        plt.title('UMAP - Leiden Clusters')
        plt.savefig(f"{output_dir}/umap_clusters.png", dpi=300, bbox_inches='tight')

        # Create t-SNE plot
        plt.figure(figsize=(8, 6))
        sc.pl.tsne(adata, color='leiden', show=False)
        plt.title('t-SNE - Leiden Clusters')
        plt.savefig(f"{output_dir}/tsne_clusters.png", dpi=300, bbox_inches='tight')

        # Save cluster statistics
        cluster_stats = adata.obs['leiden'].value_counts().sort_index()
        cluster_stats.to_csv(f"{output_dir}/cluster_sizes.csv")

        # Save top marker genes per cluster
        marker_genes = sc.get.rank_genes_groups_df(adata, None)
        marker_genes.to_csv(f"{output_dir}/marker_genes.csv", index=False)

        logger.info(f"Clustering completed. Results saved to {output_dir}")
        logger.info(f"Found {len(cluster_stats)} clusters")

        return adata

    except Exception as e:
        logger.error(f"Error in clustering: {e}")
        sys.exit(1)

def main():
    """Main function for Snakemake integration."""
    # Get input files from Snakemake
    matrix_file = snakemake.input.matrix
    barcodes_file = snakemake.input.barcodes
    genes_file = snakemake.input.genes
    output_dir = snakemake.params.output_dir
    normalization = snakemake.params.normalization
    n_pcs = snakemake.params.n_pcs
    resolution = snakemake.params.resolution
    method = snakemake.params.method

    logger.info(f"Starting single-cell clustering for {matrix_file}")

    # Load data
    adata = load_filtered_data(matrix_file, barcodes_file, genes_file)

    # Run clustering
    adata = run_clustering(adata, output_dir, normalization, n_pcs, n_pcs, resolution, method)

    logger.info("Single-cell clustering completed successfully")

if __name__ == "__main__":
    main()



