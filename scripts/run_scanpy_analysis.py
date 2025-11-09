#!/usr/bin/env python3
"""
Comprehensive single-cell analysis script using Scanpy.
Performs filtering, normalization, dimensionality reduction, clustering, and marker gene detection.
"""
import argparse
import logging
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import celltypist

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def run_annotation(adata: ad.AnnData, model_name: str, output_dir: Path):
    """
    Performs automated cell type annotation using CellTypist.
    """
    logging.info(f"Performing cell type annotation with model: {model_name}")
    try:
        # Download the model if it's not already cached
        model_path = celltypist.models.get(model=model_name)
        
        # Predict cell types
        predictions = celltypist.annotate(adata, model=model_path, majority_voting=True)
        adata.obs['predicted_labels'] = predictions.predicted_labels['majority_voting']

        # Visualization
        fig, ax = plt.subplots(figsize=(10, 10))
        sc.pl.umap(adata, color="predicted_labels", ax=ax, show=False, legend_loc="on data")
        plt.tight_layout()
        plt.savefig(output_dir / "umap_cell_types.png")

        logging.info("Cell type annotation finished successfully.")
        
    except Exception as e:
        logging.error(f"CellTypist annotation failed: {e}")

def run_analysis(
    input_dir: Path,
    output_dir: Path,
    min_genes: int,
    min_cells: int,
    max_mito_percent: float,
    n_pcs: int,
    n_neighbors: int,
    resolution: float,
    run_celltypist: bool,
    celltypist_model: str,
):
    """Run the complete Scanpy analysis workflow."""
    logging.info(f"Starting analysis for data in {input_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load Data
    adata = sc.read_10x_mtx(input_dir, var_names="gene_symbols", cache=True)
    adata.var_names_make_unique()
    logging.info(f"Loaded initial matrix: {adata.n_obs} cells x {adata.n_vars} genes")

    # 2. Quality Control and Filtering
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show=False, ax=axs)
    plt.savefig(output_dir / "qc_violin_before_filtering.png")

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs.pct_counts_mt < max_mito_percent, :]
    logging.info(f"Matrix after filtering: {adata.n_obs} cells x {adata.n_vars} genes")

    # 3. Normalization and Scaling
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)

    # 4. Dimensionality Reduction
    sc.tl.pca(adata, svd_solver="arpack", n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)

    # 5. Clustering
    sc.tl.leiden(adata, resolution=resolution)

    # 6. Find Marker Genes
    sc.tl.rank_genes_groups(adata, "leiden", method="t-test")
    marker_df = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(20)
    marker_df.to_csv(output_dir / "marker_genes_top20.tsv", sep="	", index=False)

    # 7. Visualization
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    sc.pl.umap(adata, color="leiden", ax=axs[0], show=False, legend_loc="on data")
    sc.pl.umap(adata, color="total_counts", ax=axs[1], show=False)
    plt.tight_layout()
    plt.savefig(output_dir / "umap_clusters.png")

    # 8. Automated Annotation (optional)
    if run_celltypist:
        run_annotation(adata, celltypist_model, output_dir)

    # 9. Save final object
    adata.write(output_dir / "processed_data.h5ad")
    logging.info(f"Analysis complete. Results saved to {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Run Scanpy single-cell analysis workflow.")
    parser.add_argument("--input_dir", type=Path, required=True, help="Path to 10x formatted data directory (matrix.mtx.gz, etc.)")
    parser.add_argument("--output_dir", type=Path, required=True, help="Directory to save results.")
    parser.add_argument("--min_genes", type=int, default=200, help="Minimum genes per cell.")
    parser.add_argument("--min_cells", type=int, default=3, help="Minimum cells per gene.")
    parser.add_argument("--max_mito_percent", type=float, default=15.0, help="Maximum mitochondrial percentage.")
    parser.add_argument("--n_pcs", type=int, default=30, help="Number of principal components to use.")
    parser.add_argument("--n_neighbors", type=int, default=15, help="Number of neighbors for the graph.")
    parser.add_argument("--resolution", type=float, default=0.5, help="Resolution for Leiden clustering.")
    parser.add_argument("--run_annotation", action="store_true", help="If set, run CellTypist annotation.")
    parser.add_argument("--annotation_model", type=str, default="Immune_All_Low.pkl", help="Name of the CellTypist model to use.")
    args = parser.parse_args()

    run_analysis(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        max_mito_percent=args.max_mito_percent,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        resolution=args.resolution,
        run_celltypist=args.run_annotation,
        celltypist_model=args.annotation_model,
    )

if __name__ == "__main__":
    main()



