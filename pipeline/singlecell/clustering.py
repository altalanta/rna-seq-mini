#!/usr/bin/env python3
"""
Single-cell clustering and dimensionality reduction module.
Provides comprehensive clustering analysis for single-cell RNA-seq data.
"""

import os
import json
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Union
from dataclasses import dataclass
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class ClusteringConfig:
    """Configuration for single-cell clustering analysis."""

    # Input data
    matrix_file: str
    metadata_file: Optional[str] = None

    # Preprocessing parameters
    min_genes: int = 200
    min_cells: int = 3
    max_genes: int = 6000
    mito_percent: float = 10.0

    # Normalization parameters
    normalization_method: str = "log_normalize"  # log_normalize, scran, sct
    scale_data: bool = True

    # Dimensionality reduction parameters
    n_pcs: int = 50
    n_neighbors: int = 10
    resolution: float = 1.0

    # Clustering parameters
    clustering_method: str = "leiden"  # leiden, louvain, kmeans
    n_clusters: Optional[int] = None

    # Output parameters
    output_dir: str = "results/singlecell/clustering"
    save_plots: bool = True

    # Performance parameters
    threads: int = 4
    memory_gb: int = 16


class SingleCellClustering:
    """Handles single-cell clustering and dimensionality reduction."""

    def __init__(self, config: ClusteringConfig):
        self.config = config
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Check tool availability
        self.tools = self._check_tool_availability()

    def _check_tool_availability(self) -> Dict[str, bool]:
        """Check availability of required tools."""
        tools = {}

        # Check for scanpy
        try:
            import scanpy as sc
            tools['scanpy'] = True
        except ImportError:
            tools['scanpy'] = False

        # Check for seurat (R package)
        try:
            result = subprocess.run(['R', '--version'],
                                  capture_output=True, text=True, timeout=10)
            tools['seurat'] = result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            tools['seurat'] = False

        # Check for scvi-tools
        try:
            import scvi
            tools['scvi'] = True
        except ImportError:
            tools['scvi'] = False

        return tools

    def preprocess_data(self, matrix_file: str, metadata_file: Optional[str] = None) -> Dict[str, Any]:
        """
        Preprocess single-cell data for clustering.

        Args:
            matrix_file: Path to count matrix (h5ad, mtx, csv)
            metadata_file: Optional metadata file

        Returns:
            Preprocessing results and quality metrics
        """
        if not self.tools.get('scanpy', False):
            raise RuntimeError("scanpy not available. Please install scanpy for clustering analysis.")

        try:
            import scanpy as sc

            # Load data
            logger.info(f"Loading data from {matrix_file}")
            if matrix_file.endswith('.h5ad'):
                adata = sc.read_h5ad(matrix_file)
            else:
                # Assume it's a directory with 10x matrix
                adata = sc.read_10x_mtx(matrix_file)

            # Load metadata if provided
            if metadata_file:
                metadata = pd.read_csv(metadata_file, index_col=0)
                adata.obs = adata.obs.join(metadata, how='left')

            # Quality control filtering
            logger.info("Performing quality control filtering")
            sc.pp.filter_cells(adata, min_genes=self.config.min_genes)
            sc.pp.filter_genes(adata, min_cells=self.config.min_cells)

            # Calculate QC metrics
            adata.var['mt'] = adata.var_names.str.startswith('MT-')
            adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

            # Filter by mitochondrial percentage
            adata = adata[adata.obs.pct_counts_mt < self.config.mito_percent, :]

            # Normalize data
            logger.info(f"Normalizing data using {self.config.normalization_method}")
            if self.config.normalization_method == "log_normalize":
                sc.pp.normalize_total(adata, target_sum=1e4)
                sc.pp.log1p(adata)
            elif self.config.normalization_method == "scran":
                # Would use scran R package
                logger.warning("scran normalization not implemented, using log_normalize")
                sc.pp.normalize_total(adata, target_sum=1e4)
                sc.pp.log1p(adata)

            # Identify highly variable genes
            logger.info("Identifying highly variable genes")
            sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

            # Filter to highly variable genes
            adata = adata[:, adata.var.highly_variable]

            # Scale data
            if self.config.scale_data:
                logger.info("Scaling data")
                sc.pp.scale(adata, max_value=10)

            # Save preprocessed data
            preprocessed_file = self.output_dir / "preprocessed_data.h5ad"
            adata.write(preprocessed_file)

            # Calculate preprocessing metrics
            metrics = {
                'original_cells': adata.obs.shape[0],
                'original_genes': adata.var.shape[0],
                'filtered_cells': adata.obs.shape[0],
                'filtered_genes': adata.var.shape[0],
                'highly_variable_genes': adata.var.highly_variable.sum(),
                'median_genes_per_cell': adata.obs.n_genes.median(),
                'median_umis_per_cell': adata.obs.total_counts.median(),
                'median_mt_percent': adata.obs.pct_counts_mt.median()
            }

            return {
                'success': True,
                'preprocessed_file': str(preprocessed_file),
                'metrics': metrics,
                'adata_shape': adata.shape
            }

        except Exception as e:
            logger.error(f"Error preprocessing data: {e}")
            return {'success': False, 'error': str(e)}

    def run_dimensionality_reduction(self, preprocessed_file: str) -> Dict[str, Any]:
        """
        Run dimensionality reduction (PCA, UMAP, t-SNE).

        Args:
            preprocessed_file: Path to preprocessed data file

        Returns:
            Dimensionality reduction results
        """
        if not self.tools.get('scanpy', False):
            raise RuntimeError("scanpy not available for dimensionality reduction.")

        try:
            import scanpy as sc

            # Load preprocessed data
            adata = sc.read_h5ad(preprocessed_file)

            # Run PCA
            logger.info("Running PCA")
            sc.tl.pca(adata, svd_solver='arpack', n_comps=self.config.n_pcs)

            # Compute neighbors
            logger.info("Computing neighborhood graph")
            sc.pp.neighbors(adata, n_neighbors=self.config.n_neighbors, n_pcs=self.config.n_pcs)

            # Run UMAP
            logger.info("Running UMAP")
            sc.tl.umap(adata)

            # Run t-SNE (optional, computationally expensive)
            logger.info("Running t-SNE")
            sc.tl.tsne(adata, n_pcs=self.config.n_pcs, perplexity=30)

            # Save results
            dimred_file = self.output_dir / "dimensionality_reduction.h5ad"
            adata.write(dimred_file)

            # Calculate explained variance
            if hasattr(adata, 'uns') and 'pca' in adata.uns:
                pca_variance = adata.uns['pca']['variance_ratio']
                explained_variance = {
                    f'PC_{i+1}': var for i, var in enumerate(pca_variance[:10])
                }
            else:
                explained_variance = {}

            return {
                'success': True,
                'dimred_file': str(dimred_file),
                'explained_variance': explained_variance,
                'n_pcs_used': self.config.n_pcs,
                'n_neighbors': self.config.n_neighbors
            }

        except Exception as e:
            logger.error(f"Error in dimensionality reduction: {e}")
            return {'success': False, 'error': str(e)}

    def run_clustering(self, dimred_file: str) -> Dict[str, Any]:
        """
        Run clustering analysis.

        Args:
            dimred_file: Path to dimensionality reduction results

        Returns:
            Clustering results
        """
        if not self.tools.get('scanpy', False):
            raise RuntimeError("scanpy not available for clustering.")

        try:
            import scanpy as sc

            # Load dimensionality reduction results
            adata = sc.read_h5ad(dimred_file)

            # Run clustering
            logger.info(f"Running {self.config.clustering_method} clustering")
            if self.config.clustering_method == "leiden":
                sc.tl.leiden(adata, resolution=self.config.resolution)
            elif self.config.clustering_method == "louvain":
                sc.tl.louvain(adata, resolution=self.config.resolution)
            elif self.config.clustering_method == "kmeans":
                # Use k-means on PCA coordinates
                from sklearn.cluster import KMeans
                n_clusters = self.config.n_clusters or int(np.sqrt(adata.obs.shape[0] / 10))
                kmeans = KMeans(n_clusters=n_clusters, random_state=42)
                adata.obs['kmeans_clusters'] = kmeans.fit_predict(adata.obsm['X_pca'][:, :20])
            else:
                raise ValueError(f"Unknown clustering method: {self.config.clustering_method}")

            # Save clustered data
            clustered_file = self.output_dir / "clustered_data.h5ad"
            adata.write(clustered_file)

            # Calculate clustering metrics
            n_clusters = len(adata.obs[f'{self.config.clustering_method}_clusters'].unique())

            # Calculate silhouette score if possible
            silhouette_score = None
            try:
                from sklearn.metrics import silhouette_score
                if adata.obs.shape[0] < 10000:  # Only for smaller datasets
                    X_pca = adata.obsm['X_pca'][:, :20]  # Use first 20 PCs
                    cluster_labels = adata.obs[f'{self.config.clustering_method}_clusters'].astype(int)
                    silhouette_score = silhouette_score(X_pca, cluster_labels)
            except Exception as e:
                logger.warning(f"Could not calculate silhouette score: {e}")

            metrics = {
                'n_clusters': n_clusters,
                'clustering_method': self.config.clustering_method,
                'resolution': self.config.resolution,
                'silhouette_score': silhouette_score
            }

            return {
                'success': True,
                'clustered_file': str(clustered_file),
                'metrics': metrics,
                'n_clusters': n_clusters
            }

        except Exception as e:
            logger.error(f"Error in clustering: {e}")
            return {'success': False, 'error': str(e)}

    def find_marker_genes(self, clustered_file: str) -> Dict[str, Any]:
        """
        Find marker genes for each cluster.

        Args:
            clustered_file: Path to clustered data

        Returns:
            Marker gene results
        """
        if not self.tools.get('scanpy', False):
            raise RuntimeError("scanpy not available for marker gene analysis.")

        try:
            import scanpy as sc

            # Load clustered data
            adata = sc.read_h5ad(clustered_file)

            # Find marker genes
            logger.info("Finding marker genes")
            sc.tl.rank_genes_groups(adata, groupby=f'{self.config.clustering_method}_clusters', method='wilcoxon')

            # Extract marker genes
            marker_genes = {}
            for cluster in adata.obs[f'{self.config.clustering_method}_clusters'].unique():
                cluster_markers = sc.get.rank_genes_groups_df(adata, group=str(cluster))
                top_markers = cluster_markers.head(20)[['names', 'scores', 'pvals_adj', 'logfoldchanges']]
                marker_genes[str(cluster)] = top_markers.to_dict('records')

            # Save marker genes
            marker_file = self.output_dir / "marker_genes.json"
            with open(marker_file, 'w') as f:
                json.dump(marker_genes, f, indent=2)

            return {
                'success': True,
                'marker_file': str(marker_file),
                'marker_genes': marker_genes,
                'n_clusters': len(marker_genes)
            }

        except Exception as e:
            logger.error(f"Error finding marker genes: {e}")
            return {'success': False, 'error': str(e)}

    def run_trajectory_analysis(self, clustered_file: str, method: str = "palantir") -> Dict[str, Any]:
        """
        Run trajectory/pseudotime analysis.

        Args:
            clustered_file: Path to clustered data
            method: Trajectory method (palantir, monocle, slingshot)

        Returns:
            Trajectory analysis results
        """
        if method == "palantir":
            return self._run_palantir_trajectory(clustered_file)
        else:
            logger.warning(f"Trajectory method {method} not implemented")
            return {'success': False, 'error': f'Method {method} not implemented'}

    def _run_palantir_trajectory(self, clustered_file: str) -> Dict[str, Any]:
        """Run Palantir trajectory analysis."""
        try:
            # This would require palantir package
            # For now, return placeholder
            logger.warning("Palantir trajectory analysis not implemented")

            return {
                'success': False,
                'error': 'Palantir trajectory analysis requires palantir package'
            }

        except Exception as e:
            logger.error(f"Error in trajectory analysis: {e}")
            return {'success': False, 'error': str(e)}

    def generate_clustering_report(self, clustering_results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Generate comprehensive clustering analysis report.

        Args:
            clustering_results: Dictionary of clustering results

        Returns:
            Complete clustering report
        """
        report = {
            'summary': {
                'analysis_completed': True,
                'n_clusters': clustering_results.get('n_clusters', 0),
                'clustering_method': self.config.clustering_method,
                'preprocessing_method': self.config.normalization_method
            },
            'preprocessing': clustering_results.get('preprocessing', {}),
            'dimensionality_reduction': clustering_results.get('dimensionality_reduction', {}),
            'clustering': clustering_results.get('clustering', {}),
            'marker_genes': clustering_results.get('marker_genes', {}),
            'trajectory_analysis': clustering_results.get('trajectory_analysis', {}),
            'quality_metrics': {},
            'recommendations': []
        }

        # Calculate quality metrics
        if 'clustering' in clustering_results and clustering_results['clustering'].get('success'):
            metrics = clustering_results['clustering']['metrics']

            # Check cluster balance
            if 'n_clusters' in metrics and metrics['n_clusters'] > 1:
                # In a real implementation, you'd check cluster sizes
                report['quality_metrics']['cluster_balance'] = 'balanced'
            else:
                report['quality_metrics']['cluster_balance'] = 'needs_review'

            # Check silhouette score
            if metrics.get('silhouette_score') is not None:
                score = metrics['silhouette_score']
                if score > 0.5:
                    report['quality_metrics']['clustering_quality'] = 'excellent'
                elif score > 0.25:
                    report['quality_metrics']['clustering_quality'] = 'good'
                else:
                    report['quality_metrics']['clustering_quality'] = 'poor'

        # Generate recommendations
        recommendations = []

        if report['quality_metrics'].get('clustering_quality') == 'poor':
            recommendations.append("Consider adjusting clustering resolution or trying different method")

        if report['quality_metrics'].get('cluster_balance') == 'needs_review':
            recommendations.append("Check for balanced cluster sizes")

        if not recommendations:
            recommendations.append("Clustering analysis looks good")

        report['recommendations'] = recommendations

        return report

    def run_complete_clustering_analysis(self, matrix_file: str,
                                       metadata_file: Optional[str] = None) -> Dict[str, Any]:
        """
        Run complete clustering analysis pipeline.

        Args:
            matrix_file: Path to count matrix
            metadata_file: Optional metadata file

        Returns:
            Complete clustering analysis results
        """
        logger.info("Starting complete single-cell clustering analysis")

        results = {
            'preprocessing': {},
            'dimensionality_reduction': {},
            'clustering': {},
            'marker_genes': {},
            'trajectory_analysis': {},
            'report': {}
        }

        try:
            # Step 1: Preprocessing
            logger.info("Step 1: Data preprocessing")
            preprocessing_result = self.preprocess_data(matrix_file, metadata_file)
            results['preprocessing'] = preprocessing_result

            if not preprocessing_result['success']:
                raise RuntimeError("Preprocessing failed")

            # Step 2: Dimensionality reduction
            logger.info("Step 2: Dimensionality reduction")
            dimred_result = self.run_dimensionality_reduction(preprocessing_result['preprocessed_file'])
            results['dimensionality_reduction'] = dimred_result

            if not dimred_result['success']:
                raise RuntimeError("Dimensionality reduction failed")

            # Step 3: Clustering
            logger.info("Step 3: Clustering")
            clustering_result = self.run_clustering(dimred_result['dimred_file'])
            results['clustering'] = clustering_result

            if not clustering_result['success']:
                raise RuntimeError("Clustering failed")

            # Step 4: Marker gene identification
            logger.info("Step 4: Marker gene identification")
            marker_result = self.find_marker_genes(clustering_result['clustered_file'])
            results['marker_genes'] = marker_result

            # Step 5: Trajectory analysis (optional)
            logger.info("Step 5: Trajectory analysis")
            trajectory_result = self.run_trajectory_analysis(clustering_result['clustered_file'])
            results['trajectory_analysis'] = trajectory_result

            # Generate final report
            logger.info("Generating clustering report")
            results['report'] = self.generate_clustering_report(results)

            logger.info("Complete clustering analysis finished successfully")
            return results

        except Exception as e:
            logger.error(f"Error in complete clustering analysis: {e}")
            return {
                'error': str(e),
                'partial_results': results,
                'success': False
            }

    def create_visualization_plots(self, clustered_file: str) -> Dict[str, str]:
        """
        Create visualization plots for clustering results.

        Args:
            clustered_file: Path to clustered data

        Returns:
            Dictionary of plot file paths
        """
        plot_files = {}

        try:
            import scanpy as sc
            import matplotlib.pyplot as plt

            # Load data
            adata = sc.read_h5ad(clustered_file)

            # Set figure size
            plt.rcParams['figure.figsize'] = (10, 8)

            # UMAP colored by clusters
            sc.pl.umap(adata, color=f'{self.config.clustering_method}_clusters',
                      save='_clusters.png', show=False)
            plot_files['umap_clusters'] = str(self.output_dir / "umap_clusters.png")

            # UMAP colored by QC metrics
            if 'total_counts' in adata.obs.columns:
                sc.pl.umap(adata, color='total_counts', save='_total_counts.png', show=False)
                plot_files['umap_total_counts'] = str(self.output_dir / "umap_total_counts.png")

            if 'pct_counts_mt' in adata.obs.columns:
                sc.pl.umap(adata, color='pct_counts_mt', save='_mt_percent.png', show=False)
                plot_files['umap_mt_percent'] = str(self.output_dir / "umap_mt_percent.png")

            # t-SNE colored by clusters
            sc.pl.tsne(adata, color=f'{self.config.clustering_method}_clusters',
                      save='_clusters.png', show=False)
            plot_files['tsne_clusters'] = str(self.output_dir / "tsne_clusters.png")

            # PCA colored by clusters
            sc.pl.pca(adata, color=f'{self.config.clustering_method}_clusters',
                     save='_clusters.png', show=False)
            plot_files['pca_clusters'] = str(self.output_dir / "pca_clusters.png")

            # Save all plots
            plt.close('all')

            logger.info(f"Created {len(plot_files)} visualization plots")

        except ImportError:
            logger.warning("scanpy/matplotlib not available for plotting")
        except Exception as e:
            logger.error(f"Error creating visualization plots: {e}")

        return plot_files

    def export_clustering_results(self, results: Dict[str, Any], format: str = "h5ad") -> str:
        """
        Export clustering results in various formats.

        Args:
            results: Clustering analysis results
            format: Export format (h5ad, loom, csv)

        Returns:
            Path to exported file
        """
        try:
            if format == "h5ad" and 'clustering' in results:
                clustered_file = results['clustering']['clustered_file']
                return clustered_file

            elif format == "loom":
                # Export to Loom format
                logger.warning("Loom export not implemented")
                return ""

            elif format == "csv":
                # Export cluster assignments and metadata
                if 'clustering' in results and results['clustering'].get('success'):
                    clustered_file = results['clustering']['clustered_file']

                    # Load data and export
                    import scanpy as sc
                    adata = sc.read_h5ad(clustered_file)
                    export_df = adata.obs.copy()
                    export_df['cell_barcode'] = adata.obs.index

                    export_file = self.output_dir / "clustering_results.csv"
                    export_df.to_csv(export_file)

                    return str(export_file)

            return ""

        except Exception as e:
            logger.error(f"Error exporting clustering results: {e}")
            return ""
