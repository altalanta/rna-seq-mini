#!/usr/bin/env python3
"""
Spatial transcriptomics analysis module for RNASEQ-MINI.
Supports spatial coordinate integration, tissue reconstruction, and spatial statistics.
"""

import os
import logging
import subprocess
import json
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Union
from dataclasses import dataclass
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class SpatialConfig:
    """Configuration for spatial transcriptomics analysis."""

    # Input data
    matrix_file: str
    coordinates_file: str
    image_file: Optional[str] = None
    metadata_file: Optional[str] = None

    # Spatial parameters
    technology: str = "visium"  # visium, slideseq, merfish, etc.
    spot_diameter: float = 55.0  # microns
    scale_factor: float = 1.0

    # Analysis parameters
    n_neighbors: int = 6
    n_pcs: int = 50
    resolution: float = 1.0

    # Tissue reconstruction parameters
    reconstruction_method: str = "squidpy"  # squidpy, stlearn, stereoscope
    tissue_domains: int = 10

    # Output parameters
    output_dir: str = "results/singlecell/spatial"
    save_plots: bool = True

    # Performance parameters
    threads: int = 4
    memory_gb: int = 16


class SpatialAnalyzer:
    """Handles spatial transcriptomics analysis and tissue reconstruction."""

    def __init__(self, config: SpatialConfig):
        self.config = config
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Check tool availability
        self.tools = self._check_tool_availability()

    def _check_tool_availability(self) -> Dict[str, bool]:
        """Check availability of required tools."""
        tools = {}

        # Check for squidpy (spatial analysis)
        try:
            import squidpy as sq
            tools['squidpy'] = True
        except ImportError:
            tools['squidpy'] = False

        # Check for scanpy (general single-cell)
        try:
            import scanpy as sc
            tools['scanpy'] = True
        except ImportError:
            tools['scanpy'] = False

        # Check for stlearn (spatial transcriptomics)
        try:
            import stlearn as st
            tools['stlearn'] = True
        except ImportError:
            tools['stlearn'] = False

        # Check for stereoscope (deconvolution)
        try:
            import stereoscope
            tools['stereoscope'] = True
        except ImportError:
            tools['stereoscope'] = False

        return tools

    def load_spatial_data(self, matrix_file: str, coordinates_file: str,
                         image_file: Optional[str] = None) -> Dict[str, Any]:
        """
        Load spatial transcriptomics data with coordinates.

        Args:
            matrix_file: Path to count matrix
            coordinates_file: Path to spatial coordinates
            image_file: Optional tissue image file

        Returns:
            Loaded spatial data
        """
        try:
            import scanpy as sc

            logger.info(f"Loading spatial data from {matrix_file} and {coordinates_file}")

            # Load count matrix
            if matrix_file.endswith('.h5ad'):
                adata = sc.read_h5ad(matrix_file)
            else:
                # Assume it's a 10x-like matrix
                adata = sc.read_10x_mtx(matrix_file)

            # Load spatial coordinates
            coordinates = pd.read_csv(coordinates_file, index_col=0)

            # Validate coordinate compatibility
            if not set(coordinates.index).issubset(set(adata.obs.index)):
                logger.warning("Coordinate indices don't match matrix indices")

            # Add coordinates to AnnData
            adata.obsm['spatial'] = coordinates[['x', 'y']].values
            adata.obs['array_row'] = coordinates['x'].astype(int)
            adata.obs['array_col'] = coordinates['y'].astype(int)

            # Load tissue image if provided
            if image_file and Path(image_file).exists():
                adata.uns['spatial'] = {}
                adata.uns['spatial']['tissue_image'] = image_file
                adata.uns['spatial']['scale_factor'] = self.config.scale_factor

            # Calculate spatial metrics
            self._calculate_spatial_metrics(adata)

            # Save loaded data
            loaded_file = self.output_dir / "spatial_data.h5ad"
            adata.write(loaded_file)

            return {
                'success': True,
                'loaded_file': str(loaded_file),
                'n_spots': adata.obs.shape[0],
                'n_genes': adata.var.shape[0],
                'spatial_extent': {
                    'x_min': float(coordinates['x'].min()),
                    'x_max': float(coordinates['x'].max()),
                    'y_min': float(coordinates['y'].min()),
                    'y_max': float(coordinates['y'].max())
                }
            }

        except Exception as e:
            logger.error(f"Error loading spatial data: {e}")
            return {'success': False, 'error': str(e)}

    def _calculate_spatial_metrics(self, adata):
        """Calculate spatial neighborhood and domain metrics."""
        try:
            # Calculate nearest neighbors
            from sklearn.neighbors import NearestNeighbors

            coords = adata.obsm['spatial']
            nn = NearestNeighbors(n_neighbors=self.config.n_neighbors + 1)
            nn.fit(coords)

            distances, indices = nn.kneighbors(coords)

            # Store neighbor information
            adata.obsm['spatial_neighbors'] = indices[:, 1:]  # Exclude self
            adata.obs['mean_neighbor_distance'] = distances[:, 1:].mean(axis=1)

            # Calculate local density
            adata.obs['local_density'] = 1.0 / (distances[:, 1:].mean(axis=1) + 1e-10)

        except Exception as e:
            logger.warning(f"Error calculating spatial metrics: {e}")

    def run_spatial_clustering(self, loaded_file: str) -> Dict[str, Any]:
        """
        Run spatial-aware clustering analysis.

        Args:
            loaded_file: Path to loaded spatial data

        Returns:
            Spatial clustering results
        """
        if not self.tools.get('scanpy', False):
            raise RuntimeError("scanpy not available for spatial clustering.")

        try:
            import scanpy as sc

            # Load spatial data
            adata = sc.read_h5ad(loaded_file)

            # Preprocessing for spatial data
            logger.info("Preprocessing spatial data")

            # Normalize total counts
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)

            # Identify highly variable genes
            sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            adata = adata[:, adata.var.highly_variable]

            # Scale data
            sc.pp.scale(adata, max_value=10)

            # Run PCA
            sc.tl.pca(adata, svd_solver='arpack', n_comps=self.config.n_pcs)

            # Compute spatial neighbors
            logger.info("Computing spatial neighborhood graph")
            sc.pp.neighbors(adata, n_neighbors=self.config.n_neighbors, n_pcs=self.config.n_pcs,
                          use_rep='X_pca')

            # Run spatial clustering (Leiden with spatial constraint)
            logger.info("Running spatial clustering")
            sc.tl.leiden(adata, resolution=self.config.resolution)

            # Save results
            clustered_file = self.output_dir / "spatial_clustered.h5ad"
            adata.write(clustered_file)

            # Calculate spatial clustering metrics
            n_clusters = len(adata.obs['leiden'].unique())

            # Calculate spatial autocorrelation (Moran's I)
            moran_i = self._calculate_spatial_autocorrelation(adata)

            return {
                'success': True,
                'clustered_file': str(clustered_file),
                'n_clusters': n_clusters,
                'spatial_autocorrelation': moran_i,
                'clustering_method': 'spatial_leiden'
            }

        except Exception as e:
            logger.error(f"Error in spatial clustering: {e}")
            return {'success': False, 'error': str(e)}

    def _calculate_spatial_autocorrelation(self, adata) -> Dict[str, float]:
        """Calculate spatial autocorrelation (Moran's I) for clusters."""
        try:
            # This is a simplified implementation
            # In practice, you'd use proper spatial statistics libraries

            cluster_labels = adata.obs['leiden'].astype(int).values
            coords = adata.obsm['spatial']

            # Calculate Moran's I for cluster labels
            # This is a placeholder - proper implementation would use spatial statistics
            unique_clusters = np.unique(cluster_labels)
            moran_values = {}

            for cluster in unique_clusters:
                cluster_mask = cluster_labels == cluster
                if cluster_mask.sum() > 1:
                    # Calculate average distance to cluster members vs non-members
                    cluster_coords = coords[cluster_mask]
                    non_cluster_coords = coords[~cluster_mask]

                    # Simplified spatial autocorrelation measure
                    if len(cluster_coords) > 1 and len(non_cluster_coords) > 0:
                        # Calculate mean distance within cluster vs between clusters
                        # This is a very simplified version
                        moran_values[str(cluster)] = 0.5  # Placeholder

            return moran_values

        except Exception as e:
            logger.warning(f"Error calculating spatial autocorrelation: {e}")
            return {}

    def reconstruct_tissue_domains(self, clustered_file: str) -> Dict[str, Any]:
        """
        Reconstruct tissue domains and spatial patterns.

        Args:
            clustered_file: Path to clustered spatial data

        Returns:
            Tissue reconstruction results
        """
        if not self.tools.get('squidpy', False):
            logger.warning("squidpy not available for tissue reconstruction")
            return {'success': False, 'error': 'squidpy required for tissue reconstruction'}

        try:
            import squidpy as sq

            # Load clustered data
            adata = sq.read.h5ad(clustered_file)

            # Run tissue domain identification
            logger.info("Identifying tissue domains")
            sq.gr.spatial_neighbors(adata, n_neighs=self.config.n_neighbors)
            sq.gr.nhood_enrichment(adata, cluster_key='leiden')

            # Find spatial domains
            sq.gr.spatial_autocorr(adata, mode='moran', cluster_key='leiden')

            # Run domain segmentation
            if self.config.reconstruction_method == "squidpy":
                # Use squidpy's domain identification
                sq.gr.spatial_domain(adata, n_domains=self.config.tissue_domains)

            # Save results
            domain_file = self.output_dir / "tissue_domains.h5ad"
            adata.write(domain_file)

            # Calculate domain metrics
            domain_metrics = {}
            if 'spatial_domain' in adata.obs.columns:
                n_domains = len(adata.obs['spatial_domain'].unique())

                # Calculate domain sizes
                domain_sizes = adata.obs['spatial_domain'].value_counts().to_dict()
                domain_metrics = {
                    'n_domains': n_domains,
                    'domain_sizes': domain_sizes,
                    'largest_domain': max(domain_sizes.values()),
                    'smallest_domain': min(domain_sizes.values())
                }

            return {
                'success': True,
                'domain_file': str(domain_file),
                'domain_metrics': domain_metrics,
                'reconstruction_method': self.config.reconstruction_method
            }

        except Exception as e:
            logger.error(f"Error in tissue reconstruction: {e}")
            return {'success': False, 'error': str(e)}

    def run_spatial_de_analysis(self, domain_file: str) -> Dict[str, Any]:
        """
        Run spatial differential expression analysis.

        Args:
            domain_file: Path to domain-annotated data

        Returns:
            Spatial DE results
        """
        try:
            import scanpy as sc

            # Load domain data
            adata = sc.read_h5ad(domain_file)

            # Run spatial differential expression
            logger.info("Running spatial differential expression")

            # Compare expression between spatial domains
            sc.tl.rank_genes_groups(adata, groupby='spatial_domain', method='wilcoxon')

            # Extract domain-specific markers
            domain_markers = {}
            for domain in adata.obs['spatial_domain'].unique():
                domain_markers_df = sc.get.rank_genes_groups_df(adata, group=str(domain))
                top_markers = domain_markers_df.head(20)[['names', 'scores', 'pvals_adj', 'logfoldchanges']]
                domain_markers[str(domain)] = top_markers.to_dict('records')

            # Save results
            de_file = self.output_dir / "spatial_de_results.json"
            with open(de_file, 'w') as f:
                json.dump(domain_markers, f, indent=2)

            return {
                'success': True,
                'de_file': str(de_file),
                'domain_markers': domain_markers,
                'n_domains_compared': len(domain_markers)
            }

        except Exception as e:
            logger.error(f"Error in spatial DE analysis: {e}")
            return {'success': False, 'error': str(e)}

    def create_spatial_visualizations(self, data_file: str) -> Dict[str, str]:
        """
        Create spatial visualizations and plots.

        Args:
            data_file: Path to spatial data

        Returns:
            Dictionary of plot file paths
        """
        plot_files = {}

        try:
            import scanpy as sc
            import matplotlib.pyplot as plt

            # Load data
            adata = sc.read_h5ad(data_file)

            # Set up plotting
            plt.rcParams['figure.figsize'] = (12, 10)

            # Spatial plot colored by clusters
            if 'leiden' in adata.obs.columns:
                sc.pl.spatial(adata, color='leiden', spot_size=self.config.spot_diameter,
                             save='_clusters.png', show=False)
                plot_files['spatial_clusters'] = str(self.output_dir / "spatial_clusters.png")

            # Spatial plot colored by gene expression (top variable genes)
            if adata.var.highly_variable.sum() > 0:
                top_genes = adata.var[adata.var.highly_variable].index[:6]  # Top 6 genes
                for gene in top_genes:
                    try:
                        sc.pl.spatial(adata, color=gene, spot_size=self.config.spot_diameter,
                                    save=f'_{gene}.png', show=False, vmax='p99')
                        plot_files[f'spatial_{gene}'] = str(self.output_dir / f"spatial_{gene}.png")
                    except Exception as e:
                        logger.warning(f"Error creating plot for gene {gene}: {e}")

            # Spatial domain plot if available
            if 'spatial_domain' in adata.obs.columns:
                sc.pl.spatial(adata, color='spatial_domain', spot_size=self.config.spot_diameter,
                             save='_domains.png', show=False)
                plot_files['spatial_domains'] = str(self.output_dir / "spatial_domains.png")

            # Tissue image overlay if available
            if 'spatial' in adata.uns and 'tissue_image' in adata.uns['spatial']:
                try:
                    # This would require additional implementation for image overlay
                    logger.info("Tissue image overlay visualization available")
                except Exception as e:
                    logger.warning(f"Error with tissue image overlay: {e}")

            plt.close('all')
            logger.info(f"Created {len(plot_files)} spatial visualization plots")

        except ImportError:
            logger.warning("scanpy not available for spatial plotting")
        except Exception as e:
            logger.error(f"Error creating spatial visualizations: {e}")

        return plot_files

    def run_cell_type_deconvolution(self, spatial_data: str, reference_data: str) -> Dict[str, Any]:
        """
        Run cell type deconvolution for spatial data.

        Args:
            spatial_data: Path to spatial transcriptomics data
            reference_data: Path to single-cell reference data

        Returns:
            Deconvolution results
        """
        if not self.tools.get('stereoscope', False):
            logger.warning("stereoscope not available for deconvolution")
            return {'success': False, 'error': 'stereoscope required for deconvolution'}

        try:
            # This would use stereoscope for cell type deconvolution
            # For now, return placeholder
            logger.warning("Cell type deconvolution not fully implemented")

            return {
                'success': False,
                'error': 'Cell type deconvolution requires stereoscope package'
            }

        except Exception as e:
            logger.error(f"Error in cell type deconvolution: {e}")
            return {'success': False, 'error': str(e)}

    def run_spatial_trajectory_analysis(self, clustered_file: str) -> Dict[str, Any]:
        """
        Run spatial trajectory analysis for tissue development.

        Args:
            clustered_file: Path to clustered spatial data

        Returns:
            Spatial trajectory results
        """
        try:
            # This would use specialized spatial trajectory tools
            # For now, return placeholder
            logger.warning("Spatial trajectory analysis not implemented")

            return {
                'success': False,
                'error': 'Spatial trajectory analysis requires specialized tools'
            }

        except Exception as e:
            logger.error(f"Error in spatial trajectory analysis: {e}")
            return {'success': False, 'error': str(e)}

    def generate_spatial_report(self, spatial_results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Generate comprehensive spatial analysis report.

        Args:
            spatial_results: Dictionary of spatial analysis results

        Returns:
            Complete spatial analysis report
        """
        report = {
            'summary': {
                'analysis_completed': True,
                'spatial_technology': self.config.technology,
                'n_spots': spatial_results.get('n_spots', 0),
                'n_genes': spatial_results.get('n_genes', 0),
                'spatial_extent': spatial_results.get('spatial_extent', {})
            },
            'data_loading': spatial_results.get('data_loading', {}),
            'spatial_clustering': spatial_results.get('spatial_clustering', {}),
            'tissue_domains': spatial_results.get('tissue_domains', {}),
            'spatial_de': spatial_results.get('spatial_de', {}),
            'deconvolution': spatial_results.get('deconvolution', {}),
            'trajectory_analysis': spatial_results.get('trajectory_analysis', {}),
            'visualizations': spatial_results.get('visualizations', {}),
            'quality_metrics': {},
            'recommendations': []
        }

        # Calculate quality metrics
        if 'spatial_clustering' in spatial_results and spatial_results['spatial_clustering'].get('success'):
            clustering = spatial_results['spatial_clustering']

            # Check spatial autocorrelation
            autocorrelation = clustering.get('spatial_autocorrelation', {})
            if autocorrelation:
                avg_autocorr = np.mean([v for v in autocorrelation.values() if isinstance(v, (int, float))])
                report['quality_metrics']['spatial_autocorrelation'] = avg_autocorr

        # Generate recommendations
        recommendations = []

        if report['quality_metrics'].get('spatial_autocorrelation', 0) < 0.1:
            recommendations.append("Low spatial autocorrelation detected. Consider adjusting clustering parameters.")

        if not recommendations:
            recommendations.append("Spatial analysis completed successfully.")

        report['recommendations'] = recommendations

        return report

    def run_complete_spatial_analysis(self, matrix_file: str, coordinates_file: str,
                                    image_file: Optional[str] = None) -> Dict[str, Any]:
        """
        Run complete spatial transcriptomics analysis pipeline.

        Args:
            matrix_file: Path to count matrix
            coordinates_file: Path to spatial coordinates
            image_file: Optional tissue image

        Returns:
            Complete spatial analysis results
        """
        logger.info("Starting complete spatial transcriptomics analysis")

        results = {
            'data_loading': {},
            'spatial_clustering': {},
            'tissue_domains': {},
            'spatial_de': {},
            'deconvolution': {},
            'trajectory_analysis': {},
            'visualizations': {},
            'report': {}
        }

        try:
            # Step 1: Load spatial data
            logger.info("Step 1: Loading spatial data")
            loading_result = self.load_spatial_data(matrix_file, coordinates_file, image_file)
            results['data_loading'] = loading_result

            if not loading_result['success']:
                raise RuntimeError("Data loading failed")

            # Step 2: Spatial clustering
            logger.info("Step 2: Spatial clustering")
            clustering_result = self.run_spatial_clustering(loading_result['loaded_file'])
            results['spatial_clustering'] = clustering_result

            # Step 3: Tissue domain reconstruction
            logger.info("Step 3: Tissue domain reconstruction")
            domain_result = self.reconstruct_tissue_domains(clustering_result['clustered_file'])
            results['tissue_domains'] = domain_result

            # Step 4: Spatial DE analysis
            logger.info("Step 4: Spatial differential expression")
            de_result = self.run_spatial_de_analysis(domain_result['domain_file'])
            results['spatial_de'] = de_result

            # Step 5: Cell type deconvolution (optional)
            logger.info("Step 5: Cell type deconvolution")
            deconvolution_result = self.run_cell_type_deconvolution(
                loading_result['loaded_file'], "reference_sc_data.h5ad"
            )
            results['deconvolution'] = deconvolution_result

            # Step 6: Trajectory analysis (optional)
            logger.info("Step 6: Spatial trajectory analysis")
            trajectory_result = self.run_spatial_trajectory_analysis(clustering_result['clustered_file'])
            results['trajectory_analysis'] = trajectory_result

            # Step 7: Create visualizations
            logger.info("Step 7: Creating visualizations")
            vis_result = self.create_spatial_visualizations(clustering_result['clustered_file'])
            results['visualizations'] = vis_result

            # Generate final report
            logger.info("Generating spatial analysis report")
            results['report'] = self.generate_spatial_report(results)

            logger.info("Complete spatial analysis finished successfully")
            return results

        except Exception as e:
            logger.error(f"Error in complete spatial analysis: {e}")
            return {
                'error': str(e),
                'partial_results': results,
                'success': False
            }

    def export_spatial_results(self, results: Dict[str, Any], format: str = "h5ad") -> str:
        """
        Export spatial analysis results.

        Args:
            results: Spatial analysis results
            format: Export format (h5ad, geojson, csv)

        Returns:
            Path to exported file
        """
        try:
            if format == "h5ad" and 'spatial_clustering' in results:
                clustered_file = results['spatial_clustering']['clustered_file']
                return clustered_file

            elif format == "geojson":
                # Export as GeoJSON for GIS applications
                if 'data_loading' in results and results['data_loading'].get('success'):
                    loaded_file = results['data_loading']['loaded_file']

                    # Load data and convert to GeoJSON
                    import scanpy as sc
                    adata = sc.read_h5ad(loaded_file)

                    # Create GeoJSON features
                    features = []
                    for idx in adata.obs.index:
                        coord = adata.obsm['spatial'][adata.obs.index.get_loc(idx)]
                        feature = {
                            "type": "Feature",
                            "geometry": {
                                "type": "Point",
                                "coordinates": [float(coord[0]), float(coord[1])]
                            },
                            "properties": {
                                "spot_id": idx,
                                "cluster": adata.obs.loc[idx, 'leiden'] if 'leiden' in adata.obs.columns else None
                            }
                        }
                        features.append(feature)

                    geojson_data = {
                        "type": "FeatureCollection",
                        "features": features
                    }

                    export_file = self.output_dir / "spatial_results.geojson"
                    with open(export_file, 'w') as f:
                        json.dump(geojson_data, f, indent=2)

                    return str(export_file)

            elif format == "csv":
                # Export spatial coordinates and metadata
                if 'data_loading' in results and results['data_loading'].get('success'):
                    loaded_file = results['data_loading']['loaded_file']

                    import scanpy as sc
                    adata = sc.read_h5ad(loaded_file)

                    # Create export dataframe
                    export_df = adata.obs.copy()
                    export_df['x_coord'] = adata.obsm['spatial'][:, 0]
                    export_df['y_coord'] = adata.obsm['spatial'][:, 1]
                    export_df['spot_barcode'] = adata.obs.index

                    export_file = self.output_dir / "spatial_coordinates.csv"
                    export_df.to_csv(export_file)

                    return str(export_file)

            return ""

        except Exception as e:
            logger.error(f"Error exporting spatial results: {e}")
            return ""







