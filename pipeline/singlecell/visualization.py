#!/usr/bin/env python3
"""
Single-cell visualization module for RNASEQ-MINI.
Provides comprehensive interactive visualizations for single-cell and spatial data.
"""

import os
import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Union
from dataclasses import dataclass
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class VisualizationConfig:
    """Configuration for single-cell visualizations."""

    # Input data
    data_file: str
    metadata_file: Optional[str] = None

    # Plot types
    plot_types: List[str] = None  # umap, tsne, pca, spatial, heatmap, etc.

    # Styling parameters
    color_palette: str = "default"
    figure_size: Tuple[int, int] = (10, 8)
    dpi: int = 300

    # Output parameters
    output_dir: str = "results/singlecell/visualizations"
    save_formats: List[str] = None  # png, pdf, svg, html
    show_plots: bool = False

    # Interactive parameters
    interactive: bool = True
    hover_info: bool = True

    # Performance parameters
    max_points: int = 100000
    subsample: bool = True


class SingleCellVisualizer:
    """Handles comprehensive single-cell data visualization."""

    def __init__(self, config: VisualizationConfig):
        self.config = config
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Set default plot types
        if config.plot_types is None:
            config.plot_types = ["umap", "tsne", "pca", "heatmap", "violin", "dotplot"]

        # Set default save formats
        if config.save_formats is None:
            config.save_formats = ["png", "html"] if config.interactive else ["png"]

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

        # Check for plotly
        try:
            import plotly
            tools['plotly'] = True
        except ImportError:
            tools['plotly'] = False

        # Check for matplotlib
        try:
            import matplotlib
            tools['matplotlib'] = True
        except ImportError:
            tools['matplotlib'] = False

        # Check for seaborn
        try:
            import seaborn as sns
            tools['seaborn'] = True
        except ImportError:
            tools['seaborn'] = False

        return tools

    def load_data(self, data_file: str, metadata_file: Optional[str] = None) -> Dict[str, Any]:
        """
        Load single-cell data for visualization.

        Args:
            data_file: Path to single-cell data file
            metadata_file: Optional metadata file

        Returns:
            Loaded data information
        """
        try:
            import scanpy as sc

            logger.info(f"Loading data from {data_file}")

            # Load data
            if data_file.endswith('.h5ad'):
                adata = sc.read_h5ad(data_file)
            else:
                # Assume it's a 10x-like matrix
                adata = sc.read_10x_mtx(data_file)

            # Load metadata if provided
            if metadata_file:
                metadata = pd.read_csv(metadata_file, index_col=0)
                adata.obs = adata.obs.join(metadata, how='left')

            return {
                'success': True,
                'adata': adata,
                'n_cells': adata.obs.shape[0],
                'n_genes': adata.var.shape[0],
                'available_embeddings': list(adata.obsm.keys()),
                'available_clusters': [col for col in adata.obs.columns if 'cluster' in col.lower() or 'leiden' in col.lower() or 'louvain' in col.lower()]
            }

        except Exception as e:
            logger.error(f"Error loading data: {e}")
            return {'success': False, 'error': str(e)}

    def create_umap_visualizations(self, adata) -> Dict[str, str]:
        """
        Create UMAP visualizations.

        Args:
            adata: AnnData object

        Returns:
            Dictionary of plot file paths
        """
        plot_files = {}

        try:
            import scanpy as sc
            import matplotlib.pyplot as plt

            # Set figure size
            plt.rcParams['figure.figsize'] = self.config.figure_size

            # UMAP colored by clusters
            cluster_cols = [col for col in adata.obs.columns if 'cluster' in col.lower() or 'leiden' in col.lower() or 'louvain' in col.lower()]
            for col in cluster_cols[:3]:  # Limit to first 3 cluster columns
                try:
                    sc.pl.umap(adata, color=col, save=f'_umap_{col}.png', show=self.config.show_plots)
                    plot_files[f'umap_{col}'] = str(self.output_dir / f"umap_{col}.png")
                except Exception as e:
                    logger.warning(f"Error creating UMAP for {col}: {e}")

            # UMAP colored by QC metrics
            qc_cols = ['total_counts', 'n_genes', 'pct_counts_mt', 'pct_counts_ribo']
            for col in qc_cols:
                if col in adata.obs.columns:
                    try:
                        sc.pl.umap(adata, color=col, save=f'_umap_{col}.png', show=self.config.show_plots)
                        plot_files[f'umap_{col}'] = str(self.output_dir / f"umap_{col}.png")
                    except Exception as e:
                        logger.warning(f"Error creating UMAP for {col}: {e}")

            # Interactive UMAP with plotly
            if self.tools.get('plotly', False) and 'X_umap' in adata.obsm:
                try:
                    import plotly.express as px
                    import plotly.graph_objects as go

                    # Create interactive UMAP
                    if cluster_cols:
                        fig = px.scatter(
                            x=adata.obsm['X_umap'][:, 0],
                            y=adata.obsm['X_umap'][:, 1],
                            color=adata.obs[cluster_cols[0]],
                            labels={'color': cluster_cols[0]},
                            title=f'Interactive UMAP - {cluster_cols[0]}'
                        )

                        # Save interactive plot
                        interactive_file = self.output_dir / "interactive_umap.html"
                        fig.write_html(str(interactive_file))
                        plot_files['interactive_umap'] = str(interactive_file)

                except Exception as e:
                    logger.warning(f"Error creating interactive UMAP: {e}")

            plt.close('all')

        except ImportError:
            logger.warning("scanpy not available for UMAP visualization")
        except Exception as e:
            logger.error(f"Error creating UMAP visualizations: {e}")

        return plot_files

    def create_tsne_visualizations(self, adata) -> Dict[str, str]:
        """
        Create t-SNE visualizations.

        Args:
            adata: AnnData object

        Returns:
            Dictionary of plot file paths
        """
        plot_files = {}

        try:
            import scanpy as sc
            import matplotlib.pyplot as plt

            plt.rcParams['figure.figsize'] = self.config.figure_size

            # t-SNE colored by clusters
            cluster_cols = [col for col in adata.obs.columns if 'cluster' in col.lower() or 'leiden' in col.lower() or 'louvain' in col.lower()]
            for col in cluster_cols[:3]:
                try:
                    sc.pl.tsne(adata, color=col, save=f'_tsne_{col}.png', show=self.config.show_plots)
                    plot_files[f'tsne_{col}'] = str(self.output_dir / f"tsne_{col}.png")
                except Exception as e:
                    logger.warning(f"Error creating t-SNE for {col}: {e}")

            # t-SNE colored by QC metrics
            qc_cols = ['total_counts', 'n_genes', 'pct_counts_mt']
            for col in qc_cols:
                if col in adata.obs.columns:
                    try:
                        sc.pl.tsne(adata, color=col, save=f'_tsne_{col}.png', show=self.config.show_plots)
                        plot_files[f'tsne_{col}'] = str(self.output_dir / f"tsne_{col}.png")
                    except Exception as e:
                        logger.warning(f"Error creating t-SNE for {col}: {e}")

            plt.close('all')

        except ImportError:
            logger.warning("scanpy not available for t-SNE visualization")
        except Exception as e:
            logger.error(f"Error creating t-SNE visualizations: {e}")

        return plot_files

    def create_pca_visualizations(self, adata) -> Dict[str, str]:
        """
        Create PCA visualizations.

        Args:
            adata: AnnData object

        Returns:
            Dictionary of plot file paths
        """
        plot_files = {}

        try:
            import scanpy as sc
            import matplotlib.pyplot as plt

            plt.rcParams['figure.figsize'] = self.config.figure_size

            # PCA colored by clusters
            cluster_cols = [col for col in adata.obs.columns if 'cluster' in col.lower() or 'leiden' in col.lower() or 'louvain' in col.lower()]
            for col in cluster_cols[:3]:
                try:
                    sc.pl.pca(adata, color=col, save=f'_pca_{col}.png', show=self.config.show_plots)
                    plot_files[f'pca_{col}'] = str(self.output_dir / f"pca_{col}.png")
                except Exception as e:
                    logger.warning(f"Error creating PCA for {col}: {e}")

            # PCA colored by QC metrics
            qc_cols = ['total_counts', 'n_genes', 'pct_counts_mt']
            for col in qc_cols:
                if col in adata.obs.columns:
                    try:
                        sc.pl.pca(adata, color=col, save=f'_pca_{col}.png', show=self.config.show_plots)
                        plot_files[f'pca_{col}'] = str(self.output_dir / f"pca_{col}.png")
                    except Exception as e:
                        logger.warning(f"Error creating PCA for {col}: {e}")

            plt.close('all')

        except ImportError:
            logger.warning("scanpy not available for PCA visualization")
        except Exception as e:
            logger.error(f"Error creating PCA visualizations: {e}")

        return plot_files

    def create_heatmap_visualizations(self, adata) -> Dict[str, str]:
        """
        Create heatmap visualizations.

        Args:
            adata: AnnData object

        Returns:
            Dictionary of plot file paths
        """
        plot_files = {}

        try:
            import scanpy as sc
            import matplotlib.pyplot as plt
            import seaborn as sns

            plt.rcParams['figure.figsize'] = self.config.figure_size

            # Get top variable genes
            if 'highly_variable' in adata.var.columns:
                top_genes = adata.var[adata.var.highly_variable].index[:50]
            else:
                # Use top expressed genes
                gene_means = adata.X.mean(axis=0).A1
                top_indices = np.argsort(gene_means)[-50:]
                top_genes = adata.var_names[top_indices]

            # Heatmap of top genes
            try:
                sc.pl.heatmap(adata, var_names=top_genes, save='_top_genes.png', show=self.config.show_plots)
                plot_files['heatmap_top_genes'] = str(self.output_dir / "heatmap_top_genes.png")
            except Exception as e:
                logger.warning(f"Error creating gene heatmap: {e}")

            # Heatmap colored by clusters
            cluster_cols = [col for col in adata.obs.columns if 'cluster' in col.lower() or 'leiden' in col.lower() or 'louvain' in col.lower()]
            if cluster_cols:
                try:
                    sc.pl.heatmap(adata, var_names=top_genes, groupby=cluster_cols[0],
                                save=f'_by_{cluster_cols[0]}.png', show=self.config.show_plots)
                    plot_files[f'heatmap_by_{cluster_cols[0]}'] = str(self.output_dir / f"heatmap_by_{cluster_cols[0]}.png")
                except Exception as e:
                    logger.warning(f"Error creating cluster heatmap: {e}")

            plt.close('all')

        except ImportError:
            logger.warning("scanpy/seaborn not available for heatmap visualization")
        except Exception as e:
            logger.error(f"Error creating heatmap visualizations: {e}")

        return plot_files

    def create_violin_plots(self, adata) -> Dict[str, str]:
        """
        Create violin plots for gene expression.

        Args:
            adata: AnnData object

        Returns:
            Dictionary of plot file paths
        """
        plot_files = {}

        try:
            import scanpy as sc
            import matplotlib.pyplot as plt

            plt.rcParams['figure.figsize'] = self.config.figure_size

            # Violin plots for marker genes
            marker_genes = self._get_marker_genes(adata)
            if marker_genes:
                try:
                    sc.pl.violin(adata, keys=marker_genes[:10], groupby='leiden',
                               save='_marker_genes.png', show=self.config.show_plots)
                    plot_files['violin_marker_genes'] = str(self.output_dir / "violin_marker_genes.png")
                except Exception as e:
                    logger.warning(f"Error creating marker gene violin plot: {e}")

            # Violin plots for QC metrics
            qc_cols = ['total_counts', 'n_genes', 'pct_counts_mt']
            existing_qc = [col for col in qc_cols if col in adata.obs.columns]
            if existing_qc:
                try:
                    sc.pl.violin(adata, keys=existing_qc, groupby='leiden',
                               save='_qc_metrics.png', show=self.config.show_plots)
                    plot_files['violin_qc_metrics'] = str(self.output_dir / "violin_qc_metrics.png")
                except Exception as e:
                    logger.warning(f"Error creating QC violin plot: {e}")

            plt.close('all')

        except ImportError:
            logger.warning("scanpy not available for violin plots")
        except Exception as e:
            logger.error(f"Error creating violin plots: {e}")

        return plot_files

    def create_dot_plots(self, adata) -> Dict[str, str]:
        """
        Create dot plots for gene expression.

        Args:
            adata: AnnData object

        Returns:
            Dictionary of plot file paths
        """
        plot_files = {}

        try:
            import scanpy as sc
            import matplotlib.pyplot as plt

            plt.rcParams['figure.figsize'] = self.config.figure_size

            # Dot plot for marker genes
            marker_genes = self._get_marker_genes(adata)
            if marker_genes:
                try:
                    sc.pl.dotplot(adata, var_names=marker_genes[:20], groupby='leiden',
                                save='_marker_genes.png', show=self.config.show_plots)
                    plot_files['dotplot_marker_genes'] = str(self.output_dir / "dotplot_marker_genes.png")
                except Exception as e:
                    logger.warning(f"Error creating marker gene dot plot: {e}")

            plt.close('all')

        except ImportError:
            logger.warning("scanpy not available for dot plots")
        except Exception as e:
            logger.error(f"Error creating dot plots: {e}")

        return plot_files

    def create_spatial_visualizations(self, adata) -> Dict[str, str]:
        """
        Create spatial visualizations for spatial transcriptomics.

        Args:
            adata: AnnData object with spatial coordinates

        Returns:
            Dictionary of plot file paths
        """
        plot_files = {}

        try:
            import scanpy as sc
            import matplotlib.pyplot as plt

            plt.rcParams['figure.figsize'] = self.config.figure_size

            # Check if spatial coordinates are available
            if 'spatial' in adata.obsm:
                coords = adata.obsm['spatial']

                # Spatial plot colored by clusters
                cluster_cols = [col for col in adata.obs.columns if 'cluster' in col.lower() or 'leiden' in col.lower() or 'louvain' in col.lower()]
                for col in cluster_cols[:3]:
                    try:
                        sc.pl.spatial(adata, color=col, spot_size=50,
                                    save=f'_spatial_{col}.png', show=self.config.show_plots)
                        plot_files[f'spatial_{col}'] = str(self.output_dir / f"spatial_{col}.png")
                    except Exception as e:
                        logger.warning(f"Error creating spatial plot for {col}: {e}")

                # Spatial plot colored by gene expression
                marker_genes = self._get_marker_genes(adata)
                for gene in marker_genes[:6]:  # Top 6 marker genes
                    try:
                        sc.pl.spatial(adata, color=gene, spot_size=50,
                                    save=f'_spatial_{gene}.png', show=self.config.show_plots, vmax='p99')
                        plot_files[f'spatial_{gene}'] = str(self.output_dir / f"spatial_{gene}.png")
                    except Exception as e:
                        logger.warning(f"Error creating spatial plot for {gene}: {e}")

            plt.close('all')

        except ImportError:
            logger.warning("scanpy not available for spatial visualization")
        except Exception as e:
            logger.error(f"Error creating spatial visualizations: {e}")

        return plot_files

    def create_interactive_plots(self, adata) -> Dict[str, str]:
        """
        Create interactive Plotly visualizations.

        Args:
            adata: AnnData object

        Returns:
            Dictionary of plot file paths
        """
        plot_files = {}

        try:
            import plotly.express as px
            import plotly.graph_objects as go

            # Interactive UMAP
            if 'X_umap' in adata.obsm:
                cluster_cols = [col for col in adata.obs.columns if 'cluster' in col.lower() or 'leiden' in col.lower() or 'louvain' in col.lower()]

                if cluster_cols:
                    fig = px.scatter(
                        x=adata.obsm['X_umap'][:, 0],
                        y=adata.obsm['X_umap'][:, 1],
                        color=adata.obs[cluster_cols[0]],
                        hover_data=['total_counts', 'n_genes'] if 'total_counts' in adata.obs.columns else None,
                        labels={'color': cluster_cols[0]},
                        title=f'Interactive UMAP - {cluster_cols[0]}'
                    )

                    # Save interactive plot
                    interactive_file = self.output_dir / "interactive_umap.html"
                    fig.write_html(str(interactive_file))
                    plot_files['interactive_umap'] = str(interactive_file)

            # Interactive heatmap
            if 'highly_variable' in adata.var.columns and adata.var.highly_variable.sum() > 0:
                try:
                    # Get top variable genes
                    top_genes = adata.var[adata.var.highly_variable].index[:20]

                    # Subsample cells for performance
                    n_cells = min(1000, adata.obs.shape[0])
                    subsample_idx = np.random.choice(adata.obs.shape[0], n_cells, replace=False)
                    subsample_adata = adata[subsample_idx, top_genes]

                    # Create heatmap data
                    heatmap_data = pd.DataFrame(
                        subsample_adata.X.toarray(),
                        index=subsample_adata.obs.index,
                        columns=subsample_adata.var_names
                    )

                    # Create interactive heatmap
                    fig = go.Figure(data=go.Heatmap(
                        z=heatmap_data.values,
                        x=heatmap_data.columns,
                        y=heatmap_data.index,
                        colorscale='RdBu',
                        hoverongaps=False
                    ))

                    fig.update_layout(
                        title='Interactive Expression Heatmap',
                        xaxis_title='Genes',
                        yaxis_title='Cells'
                    )

                    interactive_heatmap = self.output_dir / "interactive_heatmap.html"
                    fig.write_html(str(interactive_heatmap))
                    plot_files['interactive_heatmap'] = str(interactive_heatmap)

                except Exception as e:
                    logger.warning(f"Error creating interactive heatmap: {e}")

        except ImportError:
            logger.warning("plotly not available for interactive plots")
        except Exception as e:
            logger.error(f"Error creating interactive plots: {e}")

        return plot_files

    def _get_marker_genes(self, adata) -> List[str]:
        """Get list of marker genes for visualization."""
        # Try to get marker genes from various sources
        marker_genes = []

        # From highly variable genes
        if 'highly_variable' in adata.var.columns:
            hv_genes = adata.var[adata.var.highly_variable].index.tolist()
            marker_genes.extend(hv_genes[:20])

        # From gene names (look for common marker patterns)
        gene_names = adata.var_names.tolist()
        marker_patterns = ['GAPDH', 'ACTB', 'PTPRC', 'CD3', 'CD19', 'CD14', 'CD8', 'CD4']
        for pattern in marker_patterns:
            matching = [gene for gene in gene_names if pattern in gene.upper()]
            marker_genes.extend(matching[:5])  # Limit to 5 per pattern

        return list(set(marker_genes))[:50]  # Return up to 50 unique genes

    def create_comprehensive_visualizations(self, data_file: str, metadata_file: Optional[str] = None) -> Dict[str, Any]:
        """
        Create comprehensive set of visualizations for single-cell data.

        Args:
            data_file: Path to single-cell data file
            metadata_file: Optional metadata file

        Returns:
            Dictionary of all visualization results
        """
        logger.info("Creating comprehensive single-cell visualizations")

        # Load data
        data_info = self.load_data(data_file, metadata_file)
        if not data_info['success']:
            return {'success': False, 'error': data_info['error']}

        adata = data_info['adata']

        # Create all types of visualizations
        results = {
            'umap_plots': self.create_umap_visualizations(adata),
            'tsne_plots': self.create_tsne_visualizations(adata),
            'pca_plots': self.create_pca_visualizations(adata),
            'heatmap_plots': self.create_heatmap_visualizations(adata),
            'violin_plots': self.create_violin_plots(adata),
            'dot_plots': self.create_dot_plots(adata),
            'spatial_plots': self.create_spatial_visualizations(adata),
            'interactive_plots': self.create_interactive_plots(adata),
            'summary': {
                'total_plots': 0,
                'plot_types': [],
                'data_info': data_info
            }
        }

        # Calculate summary statistics
        total_plots = sum(len(plots) for plots in results.values() if isinstance(plots, dict))
        plot_types = [key for key in results.keys() if key != 'summary' and results[key]]

        results['summary'].update({
            'total_plots': total_plots,
            'plot_types': plot_types
        })

        logger.info(f"Created {total_plots} visualizations across {len(plot_types)} plot types")
        return results

    def generate_visualization_report(self, visualization_results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Generate report of created visualizations.

        Args:
            visualization_results: Results from create_comprehensive_visualizations

        Returns:
            Visualization report
        """
        report = {
            'summary': visualization_results.get('summary', {}),
            'plot_inventory': {},
            'recommendations': []
        }

        # Inventory all created plots
        for plot_type, plots in visualization_results.items():
            if plot_type != 'summary' and isinstance(plots, dict):
                report['plot_inventory'][plot_type] = {
                    'count': len(plots),
                    'plots': list(plots.keys()),
                    'files': list(plots.values())
                }

        # Generate recommendations
        recommendations = []

        if report['summary'].get('total_plots', 0) == 0:
            recommendations.append("No visualizations were created. Check data format and required packages.")
        else:
            recommendations.append(f"Successfully created {report['summary']['total_plots']} visualizations.")

        if 'interactive_plots' not in report['plot_inventory'] or report['plot_inventory']['interactive_plots']['count'] == 0:
            recommendations.append("Consider installing plotly for interactive visualizations.")

        report['recommendations'] = recommendations

        return report

    def export_visualizations(self, visualization_results: Dict[str, Any], format: str = "html") -> str:
        """
        Export visualizations in specified format.

        Args:
            visualization_results: Visualization results
            format: Export format (html, pdf, png)

        Returns:
            Path to exported file
        """
        try:
            if format == "html":
                # Create HTML gallery
                html_content = self._create_html_gallery(visualization_results)
                export_file = self.output_dir / "visualization_gallery.html"

                with open(export_file, 'w') as f:
                    f.write(html_content)

                return str(export_file)

            elif format == "pdf":
                # Combine plots into PDF
                logger.warning("PDF export not fully implemented")
                return ""

            return ""

        except Exception as e:
            logger.error(f"Error exporting visualizations: {e}")
            return ""

    def _create_html_gallery(self, visualization_results: Dict[str, Any]) -> str:
        """Create HTML gallery of visualizations."""
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Single-Cell Visualization Gallery</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                .plot-section {{ margin-bottom: 40px; border: 1px solid #ddd; padding: 20px; }}
                .plot-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; }}
                .plot-item {{ border: 1px solid #eee; padding: 10px; }}
                img {{ max-width: 100%; height: auto; }}
                iframe {{ width: 100%; height: 500px; border: none; }}
            </style>
        </head>
        <body>
            <h1>Single-Cell Analysis Visualizations</h1>
            <p>Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        """

        # Add each plot type section
        for plot_type, plots in visualization_results.items():
            if plot_type != 'summary' and isinstance(plots, dict) and plots:
                html += f"""
                <div class="plot-section">
                    <h2>{plot_type.replace('_', ' ').title()}</h2>
                    <div class="plot-grid">
                """

                for plot_name, plot_file in plots.items():
                    if plot_file.endswith('.html'):
                        html += f'<div class="plot-item"><iframe src="{plot_file}"></iframe></div>'
                    else:
                        html += f'<div class="plot-item"><img src="{plot_file}" alt="{plot_name}"></div>'

                html += "</div></div>"

        html += """
        </body>
        </html>
        """

        return html

















