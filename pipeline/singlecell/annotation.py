#!/usr/bin/env python3
"""
Cell type annotation module for single-cell RNA-seq data.
Provides automated cell type identification using marker genes and reference datasets.
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
class AnnotationConfig:
    """Configuration for cell type annotation."""

    # Input data
    clustered_file: str
    marker_genes: Optional[str] = None

    # Annotation parameters
    annotation_method: str = "marker_based"  # marker_based, reference_mapping, supervised
    min_confidence: float = 0.7
    max_annotations: int = 3

    # Reference parameters
    reference_dataset: Optional[str] = None
    species: str = "human"  # human, mouse, other

    # Output parameters
    output_dir: str = "results/singlecell/annotation"
    save_predictions: bool = True

    # Advanced parameters
    use_ontology: bool = True
    hierarchical_annotation: bool = True
    cross_validation: bool = True


class CellTypeAnnotation:
    """Handles automated cell type annotation for single-cell data."""

    def __init__(self, config: AnnotationConfig):
        self.config = config
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Load reference data
        self.reference_markers = {}
        self.cell_ontology = {}
        self._load_reference_data()

    def _load_reference_data(self):
        """Load reference marker genes and cell ontology."""
        try:
            # Load built-in marker genes
            builtin_markers_file = Path(__file__).parent / "reference" / "marker_genes.json"
            if builtin_markers_file.exists():
                with open(builtin_markers_file, 'r') as f:
                    self.reference_markers = json.load(f)

            # Load cell ontology
            ontology_file = Path(__file__).parent / "reference" / "cell_ontology.json"
            if ontology_file.exists():
                with open(ontology_file, 'r') as f:
                    self.cell_ontology = json.load(f)

            logger.info(f"Loaded {len(self.reference_markers)} reference cell types")

        except Exception as e:
            logger.error(f"Error loading reference data: {e}")

    def annotate_by_markers(self, clustered_file: str, marker_genes: Dict[str, List[str]] = None) -> Dict[str, Any]:
        """
        Annotate cell types using marker gene expression.

        Args:
            clustered_file: Path to clustered single-cell data
            marker_genes: Optional custom marker gene dictionary

        Returns:
            Annotation results
        """
        try:
            import scanpy as sc

            # Load clustered data
            adata = sc.read_h5ad(clustered_file)

            # Use provided marker genes or reference markers
            if marker_genes is None:
                marker_genes = self.reference_markers

            if not marker_genes:
                raise ValueError("No marker genes provided for annotation")

            # Calculate marker gene scores for each cluster
            cluster_annotations = {}

            for cluster in adata.obs['leiden'].unique():
                cluster_mask = adata.obs[f'{leiden}_clusters'] == cluster

                # Calculate average expression of marker genes for this cluster
                cluster_expr = adata[cluster_mask, :].X.mean(axis=0).A1

                # Score each cell type based on marker gene expression
                cell_type_scores = {}

                for cell_type, markers in marker_genes.items():
                    # Find markers that exist in the data
                    existing_markers = [m for m in markers if m in adata.var_names]

                    if existing_markers:
                        # Calculate average expression of existing markers
                        marker_indices = [adata.var_names.get_loc(m) for m in existing_markers]
                        marker_scores = cluster_expr[marker_indices]

                        # Calculate score as mean expression of marker genes
                        score = np.mean(marker_scores) if len(marker_scores) > 0 else 0.0

                        if score > 0:
                            cell_type_scores[cell_type] = {
                                'score': float(score),
                                'markers_found': len(existing_markers),
                                'total_markers': len(markers),
                                'confidence': min(1.0, score * len(existing_markers) / len(markers))
                            }

                # Sort cell types by score and select top candidates
                sorted_scores = sorted(cell_type_scores.items(), key=lambda x: x[1]['score'], reverse=True)

                cluster_annotations[str(cluster)] = {
                    'top_predictions': sorted_scores[:self.config.max_annotations],
                    'all_scores': cell_type_scores,
                    'n_cell_types_scored': len(cell_type_scores)
                }

            # Save annotations
            annotation_file = self.output_dir / "marker_annotations.json"
            with open(annotation_file, 'w') as f:
                json.dump(cluster_annotations, f, indent=2)

            return {
                'success': True,
                'annotation_file': str(annotation_file),
                'cluster_annotations': cluster_annotations,
                'method': 'marker_based'
            }

        except Exception as e:
            logger.error(f"Error in marker-based annotation: {e}")
            return {'success': False, 'error': str(e)}

    def annotate_by_reference_mapping(self, clustered_file: str, reference_file: str) -> Dict[str, Any]:
        """
        Annotate cell types by mapping to reference dataset.

        Args:
            clustered_file: Path to clustered single-cell data
            reference_file: Path to reference annotated dataset

        Returns:
            Reference mapping results
        """
        try:
            import scanpy as sc

            # Load query data
            query_adata = sc.read_h5ad(clustered_file)

            # Load reference data
            if reference_file.endswith('.h5ad'):
                ref_adata = sc.read_h5ad(reference_file)
            else:
                raise ValueError("Reference file must be in h5ad format")

            # Check for common genes
            common_genes = list(set(query_adata.var_names) & set(ref_adata.var_names))

            if len(common_genes) < 100:
                logger.warning(f"Only {len(common_genes)} common genes found between datasets")

            # Subset to common genes
            query_subset = query_adata[:, common_genes]
            ref_subset = ref_adata[:, common_genes]

            # Normalize both datasets
            sc.pp.normalize_total(query_subset, target_sum=1e4)
            sc.pp.log1p(query_subset)

            sc.pp.normalize_total(ref_subset, target_sum=1e4)
            sc.pp.log1p(ref_subset)

            # Find nearest neighbors in reference space
            from sklearn.neighbors import NearestNeighbors

            # Use PCA for dimensionality reduction
            sc.tl.pca(query_subset, n_comps=50)
            sc.tl.pca(ref_subset, n_comps=50)

            # Fit nearest neighbors on reference data
            nn = NearestNeighbors(n_neighbors=10)
            nn.fit(ref_subset.obsm['X_pca'])

            # Find nearest neighbors for query cells
            distances, indices = nn.kneighbors(query_subset.obsm['X_pca'])

            # Map predictions
            predictions = {}
            for i, (dist, idx) in enumerate(zip(distances, indices)):
                # Get cell types of nearest neighbors
                neighbor_types = ref_subset.obs.iloc[idx]['cell_type'].values

                # Count votes for each cell type
                type_counts = {}
                for cell_type in neighbor_types:
                    type_counts[cell_type] = type_counts.get(cell_type, 0) + 1

                # Select most frequent cell type
                predicted_type = max(type_counts.items(), key=lambda x: x[1])[0]

                # Calculate confidence based on distance and agreement
                confidence = min(1.0, 1.0 / (1.0 + np.mean(dist)))

                predictions[query_subset.obs.index[i]] = {
                    'predicted_cell_type': predicted_type,
                    'confidence': float(confidence),
                    'neighbor_votes': type_counts,
                    'mean_distance': float(np.mean(dist))
                }

            # Aggregate predictions by cluster
            cluster_predictions = {}
            for cluster in query_adata.obs[f'{leiden}_clusters'].unique():
                cluster_cells = query_adata.obs[query_adata.obs[f'{leiden}_clusters'] == cluster].index
                cluster_preds = [predictions[cell] for cell in cluster_cells if cell in predictions]

                if cluster_preds:
                    # Find most common prediction
                    cell_types = [pred['predicted_cell_type'] for pred in cluster_preds]
                    type_counts = {}
                    for cell_type in cell_types:
                        type_counts[cell_type] = type_counts.get(cell_type, 0) + 1

                    predicted_type = max(type_counts.items(), key=lambda x: x[1])[0]
                    confidence = np.mean([pred['confidence'] for pred in cluster_preds])

                    cluster_predictions[str(cluster)] = {
                        'predicted_cell_type': predicted_type,
                        'confidence': float(confidence),
                        'n_cells_annotated': len(cluster_preds),
                        'total_cells': len(cluster_cells),
                        'type_distribution': type_counts
                    }

            # Save predictions
            prediction_file = self.output_dir / "reference_predictions.json"
            with open(prediction_file, 'w') as f:
                json.dump(cluster_predictions, f, indent=2)

            return {
                'success': True,
                'prediction_file': str(prediction_file),
                'cluster_predictions': cluster_predictions,
                'method': 'reference_mapping',
                'common_genes': len(common_genes)
            }

        except Exception as e:
            logger.error(f"Error in reference mapping annotation: {e}")
            return {'success': False, 'error': str(e)}

    def annotate_hierarchically(self, base_annotations: Dict[str, Any],
                               ontology_tree: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Perform hierarchical annotation using cell ontology.

        Args:
            base_annotations: Base cell type annotations
            ontology_tree: Optional custom ontology tree

        Returns:
            Hierarchical annotation results
        """
        try:
            # Use provided ontology or default
            if ontology_tree is None:
                ontology_tree = self.cell_ontology

            if not ontology_tree:
                logger.warning("No cell ontology available for hierarchical annotation")
                return base_annotations

            hierarchical_annotations = {}

            for cluster, annotation in base_annotations.items():
                base_type = annotation.get('predicted_cell_type', '')

                # Find hierarchical relationships
                if base_type in ontology_tree:
                    hierarchy = ontology_tree[base_type]
                    hierarchical_annotations[cluster] = {
                        **annotation,
                        'hierarchical_path': hierarchy.get('path', []),
                        'parent_types': hierarchy.get('parents', []),
                        'child_types': hierarchy.get('children', [])
                    }
                else:
                    hierarchical_annotations[cluster] = annotation

            # Save hierarchical annotations
            hierarchy_file = self.output_dir / "hierarchical_annotations.json"
            with open(hierarchy_file, 'w') as f:
                json.dump(hierarchical_annotations, f, indent=2)

            return {
                'success': True,
                'hierarchy_file': str(hierarchy_file),
                'hierarchical_annotations': hierarchical_annotations,
                'ontology_used': bool(ontology_tree)
            }

        except Exception as e:
            logger.error(f"Error in hierarchical annotation: {e}")
            return {'success': False, 'error': str(e)}

    def validate_annotations(self, annotations: Dict[str, Any], validation_method: str = "cross_validation") -> Dict[str, Any]:
        """
        Validate cell type annotations using cross-validation or reference comparison.

        Args:
            annotations: Cell type annotation results
            validation_method: Validation approach

        Returns:
            Validation results
        """
        try:
            if validation_method == "cross_validation":
                return self._cross_validate_annotations(annotations)
            else:
                logger.warning(f"Validation method {validation_method} not implemented")
                return {'success': False, 'error': f'Method {validation_method} not implemented'}

        except Exception as e:
            logger.error(f"Error validating annotations: {e}")
            return {'success': False, 'error': str(e)}

    def _cross_validate_annotations(self, annotations: Dict[str, Any]) -> Dict[str, Any]:
        """Perform cross-validation on annotation accuracy."""
        try:
            # This would require a labeled dataset for validation
            # For now, return placeholder
            logger.warning("Cross-validation requires labeled reference data")

            return {
                'success': False,
                'error': 'Cross-validation requires labeled reference dataset'
            }

        except Exception as e:
            logger.error(f"Error in cross-validation: {e}")
            return {'success': False, 'error': str(e)}

    def generate_annotation_report(self, annotation_results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Generate comprehensive annotation report.

        Args:
            annotation_results: Annotation analysis results

        Returns:
            Complete annotation report
        """
        report = {
            'summary': {
                'annotation_completed': True,
                'method': annotation_results.get('method', 'unknown'),
                'n_clusters_annotated': len(annotation_results.get('cluster_annotations', {})),
                'min_confidence_threshold': self.config.min_confidence
            },
            'annotations': annotation_results,
            'quality_metrics': {},
            'recommendations': []
        }

        # Calculate quality metrics
        if 'cluster_annotations' in annotation_results:
            annotations = annotation_results['cluster_annotations']

            # Calculate confidence statistics
            confidences = []
            for cluster_data in annotations.values():
                if 'top_predictions' in cluster_data:
                    for pred in cluster_data['top_predictions']:
                        if isinstance(pred, dict) and 'confidence' in pred:
                            confidences.append(pred['confidence'])

            if confidences:
                report['quality_metrics'] = {
                    'mean_confidence': np.mean(confidences),
                    'median_confidence': np.median(confidences),
                    'high_confidence_annotations': len([c for c in confidences if c >= 0.8]),
                    'low_confidence_annotations': len([c for c in confidences if c < 0.5])
                }

        # Generate recommendations
        recommendations = []

        if report['quality_metrics'].get('low_confidence_annotations', 0) > 0:
            recommendations.append("Some annotations have low confidence. Consider using reference mapping or manual curation.")

        if report['quality_metrics'].get('mean_confidence', 1.0) < 0.6:
            recommendations.append("Overall annotation confidence is low. Consider improving marker gene sets or using multiple annotation methods.")

        if not recommendations:
            recommendations.append("Cell type annotations look reliable.")

        report['recommendations'] = recommendations

        return report

    def run_complete_annotation(self, clustered_file: str,
                              annotation_method: str = "marker_based") -> Dict[str, Any]:
        """
        Run complete cell type annotation pipeline.

        Args:
            clustered_file: Path to clustered single-cell data
            annotation_method: Primary annotation method

        Returns:
            Complete annotation results
        """
        logger.info(f"Starting cell type annotation using {annotation_method} method")

        results = {
            'method': annotation_method,
            'base_annotations': {},
            'hierarchical_annotations': {},
            'validation': {},
            'report': {}
        }

        try:
            # Step 1: Base annotation
            logger.info("Step 1: Running base annotation")
            if annotation_method == "marker_based":
                base_result = self.annotate_by_markers(clustered_file)
            elif annotation_method == "reference_mapping":
                base_result = self.annotate_by_reference_mapping(clustered_file, self.config.reference_dataset)
            else:
                raise ValueError(f"Unknown annotation method: {annotation_method}")

            results['base_annotations'] = base_result

            if not base_result['success']:
                raise RuntimeError("Base annotation failed")

            # Step 2: Hierarchical annotation
            logger.info("Step 2: Running hierarchical annotation")
            hierarchical_result = self.annotate_hierarchically(base_result['cluster_annotations'])
            results['hierarchical_annotations'] = hierarchical_result

            # Step 3: Validation
            logger.info("Step 3: Validating annotations")
            validation_result = self.validate_annotations(base_result['cluster_annotations'])
            results['validation'] = validation_result

            # Step 4: Generate report
            logger.info("Step 4: Generating annotation report")
            results['report'] = self.generate_annotation_report(results)

            logger.info("Complete annotation analysis finished successfully")
            return results

        except Exception as e:
            logger.error(f"Error in complete annotation analysis: {e}")
            return {
                'error': str(e),
                'partial_results': results,
                'success': False
            }

    def export_annotation_results(self, results: Dict[str, Any], format: str = "json") -> str:
        """
        Export annotation results in various formats.

        Args:
            results: Annotation analysis results
            format: Export format (json, csv, tsv)

        Returns:
            Path to exported file
        """
        try:
            if format == "json" and 'base_annotations' in results:
                export_file = self.output_dir / "annotation_results.json"
                with open(export_file, 'w') as f:
                    json.dump(results, f, indent=2)
                return str(export_file)

            elif format in ["csv", "tsv"]:
                # Export cluster annotations as table
                if 'base_annotations' in results and results['base_annotations'].get('success'):
                    annotations = results['base_annotations']['cluster_annotations']

                    # Create dataframe with annotations
                    rows = []
                    for cluster, cluster_data in annotations.items():
                        top_pred = cluster_data.get('top_predictions', [{}])[0]
                        rows.append({
                            'cluster': cluster,
                            'predicted_cell_type': top_pred.get('cell_type', 'Unknown'),
                            'confidence': top_pred.get('confidence', 0.0),
                            'n_markers_found': top_pred.get('markers_found', 0),
                            'method': results.get('method', 'unknown')
                        })

                    df = pd.DataFrame(rows)

                    if format == "csv":
                        export_file = self.output_dir / "annotation_results.csv"
                        df.to_csv(export_file, index=False)
                    else:
                        export_file = self.output_dir / "annotation_results.tsv"
                        df.to_csv(export_file, sep='\t', index=False)

                    return str(export_file)

            return ""

        except Exception as e:
            logger.error(f"Error exporting annotation results: {e}")
            return ""

    def create_annotation_visualizations(self, clustered_file: str, annotations: Dict[str, Any]) -> Dict[str, str]:
        """
        Create visualizations for annotation results.

        Args:
            clustered_file: Path to clustered data
            annotations: Annotation results

        Returns:
            Dictionary of plot file paths
        """
        plot_files = {}

        try:
            import scanpy as sc
            import matplotlib.pyplot as plt

            # Load data
            adata = sc.read_h5ad(clustered_file)

            # Add annotations to data
            if 'cluster_annotations' in annotations:
                # Create cell type column
                cell_types = {}
                for cluster, cluster_data in annotations['cluster_annotations'].items():
                    top_pred = cluster_data.get('top_predictions', [{}])[0]
                    cell_type = top_pred.get('cell_type', f'Cluster_{cluster}')
                    cluster_mask = adata.obs[f'{leiden}_clusters'] == cluster
                    cell_types.update({idx: cell_type for idx in adata.obs[cluster_mask].index})

                adata.obs['predicted_cell_type'] = pd.Series(cell_types)

                # Plot UMAP colored by predicted cell types
                sc.pl.umap(adata, color='predicted_cell_type', save='_cell_types.png', show=False)
                plot_files['umap_cell_types'] = str(self.output_dir / "umap_cell_types.png")

                # Plot confidence scores if available
                if 'confidence' in adata.obs.columns:
                    sc.pl.umap(adata, color='confidence', save='_confidence.png', show=False)
                    plot_files['umap_confidence'] = str(self.output_dir / "umap_confidence.png")

            plt.close('all')
            logger.info(f"Created {len(plot_files)} annotation visualization plots")

        except ImportError:
            logger.warning("scanpy not available for annotation plotting")
        except Exception as e:
            logger.error(f"Error creating annotation visualizations: {e}")

        return plot_files
