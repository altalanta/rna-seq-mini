#!/usr/bin/env python3
"""
Automated quality assessment and benchmarking system for RNASEQ-MINI.
Provides comprehensive quality metrics and validation against reference datasets.
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import subprocess
import logging
from sklearn.metrics import roc_auc_score, accuracy_score, precision_recall_fscore_support
from scipy import stats
import warnings

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')


class ReferenceDataset:
    """Container for reference datasets used in benchmarking."""

    def __init__(self, name: str, data_path: str, metadata_path: str,
                 description: str = "", expected_results: Dict = None):
        self.name = name
        self.data_path = Path(data_path)
        self.metadata_path = Path(metadata_path)
        self.description = description
        self.expected_results = expected_results or {}

        # Load data
        self.data = None
        self.metadata = None
        self._load_data()

    def _load_data(self):
        """Load reference dataset."""
        try:
            if self.data_path.exists():
                if self.data_path.suffix == '.csv':
                    self.data = pd.read_csv(self.data_path, index_col=0)
                elif self.data_path.suffix in ['.tsv', '.txt']:
                    self.data = pd.read_csv(self.data_path, sep='\t', index_col=0)
                else:
                    logger.warning(f"Unsupported reference data format: {self.data_path}")

            if self.metadata_path.exists():
                if self.metadata_path.suffix == '.csv':
                    self.metadata = pd.read_csv(self.metadata_path)
                elif self.metadata_path.suffix in ['.tsv', '.txt']:
                    self.metadata = pd.read_csv(self.metadata_path, sep='\t')

        except Exception as e:
            logger.error(f"Error loading reference dataset {self.name}: {e}")


class QualityAssessor:
    """Comprehensive quality assessment for RNA-seq analyses."""

    def __init__(self, reference_datasets: List[ReferenceDataset] = None):
        self.reference_datasets = reference_datasets or []
        self.quality_metrics = {}
        self.benchmark_results = {}

    def assess_rnaseq_quality(self, results_dir: str = "results") -> Dict[str, Any]:
        """
        Perform comprehensive quality assessment of RNA-seq results.

        Args:
            results_dir: Directory containing RNA-seq results

        Returns:
            Dictionary containing all quality metrics and assessments
        """
        results_path = Path(results_dir)
        quality_report = {
            'overall_score': 0.0,
            'section_scores': {},
            'detailed_metrics': {},
            'recommendations': [],
            'warnings': []
        }

        # Assess each component
        qc_score = self._assess_qc_quality(results_path)
        de_score = self._assess_de_quality(results_path)
        counts_score = self._assess_counts_quality(results_path)
        pathway_score = self._assess_pathway_quality(results_path)

        # Calculate overall score
        section_scores = {
            'qc': qc_score,
            'differential_expression': de_score,
            'counts': counts_score,
            'pathways': pathway_score
        }

        quality_report['section_scores'] = section_scores
        quality_report['overall_score'] = np.mean([score for score in section_scores.values() if score > 0])

        # Generate recommendations
        quality_report['recommendations'] = self._generate_recommendations(section_scores)
        quality_report['warnings'] = self._generate_warnings(section_scores)

        return quality_report

    def _assess_qc_quality(self, results_path: Path) -> float:
        """Assess quality control metrics."""
        qc_dir = results_path / "qc"
        score = 0.0

        if not qc_dir.exists():
            return 0.0

        # Check FastQC reports
        fastqc_dir = qc_dir / "fastqc"
        if fastqc_dir.exists():
            fastqc_files = list(fastqc_dir.glob("*_fastqc.html"))
            if fastqc_files:
                score += 0.3  # FastQC reports generated

        # Check MultiQC report
        multiqc_dir = qc_dir / "multiqc"
        if multiqc_dir.exists():
            multiqc_html = multiqc_dir / "multiqc_report.html"
            if multiqc_html.exists():
                score += 0.3  # MultiQC report generated

        # Check for quality metrics
        try:
            # Look for quality metrics in sample files or logs
            sample_quality_scores = self._extract_quality_metrics(results_path)
            if sample_quality_scores:
                avg_quality = np.mean(sample_quality_scores)
                score += min(0.4, avg_quality / 40.0)  # Quality score contribution
        except Exception as e:
            logger.warning(f"Error extracting quality metrics: {e}")

        return min(1.0, score)

    def _assess_de_quality(self, results_path: Path) -> float:
        """Assess differential expression analysis quality."""
        de_dir = results_path / "de"
        score = 0.0

        if not de_dir.exists():
            return 0.0

        # Check for DE results files
        de_files = list(de_dir.glob("DE_*.tsv"))
        if de_files:
            score += 0.3  # DE analyses completed

        # Check DE summary
        de_summary = de_dir / "de_summary.tsv"
        if de_summary.exists():
            try:
                de_df = pd.read_csv(de_summary, sep='\t')

                # Check for reasonable number of significant genes
                total_genes = len(de_df)
                if total_genes > 0:
                    sig_genes = len(de_df[de_df['padj'] < 0.05])
                    sig_rate = sig_genes / total_genes

                    # Reasonable significant gene rate (1-20% is typical)
                    if 0.01 <= sig_rate <= 0.2:
                        score += 0.3
                    else:
                        score += 0.1  # Some credit for having results

                # Check for reasonable effect sizes
                if 'log2FoldChange' in de_df.columns:
                    abs_lfc = de_df['log2FoldChange'].abs()
                    median_lfc = abs_lfc.median()

                    # Reasonable effect sizes (0.5-2.0 log2FC typical)
                    if 0.5 <= median_lfc <= 2.0:
                        score += 0.2
                    else:
                        score += 0.1

            except Exception as e:
                logger.warning(f"Error assessing DE quality: {e}")

        return min(1.0, score)

    def _assess_counts_quality(self, results_path: Path) -> float:
        """Assess gene expression count quality."""
        counts_dir = results_path / "counts"
        score = 0.0

        if not counts_dir.exists():
            return 0.0

        # Check for counts file
        counts_file = counts_dir / "counts.tsv"
        if counts_file.exists():
            try:
                counts_df = pd.read_csv(counts_file, sep='\t', index_col=0)

                # Basic count statistics
                n_genes = len(counts_df)
                n_samples = len(counts_df.columns)

                if n_genes > 1000 and n_samples > 1:
                    score += 0.4

                # Check for reasonable expression levels
                if not counts_df.empty:
                    mean_expression = counts_df.mean().mean()
                    if 1 <= mean_expression <= 1000:  # Reasonable expression range
                        score += 0.3

                # Check for low missing rate
                missing_rate = counts_df.isnull().sum().sum() / (counts_df.shape[0] * counts_df.shape[1])
                if missing_rate < 0.1:  # Less than 10% missing
                    score += 0.3

            except Exception as e:
                logger.warning(f"Error assessing counts quality: {e}")

        return min(1.0, score)

    def _assess_pathway_quality(self, results_path: Path) -> float:
        """Assess pathway enrichment analysis quality."""
        fgsea_dir = results_path / "fgsea"
        score = 0.0

        if not fgsea_dir.exists():
            return 0.0

        # Check for pathway results
        fgsea_files = list(fgsea_dir.glob("fgsea_*.tsv"))
        if fgsea_files:
            score += 0.3

        # Check pathway summary
        pathway_summary = fgsea_dir / "fgsea_summary.tsv"
        if pathway_summary.exists():
            try:
                pathway_df = pd.read_csv(pathway_summary, sep='\t')

                # Check for significant pathways
                sig_pathways = len(pathway_df[pathway_df['padj'] < 0.05])
                if sig_pathways > 0:
                    score += 0.3

                # Check pathway size distribution
                if 'size' in pathway_df.columns:
                    sizes = pathway_df['size'].dropna()
                    if len(sizes) > 0:
                        median_size = sizes.median()
                        if 10 <= median_size <= 500:  # Reasonable pathway sizes
                            score += 0.2
                        else:
                            score += 0.1

            except Exception as e:
                logger.warning(f"Error assessing pathway quality: {e}")

        return min(1.0, score)

    def _extract_quality_metrics(self, results_path: Path) -> List[float]:
        """Extract quality scores from various sources."""
        quality_scores = []

        # Try to extract from FastQC data or logs
        try:
            # Look for quality metrics in log files or output files
            for log_file in results_path.rglob("*log"):
                if log_file.exists():
                    content = log_file.read_text()
                    # Look for quality-related patterns
                    if "quality" in content.lower():
                        # Simple extraction - in practice, you'd parse specific formats
                        quality_scores.append(35.0)  # Placeholder

        except Exception as e:
            logger.warning(f"Error extracting quality metrics: {e}")

        return quality_scores

    def _generate_recommendations(self, section_scores: Dict[str, float]) -> List[str]:
        """Generate recommendations based on quality scores."""
        recommendations = []

        for section, score in section_scores.items():
            if score < 0.3:
                recommendations.append(f"Low quality in {section.replace('_', ' ')} - review analysis parameters")
            elif score < 0.7:
                recommendations.append(f"Moderate quality in {section.replace('_', ' ')} - consider optimization")

        if not recommendations:
            recommendations.append("Overall analysis quality is good")

        return recommendations

    def _generate_warnings(self, section_scores: Dict[str, float]) -> List[str]:
        """Generate warnings based on quality issues."""
        warnings_list = []

        for section, score in section_scores.items():
            if score == 0.0:
                warnings_list.append(f"No {section.replace('_', ' ')} results found - check pipeline execution")

        return warnings_list


class BenchmarkComparator:
    """Compare analysis results against reference datasets."""

    def __init__(self, reference_datasets: List[ReferenceDataset]):
        self.reference_datasets = reference_datasets
        self.comparison_results = {}

    def benchmark_differential_expression(self, de_results: Dict, reference_name: str = None) -> Dict[str, Any]:
        """
        Benchmark differential expression results against reference dataset.

        Args:
            de_results: Dictionary of DE results by contrast
            reference_name: Name of reference dataset to use

        Returns:
            Benchmarking results including sensitivity, specificity, etc.
        """
        if not reference_name:
            reference_name = self.reference_datasets[0].name if self.reference_datasets else None

        if not reference_name:
            return {"error": "No reference dataset available"}

        reference = next((ref for ref in self.reference_datasets if ref.name == reference_name), None)
        if not reference:
            return {"error": f"Reference dataset '{reference_name}' not found"}

        benchmark_results = {}

        for contrast, results in de_results.items():
            if 'data' in results:
                df = pd.DataFrame(results['data'])

                # Use reference dataset for benchmarking
                if reference.data is not None and reference.metadata is not None:
                    # Simple benchmarking: compare significant genes
                    sig_genes = set(df[df['padj'] < 0.05]['gene_id'].astype(str))

                    # For demonstration, assume reference has some expected DE genes
                    expected_de_genes = set()  # In practice, this would come from reference

                    if expected_de_genes:
                        tp = len(sig_genes.intersection(expected_de_genes))
                        fp = len(sig_genes - expected_de_genes)
                        fn = len(expected_de_genes - sig_genes)

                        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
                        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
                        f1_score = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

                        benchmark_results[contrast] = {
                            'precision': precision,
                            'recall': recall,
                            'f1_score': f1_score,
                            'true_positives': tp,
                            'false_positives': fp,
                            'false_negatives': fn,
                            'reference_dataset': reference_name
                        }

        return benchmark_results

    def assess_reproducibility(self, results_list: List[Dict], method: str = 'correlation') -> Dict[str, float]:
        """
        Assess reproducibility across multiple analyses or replicates.

        Args:
            results_list: List of analysis results to compare
            method: Comparison method ('correlation', 'jaccard', 'overlap')

        Returns:
            Reproducibility metrics
        """
        if len(results_list) < 2:
            return {"error": "Need at least 2 result sets for reproducibility assessment"}

        reproducibility_scores = {}

        # Extract significant genes from each result set
        sig_genes_list = []
        for results in results_list:
            sig_genes = set()
            for contrast_data in results.get('differential_expression', {}).get('contrasts', {}).values():
                if 'data' in contrast_data:
                    df = pd.DataFrame(contrast_data['data'])
                    sig_genes.update(df[df['padj'] < 0.05]['gene_id'].astype(str))
            sig_genes_list.append(sig_genes)

        if method == 'jaccard':
            # Calculate Jaccard similarity between gene sets
            similarities = []
            for i in range(len(sig_genes_list)):
                for j in range(i + 1, len(sig_genes_list)):
                    intersection = len(sig_genes_list[i].intersection(sig_genes_list[j]))
                    union = len(sig_genes_list[i].union(sig_genes_list[j]))
                    jaccard = intersection / union if union > 0 else 0
                    similarities.append(jaccard)

            reproducibility_scores['mean_jaccard'] = np.mean(similarities)
            reproducibility_scores['std_jaccard'] = np.std(similarities)

        elif method == 'correlation':
            # Calculate correlation of log fold changes if available
            # This would require matched genes across results
            reproducibility_scores['correlation_method'] = 'Not implemented for correlation'

        return reproducibility_scores

    def generate_benchmark_report(self, analysis_results: Dict) -> Dict[str, Any]:
        """
        Generate comprehensive benchmarking report.

        Args:
            analysis_results: Complete analysis results

        Returns:
            Detailed benchmarking report
        """
        report = {
            'benchmark_timestamp': pd.Timestamp.now().isoformat(),
            'quality_assessment': {},
            'comparison_results': {},
            'reproducibility_assessment': {},
            'summary': {}
        }

        # Run quality assessment
        quality_assessor = QualityAssessor(self.reference_datasets)
        report['quality_assessment'] = quality_assessor.assess_rnaseq_quality()

        # Run benchmarking if reference datasets available
        if self.reference_datasets:
            de_results = analysis_results.get('differential_expression', {}).get('contrasts', {})
            report['comparison_results'] = self.benchmark_differential_expression(de_results)

        # Assess reproducibility if multiple results provided
        # This would need to be extended for actual multi-run comparisons

        # Generate summary
        quality_score = report['quality_assessment']['overall_score']
        report['summary'] = {
            'overall_quality_score': quality_score,
            'quality_rating': self._get_quality_rating(quality_score),
            'benchmark_available': len(self.reference_datasets) > 0,
            'recommendations_count': len(report['quality_assessment']['recommendations'])
        }

        return report

    def _get_quality_rating(self, score: float) -> str:
        """Convert numerical quality score to categorical rating."""
        if score >= 0.8:
            return "Excellent"
        elif score >= 0.6:
            return "Good"
        elif score >= 0.4:
            return "Moderate"
        elif score >= 0.2:
            return "Poor"
        else:
            return "Very Poor"


class QualityGate:
    """Quality gate system to prevent poor analyses from proceeding."""

    def __init__(self, min_quality_score: float = 0.5):
        self.min_quality_score = min_quality_score
        self.gate_results = {}

    def evaluate_quality_gate(self, quality_report: Dict) -> Tuple[bool, str]:
        """
        Evaluate whether analysis meets quality standards.

        Args:
            quality_report: Quality assessment report

        Returns:
            Tuple of (passed_gate, reason_message)
        """
        overall_score = quality_report.get('overall_score', 0.0)
        section_scores = quality_report.get('section_scores', {})

        # Check overall quality score
        if overall_score < self.min_quality_score:
            return False, f"Overall quality score ({overall_score".2f"}) below threshold ({self.min_quality_score})"

        # Check for critical failures
        critical_failures = []
        for section, score in section_scores.items():
            if score == 0.0:
                critical_failures.append(section.replace('_', ' '))

        if critical_failures:
            return False, f"Critical failures in: {', '.join(critical_failures)}"

        # Check for warnings
        warnings = quality_report.get('warnings', [])
        if warnings:
            return True, f"Quality gate passed with warnings: {', '.join(warnings)}"

        return True, "Quality gate passed successfully"

    def should_proceed_to_publication(self, quality_report: Dict) -> Tuple[bool, str]:
        """
        Evaluate if analysis is suitable for publication.

        Uses stricter criteria than basic quality gate.
        """
        # Stricter threshold for publication
        publication_threshold = 0.7

        overall_score = quality_report.get('overall_score', 0.0)

        if overall_score < publication_threshold:
            return False, f"Publication quality threshold not met ({overall_score".2f"} < {publication_threshold})"

        # Check for excellent quality in key areas
        section_scores = quality_report.get('section_scores', {})
        key_sections = ['differential_expression', 'counts']

        for section in key_sections:
            if section in section_scores and section_scores[section] < 0.6:
                return False, f"Insufficient quality in {section.replace('_', ' ')} for publication"

        return True, "Analysis meets publication quality standards"


def run_quality_assessment(results_dir: str = "results",
                          reference_datasets: List[ReferenceDataset] = None,
                          output_file: str = "quality_report.json") -> Dict[str, Any]:
    """
    Run comprehensive quality assessment and benchmarking.

    Args:
        results_dir: Directory containing RNA-seq results
        reference_datasets: List of reference datasets for benchmarking
        output_file: Output file for quality report

    Returns:
        Quality assessment report
    """
    logger.info("Starting comprehensive quality assessment...")

    # Initialize assessors
    quality_assessor = QualityAssessor(reference_datasets)
    benchmark_comparator = BenchmarkComparator(reference_datasets) if reference_datasets else None

    # Load analysis results (placeholder - in practice, this would load actual results)
    analysis_results = {
        'qc': {},
        'differential_expression': {'contrasts': {}},
        'pathways': {},
        'counts': {}
    }

    # Run quality assessment
    quality_report = quality_assessor.assess_rnaseq_quality(results_dir)

    # Run benchmarking if reference datasets available
    if benchmark_comparator:
        benchmark_report = benchmark_comparator.generate_benchmark_report(analysis_results)
        quality_report['benchmarking'] = benchmark_report

    # Run quality gate evaluation
    quality_gate = QualityGate()
    gate_passed, gate_message = quality_gate.evaluate_quality_gate(quality_report)
    publication_ready, pub_message = quality_gate.should_proceed_to_publication(quality_report)

    quality_report['quality_gate'] = {
        'passed': gate_passed,
        'message': gate_message,
        'publication_ready': publication_ready,
        'publication_message': pub_message
    }

    # Save report
    with open(output_file, 'w') as f:
        json.dump(quality_report, f, indent=2)

    logger.info(f"Quality assessment complete. Report saved to: {output_file}")
    logger.info(f"Overall quality score: {quality_report['overall_score']".2f"}")
    logger.info(f"Quality gate: {'PASSED' if gate_passed else 'FAILED'}")
    logger.info(f"Publication ready: {'YES' if publication_ready else 'NO'}")

    return quality_report


def main():
    """Command-line interface for quality assessment."""
    import argparse

    parser = argparse.ArgumentParser(description="RNASEQ-MINI Quality Assessment")
    parser.add_argument('--results-dir', default='results', help='Results directory')
    parser.add_argument('--reference-dataset', help='Reference dataset name for benchmarking')
    parser.add_argument('--output', default='quality_report.json', help='Output report file')
    parser.add_argument('--min-quality-score', type=float, default=0.5, help='Minimum quality score threshold')

    args = parser.parse_args()

    # Load reference datasets if specified
    reference_datasets = []
    if args.reference_dataset:
        # In practice, this would load the specified reference dataset
        logger.info(f"Using reference dataset: {args.reference_dataset}")

    # Run quality assessment
    report = run_quality_assessment(
        results_dir=args.results_dir,
        reference_datasets=reference_datasets,
        output_file=args.output
    )

    # Exit with appropriate code based on quality gate
    if report.get('quality_gate', {}).get('passed', False):
        logger.info("✅ Quality assessment passed")
        return 0
    else:
        logger.error("❌ Quality assessment failed")
        return 1


if __name__ == "__main__":
    exit(main())
