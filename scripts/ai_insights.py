#!/usr/bin/env python3
"""
AI-Powered Result Interpretation & Insights for RNASEQ-MINI

This module provides automated interpretation of RNA-seq results using
machine learning and statistical methods to identify biologically meaningful
patterns and provide actionable insights.
"""

import pandas as pd
import numpy as np
import json
import os
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
import warnings
import logging
from dataclasses import dataclass, asdict
import pickle
import joblib

# ML and statistical libraries
try:
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.model_selection import cross_val_score, train_test_split
    from sklearn.metrics import classification_report, roc_auc_score
    from sklearn.preprocessing import StandardScaler
    from sklearn.feature_selection import SelectKBest, f_classif
    from sklearn.decomposition import PCA
    import xgboost as xgb
    from scipy import stats
    from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
    from scipy.spatial.distance import pdist
    import networkx as nx
    from statsmodels.stats.multitest import multipletests
    ML_AVAILABLE = True
except ImportError as e:
    print(f"⚠️ Some ML libraries not available: {e}")
    ML_AVAILABLE = False


@dataclass
class InsightResult:
    """Container for AI-generated insights."""
    insight_type: str
    title: str
    description: str
    confidence: float
    data: Dict[str, Any]
    recommendations: List[str] = None
    supporting_evidence: List[str] = None
    generated_at: datetime = None

    def __post_init__(self):
        if self.recommendations is None:
            self.recommendations = []
        if self.supporting_evidence is None:
            self.supporting_evidence = []
        if self.generated_at is None:
            self.generated_at = datetime.now()


class AIInsightsEngine:
    """Main engine for AI-powered RNA-seq result interpretation."""

    def __init__(self, results_dir: str = "results", cache_dir: str = ".ai_cache"):
        self.results_dir = Path(results_dir)
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        self.models_dir = self.cache_dir / "models"
        self.models_dir.mkdir(exist_ok=True)

        # Initialize logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)

        # Load pre-trained models if available
        self.models = {}
        self._load_pretrained_models()

    def _load_pretrained_models(self):
        """Load pre-trained ML models for insights."""
        if not ML_AVAILABLE:
            return

        # Try to load pre-trained models for common tasks
        model_files = {
            "pathway_predictor": "pathway_impact_model.pkl",
            "biomarker_classifier": "biomarker_classifier.pkl",
            "gene_network_model": "gene_network_model.pkl"
        }

        for model_name, filename in model_files.items():
            model_path = self.models_dir / filename
            if model_path.exists():
                try:
                    self.models[model_name] = joblib.load(model_path)
                    self.logger.info(f"✅ Loaded {model_name} from {model_path}")
                except Exception as e:
                    self.logger.warning(f"⚠️ Failed to load {model_name}: {e}")

    def analyze_differential_expression(self, de_results_path: str = None) -> List[InsightResult]:
        """Analyze differential expression results for automated insights."""
        if de_results_path is None:
            de_results_path = self.results_dir / "de" / "deseq2_results.tsv"

        if not Path(de_results_path).exists():
            return [InsightResult(
                insight_type="error",
                title="No DE Results Found",
                description="Differential expression results not found. Please run DESeq2 analysis first.",
                confidence=1.0,
                data={"error": "missing_de_results"}
            )]

        try:
            de_results = pd.read_csv(de_results_path, sep='\t')
        except Exception as e:
            return [InsightResult(
                insight_type="error",
                title="Failed to Load DE Results",
                description=f"Could not load differential expression results: {e}",
                confidence=1.0,
                data={"error": str(e)}
            )]

        insights = []

        # 1. Effect size distribution analysis
        insights.extend(self._analyze_effect_sizes(de_results))

        # 2. Gene set enrichment pattern analysis
        insights.extend(self._analyze_enrichment_patterns(de_results))

        # 3. Biological coherence assessment
        insights.extend(self._assess_biological_coherence(de_results))

        # 4. Statistical power analysis
        insights.extend(self._analyze_statistical_power(de_results))

        # 5. Potential biomarker identification
        if ML_AVAILABLE:
            insights.extend(self._identify_potential_biomarkers(de_results))

        return insights

    def _analyze_effect_sizes(self, de_results: pd.DataFrame) -> List[InsightResult]:
        """Analyze distribution of effect sizes for insights."""
        insights = []

        # Filter significant genes
        sig_genes = de_results[(de_results['padj'] < 0.05) & (abs(de_results['log2FoldChange']) > 1)]

        if len(sig_genes) == 0:
            return [InsightResult(
                insight_type="warning",
                title="No Significant Differential Expression",
                description="No genes meet the significance threshold (padj < 0.05, |LFC| > 1). This could indicate insufficient power or biological effect.",
                confidence=0.8,
                data={"significant_genes": 0},
                recommendations=[
                    "Consider relaxing significance thresholds for exploratory analysis",
                    "Check if sample size provides adequate statistical power",
                    "Verify that contrasts are correctly specified"
                ]
            )]

        # Analyze effect size distribution
        effect_sizes = sig_genes['log2FoldChange'].abs()
        median_effect = effect_sizes.median()
        q75_effect = effect_sizes.quantile(0.75)

        # Classify effect size strength
        if median_effect > 2:
            strength = "strong"
            confidence = 0.9
        elif median_effect > 1:
            strength = "moderate"
            confidence = 0.8
        else:
            strength = "weak"
            confidence = 0.7

        insights.append(InsightResult(
            insight_type="effect_size_analysis",
            title=f"{strength.title()} Effect Sizes Detected",
            description=f"Median absolute log2 fold change: {median_effect".2f"}. {len(sig_genes)} significant genes identified.",
            confidence=confidence,
            data={
                "median_effect_size": median_effect,
                "significant_genes": len(sig_genes),
                "effect_strength": strength,
                "q75_effect_size": q75_effect
            },
            recommendations=self._get_effect_size_recommendations(strength, len(sig_genes))
        ))

        return insights

    def _get_effect_size_recommendations(self, strength: str, n_genes: int) -> List[str]:
        """Generate recommendations based on effect size analysis."""
        recommendations = []

        if strength == "weak":
            recommendations.extend([
                "Consider if the biological effect is subtle or if more samples are needed",
                "Review experimental design for potential confounding factors"
            ])
        elif strength == "moderate":
            recommendations.append("Effect sizes are biologically plausible - proceed with pathway analysis")
        elif strength == "strong":
            recommendations.extend([
                "Strong effects detected - prioritize these genes for validation",
                "Consider dose-response or time-course experiments to confirm findings"
            ])

        if n_genes > 1000:
            recommendations.append("Large number of DE genes suggests broad biological response")
        elif n_genes < 10:
            recommendations.append("Small number of DE genes may indicate subtle or specific effects")

        return recommendations

    def _analyze_enrichment_patterns(self, de_results: pd.DataFrame) -> List[InsightResult]:
        """Analyze patterns in gene set enrichment for insights."""
        insights = []

        # Load pathway results if available
        fgsea_path = self.results_dir / "fgsea" / "fgsea_results.tsv"
        if fgsea_path.exists():
            try:
                fgsea_results = pd.read_csv(fgsea_path, sep='\t')
                sig_pathways = fgsea_results[fgsea_results['padj'] < 0.05]

                if len(sig_pathways) > 0:
                    # Analyze pathway characteristics
                    top_pathways = sig_pathways.head(10)

                    insights.append(InsightResult(
                        insight_type="pathway_enrichment",
                        title=f"Significant Pathway Enrichment Detected",
                        description=f"Found {len(sig_pathways)} significantly enriched pathways. Top pathway: {top_pathways.iloc[0]['pathway']} (NES: {top_pathways.iloc[0]['NES']".2f"})",
                        confidence=0.85,
                        data={
                            "significant_pathways": len(sig_pathways),
                            "top_pathway": top_pathways.iloc[0]['pathway'],
                            "top_nes": top_pathways.iloc[0]['NES'],
                            "pathway_categories": self._categorize_pathways(top_pathways)
                        },
                        recommendations=self._get_pathway_recommendations(sig_pathways)
                    ))
            except Exception as e:
                self.logger.warning(f"Failed to analyze pathway results: {e}")

        return insights

    def _categorize_pathways(self, pathways: pd.DataFrame) -> Dict[str, int]:
        """Categorize pathways by functional themes."""
        categories = {
            "metabolic": 0,
            "signaling": 0,
            "immune": 0,
            "cell_cycle": 0,
            "development": 0,
            "stress_response": 0,
            "other": 0
        }

        pathway_names = pathways['pathway'].str.lower()

        for name in pathway_names:
            if any(word in name for word in ['metabol', 'biosynth', 'catabol']):
                categories["metabolic"] += 1
            elif any(word in name for word in ['signal', 'kinase', 'receptor']):
                categories["signaling"] += 1
            elif any(word in name for word in ['immune', 'inflam', 'cytokine']):
                categories["immune"] += 1
            elif any(word in name for word in ['cell.cycle', 'mitosis', 'division']):
                categories["cell_cycle"] += 1
            elif any(word in name for word in ['develop', 'differentiation']):
                categories["development"] += 1
            elif any(word in name for word in ['stress', 'apoptosis', 'response']):
                categories["stress_response"] += 1
            else:
                categories["other"] += 1

        return categories

    def _get_pathway_recommendations(self, pathways: pd.DataFrame) -> List[str]:
        """Generate recommendations based on pathway analysis."""
        recommendations = []

        categories = self._categorize_pathways(pathways)
        dominant_category = max(categories, key=categories.get)

        if categories[dominant_category] > len(pathways) * 0.3:
            recommendations.append(f"Strong enrichment in {dominant_category.replace('_', ' ')} pathways suggests this is a key biological theme")

        if len(pathways) > 20:
            recommendations.append("Large number of enriched pathways suggests broad biological response")

        # Check for pathway NES patterns
        high_nes = pathways[abs(pathways['NES']) > 2]
        if len(high_nes) > 0:
            recommendations.append(f"{len(high_nes)} pathways show very strong enrichment - prioritize these for follow-up")

        return recommendations

    def _assess_biological_coherence(self, de_results: pd.DataFrame) -> List[InsightResult]:
        """Assess biological coherence of DE results."""
        insights = []

        if not ML_AVAILABLE:
            return insights

        try:
            # Simple coherence assessment based on gene clustering
            sig_genes = de_results[de_results['padj'] < 0.05].copy()

            if len(sig_genes) > 50:  # Need enough genes for clustering
                # Use log fold changes for clustering
                lfc_matrix = sig_genes[['log2FoldChange']].T

                # Perform hierarchical clustering
                if len(sig_genes) > 1:
                    distance_matrix = pdist(lfc_matrix)
                    linkage_matrix = linkage(distance_matrix, method='ward')

                    # Cut dendrogram to get clusters
                    clusters = fcluster(linkage_matrix, t=3, criterion='maxclust')

                    # Analyze cluster characteristics
                    cluster_sizes = pd.Series(clusters).value_counts()

                    if len(cluster_sizes) > 1:
                        insights.append(InsightResult(
                            insight_type="biological_coherence",
                            title="Gene Expression Patterns Suggest Biological Coherence",
                            description=f"Identified {len(cluster_sizes)} distinct patterns in gene expression, suggesting coordinated biological responses.",
                            confidence=0.75,
                            data={
                                "n_clusters": len(cluster_sizes),
                                "cluster_sizes": cluster_sizes.to_dict(),
                                "largest_cluster": cluster_sizes.max()
                            },
                            recommendations=[
                                "Investigate the biological functions of genes in the largest cluster",
                                "Consider these patterns in downstream functional analysis"
                            ]
                        ))

        except Exception as e:
            self.logger.warning(f"Biological coherence analysis failed: {e}")

        return insights

    def _analyze_statistical_power(self, de_results: pd.DataFrame) -> List[InsightResult]:
        """Analyze statistical power of the experiment."""
        insights = []

        # Calculate power metrics
        total_genes = len(de_results)
        sig_genes = len(de_results[de_results['padj'] < 0.05])
        power_estimate = sig_genes / total_genes if total_genes > 0 else 0

        # Effect size distribution for power assessment
        effect_sizes = de_results['log2FoldChange'].abs()
        detectable_effects = len(effect_sizes[effect_sizes > 1])

        insights.append(InsightResult(
            insight_type="power_analysis",
            title="Statistical Power Assessment",
            description=f"Power estimate: {power_estimate".1%"}. {detectable_effects} genes show large effect sizes (>|1| LFC).",
            confidence=0.7,
            data={
                "total_genes": total_genes,
                "significant_genes": sig_genes,
                "power_estimate": power_estimate,
                "detectable_effects": detectable_effects
            },
            recommendations=self._get_power_recommendations(power_estimate, detectable_effects, total_genes)
        ))

        return insights

    def _get_power_recommendations(self, power: float, detectable: int, total: int) -> List[str]:
        """Generate recommendations based on power analysis."""
        recommendations = []

        if power < 0.1:
            recommendations.extend([
                "Low statistical power detected - consider increasing sample size",
                "May need more biological replicates for reliable results"
            ])
        elif power > 0.5:
            recommendations.append("Good statistical power - results likely reliable")

        if detectable < 10:
            recommendations.append("Few large effect genes detected - may indicate subtle biological effects")

        return recommendations

    def _identify_potential_biomarkers(self, de_results: pd.DataFrame) -> List[InsightResult]:
        """Identify potential biomarkers using ML approaches."""
        insights = []

        if not ML_AVAILABLE or len(de_results) < 100:
            return insights

        try:
            # Prepare features for biomarker prediction
            features = self._extract_biomarker_features(de_results)

            if features is None:
                return insights

            # Use pre-trained model if available, otherwise use rule-based approach
            if "biomarker_classifier" in self.models:
                biomarker_scores = self.models["biomarker_classifier"].predict_proba(features)[:, 1]
            else:
                # Rule-based biomarker scoring
                biomarker_scores = self._calculate_biomarker_scores(de_results)

            # Identify top biomarker candidates
            top_indices = np.argsort(biomarker_scores)[-20:]  # Top 20 candidates
            biomarker_candidates = de_results.iloc[top_indices]

            insights.append(InsightResult(
                insight_type="biomarker_discovery",
                title="Potential Biomarkers Identified",
                description=f"Identified {len(biomarker_candidates)} high-confidence biomarker candidates using ML analysis.",
                confidence=0.8,
                data={
                    "n_candidates": len(biomarker_candidates),
                    "top_genes": biomarker_candidates['gene_id'].tolist()[:5],
                    "avg_confidence": np.mean(biomarker_scores[top_indices])
                },
                recommendations=[
                    "Validate top biomarker candidates in independent cohorts",
                    "Consider these genes for diagnostic or prognostic applications",
                    "Assess biomarker performance using ROC analysis"
                ]
            ))

        except Exception as e:
            self.logger.warning(f"Biomarker identification failed: {e}")

        return insights

    def _extract_biomarker_features(self, de_results: pd.DataFrame) -> Optional[np.ndarray]:
        """Extract features for biomarker prediction."""
        try:
            # Select relevant features for biomarker prediction
            features = de_results[[
                'log2FoldChange', 'padj', 'baseMean'
            ]].copy()

            # Add derived features
            features['abs_lfc'] = features['log2FoldChange'].abs()
            features['lfc_rank'] = features['log2FoldChange'].abs().rank(pct=True)
            features['expression_level'] = np.log1p(features['baseMean'])

            # Handle missing values
            features = features.fillna(0)

            return features.values

        except Exception as e:
            self.logger.warning(f"Feature extraction failed: {e}")
            return None

    def _calculate_biomarker_scores(self, de_results: pd.DataFrame) -> np.ndarray:
        """Calculate biomarker scores using rule-based approach."""
        scores = np.zeros(len(de_results))

        # Effect size component (30% weight)
        scores += 0.3 * (de_results['log2FoldChange'].abs() / de_results['log2FoldChange'].abs().max())

        # Significance component (25% weight)
        scores += 0.25 * (1 - de_results['padj'].fillna(1))

        # Expression level component (20% weight)
        scores += 0.2 * (de_results['baseMean'] / de_results['baseMean'].max())

        # Consistency component (25% weight) - higher for moderate, consistent effects
        consistency = 1 / (1 + np.exp(-abs(de_results['log2FoldChange'])))  # Sigmoid around 0
        scores += 0.25 * consistency

        return scores

    def generate_result_summary(self, insights: List[InsightResult]) -> str:
        """Generate a natural language summary of insights."""
        if not insights:
            return "No significant insights generated from the analysis."

        summary_parts = []

        # Group insights by type
        insight_types = {}
        for insight in insights:
            if insight.insight_type not in insight_types:
                insight_types[insight.insight_type] = []
            insight_types[insight.insight_type].append(insight)

        # Generate summary sections
        if "effect_size_analysis" in insight_types:
            effect_insights = insight_types["effect_size_analysis"]
            effect = effect_insights[0]
            summary_parts.append(
                f"The analysis revealed {effect.data['effect_strength']} effect sizes with "
                f"{effect.data['significant_genes']} significantly differentially expressed genes."
            )

        if "pathway_enrichment" in insight_types:
            pathway_insights = insight_types["pathway_enrichment"]
            pathway = pathway_insights[0]
            summary_parts.append(
                f"Pathway analysis identified {pathway.data['significant_pathways']} enriched pathways, "
                f"with the strongest enrichment in {pathway.data['top_pathway']}."
            )

        if "biomarker_discovery" in insight_types:
            biomarker_insights = insight_types["biomarker_discovery"]
            biomarker = biomarker_insights[0]
            summary_parts.append(
                f"Machine learning analysis identified {biomarker.data['n_candidates']} potential biomarkers "
                f"for further validation."
            )

        if "power_analysis" in insight_types:
            power_insights = insight_types["power_analysis"]
            power = power_insights[0]
            summary_parts.append(
                f"Statistical power analysis suggests {power.data['power_estimate']".1%"} power "
                f"to detect differential expression."
            )

        # Combine into coherent summary
        if summary_parts:
            summary = " ".join(summary_parts)

            # Add recommendations if available
            all_recommendations = []
            for insight in insights:
                all_recommendations.extend(insight.recommendations)

            if all_recommendations:
                summary += f" Key recommendations include: {', '.join(set(all_recommendations[:3]))}."

            return summary

        return "Analysis completed successfully. Review the detailed insights for specific findings."

    def export_insights(self, insights: List[InsightResult], format: str = "json") -> str:
        """Export insights in various formats."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"ai_insights_{timestamp}"

        if format == "json":
            export_data = {
                "generated_at": datetime.now().isoformat(),
                "insights": [asdict(insight) for insight in insights],
                "summary": self.generate_result_summary(insights)
            }

            output_file = self.cache_dir / f"{filename}.json"
            with open(output_file, 'w') as f:
                json.dump(export_data, f, indent=2, default=str)

        elif format == "markdown":
            output_file = self.cache_dir / f"{filename}.md"
            with open(output_file, 'w') as f:
                f.write("# 🔬 AI-Powered RNA-seq Insights\n\n")
                f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

                # Add summary
                summary = self.generate_result_summary(insights)
                f.write(f"## Summary\n\n{summary}\n\n")

                # Add detailed insights
                f.write("## Detailed Insights\n\n")
                for insight in insights:
                    f.write(f"### {insight.title}\n\n")
                    f.write(f"**Confidence:** {insight.confidence".1%"}  \n")
                    f.write(f"**Type:** {insight.insight_type}  \n")
                    f.write(f"**Description:** {insight.description}\n\n")

                    if insight.recommendations:
                        f.write("**Recommendations:**\n")
                        for rec in insight.recommendations:
                            f.write(f"- {rec}\n")
                        f.write("\n")

                    if insight.supporting_evidence:
                        f.write("**Supporting Evidence:**\n")
                        for evidence in insight.supporting_evidence:
                            f.write(f"- {evidence}\n")
                        f.write("\n")

        return str(output_file)


def run_ai_insights(results_dir: str = "results", output_format: str = "both") -> str:
    """Run AI-powered insights analysis."""
    engine = AIInsightsEngine(results_dir)

    print("🧠 Running AI-powered RNA-seq insights analysis...")

    # Analyze differential expression results
    insights = engine.analyze_differential_expression()

    if not insights or all(i.insight_type == "error" for i in insights):
        print("❌ No analyzable results found")
        return None

    print(f"✅ Generated {len(insights)} insights")

    # Generate summary
    summary = engine.generate_result_summary(insights)
    print(f"\n📋 Summary: {summary}")

    # Export results
    if output_format in ["json", "both"]:
        json_file = engine.export_insights(insights, "json")
        print(f"📄 Exported JSON insights to: {json_file}")

    if output_format in ["markdown", "both"]:
        md_file = engine.export_insights(insights, "markdown")
        print(f"📄 Exported Markdown report to: {md_file}")

    return summary


def train_insights_models(training_data_dir: str = "training_data"):
    """Train ML models for improved insights (requires training data)."""
    if not ML_AVAILABLE:
        print("❌ ML libraries not available for model training")
        return False

    print("🧠 Training AI models for RNA-seq insights...")

    # This would implement model training using provided training data
    # For now, create placeholder models

    models_dir = Path("scripts") / ".ai_cache" / "models"
    models_dir.mkdir(parents=True, exist_ok=True)

    # Create simple placeholder models
    try:
        # Pathway impact model (placeholder)
        pathway_model = {"type": "placeholder", "version": "1.0"}
        with open(models_dir / "pathway_impact_model.pkl", 'wb') as f:
            pickle.dump(pathway_model, f)

        # Biomarker classifier (placeholder)
        biomarker_model = {"type": "placeholder", "version": "1.0"}
        with open(models_dir / "biomarker_classifier.pkl", 'wb') as f:
            pickle.dump(biomarker_model, f)

        print("✅ Placeholder models created (replace with real training data)")
        return True

    except Exception as e:
        print(f"❌ Model training failed: {e}")
        return False


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="AI-Powered RNA-seq Insights")
    subparsers = parser.add_subparsers(dest="command")

    # Analyze command
    analyze_parser = subparsers.add_parser("analyze", help="Run AI insights analysis")
    analyze_parser.add_argument("--results-dir", default="results", help="Results directory")
    analyze_parser.add_argument("--format", choices=["json", "markdown", "both"], default="both",
                               help="Output format")

    # Train command
    train_parser = subparsers.add_parser("train", help="Train AI models")
    train_parser.add_argument("--data-dir", default="training_data", help="Training data directory")

    args = parser.parse_args()

    if args.command == "analyze":
        summary = run_ai_insights(args.results_dir, args.format)
        if summary:
            print(f"\n🎯 Analysis complete: {summary}")
    elif args.command == "train":
        success = train_insights_models(args.data_dir)
        if success:
            print("✅ Model training complete")
    else:
        parser.print_help()
