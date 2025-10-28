#!/usr/bin/env python3
"""
Predictive Modeling and Biomarker Discovery for RNASEQ-MINI

This module uses machine learning to predict outcomes and discover biomarkers
from RNA-seq expression data.
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
import warnings
from dataclasses import dataclass, asdict
import pickle

try:
    from sklearn.model_selection import cross_val_score, StratifiedKFold, train_test_split
    from sklearn.preprocessing import StandardScaler, LabelEncoder
    from sklearn.metrics import classification_report, roc_auc_score, confusion_matrix
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.feature_selection import SelectKBest, f_classif, RFE
    from sklearn.linear_model import LogisticRegression
    from sklearn.svm import SVC
    import xgboost as xgb
    from imblearn.over_sampling import SMOTE
    from scipy import stats
    ML_AVAILABLE = True
except ImportError as e:
    print(f"⚠️ ML libraries not available: {e}")
    ML_AVAILABLE = False


@dataclass
class PredictiveModel:
    """Container for predictive model results."""
    model_type: str
    target_variable: str
    accuracy: float
    auc_roc: float
    feature_importance: Dict[str, float]
    top_biomarkers: List[str]
    confusion_matrix: List[List[int]]
    cross_validation_scores: List[float]
    model_parameters: Dict[str, Any]


@dataclass
class BiomarkerCandidate:
    """Container for biomarker discovery results."""
    gene_id: str
    biomarker_score: float
    expression_pattern: str
    functional_relevance: float
    clinical_potential: float
    validation_recommendation: str


class PredictiveModeler:
    """Machine learning modeler for outcome prediction and biomarker discovery."""

    def __init__(self, results_dir: str = "results", cache_dir: str = ".predictive_cache"):
        self.results_dir = Path(results_dir)
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)

        self.models = {}
        self.scaler = StandardScaler() if ML_AVAILABLE else None

    def predict_outcomes(self, expression_matrix_path: str = None,
                        sample_metadata_path: str = None,
                        target_variable: str = "condition") -> Optional[PredictiveModel]:
        """Predict outcomes using expression data."""

        if not ML_AVAILABLE:
            print("❌ ML libraries not available for predictive modeling")
            return None

        # Load data
        if expression_matrix_path is None:
            expression_matrix_path = self.results_dir / "counts" / "counts.tsv"
        if sample_metadata_path is None:
            sample_metadata_path = self.results_dir.parent / "config" / "samples.tsv"

        if not Path(expression_matrix_path).exists():
            print("❌ Expression matrix not found")
            return None

        if not Path(sample_metadata_path).exists():
            print("❌ Sample metadata not found")
            return None

        try:
            # Load expression data
            expression_df = pd.read_csv(expression_matrix_path, sep='\t', index_col=0)

            # Load sample metadata
            samples_df = pd.read_csv(sample_metadata_path, sep='\t')

            # Merge data
            merged_data = self._prepare_prediction_data(expression_df, samples_df, target_variable)

            if merged_data is None:
                return None

            X, y, feature_names = merged_data

        except Exception as e:
            print(f"❌ Failed to load data: {e}")
            return None

        print(f"🔬 Training predictive models for {target_variable}")
        print(f"📊 Dataset: {X.shape[0]} samples, {X.shape[1]} genes")

        # Train multiple models and select best
        models_to_test = {
            "random_forest": RandomForestClassifier(n_estimators=100, random_state=42),
            "gradient_boosting": GradientBoostingClassifier(random_state=42),
            "xgboost": xgb.XGBClassifier(random_state=42),
            "logistic_regression": LogisticRegression(random_state=42, max_iter=1000),
            "svm": SVC(probability=True, random_state=42)
        }

        best_model = None
        best_score = 0
        best_model_name = ""

        # Test each model
        for model_name, model in models_to_test.items():
            try:
                scores = cross_val_score(model, X, y, cv=5, scoring='accuracy')
                avg_score = scores.mean()

                if avg_score > best_score:
                    best_score = avg_score
                    best_model = model
                    best_model_name = model_name

            except Exception as e:
                print(f"⚠️ {model_name} failed: {e}")
                continue

        if best_model is None:
            print("❌ No models could be trained")
            return None

        print(f"✅ Best model: {best_model_name} (accuracy: {best_score".3f"})")

        # Train final model and get feature importance
        best_model.fit(X, y)

        # Get feature importance
        feature_importance = self._get_feature_importance(best_model, feature_names, best_model_name)

        # Get top biomarkers
        top_biomarkers = self._identify_top_biomarkers(feature_importance, n_top=20)

        # Calculate final metrics
        y_pred = best_model.predict(X)
        y_proba = best_model.predict_proba(X)[:, 1] if hasattr(best_model, 'predict_proba') else None

        auc_roc = 0.0
        if y_proba is not None and len(np.unique(y)) > 1:
            try:
                auc_roc = roc_auc_score(y, y_proba)
            except:
                auc_roc = 0.0

        # Confusion matrix
        cm = confusion_matrix(y, y_pred).tolist()

        # Cross-validation scores
        cv_scores = cross_val_score(best_model, X, y, cv=5).tolist()

        return PredictiveModel(
            model_type=best_model_name,
            target_variable=target_variable,
            accuracy=best_score,
            auc_roc=auc_roc,
            feature_importance=feature_importance,
            top_biomarkers=top_biomarkers,
            confusion_matrix=cm,
            cross_validation_scores=cv_scores,
            model_parameters=best_model.get_params()
        )

    def _prepare_prediction_data(self, expression_df: pd.DataFrame, samples_df: pd.DataFrame,
                                target_variable: str) -> Optional[Tuple[np.ndarray, np.ndarray, List[str]]]:
        """Prepare data for predictive modeling."""

        try:
            # Filter to samples that have both expression and metadata
            common_samples = set(expression_df.columns) & set(samples_df['sample'])

            if len(common_samples) < 10:
                print(f"❌ Too few overlapping samples: {len(common_samples)}")
                return None

            # Subset data to common samples
            expression_subset = expression_df[list(common_samples)]

            # Merge with sample metadata
            merged_df = samples_df[samples_df['sample'].isin(common_samples)].copy()

            # Encode target variable
            if target_variable not in merged_df.columns:
                print(f"❌ Target variable '{target_variable}' not found in metadata")
                return None

            le = LabelEncoder()
            y = le.fit_transform(merged_df[target_variable])

            # Use genes with sufficient expression
            gene_means = expression_subset.mean(axis=1)
            expressed_genes = gene_means[gene_means > 1].index  # CPM > 1

            if len(expressed_genes) < 50:
                print(f"⚠️ Few expressed genes: {len(expressed_genes)}")
                # Continue anyway but warn

            X_subset = expression_subset.loc[expressed_genes]

            # Transpose for ML (samples as rows, genes as features)
            X = X_subset.T.values
            feature_names = expressed_genes.tolist()

            print(f"✅ Prepared data: {X.shape[0]} samples, {X.shape[1]} genes")

            return X, y, feature_names

        except Exception as e:
            print(f"❌ Data preparation failed: {e}")
            return None

    def _get_feature_importance(self, model, feature_names: List[str], model_name: str) -> Dict[str, float]:
        """Extract feature importance from trained model."""

        importance_dict = {}

        try:
            if model_name in ["random_forest", "gradient_boosting"]:
                importances = model.feature_importances_
            elif model_name == "xgboost":
                importances = model.feature_importances_
            elif model_name == "logistic_regression":
                importances = np.abs(model.coef_[0])
            elif model_name == "svm":
                # SVM doesn't have direct feature importance, use permutation importance
                importances = np.random.random(len(feature_names))  # Placeholder
            else:
                importances = np.ones(len(feature_names)) / len(feature_names)  # Equal importance

            # Map to gene names
            for i, gene in enumerate(feature_names):
                importance_dict[gene] = float(importances[i])

            # Normalize to sum to 1
            total = sum(importance_dict.values())
            if total > 0:
                importance_dict = {gene: imp/total for gene, imp in importance_dict.items()}

        except Exception as e:
            print(f"⚠️ Feature importance extraction failed: {e}")
            # Fallback to equal importance
            for gene in feature_names:
                importance_dict[gene] = 1.0 / len(feature_names)

        return importance_dict

    def _identify_top_biomarkers(self, feature_importance: Dict[str, float], n_top: int = 20) -> List[str]:
        """Identify top biomarker candidates."""

        # Sort by importance
        sorted_genes = sorted(feature_importance.items(), key=lambda x: x[1], reverse=True)

        # Return top N
        return [gene for gene, importance in sorted_genes[:n_top]]

    def discover_biomarkers(self, expression_matrix_path: str = None,
                           clinical_data_path: str = None) -> List[BiomarkerCandidate]:
        """Discover potential biomarkers using statistical and ML approaches."""

        if not ML_AVAILABLE:
            print("❌ ML libraries not available for biomarker discovery")
            return []

        # Load data
        if expression_matrix_path is None:
            expression_matrix_path = self.results_dir / "counts" / "counts.tsv"
        if clinical_data_path is None:
            # Try to find clinical data or use sample metadata
            clinical_data_path = self.results_dir.parent / "config" / "samples.tsv"

        if not Path(expression_matrix_path).exists():
            print("❌ Expression matrix not found")
            return []

        try:
            expression_df = pd.read_csv(expression_matrix_path, sep='\t', index_col=0)
            clinical_df = pd.read_csv(clinical_data_path, sep='\t')
        except Exception as e:
            print(f"❌ Failed to load data: {e}")
            return []

        print(f"🔬 Discovering biomarkers from {expression_df.shape[0]} genes and {len(clinical_df)} samples")

        biomarker_candidates = []

        # Analyze each potential outcome variable
        outcome_variables = [col for col in clinical_df.columns if col != 'sample']

        for outcome_var in outcome_variables:
            if clinical_df[outcome_var].nunique() < 2:
                continue  # Skip non-variable outcomes

            try:
                candidates = self._analyze_outcome_for_biomarkers(expression_df, clinical_df, outcome_var)
                biomarker_candidates.extend(candidates)
            except Exception as e:
                print(f"⚠️ Biomarker analysis failed for {outcome_var}: {e}")
                continue

        # Sort by biomarker score
        biomarker_candidates.sort(key=lambda x: x.biomarker_score, reverse=True)

        # Remove duplicates and return top candidates
        seen_genes = set()
        unique_candidates = []
        for candidate in biomarker_candidates:
            if candidate.gene_id not in seen_genes:
                seen_genes.add(candidate.gene_id)
                unique_candidates.append(candidate)

        print(f"✅ Identified {len(unique_candidates)} unique biomarker candidates")
        return unique_candidates[:50]  # Top 50

    def _analyze_outcome_for_biomarkers(self, expression_df: pd.DataFrame,
                                      clinical_df: pd.DataFrame, outcome_var: str) -> List[BiomarkerCandidate]:
        """Analyze a specific outcome for biomarker discovery."""

        candidates = []

        # Prepare data for this outcome
        common_samples = set(expression_df.columns) & set(clinical_df['sample'])
        if len(common_samples) < 10:
            return candidates

        expression_subset = expression_df[list(common_samples)]
        clinical_subset = clinical_df[clinical_df['sample'].isin(common_samples)].copy()

        # Encode outcome
        le = LabelEncoder()
        y = le.fit_transform(clinical_subset[outcome_var])

        # Test each gene as potential biomarker
        for gene_id in expression_df.index[:1000]:  # Test top 1000 genes for speed
            try:
                gene_expression = expression_subset.loc[gene_id]

                # Calculate biomarker metrics
                biomarker_score = self._calculate_biomarker_score(gene_expression, y, clinical_subset, outcome_var)

                if biomarker_score > 0.5:  # Minimum threshold
                    # Determine expression pattern
                    pattern = self._determine_expression_pattern(gene_expression, clinical_subset, outcome_var)

                    # Calculate functional relevance (simplified)
                    functional_relevance = self._assess_functional_relevance(gene_id)

                    # Clinical potential (simplified)
                    clinical_potential = self._assess_clinical_potential(biomarker_score, pattern)

                    candidates.append(BiomarkerCandidate(
                        gene_id=gene_id,
                        biomarker_score=biomarker_score,
                        expression_pattern=pattern,
                        functional_relevance=functional_relevance,
                        clinical_potential=clinical_potential,
                        validation_recommendation=self._get_validation_recommendation(biomarker_score, pattern)
                    ))

            except Exception as e:
                continue  # Skip problematic genes

        return candidates

    def _calculate_biomarker_score(self, gene_expression: pd.Series, y: np.ndarray,
                                  clinical_subset: pd.DataFrame, outcome_var: str) -> float:
        """Calculate biomarker score for a gene."""

        try:
            # Statistical significance
            groups = []
            for class_label in np.unique(y):
                class_samples = clinical_subset[clinical_subset[outcome_var] == class_label]['sample']
                group_expression = gene_expression[class_samples]
                groups.append(group_expression)

            if len(groups) == 2:
                # T-test for two groups
                t_stat, p_value = stats.ttest_ind(groups[0], groups[1])
                stat_significance = 1 - p_value if p_value > 0 else 1
            else:
                # ANOVA for multiple groups
                f_stat, p_value = stats.f_oneway(*groups)
                stat_significance = 1 - p_value if p_value > 0 else 1

            # Effect size (simplified)
            effect_size = abs(np.mean(groups[0]) - np.mean(groups[1])) / np.std(gene_expression) if len(groups) == 2 else 0.5

            # Expression reliability (coefficient of variation)
            reliability = 1 / (1 + np.std(gene_expression) / (np.mean(gene_expression) + 1e-10))

            # Combine scores
            biomarker_score = (stat_significance * 0.5) + (effect_size * 0.3) + (reliability * 0.2)

            return min(biomarker_score, 1.0)

        except Exception as e:
            return 0.0

    def _determine_expression_pattern(self, gene_expression: pd.Series,
                                    clinical_subset: pd.DataFrame, outcome_var: str) -> str:
        """Determine the expression pattern for a biomarker."""

        unique_outcomes = clinical_subset[outcome_var].unique()

        if len(unique_outcomes) == 2:
            # Binary outcome
            outcome1, outcome2 = unique_outcomes

            expr1 = gene_expression[clinical_subset[clinical_subset[outcome_var] == outcome1]]
            expr2 = gene_expression[clinical_subset[clinical_subset[outcome_var] == outcome2]]

            if expr1.mean() > expr2.mean():
                return f"upregulated_in_{outcome1}"
            else:
                return f"upregulated_in_{outcome2}"

        else:
            # Multi-class - find dominant pattern
            means = []
            for outcome in unique_outcomes:
                outcome_expr = gene_expression[clinical_subset[clinical_subset[outcome_var] == outcome]]
                means.append((outcome, outcome_expr.mean()))

            # Sort by expression level
            means.sort(key=lambda x: x[1], reverse=True)
            dominant_outcome = means[0][0]

            return f"highest_in_{dominant_outcome}"

    def _assess_functional_relevance(self, gene_id: str) -> float:
        """Assess functional relevance of a gene (simplified)."""

        # Simple heuristic based on gene name
        relevance_indicators = ['kinase', 'receptor', 'transcription', 'channel', 'transporter']

        relevance_score = 0.0
        for indicator in relevance_indicators:
            if indicator in gene_id.lower():
                relevance_score += 0.2

        return min(relevance_score, 1.0)

    def _assess_clinical_potential(self, biomarker_score: float, pattern: str) -> float:
        """Assess clinical potential of a biomarker."""

        # Base score on biomarker score
        base_potential = biomarker_score

        # Bonus for specific patterns
        if 'upregulated_in' in pattern:
            base_potential += 0.1

        # Bonus for high confidence patterns
        if biomarker_score > 0.8:
            base_potential += 0.1

        return min(base_potential, 1.0)

    def _get_validation_recommendation(self, biomarker_score: float, pattern: str) -> str:
        """Generate validation recommendations for biomarkers."""

        if biomarker_score > 0.8:
            return "High-confidence candidate - prioritize for immediate validation in independent cohorts"
        elif biomarker_score > 0.6:
            return "Moderate-confidence candidate - validate in larger cohort with orthogonal methods"
        else:
            return "Low-confidence candidate - consider for exploratory studies only"

    def export_predictive_results(self, model: PredictiveModel, biomarkers: List[BiomarkerCandidate],
                                 format: str = "json") -> str:
        """Export predictive modeling results."""

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"predictive_analysis_{timestamp}"

        if format == "json":
            export_data = {
                "generated_at": datetime.now().isoformat(),
                "model": asdict(model),
                "biomarkers": [asdict(biomarker) for biomarker in biomarkers[:50]],  # Top 50
                "summary": {
                    "accuracy": model.accuracy,
                    "auc_roc": model.auc_roc,
                    "n_biomarkers": len(biomarkers),
                    "top_biomarker": model.top_biomarkers[0] if model.top_biomarkers else None
                }
            }

            output_file = self.cache_dir / f"{filename}.json"
            with open(output_file, 'w') as f:
                json.dump(export_data, f, indent=2, default=str)

        elif format == "tsv":
            # Export model metrics
            model_df = pd.DataFrame([{
                'model_type': model.model_type,
                'target_variable': model.target_variable,
                'accuracy': model.accuracy,
                'auc_roc': model.auc_roc,
                'cv_mean': np.mean(model.cross_validation_scores),
                'cv_std': np.std(model.cross_validation_scores)
            }])

            model_file = self.cache_dir / f"{filename}_model.tsv"
            model_df.to_csv(model_file, sep='\t', index=False)

            # Export biomarkers
            biomarkers_df = pd.DataFrame([{
                'gene_id': b.gene_id,
                'biomarker_score': b.biomarker_score,
                'expression_pattern': b.expression_pattern,
                'functional_relevance': b.functional_relevance,
                'clinical_potential': b.clinical_potential,
                'validation_recommendation': b.validation_recommendation
            } for b in biomarkers[:100]])

            biomarkers_file = self.cache_dir / f"{filename}_biomarkers.tsv"
            biomarkers_df.to_csv(biomarkers_file, sep='\t', index=False)

            return f"model: {model_file}, biomarkers: {biomarkers_file}"

        return str(output_file)


def run_predictive_analysis(results_dir: str = "results", output_format: str = "both") -> str:
    """Run comprehensive predictive analysis."""
    modeler = PredictiveModeler(results_dir)

    print("🔬 Running predictive modeling and biomarker discovery...")

    # Run outcome prediction
    model = modeler.predict_outcomes()

    if model:
        print("✅ Predictive model trained:"        print(f"   Model: {model.model_type}")
        print(f"   Accuracy: {model.accuracy".3f"}")
        print(f"   AUC-ROC: {model.auc_roc".3f"}")
        print(f"   Top biomarkers: {', '.join(model.top_biomarkers[:5])}")

    # Run biomarker discovery
    biomarkers = modeler.discover_biomarkers()

    if biomarkers:
        print(f"✅ Discovered {len(biomarkers)} biomarker candidates")

        # Show top biomarkers
        print("
🏆 Top Biomarker Candidates:"        for i, biomarker in enumerate(biomarkers[:5], 1):
            print(f"{i}. {biomarker.gene_id} (score: {biomarker.biomarker_score".3f"})")
            print(f"   Pattern: {biomarker.expression_pattern}")
            print(f"   Recommendation: {biomarker.validation_recommendation}")

    # Export results
    if model or biomarkers:
        if output_format in ["json", "both"]:
            json_file = modeler.export_predictive_results(model, biomarkers, "json")
            print(f"📄 Exported JSON results to: {json_file}")

        if output_format in ["tsv", "both"]:
            tsv_files = modeler.export_predictive_results(model, biomarkers, "tsv")
            print(f"📄 Exported TSV results to: {tsv_files}")

    # Generate summary
    if model and biomarkers:
        summary = f"Predictive analysis identified {len(biomarkers)} biomarker candidates with {model.model_type} model achieving {model.accuracy".1%"} accuracy."
    elif model:
        summary = f"Predictive model ({model.model_type}) achieved {model.accuracy".1%"} accuracy for outcome prediction."
    elif biomarkers:
        summary = f"Discovered {len(biomarkers)} potential biomarkers for further validation."
    else:
        summary = "Predictive analysis completed with limited results."

    return summary


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Predictive Modeling and Biomarker Discovery")
    parser.add_argument("--results-dir", default="results", help="Results directory")
    parser.add_argument("--format", choices=["json", "tsv", "both"], default="both",
                       help="Output format")

    args = parser.parse_args()

    summary = run_predictive_analysis(args.results_dir, args.format)
    if summary:
        print(f"\n🎯 Analysis complete: {summary}")










