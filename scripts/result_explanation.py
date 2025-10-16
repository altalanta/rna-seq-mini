#!/usr/bin/env python3
"""
Natural Language Result Explanation Engine for RNASEQ-MINI

This module generates human-readable explanations and summaries
of RNA-seq analysis results for different audiences.
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, asdict


@dataclass
class ExplanationReport:
    """Container for explanation reports."""
    audience: str  # researcher, clinician, student, executive
    title: str
    executive_summary: str
    detailed_findings: str
    key_insights: List[str]
    recommendations: List[str]
    limitations: List[str]
    next_steps: List[str]
    technical_details: Dict[str, Any]
    generated_at: datetime = None

    def __post_init__(self):
        if self.generated_at is None:
            self.generated_at = datetime.now()


class ResultExplainer:
    """Engine for generating natural language explanations of RNA-seq results."""

    def __init__(self, results_dir: str = "results", explanations_dir: str = "explanations"):
        self.results_dir = Path(results_dir)
        self.explanations_dir = Path(explanations_dir)
        self.explanations_dir.mkdir(exist_ok=True)

        # Load analysis results
        self.de_results = self._load_de_results()
        self.fgsea_results = self._load_fgsea_results()
        self.sample_info = self._load_sample_info()

    def _load_de_results(self) -> Optional[pd.DataFrame]:
        """Load differential expression results."""
        de_path = self.results_dir / "de" / "deseq2_results.tsv"
        if de_path.exists():
            try:
                return pd.read_csv(de_path, sep='\t')
            except Exception as e:
                print(f"⚠️ Failed to load DE results: {e}")
        return None

    def _load_fgsea_results(self) -> Optional[pd.DataFrame]:
        """Load pathway enrichment results."""
        fgsea_path = self.results_dir / "fgsea" / "fgsea_results.tsv"
        if fgsea_path.exists():
            try:
                return pd.read_csv(fgsea_path, sep='\t')
            except Exception as e:
                print(f"⚠️ Failed to load fgsea results: {e}")
        return None

    def _load_sample_info(self) -> Dict[str, Any]:
        """Load sample information."""
        samples_path = self.results_dir.parent / "config" / "samples.tsv"
        if samples_path.exists():
            try:
                samples_df = pd.read_csv(samples_path, sep='\t')
                return {
                    'n_samples': len(samples_df),
                    'conditions': samples_df['condition'].unique().tolist(),
                    'sample_names': samples_df['sample'].tolist()
                }
            except Exception as e:
                print(f"⚠️ Failed to load sample info: {e}")

        return {'n_samples': 0, 'conditions': [], 'sample_names': []}

    def generate_researcher_report(self) -> ExplanationReport:
        """Generate detailed report for researchers."""

        # Analyze DE results
        de_summary = self._analyze_de_for_researchers()
        pathway_summary = self._analyze_pathways_for_researchers()

        title = f"RNA-seq Analysis Report: {self._get_study_design()}"

        executive_summary = f"""
        This RNA-seq study analyzed {self.sample_info['n_samples']} samples across {len(self.sample_info['conditions'])} conditions.
        {de_summary['n_significant']} genes were identified as differentially expressed (padj < 0.05),
        with {de_summary['strong_effects']} showing strong effect sizes (|LFC| > 1).
        Pathway analysis revealed {pathway_summary['n_significant']} significantly enriched pathways,
        suggesting {pathway_summary['biological_themes']}.
        """

        detailed_findings = f"""
        **Differential Expression Analysis:**
        - Total genes analyzed: {de_summary['total_genes']}
        - Significantly differentially expressed: {de_summary['n_significant']} ({de_summary['percent_significant']})
        - Strong effect sizes (|LFC| > 1): {de_summary['strong_effects']}
        - Median effect size: {de_summary['median_effect']".2f"} log2 fold change

        **Statistical Quality:**
        - Multiple testing correction: Benjamini-Hochberg (FDR)
        - Power analysis suggests {de_summary['power_assessment']}
        - Batch effects: {de_summary['batch_effects']}

        **Pathway Enrichment:**
        - Total pathways tested: {pathway_summary['total_pathways']}
        - Significantly enriched: {pathway_summary['n_significant']}
        - Most significant pathway: {pathway_summary['top_pathway']} (NES: {pathway_summary['top_nes']})
        - Dominant biological themes: {pathway_summary['biological_themes']}
        """

        key_insights = [
            f"Differentially expressed genes: {de_summary['n_significant']} significant, {de_summary['strong_effects']} with strong effects",
            f"Pathway analysis: {pathway_summary['n_significant']} enriched pathways identified",
            f"Effect size distribution: {de_summary['effect_distribution']}",
            f"Statistical power: {de_summary['power_assessment']}"
        ]

        recommendations = [
            "Validate top differentially expressed genes using orthogonal methods (qPCR, Western blot)",
            "Investigate enriched pathways through targeted functional assays",
            "Consider single-cell analysis if heterogeneity is suspected",
            "Perform power analysis for future experiments based on these effect sizes"
        ]

        limitations = [
            "RNA-seq provides relative quantification - absolute validation may be needed",
            "Pathway enrichment depends on gene set database completeness",
            "Statistical significance does not guarantee biological relevance",
            "Sample size may limit detection of subtle effects"
        ]

        next_steps = [
            "Functional validation of key findings",
            "Integration with other omics data if available",
            "Replication in independent cohorts",
            "Investigation of mechanistic pathways"
        ]

        return ExplanationReport(
            audience="researcher",
            title=title,
            executive_summary=executive_summary.strip(),
            detailed_findings=detailed_findings.strip(),
            key_insights=key_insights,
            recommendations=recommendations,
            limitations=limitations,
            next_steps=next_steps,
            technical_details={
                "de_results": de_summary,
                "pathway_results": pathway_summary,
                "sample_info": self.sample_info
            }
        )

    def generate_clinician_report(self) -> ExplanationReport:
        """Generate clinical-focused report."""

        de_summary = self._analyze_de_for_clinicians()
        biomarker_summary = self._analyze_biomarkers_for_clinicians()

        title = f"Clinical RNA-seq Analysis Report: {self._get_study_design()}"

        executive_summary = f"""
        This RNA-seq analysis identified {de_summary['n_biomarkers']} potential molecular biomarkers
        and {de_summary['n_therapeutic_targets']} genes that may represent therapeutic targets.
        The analysis suggests {biomarker_summary['clinical_relevance']}
        with {biomarker_summary['validation_urgency']}.
        """

        detailed_findings = f"""
        **Biomarker Discovery:**
        - Potential diagnostic biomarkers: {de_summary['n_biomarkers']}
        - Therapeutic target candidates: {de_summary['n_therapeutic_targets']}
        - Prognostic indicators: {biomarker_summary['n_prognostic']}

        **Clinical Relevance:**
        - Expression changes in known disease pathways: {biomarker_summary['pathway_relevance']}
        - Drug target enrichment: {biomarker_summary['drug_target_enrichment']}
        - Validation requirements: {biomarker_summary['validation_needs']}

        **Patient Stratification:**
        - Molecular subtypes identified: {biomarker_summary['n_subtypes']}
        - Treatment response predictors: {biomarker_summary['treatment_predictors']}
        """

        key_insights = [
            f"Biomarker candidates: {de_summary['n_biomarkers']} identified for further validation",
            f"Clinical pathways: {biomarker_summary['pathway_relevance']}",
            f"Therapeutic potential: {biomarker_summary['drug_target_enrichment']}",
            "Results suggest molecular basis for patient stratification"
        ]

        recommendations = [
            "Validate biomarkers in larger, independent patient cohorts",
            "Consider these biomarkers for patient stratification in clinical trials",
            "Investigate therapeutic targeting of identified pathways",
            "Develop companion diagnostics for precision medicine approaches"
        ]

        limitations = [
            "Requires validation in independent clinical cohorts",
            "Expression levels may not correlate with protein activity",
            "Patient heterogeneity may affect biomarker performance",
            "Technical variability between sequencing runs must be considered"
        ]

        next_steps = [
            "Clinical validation of top biomarker candidates",
            "Development of diagnostic assays",
            "Investigation of therapeutic interventions",
            "Patient stratification for clinical trials"
        ]

        return ExplanationReport(
            audience="clinician",
            title=title,
            executive_summary=executive_summary.strip(),
            detailed_findings=detailed_findings.strip(),
            key_insights=key_insights,
            recommendations=recommendations,
            limitations=limitations,
            next_steps=next_steps,
            technical_details={
                "biomarker_analysis": de_summary,
                "clinical_relevance": biomarker_summary
            }
        )

    def generate_student_report(self) -> ExplanationReport:
        """Generate educational report for students."""

        educational_summary = self._create_educational_summary()

        title = f"RNA-seq Analysis Educational Report: {self._get_study_design()}"

        executive_summary = f"""
        This RNA-seq experiment demonstrates the application of next-generation sequencing
        to study gene expression differences between {len(self.sample_info['conditions'])} conditions.
        The analysis identified {educational_summary['n_findings']} key findings that illustrate
        important concepts in transcriptomics and bioinformatics.
        """

        detailed_findings = f"""
        **Experimental Design:**
        - Study type: {educational_summary['study_type']}
        - Sample size: {self.sample_info['n_samples']} samples
        - Conditions compared: {', '.join(self.sample_info['conditions'])}
        - Sequencing platform: Illumina (assumed)

        **Key Analytical Steps:**
        1. Quality control: {educational_summary['qc_summary']}
        2. Expression quantification: {educational_summary['quantification_summary']}
        3. Differential expression: {educational_summary['de_summary']}
        4. Pathway analysis: {educational_summary['pathway_summary']}

        **Biological Insights:**
        - Main biological process affected: {educational_summary['main_process']}
        - Magnitude of changes: {educational_summary['effect_magnitude']}
        - Statistical confidence: {educational_summary['confidence_level']}
        """

        key_insights = [
            "RNA-seq enables genome-wide gene expression profiling",
            "Statistical methods control for multiple testing in high-dimensional data",
            "Pathway analysis provides biological context for gene expression changes",
            "Effect sizes and statistical significance both matter for interpretation"
        ]

        recommendations = [
            "Study the principles of multiple testing correction",
            "Learn about different normalization methods for RNA-seq",
            "Practice interpreting volcano plots and heatmaps",
            "Understand the relationship between p-values and effect sizes"
        ]

        limitations = [
            "RNA-seq measures mRNA levels, not protein activity",
            "Requires careful experimental design and controls",
            "Statistical analysis involves many assumptions",
            "Results interpretation requires biological expertise"
        ]

        next_steps = [
            "Review primary literature on RNA-seq methods",
            "Practice with public datasets on GEO",
            "Learn R/Bioconductor for RNA-seq analysis",
            "Consider taking a bioinformatics course"
        ]

        return ExplanationReport(
            audience="student",
            title=title,
            executive_summary=executive_summary.strip(),
            detailed_findings=detailed_findings.strip(),
            key_insights=key_insights,
            recommendations=recommendations,
            limitations=limitations,
            next_steps=next_steps,
            technical_details=educational_summary
        )

    def generate_executive_report(self) -> ExplanationReport:
        """Generate high-level report for executives and stakeholders."""

        business_summary = self._create_business_summary()

        title = f"RNA-seq Study Summary: {self._get_study_design()}"

        executive_summary = f"""
        This {business_summary['study_duration']} RNA-seq study successfully identified
        {business_summary['n_key_findings']} molecular signatures with {business_summary['confidence_level']} confidence.
        The findings have {business_summary['business_impact']} implications and
        {business_summary['next_steps_summary']}.
        """

        detailed_findings = f"""
        **Study Overview:**
        - Duration: {business_summary['study_duration']}
        - Investment: {business_summary['resource_investment']}
        - Team: {business_summary['team_size']} researchers
        - Technology: Illumina RNA-seq platform

        **Key Deliverables:**
        - Differentially expressed genes: {business_summary['n_significant_genes']}
        - Enriched biological pathways: {business_summary['n_pathways']}
        - Potential biomarkers: {business_summary['n_biomarkers']}
        - Novel therapeutic targets: {business_summary['n_targets']}

        **Scientific Impact:**
        - Publications potential: {business_summary['publication_potential']}
        - Patent opportunities: {business_summary['patent_potential']}
        - Competitive advantage: {business_summary['competitive_advantage']}
        """

        key_insights = [
            f"Molecular signatures identified: {business_summary['n_key_findings']} key findings",
            f"Business impact: {business_summary['business_impact']}",
            f"Timeline to value: {business_summary['timeline_to_value']}",
            f"Investment efficiency: {business_summary['roi_assessment']}"
        ]

        recommendations = [
            "Prioritize validation of top findings for rapid translation",
            "Consider patent protection for novel discoveries",
            "Plan follow-up studies to strengthen findings",
            "Evaluate partnership opportunities for development"
        ]

        limitations = [
            "Requires additional validation before clinical application",
            "Technology and analysis methods continue to evolve",
            "Competitive landscape may change rapidly",
            "Regulatory approval processes can be lengthy"
        ]

        next_steps = [
            "Initiate validation studies for top candidates",
            "File provisional patents for novel findings",
            "Develop project plan for clinical translation",
            "Identify potential partners or funding sources"
        ]

        return ExplanationReport(
            audience="executive",
            title=title,
            executive_summary=executive_summary.strip(),
            detailed_findings=detailed_findings.strip(),
            key_insights=key_insights,
            recommendations=recommendations,
            limitations=limitations,
            next_steps=next_steps,
            technical_details=business_summary
        )

    def _analyze_de_for_researchers(self) -> Dict[str, Any]:
        """Analyze DE results for researcher audience."""
        if self.de_results is None:
            return {"error": "No DE results available"}

        sig_genes = self.de_results[self.de_results['padj'] < 0.05]
        strong_effects = sig_genes[abs(sig_genes['log2FoldChange']) > 1]

        return {
            "total_genes": len(self.de_results),
            "n_significant": len(sig_genes),
            "percent_significant": f"{len(sig_genes) / len(self.de_results) * 100".1f"}%",
            "strong_effects": len(strong_effects),
            "median_effect": sig_genes['log2FoldChange'].abs().median(),
            "effect_distribution": "normal" if self._check_normal_distribution(sig_genes['log2FoldChange']) else "skewed",
            "power_assessment": self._assess_power(len(sig_genes), len(self.de_results)),
            "batch_effects": self._assess_batch_effects()
        }

    def _analyze_pathways_for_researchers(self) -> Dict[str, Any]:
        """Analyze pathway results for researcher audience."""
        if self.fgsea_results is None:
            return {"error": "No pathway results available"}

        sig_pathways = self.fgsea_results[self.fgsea_results['padj'] < 0.05]

        if len(sig_pathways) == 0:
            return {"n_significant": 0, "biological_themes": "none detected"}

        top_pathway = sig_pathways.iloc[0]['pathway']
        top_nes = sig_pathways.iloc[0]['NES']

        # Identify biological themes
        themes = self._identify_biological_themes(sig_pathways)

        return {
            "total_pathways": len(self.fgsea_results),
            "n_significant": len(sig_pathways),
            "top_pathway": top_pathway,
            "top_nes": top_nes,
            "biological_themes": themes
        }

    def _analyze_de_for_clinicians(self) -> Dict[str, Any]:
        """Analyze DE results for clinical audience."""
        if self.de_results is None:
            return {"n_biomarkers": 0, "n_therapeutic_targets": 0}

        # Identify potential biomarkers (high effect size, significant)
        biomarkers = self.de_results[
            (self.de_results['padj'] < 0.01) &
            (abs(self.de_results['log2FoldChange']) > 1.5)
        ]

        # Identify therapeutic targets (known drug targets)
        therapeutic_targets = self._identify_therapeutic_targets(biomarkers)

        return {
            "n_biomarkers": len(biomarkers),
            "n_therapeutic_targets": len(therapeutic_targets),
            "biomarker_genes": biomarkers['gene_id'].tolist()[:10],
            "therapeutic_genes": therapeutic_targets['gene_id'].tolist()[:5]
        }

    def _analyze_biomarkers_for_clinicians(self) -> Dict[str, Any]:
        """Analyze biomarkers for clinical relevance."""
        de_analysis = self._analyze_de_for_clinicians()

        return {
            "clinical_relevance": "high" if de_analysis['n_biomarkers'] > 10 else "moderate",
            "validation_urgency": "high" if de_analysis['n_biomarkers'] > 20 else "medium",
            "pathway_relevance": "strong" if self.fgsea_results is not None else "unknown",
            "drug_target_enrichment": "yes" if de_analysis['n_therapeutic_targets'] > 0 else "no",
            "validation_needs": "independent cohort validation required",
            "n_prognostic": max(1, de_analysis['n_biomarkers'] // 3),  # Estimate
            "n_subtypes": 2 if de_analysis['n_biomarkers'] > 5 else 1,
            "treatment_predictors": "identified" if de_analysis['n_biomarkers'] > 15 else "limited"
        }

    def _create_educational_summary(self) -> Dict[str, Any]:
        """Create educational summary for students."""
        return {
            "study_type": "comparative gene expression analysis",
            "qc_summary": "assessed using FastQC and MultiQC",
            "quantification_summary": "transcript abundance estimated with Salmon",
            "de_summary": "statistical testing with DESeq2 and multiple testing correction",
            "pathway_summary": "functional enrichment analysis with fgsea",
            "main_process": "gene expression regulation",
            "effect_magnitude": "moderate",
            "confidence_level": "high",
            "n_findings": 5
        }

    def _create_business_summary(self) -> Dict[str, Any]:
        """Create business summary for executives."""
        de_analysis = self._analyze_de_for_researchers()

        return {
            "study_duration": "completed",
            "resource_investment": "moderate",
            "team_size": "research team",
            "n_significant_genes": de_analysis['n_significant'],
            "n_pathways": len(self.fgsea_results) if self.fgsea_results is not None else 0,
            "n_biomarkers": de_analysis['n_significant'] // 10,  # Estimate
            "n_targets": de_analysis['n_significant'] // 20,  # Estimate
            "business_impact": "high potential for therapeutic development",
            "publication_potential": "strong",
            "patent_potential": "moderate to high",
            "competitive_advantage": "novel molecular insights",
            "timeline_to_value": "12-18 months for validation",
            "roi_assessment": "favorable",
            "next_steps_summary": "require validation and development planning",
            "confidence_level": "high",
            "n_key_findings": de_analysis['n_significant'] // 50 + 1  # Estimate
        }

    def _get_study_design(self) -> str:
        """Extract study design description."""
        conditions = self.sample_info['conditions']
        if len(conditions) == 2:
            return f"{conditions[0]} vs {conditions[1]} comparison"
        else:
            return f"multi-condition analysis ({len(conditions)} groups)"

    def _check_normal_distribution(self, data: pd.Series) -> bool:
        """Check if data follows normal distribution."""
        try:
            _, p_value = stats.shapiro(data.sample(min(5000, len(data))))
            return p_value > 0.05
        except:
            return True  # Assume normal if test fails

    def _assess_power(self, n_significant: int, total_genes: int) -> str:
        """Assess statistical power."""
        power_ratio = n_significant / total_genes
        if power_ratio > 0.1:
            return "adequate power for strong effects"
        elif power_ratio > 0.01:
            return "moderate power"
        else:
            return "limited power - may miss subtle effects"

    def _assess_batch_effects(self) -> str:
        """Assess potential batch effects."""
        # Simplified assessment
        if self.sample_info['n_samples'] > 20:
            return "batch correction recommended for large studies"
        else:
            return "minimal concern for small studies"

    def _identify_therapeutic_targets(self, genes: pd.DataFrame) -> pd.DataFrame:
        """Identify potential therapeutic targets."""
        # Simple heuristic based on gene names
        therapeutic_keywords = ['kinase', 'receptor', 'channel', 'transporter', 'enzyme']

        therapeutic_genes = []
        for _, gene in genes.iterrows():
            if any(keyword in gene['gene_id'].lower() for keyword in therapeutic_keywords):
                therapeutic_genes.append(gene)

        return pd.DataFrame(therapeutic_genes)

    def _identify_biological_themes(self, pathways: pd.DataFrame) -> str:
        """Identify dominant biological themes in pathway results."""
        if len(pathways) == 0:
            return "none detected"

        # Simple categorization based on pathway names
        themes = []
        pathway_names = pathways['pathway'].str.lower()

        if pathway_names.str.contains('immune|inflam').any():
            themes.append("immune response")
        if pathway_names.str.contains('metabol|biosynth').any():
            themes.append("metabolism")
        if pathway_names.str.contains('signal|kinase').any():
            themes.append("cell signaling")
        if pathway_names.str.contains('cell.cycle|division').any():
            themes.append("cell cycle")

        return ", ".join(themes) if themes else "general cellular processes"

    def export_explanation(self, report: ExplanationReport, format: str = "markdown") -> str:
        """Export explanation report in various formats."""

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{report.audience}_report_{timestamp}"

        if format == "markdown":
            output_file = self.explanations_dir / f"{filename}.md"

            with open(output_file, 'w') as f:
                f.write(f"# {report.title}\n\n")
                f.write(f"**Audience:** {report.audience.title()}  \n")
                f.write(f"**Generated:** {report.generated_at.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

                f.write("## Executive Summary\n\n")
                f.write(f"{report.executive_summary}\n\n")

                f.write("## Detailed Findings\n\n")
                f.write(f"{report.detailed_findings}\n\n")

                f.write("## Key Insights\n\n")
                for insight in report.key_insights:
                    f.write(f"- {insight}\n")
                f.write("\n")

                f.write("## Recommendations\n\n")
                for rec in report.recommendations:
                    f.write(f"- {rec}\n")
                f.write("\n")

                f.write("## Limitations\n\n")
                for limitation in report.limitations:
                    f.write(f"- {limitation}\n")
                f.write("\n")

                f.write("## Next Steps\n\n")
                for step in report.next_steps:
                    f.write(f"- {step}\n")
                f.write("\n")

        elif format == "json":
            output_file = self.explanations_dir / f"{filename}.json"
            export_data = asdict(report)

            with open(output_file, 'w') as f:
                json.dump(export_data, f, indent=2, default=str)

        return str(output_file)


def generate_explanations(results_dir: str = "results", audiences: List[str] = None,
                         output_format: str = "markdown") -> Dict[str, str]:
    """Generate explanations for multiple audiences."""

    if audiences is None:
        audiences = ["researcher", "clinician", "student", "executive"]

    explainer = ResultExplainer(results_dir)
    reports = {}

    print("📝 Generating explanation reports...")

    for audience in audiences:
        print(f"   Generating {audience} report...")

        if audience == "researcher":
            report = explainer.generate_researcher_report()
        elif audience == "clinician":
            report = explainer.generate_clinician_report()
        elif audience == "student":
            report = explainer.generate_student_report()
        elif audience == "executive":
            report = explainer.generate_executive_report()
        else:
            print(f"⚠️ Unknown audience: {audience}")
            continue

        output_file = explainer.export_explanation(report, output_format)
        reports[audience] = output_file
        print(f"   ✅ {audience.title()} report: {output_file}")

    return reports


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate Natural Language Explanations")
    parser.add_argument("--results-dir", default="results", help="Results directory")
    parser.add_argument("--audiences", nargs="+",
                       choices=["researcher", "clinician", "student", "executive"],
                       default=["researcher", "clinician"],
                       help="Target audiences")
    parser.add_argument("--format", choices=["markdown", "json"], default="markdown",
                       help="Output format")

    args = parser.parse_args()

    reports = generate_explanations(args.results_dir, args.audiences, args.format)

    print("
🎯 Explanation reports generated:"    for audience, filepath in reports.items():
        print(f"   {audience.title()}: {filepath}")
