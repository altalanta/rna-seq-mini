#!/usr/bin/env python3
"""
Advanced Pathway Impact Analysis for RNASEQ-MINI

This module provides sophisticated pathway impact scoring that goes beyond
simple enrichment p-values to assess functional significance and biological relevance.
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
import warnings
from dataclasses import dataclass, asdict
import networkx as nx
from scipy import stats
from collections import defaultdict

try:
    from sklearn.preprocessing import StandardScaler
    from sklearn.cluster import KMeans
    ML_AVAILABLE = True
except ImportError:
    ML_AVAILABLE = False


@dataclass
class PathwayImpact:
    """Container for pathway impact analysis results."""
    pathway_name: str
    pathway_id: str
    enrichment_score: float  # NES from fgsea
    adjusted_pvalue: float
    impact_score: float  # Our calculated impact score
    confidence_score: float
    functional_relevance: float
    expression_impact: float
    network_centrality: float
    gene_impact_distribution: Dict[str, float]
    key_drivers: List[str]
    biological_interpretation: str
    therapeutic_relevance: float = 0.0
    biomarker_potential: float = 0.0


class PathwayImpactAnalyzer:
    """Advanced analyzer for pathway impact assessment."""

    def __init__(self, results_dir: str = "results", cache_dir: str = ".pathway_cache"):
        self.results_dir = Path(results_dir)
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)

        # Load pathway databases
        self.pathway_genes = self._load_pathway_databases()

    def _load_pathway_databases(self) -> Dict[str, Dict[str, Any]]:
        """Load pathway gene set databases."""
        pathway_genes = {}

        # Load common pathway databases
        databases = {
            "kegg": "pathway_databases/kegg_pathways.json",
            "reactome": "pathway_databases/reactome_pathways.json",
            "hallmark": "pathway_databases/hallmark_pathways.json",
            "go_bp": "pathway_databases/go_biological_process.json"
        }

        for db_name, db_file in databases.items():
            db_path = Path(__file__).parent / db_file
            if db_path.exists():
                try:
                    with open(db_path) as f:
                        pathway_genes[db_name] = json.load(f)
                except Exception as e:
                    print(f"⚠️ Failed to load {db_name} database: {e}")

        return pathway_genes

    def analyze_pathway_impact(self, de_results_path: str = None,
                              fgsea_results_path: str = None) -> List[PathwayImpact]:
        """Perform comprehensive pathway impact analysis."""
        # Load input files
        if de_results_path is None:
            de_results_path = self.results_dir / "de" / "deseq2_results.tsv"
        if fgsea_results_path is None:
            fgsea_results_path = self.results_dir / "fgsea" / "fgsea_results.tsv"

        if not Path(de_results_path).exists():
            print("❌ DE results not found")
            return []

        if not Path(fgsea_results_path).exists():
            print("❌ fgsea results not found")
            return []

        try:
            de_results = pd.read_csv(de_results_path, sep='\t')
            fgsea_results = pd.read_csv(fgsea_results_path, sep='\t')
        except Exception as e:
            print(f"❌ Failed to load results: {e}")
            return []

        # Filter significant pathways
        sig_pathways = fgsea_results[fgsea_results['padj'] < 0.05].copy()

        if len(sig_pathways) == 0:
            print("⚠️ No significant pathways found")
            return []

        pathway_impacts = []

        for _, pathway_row in sig_pathways.iterrows():
            pathway_name = pathway_row['pathway']

            try:
                impact = self._calculate_pathway_impact(pathway_name, pathway_row, de_results)
                if impact:
                    pathway_impacts.append(impact)
            except Exception as e:
                print(f"⚠️ Failed to analyze pathway {pathway_name}: {e}")
                continue

        # Sort by impact score
        pathway_impacts.sort(key=lambda x: x.impact_score, reverse=True)

        print(f"✅ Analyzed {len(pathway_impacts)} pathways")
        return pathway_impacts

    def _calculate_pathway_impact(self, pathway_name: str, pathway_row: pd.Series,
                                 de_results: pd.DataFrame) -> Optional[PathwayImpact]:
        """Calculate comprehensive impact score for a pathway."""

        # Get genes in this pathway
        pathway_genes = self._get_pathway_genes(pathway_name)
        if not pathway_genes:
            return None

        # Get DE results for pathway genes
        pathway_de = de_results[de_results['gene_id'].isin(pathway_genes)].copy()

        if len(pathway_de) == 0:
            return None

        # Calculate various impact components
        functional_relevance = self._calculate_functional_relevance(pathway_genes, pathway_de)
        expression_impact = self._calculate_expression_impact(pathway_de)
        network_centrality = self._calculate_network_centrality(pathway_genes, pathway_de)

        # Calculate gene impact distribution
        gene_impacts = self._calculate_gene_impacts(pathway_de)

        # Calculate overall impact score
        impact_score = self._calculate_overall_impact(
            pathway_row['NES'], pathway_row['padj'],
            functional_relevance, expression_impact, network_centrality
        )

        # Identify key driver genes
        key_drivers = self._identify_key_drivers(pathway_de, gene_impacts)

        # Generate biological interpretation
        interpretation = self._generate_interpretation(
            pathway_name, pathway_row['NES'], functional_relevance,
            expression_impact, key_drivers
        )

        # Calculate confidence
        confidence = self._calculate_confidence(
            len(pathway_genes), len(pathway_de), pathway_row['padj']
        )

        return PathwayImpact(
            pathway_name=pathway_name,
            pathway_id=pathway_row.get('pathway_id', pathway_name),
            enrichment_score=pathway_row['NES'],
            adjusted_pvalue=pathway_row['padj'],
            impact_score=impact_score,
            confidence_score=confidence,
            functional_relevance=functional_relevance,
            expression_impact=expression_impact,
            network_centrality=network_centrality,
            gene_impact_distribution=gene_impacts,
            key_drivers=key_drivers,
            biological_interpretation=interpretation
        )

    def _get_pathway_genes(self, pathway_name: str) -> List[str]:
        """Get list of genes in a pathway."""
        # Try to find pathway in loaded databases
        for db_name, db_data in self.pathway_genes.items():
            if pathway_name in db_data:
                return db_data[pathway_name].get('genes', [])

        # Fallback: try to extract from pathway name or use placeholder
        return []

    def _calculate_functional_relevance(self, pathway_genes: List[str], pathway_de: pd.DataFrame) -> float:
        """Calculate functional relevance score based on gene annotations and interactions."""
        if len(pathway_genes) == 0:
            return 0.0

        # Calculate proportion of pathway genes that are differentially expressed
        de_genes_in_pathway = len(pathway_de)
        total_pathway_genes = len(pathway_genes)

        proportion_de = de_genes_in_pathway / total_pathway_genes if total_pathway_genes > 0 else 0

        # Bonus for having key regulatory genes
        regulatory_bonus = 0.0
        regulatory_keywords = ['kinase', 'phosphatase', 'transcription', 'receptor', 'ligand']

        for gene_id in pathway_de['gene_id']:
            if any(keyword in gene_id.lower() for keyword in regulatory_keywords):
                regulatory_bonus += 0.1

        relevance_score = min(proportion_de + regulatory_bonus, 1.0)
        return relevance_score

    def _calculate_expression_impact(self, pathway_de: pd.DataFrame) -> float:
        """Calculate expression impact based on effect sizes and consistency."""
        if len(pathway_de) == 0:
            return 0.0

        # Calculate average absolute effect size
        avg_effect = pathway_de['log2FoldChange'].abs().mean()

        # Calculate effect size consistency (lower variance = higher consistency)
        effect_variance = pathway_de['log2FoldChange'].var()
        max_possible_variance = (pathway_de['log2FoldChange'].max() - pathway_de['log2FoldChange'].min()) ** 2 / 4
        consistency = 1 - (effect_variance / max_possible_variance) if max_possible_variance > 0 else 0

        # Combine effect size and consistency
        impact_score = (avg_effect * 0.7) + (consistency * 0.3)

        return min(impact_score / 5, 1.0)  # Normalize to 0-1 scale

    def _calculate_network_centrality(self, pathway_genes: List[str], pathway_de: pd.DataFrame) -> float:
        """Calculate network centrality based on gene interactions and connectivity."""
        if len(pathway_genes) == 0:
            return 0.0

        # Simple centrality measure based on expression correlation
        if len(pathway_de) > 5:
            # Calculate how interconnected the DE genes are
            # This is a simplified measure - in practice would use PPI networks
            centrality_score = min(len(pathway_de) / 20, 1.0)  # Simple scaling
        else:
            centrality_score = 0.0

        return centrality_score

    def _calculate_gene_impacts(self, pathway_de: pd.DataFrame) -> Dict[str, float]:
        """Calculate individual gene impact scores within the pathway."""
        if len(pathway_de) == 0:
            return {}

        gene_impacts = {}

        for _, gene in pathway_de.iterrows():
            gene_id = gene['gene_id']

            # Combine effect size and significance
            effect_component = abs(gene['log2FoldChange']) / 5  # Normalize
            significance_component = 1 - gene['padj'] if gene['padj'] > 0 else 1

            # Expression level component
            expression_level = min(gene['baseMean'] / 10000, 1.0)  # Normalize

            # Combined score
            impact = (effect_component * 0.4) + (significance_component * 0.4) + (expression_level * 0.2)
            gene_impacts[gene_id] = min(impact, 1.0)

        return gene_impacts

    def _calculate_overall_impact(self, nes: float, padj: float,
                                 functional_relevance: float, expression_impact: float,
                                 network_centrality: float) -> float:
        """Calculate overall pathway impact score."""

        # Weight components
        weights = {
            'enrichment': 0.3,  # NES from fgsea
            'significance': 0.2,  # Adjusted p-value
            'functional': 0.25,  # Functional relevance
            'expression': 0.15,  # Expression impact
            'network': 0.1  # Network centrality
        }

        # Normalize NES (assuming typical range -3 to 3)
        normalized_nes = max(0, nes) / 3

        # Convert p-value to score (lower p-value = higher score)
        pvalue_score = 1 - padj if padj > 0 else 1

        # Calculate weighted score
        impact_score = (
            weights['enrichment'] * normalized_nes +
            weights['significance'] * pvalue_score +
            weights['functional'] * functional_relevance +
            weights['expression'] * expression_impact +
            weights['network'] * network_centrality
        )

        return min(impact_score, 1.0)

    def _identify_key_drivers(self, pathway_de: pd.DataFrame, gene_impacts: Dict[str, float]) -> List[str]:
        """Identify key driver genes within the pathway."""
        if len(pathway_de) == 0:
            return []

        # Get top genes by impact score
        top_genes = sorted(gene_impacts.items(), key=lambda x: x[1], reverse=True)

        # Return top 5 or top 25% of genes, whichever is smaller
        n_genes = min(5, len(top_genes))
        if len(top_genes) > 4:
            n_genes = max(n_genes, int(len(top_genes) * 0.25))

        return [gene_id for gene_id, _ in top_genes[:n_genes]]

    def _generate_interpretation(self, pathway_name: str, nes: float, functional_relevance: float,
                                expression_impact: float, key_drivers: List[str]) -> str:
        """Generate biological interpretation of pathway impact."""

        # Determine impact level
        if nes > 2:
            impact_level = "very strong"
        elif nes > 1.5:
            impact_level = "strong"
        elif nes > 1:
            impact_level = "moderate"
        else:
            impact_level = "weak"

        # Determine consistency
        if functional_relevance > 0.7:
            consistency = "highly consistent"
        elif functional_relevance > 0.5:
            consistency = "moderately consistent"
        else:
            consistency = "variable"

        # Generate interpretation
        interpretation = (
            f"This pathway shows {impact_level} enrichment (NES: {nes".2f"}) with "
            f"{consistency} functional representation of differentially expressed genes. "
        )

        if key_drivers:
            interpretation += f"Key drivers include: {', '.join(key_drivers[:3])}."

        return interpretation

    def _calculate_confidence(self, total_genes: int, de_genes: int, padj: float) -> float:
        """Calculate confidence score for the pathway impact assessment."""

        # Base confidence on multiple factors
        confidence_components = []

        # Confidence increases with more genes in pathway
        gene_ratio = de_genes / total_genes if total_genes > 0 else 0
        confidence_components.append(min(gene_ratio * 2, 1.0))

        # Confidence increases with better p-value
        pvalue_confidence = 1 - padj if padj > 0 else 1
        confidence_components.append(pvalue_confidence)

        # Confidence decreases with very small or very large pathways
        if total_genes < 10:
            size_penalty = 0.8
        elif total_genes > 500:
            size_penalty = 0.9
        else:
            size_penalty = 1.0

        # Combine components
        overall_confidence = np.mean(confidence_components) * size_penalty

        return min(overall_confidence, 1.0)

    def export_impact_analysis(self, pathway_impacts: List[PathwayImpact],
                              format: str = "json") -> str:
        """Export pathway impact analysis results."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"pathway_impact_analysis_{timestamp}"

        if format == "json":
            export_data = {
                "generated_at": datetime.now().isoformat(),
                "n_pathways": len(pathway_impacts),
                "pathways": [asdict(impact) for impact in pathway_impacts[:50]]  # Top 50
            }

            output_file = self.cache_dir / f"{filename}.json"
            with open(output_file, 'w') as f:
                json.dump(export_data, f, indent=2, default=str)

        elif format == "tsv":
            # Create detailed TSV export
            rows = []
            for impact in pathway_impacts:
                rows.append({
                    'pathway_name': impact.pathway_name,
                    'enrichment_score': impact.enrichment_score,
                    'adjusted_pvalue': impact.adjusted_pvalue,
                    'impact_score': impact.impact_score,
                    'confidence': impact.confidence_score,
                    'functional_relevance': impact.functional_relevance,
                    'expression_impact': impact.expression_impact,
                    'network_centrality': impact.network_centrality,
                    'n_de_genes': len(impact.gene_impact_distribution),
                    'key_drivers': ';'.join(impact.key_drivers[:5]),
                    'interpretation': impact.biological_interpretation
                })

            df = pd.DataFrame(rows)
            output_file = self.cache_dir / f"{filename}.tsv"
            df.to_csv(output_file, sep='\t', index=False)

        return str(output_file)


def run_pathway_impact_analysis(results_dir: str = "results", output_format: str = "both") -> str:
    """Run comprehensive pathway impact analysis."""
    analyzer = PathwayImpactAnalyzer(results_dir)

    print("🔬 Running advanced pathway impact analysis...")

    # Perform analysis
    pathway_impacts = analyzer.analyze_pathway_impact()

    if not pathway_impacts:
        print("❌ No pathway impacts to analyze")
        return None

    print(f"✅ Analyzed {len(pathway_impacts)} pathways")

    # Show top insights
    top_pathways = pathway_impacts[:5]
    print("
🏆 Top Pathways by Impact Score:"    for i, impact in enumerate(top_pathways, 1):
        print(f"{i}. {impact.pathway_name} (Impact: {impact.impact_score".3f"})")
        print(f"   NES: {impact.enrichment_score".2f"}, Confidence: {impact.confidence_score".1%"}")

    # Export results
    if output_format in ["json", "both"]:
        json_file = analyzer.export_impact_analysis(pathway_impacts, "json")
        print(f"📄 Exported JSON results to: {json_file}")

    if output_format in ["tsv", "both"]:
        tsv_file = analyzer.export_impact_analysis(pathway_impacts, "tsv")
        print(f"📄 Exported TSV results to: {tsv_file}")

    # Generate summary
    summary = f"Analyzed {len(pathway_impacts)} pathways with impact scores ranging from {pathway_impacts[-1].impact_score".3f"} to {pathway_impacts[0].impact_score".3f"}."

    return summary


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Advanced Pathway Impact Analysis")
    parser.add_argument("--results-dir", default="results", help="Results directory")
    parser.add_argument("--format", choices=["json", "tsv", "both"], default="both",
                       help="Output format")

    args = parser.parse_args()

    summary = run_pathway_impact_analysis(args.results_dir, args.format)
    if summary:
        print(f"\n🎯 Analysis complete: {summary}")











