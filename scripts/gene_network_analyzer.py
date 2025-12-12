#!/usr/bin/env python3
"""
Gene Interaction Network Analysis for RNASEQ-MINI

This module analyzes gene-gene interactions and regulatory networks
to identify key regulatory patterns and biological modules.
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
import networkx as nx
from scipy import stats
from collections import defaultdict
import matplotlib.pyplot as plt

try:
    import plotly.graph_objects as go
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False


@dataclass
class NetworkModule:
    """Container for gene network modules."""
    module_id: str
    genes: List[str]
    hub_genes: List[str]
    module_score: float
    functional_enrichment: Dict[str, float]
    regulatory_patterns: List[str]
    biological_interpretation: str


@dataclass
class GeneInteraction:
    """Container for gene-gene interactions."""
    gene1: str
    gene2: str
    interaction_type: str  # coexpression, regulatory, protein_protein, etc.
    strength: float
    direction: str  # positive, negative, bidirectional
    confidence: float


class GeneNetworkAnalyzer:
    """Analyzer for gene interaction networks and modules."""

    def __init__(self, results_dir: str = "results", cache_dir: str = ".network_cache"):
        self.results_dir = Path(results_dir)
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)

        # Load interaction databases
        self.interaction_db = self._load_interaction_databases()

    def _load_interaction_databases(self) -> Dict[str, Any]:
        """Load gene interaction databases."""
        interaction_db = {}

        # Try to load common interaction databases
        db_files = {
            "string_ppi": "interaction_databases/string_ppi.json",
            "regulatory": "interaction_databases/regulatory_interactions.json",
            "pathway_interactions": "interaction_databases/pathway_interactions.json"
        }

        for db_name, db_file in db_files.items():
            db_path = Path(__file__).parent / db_file
            if db_path.exists():
                try:
                    with open(db_path) as f:
                        interaction_db[db_name] = json.load(f)
                except Exception as e:
                    print(f"⚠️ Failed to load {db_name}: {e}")

        return interaction_db

    def analyze_gene_networks(self, de_results_path: str = None,
                             expression_matrix_path: str = None) -> Tuple[List[NetworkModule], List[GeneInteraction]]:
        """Analyze gene interaction networks from DE results."""

        # Load input files
        if de_results_path is None:
            de_results_path = self.results_dir / "de" / "deseq2_results.tsv"
        if expression_matrix_path is None:
            expression_matrix_path = self.results_dir / "counts" / "counts.tsv"

        if not Path(de_results_path).exists():
            print("❌ DE results not found")
            return [], []

        try:
            de_results = pd.read_csv(de_results_path, sep='\t')
        except Exception as e:
            print(f"❌ Failed to load DE results: {e}")
            return [], []

        # Filter significant genes
        sig_genes = de_results[de_results['padj'] < 0.05].copy()

        if len(sig_genes) < 10:
            print("⚠️ Too few significant genes for network analysis")
            return [], []

        print(f"🔬 Analyzing networks for {len(sig_genes)} significant genes")

        # Build different types of networks
        networks = {}

        # 1. Co-expression network
        if Path(expression_matrix_path).exists():
            networks['coexpression'] = self._build_coexpression_network(sig_genes, expression_matrix_path)

        # 2. Regulatory network
        networks['regulatory'] = self._build_regulatory_network(sig_genes)

        # 3. Protein-protein interaction network
        networks['ppi'] = self._build_ppi_network(sig_genes)

        # Analyze network modules
        modules = self._identify_network_modules(networks)

        # Extract key interactions
        interactions = self._extract_key_interactions(networks)

        print(f"✅ Identified {len(modules)} network modules")
        return modules, interactions

    def _build_coexpression_network(self, sig_genes: pd.DataFrame,
                                   expression_matrix_path: str) -> nx.Graph:
        """Build co-expression network based on correlation."""

        try:
            # Load normalized expression data
            expression_df = pd.read_csv(expression_matrix_path, sep='\t', index_col=0)

            # Get expression for significant genes
            common_genes = [g for g in sig_genes['gene_id'] if g in expression_df.index]
            if len(common_genes) < 5:
                return nx.Graph()

            gene_expression = expression_df.loc[common_genes]

            # Calculate correlations
            correlation_matrix = gene_expression.T.corr()

            # Build network
            G = nx.Graph()

            for i, gene1 in enumerate(common_genes):
                for j, gene2 in enumerate(common_genes):
                    if i < j:  # Avoid duplicate edges
                        corr = correlation_matrix.loc[gene1, gene2]
                        if abs(corr) > 0.7:  # Strong correlation threshold
                            G.add_edge(gene1, gene2,
                                     weight=abs(corr),
                                     correlation=corr,
                                     interaction_type='coexpression')

            return G

        except Exception as e:
            print(f"⚠️ Co-expression network failed: {e}")
            return nx.Graph()

    def _build_regulatory_network(self, sig_genes: pd.DataFrame) -> nx.Graph:
        """Build regulatory interaction network."""

        G = nx.DiGraph()  # Directed graph for regulatory interactions

        # Use simplified regulatory logic based on gene names and expression patterns
        regulatory_genes = []
        target_genes = []

        # Identify potential regulators (kinases, transcription factors, etc.)
        for gene_id in sig_genes['gene_id']:
            if any(keyword in gene_id.lower() for keyword in ['kinase', 'tf', 'transcription']):
                regulatory_genes.append(gene_id)
            else:
                target_genes.append(gene_id)

        # Create regulatory edges based on effect size patterns
        for regulator in regulatory_genes:
            regulator_data = sig_genes[sig_genes['gene_id'] == regulator]
            if len(regulator_data) == 0:
                continue

            regulator_lfc = regulator_data['log2FoldChange'].iloc[0]

            for target in target_genes:
                target_data = sig_genes[sig_genes['gene_id'] == target]
                if len(target_data) == 0:
                    continue

                target_lfc = target_data['log2FoldChange'].iloc[0]

                # Simple regulatory logic: opposite direction suggests regulation
                if (regulator_lfc > 0 and target_lfc < 0) or (regulator_lfc < 0 and target_lfc > 0):
                    correlation = -1 if regulator_lfc * target_lfc < 0 else 1
                    G.add_edge(regulator, target,
                             weight=abs(regulator_lfc) * abs(target_lfc) / 10,
                             correlation=correlation,
                             interaction_type='putative_regulation')

        return G

    def _build_ppi_network(self, sig_genes: pd.DataFrame) -> nx.Graph:
        """Build protein-protein interaction network."""

        G = nx.Graph()

        # Use STRING database if available, otherwise use simple heuristics
        if "string_ppi" in self.interaction_db:
            string_data = self.interaction_db["string_ppi"]

            for gene1 in sig_genes['gene_id']:
                for gene2 in sig_genes['gene_id']:
                    if gene1 != gene2:
                        interaction_key = f"{gene1}_{gene2}"
                        if interaction_key in string_data:
                            interaction_score = string_data[interaction_key]
                            if interaction_score > 0.7:  # High confidence interactions
                                G.add_edge(gene1, gene2,
                                         weight=interaction_score,
                                         interaction_type='protein_protein')

        return G

    def _identify_network_modules(self, networks: Dict[str, nx.Graph]) -> List[NetworkModule]:
        """Identify functional modules in gene networks."""

        modules = []

        for network_type, G in networks.items():
            if G.number_of_nodes() < 5:
                continue

            try:
                # Use community detection if available
                if hasattr(nx, 'community'):
                    try:
                        communities = nx.community.greedy_modularity_communities(G)
                        community_list = list(communities)
                    except Exception:
                        # Fallback to connected components
                        community_list = list(nx.connected_components(G))
                else:
                    # Fallback to connected components
                    community_list = list(nx.connected_components(G))

                for i, community in enumerate(community_list):
                    if len(community) >= 3:  # Minimum module size
                        module = self._analyze_module(community, G, network_type, i)
                        if module:
                            modules.append(module)

            except Exception as e:
                print(f"⚠️ Module identification failed for {network_type}: {e}")

        # Sort by module score
        modules.sort(key=lambda x: x.module_score, reverse=True)

        return modules

    def _analyze_module(self, genes: List[str], G: nx.Graph,
                       network_type: str, module_id: int) -> Optional[NetworkModule]:
        """Analyze a single network module."""

        try:
            # Extract subgraph for this module
            module_subgraph = G.subgraph(genes)

            # Calculate module score based on connectivity
            density = nx.density(module_subgraph)
            avg_degree = sum(dict(module_subgraph.degree()).values()) / len(genes)

            # Simple module score
            module_score = (density + avg_degree / 10) / 2

            # Identify hub genes (high degree)
            degrees = dict(module_subgraph.degree())
            hub_genes = [gene for gene, degree in degrees.items() if degree > avg_degree]

            # Generate functional interpretation
            interpretation = self._interpret_module_function(genes, network_type)

            # Identify regulatory patterns
            patterns = self._identify_regulatory_patterns(module_subgraph)

            return NetworkModule(
                module_id=f"{network_type}_module_{module_id}",
                genes=list(genes),
                hub_genes=hub_genes,
                module_score=module_score,
                functional_enrichment={},  # Would implement pathway enrichment
                regulatory_patterns=patterns,
                biological_interpretation=interpretation
            )

        except Exception as e:
            print(f"⚠️ Module analysis failed: {e}")
            return None

    def _interpret_module_function(self, genes: List[str], network_type: str) -> str:
        """Generate biological interpretation of module function."""

        # Simple interpretation based on network type and gene patterns
        if network_type == "coexpression":
            interpretation = "Genes in this module show coordinated expression patterns, suggesting shared regulatory control."
        elif network_type == "regulatory":
            interpretation = "This module contains regulatory relationships, indicating potential signaling cascades."
        elif network_type == "ppi":
            interpretation = "Protein-protein interactions suggest functional complexes or pathways."
        else:
            interpretation = "This module represents a functionally related group of genes."

        # Add size information
        interpretation += f" The module contains {len(genes)} genes."

        return interpretation

    def _identify_regulatory_patterns(self, subgraph: nx.Graph) -> List[str]:
        """Identify regulatory patterns in the module."""

        patterns = []

        # Check for hub-and-spoke pattern
        degrees = dict(subgraph.degree())
        if degrees:
            max_degree = max(degrees.values())
            avg_degree = sum(degrees.values()) / len(degrees)

            if max_degree > avg_degree * 2:
                hub_gene = max(degrees, key=degrees.get)
                patterns.append(f"Hub-and-spoke pattern with {hub_gene} as central regulator")

        # Check for cascade patterns (would need directed graph for this)
        if subgraph.is_directed():
            # Analyze for cascade patterns
            pass

        if not patterns:
            patterns.append("General connectivity pattern")

        return patterns

    def _extract_key_interactions(self, networks: Dict[str, nx.Graph]) -> List[GeneInteraction]:
        """Extract key interactions from all networks."""

        interactions = []

        for network_type, G in networks.items():
            for edge in G.edges(data=True):
                gene1, gene2, edge_data = edge

                interaction = GeneInteraction(
                    gene1=gene1,
                    gene2=gene2,
                    interaction_type=edge_data.get('interaction_type', network_type),
                    strength=edge_data.get('weight', 1.0),
                    direction="bidirectional",  # Default
                    confidence=0.8  # Default confidence
                )

                # Set direction based on correlation or regulatory logic
                if 'correlation' in edge_data:
                    correlation = edge_data['correlation']
                    if correlation > 0:
                        interaction.direction = "positive"
                    else:
                        interaction.direction = "negative"

                interactions.append(interaction)

        # Sort by strength and remove duplicates
        interactions.sort(key=lambda x: x.strength, reverse=True)

        # Remove duplicates (same gene pair with different types)
        seen_pairs = set()
        unique_interactions = []
        for interaction in interactions:
            pair = tuple(sorted([interaction.gene1, interaction.gene2]))
            if pair not in seen_pairs:
                seen_pairs.add(pair)
                unique_interactions.append(interaction)

        return unique_interactions[:100]  # Top 100 interactions

    def visualize_network(self, modules: List[NetworkModule], interactions: List[GeneInteraction],
                         output_dir: str = "network_visualizations"):
        """Create network visualizations."""

        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)

        # Create module summary visualization
        if PLOTLY_AVAILABLE:
            self._create_interactive_visualization(modules, interactions, output_path)
        else:
            self._create_static_visualization(modules, interactions, output_path)

    def _create_interactive_visualization(self, modules: List[NetworkModule],
                                        interactions: List[GeneInteraction], output_path: Path):
        """Create interactive network visualization using Plotly."""

        # Create network graph for visualization
        G = nx.Graph()

        # Add top interactions as edges
        for interaction in interactions[:50]:  # Top 50 for clarity
            G.add_edge(interaction.gene1, interaction.gene2,
                      weight=interaction.strength,
                      interaction_type=interaction.interaction_type)

        # Get positions
        pos = nx.spring_layout(G, k=1, iterations=50)

        # Create edge traces
        edge_traces = []
        for edge in G.edges(data=True):
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]

            edge_trace = go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                line=dict(width=edge[2]['weight']*3, color='rgba(100,100,100,0.5)'),
                hoverinfo='none',
                mode='lines'
            )
            edge_traces.append(edge_trace)

        # Create node traces
        node_x = []
        node_y = []
        node_info = []
        node_colors = []

        for node in G.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            node_info.append(node)

            # Color by module membership (if applicable)
            node_colors.append('blue')

        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers+text',
            hoverinfo='text',
            text=node_info,
            textposition="middle right",
            marker=dict(size=10, color=node_colors),
            textfont=dict(size=8)
        )

        # Create figure
        fig = go.Figure(data=edge_traces + [node_trace],
                       layout=go.Layout(
                           title='Gene Interaction Network',
                           titlefont_size=16,
                           showlegend=False,
                           hovermode='closest',
                           margin=dict(b=20,l=5,r=5,t=40),
                           xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                           yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
                       ))

        # Save interactive plot
        output_file = output_path / "network_interactive.html"
        fig.write_html(str(output_file))
        print(f"📊 Interactive network visualization: {output_file}")

    def _create_static_visualization(self, modules: List[NetworkModule],
                                    interactions: List[GeneInteraction], output_path: Path):
        """Create static network visualization using matplotlib."""

        # Create a simple network plot
        plt.figure(figsize=(12, 8))

        # Plot module sizes
        module_sizes = [len(module.genes) for module in modules[:10]]
        module_names = [module.module_id for module in modules[:10]]

        plt.bar(range(len(module_sizes)), module_sizes)
        plt.xticks(range(len(module_names)), module_names, rotation=45, ha='right')
        plt.title('Network Module Sizes')
        plt.ylabel('Number of Genes')
        plt.tight_layout()

        # Save plot
        output_file = output_path / "network_modules.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 Static network visualization: {output_file}")

    def export_network_analysis(self, modules: List[NetworkModule],
                               interactions: List[GeneInteraction], format: str = "json") -> str:
        """Export network analysis results."""

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"gene_network_analysis_{timestamp}"

        if format == "json":
            export_data = {
                "generated_at": datetime.now().isoformat(),
                "n_modules": len(modules),
                "n_interactions": len(interactions),
                "modules": [asdict(module) for module in modules],
                "interactions": [{
                    "gene1": interaction.gene1,
                    "gene2": interaction.gene2,
                    "interaction_type": interaction.interaction_type,
                    "strength": interaction.strength,
                    "direction": interaction.direction,
                    "confidence": interaction.confidence
                } for interaction in interactions[:100]]  # Top 100
            }

            output_file = self.cache_dir / f"{filename}.json"
            with open(output_file, 'w') as f:
                json.dump(export_data, f, indent=2, default=str)

        elif format == "tsv":
            # Export modules
            modules_data = []
            for module in modules:
                modules_data.append({
                    'module_id': module.module_id,
                    'n_genes': len(module.genes),
                    'hub_genes': ';'.join(module.hub_genes),
                    'module_score': module.module_score,
                    'interpretation': module.biological_interpretation
                })

            modules_df = pd.DataFrame(modules_data)
            modules_file = self.cache_dir / f"{filename}_modules.tsv"
            modules_df.to_csv(modules_file, sep='\t', index=False)

            # Export interactions
            interactions_data = []
            for interaction in interactions[:200]:  # Top 200
                interactions_data.append({
                    'gene1': interaction.gene1,
                    'gene2': interaction.gene2,
                    'interaction_type': interaction.interaction_type,
                    'strength': interaction.strength,
                    'direction': interaction.direction,
                    'confidence': interaction.confidence
                })

            interactions_df = pd.DataFrame(interactions_data)
            interactions_file = self.cache_dir / f"{filename}_interactions.tsv"
            interactions_df.to_csv(interactions_file, sep='\t', index=False)

            return f"modules: {modules_file}, interactions: {interactions_file}"

        return str(output_file)


def run_gene_network_analysis(results_dir: str = "results", output_format: str = "both") -> str:
    """Run gene network analysis."""
    analyzer = GeneNetworkAnalyzer(results_dir)

    print("🔬 Running gene network analysis...")

    # Perform analysis
    modules, interactions = analyzer.analyze_gene_networks()

    if not modules and not interactions:
        print("❌ No network analysis results generated")
        return None

    print(f"✅ Identified {len(modules)} network modules and {len(interactions)} key interactions")

    # Show top modules
    if modules:
        print("\n🏆 Top Network Modules:")
        for i, module in enumerate(modules[:3], 1):
            print(f"{i}. {module.module_id}: {len(module.genes)} genes (score: {module.module_score:.3f})")
            if module.hub_genes:
                print(f"   Hub genes: {', '.join(module.hub_genes[:3])}")

    # Show top interactions
    if interactions:
        print("\n🔗 Top Gene Interactions:")
        for i, interaction in enumerate(interactions[:5], 1):
            direction_symbol = "→" if interaction.direction == "positive" else "←"
            print(f"{i}. {interaction.gene1} {direction_symbol} {interaction.gene2} "
                  f"({interaction.interaction_type}, strength: {interaction.strength:.2f})")

    # Create visualizations
    analyzer.visualize_network(modules, interactions)

    # Export results
    if output_format in ["json", "both"]:
        json_file = analyzer.export_network_analysis(modules, interactions, "json")
        print(f"📄 Exported JSON results to: {json_file}")

    if output_format in ["tsv", "both"]:
        tsv_files = analyzer.export_network_analysis(modules, interactions, "tsv")
        print(f"📄 Exported TSV results to: {tsv_files}")

    # Generate summary
    summary = f"Network analysis identified {len(modules)} functional modules and {len(interactions)} key interactions among differentially expressed genes."

    return summary


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Gene Network Analysis")
    parser.add_argument("--results-dir", default="results", help="Results directory")
    parser.add_argument("--format", choices=["json", "tsv", "both"], default="both",
                       help="Output format")

    args = parser.parse_args()

    summary = run_gene_network_analysis(args.results_dir, args.format)
    if summary:
        print(f"\n🎯 Analysis complete: {summary}")












