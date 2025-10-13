#!/usr/bin/env python3
"""
Interactive web-based analysis environment for RNASEQ-MINI results.
Provides dynamic exploration of QC, DE, and pathway analysis results.
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Any
import logging
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.responses import HTMLResponse, JSONResponse
import plotly
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Initialize FastAPI app
app = FastAPI(
    title="RNASEQ-MINI Interactive Analysis",
    description="Interactive web interface for exploring RNA-seq analysis results",
    version="1.0.0"
)

# Mount static files and templates
app.mount("/static", StaticFiles(directory="web_app/static"), name="static")
templates = Jinja2Templates(directory="web_app/templates")


class ResultsExplorer:
    """Handles loading and serving RNA-seq results data."""

    def __init__(self, results_dir: str = "results"):
        self.results_dir = Path(results_dir)
        self.data_cache = {}

    def load_qc_data(self) -> Dict[str, Any]:
        """Load QC data for visualization."""
        qc_data = {}

        # Load MultiQC data if available
        multiqc_path = self.results_dir / "qc" / "multiqc" / "multiqc_data.json"
        if multiqc_path.exists():
            try:
                with open(multiqc_path, 'r') as f:
                    multiqc_data = json.load(f)

                # Extract key metrics for plotting
                qc_data['multiqc'] = {
                    'samples': list(multiqc_data.get('report_data_sources', {}).keys()),
                    'metrics': multiqc_data.get('report_general_stats_data', [])
                }
            except Exception as e:
                logger.warning(f"Error loading MultiQC data: {e}")

        # Load FastQC data
        fastqc_dir = self.results_dir / "qc" / "fastqc"
        if fastqc_dir.exists():
            qc_data['fastqc'] = {}
            for sample_dir in fastqc_dir.iterdir():
                if sample_dir.is_dir():
                    sample_name = sample_dir.name.replace('_fastqc', '')
                    qc_data['fastqc'][sample_name] = {
                        'path': str(sample_dir),
                        'exists': True
                    }

        return qc_data

    def load_differential_expression_data(self) -> Dict[str, Any]:
        """Load differential expression results."""
        de_data = {}

        # Load DE summary
        de_summary_path = self.results_dir / "de" / "de_summary.tsv"
        if de_summary_path.exists():
            try:
                de_summary = pd.read_csv(de_summary_path, sep='\t')
                de_data['summary'] = de_summary.to_dict('records')
            except Exception as e:
                logger.warning(f"Error loading DE summary: {e}")

        # Load individual contrast results
        de_dir = self.results_dir / "de"
        if de_dir.exists():
            de_data['contrasts'] = {}
            for contrast_file in de_dir.glob("DE_*.tsv"):
                contrast_name = contrast_file.stem.replace('DE_', '')
                try:
                    df = pd.read_csv(contrast_file, sep='\t')
                    de_data['contrasts'][contrast_name] = {
                        'data': df.to_dict('records'),
                        'columns': list(df.columns),
                        'significant_genes': len(df[df['padj'] < 0.05])
                    }
                except Exception as e:
                    logger.warning(f"Error loading contrast {contrast_name}: {e}")

        return de_data

    def load_pathway_data(self) -> Dict[str, Any]:
        """Load pathway enrichment results."""
        pathway_data = {}

        # Load pathway summary
        pathway_summary_path = self.results_dir / "fgsea" / "fgsea_summary.tsv"
        if pathway_summary_path.exists():
            try:
                pathway_summary = pd.read_csv(pathway_summary_path, sep='\t')
                pathway_data['summary'] = pathway_summary.to_dict('records')
            except Exception as e:
                logger.warning(f"Error loading pathway summary: {e}")

        # Load individual contrast pathways
        fgsea_dir = self.results_dir / "fgsea"
        if fgsea_dir.exists():
            pathway_data['contrasts'] = {}
            for contrast_file in fgsea_dir.glob("fgsea_*.tsv"):
                contrast_name = contrast_file.stem.replace('fgsea_', '')
                try:
                    df = pd.read_csv(contrast_file, sep='\t')
                    pathway_data['contrasts'][contrast_name] = {
                        'data': df.to_dict('records'),
                        'columns': list(df.columns),
                        'significant_pathways': len(df[df['padj'] < 0.05])
                    }
                except Exception as e:
                    logger.warning(f"Error loading pathway contrast {contrast_name}: {e}")

        return pathway_data

    def load_counts_data(self) -> Dict[str, Any]:
        """Load gene expression count data."""
        counts_data = {}

        # Load counts matrix
        counts_path = self.results_dir / "counts" / "counts.tsv"
        if counts_path.exists():
            try:
                counts = pd.read_csv(counts_path, sep='\t', index_col=0)
                counts_data['counts'] = {
                    'data': counts.to_dict('records'),
                    'columns': list(counts.columns),
                    'genes': list(counts.index),
                    'samples': list(counts.columns)
                }

                # Calculate basic statistics
                counts_data['stats'] = {
                    'total_genes': len(counts),
                    'total_samples': len(counts.columns),
                    'mean_expression': counts.mean().mean(),
                    'median_expression': counts.median().median()
                }

            except Exception as e:
                logger.warning(f"Error loading counts data: {e}")

        return counts_data

    def get_all_results(self) -> Dict[str, Any]:
        """Load all available results data."""
        if 'all_results' in self.data_cache:
            return self.data_cache['all_results']

        results = {
            'qc': self.load_qc_data(),
            'differential_expression': self.load_differential_expression_data(),
            'pathways': self.load_pathway_data(),
            'counts': self.load_counts_data(),
            'metadata': {
                'results_dir': str(self.results_dir),
                'last_updated': pd.Timestamp.now().isoformat(),
                'available_sections': []
            }
        }

        # Mark available sections
        if results['qc']:
            results['metadata']['available_sections'].append('qc')
        if results['differential_expression']:
            results['metadata']['available_sections'].append('differential_expression')
        if results['pathways']:
            results['metadata']['available_sections'].append('pathways')
        if results['counts']:
            results['metadata']['available_sections'].append('counts')

        self.data_cache['all_results'] = results
        return results


# Initialize results explorer
results_explorer = ResultsExplorer()


@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    """Main dashboard page."""
    results = results_explorer.get_all_results()
    return templates.TemplateResponse("dashboard.html", {
        "request": request,
        "results": results,
        "title": "RNASEQ-MINI Analysis Dashboard"
    })


@app.get("/api/results")
async def get_results():
    """API endpoint to get all results data."""
    return results_explorer.get_all_results()


@app.get("/api/qc/summary")
async def get_qc_summary():
    """Get QC summary data for visualization."""
    qc_data = results_explorer.load_qc_data()

    # Create sample quality overview plot
    if qc_data.get('multiqc', {}).get('metrics'):
        metrics = qc_data['multiqc']['metrics']
        samples = qc_data['multiqc']['samples']

        # Create a simple quality score overview
        quality_data = []
        for sample in samples:
            sample_metrics = next((m for m in metrics if m.get('sample') == sample), {})
            quality_data.append({
                'sample': sample,
                'total_reads': sample_metrics.get('fastqc_total_sequences', 0),
                'avg_quality': sample_metrics.get('fastqc_percent_gc', 50),  # Placeholder
                'gc_content': sample_metrics.get('fastqc_percent_gc', 50)
            })

        return {
            'qc_overview': quality_data,
            'sample_count': len(samples),
            'fastqc_reports': qc_data.get('fastqc', {})
        }

    return {"message": "No QC data available"}


@app.get("/api/de/volcano")
async def get_volcano_plot_data(
    contrast: str = Query(..., description="Contrast name for volcano plot"),
    pvalue_threshold: float = Query(0.05, description="P-value threshold"),
    logfc_threshold: float = Query(1.0, description="Log fold change threshold")
):
    """Generate volcano plot data for differential expression results."""
    de_data = results_explorer.load_differential_expression_data()

    if contrast not in de_data.get('contrasts', {}):
        raise HTTPException(status_code=404, detail=f"Contrast '{contrast}' not found")

    contrast_data = de_data['contrasts'][contrast]['data']

    # Convert to DataFrame for easier processing
    df = pd.DataFrame(contrast_data)

    # Add significance and fold change categories
    df['significant'] = (df['padj'] < pvalue_threshold) & (abs(df['log2FoldChange']) > logfc_threshold)
    df['direction'] = np.where(df['log2FoldChange'] > 0, 'up', 'down')

    # Create volcano plot data
    volcano_data = []
    for _, row in df.iterrows():
        volcano_data.append({
            'gene': row.get('gene_name', row.get('gene_id', 'Unknown')),
            'log2FoldChange': float(row['log2FoldChange']),
            'padj': float(row['padj']),
            'significant': bool(row['significant']),
            'direction': row['direction'],
            'baseMean': float(row.get('baseMean', 0))
        })

    return {
        'volcano_data': volcano_data,
        'thresholds': {
            'pvalue': pvalue_threshold,
            'logfc': logfc_threshold
        },
        'summary': {
            'total_genes': len(df),
            'significant_genes': len(df[df['significant']]),
            'upregulated': len(df[(df['significant']) & (df['log2FoldChange'] > 0)]),
            'downregulated': len(df[(df['significant']) & (df['log2FoldChange'] < 0)])
        }
    }


@app.get("/api/de/heatmap")
async def get_heatmap_data(
    contrast: str = Query(..., description="Contrast name for heatmap"),
    top_n: int = Query(50, description="Number of top genes to show"),
    method: str = Query("padj", description="Sorting method (padj, log2FoldChange)")
):
    """Generate heatmap data for top differentially expressed genes."""
    de_data = results_explorer.load_differential_expression_data()

    if contrast not in de_data.get('contrasts', {}):
        raise HTTPException(status_code=404, detail=f"Contrast '{contrast}' not found")

    contrast_data = de_data['contrasts'][contrast]['data']
    df = pd.DataFrame(contrast_data)

    # Get counts data for expression values
    counts_data = results_explorer.load_counts_data()
    if not counts_data.get('counts'):
        raise HTTPException(status_code=404, detail="No counts data available for heatmap")

    counts_df = pd.DataFrame(counts_data['counts']['data'])
    counts_df = counts_df.set_index(counts_data['counts']['genes'])

    # Select top genes based on specified method
    if method == 'padj':
        top_genes = df.nsmallest(top_n, 'padj')['gene_id'].tolist()
    elif method == 'log2FoldChange':
        top_genes = df.reindex(df['log2FoldChange'].abs().nlargest(top_n).index)['gene_id'].tolist()
    else:
        top_genes = df.head(top_n)['gene_id'].tolist()

    # Filter to genes that exist in counts data
    available_genes = [g for g in top_genes if g in counts_df.index]

    if not available_genes:
        raise HTTPException(status_code=404, detail="No matching genes found in counts data")

    # Create heatmap data
    heatmap_data = counts_df.loc[available_genes].T.to_dict('records')

    return {
        'heatmap_data': heatmap_data,
        'genes': available_genes,
        'samples': counts_data['counts']['samples'],
        'method': method,
        'gene_count': len(available_genes)
    }


@app.get("/api/pathways/enrichment")
async def get_pathway_enrichment(
    contrast: str = Query(..., description="Contrast name for pathway analysis"),
    top_n: int = Query(20, description="Number of top pathways to show"),
    pvalue_threshold: float = Query(0.05, description="P-value threshold")
):
    """Get pathway enrichment analysis results."""
    pathway_data = results_explorer.load_pathway_data()

    if contrast not in pathway_data.get('contrasts', {}):
        raise HTTPException(status_code=404, detail=f"Contrast '{contrast}' not found")

    contrast_data = pathway_data['contrasts'][contrast]['data']
    df = pd.DataFrame(contrast_data)

    # Filter significant pathways
    significant_pathways = df[df['padj'] < pvalue_threshold]

    # Sort by NES (Normalized Enrichment Score) or p-value
    top_pathways = significant_pathways.nlargest(top_n, 'NES')

    pathway_results = []
    for _, row in top_pathways.iterrows():
        pathway_results.append({
            'pathway': row.get('pathway', 'Unknown'),
            'NES': float(row['NES']),
            'padj': float(row['padj']),
            'size': int(row.get('size', 0)),
            'leading_edge_genes': row.get('leadingEdge', '').split(',')[:10]  # Top 10 genes
        })

    return {
        'pathway_data': pathway_results,
        'threshold': pvalue_threshold,
        'total_significant': len(significant_pathways),
        'shown_pathways': len(pathway_results)
    }


@app.get("/api/export/data")
async def export_data(
    data_type: str = Query(..., description="Type of data to export"),
    contrast: Optional[str] = Query(None, description="Contrast name (for DE/pathway data)"),
    format: str = Query("json", description="Export format (json, csv)")
):
    """Export data in various formats."""
    if data_type == "de" and contrast:
        de_data = results_explorer.load_differential_expression_data()
        if contrast in de_data.get('contrasts', {}):
            df = pd.DataFrame(de_data['contrasts'][contrast]['data'])

            if format == "csv":
                return df.to_csv(index=False)
            else:
                return df.to_dict('records')

    elif data_type == "pathways" and contrast:
        pathway_data = results_explorer.load_pathway_data()
        if contrast in pathway_data.get('contrasts', {}):
            df = pd.DataFrame(pathway_data['contrasts'][contrast]['data'])

            if format == "csv":
                return df.to_csv(index=False)
            else:
                return df.to_dict('records')

    elif data_type == "counts":
        counts_data = results_explorer.load_counts_data()
        if counts_data.get('counts'):
            df = pd.DataFrame(counts_data['counts']['data'])

            if format == "csv":
                return df.to_csv(index=False)
            else:
                return df.to_dict('records')

    raise HTTPException(status_code=404, detail=f"Data type '{data_type}' not found or contrast not specified")


@app.get("/api/stats/overview")
async def get_stats_overview():
    """Get overall statistics for the analysis."""
    results = results_explorer.get_all_results()

    stats = {
        'qc': {
            'samples_analyzed': len(results.get('qc', {}).get('fastqc', {})),
            'multiqc_available': bool(results.get('qc', {}).get('multiqc'))
        },
        'differential_expression': {
            'contrasts_analyzed': len(results.get('differential_expression', {}).get('contrasts', {})),
            'total_genes_tested': 0,
            'significant_genes': 0
        },
        'pathways': {
            'contrasts_analyzed': len(results.get('pathways', {}).get('contrasts', {})),
            'total_pathways_tested': 0,
            'significant_pathways': 0
        },
        'counts': {
            'genes_measured': results.get('counts', {}).get('stats', {}).get('total_genes', 0),
            'samples_processed': results.get('counts', {}).get('stats', {}).get('total_samples', 0)
        }
    }

    # Calculate totals
    for contrast_data in results.get('differential_expression', {}).get('contrasts', {}).values():
        stats['differential_expression']['total_genes_tested'] += len(contrast_data.get('data', []))
        stats['differential_expression']['significant_genes'] += contrast_data.get('significant_genes', 0)

    for contrast_data in results.get('pathways', {}).get('contrasts', {}).values():
        stats['pathways']['total_pathways_tested'] += len(contrast_data.get('data', []))
        stats['pathways']['significant_pathways'] += contrast_data.get('significant_pathways', 0)

    return stats


if __name__ == "__main__":
    import uvicorn

    print("🚀 Starting RNASEQ-MINI Interactive Analysis Server")
    print("📊 Open http://localhost:8000 in your browser to explore results")
    print("🔄 Server will reload automatically when files change")

    uvicorn.run(
        "app:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
        reload_dirs=["web_app"]
    )


