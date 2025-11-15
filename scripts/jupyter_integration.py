#!/usr/bin/env python3
"""
Jupyter Notebook Integration for RNASEQ-MINI Collaborative Analysis

This module provides seamless integration between RNASEQ-MINI and Jupyter notebooks,
enabling interactive analysis within collaborative sessions.
"""

import json
import base64
import os
import sys
from pathlib import Path
from typing import Dict, Any, Optional
import subprocess
import tempfile
import webbrowser
from datetime import datetime


class JupyterIntegration:
    """Jupyter notebook integration for collaborative analysis."""

    def __init__(self, collaboration_manager=None):
        self.collaboration_manager = collaboration_manager
        self.notebooks_dir = Path("collaborative_notebooks")
        self.notebooks_dir.mkdir(exist_ok=True)

    def create_analysis_notebook(self, session_id: str, template: str = "standard") -> str:
        """Create a Jupyter notebook for interactive analysis."""
        session = self.collaboration_manager.db.get_session(session_id) if self.collaboration_manager else None

        if not session:
            print("❌ Session not found. Please create or join a session first.")
            return None

        notebook_path = self.notebooks_dir / f"analysis_{session_id}.ipynb"

        # Create notebook content based on template
        notebook_content = self._create_notebook_content(session, template)

        with open(notebook_path, 'w') as f:
            json.dump(notebook_content, f, indent=2)

        print(f"📓 Created analysis notebook: {notebook_path}")
        print(f"🔗 Open in Jupyter: jupyter notebook {notebook_path}")

        return str(notebook_path)

    def _create_notebook_content(self, session: Any, template: str) -> Dict[str, Any]:
        """Create notebook content based on session and template."""
        base_content = {
            "cells": [],
            "metadata": {
                "kernelspec": {
                    "display_name": "RNASEQ-MINI",
                    "language": "python",
                    "name": "rnaseq-mini"
                },
                "language_info": {
                    "name": "python",
                    "version": "3.8+"
                }
            },
            "nbformat": 4,
            "nbformat_minor": 4
        }

        # Add cells based on template
        if template == "standard":
            base_content["cells"] = self._get_standard_cells(session)
        elif template == "qc_focused":
            base_content["cells"] = self._get_qc_cells(session)
        elif template == "de_focused":
            base_content["cells"] = self._get_de_cells(session)

        return base_content

    def _get_standard_cells(self, session: Any) -> List[Dict[str, Any]]:
        """Get cells for standard analysis workflow."""
        return [
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "# 🔬 RNASEQ-MINI Collaborative Analysis\n",
                    f"\n**Session:** {session.name}\n",
                    f"**Project:** {session.project_config.get('project', 'N/A')}\n",
                    f"**Participants:** {len(session.participants)}\n",
                    "\nThis notebook provides interactive analysis capabilities for your collaborative RNA-seq study."
                ]
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## 📊 Session Overview\n",
                    "\nLet's start by examining the current analysis session and available data."
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# Import RNASEQ-MINI collaboration tools\n",
                    "import sys\n",
                    "sys.path.append('../scripts')\n",
                    "from collaboration_manager import CollaborationManager\n",
                    "\n",
                    "# Connect to collaboration session\n",
                    "collab = CollaborationManager()\n",
                    "session = collab.db.get_session('" + (session.id if hasattr(session, 'id') else 'SESSION_ID') + "')\n",
                    "\n",
                    "print(f\"Session: {session.name}\")\n",
                    "print(f\"Current step: {session.current_step}\")\n",
                    "print(f\"Participants: {len(session.participants)}\")"
                ]
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## 🔍 Data Exploration\n",
                    "\nExamine your input data and quality metrics."
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# Load sample information\n",
                    "import pandas as pd\n",
                    "import yaml\n",
                    "\n",
                    "# Load configuration\n",
                    "with open('../config/params.yaml') as f:\n",
                    "    config = yaml.safe_load(f)\n",
                    "\n",
                    "print(\"Project configuration:\")\n",
                    "for key, value in config.items():\n",
                    "    if not key.startswith('_'):\n",
                    "        print(f\"  {key}: {value}\")"
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# Check for existing results\n",
                    "import os\n",
                    "from pathlib import Path\n",
                    "\n",
                    "results_dir = Path('../results')\n",
                    "if results_dir.exists():\n",
                    "    print(\"Available result files:\")\n",
                    "    for file in results_dir.rglob('*'):\n",
                    "        if file.is_file():\n",
                    "            print(f\"  {file.relative_to(results_dir)}\")\n",
                    "else:\n",
                    "    print(\"No results found. Run analysis first with: make run\")"
                ]
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## 📈 Quality Control Analysis\n",
                    "\nInteractive exploration of quality control metrics."
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# Load and display MultiQC report if available\n",
                    "import IPython.display as display\n",
                    "\n",
                    "multiqc_path = Path('../results/qc/multiqc_report.html')\n",
                    "if multiqc_path.exists():\n",
                    "    print(\"MultiQC Report:\")\n",
                    "    display.IFrame(str(multiqc_path), width='100%', height=600)\n",
                    "else:\n",
                    "    print(\"MultiQC report not found. Run quality control first.\")"
                ]
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## 🧬 Differential Expression Analysis\n",
                    "\nInteractive exploration of differential expression results."
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# Load DESeq2 results\n",
                    "import pandas as pd\n",
                    "import plotly.express as px\n",
                    "import plotly.graph_objects as go\n",
                    "\n",
                    "de_results_path = Path('../results/de/deseq2_results.tsv')\n",
                    "if de_results_path.exists():\n",
                    "    de_results = pd.read_csv(de_results_path, sep='\\t')\n",
                    "    \n",
                    "    # Create volcano plot\n",
                    "    fig = go.Figure()\n",
                    "    \n",
                    "    # Color by significance\n",
                    "    colors = ['red' if (p < 0.05 and abs(lfc) > 1) else 'blue' \n",
                    "              for p, lfc in zip(de_results['padj'], de_results['log2FoldChange'])]\n",
                    "    \n",
                    "    fig.add_trace(go.Scatter(\n",
                    "        x=de_results['log2FoldChange'],\n",
                    "        y=-de_results['padj'].apply(lambda x: -math.log10(x) if x > 0 else 0),\n",
                    "        mode='markers',\n",
                    "        marker=dict(color=colors, size=4),\n",
                    "        text=de_results['gene_id'],\n",
                    "        hovertemplate='%{text}<br>log2FC: %{x:.2f}<br>-log10(padj): %{y:.2f}'\n",
                    "    ))\n",
                    "    \n",
                    "    fig.update_layout(\n",
                    "        title='Volcano Plot - Differential Expression',\n",
                    "        xaxis_title='log2 Fold Change',\n",
                    "        yaxis_title='-log10(adjusted p-value)',\n",
                    "        showlegend=False\n",
                    "    )\n",
                    "    \n",
                    "    fig.show()\n",
                    "    \n",
                    "    print(f\"Total genes: {len(de_results)}\")\n",
                    "    print(f\"Significant genes (padj < 0.05, |LFC| > 1): {sum((de_results['padj'] < 0.05) & (abs(de_results['log2FoldChange']) > 1))}\")\n",
                    "else:\n",
                    "    print(\"DESeq2 results not found. Run differential expression analysis first.\")"
                ]
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## 🛣️ Pathway Enrichment Analysis\n",
                    "\nInteractive exploration of pathway enrichment results."
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# Load fgsea results\n",
                    "fgsea_path = Path('../results/fgsea/fgsea_results.tsv')\n",
                    "if fgsea_path.exists():\n",
                    "    fgsea_results = pd.read_csv(fgsea_path, sep='\\t')\n",
                    "    \n",
                    "    # Filter significant pathways\n",
                    "    sig_pathways = fgsea_results[fgsea_results['padj'] < 0.05].head(20)\n",
                    "    \n",
                    "    # Create enrichment plot\n",
                    "    fig = px.bar(\n",
                    "        sig_pathways,\n",
                    "        x='NES',\n",
                    "        y='pathway',\n",
                    "        orientation='h',\n",
                    "        color='padj',\n",
                    "        color_continuous_scale='RdBu_r',\n",
                    "        title='Top Enriched Pathways'\n",
                    "    )\n",
                    "    \n",
                    "    fig.update_layout(yaxis={'categoryorder': 'total ascending'})\n",
                    "    fig.show()\n",
                    "    \n",
                    "    print(f\"Total pathways tested: {len(fgsea_results)}\")\n",
                    "    print(f\"Significant pathways (padj < 0.05): {len(sig_pathways)}\")\n",
                    "else:\n",
                    "    print(\"fgsea results not found. Run pathway analysis first.\")"
                ]
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## 🤝 Collaboration Features\n",
                    "\nUse these cells to interact with your collaborative session."
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# Update analysis step\n",
                    "collab.update_session_step('" + (session.id if hasattr(session, 'id') else 'SESSION_ID') + "', 'notebook_analysis', 'current_user')\n",
                    "\n",
                    "print(\"Updated session step to 'notebook_analysis'\")"
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# Log a custom analysis event\n",
                    "collab.log_analysis_event(\n",
                    "    '" + (session.id if hasattr(session, 'id') else 'SESSION_ID') + "',\n",
                    "    'current_user',\n",
                    "    'custom_analysis',\n",
                    "    'Performed custom analysis in Jupyter notebook',\n",
                    "    {'notebook_cell': 'custom_analysis', 'timestamp': str(pd.Timestamp.now())}\n",
                    ")\n",
                    "\n",
                    "print(\"Logged custom analysis event\")"
                ]
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## 💾 Export Results\n",
                    "\nExport your analysis results for sharing or publication."
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# Export filtered results\n",
                    "significant_genes = de_results[(de_results['padj'] < 0.05) & (abs(de_results['log2FoldChange']) > 1)]\n",
                    "\n",
                    "# Save to CSV\n",
                    "significant_genes.to_csv('../results/significant_genes.csv', index=False)\n",
                    "\n",
                    "# Export to Excel with multiple sheets\n",
                    "with pd.ExcelWriter('../results/analysis_results.xlsx') as writer:\n",
                    "    de_results.to_excel(writer, sheet_name='All_Results', index=False)\n",
                    "    significant_genes.to_excel(writer, sheet_name='Significant_Genes', index=False)\n",
                    "    if 'fgsea_results' in locals():\n",
                    "        sig_pathways.to_excel(writer, sheet_name='Enriched_Pathways', index=False)\n",
                    "\n",
                    "print(\"Exported results to CSV and Excel formats\")"
                ]
            }
        ]

    def _get_qc_cells(self, session: Any) -> List[Dict[str, Any]]:
        """Get cells focused on quality control analysis."""
        return [
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "# 🔍 Quality Control Deep Dive\n",
                    "\nDetailed exploration of RNA-seq quality control metrics."
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# Load and analyze FastQC results\n",
                    "import json\n",
                    "import glob\n",
                    "\n",
                    "fastqc_dir = Path('../results/qc/fastqc')\n",
                    "if fastqc_dir.exists():\n",
                    "    fastqc_files = list(fastqc_dir.glob('*_fastqc.zip'))\n",
                    "    print(f\"Found {len(fastqc_files)} FastQC reports\")\n",
                    "    \n",
                    "    # Extract summary data (simplified example)\n",
                    "    for file in fastqc_files[:3]:  # Show first 3\n",
                    "        print(f\"\\n📊 {file.stem}:\")\n",
                    "        print(f\"   File: {file}\")\n",
                    "        print(f\"   Size: {file.stat().st_size / 1024:.1f} KB\")\n",
                    "else:\n",
                    "    print(\"FastQC reports not found. Run quality control first.\")"
                ]
            }
        ]

    def _get_de_cells(self, session: Any) -> List[Dict[str, Any]]:
        """Get cells focused on differential expression analysis."""
        return [
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "# 🧬 Differential Expression Deep Dive\n",
                    "\nAdvanced exploration of differential expression results."
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# Advanced DE analysis\n",
                    "# Load normalized counts and perform custom statistical tests\n",
                    "import numpy as np\n",
                    "from scipy import stats\n",
                    "\n",
                    "counts_path = Path('../results/counts/counts.tsv')\n",
                    "if counts_path.exists():\n",
                    "    counts = pd.read_csv(counts_path, sep='\\t', index_col=0)\n",
                    "    print(f\"Loaded counts matrix: {counts.shape}\")\n",
                    "    \n",
                    "    # Example: Calculate coefficient of variation\n",
                    "    cv = counts.std() / counts.mean()\n",
                    "    high_variation_genes = cv[cv > cv.quantile(0.9)].index\n",
                    "    \n",
                    "    print(f\"Genes with high variation (>90th percentile): {len(high_variation_genes)}\")\n",
                    "else:\n",
                    "    print(\"Counts matrix not found.\")"
                ]
            }
        ]

    def launch_notebook_server(self, notebook_path: str, port: int = 8888):
        """Launch Jupyter notebook server for collaborative editing."""
        try:
            # Check if notebook exists
            if not Path(notebook_path).exists():
                print(f"❌ Notebook not found: {notebook_path}")
                return False

            # Start Jupyter server
            cmd = [
                "jupyter", "notebook",
                str(notebook_path),
                "--port", str(port),
                "--no-browser",
                "--NotebookApp.token=''",  # No token for easy access
                "--NotebookApp.password=''",  # No password for easy access
                "--NotebookApp.allow_origin='*'"
            ]

            print(f"🚀 Starting Jupyter server on port {port}")
            print(f"📖 Notebook: {notebook_path}")
            print(f"🔗 Access at: http://localhost:{port}")

            # Open in browser
            webbrowser.open(f"http://localhost:{port}")

            # Start server
            subprocess.run(cmd, cwd=self.notebooks_dir.parent)

        except Exception as e:
            print(f"❌ Failed to start Jupyter server: {e}")
            print("💡 Try: pip install jupyter notebook")
            return False

    def embed_collaborative_widgets(self, session_id: str):
        """Embed collaborative widgets in the notebook."""
        session = self.collaboration_manager.db.get_session(session_id) if self.collaboration_manager else None

        if not session:
            print("❌ Session not found")
            return

        # Create collaborative status widget
        widget_code = f"""
# Collaborative Session Widget
print("🛜 Connected to session: {session.name}")
print(f"👥 Participants: {{len(session.participants)}}")
print(f"📍 Current step: {session.current_step}")

# Add real-time updates (JavaScript integration would go here)
from IPython.display import HTML
HTML('''
<div id="collab-status" style="border: 1px solid #ddd; padding: 10px; margin: 10px 0;">
    <h4>🛜 Live Collaboration Status</h4>
    <p>Session: {session.name}</p>
    <p>Participants: {len(session.participants)}</p>
    <p>Current Step: {session.current_step}</p>
    <div id="live-updates" style="max-height: 200px; overflow-y: auto;">
        <!-- Real-time updates would appear here -->
    </div>
</div>
''')
"""

        print("📋 Collaborative widget code generated")
        print("💡 Copy this code into your Jupyter notebook:")
        print("=" * 50)
        print(widget_code)
        print("=" * 50)


def create_collaborative_notebook(session_id: str, template: str = "standard") -> Optional[str]:
    """Create a collaborative Jupyter notebook for a session."""
    integration = JupyterIntegration()
    return integration.create_analysis_notebook(session_id, template)


def launch_collaborative_notebook(notebook_path: str, port: int = 8888) -> bool:
    """Launch a collaborative Jupyter notebook server."""
    integration = JupyterIntegration()
    return integration.launch_notebook_server(notebook_path, port)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Jupyter Integration for RNASEQ-MINI")
    subparsers = parser.add_subparsers(dest="command")

    # Create notebook command
    create_parser = subparsers.add_parser("create", help="Create a collaborative notebook")
    create_parser.add_argument("session_id", help="Session ID")
    create_parser.add_argument("--template", choices=["standard", "qc_focused", "de_focused"],
                              default="standard", help="Notebook template")

    # Launch notebook command
    launch_parser = subparsers.add_parser("launch", help="Launch Jupyter server")
    launch_parser.add_argument("notebook_path", help="Path to notebook file")
    launch_parser.add_argument("--port", type=int, default=8888, help="Server port")

    # Embed widgets command
    embed_parser = subparsers.add_parser("embed", help="Generate collaborative widgets")
    embed_parser.add_argument("session_id", help="Session ID")

    args = parser.parse_args()

    if args.command == "create":
        notebook_path = create_collaborative_notebook(args.session_id, args.template)
        if notebook_path:
            print(f"\n🚀 Launch with: python scripts/jupyter_integration.py launch {notebook_path}")
    elif args.command == "launch":
        success = launch_collaborative_notebook(args.notebook_path, args.port)
        if not success:
            sys.exit(1)
    elif args.command == "embed":
        integration = JupyterIntegration()
        integration.embed_collaborative_widgets(args.session_id)
    else:
        parser.print_help()









