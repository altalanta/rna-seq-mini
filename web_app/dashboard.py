import dash
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
import glob
from pathlib import Path
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import numpy as np
import json
from api.server import get_job_data
import yaml

# --- Data Loading ---
def load_multiqc_stats(results_dir):
    """Load general stats from the MultiQC report."""
    try:
        stats_file = results_dir / "qc" / "multiqc" / "multiqc_data" / "multiqc_general_stats.txt"
        df = pd.read_csv(stats_file, sep='	')
        # We need to pivot the table to get samples as rows and stats as columns
        df_pivot = df.set_index('Sample').T
        return df_pivot
    except FileNotFoundError:
        return pd.DataFrame()

def load_normalized_counts(results_dir):
    """Load normalized counts (TPM)."""
    try:
        counts_file = results_dir / "counts" / "tpm.tsv"
        df = pd.read_csv(counts_file, sep='	', index_col=0)
        return df
    except FileNotFoundError:
        return pd.DataFrame()

def load_de_results(results_dir):
    """Load differential expression results from a specified results directory."""
    # This function is now superseded by get_job_data, but we keep it for structure
    # The actual data loading for DE will happen inside the callbacks via get_job_data
    de_files = glob.glob(str(results_dir / "de" / "DE_*.tsv"))
    if not de_files:
        return {}, {}
    
    contrasts = [Path(f).stem.replace("DE_", "") for f in de_files]
    return {c: None for c in contrasts}, {} # Return dummy data

def load_samples_info(results_dir):
    """Load sample metadata to map samples to conditions."""
    try:
        # A bit of a hack: find the params file used for the run to locate the samples file
        params_file = results_dir / "params.yaml"
        with open(params_file, 'r') as f:
            config = yaml.safe_load(f)
        samples_path = Path(config["paths"]["samples"])
        df = pd.read_csv(samples_path, sep='	')
        return df.set_index('sample')['condition'].to_dict()
    except Exception:
        return {}

def load_gsea_results(results_dir, contrast):
    """Load GSEA results for a specific contrast."""
    try:
        gsea_file = results_dir / "fgsea" / contrast / "fgsea_results.tsv"
        df = pd.read_csv(gsea_file, sep='	')
        return df
    except FileNotFoundError:
        return pd.DataFrame()

# --- App Initialization ---
def create_dashboard(server, results_dir_path="results"):
    """Create the Dash dashboard and mount it on the FastAPI server."""
    results_dir = Path(results_dir_path)

    summary_df, all_de_data = load_de_results(results_dir)
    multiqc_stats = load_multiqc_stats(results_dir)
    tpm_counts = load_normalized_counts(results_dir)
    contrasts = list(all_de_data.keys())
    samples_map = load_samples_info(results_dir)

    dash_app = dash.Dash(
        server=server,
        routes_pathname_prefix="/interactive/",
        external_stylesheets=[dbc.themes.BOOTSTRAP]
    )

    # --- Layout ---
    # Helper function to create summary cards
    def make_summary_card(title, value, icon):
        return dbc.Card(
            dbc.CardBody([
                html.H4(title, className="card-title"),
                html.P(f"{value}", className="card-text fs-3"),
                html.I(className=f"fas {icon} fa-2x position-absolute top-50 end-0 translate-middle-y me-3 text-muted")
            ]),
            className="mb-4 position-relative"
        )

    # Prepare summary stats
    total_reads = multiqc_stats.get('total_sequences', pd.Series([0])).sum() * 2 if not multiqc_stats.empty else 0
    total_genes = len(tpm_counts) if not tpm_counts.empty else 0
    num_contrasts = len(contrasts)

    dash_app.layout = dbc.Container(fluid=True, children=[
        # Header
        dbc.Row([
            dbc.Col(html.H1("RNA-seq Analysis Dashboard", className="text-white"), width=10),
            dbc.Col(
                dbc.Button("Start Tour", id="start-tour-button", color="info", className="mt-2"),
                width=2,
                className="text-end"
            )
        ], className="bg-primary p-3 align-items-center"),

        # Main content
        dbc.Tabs([
            dbc.Tab(label="Differential Expression", children=[
                dbc.Row([
                    dbc.Col([
                        html.Label("Select Contrast:"),
                        dcc.Dropdown(
                            id='contrast-dropdown',
                            options=[{'label': c, 'value': c} for c in contrasts],
                            value=contrasts[0] if contrasts else None
                        ),
                    ], width=6, className="mt-3")
                ]),
                dbc.Row([
                    dbc.Col(dcc.Graph(id='volcano-plot'), width=7),
                    dbc.Col(dcc.Graph(id='heatmap'), width=5),
                ], className="mt-4"),
                dbc.Row([
                    dbc.Col([
                        html.H3("Differential Expression Results"),
                        dash_table.DataTable(
                            id='de-table',
                            columns=[], data=[], sort_action="native", filter_action="native", page_action="native", page_size=15,
                            style_table={'overflowX': 'auto'}
                        )
                    ], width=12, className="mt-4")
                ])
            ]),
            dbc.Tab(label="Pathway Analysis", children=[
                dbc.Row([
                    dbc.Col([
                        html.H4("GSEA Results"),
                        dash_table.DataTable(
                            id='gsea-table',
                            sort_action="native",
                            filter_action="native",
                            page_action="native",
                            page_size=15,
                            style_table={'overflowX': 'auto'},
                            style_cell={'textAlign': 'left', 'minWidth': '120px'},
                        )
                    ], width=12, className="mt-3")
                ]),
                dbc.Row([
                    dbc.Col([
                        html.H4("Top Pathways Plot"),
                        html.Img(id='gsea-plot-table', style={'width': '100%'})
                    ], width=12, className="mt-4")
                ])
            ]),
            dbc.Tab(label="QC Summary", children=[
                html.H3("MultiQC Summary", className="mt-3"),
                dash_table.DataTable(
                    id='multiqc-table',
                    columns=[{"name": i, "id": i} for i in multiqc_stats.reset_index().columns],
                    data=multiqc_stats.reset_index().to_dict('records'),
                    sort_action="native", filter_action="native", page_action="native", page_size=15,
                    style_table={'overflowX': 'auto'}
                )
            ]),
            dbc.Tab(label="Gene View", children=[
                dbc.Row([
                    dbc.Col(dcc.Input(id='gene-input', type='text', placeholder='Enter Gene ID...', debounce=True), width=4),
                    dbc.Col(html.Button('Search', id='gene-search-button', n_clicks=0), width=2),
                ], className="mt-3"),
                dbc.Row([
                    dbc.Col(dcc.Graph(id='gene-expression-plot'), width=6),
                    dbc.Col([
                        html.H4("Differential Expression Stats"),
                        dash_table.DataTable(
                            id='gene-de-table',
                            style_table={'overflowX': 'auto'}
                        )
                    ], width=6)
                ], className="mt-4")
            ]),
        ])
    ], fluid=True)


    # --- Callbacks ---
    @dash_app.callback(
        [Output('volcano-plot', 'figure'),
         Output('de-table', 'columns'),
         Output('de-table', 'data'),
         Output('heatmap', 'figure')],
        [Input('contrast-dropdown', 'value')]
    )
    def update_de_visuals(selected_contrast):
        if not selected_contrast or not all_de_data:
            empty_fig = go.Figure().update_layout(title="No data available", xaxis_visible=False, yaxis_visible=False)
            return empty_fig, [], [], empty_fig

        df = all_de_data[selected_contrast].copy()
        
        # Volcano plot (code from previous step)
        df['color'] = np.where((df['padj'] < 0.05) & (abs(df['log2FoldChange']) > 1), 'significant', 'non-significant')
        df['-log10(padj)'] = -np.log10(df['padj'].replace(0, 1e-300))
        volcano_fig = px.scatter(
            df, x='log2FoldChange', y='-log10(padj)', hover_data=['gene'], color='color',
            color_discrete_map={'significant': 'red', 'non-significant': 'gray'},
            title=f"Volcano Plot: {selected_contrast}",
            labels={'log2FoldChange': 'Log2 Fold Change', '-log10(padj)': '-Log10 Adjusted P-value'}
        )
        volcano_fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="green")
        volcano_fig.add_vline(x=1, line_dash="dash", line_color="blue")
        volcano_fig.add_vline(x=-1, line_dash="dash", line_color="blue")

        # Heatmap of top 20 significant genes
        sig_genes = df[(df['padj'] < 0.05) & (abs(df['log2FoldChange']) > 1)].sort_values('padj').head(20)
        if not sig_genes.empty and not tpm_counts.empty:
            # Check if gene names exist in the TPM index before trying to locate them
            valid_genes = sig_genes['gene'][sig_genes['gene'].isin(tpm_counts.index)]
            if not valid_genes.empty:
                heatmap_data = tpm_counts.loc[valid_genes].dropna()
                # Log transform and center the data for better visualization
                heatmap_data = np.log2(heatmap_data + 1)
                heatmap_data = heatmap_data.sub(heatmap_data.mean(axis=1), axis=0)
                
                heatmap_fig = px.imshow(
                    heatmap_data,
                    labels=dict(x="Sample", y="Gene", color="Log2(TPM)"),
                    title="Top 20 DE Genes"
                )
            else:
                heatmap_fig = go.Figure().update_layout(title="No significant genes to display in heatmap", xaxis_visible=False, yaxis_visible=False)
        else:
            heatmap_fig = go.Figure().update_layout(title="No significant genes to display in heatmap", xaxis_visible=False, yaxis_visible=False)

        # DE Table
        columns = [{"name": i, "id": i} for i in df.columns if i not in ['color', '-log10(padj)']]
        data = df.to_dict('records')

        return volcano_fig, columns, data, heatmap_fig

    @dash_app.callback(
        [Output('gene-expression-plot', 'figure'),
         Output('gene-de-table', 'data'),
         Output('gene-de-table', 'columns')],
        [Input('gene-search-button', 'n_clicks')],
        [State('gene-input', 'value')]
    )
    def update_gene_view(n_clicks, gene_id):
        if not n_clicks or not gene_id:
            empty_fig = go.Figure().update_layout(title="Enter a gene ID and click Search", xaxis_visible=False, yaxis_visible=False)
            return empty_fig, [], []

        job_data = get_job_data(str(results_dir))
        tpm_df = job_data["tpm"]
        de_data = job_data["de"]
        
        # Expression Plot
        if not tpm_df.empty and gene_id in tpm_df.index:
            expression_data = tpm_df.loc[gene_id]
            exp_df = pd.DataFrame({'sample': expression_data.index, 'tpm': expression_data.values})
            exp_df['condition'] = exp_df['sample'].map(samples_map)
            
            expression_fig = px.box(exp_df, x='condition', y='tpm', title=f"TPM Expression for {gene_id}", points="all")
        else:
            expression_fig = go.Figure().update_layout(title=f"Gene {gene_id} not found in TPM data", xaxis_visible=False, yaxis_visible=False)

        # DE Stats Table
        de_stats = []
        for contrast, df in de_data.items():
            gene_row = df[df['gene'] == gene_id]
            if not gene_row.empty:
                stats = gene_row.iloc[0].to_dict()
                stats['contrast'] = contrast
                de_stats.append(stats)
        
        de_columns = [
            {"name": "Contrast", "id": "contrast"},
            {"name": "log2FC", "id": "log2FoldChange"},
            {"name": "p-value", "id": "pvalue"},
            {"name": "adj p-value", "id": "padj"}
        ]
        
        return expression_fig, de_stats, de_columns

    @dash_app.callback(
        [Output('gsea-table', 'columns'),
         Output('gsea-table', 'data'),
         Output('gsea-plot-table', 'src')],
        [Input('contrast-dropdown', 'value')]
    )
    def update_gsea_view(selected_contrast):
        if not selected_contrast:
            return [], [], ""

        gsea_df = load_gsea_results(results_dir, selected_contrast)
        
        if gsea_df.empty:
            return [], [], ""

        # Prepare table data
        gsea_df['leadingEdge'] = gsea_df['leadingEdge'].apply(lambda x: ', '.join(x.split(' ')))
        columns = [{"name": i, "id": i} for i in gsea_df.columns]
        data = gsea_df.to_dict('records')

        # Prepare plot path
        plot_path = results_dir / "fgsea" / selected_contrast / "top_pathways_table.pdf"
        # Dash can't serve PDFs directly, so for now we will return an empty string.
        # A more advanced implementation could convert the PDF to a PNG.
        plot_src = "" # Placeholder

        return columns, data, plot_src

    return dash_app.server


