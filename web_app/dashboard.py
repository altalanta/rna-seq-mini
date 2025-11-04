import dash
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output
import plotly.express as px
import pandas as pd
import glob
from pathlib import Path
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import numpy as np

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
    de_files = glob.glob(str(results_dir / "de" / "DE_*.tsv"))

    if not de_files:
        return pd.DataFrame(), {}

    all_de = {}
    for f in de_files:
        contrast = Path(f).stem.replace("DE_", "")
        df = pd.read_csv(f, sep='\t')
        df.rename(columns={'Unnamed: 0': 'gene'}, inplace=True)
        all_de[contrast] = df

    # Create a summary of significant genes
    summary = []
    for contrast, df in all_de.items():
        significant = df[df['padj'] < 0.05]
        summary.append({
            "contrast": contrast,
            "upregulated": (significant['log2FoldChange'] > 0).sum(),
            "downregulated": (significant['log2FoldChange'] < 0).sum()
        })

    summary_df = pd.DataFrame(summary)
    return summary_df, all_de

# --- App Initialization ---
def create_dashboard(server, results_dir_path="results"):
    """Create the Dash dashboard and mount it on the FastAPI server."""
    results_dir = Path(results_dir_path)

    summary_df, all_de_data = load_de_results(results_dir)
    multiqc_stats = load_multiqc_stats(results_dir)
    tpm_counts = load_normalized_counts(results_dir)
    contrasts = list(all_de_data.keys())

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

    dash_app.layout = dbc.Container([
        dbc.Row(dbc.Col(html.H1("Interactive Analysis Dashboard"), width=12, className="mb-4")),

        # Summary Cards
        dbc.Row([
            dbc.Col(make_summary_card("Total Reads", f"{total_reads/1e6:.1f}M", "fa-book-open"), md=4),
            dbc.Col(make_summary_card("Total Genes", f"{total_genes:,}", "fa-dna"), md=4),
            dbc.Col(make_summary_card("Contrasts", num_contrasts, "fa-balance-scale"), md=4),
        ]),

        # Main Tabs
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
            heatmap_data = tpm_counts.loc[sig_genes['gene']].dropna()
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

        # DE Table
        columns = [{"name": i, "id": i} for i in df.columns if i not in ['color', '-log10(padj)']]
        data = df.to_dict('records')

        return volcano_fig, columns, data, heatmap_fig

    return dash_app


