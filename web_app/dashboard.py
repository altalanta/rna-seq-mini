import dash
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output
import plotly.express as px
import pandas as pd
import glob
from pathlib import Path
import dash_bootstrap_components as dbc

# --- Data Loading ---
def load_de_results():
    """Load differential expression results from the results directory."""
    results_dir = Path.cwd() / "results"
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
def create_dashboard(server):
    """Create the Dash dashboard and mount it on the FastAPI server."""
    
    summary_df, all_de_data = load_de_results()
    contrasts = list(all_de_data.keys())
    
    dash_app = dash.Dash(
        server=server,
        routes_pathname_prefix="/interactive/",
        external_stylesheets=[dbc.themes.BOOTSTRAP]
    )

    # --- Layout ---
    dash_app.layout = dbc.Container([
        dbc.Row(dbc.Col(html.H1("Interactive Analysis Dashboard"), width=12)),
        
        dbc.Row([
            dbc.Col([
                html.Label("Select Contrast:"),
                dcc.Dropdown(
                    id='contrast-dropdown',
                    options=[{'label': c, 'value': c} for c in contrasts],
                    value=contrasts[0] if contrasts else None
                ),
            ], width=6)
        ], className="mb-4"),

        dbc.Row([
            dbc.Col(dcc.Graph(id='volcano-plot'), width=12)
        ], className="mb-4"),

        dbc.Row([
            dbc.Col([
                html.H3("Differential Expression Results"),
                dash_table.DataTable(
                    id='de-table',
                    columns=[],
                    data=[],
                    sort_action="native",
                    filter_action="native",
                    page_action="native",
                    page_size=15,
                    style_table={'overflowX': 'auto'}
                )
            ], width=12)
        ])
    ], fluid=True)
    
    # --- Callbacks ---
    @dash_app.callback(
        [Output('volcano-plot', 'figure'),
         Output('de-table', 'columns'),
         Output('de-table', 'data')],
        [Input('contrast-dropdown', 'value')]
    )
    def update_dashboard(selected_contrast):
        if not selected_contrast or not all_de_data:
            return px.scatter(title="No data available"), [], []

        df = all_de_data[selected_contrast]
        
        # Volcano plot
        fig = px.scatter(
            df,
            x='log2FoldChange',
            y='padj',
            hover_data=['gene'],
            title=f"Volcano Plot: {selected_contrast}"
        )
        fig.update_yaxes(type="log") # Use log scale for p-values
        
        # DE Table
        columns = [{"name": i, "id": i} for i in df.columns]
        data = df.to_dict('records')
        
        return fig, columns, data

    return dash_app


