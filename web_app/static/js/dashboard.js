// RNASEQ-MINI Interactive Dashboard JavaScript

// Global variables
let currentResults = {};
let currentContrast = '';
let volcanoPlot = null;
let heatmapPlot = null;
let pathwayPlot = null;

// Initialize dashboard when page loads
document.addEventListener('DOMContentLoaded', function() {
    loadInitialData();
    setupEventListeners();
});

// Load initial data and populate UI
async function loadInitialData() {
    try {
        showLoading();

        // Load all results data
        const response = await fetch('/api/results');
        currentResults = await response.json();

        // Populate dropdowns with available contrasts
        populateContrastDropdowns();

        // Load statistics overview
        await loadStatsOverview();

        // Load QC data
        await loadQCData();

        hideLoading();

    } catch (error) {
        console.error('Error loading initial data:', error);
        hideLoading();
    }
}

// Setup event listeners for interactive elements
function setupEventListeners() {
    // Contrast selection dropdowns
    document.getElementById('de-contrast-select').addEventListener('change', function() {
        currentContrast = this.value;
        if (currentContrast) {
            loadVolcanoPlot();
            loadHeatmap();
        }
    });

    document.getElementById('pathway-contrast-select').addEventListener('change', function() {
        currentContrast = this.value;
        if (currentContrast) {
            loadPathwayPlot();
        }
    });

    // Export functionality
    document.getElementById('export-type').addEventListener('change', function() {
        const contrastSelect = document.getElementById('export-contrast');
        if (this.value === 'counts') {
            contrastSelect.disabled = true;
            contrastSelect.value = '';
        } else {
            contrastSelect.disabled = false;
        }
    });
}

// Populate contrast dropdowns with available contrasts
function populateContrastDropdowns() {
    const deContrasts = Object.keys(currentResults.differential_expression?.contrasts || {});
    const pathwayContrasts = Object.keys(currentResults.pathways?.contrasts || {});

    // Populate DE contrast dropdown
    const deSelect = document.getElementById('de-contrast-select');
    deContrasts.forEach(contrast => {
        const option = document.createElement('option');
        option.value = contrast;
        option.textContent = contrast.replace('_vs_', ' vs ');
        deSelect.appendChild(option);
    });

    // Populate pathway contrast dropdown
    const pathwaySelect = document.getElementById('pathway-contrast-select');
    pathwayContrasts.forEach(contrast => {
        const option = document.createElement('option');
        option.value = contrast;
        option.textContent = contrast.replace('_vs_', ' vs ');
        pathwaySelect.appendChild(option);
    });

    // Populate export contrast dropdown
    const exportSelect = document.getElementById('export-contrast');
    const allContrasts = [...new Set([...deContrasts, ...pathwayContrasts])];
    allContrasts.forEach(contrast => {
        const option = document.createElement('option');
        option.value = contrast;
        option.textContent = contrast.replace('_vs_', ' vs ');
        exportSelect.appendChild(option);
    });
}

// Load statistics overview
async function loadStatsOverview() {
    try {
        const response = await fetch('/api/stats/overview');
        const stats = await response.json();

        const statsContainer = document.getElementById('stats-cards');
        statsContainer.innerHTML = `
            <div class="col-md-3">
                <div class="card bg-primary text-white">
                    <div class="card-body">
                        <h5 class="card-title">${stats.qc.samples_analyzed}</h5>
                        <p class="card-text">Samples Analyzed</p>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card bg-success text-white">
                    <div class="card-body">
                        <h5 class="card-title">${stats.differential_expression.significant_genes}</h5>
                        <p class="card-text">Significant Genes</p>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card bg-warning text-white">
                    <div class="card-body">
                        <h5 class="card-title">${stats.pathways.significant_pathways}</h5>
                        <p class="card-text">Significant Pathways</p>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card bg-info text-white">
                    <div class="card-body">
                        <h5 class="card-title">${stats.counts.genes_measured.toLocaleString()}</h5>
                        <p class="card-text">Genes Measured</p>
                    </div>
                </div>
            </div>
        `;

    } catch (error) {
        console.error('Error loading stats overview:', error);
    }
}

// Load QC data and create visualizations
async function loadQCData() {
    try {
        const response = await fetch('/api/qc/summary');
        const qcData = await response.json();

        if (qcData.qc_overview) {
            // Create QC overview plot
            const samples = qcData.qc_overview.map(d => d.sample);
            const totalReads = qcData.qc_overview.map(d => d.total_reads);
            const avgQuality = qcData.qc_overview.map(d => d.avg_quality);

            const trace1 = {
                x: samples,
                y: totalReads,
                type: 'bar',
                name: 'Total Reads',
                marker: { color: '#1f77b4' }
            };

            const trace2 = {
                x: samples,
                y: avgQuality,
                type: 'scatter',
                mode: 'lines+markers',
                name: 'Avg Quality',
                yaxis: 'y2',
                marker: { color: '#ff7f0e' }
            };

            const layout = {
                title: 'Sample Quality Overview',
                xaxis: { title: 'Sample' },
                yaxis: { title: 'Total Reads' },
                yaxis2: {
                    title: 'Average Quality',
                    overlaying: 'y',
                    side: 'right'
                },
                height: 350
            };

            Plotly.newPlot('qc-plot', [trace1, trace2], layout);
        }

        // Populate FastQC links
        const fastqcContainer = document.getElementById('fastqc-links');
        if (qcData.fastqc_reports) {
            let linksHtml = '<ul class="list-group">';
            for (const [sample, info] of Object.entries(qcData.fastqc_reports)) {
                linksHtml += `
                    <li class="list-group-item d-flex justify-content-between align-items-center">
                        ${sample}
                        <span class="badge bg-primary rounded-pill">
                            <i class="fas fa-file-alt me-1"></i>Report
                        </span>
                    </li>
                `;
            }
            linksHtml += '</ul>';
            fastqcContainer.innerHTML = linksHtml;
        }

    } catch (error) {
        console.error('Error loading QC data:', error);
    }
}

// Load and display volcano plot
async function loadVolcanoPlot() {
    if (!currentContrast) return;

    try {
        showLoading();

        const pvalueThreshold = parseFloat(document.getElementById('pvalue-threshold').value);
        const logfcThreshold = parseFloat(document.getElementById('logfc-threshold').value);

        const response = await fetch(`/api/de/volcano?contrast=${currentContrast}&pvalue_threshold=${pvalueThreshold}&logfc_threshold=${logfcThreshold}`);
        const data = await response.json();

        createVolcanoPlot(data);

        hideLoading();

    } catch (error) {
        console.error('Error loading volcano plot:', error);
        hideLoading();
    }
}

// Create interactive volcano plot
function createVolcanoPlot(data) {
    const volcanoData = data.volcano_data;

    // Separate significant and non-significant points
    const significant = volcanoData.filter(d => d.significant);
    const nonSignificant = volcanoData.filter(d => !d.significant);

    const traces = [];

    // Non-significant points (gray, smaller)
    if (nonSignificant.length > 0) {
        traces.push({
            x: nonSignificant.map(d => d.log2FoldChange),
            y: nonSignificant.map(d => -Math.log10(d.padj)),
            mode: 'markers',
            type: 'scatter',
            name: 'Not Significant',
            marker: {
                color: '#cccccc',
                size: 4,
                opacity: 0.6
            },
            hovertemplate: '<b>%{text}</b><br>Log2FC: %{x:.2f}<br>-Log10(P): %{y:.2f}<extra></extra>',
            text: nonSignificant.map(d => d.gene)
        });
    }

    // Significant upregulated (red)
    const upregulated = significant.filter(d => d.direction === 'up');
    if (upregulated.length > 0) {
        traces.push({
            x: upregulated.map(d => d.log2FoldChange),
            y: upregulated.map(d => -Math.log10(d.padj)),
            mode: 'markers',
            type: 'scatter',
            name: 'Upregulated',
            marker: {
                color: '#d62728',
                size: 6
            },
            hovertemplate: '<b>%{text}</b><br>Log2FC: %{x:.2f}<br>-Log10(P): %{y:.2f}<extra></extra>',
            text: upregulated.map(d => d.gene)
        });
    }

    // Significant downregulated (blue)
    const downregulated = significant.filter(d => d.direction === 'down');
    if (downregulated.length > 0) {
        traces.push({
            x: downregulated.map(d => d.log2FoldChange),
            y: downregulated.map(d => -Math.log10(d.padj)),
            mode: 'markers',
            type: 'scatter',
            name: 'Downregulated',
            marker: {
                color: '#1f77b4',
                size: 6
            },
            hovertemplate: '<b>%{text}</b><br>Log2FC: %{x:.2f}<br>-Log10(P): %{y:.2f}<extra></extra>',
            text: downregulated.map(d => d.gene)
        });
    }

    // Threshold lines
    traces.push({
        x: [data.thresholds.logfc, data.thresholds.logfc],
        y: [0, 100],
        mode: 'lines',
        name: 'LogFC Threshold',
        line: { color: 'red', dash: 'dash' },
        showlegend: false
    });

    traces.push({
        x: [-data.thresholds.logfc, -data.thresholds.logfc],
        y: [0, 100],
        mode: 'lines',
        name: 'LogFC Threshold',
        line: { color: 'red', dash: 'dash' },
        showlegend: false
    });

    traces.push({
        x: [-10, 10],
        y: [-Math.log10(data.thresholds.pvalue), -Math.log10(data.thresholds.pvalue)],
        mode: 'lines',
        name: 'P-value Threshold',
        line: { color: 'red', dash: 'dot' },
        showlegend: false
    });

    const layout = {
        title: `Volcano Plot - ${currentContrast.replace('_vs_', ' vs ')}`,
        xaxis: { title: 'Log2 Fold Change' },
        yaxis: { title: '-Log10(Adjusted P-value)' },
        height: 450,
        hovermode: 'closest',
        legend: { x: 1, y: 1 }
    };

    volcanoPlot = Plotly.newPlot('volcano-plot', traces, layout);
}

// Load and display expression heatmap
async function loadHeatmap(method = 'padj') {
    if (!currentContrast) return;

    try {
        showLoading();

        const topN = 50; // Default value
        const response = await fetch(`/api/de/heatmap?contrast=${currentContrast}&top_n=${topN}&method=${method}`);
        const data = await response.json();

        createHeatmap(data);

        hideLoading();

    } catch (error) {
        console.error('Error loading heatmap:', error);
        hideLoading();
    }
}

// Create interactive heatmap
function createHeatmap(data) {
    const heatmapData = data.heatmap_data;
    const genes = data.genes;
    const samples = data.samples;

    // Create z-values for heatmap
    const z = genes.map(gene => samples.map(sample => heatmapData.find(d => d.gene_id === gene)?.[sample] || 0));

    // Create annotations for gene labels
    const annotations = genes.map((gene, i) =>
        Object.assign({}, {
            x: samples.length + 0.5,
            y: i,
            text: gene,
            showarrow: false,
            font: { size: 8 },
            align: 'left'
        })
    );

    const layout = {
        title: `Expression Heatmap - Top ${data.gene_count} Genes (${data.method})`,
        xaxis: { title: 'Sample' },
        yaxis: { title: 'Gene' },
        height: 550,
        annotations: annotations
    };

    const trace = {
        z: z,
        x: samples,
        y: genes,
        type: 'heatmap',
        colorscale: 'RdBu',
        reversescale: true,
        hovertemplate: 'Sample: %{x}<br>Gene: %{y}<br>Expression: %{z}<extra></extra>'
    };

    heatmapPlot = Plotly.newPlot('heatmap-plot', [trace], layout);
}

// Load and display pathway enrichment plot
async function loadPathwayPlot() {
    if (!currentContrast) return;

    try {
        showLoading();

        const pvalueThreshold = parseFloat(document.getElementById('pathway-threshold').value);
        const topN = parseInt(document.getElementById('top-pathways').value);

        const response = await fetch(`/api/pathways/enrichment?contrast=${currentContrast}&top_n=${topN}&pvalue_threshold=${pvalueThreshold}`);
        const data = await response.json();

        createPathwayPlot(data);

        hideLoading();

    } catch (error) {
        console.error('Error loading pathway plot:', error);
        hideLoading();
    }
}

// Create interactive pathway enrichment plot
function createPathwayPlot(data) {
    const pathwayData = data.pathway_data;

    const pathways = pathwayData.map(d => d.pathway);
    const nesValues = pathwayData.map(d => d.NES);
    const padjValues = pathwayData.map(d => -Math.log10(d.padj));

    const trace1 = {
        x: pathways,
        y: nesValues,
        type: 'bar',
        name: 'NES',
        marker: { color: nesValues.map(v => v > 0 ? '#d62728' : '#1f77b4') }
    };

    const trace2 = {
        x: pathways,
        y: padjValues,
        type: 'scatter',
        mode: 'markers',
        name: '-Log10(P-adj)',
        yaxis: 'y2',
        marker: { color: '#ff7f0e', size: 8 }
    };

    const layout = {
        title: `Pathway Enrichment - ${currentContrast.replace('_vs_', ' vs ')}`,
        xaxis: { title: 'Pathway' },
        yaxis: { title: 'Normalized Enrichment Score (NES)' },
        yaxis2: {
            title: '-Log10(Adjusted P-value)',
            overlaying: 'y',
            side: 'right'
        },
        height: 450,
        margin: { l: 200 },
        showlegend: false
    };

    pathwayPlot = Plotly.newPlot('pathway-plot', [trace1, trace2], layout);

    // Populate pathway details
    populatePathwayDetails(pathwayData);
}

// Populate pathway details panel
function populatePathwayDetails(pathwayData) {
    const container = document.getElementById('pathway-details');
    let html = '<div class="list-group">';

    pathwayData.forEach((pathway, index) => {
        html += `
            <div class="list-group-item">
                <div class="d-flex w-100 justify-content-between">
                    <h6 class="mb-1">${index + 1}. ${pathway.pathway}</h6>
                    <small>NES: ${pathway.NES.toFixed(2)}, P: ${pathway.padj.toExponential(2)}</small>
                </div>
                <p class="mb-1">Size: ${pathway.size} genes</p>
                <small class="text-muted">Leading genes: ${pathway.leading_edge_genes.join(', ')}</small>
            </div>
        `;
    });

    html += '</div>';
    container.innerHTML = html;
}

// Update functions for interactive elements
function updateVolcanoPlot() {
    loadVolcanoPlot();
}

function updateHeatmap(method) {
    loadHeatmap(method);
}

function updatePathwayPlot() {
    loadPathwayPlot();
}

// Export data functionality
async function exportData() {
    const dataType = document.getElementById('export-type').value;
    const contrast = document.getElementById('export-contrast').value;
    const format = document.getElementById('export-format').value;

    if (!dataType) {
        alert('Please select a data type to export');
        return;
    }

    if ((dataType === 'de' || dataType === 'pathways') && !contrast) {
        alert('Please select a contrast for differential expression or pathway data');
        return;
    }

    try {
        showLoading();

        const response = await fetch(`/api/export/data?data_type=${dataType}&contrast=${contrast}&format=${format}`);

        if (format === 'json') {
            const data = await response.json();
            downloadJSON(data, `${dataType}_${contrast || 'all'}.${format}`);
        } else {
            const csvData = await response.text();
            downloadCSV(csvData, `${dataType}_${contrast || 'all'}.${format}`);
        }

        hideLoading();

    } catch (error) {
        console.error('Error exporting data:', error);
        hideLoading();
        alert('Error exporting data. Please try again.');
    }
}

// Download functions
function downloadJSON(data, filename) {
    const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
}

function downloadCSV(data, filename) {
    const blob = new Blob([data], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
}

// Refresh all data
async function refreshAllData() {
    try {
        showLoading();
        await loadInitialData();
        hideLoading();
    } catch (error) {
        console.error('Error refreshing data:', error);
        hideLoading();
    }
}

// Utility functions
function showLoading() {
    const modal = new bootstrap.Modal(document.getElementById('loadingModal'));
    modal.show();
}

function hideLoading() {
    const modal = bootstrap.Modal.getInstance(document.getElementById('loadingModal'));
    if (modal) {
        modal.hide();
    }
}

// Auto-refresh functionality (optional)
setInterval(async function() {
    // Check if results have been updated (could be improved with proper file watching)
    try {
        const response = await fetch('/api/results');
        const newResults = await response.json();

        // Simple check - if the timestamp changed, refresh overview stats
        if (newResults.metadata?.last_updated !== currentResults.metadata?.last_updated) {
            currentResults = newResults;
            loadStatsOverview();
        }
    } catch (error) {
        // Silently fail - this is just for auto-refresh
    }
}, 30000); // Check every 30 seconds

