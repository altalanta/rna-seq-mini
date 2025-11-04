// RNASEQ-MINI Interactive Dashboard JavaScript

// Global variables
let currentResults = {};
let currentContrast = '';
let volcanoPlot = null;
let heatmapPlot = null;
let pathwayPlot = null;
let currentDeData = [];
let currentPathwayData = [];
let filteredDeData = [];
let filteredPathwayData = [];

// Initialize dashboard when page loads
document.addEventListener('DOMContentLoaded', function() {
    loadInitialData();
    setupEventListeners();
    setupMobileInteractions();
    loadConfigurationsFromStorage();
    restoreStateFromUrl();
    initializeProgressTracking();
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

        // Initialize comparison view
        initializeComparisonView();

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

    // Table search and filter functionality
    document.getElementById('de-search').addEventListener('input', debounce(filterDeTable, 300));
    document.getElementById('pathway-search').addEventListener('input', debounce(filterPathwayTable, 300));

    document.getElementById('de-significance-filter').addEventListener('change', filterDeTable);
    document.getElementById('pathway-significance-filter').addEventListener('change', filterPathwayTable);

    document.getElementById('de-logfc-filter').addEventListener('input', debounce(filterDeTable, 300));
    document.getElementById('pathway-nes-filter').addEventListener('input', debounce(filterPathwayTable, 300));
}

// Debounce function for search inputs
function debounce(func, wait) {
    let timeout;
    return function executedFunction(...args) {
        const later = () => {
            clearTimeout(timeout);
            func(...args);
        };
        clearTimeout(timeout);
        timeout = setTimeout(later, wait);
    };
}

// Mobile Interactions Setup
function setupMobileInteractions() {
    // Add touch feedback for clickable elements
    addTouchFeedback();

    // Setup responsive navigation
    setupResponsiveNavigation();

    // Setup touch gestures for plots
    setupPlotGestures();

    // Setup keyboard navigation for accessibility
    setupKeyboardNavigation();

    // Setup viewport detection for responsive adjustments
    handleViewportChanges();
}

// Add touch feedback for clickable elements
function addTouchFeedback() {
    // Add touch feedback to clickable genes
    document.querySelectorAll('.clickable-gene').forEach(element => {
        element.addEventListener('touchstart', function() {
            this.style.backgroundColor = 'rgba(13, 110, 253, 0.2)';
        });

        element.addEventListener('touchend', function() {
            setTimeout(() => {
                this.style.backgroundColor = '';
            }, 150);
        });
    });

    // Add touch feedback to buttons
    document.querySelectorAll('.btn').forEach(button => {
        button.addEventListener('touchstart', function() {
            this.style.transform = 'scale(0.95)';
        });

        button.addEventListener('touchend', function() {
            setTimeout(() => {
                this.style.transform = '';
            }, 150);
        });
    });
}

// Setup responsive navigation
function setupResponsiveNavigation() {
    // Smooth scroll for navigation links
    document.querySelectorAll('a[href^="#"]').forEach(anchor => {
        anchor.addEventListener('click', function(e) {
            e.preventDefault();
            const target = document.querySelector(this.getAttribute('href'));

            if (target) {
                // Close mobile navbar if open
                const navbarCollapse = document.querySelector('.navbar-collapse');
                if (navbarCollapse && navbarCollapse.classList.contains('show')) {
                    const navbarToggler = document.querySelector('.navbar-toggler');
                    if (navbarToggler) {
                        navbarToggler.click();
                    }
                }

                // Smooth scroll to target
                target.scrollIntoView({
                    behavior: 'smooth',
                    block: 'start'
                });
            }
        });
    });

    // Auto-close mobile navbar when clicking outside
    document.addEventListener('click', function(e) {
        const navbarCollapse = document.querySelector('.navbar-collapse');
        const navbarToggler = document.querySelector('.navbar-toggler');

        if (navbarCollapse && navbarCollapse.classList.contains('show') &&
            !navbarCollapse.contains(e.target) && !navbarToggler.contains(e.target)) {
            navbarToggler.click();
        }
    });
}

// Setup touch gestures for plots
function setupPlotGestures() {
    // Enable touch gestures for Plotly plots
    const plotElements = ['volcano-plot', 'heatmap-plot', 'pathway-plot', 'qc-plot'];

    plotElements.forEach(plotId => {
        const plotElement = document.getElementById(plotId);
        if (plotElement) {
            // Add double-tap to zoom functionality for mobile
            let lastTouchTime = 0;

            plotElement.addEventListener('touchend', function(e) {
                const currentTime = new Date().getTime();
                const timeDiff = currentTime - lastTouchTime;

                if (timeDiff < 300 && timeDiff > 0) {
                    // Double tap detected - could implement zoom here
                    e.preventDefault();
                }

                lastTouchTime = currentTime;
            });
        }
    });
}

// Setup keyboard navigation for accessibility
function setupKeyboardNavigation() {
    // Add keyboard support for modal interactions
    document.addEventListener('keydown', function(e) {
        // ESC key to close modals
        if (e.key === 'Escape') {
            const openModal = document.querySelector('.modal.show');
            if (openModal) {
                const closeButton = openModal.querySelector('.btn-close');
                if (closeButton) {
                    closeButton.click();
                }
            }
        }

        // Enter key for buttons and links
        if (e.key === 'Enter' && e.target.tagName === 'BUTTON') {
            e.target.click();
        }
    });

    // Add focus indicators for keyboard navigation
    document.querySelectorAll('button, a, input, select').forEach(element => {
        element.addEventListener('focus', function() {
            this.setAttribute('data-focused', 'true');
        });

        element.addEventListener('blur', function() {
            this.removeAttribute('data-focused');
        });
    });
}

// Handle viewport changes for responsive adjustments
function handleViewportChanges() {
    // Detect orientation changes
    window.addEventListener('orientationchange', function() {
        setTimeout(function() {
            // Force plot resize after orientation change
            window.dispatchEvent(new Event('resize'));

            // Update modal positioning if open
            const openModal = document.querySelector('.modal.show');
            if (openModal && window.innerWidth <= 576) {
                // Adjust modal for mobile landscape
                openModal.style.margin = '0.5rem';
            }
        }, 100);
    });

    // Handle window resize for responsive adjustments
    window.addEventListener('resize', debounce(function() {
        // Update plot sizes if needed
        const isMobile = window.innerWidth <= 768;

        if (isMobile) {
            // Add mobile-specific classes or adjustments
            document.body.classList.add('mobile-view');
        } else {
            document.body.classList.remove('mobile-view');
        }
    }, 250));
}

// Sharing and Configuration Functions
// Save current configuration to local storage
function saveCurrentConfiguration() {
    const configName = document.getElementById('config-name').value.trim();

    if (!configName) {
        alert('Please enter a name for this configuration');
        return;
    }

    // Get current state
    const currentState = getCurrentState();

    // Get existing configurations
    let configurations = JSON.parse(localStorage.getItem('rnaseq-configs') || '{}');

    // Add timestamp and metadata
    currentState.metadata = {
        name: configName,
        created: new Date().toISOString(),
        version: '1.0'
    };

    // Save configuration
    configurations[configName] = currentState;

    // Store in localStorage
    localStorage.setItem('rnaseq-configs', JSON.stringify(configurations));

    // Update UI
    document.getElementById('config-name').value = '';
    loadConfigurationsFromStorage();

    // Show success message
    showAlert('Configuration saved successfully!', 'success');
}

// Load configurations from localStorage
function loadConfigurationsFromStorage() {
    const configurations = JSON.parse(localStorage.getItem('rnaseq-configs') || '{}');

    if (Object.keys(configurations).length === 0) {
        document.getElementById('saved-configurations').style.display = 'none';
        return;
    }

    document.getElementById('saved-configurations').style.display = 'block';

    const configsList = document.getElementById('saved-configs-list');
    configsList.innerHTML = '';

    Object.entries(configurations).forEach(([name, config]) => {
        const configItem = document.createElement('div');
        configItem.className = 'list-group-item d-flex justify-content-between align-items-center';
        configItem.innerHTML = `
            <div>
                <strong>${name}</strong><br>
                <small class="text-muted">Created: ${new Date(config.metadata.created).toLocaleDateString()}</small>
            </div>
            <button class="btn btn-sm btn-outline-primary" onclick="loadConfigurationFromName('${name}')">
                <i class="fas fa-upload me-1"></i>Load
            </button>
        `;
        configsList.appendChild(configItem);
    });

    // Update load dropdown
    const loadSelect = document.getElementById('load-config-select');
    loadSelect.innerHTML = '<option value="">Select a saved configuration...</option>';

    Object.keys(configurations).forEach(name => {
        const option = document.createElement('option');
        option.value = name;
        option.textContent = name;
        loadSelect.appendChild(option);
    });
}

// Get current application state
function getCurrentState() {
    return {
        // Current contrast selections
        currentContrast: currentContrast,
        deContrast: document.getElementById('de-contrast-select').value,
        pathwayContrast: document.getElementById('pathway-contrast-select').value,

        // Thresholds and filters
        pvalueThreshold: document.getElementById('pvalue-threshold').value,
        logfcThreshold: document.getElementById('logfc-threshold').value,
        pathwayThreshold: document.getElementById('pathway-threshold').value,
        topPathways: document.getElementById('top-pathways').value,

        // Table filters
        deSearch: document.getElementById('de-search').value,
        deSignificanceFilter: document.getElementById('de-significance-filter').value,
        deLogfcFilter: document.getElementById('de-logfc-filter').value,
        pathwaySearch: document.getElementById('pathway-search').value,
        pathwaySignificanceFilter: document.getElementById('pathway-significance-filter').value,
        pathwayNesFilter: document.getElementById('pathway-nes-filter').value,

        // Table visibility
        deTableVisible: document.getElementById('de-table-container').style.display !== 'none',
        pathwayTableVisible: document.getElementById('pathway-table-card').style.display !== 'none',

        // Gene search
        geneSearchQuery: document.getElementById('gene-search-input').value,

        // Active section (for navigation)
        activeSection: window.location.hash || '#overview'
    };
}

// Restore state from configuration
function restoreState(state) {
    // Restore contrast selections
    if (state.deContrast) {
        document.getElementById('de-contrast-select').value = state.deContrast;
        currentContrast = state.deContrast;
        if (currentContrast) {
            loadVolcanoPlot();
            loadHeatmap();
        }
    }

    if (state.pathwayContrast) {
        document.getElementById('pathway-contrast-select').value = state.pathwayContrast;
        if (state.pathwayContrast) {
            loadPathwayPlot();
        }
    }

    // Restore thresholds
    if (state.pvalueThreshold) {
        document.getElementById('pvalue-threshold').value = state.pvalueThreshold;
    }
    if (state.logfcThreshold) {
        document.getElementById('logfc-threshold').value = state.logfcThreshold;
    }
    if (state.pathwayThreshold) {
        document.getElementById('pathway-threshold').value = state.pathwayThreshold;
    }
    if (state.topPathways) {
        document.getElementById('top-pathways').value = state.topPathways;
    }

    // Restore table filters
    if (state.deSearch) {
        document.getElementById('de-search').value = state.deSearch;
    }
    if (state.deSignificanceFilter) {
        document.getElementById('de-significance-filter').value = state.deSignificanceFilter;
    }
    if (state.deLogfcFilter) {
        document.getElementById('de-logfc-filter').value = state.deLogfcFilter;
    }
    if (state.pathwaySearch) {
        document.getElementById('pathway-search').value = state.pathwaySearch;
    }
    if (state.pathwaySignificanceFilter) {
        document.getElementById('pathway-significance-filter').value = state.pathwaySignificanceFilter;
    }
    if (state.pathwayNesFilter) {
        document.getElementById('pathway-nes-filter').value = state.pathwayNesFilter;
    }

    // Restore table visibility
    if (state.deTableVisible) {
        toggleDeTable();
    }
    if (state.pathwayTableVisible) {
        togglePathwayTable();
    }

    // Restore gene search
    if (state.geneSearchQuery) {
        document.getElementById('gene-search-input').value = state.geneSearchQuery;
    }

    // Navigate to active section
    if (state.activeSection && state.activeSection !== '#overview') {
        setTimeout(() => {
            document.querySelector(`a[href="${state.activeSection}"]`).click();
        }, 500);
    }

    // Apply filters if data is loaded
    if (state.deSearch || state.deSignificanceFilter || state.deLogfcFilter) {
        filterDeTable();
    }
    if (state.pathwaySearch || state.pathwaySignificanceFilter || state.pathwayNesFilter) {
        filterPathwayTable();
    }
}

// Load configuration by name
function loadConfigurationFromName(configName) {
    const configurations = JSON.parse(localStorage.getItem('rnaseq-configs') || '{}');
    const config = configurations[configName];

    if (config) {
        restoreState(config);
        showAlert(`Configuration "${configName}" loaded successfully!`, 'success');
    }
}

// Load configuration from dropdown
function loadConfiguration() {
    const configName = document.getElementById('load-config-select').value;
    if (configName) {
        loadConfigurationFromName(configName);
    }
}

// Delete configuration
function deleteConfiguration() {
    const configName = document.getElementById('load-config-select').value;
    if (!configName) {
        alert('Please select a configuration to delete');
        return;
    }

    if (confirm(`Are you sure you want to delete the configuration "${configName}"?`)) {
        let configurations = JSON.parse(localStorage.getItem('rnaseq-configs') || '{}');
        delete configurations[configName];
        localStorage.setItem('rnaseq-configs', JSON.stringify(configurations));

        loadConfigurationsFromStorage();
        showAlert(`Configuration "${configName}" deleted successfully!`, 'info');
    }
}

// Generate shareable URL with current state
function generateShareUrl() {
    const state = getCurrentState();
    const stateString = btoa(JSON.stringify(state));
    const baseUrl = window.location.origin + window.location.pathname;
    const shareUrl = `${baseUrl}?state=${stateString}`;

    document.getElementById('share-url').value = shareUrl;

    // Show success message
    showAlert('Share URL generated! You can now copy and share this link.', 'success');
}

// Copy share URL to clipboard
function copyShareUrl() {
    const shareUrlInput = document.getElementById('share-url');

    if (!shareUrlInput.value) {
        alert('Please generate a share URL first');
        return;
    }

    shareUrlInput.select();
    shareUrlInput.setSelectionRange(0, 99999); // For mobile devices

    try {
        document.execCommand('copy');
        showAlert('Share URL copied to clipboard!', 'success');
    } catch (err) {
        // Fallback for modern browsers
        navigator.clipboard.writeText(shareUrlInput.value).then(() => {
            showAlert('Share URL copied to clipboard!', 'success');
        }).catch(() => {
            alert('Failed to copy URL. Please copy manually.');
        });
    }
}

// Restore state from URL parameters
function restoreStateFromUrl() {
    const urlParams = new URLSearchParams(window.location.search);
    const stateParam = urlParams.get('state');

    if (stateParam) {
        try {
            const state = JSON.parse(atob(stateParam));
            restoreState(state);

            // Show that URL state was restored
            showAlert('Analysis state restored from shared URL!', 'info');

            // Clean URL after loading
            const cleanUrl = window.location.origin + window.location.pathname + window.location.hash;
            window.history.replaceState({}, document.title, cleanUrl);

        } catch (error) {
            console.error('Error restoring state from URL:', error);
            showAlert('Error loading shared configuration. The link may be corrupted.', 'danger');
        }
    }
}

// Show alert messages
function showAlert(message, type = 'info') {
    // Create alert element
    const alertDiv = document.createElement('div');
    alertDiv.className = `alert alert-${type} alert-dismissible fade show position-fixed`;
    alertDiv.style.cssText = 'top: 20px; right: 20px; z-index: 9999; min-width: 300px;';
    alertDiv.innerHTML = `
        ${message}
        <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
    `;

    document.body.appendChild(alertDiv);

    // Auto-remove after 5 seconds
    setTimeout(() => {
        if (alertDiv.parentNode) {
            alertDiv.remove();
        }
    }, 5000);
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

        // Also load table data if table is visible
        if (document.getElementById('de-table-container').style.display !== 'none') {
            currentDeData = data.volcano_data;
            filteredDeData = [...currentDeData];
            populateDeTable(filteredDeData);
            updateDeTableInfo(filteredDeData.length, currentDeData.length);
        }

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

    // Non-significant points (gray, smaller, clickable)
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
            hovertemplate: '<b>%{text}</b><br>Log2FC: %{x:.2f}<br>-Log10(P): %{y:.2f}<br><i>Click for details</i><extra></extra>',
            text: nonSignificant.map(d => d.gene),
            customdata: nonSignificant.map(d => ({
                gene: d.gene,
                gene_id: d.gene_id || d.gene,
                log2FoldChange: d.log2FoldChange,
                padj: d.padj,
                baseMean: d.baseMean
            }))
        });
    }

    // Significant upregulated (red, clickable)
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
            hovertemplate: '<b>%{text}</b><br>Log2FC: %{x:.2f}<br>-Log10(P): %{y:.2f}<br><i>Click for details</i><extra></extra>',
            text: upregulated.map(d => d.gene),
            customdata: upregulated.map(d => ({
                gene: d.gene,
                gene_id: d.gene_id || d.gene,
                log2FoldChange: d.log2FoldChange,
                padj: d.padj,
                baseMean: d.baseMean
            }))
        });
    }

    // Significant downregulated (blue, clickable)
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
            hovertemplate: '<b>%{text}</b><br>Log2FC: %{x:.2f}<br>-Log10(P): %{y:.2f}<br><i>Click for details</i><extra></extra>',
            text: downregulated.map(d => d.gene),
            customdata: downregulated.map(d => ({
                gene: d.gene,
                gene_id: d.gene_id || d.gene,
                log2FoldChange: d.log2FoldChange,
                padj: d.padj,
                baseMean: d.baseMean
            }))
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

    // Add click event listener to the plot
    volcanoPlot.then(function() {
        const plotDiv = document.getElementById('volcano-plot');
        plotDiv.on('plotly_click', function(data) {
            if (data.points && data.points.length > 0) {
                const point = data.points[0];
                const customData = point.customdata;
                if (customData) {
                    showGeneDetails(customData.gene_id, customData.gene, currentContrast);
                }
            }
        });
    });
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

    // Create annotations for gene labels (clickable)
    const annotations = genes.map((gene, i) =>
        Object.assign({}, {
            x: samples.length + 0.5,
            y: i,
            text: gene,
            showarrow: false,
            font: { size: 8 },
            align: 'left',
            clicktoshow: 'on'
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
        hovertemplate: 'Sample: %{x}<br>Gene: %{y}<br>Expression: %{z}<br><i>Click for details</i><extra></extra>',
        customdata: genes.map(gene => ({ gene_id: gene, gene_name: gene }))
    };

    heatmapPlot = Plotly.newPlot('heatmap-plot', [trace], layout);

    // Add click event listener to the heatmap
    heatmapPlot.then(function() {
        const plotDiv = document.getElementById('heatmap-plot');
        plotDiv.on('plotly_click', function(data) {
            if (data.points && data.points.length > 0) {
                const point = data.points[0];
                const customData = point.customdata;
                if (customData) {
                    showGeneDetails(customData.gene_id, customData.gene_name, currentContrast);
                }
            }
        });
    });
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

        // Also load table data if table is visible
        if (document.getElementById('pathway-table-card').style.display !== 'none') {
            currentPathwayData = data.pathway_data;
            filteredPathwayData = [...currentPathwayData];
            populatePathwayTable(filteredPathwayData);
            updatePathwayTableInfo(filteredPathwayData.length, currentPathwayData.length);
        }

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

// Gene Search Functions
// Search for a specific gene across all contrasts
async function searchForGene() {
    // Check both the main search input and navbar search input
    let geneQuery = document.getElementById('gene-search-input').value.trim() ||
                   document.getElementById('navbar-gene-search').value.trim();

    if (!geneQuery) {
        alert('Please enter a gene name or ID to search for');
        return;
    }

    try {
        showLoading();

        const response = await fetch(`/api/genes/search?gene_query=${encodeURIComponent(geneQuery)}&exact_match=false`);
        const data = await response.json();

        displayGeneSearchResults(data);

        // Scroll to results
        document.getElementById('gene-search-results').scrollIntoView({ behavior: 'smooth' });

        hideLoading();

    } catch (error) {
        console.error('Error searching for gene:', error);
        hideLoading();
        alert('Error searching for gene. Please try again.');
    }
}

// Display gene search results
function displayGeneSearchResults(data) {
    const resultsContainer = document.getElementById('gene-search-results');
    const tableBody = document.querySelector('#gene-search-table tbody');

    if (data.results.length === 0) {
        resultsContainer.style.display = 'block';
        tableBody.innerHTML = '<tr><td colspan="6" class="text-center text-muted">No genes found matching your search</td></tr>';
        return;
    }

    resultsContainer.style.display = 'block';
    tableBody.innerHTML = '';

    data.results.forEach(result => {
        const row = document.createElement('tr');

        const significantBadge = result.significant
            ? '<span class="badge bg-success">Yes</span>'
            : '<span class="badge bg-secondary">No</span>';

        row.innerHTML = `
            <td><strong>${result.gene}</strong></td>
            <td>${result.contrast.replace('_vs_', ' vs ')}</td>
            <td>${result.log2FoldChange.toFixed(3)}</td>
            <td>${result.padj.toExponential(2)}</td>
            <td>${significantBadge}</td>
            <td>
                <button class="btn btn-sm btn-outline-primary" onclick="viewGeneInVolcano('${result.gene}', '${result.contrast}', ${result.log2FoldChange}, ${-Math.log10(result.padj)})">
                    <i class="fas fa-volcano me-1"></i>View in Plot
                </button>
            </td>
        `;

        tableBody.appendChild(row);
    });
}

// View gene in volcano plot (placeholder for future implementation)
function viewGeneInVolcano(geneName, contrast, logFC, negLogP) {
    // Switch to the DE section and select the appropriate contrast
    document.querySelector('a[href="#de"]').click();

    // Set the contrast dropdown to the correct value
    const deSelect = document.getElementById('de-contrast-select');
    deSelect.value = contrast;
    currentContrast = contrast;

    // Trigger the change event to load the data
    deSelect.dispatchEvent(new Event('change'));

    console.log(`Would highlight gene: ${geneName} at (${logFC}, ${negLogP}) in contrast: ${contrast}`);
    // This would be implemented to highlight specific points in the volcano plot
}

// Clear gene search results
function clearGeneSearch() {
    document.getElementById('gene-search-input').value = '';
    document.getElementById('navbar-gene-search').value = '';
    document.getElementById('gene-search-results').style.display = 'none';
}

// Table Management Functions
// Toggle DE table visibility
function toggleDeTable() {
    const container = document.getElementById('de-table-container');
    const button = document.querySelector('#de-table-container').previousElementSibling.querySelector('.btn');

    if (container.style.display === 'none') {
        container.style.display = 'block';
        button.innerHTML = '<i class="fas fa-eye-slash me-1"></i>Hide Table';
        if (currentDeData.length === 0 && currentContrast) {
            loadDeTableData();
        }
    } else {
        container.style.display = 'none';
        button.innerHTML = '<i class="fas fa-eye me-1"></i>Show Table';
    }
}

// Toggle pathway table visibility
function togglePathwayTable() {
    const tableCard = document.getElementById('pathway-table-card');
    const button = document.querySelector('[onclick="togglePathwayTable()"]');

    if (tableCard.style.display === 'none') {
        tableCard.style.display = 'block';
        button.innerHTML = '<i class="fas fa-eye-slash me-1"></i>Hide Table';
        if (currentPathwayData.length === 0 && currentContrast) {
            loadPathwayTableData();
        }
    } else {
        tableCard.style.display = 'none';
        button.innerHTML = '<i class="fas fa-table me-1"></i>View Table';
    }
}

// Load DE data for table
async function loadDeTableData() {
    if (!currentContrast) return;

    try {
        const pvalueThreshold = parseFloat(document.getElementById('pvalue-threshold').value);
        const logfcThreshold = parseFloat(document.getElementById('logfc-threshold').value);

        const response = await fetch(`/api/de/volcano?contrast=${currentContrast}&pvalue_threshold=${pvalueThreshold}&logfc_threshold=${logfcThreshold}`);
        const data = await response.json();

        currentDeData = data.volcano_data;
        filteredDeData = [...currentDeData];
        populateDeTable(filteredDeData);

        updateDeTableInfo(filteredDeData.length, currentDeData.length);

    } catch (error) {
        console.error('Error loading DE table data:', error);
    }
}

// Load pathway data for table
async function loadPathwayTableData() {
    if (!currentContrast) return;

    try {
        const pvalueThreshold = parseFloat(document.getElementById('pathway-threshold').value);
        const topN = parseInt(document.getElementById('top-pathways').value);

        const response = await fetch(`/api/pathways/enrichment?contrast=${currentContrast}&top_n=${topN}&pvalue_threshold=${pvalueThreshold}`);
        const data = await response.json();

        currentPathwayData = data.pathway_data;
        filteredPathwayData = [...currentPathwayData];
        populatePathwayTable(filteredPathwayData);

        updatePathwayTableInfo(filteredPathwayData.length, currentPathwayData.length);

    } catch (error) {
        console.error('Error loading pathway table data:', error);
    }
}

// Populate DE table with data
function populateDeTable(data) {
    const tbody = document.querySelector('#de-results-table tbody');
    tbody.innerHTML = '';

    data.forEach((row, index) => {
        const tr = document.createElement('tr');

        // Make gene name clickable to show gene details
        const geneCell = document.createElement('td');
        geneCell.innerHTML = `<a href="#" onclick="showGeneDetails('${row.gene_id || row.gene}', '${row.gene}', '${currentContrast}')" style="color: #0d6efd; text-decoration: none;" title="Click for detailed gene information">${row.gene}</a>`;
        tr.appendChild(geneCell);

        tr.innerHTML += `
            <td>${row.baseMean.toFixed(1)}</td>
            <td>${row.log2FoldChange.toFixed(3)}</td>
            <td>${row.padj.toExponential(2)}</td>
            <td>${row.padj.toExponential(2)}</td>
            <td><span class="badge ${row.significant ? 'bg-success' : 'bg-secondary'}">${row.significant ? 'Yes' : 'No'}</span></td>
            <td><span class="badge ${row.direction === 'up' ? 'bg-danger' : 'bg-primary'}">${row.direction === 'up' ? 'Up' : 'Down'}</span></td>
        `;

        tbody.appendChild(tr);
    });
}

// Populate pathway table with data
function populatePathwayTable(data) {
    const tbody = document.querySelector('#pathway-results-table tbody');
    tbody.innerHTML = '';

    data.forEach((row, index) => {
        const tr = document.createElement('tr');

        // Make pathway name clickable
        const pathwayCell = document.createElement('td');
        pathwayCell.innerHTML = `<strong>${row.pathway}</strong>`;
        tr.appendChild(pathwayCell);

        tr.innerHTML += `
            <td>${row.NES.toFixed(3)}</td>
            <td>${row.padj.toExponential(2)}</td>
            <td>${row.padj.toExponential(2)}</td>
            <td>${row.size}</td>
            <td><small class="text-muted">${row.leading_edge_genes.slice(0, 3).join(', ')}${row.leading_edge_genes.length > 3 ? '...' : ''}</small></td>
        `;

        tbody.appendChild(tr);
    });
}

// Filter DE table
function filterDeTable() {
    const searchTerm = document.getElementById('de-search').value.toLowerCase();
    const significanceFilter = document.getElementById('de-significance-filter').value;
    const logfcFilter = parseFloat(document.getElementById('de-logfc-filter').value) || 0;

    filteredDeData = currentDeData.filter(row => {
        // Search filter
        if (searchTerm && !row.gene.toLowerCase().includes(searchTerm)) {
            return false;
        }

        // Significance filter
        if (significanceFilter === 'significant' && !row.significant) {
            return false;
        }
        if (significanceFilter === 'up' && (!row.significant || row.direction !== 'up')) {
            return false;
        }
        if (significanceFilter === 'down' && (!row.significant || row.direction !== 'down')) {
            return false;
        }

        // LogFC filter
        if (logfcFilter > 0 && Math.abs(row.log2FoldChange) < logfcFilter) {
            return false;
        }

        return true;
    });

    populateDeTable(filteredDeData);
    updateDeTableInfo(filteredDeData.length, currentDeData.length);
}

// Filter pathway table
function filterPathwayTable() {
    const searchTerm = document.getElementById('pathway-search').value.toLowerCase();
    const significanceFilter = document.getElementById('pathway-significance-filter').value;
    const nesFilter = parseFloat(document.getElementById('pathway-nes-filter').value) || 0;

    filteredPathwayData = currentPathwayData.filter(row => {
        // Search filter
        if (searchTerm && !row.pathway.toLowerCase().includes(searchTerm)) {
            return false;
        }

        // Significance filter
        if (significanceFilter === 'significant' && row.padj >= 0.05) {
            return false;
        }

        // NES filter
        if (nesFilter > 0 && Math.abs(row.NES) < nesFilter) {
            return false;
        }

        return true;
    });

    populatePathwayTable(filteredPathwayData);
    updatePathwayTableInfo(filteredPathwayData.length, currentPathwayData.length);
}

// Update DE table info
function updateDeTableInfo(filteredCount, totalCount) {
    const info = document.getElementById('de-table-info');
    info.textContent = `Showing ${filteredCount} of ${totalCount} genes`;
}

// Update pathway table info
function updatePathwayTableInfo(filteredCount, totalCount) {
    const info = document.getElementById('pathway-table-info');
    info.textContent = `Showing ${filteredCount} of ${totalCount} pathways`;
}

// Export DE table data
async function exportDeTable() {
    if (filteredDeData.length === 0) {
        alert('No data to export. Please load data first.');
        return;
    }

    try {
        const csvContent = convertToCSV(filteredDeData);
        downloadCSV(csvContent, `de_results_${currentContrast}_${new Date().toISOString().split('T')[0]}.csv`);
    } catch (error) {
        console.error('Error exporting DE table:', error);
        alert('Error exporting data. Please try again.');
    }
}

// Export pathway table data
async function exportPathwayTable() {
    if (filteredPathwayData.length === 0) {
        alert('No data to export. Please load data first.');
        return;
    }

    try {
        const csvContent = convertToCSV(filteredPathwayData);
        downloadCSV(csvContent, `pathway_results_${currentContrast}_${new Date().toISOString().split('T')[0]}.csv`);
    } catch (error) {
        console.error('Error exporting pathway table:', error);
        alert('Error exporting data. Please try again.');
    }
}

// Convert data to CSV format
function convertToCSV(data) {
    if (data.length === 0) return '';

    const headers = Object.keys(data[0]);
    const csvRows = [
        headers.join(','),
        ...data.map(row =>
            headers.map(header => {
                const value = row[header];
                // Escape commas and quotes in CSV values
                if (typeof value === 'string' && (value.includes(',') || value.includes('"'))) {
                    return `"${value.replace(/"/g, '""')}"`;
                }
                return value;
            }).join(',')
        )
    ];

    return csvRows.join('\n');
}

// Gene Details Functions
// Show gene details modal
async function showGeneDetails(geneId, geneName, contrast = null) {
    try {
        showLoading();

        const response = await fetch(`/api/genes/details?gene_id=${encodeURIComponent(geneId)}&contrast=${encodeURIComponent(contrast || '')}`);
        const geneDetails = await response.json();

        displayGeneDetails(geneDetails);

        hideLoading();

    } catch (error) {
        console.error('Error loading gene details:', error);
        hideLoading();
        alert('Error loading gene details. Please try again.');
    }
}

// Display gene details in modal
function displayGeneDetails(geneDetails) {
    // Update modal title
    document.getElementById('gene-name-title').textContent = geneDetails.gene_name || geneDetails.gene_id;

    const content = document.getElementById('gene-details-content');

    let html = `
        <div class="row mb-4">
            <div class="col-12">
                <h5><i class="fas fa-info-circle me-2"></i>Gene Summary</h5>
                <div class="row">
                    <div class="col-md-3">
                        <div class="card bg-light">
                            <div class="card-body text-center">
                                <h6>${geneDetails.summary.contrasts_analyzed}</h6>
                                <small class="text-muted">Contrasts Analyzed</small>
                            </div>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="card bg-light">
                            <div class="card-body text-center">
                                <h6>${geneDetails.summary.significant_in_contrasts}</h6>
                                <small class="text-muted">Significant Contrasts</small>
                            </div>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="card bg-light">
                            <div class="card-body text-center">
                                <h6>${geneDetails.summary.max_logfc.toFixed(2)}</h6>
                                <small class="text-muted">Max |LogFC|</small>
                            </div>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="card bg-light">
                            <div class="card-body text-center">
                                <h6>${geneDetails.expression_data ? 'Yes' : 'No'}</h6>
                                <small class="text-muted">Expression Data</small>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    `;

    // Differential Expression Section
    if (Object.keys(geneDetails.differential_expression).length > 0) {
        html += `
            <div class="row mb-4">
                <div class="col-12">
                    <h5><i class="fas fa-chart-line me-2"></i>Differential Expression</h5>
                    <div class="table-responsive">
                        <table class="table table-sm table-striped">
                            <thead class="table-dark">
                                <tr>
                                    <th>Contrast</th>
                                    <th>Log2 FC</th>
                                    <th>P-value</th>
                                    <th>Adj P-value</th>
                                    <th>Base Mean</th>
                                    <th>Significant</th>
                                    <th>Direction</th>
                                </tr>
                            </thead>
                            <tbody>
        `;

        Object.entries(geneDetails.differential_expression).forEach(([contrast, data]) => {
            const directionBadge = data.direction === 'up'
                ? '<span class="badge bg-danger">Up</span>'
                : '<span class="badge bg-primary">Down</span>';

            const significantBadge = data.significant
                ? '<span class="badge bg-success">Yes</span>'
                : '<span class="badge bg-secondary">No</span>';

            html += `
                <tr>
                    <td>${contrast.replace('_vs_', ' vs ')}</td>
                    <td>${data.log2FoldChange.toFixed(3)}</td>
                    <td>${data.pvalue.toExponential(2)}</td>
                    <td>${data.padj.toExponential(2)}</td>
                    <td>${data.baseMean.toFixed(1)}</td>
                    <td>${significantBadge}</td>
                    <td>${directionBadge}</td>
                </tr>
            `;
        });

        html += `
                            </tbody>
                        </table>
                    </div>
                </div>
            </div>
        `;
    }

    // Expression Data Section
    if (geneDetails.expression_data && geneDetails.expression_data.samples) {
        html += `
            <div class="row mb-4">
                <div class="col-12">
                    <h5><i class="fas fa-chart-bar me-2"></i>Expression Profile</h5>
                    <div class="row">
                        <div class="col-md-6">
                            <div class="card">
                                <div class="card-header">
                                    <h6>Expression Statistics</h6>
                                </div>
                                <div class="card-body">
                                    <p><strong>Mean Expression:</strong> ${geneDetails.expression_data.mean_expression.toFixed(2)}</p>
                                    <p><strong>Max Expression:</strong> ${geneDetails.expression_data.max_expression.toFixed(2)}</p>
                                    <p><strong>Min Expression:</strong> ${geneDetails.expression_data.min_expression.toFixed(2)}</p>
                                </div>
                            </div>
                        </div>
                        <div class="col-md-6">
                            <div id="gene-expression-plot" style="height: 200px;"></div>
                        </div>
                    </div>
                </div>
            </div>
        `;

        // Store data for plotting after modal is shown
        setTimeout(() => {
            createGeneExpressionPlot(geneDetails.expression_data);
        }, 100);
    }

    // Pathway Information Section
    if (geneDetails.pathway_info && geneDetails.pathway_info.length > 0) {
        html += `
            <div class="row mb-4">
                <div class="col-12">
                    <h5><i class="fas fa-project-diagram me-2"></i>Pathway Information</h5>
                    <p class="text-muted">This gene appears in ${geneDetails.pathway_info.length} pathway(s) as a leading edge gene</p>
                    <div class="table-responsive">
                        <table class="table table-sm table-striped">
                            <thead class="table-dark">
                                <tr>
                                    <th>Pathway</th>
                                    <th>Contrast</th>
                                    <th>NES</th>
                                    <th>P-value</th>
                                    <th>Pathway Size</th>
                                </tr>
                            </thead>
                            <tbody>
        `;

        geneDetails.pathway_info.forEach(pathway => {
            html += `
                <tr>
                    <td><strong>${pathway.pathway}</strong></td>
                    <td>${pathway.contrast.replace('_vs_', ' vs ')}</td>
                    <td>${pathway.NES.toFixed(3)}</td>
                    <td>${pathway.padj.toExponential(2)}</td>
                    <td>${pathway.size}</td>
                </tr>
            `;
        });

        html += `
                            </tbody>
                        </table>
                    </div>
                </div>
            </div>
        `;
    }

    content.innerHTML = html;

    // Show the modal
    const modal = new bootstrap.Modal(document.getElementById('geneDetailsModal'));
    modal.show();
}

// Create expression profile plot for gene details
function createGeneExpressionPlot(expressionData) {
    const samples = expressionData.samples;
    const values = expressionData.expression_values;

    const trace = {
        x: samples,
        y: values,
        type: 'scatter',
        mode: 'lines+markers',
        name: 'Expression',
        marker: { color: '#1f77b4' },
        line: { color: '#1f77b4' }
    };

    const layout = {
        title: 'Gene Expression Across Samples',
        xaxis: { title: 'Sample' },
        yaxis: { title: 'Expression Level' },
        height: 180,
        margin: { t: 30, r: 20, b: 40, l: 50 }
    };

    Plotly.newPlot('gene-expression-plot', [trace], layout);
}

// Export gene data
async function exportGeneData() {
    const geneName = document.getElementById('gene-name-title').textContent;
    const geneId = geneName; // For now, use name as ID

    try {
        const response = await fetch(`/api/genes/details?gene_id=${encodeURIComponent(geneId)}`);
        const geneDetails = await response.json();

        // Convert to CSV format
        const csvData = convertGeneDetailsToCSV(geneDetails);
        downloadCSV(csvData, `gene_details_${geneId}_${new Date().toISOString().split('T')[0]}.csv`);

    } catch (error) {
        console.error('Error exporting gene data:', error);
        alert('Error exporting gene data. Please try again.');
    }
}

// Convert gene details to CSV format
function convertGeneDetailsToCSV(geneDetails) {
    const rows = [];

    // Add header
    rows.push('Gene Details Export');
    rows.push(`Gene ID,Gene Name,Contrasts Analyzed,Significant Contrasts,Max LogFC`);
    rows.push(`${geneDetails.gene_id},${geneDetails.gene_name},${geneDetails.summary.contrasts_analyzed},${geneDetails.summary.significant_in_contrasts},${geneDetails.summary.max_logfc.toFixed(3)}`);
    rows.push('');

    // Differential Expression data
    if (geneDetails.differential_expression) {
        rows.push('Differential Expression Results');
        rows.push('Contrast,Log2 Fold Change,P-value,Adj P-value,Base Mean,Significant,Direction');

        Object.entries(geneDetails.differential_expression).forEach(([contrast, data]) => {
            rows.push(`${contrast},${data.log2FoldChange.toFixed(3)},${data.pvalue.toExponential(2)},${data.padj.toExponential(2)},${data.baseMean.toFixed(1)},${data.significant},${data.direction}`);
        });
        rows.push('');
    }

    // Expression data
    if (geneDetails.expression_data && geneDetails.expression_data.samples) {
        rows.push('Expression Data');
        rows.push('Sample,Expression Level');
        geneDetails.expression_data.samples.forEach((sample, index) => {
            rows.push(`${sample},${geneDetails.expression_data.expression_values[index].toFixed(2)}`);
        });
        rows.push('');
    }

    // Pathway data
    if (geneDetails.pathway_info && geneDetails.pathway_info.length > 0) {
        rows.push('Pathway Information');
        rows.push('Pathway,Contrast,NES,P-value,Pathway Size');

        geneDetails.pathway_info.forEach(pathway => {
            rows.push(`${pathway.pathway},${pathway.contrast},${pathway.NES.toFixed(3)},${pathway.padj.toExponential(2)},${pathway.size}`);
        });
    }

    return rows.join('\n');
}

// Highlight gene in volcano plot (placeholder for future implementation)
function highlightGeneInVolcano(geneName, logFC, negLogP) {
    console.log(`Would highlight gene: ${geneName} at (${logFC}, ${negLogP})`);
    // This would be implemented to highlight specific points in the volcano plot
}

// Comparison View Functions
// Initialize comparison view by populating contrast checkboxes
function initializeComparisonView() {
    const container = document.getElementById('comparison-contrast-selection');
    const deContrasts = Object.keys(currentResults.differential_expression?.contrasts || {});

    if (deContrasts.length === 0) {
        container.innerHTML = '<div class="col-12"><p class="text-muted">No contrasts available for comparison</p></div>';
        return;
    }

    let checkboxesHtml = '';
    deContrasts.forEach((contrast, index) => {
        if (index % 3 === 0) checkboxesHtml += '<div class="row">';
        checkboxesHtml += `
            <div class="col-md-4">
                <div class="form-check">
                    <input class="form-check-input" type="checkbox" value="${contrast}" id="compare-${contrast}">
                    <label class="form-check-label" for="compare-${contrast}">
                        ${contrast.replace('_vs_', ' vs ')}
                    </label>
                </div>
            </div>
        `;
        if ((index + 1) % 3 === 0 || index === deContrasts.length - 1) checkboxesHtml += '</div>';
    });

    container.innerHTML = checkboxesHtml;
}

// Load comparison view with selected contrasts
async function loadComparisonView() {
    const selectedContrasts = Array.from(document.querySelectorAll('#comparison-contrast-selection input[type="checkbox"]:checked'))
        .map(cb => cb.value);

    if (selectedContrasts.length < 2) {
        alert('Please select at least 2 contrasts to compare');
        return;
    }

    try {
        showLoading();

        // Load data for all selected contrasts
        const contrastData = {};
        for (const contrast of selectedContrasts) {
            const response = await fetch(`/api/de/volcano?contrast=${contrast}&pvalue_threshold=0.05&logfc_threshold=1.0`);
            const data = await response.json();
            contrastData[contrast] = data.volcano_data;
        }

        // Generate comparison analysis
        const comparisonResults = analyzeContrastOverlap(contrastData);

        // Display results
        displayComparisonSummary(comparisonResults, selectedContrasts);
        displayOverlapPlot(comparisonResults);
        displayComparisonVolcanoPlots(contrastData, selectedContrasts);

        document.getElementById('comparison-results').style.display = 'block';

        hideLoading();

    } catch (error) {
        console.error('Error loading comparison view:', error);
        hideLoading();
        alert('Error loading comparison data. Please try again.');
    }
}

// Analyze overlap between contrasts
function analyzeContrastOverlap(contrastData) {
    const allGenes = new Set();
    const geneContrastMap = {};

    // Collect all genes and their data across contrasts
    Object.entries(contrastData).forEach(([contrast, genes]) => {
        genes.forEach(gene => {
            const geneId = gene.gene_id || gene.gene;
            if (!geneContrastMap[geneId]) {
                geneContrastMap[geneId] = {
                    gene: gene.gene,
                    gene_id: geneId,
                    contrasts: {},
                    directions: [],
                    logFCs: [],
                    significant: []
                };
            }

            geneContrastMap[geneId].contrasts[contrast] = {
                log2FoldChange: gene.log2FoldChange,
                padj: gene.padj,
                significant: gene.significant
            };

            if (gene.significant) {
                geneContrastMap[geneId].directions.push(gene.log2FoldChange > 0 ? 'up' : 'down');
                geneContrastMap[geneId].logFCs.push(gene.log2FoldChange);
                geneContrastMap[geneId].significant.push(contrast);
            }
        });
    });

    // Calculate overlap statistics
    const significantGenes = Object.values(geneContrastMap).filter(gene => gene.significant.length > 0);

    // Group genes by number of contrasts they're significant in
    const overlapGroups = {};
    significantGenes.forEach(gene => {
        const count = gene.significant.length;
        if (!overlapGroups[count]) {
            overlapGroups[count] = [];
        }
        overlapGroups[count].push(gene);
    });

    return {
        totalGenes: Object.keys(geneContrastMap).length,
        significantGenes: significantGenes.length,
        overlapGroups,
        geneContrastMap
    };
}

// Display comparison summary
function displayComparisonSummary(results, selectedContrasts) {
    const summaryContainer = document.getElementById('comparison-summary');

    const summaryHtml = `
        <div class="col-md-3">
            <div class="card bg-primary text-white">
                <div class="card-body">
                    <h5 class="card-title">${selectedContrasts.length}</h5>
                    <p class="card-text">Contrasts Compared</p>
                </div>
            </div>
        </div>
        <div class="col-md-3">
            <div class="card bg-success text-white">
                <div class="card-body">
                    <h5 class="card-title">${results.significantGenes}</h5>
                    <p class="card-text">Total Significant Genes</p>
                </div>
            </div>
        </div>
        <div class="col-md-3">
            <div class="card bg-warning text-white">
                <div class="card-body">
                    <h5 class="card-title">${Object.keys(results.overlapGroups).length}</h5>
                    <p class="card-text">Overlap Categories</p>
                </div>
            </div>
        </div>
        <div class="col-md-3">
            <div class="card bg-info text-white">
                <div class="card-body">
                    <h5 class="card-title">${results.overlapGroups[Object.keys(results.overlapGroups).length]?.length || 0}</h5>
                    <p class="card-text">Shared Across All</p>
                </div>
            </div>
        </div>
    `;

    summaryContainer.innerHTML = summaryHtml;
}

// Display overlap visualization
function displayOverlapPlot(results) {
    const overlapData = Object.entries(results.overlapGroups).map(([count, genes]) => ({
        count: parseInt(count),
        genes: genes.length
    })).sort((a, b) => b.count - a.count);

    const trace = {
        x: overlapData.map(d => `${d.count} contrast${d.count > 1 ? 's' : ''}`),
        y: overlapData.map(d => d.genes),
        type: 'bar',
        marker: { color: '#1f77b4' }
    };

    const layout = {
        title: 'Gene Overlap Across Contrasts',
        xaxis: { title: 'Number of Contrasts' },
        yaxis: { title: 'Number of Genes' },
        height: 250
    };

    Plotly.newPlot('overlap-summary-plot', [trace], layout);
}

// Display individual volcano plots for comparison
function displayComparisonVolcanoPlots(contrastData, selectedContrasts) {
    const plotsContainer = document.getElementById('comparison-volcano-plots');
    plotsContainer.innerHTML = '';

    selectedContrasts.forEach((contrast, index) => {
        const genes = contrastData[contrast];

        const significant = genes.filter(g => g.significant);
        const upregulated = significant.filter(g => g.direction === 'up');
        const downregulated = significant.filter(g => g.direction === 'down');

        const traces = [];

        // Non-significant points
        const nonSignificant = genes.filter(g => !g.significant);
        if (nonSignificant.length > 0) {
            traces.push({
                x: nonSignificant.map(g => g.log2FoldChange),
                y: nonSignificant.map(g => -Math.log10(g.padj)),
                mode: 'markers',
                type: 'scatter',
                name: 'Not Significant',
                marker: { color: '#cccccc', size: 3, opacity: 0.5 },
                showlegend: index === 0 // Only show legend for first plot
            });
        }

        // Significant points
        if (upregulated.length > 0) {
            traces.push({
                x: upregulated.map(g => g.log2FoldChange),
                y: upregulated.map(g => -Math.log10(g.padj)),
                mode: 'markers',
                type: 'scatter',
                name: 'Upregulated',
                marker: { color: '#d62728', size: 5 },
                showlegend: index === 0
            });
        }

        if (downregulated.length > 0) {
            traces.push({
                x: downregulated.map(g => g.log2FoldChange),
                y: downregulated.map(g => -Math.log10(g.padj)),
                mode: 'markers',
                type: 'scatter',
                name: 'Downregulated',
                marker: { color: '#1f77b4', size: 5 },
                showlegend: index === 0
            });
        }

        const layout = {
            title: contrast.replace('_vs_', ' vs '),
            xaxis: { title: 'Log2 Fold Change' },
            yaxis: { title: '-Log10(Adjusted P-value)' },
            height: 300,
            margin: { t: 50, r: 20, b: 40, l: 50 }
        };

        const plotDiv = document.createElement('div');
        plotDiv.className = 'col-md-6 mb-4';
        plotDiv.innerHTML = `<div class="card"><div class="card-body"><div id="comparison-volcano-${index}" style="height: 250px;"></div></div></div>`;
        plotsContainer.appendChild(plotDiv);

        Plotly.newPlot(`comparison-volcano-${index}`, traces, layout);
    });
}

// Toggle overlap table visibility
function toggleOverlapTable() {
    const tableContainer = document.getElementById('overlap-table-container');
    const button = document.querySelector('[onclick="toggleOverlapTable()"]');

    if (tableContainer.style.display === 'none') {
        tableContainer.style.display = 'block';
        button.innerHTML = '<i class="fas fa-table me-1"></i>Hide Table';
    } else {
        tableContainer.style.display = 'none';
        button.innerHTML = '<i class="fas fa-table me-1"></i>View Table';
    }
}

// Clear comparison view
function clearComparisonView() {
    document.getElementById('comparison-results').style.display = 'none';

    // Uncheck all contrast checkboxes
    document.querySelectorAll('#comparison-contrast-selection input[type="checkbox"]').forEach(cb => {
        cb.checked = false;
    });
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

// Progress Tracking Functions
function initializeProgressTracking() {
    // Check if we're on the dashboard page (has progress tracking section)
    const progressSection = document.getElementById('progress-tracking');
    if (!progressSection) return;

    // Load progress from localStorage or API
    const progress = loadSetupProgress();

    // Update progress circles
    updateProgressCircles(progress);

    // Set up periodic progress checks
    setInterval(checkProgressUpdates, 5000);
}

function loadSetupProgress() {
    // Check localStorage for saved progress
    const savedProgress = localStorage.getItem('rnaseq-setup-progress');
    if (savedProgress) {
        try {
            return JSON.parse(savedProgress);
        } catch (e) {
            console.warn('Invalid progress data in localStorage');
        }
    }

    // Default progress based on configuration files
    const progress = {
        environment: checkEnvironmentSetup() ? 100 : 0,
        config: checkConfigFiles() ? 100 : 0,
        samples: checkSampleFiles() ? 100 : 0,
        analysis: checkAnalysisReady() ? 100 : 0
    };

    return progress;
}

function checkEnvironmentSetup() {
    // Check if conda environments exist or if Docker is available
    return true; // Simplified for now - would check actual environment status
}

function checkConfigFiles() {
    // Check if config files exist and are valid
    const configFiles = ['config/params.yaml', 'config/samples.tsv', 'config/contrasts.tsv'];
    return configFiles.every(file => {
        // This would check if files exist and are valid YAML/TSV
        return true; // Simplified
    });
}

function checkSampleFiles() {
    // Check if sample files are configured and exist
    return true; // Simplified
}

function checkAnalysisReady() {
    // Check if analysis can be run (all prerequisites met)
    return false; // Simplified
}

function updateProgressCircles(progress) {
    const circles = document.querySelectorAll('.progress-circle');
    circles.forEach((circle, index) => {
        const stepKeys = ['environment', 'config', 'samples', 'analysis'];
        const stepProgress = progress[stepKeys[index]] || 0;

        circle.setAttribute('data-progress', stepProgress);
        circle.style.setProperty('--progress', `${stepProgress}%`);

        // Update visual state
        if (stepProgress === 100) {
            circle.classList.add('completed');
        } else if (stepProgress > 0) {
            circle.classList.add('in-progress');
        }
    });
}

function checkProgressUpdates() {
    // Periodically check for progress updates (could be from API)
    const updatedProgress = loadSetupProgress();
    updateProgressCircles(updatedProgress);
    saveSetupProgress(updatedProgress);
}

function saveSetupProgress(progress) {
    localStorage.setItem('rnaseq-setup-progress', JSON.stringify(progress));
}

function showSetupWizard() {
    // Open setup wizard in new tab or modal
    window.open('/tutorial', '_blank');
}

// Tutorial Functions
function startTutorial() {
    // Redirect to tutorial page
    window.location.href = '/tutorial';
}

function completeTutorialStep(step) {
    // Mark tutorial step as completed
    const tutorialProgress = JSON.parse(localStorage.getItem('tutorial-progress') || '{}');
    tutorialProgress[step] = true;
    localStorage.setItem('tutorial-progress', JSON.stringify(tutorialProgress));

    // Update UI if needed
    updateTutorialProgress();
}

function updateTutorialProgress() {
    // Update tutorial progress indicators
    const tutorialProgress = JSON.parse(localStorage.getItem('tutorial-progress') || '{}');
    const dots = document.querySelectorAll('.tutorial-progress-dot');

    dots.forEach((dot, index) => {
        const step = index + 1;
        if (tutorialProgress[step]) {
            dot.classList.add('completed');
        } else {
            dot.classList.remove('completed');
        }
    });
}

// Configuration Management Functions
function saveCurrentConfiguration() {
    const configName = document.getElementById('config-name').value.trim();
    if (!configName) {
        alert('Please enter a name for the configuration');
        return;
    }

    // Save current dashboard state
    const config = {
        name: configName,
        timestamp: new Date().toISOString(),
        state: {
            currentContrast: currentContrast,
            filters: getCurrentFilters(),
            plots: getCurrentPlotSettings()
        }
    };

    // Save to localStorage or send to API
    saveConfigurationToStorage(config);

    // Update UI
    loadConfigurationsFromStorage();
    document.getElementById('config-name').value = '';
}

function getCurrentFilters() {
    return {
        pvalueThreshold: document.getElementById('pvalue-threshold').value,
        logfcThreshold: document.getElementById('logfc-threshold').value,
        pathwayThreshold: document.getElementById('pathway-threshold').value
    };
}

function getCurrentPlotSettings() {
    return {
        volcanoPlot: volcanoPlot ? volcanoPlot.layout : null,
        heatmapPlot: heatmapPlot ? heatmapPlot.layout : null,
        pathwayPlot: pathwayPlot ? pathwayPlot.layout : null
    };
}

function saveConfigurationToStorage(config) {
    const configs = JSON.parse(localStorage.getItem('saved-configurations') || '[]');
    configs.push(config);
    localStorage.setItem('saved-configurations', JSON.stringify(configs));
}

function loadConfigurationsFromStorage() {
    const configs = JSON.parse(localStorage.getItem('saved-configurations') || '[]');
    const select = document.getElementById('load-config-select');

    // Clear existing options
    select.innerHTML = '<option value="">Select a saved configuration...</option>';

    // Add saved configurations
    configs.forEach((config, index) => {
        const option = document.createElement('option');
        option.value = index;
        option.textContent = `${config.name} (${new Date(config.timestamp).toLocaleDateString()})`;
        select.appendChild(option);
    });

    // Show configurations section if we have any
    if (configs.length > 0) {
        document.getElementById('saved-configurations').style.display = 'block';
    }
}

function loadConfiguration() {
    const select = document.getElementById('load-config-select');
    const configIndex = select.value;

    if (!configIndex) return;

    const configs = JSON.parse(localStorage.getItem('saved-configurations') || '[]');
    const config = configs[configIndex];

    if (config) {
        // Restore state
        restoreConfigurationState(config.state);
    }
}

function restoreConfigurationState(state) {
    // Restore filters
    if (state.filters) {
        if (state.filters.pvalueThreshold) {
            document.getElementById('pvalue-threshold').value = state.filters.pvalueThreshold;
        }
        if (state.filters.logfcThreshold) {
            document.getElementById('logfc-threshold').value = state.filters.logfcThreshold;
        }
        if (state.filters.pathwayThreshold) {
            document.getElementById('pathway-threshold').value = state.filters.pathwayThreshold;
        }
    }

    // Restore contrast selection
    if (state.currentContrast) {
        currentContrast = state.currentContrast;
        const deSelect = document.getElementById('de-contrast-select');
        const pathwaySelect = document.getElementById('pathway-contrast-select');

        if (deSelect.querySelector(`option[value="${state.currentContrast}"]`)) {
            deSelect.value = state.currentContrast;
        }
        if (pathwaySelect.querySelector(`option[value="${state.currentContrast}"]`)) {
            pathwaySelect.value = state.currentContrast;
        }
    }

    // Update plots
    if (state.currentContrast) {
        updateVolcanoPlot();
        updateHeatmap();
        updatePathwayPlot();
    }
}

function deleteConfiguration() {
    const select = document.getElementById('load-config-select');
    const configIndex = select.value;

    if (!configIndex) return;

    if (confirm('Are you sure you want to delete this configuration?')) {
        const configs = JSON.parse(localStorage.getItem('saved-configurations') || '[]');
        configs.splice(configIndex, 1);
        localStorage.setItem('saved-configurations', JSON.stringify(configs));

        loadConfigurationsFromStorage();
    }
}

function generateShareUrl() {
    // Generate shareable URL with current state
    const state = {
        contrast: currentContrast,
        filters: getCurrentFilters(),
        plots: getCurrentPlotSettings()
    };

    const params = new URLSearchParams();
    params.set('state', btoa(JSON.stringify(state)));

    const shareUrl = `${window.location.origin}${window.location.pathname}?${params.toString()}`;
    document.getElementById('share-url').value = shareUrl;
}

function copyShareUrl() {
    const shareUrlInput = document.getElementById('share-url');
    shareUrlInput.select();
    document.execCommand('copy');

    // Show feedback
    const originalPlaceholder = shareUrlInput.placeholder;
    shareUrlInput.placeholder = 'URL copied to clipboard!';
    setTimeout(() => {
        shareUrlInput.placeholder = originalPlaceholder;
    }, 2000);
}

function restoreStateFromUrl() {
    // Restore state from URL parameters
    const urlParams = new URLSearchParams(window.location.search);
    const stateParam = urlParams.get('state');

    if (stateParam) {
        try {
            const state = JSON.parse(atob(stateParam));
            restoreConfigurationState(state);
        } catch (error) {
            console.warn('Failed to restore state from URL:', error);
        }
    }
}
