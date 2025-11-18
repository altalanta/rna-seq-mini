document.addEventListener('DOMContentLoaded', function() {

    const tour = new Shepherd.Tour({
        useModalOverlay: true,
        defaultStepOptions: {
            classes: 'shadow-md bg-purple-dark',
            scrollTo: true
        }
    });

    tour.addStep({
        id: 'welcome',
        text: 'Welcome to the RNA-seq Mini Dashboard! This tour will guide you through the key features of this interactive report.',
        buttons: [{
            text: 'Next',
            action: tour.next
        }]
    });

    tour.addStep({
        id: 'contrast-dropdown',
        text: 'First, select the experimental comparison (contrast) you want to explore from this dropdown menu. All plots and tables will update based on your selection.',
        attachTo: {
            element: '#contrast-dropdown',
            on: 'bottom'
        },
        buttons: [{
            text: 'Back',
            action: tour.back
        }, {
            text: 'Next',
            action: tour.next
        }]
    });

    tour.addStep({
        id: 'volcano-plot',
        text: 'This is the volcano plot. It shows the relationship between fold change (x-axis) and statistical significance (y-axis) for every gene. Genes in red are significantly up- or down-regulated.',
        attachTo: {
            element: '#volcano-plot',
            on: 'bottom'
        },
        buttons: [{
            text: 'Back',
            action: tour.back
        }, {
            text: 'Next',
            action: tour.next
        }]
    });

    tour.addStep({
        id: 'de-table',
        text: 'The table below the plot contains the detailed differential expression results. You can sort by any column, or filter by gene name.',
        attachTo: {
            element: '#de-table',
            on: 'top'
        },
        buttons: [{
            text: 'Back',
            action: tour.back
        }, {
            text: 'Next',
            action: tour.next
        }]
    });

    tour.addStep({
        id: 'gene-view-tab',
        text: 'Now, let\'s explore a specific gene. Click on the "Gene View" tab.',
        attachTo: {
            // Dash dynamically creates tab IDs, so we target it by its value
            element: 'li.nav-item a[data-bs-target="#tabs-content .tab-pane:nth-child(2)"]',
            on: 'bottom'
        },
        buttons: [{
            text: 'Back',
            action: tour.back
        }, {
            text: 'Next',
            action: tour.next
        }]
    });
    
    tour.addStep({
        id: 'gene-search',
        text: 'In the "Gene View" tab, you can type a gene ID into this search box and click "Search" to see its expression across all samples and its stats in all comparisons.',
        attachTo: {
            element: '#gene-input',
            on: 'bottom'
        },
        buttons: [{
            text: 'Back',
            action: tour.back
        }, {
            text: 'Next',
            action: tour.next
        }]
    });

    tour.addStep({
        id: 'pathway-tab',
        text: 'The "Pathway Analysis" tab shows the results of Gene Set Enrichment Analysis (GSEA). This helps you understand the biological processes that are impacted by the gene expression changes.',
        attachTo: {
             element: 'li.nav-item a[data-bs-target="#tabs-content .tab-pane:nth-child(3)"]',
            on: 'bottom'
        },
        buttons: [{
            text: 'Back',
            action: tour.back
        }, {
            text: 'Next',
            action: tour.next
        }]
    });

    tour.addStep({
        id: 'qc-tab',
        text: 'Finally, the "QC Summary" tab contains key quality control metrics from the pipeline run, helping you assess the quality of the sequencing data.',
        attachTo: {
             element: 'li.nav-item a[data-bs-target="#tabs-content .tab-pane:nth-child(4)"]',
            on: 'bottom'
        },
        buttons: [{
            text: 'Back',
            action: tour.back
        }, {
            text: 'Finish',
            action: tour.complete
        }]
    });

    const startTourButton = document.getElementById('start-tour-button');
    if (startTourButton) {
        startTourButton.addEventListener('click', () => {
            tour.start();
        });
    }

});


