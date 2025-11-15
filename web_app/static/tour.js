document.addEventListener('DOMContentLoaded', function() {
  const tour = new Shepherd.Tour({
    useModalOverlay: true,
    defaultStepOptions: {
      classes: 'shadow-md bg-purple-dark',
      scrollTo: true,
      cancelIcon: {
        enabled: true
      }
    }
  });

  tour.addStep({
    id: 'welcome',
    text: 'Welcome to the RNA-seq Mini dashboard! This tour will guide you through the main features of this interface.',
    buttons: [
      {
        text: 'Next',
        action: tour.next
      }
    ]
  });

  tour.addStep({
    id: 'contrast-selection',
    text: 'Start by selecting an experimental contrast from this dropdown. All plots and tables on this page will update to reflect the selected comparison.',
    attachTo: {
      element: '#contrast-dropdown',
      on: 'bottom'
    },
    buttons: [
      {
        text: 'Back',
        action: tour.back
      },
      {
        text: 'Next',
        action: tour.next
      }
    ]
  });
  
  tour.addStep({
    id: 'volcano-plot',
    text: 'This volcano plot visualizes differential expression. Genes on the right are upregulated, genes on the left are downregulated. The higher a gene is plotted, the more statistically significant the change.',
    attachTo: {
      element: '#volcano-plot',
      on: 'bottom'
    },
    buttons: [
        {
          text: 'Back',
          action: tour.back
        },
        {
          text: 'Next',
          action: tour.next
        }
    ]
  });

  tour.addStep({
    id: 'gene-view',
    text: 'You can search for a specific gene here. Try typing a gene name (e.g., "gene1") and clicking search to see its expression and statistics.',
    attachTo: {
        element: '#gene-input',
        on: 'bottom'
    },
    buttons: [
        {
          text: 'Back',
          action: tour.back
        },
        {
          text: 'Next',
          action: tour.next
        }
    ]
  });

  tour.addStep({
    id: 'pathway-analysis',
    text: 'This tab contains the Gene Set Enrichment Analysis (GSEA) results, which show which biological pathways are enriched in your differentially expressed genes.',
    attachTo: {
      element: 'a[href="#pathway-analysis-tab"]',
      on: 'bottom'
    },
    buttons: [
      {
        text: 'Back',
        action: tour.back
      },
      {
        text: 'Finish',
        action: tour.complete
      }
    ]
  });

  const startTourButton = document.getElementById('start-tour-button');
  if (startTourButton) {
    startTourButton.addEventListener('click', () => {
      tour.start();
    });
  }
});
