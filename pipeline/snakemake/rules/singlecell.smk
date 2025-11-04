"""
Single-cell RNA-seq analysis rules for RNASEQ-MINI.
Supports multiple quantification methods and analysis workflows.
"""

from pathlib import Path

# Single-cell configuration
SINGLECELL_ENABLED = config.get("singlecell", {}).get("enabled", False)
SINGLECELL_INPUT_DIR = Path(config.get("singlecell", {}).get("input_dir", ""))
SINGLECELL_OUTPUT_DIR = Path(config.get("singlecell", {}).get("paths", {}).get("output_dir", "results/singlecell/analysis"))

# Create output directory
if SINGLECELL_ENABLED:
    SINGLECELL_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Main analysis rule that runs the complete Scanpy workflow
if SINGLECELL_ENABLED:
    rule singlecell_scanpy_analysis:
        input:
            # The input is the directory, but we need a file to act as a trigger.
            # We'll use the matrix file as the primary input dependency.
            matrix=SINGLECELL_INPUT_DIR / "matrix.mtx.gz",
            barcodes=SINGLECELL_INPUT_DIR / "barcodes.tsv.gz",
            features=SINGLECELL_INPUT_DIR / "features.tsv.gz"
        output:
            adata=SINGLECELL_OUTPUT_DIR / "processed_data.h5ad",
            umap_plot=SINGLECELL_OUTPUT_DIR / "umap_clusters.png",
            marker_genes=SINGLECELL_OUTPUT_DIR / "marker_genes_top20.tsv",
            qc_plot=SINGLECELL_OUTPUT_DIR / "qc_violin_before_filtering.png"
        params:
            input_dir=str(SINGLECELL_INPUT_DIR),
            output_dir=str(SINGLECELL_OUTPUT_DIR),
            min_genes=config["singlecell"]["analysis"]["min_genes_per_cell"],
            min_cells=config["singlecell"]["analysis"]["min_cells_per_gene"],
            max_mito=config["singlecell"]["analysis"]["max_mito_percent"],
            n_pcs=config["singlecell"]["analysis"]["n_pcs"],
            n_neighbors=config["singlecell"]["analysis"]["n_neighbors"],
            resolution=config["singlecell"]["analysis"]["clustering_resolution"]
        log:
            SINGLECELL_OUTPUT_DIR / "scanpy_analysis.log"
        conda:
            "../../envs/rnaseq-analysis.yml"
        threads: config["threads"]
        shell:
            """
            python scripts/run_scanpy_analysis.py \\
                --input_dir {params.input_dir} \\
                --output_dir {params.output_dir} \\
                --min_genes {params.min_genes} \\
                --min_cells {params.min_cells} \\
                --max_mito_percent {params.max_mito} \\
                --n_pcs {params.n_pcs} \\
                --n_neighbors {params.n_neighbors} \\
                --resolution {params.resolution} > {log} 2>&1
            """

# Single-cell analysis completion marker
rule singlecell_complete:
    """Mark single-cell analysis as complete."""
    input:
        # This rule will run if single-cell is enabled, and it depends on the final output file.
        analysis_outputs=rules.singlecell_scanpy_analysis.output if SINGLECELL_ENABLED else []
    output:
        complete=touch(SINGLECELL_OUTPUT_DIR / ".singlecell_complete")
    run:
        if SINGLECELL_ENABLED:
            print(f"Single-cell analysis completed. Outputs in {SINGLECELL_OUTPUT_DIR}")
        else:
            shell("echo 'Single-cell analysis not enabled.' > {output.complete}")


