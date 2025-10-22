/*
Single-cell RNA-seq analysis module for RNASEQ-MINI Nextflow pipeline.
Supports multiple quantification methods and analysis workflows.
*/

params.singlecell_enabled = false
params.singlecell_technology = "auto"
params.singlecell_method = "auto"
params.singlecell_min_genes = 200
params.singlecell_max_genes = 6000
params.singlecell_mito_threshold = 10.0

// Single-cell quantification processes
process SINGLECELL_CELLRANGER_COUNT {
    tag "${sample_id}"
    label 'process_high'

    conda "${projectDir}/envs/singlecell.yml"

    input:
    tuple val(sample_id), path(reads)
    path(reference)

    output:
    path("${sample_id}/outs"), emit: cellranger_outs
    path("${sample_id}/outs/filtered_feature_bc_matrix.h5"), emit: matrix_h5
    path("${sample_id}/outs/metrics_summary.csv"), emit: metrics

    script:
    def chemistry = params.singlecell_chemistry ?: "auto"
    def memory_gb = params.singlecell_memory_gb ?: 32
    """
    cellranger count \\
        --id=${sample_id} \\
        --transcriptome=${reference} \\
        --fastqs=. \\
        --sample=${sample_id} \\
        --chemistry=${chemistry} \\
        --localcores=${task.cpus} \\
        --localmem=${memory_gb} \\
        --output-dir=${sample_id}
    """
}

process SINGLECELL_KALLISTO_BUS {
    tag "${sample_id}"
    label 'process_medium'

    conda "${projectDir}/envs/singlecell.yml"

    input:
    tuple val(sample_id), path(reads)
    path(index)
    path(t2g)

    output:
    path("${sample_id}/output.bus"), emit: busfile
    path("${sample_id}"), emit: kallisto_dir

    script:
    def technology = params.singlecell_technology ?: "10xV3"
    """
    kallisto bus \\
        -i ${index} \\
        -o ${sample_id} \\
        -x ${technology} \\
        -t ${task.cpus} \\
        ${reads}
    """
}

process SINGLECELL_BUSTOOLS_COUNT {
    tag "${sample_id}"
    label 'process_medium'

    conda "${projectDir}/envs/singlecell.yml"

    input:
    path(busfile)
    path(t2g)
    val(sample_id)

    output:
    path("${sample_id}/counts.mtx"), emit: matrix
    path("${sample_id}/barcodes.txt"), emit: barcodes
    path("${sample_id}/genes.txt"), emit: genes

    script:
    """
    mkdir -p ${sample_id}
    bustools count \\
        -o ${sample_id}/counts \\
        -g ${t2g} \\
        -e ${sample_id}/matrix.ec \\
        -t ${sample_id}/transcripts.txt \\
        ${busfile}

    # Move outputs to expected locations
    mv ${sample_id}/counts.mtx ${sample_id}/counts.mtx
    mv ${sample_id}/counts.barcodes.txt ${sample_id}/barcodes.txt
    mv ${sample_id}/counts.genes.txt ${sample_id}/genes.txt
    """
}

process SINGLECELL_QC {
    tag "${sample_id}"
    label 'process_medium'

    conda "${projectDir}/envs/singlecell.yml"

    input:
    path(matrix)
    path(barcodes)
    path(genes)
    val(sample_id)

    output:
    path("${sample_id}/qc_report.html"), emit: qc_report
    path("${sample_id}/filtered_counts.mtx"), emit: filtered_matrix
    path("${sample_id}"), emit: qc_dir

    script:
    def min_genes = params.singlecell_min_genes ?: 200
    def max_genes = params.singlecell_max_genes ?: 6000
    def mito_threshold = params.singlecell_mito_threshold ?: 10.0
    """
    # Run single-cell QC script
    singlecell_qc.py \\
        --matrix ${matrix} \\
        --barcodes ${barcodes} \\
        --genes ${genes} \\
        --output-dir ${sample_id} \\
        --min-genes ${min_genes} \\
        --max-genes ${max_genes} \\
        --mito-threshold ${mito_threshold}
    """
}

process SINGLECELL_CLUSTERING {
    tag "${sample_id}"
    label 'process_medium'

    conda "${projectDir}/envs/singlecell.yml"

    input:
    path(matrix)
    path(barcodes)
    path(genes)
    val(sample_id)

    output:
    path("${sample_id}/clustering_results.rds"), emit: clustering_results
    path("${sample_id}/umap_clusters.png"), emit: umap_plot
    path("${sample_id}/tsne_clusters.png"), emit: tsne_plot
    path("${sample_id}"), emit: clustering_dir

    script:
    def normalization = params.singlecell_clustering_normalization ?: "log_normalize"
    def n_pcs = params.singlecell_clustering_n_pcs ?: 50
    def resolution = params.singlecell_clustering_resolution ?: 1.0
    def method = params.singlecell_clustering_method ?: "leiden"
    """
    # Run single-cell clustering script
    singlecell_clustering.py \\
        --matrix ${matrix} \\
        --barcodes ${barcodes} \\
        --genes ${genes} \\
        --output-dir ${sample_id} \\
        --normalization ${normalization} \\
        --n-pcs ${n_pcs} \\
        --resolution ${resolution} \\
        --method ${method}
    """
}

// Single-cell workflow
workflow SINGLECELL_ANALYSIS {
    take:
    samples_ch
    reference_index
    t2g_file

    main:
    // Detect single-cell samples and prepare reads
    samples_ch
        .filter { sample_id, reads ->
            params.singlecell_enabled &&
            (sample_id.contains('10x') || sample_id.contains('chromium') ||
             params.singlecell_technology != "auto")
        }
        .set { singlecell_samples }

    if (params.singlecell_method == "cellranger" || params.singlecell_method == "auto") {
        // CellRanger workflow
        reference_cellranger = file("${reference_index}/cellranger_reference")

        SINGLECELL_CELLRANGER_COUNT(singlecell_samples, reference_cellranger)

        // Convert CellRanger output to standard format for downstream analysis
        SINGLECELL_CELLRANGER_COUNT.out.matrix_h5
            .map { matrix_h5 ->
                def sample_id = matrix_h5.getParent().getName()
                return tuple(sample_id, matrix_h5, file("${matrix_h5.getParent()}/outs/metrics_summary.csv"))
            }
            .set { cellranger_results }

        // Convert CellRanger H5 to MTX format for clustering
        CELLRANGER_TO_MTX(cellranger_results)

        // Run QC and clustering
        SINGLECELL_QC(CELLRANGER_TO_MTX.out.matrix,
                     CELLRANGER_TO_MTX.out.barcodes,
                     CELLRANGER_TO_MTX.out.genes,
                     CELLRANGER_TO_MTX.out.sample_id)

        SINGLECELL_CLUSTERING(SINGLECELL_QC.out.filtered_matrix,
                             CELLRANGER_TO_MTX.out.barcodes,
                             CELLRANGER_TO_MTX.out.genes,
                             SINGLECELL_QC.out.sample_id)

    } else if (params.singlecell_method == "kallisto" || params.singlecell_method == "auto") {
        // Kallisto|bustools workflow
        SINGLECELL_KALLISTO_BUS(singlecell_samples, reference_index, t2g_file)

        SINGLECELL_BUSTOOLS_COUNT(SINGLECELL_KALLISTO_BUS.out.busfile,
                                 t2g_file,
                                 SINGLECELL_KALLISTO_BUS.out.sample_id)

        // Run QC and clustering
        SINGLECELL_QC(SINGLECELL_BUSTOOLS_COUNT.out.matrix,
                     SINGLECELL_BUSTOOLS_COUNT.out.barcodes,
                     SINGLECELL_BUSTOOLS_COUNT.out.genes,
                     SINGLECELL_BUSTOOLS_COUNT.out.sample_id)

        SINGLECELL_CLUSTERING(SINGLECELL_QC.out.filtered_matrix,
                             SINGLECELL_BUSTOOLS_COUNT.out.barcodes,
                             SINGLECELL_BUSTOOLS_COUNT.out.genes,
                             SINGLECELL_QC.out.sample_id)
    }

    emit:
    qc_reports = SINGLECELL_QC.out.qc_report
    clustering_results = SINGLECELL_CLUSTERING.out.clustering_results
    umap_plots = SINGLECELL_CLUSTERING.out.umap_plot
    tsne_plots = SINGLECELL_CLUSTERING.out.tsne_plot
}

// Helper process to convert CellRanger H5 to MTX format
process CELLRANGER_TO_MTX {
    tag "${sample_id}"
    label 'process_medium'

    conda "${projectDir}/envs/singlecell.yml"

    input:
    tuple val(sample_id), path(matrix_h5), path(metrics)

    output:
    path("${sample_id}/counts.mtx"), emit: matrix
    path("${sample_id}/barcodes.txt"), emit: barcodes
    path("${sample_id}/genes.txt"), emit: genes
    val(sample_id), emit: sample_id

    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    from scipy.io import mmwrite
    from pathlib import Path

    # Read CellRanger H5 file
    adata = sc.read_10x_h5('${matrix_h5}')

    # Convert to MTX format
    Path('${sample_id}').mkdir(exist_ok=True)

    # Write matrix
    mmwrite('${sample_id}/counts.mtx', adata.X.T)

    # Write barcodes
    adata.obs_names.to_series().to_csv('${sample_id}/barcodes.txt', index=False, header=False)

    # Write genes
    genes_df = pd.DataFrame({
        'gene_id': adata.var_names,
        'gene_symbol': adata.var_names
    })
    genes_df.to_csv('${sample_id}/genes.txt', index=False, header=False)

    print(f"Converted CellRanger data for {sample_id}: {adata.n_obs} cells, {adata.n_vars} genes")
    """
}




