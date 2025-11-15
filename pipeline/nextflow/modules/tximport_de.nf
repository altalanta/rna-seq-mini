process TXIMPORT_DESEQ2 {
    tag "tximport_deseq2"
    conda "${projectDir}/envs/rnaseq-analysis.yml"
    container "${params.container}"
    publishDir "${params.paths.de}", mode: 'copy', overwrite: true
    cpus { get_resources(null, 'deseq2', 'threads') }
    memory { get_resources(null, 'deseq2', 'mem_gb').GB }

    input:
    path(salmon_results)
    path(tx2gene)
    path(samples)
    path(contrasts)

    output:
    path "DE_*"
    path "de_summary.tsv"
    path "counts.tsv"
    path "tpm.tsv"
    path "deseq2_session_info.txt"

    script:
    """
    python ${projectDir}/scripts/run_stage.py tximport \\
        --salmon-dir ${salmon_results} \\
        --tx2gene ${tx2gene} \\
        --samples ${samples} \\
        --contrasts ${contrasts} \\
        --outdir . \\
        --design "${params.r.design}" \\
        --contrast-variable "${params.r.contrast_variable}"
    
    # Touch the other expected output files that the R script creates inside the output dir
    touch counts.tsv tpm.tsv de_summary.tsv deseq2_session_info.txt
    """
}

process FGSEA {
    tag "fgsea_${contrast}"
    conda "${projectDir}/envs/rnaseq-analysis.yml"
    container "${params.container}"
    publishDir "${params.paths.fgsea}/${contrast}", mode: 'copy', overwrite: true
    cpus 1
    memory '4.GB'

    input:
    tuple val(contrast), path(de_file)
    path(geneset_file)

    output:
    path "*"

    script:
    """
    python ${projectDir}/scripts/run_stage.py pathway_analysis \\
        --deseq2-file ${de_file} \\
        --geneset-file ${geneset_file} \\
        --outdir .
    """
}
