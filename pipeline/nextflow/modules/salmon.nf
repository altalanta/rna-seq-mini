process SALMON_INDEX {
    tag "salmon_index"
    conda "${projectDir}/envs/rnaseq-core.yml"
    container "${params.container}"
    publishDir "${params.reference.salmon_index}", mode: 'copy', overwrite: true
    cpus 8
    memory '32.GB'

    input:
    path(transcripts)
    path(annotation)
    path(decoy)

    output:
    path "."

    script:
    def decoy_arg = decoy ? "-d ${decoy}" : ""
    """
    bash ${projectDir}/scripts/build_salmon_index.sh \\
        -t ${transcripts} \\
        -a ${annotation} \\
        \${decoy_arg} \\
        -o . \\
        -p ${task.cpus}
    """
}

process SALMON_QUANT {
    tag "${sample_id}"
    conda "${projectDir}/envs/rnaseq-core.yml"
    container "${params.container}"
    publishDir "${params.paths.salmon}/${sample_id}", mode: 'copy', overwrite: true
    cpus { get_resources(reads[0], 'salmon_quant', 'threads') }
    memory { get_resources(reads[0], 'salmon_quant', 'mem_gb').GB }

    input:
    tuple val(sample_id), path(reads)
    path(index)

    output:
    path "."

    script:
    def fastq2_arg = reads.size() > 1 ? "--fastq2 ${reads[1]}" : ""
    """
    python ${projectDir}/scripts/run_stage.py salmon_quant \\
        --index ${index} \\
        --outdir . \\
        --libtype "${params.salmon.libtype}" \\
        --threads ${task.cpus} \\
        --fastq1 ${reads[0]} \\
        \${fastq2_arg} \\
        --extra-opts "${params.salmon.extra}"
    """
}
