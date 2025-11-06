process FASTQC {
    tag "${sample}_${read}"
    conda "${projectDir}/envs/rnaseq-core.yml"
    container "${params.container}"
    publishDir "${params.paths.qc}/fastqc", mode: 'copy', overwrite: true
    cpus 2

    input:
        tuple val(sample_id), path(reads)

    output:
        path "fastqc"

    script:
        def read_files = reads.join(' ')
        """
        mkdir -p fastqc
        python ${projectDir}/scripts/run_stage.py fastqc \\
            --fastq ${reads[0]} \\
            --outdir fastqc
        
        if [ -f "${reads[1]}" ]; then
            python ${projectDir}/scripts/run_stage.py fastqc \\
                --fastq ${reads[1]} \\
                --outdir fastqc
        fi
        """
}

process MULTIQC {
    tag "multiqc"
    conda "${projectDir}/envs/rnaseq-core.yml"
    container "${params.container}"
    publishDir "${params.paths.qc}/multiqc", mode: 'copy', overwrite: true
    cpus 2

    input:
        path(fastqc_results)

    output:
        path "multiqc"

    script:
        """
        python ${projectDir}/scripts/run_stage.py multiqc \\
            --analysis-dir . \\
            --title "${params.multiqc_title}" \\
            --outdir "multiqc"
        """
}
