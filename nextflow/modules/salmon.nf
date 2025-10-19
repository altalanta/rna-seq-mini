process SALMON_INDEX {
    tag "salmon_index"
    conda "${projectDir}/envs/salmon.yml"
    publishDir "${params.reference.salmon_index}", mode: 'copy', overwrite: true
    cpus params.salmon?.threads ?: 4

    input:
        tuple path(transcripts), path(annotation), path(decoy, optional: true)

    output:
        path "salmon_index", emit: index
        path "tx2gene.tsv", emit: tx2gene

    script:
        """
        set -euo pipefail
        mkdir -p salmon_index
        decoy_flag=""
        if [ -s "${decoy}" ]; then
          decoy_flag="-d ${decoy}"
        fi
        bash ${projectDir}/scripts/build_salmon_index.sh \
          -t ${transcripts} \
          -a ${annotation} \
          ${decoy_flag} \
          -o salmon_index \
          -p ${task.cpus}
        cp salmon_index/tx2gene.tsv tx2gene.tsv
        """
}

process SALMON_QUANT {
    tag "${sample}"
    conda "${projectDir}/envs/salmon.yml"
    publishDir { "${params.paths.salmon}/${sample}" }, mode: 'copy', overwrite: true
    cpus params.salmon?.threads ?: 4

    input:
        tuple val(sample), val(meta), path(fastq1), path(fastq2, optional: true)
        path index_dir

    output:
        tuple val(sample), path("quant.sf"), path("lib_format_counts.json")

    script:
        """
        set -euo pipefail
        mkdir -p quant
        if [ -s "${fastq2}" ] && [ "${params.se}" != "true" ]; then
          salmon quant -i ${index_dir} --libType ${params.salmon?.libtype ?: 'A'} \
            -1 ${fastq1} -2 ${fastq2} --threads ${task.cpus} ${params.salmon?.extra ?: ''} \
            -o quant
        else
          salmon quant -i ${index_dir} --libType ${params.salmon?.libtype ?: 'A'} \
            -r ${fastq1} --threads ${task.cpus} ${params.salmon?.extra ?: ''} \
            -o quant
        fi
        mv quant/quant.sf .
        mv quant/lib_format_counts.json .
        """
}
