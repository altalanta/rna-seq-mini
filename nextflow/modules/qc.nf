process FASTQC {
    tag "${sample}_${read}"
    conda "${projectDir}/envs/qc.yml"
    publishDir "${params.paths.qc}/fastqc", mode: 'copy', overwrite: true
    cpus 2

    input:
        tuple val(sample), val(read), path(fastq)

    output:
        tuple val(sample), val(read), path("${sample}_${read}_fastqc.zip"), path("${sample}_${read}_fastqc.html")

    script:
        """
        set -euo pipefail
        mkdir -p fastqc
        extra="${params.fastqc?.extra ?: ''}"
        fastqc ${extra} --threads ${task.cpus} --outdir fastqc ${fastq}
        base=$(basename ${fastq})
        base=${base%%.*}
        mv fastqc/${base}_fastqc.zip ${sample}_${read}_fastqc.zip
        mv fastqc/${base}_fastqc.html ${sample}_${read}_fastqc.html
        """
}

process MULTIQC {
    tag "multiqc"
    conda "${projectDir}/envs/qc.yml"
    publishDir "${params.paths.qc}/multiqc", mode: 'copy', overwrite: true
    cpus 2

    input:
        val trigger

    output:
        path "multiqc_report.html"

    script:
        """
        set -euo pipefail
        mkdir -p output
        multiqc ${params.paths.qc}/fastqc --outdir output --filename multiqc_report.html --title "${params.multiqc?.title ?: 'RNA-seq QC'}"
        mv output/multiqc_report.html .
        """
}
