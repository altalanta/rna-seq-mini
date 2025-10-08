process TXIMPORT_DESEQ2 {
    tag "tximport_deseq2"
    conda "${projectDir}/envs/r.yml"
    publishDir "${params.paths.outdir}", mode: 'copy', overwrite: true
    cpus params.threads ?: 4

    input:
        val trigger
        path sample_sheet
        path contrasts
        path tx2gene

    output:
        path "results/counts", emit: counts_dir
        path "results/de", emit: de_dir

    script:
        """
        set -euo pipefail
        mkdir -p results/counts results/de
        ${projectDir}/scripts/tximport_deseq2.R \
          --params ${projectDir}/config/params.yaml \
          --sample-sheet ${sample_sheet} \
          --quant-dir ${params.paths.salmon} \
          --tx2gene ${tx2gene} \
          --out-counts results/counts/counts.tsv \
          --out-tpm results/counts/tpm.tsv \
          --out-txi results/counts/txi.rds \
          --out-dds results/de/dds.rds \
          --out-summary results/de/de_summary.tsv \
          --contrast-file ${contrasts} \
          --figdir results/de/figures ${params.se ? '--se' : ''}
        """
}

process FGSEA {
    tag "fgsea"
    conda "${projectDir}/envs/r.yml"
    publishDir "${params.paths.outdir}", mode: 'copy', overwrite: true
    cpus 2

    input:
        path de_dir
        path contrast_file

    output:
        path "results/fgsea", emit: fgsea_dir

    script:
        """
        set -euo pipefail
        mkdir -p results
        cp -r ${de_dir} results/
        ${projectDir}/scripts/fgsea_pathways.R \
          --params ${projectDir}/config/params.yaml \
          --de-dir results/de \
          --contrast-file ${contrast_file} \
          --outdir results/fgsea \
          --figdir results/fgsea/figures
        """
}
