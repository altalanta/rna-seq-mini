process REPORT {
    tag "report"
    conda "${projectDir}/envs/r.yml"
    publishDir "${params.paths.outdir}", mode: 'copy', overwrite: true
    cpus 2

    input:
        path multiqc_html
        path counts_dir
        path de_dir
        path fgsea_dir

    output:
        path "results/report.html"

    script:
        """
        set -euo pipefail
        mkdir -p results/qc/multiqc results/counts results/de results/fgsea
        cp ${multiqc_html} results/qc/multiqc/
        cp -r ${counts_dir}/* results/counts/ || true
        cp -r ${de_dir}/* results/de/ || true
        cp -r ${fgsea_dir}/* results/fgsea/ || true
        ${projectDir}/scripts/render_report.R \
          --params ${projectDir}/config/params.yaml \
          --template ${projectDir}/report/rnaseq_report.Rmd \
          --output results/report.html \
          --results-dir results
        """
}
