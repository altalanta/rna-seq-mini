process RENDER_REPORT {
    tag "render_report"
    conda "${projectDir}/envs/rnaseq-analysis.yml"
    container "${params.container}"
    publishDir "${params.paths.report_dir}", mode: 'copy', overwrite: true

    input:
    path(de_results)
    path(samples)
    path(r_env)

    output:
    path "report.html"

    script:
    """
    python ${projectDir}/scripts/run_stage.py render_report \\
        --config "${projectDir}/config/params.yaml" \\
        --samples ${samples} \\
        --r-env ${r_env}
    """
}
