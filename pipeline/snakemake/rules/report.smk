rule report_html:
    input:
        multiqc=QC_MULTIQC_DIR / "multiqc_report.html",
        counts=COUNTS_DIR / "counts.tsv",
        de_summary=DE_DIR / "de_summary.tsv",
        fgsea=FGSEA_DIR / "fgsea_summary.tsv",
        de_tables=DE_TABLE_PATHS,
        fgsea_tables=FGSEA_TABLE_PATHS
    output:
        html=REPORT_HTML
    log:
        LOG_DIR / "r" / "report.log"
    params:
        template="report/rnaseq_report.Rmd",
        results_dir=RESULTS_DIR
    threads: 2
    resources:
        mem_mb=4000
    conda: "../../envs/r.yml"
    shell:
        """
        mkdir -p {LOG_DIR / "r"}
        scripts/render_report.R \
          --params config/params.yaml \
          --template {params.template} \
          --output {output.html} \
          --results-dir {params.results_dir} > {log} 2>&1
        """
