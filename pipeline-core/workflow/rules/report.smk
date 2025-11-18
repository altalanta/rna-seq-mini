rule render_report:
    input:
        de_summary=rules.run_deseq2.output.de_summary,
        fgsea_summary=rules.run_fgsea.output.fgsea_summary,
        samples=config["paths"]["samples"],
        r_env=rules.run_deseq2.output.r_env
    output:
        report=REPORT_HTML
    log:
        LOG_DIR / "report" / "render_report.log"
    threads: 1
    resources:
        mem_mb=4000
    conda:
        "../../envs/rnaseq-analysis.yml"
    shell:
        """
        python scripts/run_stage.py render_report \\
            --config "config/params.yaml" \\
            --samples {input.samples} \\
            --r-env {input.r_env} > {log} 2>&1
        """
