rule run_deseq2:
    input:
        salmon_dir=SALMON_DIR,
        samples=config["paths"]["samples"],
        contrasts=config["r"]["contrasts_file"],
        tx2gene=rules.salmon_index.output.tx2gene
    output:
        counts_table=COUNTS_DIR / "counts.tsv",
        tpm_table=COUNTS_DIR / "tpm.tsv",
        de_summary=DE_DIR / "de_summary.tsv",
        r_env=DE_DIR / "deseq2_session_info.txt",
        de_tables=DE_TABLE_PATHS
    log:
        LOG_DIR / "deseq2" / "deseq2.log"
    threads: lambda wildcards: get_resources(wildcards, "deseq2", "threads")
    resources:
        mem_gb=lambda wildcards: get_resources(wildcards, "deseq2", "mem_gb")
    conda:
        "../../envs/rnaseq-analysis.yml"
    shell:
        """
        python scripts/run_stage.py tximport \\
            --salmon-dir {input.salmon_dir} \\
            --tx2gene {input.tx2gene} \\
            --samples {input.samples} \\
            --contrasts {input.contrasts} \\
            --outdir {DE_DIR} \\
            --design "{config[r][design]}" \\
            --contrast-variable "{config[r][contrast_variable]}" > {log} 2>&1
        touch {output.counts_table} {output.tpm_table} {output.de_summary} {output.r_env}
        """


rule fgsea:
    input:
        de_tables=DE_TABLE_PATHS,
        contrasts=config["r"]["contrasts_file"]
    output:
        summary=FGSEA_DIR / "fgsea_summary.tsv",
        fgsea_tables=FGSEA_TABLE_PATHS
    log:
        LOG_DIR / "r" / "fgsea.log"
    params:
        figdir=FGSEA_FIG_DIR
    threads: 2
    resources:
        mem_mb=4000
    conda: "../../envs/r.yml"
    shell:
        """
        mkdir -p {LOG_DIR / "r"}
        scripts/fgsea_pathways.R \
          --params config/params.yaml \
          --de-dir {DE_DIR} \
          --contrast-file {input.contrasts} \
          --outdir {FGSEA_DIR} \
          --figdir {params.figdir} > {log} 2>&1
        """
