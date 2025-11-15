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

rule run_fgsea:
    input:
        de_file=DE_DIR / "DE_{contrast}.tsv",
        geneset_file=config["paths"]["genesets"]
    output:
        results_table=FGSEA_DIR / "{contrast}" / "fgsea_results.tsv",
        plot_table=FGSEA_DIR / "{contrast}" / "top_pathways_table.pdf",
        plot_enrichment=FGSEA_DIR / "{contrast}" / "top_enrichment_plot.pdf"
    params:
        outdir=lambda wildcards: FGSEA_DIR / wildcards.contrast
    log:
        LOG_DIR / "fgsea" / "{contrast}.log"
    threads: 1
    resources:
        mem_gb=4
    conda:
        "../../envs/rnaseq-analysis.yml"
    run:
        run_command(f"""
            python scripts/run_stage.py pathway_analysis \\
                --deseq2-file {input.de_file} \\
                --geneset-file {input.geneset_file} \\
                --outdir {params.outdir}
        """)

# --- Aggregation Rule ---
rule all_de:
    input:
        de_tables=DE_TABLE_PATHS
    output:
        de_summary=DE_DIR / "de_summary.tsv"
    log:
        LOG_DIR / "de_summary.log"
    run:
        run_command(f"""
            python scripts/run_stage.py aggregate_de \\
                --de-dir {DE_DIR} \\
                --outdir {DE_DIR}
        """)
