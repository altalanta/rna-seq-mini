rule tximport_deseq2:
    input:
        quant=expand(SALMON_DIR / "{sample}" / "quant.sf", sample=SAMPLES),
        tx2gene=rules.salmon_index.output.tx2gene,
        samples=config["paths"]["samples"],
        contrasts=config["r"]["contrasts_file"]
    output:
        counts=COUNTS_DIR / "counts.tsv",
        tpm=COUNTS_DIR / "tpm.tsv",
        txi=COUNTS_DIR / "txi.rds",
        dds=DE_DIR / "dds.rds",
        summary=DE_DIR / "de_summary.tsv",
        normalized=DE_DIR / "normalized_counts.tsv",
        de_tables=DE_TABLE_PATHS
    log:
        LOG_DIR / "r" / "tximport_deseq2.log"
    params:
        figdir=DE_FIG_DIR,
        quant_dir=SALMON_DIR,
        se_flag="--se" if SEQUENCING_MODE == "single" else ""
    threads: config["threads"]
    resources:
        mem_mb=config["memory_gb"] * 1024
    conda: "../../envs/r.yml"
    shell:
        """
        mkdir -p {LOG_DIR / "r"}
        scripts/tximport_deseq2.R \
          --params config/params.yaml \
          --sample-sheet {input.samples} \
          --quant-dir {params.quant_dir} \
          --tx2gene {input.tx2gene} \
          --out-counts {output.counts} \
          --out-tpm {output.tpm} \
          --out-txi {output.txi} \
          --out-dds {output.dds} \
          --out-summary {output.summary} \
          --contrast-file {input.contrasts} \
          --figdir {params.figdir} \
          {params.se_flag} > {log} 2>&1
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
