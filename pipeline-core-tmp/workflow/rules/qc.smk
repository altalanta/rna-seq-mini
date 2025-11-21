rule fastqc:
    input:
        fastq=lambda wildcards: FASTQ_LOOKUP[wildcards.sample][wildcards.read]
    output:
        # The output of the script runner is the directory. We just need to touch a sentinel file.
        outdir=directory(QC_FASTQC_DIR),
        html=lambda wildcards: QC_FASTQC_DIR / f"{wildcards.sample}_{wildcards.read}_fastqc.html"
    log:
        lambda wildcards: LOG_DIR / "fastqc" / f"{wildcards.sample}_{wildcards.read}.log"
    threads: lambda wildcards: get_resources(wildcards, "fastqc", "threads")
    resources:
        mem_gb=lambda wildcards: get_resources(wildcards, "fastqc", "mem_gb")
    conda: "../../envs/rnaseq-core.yml"
    shell:
        """
        python scripts/run_stage.py fastqc \\
            --fastq {input.fastq} \\
            --outdir {output.outdir} \\
            --threads {threads} > {log} 2>&1
        """


rule multiqc:
    input:
        fastqc_reports=expand(QC_FASTQC_DIR / "{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"] if SEQUENCING_MODE == "paired" else ["R1"])
    output:
        report=QC_MULTIQC_DIR / "multiqc_report.html"
    log:
        LOG_DIR / "multiqc" / "multiqc.log"
    threads: 2
    resources:
        mem_mb=2000
    conda: "../../envs/rnaseq-core.yml"
    shell:
        """
        python scripts/run_stage.py multiqc \\
            --analysis-dir {QC_FASTQC_DIR} \\
            --title "{config[multiqc][title]}" \\
            --outdir {QC_MULTIQC_DIR} > {log} 2>&1
        """
