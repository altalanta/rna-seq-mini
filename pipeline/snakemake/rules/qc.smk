rule fastqc:
    input:
        fastq=lambda wildcards: FASTQ_LOOKUP[wildcards.sample][wildcards.read]
    output:
        html=lambda wildcards: QC_FASTQC_DIR / f"{wildcards.sample}_{wildcards.read}_fastqc.html",
        zip=lambda wildcards: QC_FASTQC_DIR / f"{wildcards.sample}_{wildcards.read}_fastqc.zip"
    log:
        lambda wildcards: LOG_DIR / "fastqc" / f"{wildcards.sample}_{wildcards.read}.log"
    threads: 2
    resources:
        mem_mb=2000
    conda: "../../envs/qc.yml"
    shell:
        """
        mkdir -p {LOG_DIR / "fastqc"}
        fastqc {config["fastqc"]["extra"]} --threads {threads} --outdir {QC_FASTQC_DIR} {input.fastq} > {log} 2>&1
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
    conda: "../../envs/qc.yml"
    shell:
        """
        mkdir -p {LOG_DIR / "multiqc"}
        multiqc {QC_FASTQC_DIR} --filename multiqc_report.html --outdir {QC_MULTIQC_DIR} > {log} 2>&1
        """
