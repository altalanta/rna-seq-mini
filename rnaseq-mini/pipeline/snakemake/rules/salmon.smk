rule salmon_index:
    input:
        transcripts=config["reference"]["transcripts_fa"],
        annotation=config["reference"]["annotation_gtf"]
    output:
        index=directory(str(SALMON_INDEX_DIR)),
        tx2gene=SALMON_INDEX_DIR / "tx2gene.tsv"
    params:
        decoy=config["reference"].get("decoy_fasta")
    log:
        LOG_DIR / "salmon" / "index.log"
    threads: config["salmon"]["threads"]
    resources:
        mem_mb=config["memory_gb"] * 1024
    conda: "../../envs/salmon.yml"
    shell:
        """
        mkdir -p {LOG_DIR / "salmon"}
        DEC=""
        if [ -n "{params.decoy}" ] && [ "{params.decoy}" != "None" ]; then
          DEC="-d {params.decoy}"
        fi
        bash scripts/build_salmon_index.sh \
          -t {input.transcripts} \
          -a {input.annotation} \
          ${DEC} \
          -o {output.index} \
          -p {threads} > {log} 2>&1
        """


rule salmon_quant:
    input:
        fastq1=lambda wildcards: FASTQ_LOOKUP[wildcards.sample]["R1"],
        fastq2=lambda wildcards: FASTQ_LOOKUP[wildcards.sample].get("R2"),
        index=rules.salmon_index.output.index
    output:
        quant=SALMON_DIR / "{sample}" / "quant.sf",
        libjson=SALMON_DIR / "{sample}" / "lib_format_counts.json"
    log:
        lambda wildcards: LOG_DIR / "salmon" / f"{wildcards.sample}.log"
    threads: config["salmon"]["threads"]
    resources:
        mem_mb=config["memory_gb"] * 1024
    conda: "../../envs/salmon.yml"
    shell:
        """
        mkdir -p {SALMON_DIR / wildcards.sample}
        mkdir -p {LOG_DIR / "salmon"}
        if [ -n "{input.fastq2}" ] && [ "{input.fastq2}" != "None" ] && [ "{SEQUENCING_MODE}" = "paired" ]; then \
          salmon quant -i {input.index} --libType {config["salmon"]["libtype"]} \
            -1 {input.fastq1} -2 {input.fastq2} \
            --threads {threads} {config["salmon"]["extra"]} \
            -o {SALMON_DIR / wildcards.sample} > {log} 2>&1; \
        else \
          salmon quant -i {input.index} --libType {config["salmon"]["libtype"]} \
            -r {input.fastq1} --threads {threads} {config["salmon"]["extra"]} \
            -o {SALMON_DIR / wildcards.sample} > {log} 2>&1; \
        fi
        """
