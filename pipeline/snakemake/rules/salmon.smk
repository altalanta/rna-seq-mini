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
    run:
        # Check cache before running
        if cache_manager and CACHE_ENABLED:
            input_files = [Path(input.transcripts), Path(input.annotation)]
            if params.decoy and params.decoy != "None":
                input_files.append(Path(params.decoy))

            output_files = [Path(output.index), Path(output.tx2gene)]
            parameters = {
                'decoy': params.decoy,
                'threads': threads,
                'salmon_version': '1.9.0'  # Could be made configurable
            }

            should_skip, cache_key = should_skip_stage(
                cache_manager, "salmon_index", input_files, output_files, parameters
            )

            if should_skip:
                print(f"Using cached Salmon index (key: {cache_key})")
                shell("touch {output.index} {output.tx2gene}")  # Ensure outputs exist
            else:
                print(f"Building Salmon index (key: {cache_key})")
                shell("""
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
                """)

                # Mark as complete in cache
                if cache_key:
                    mark_stage_complete(
                        cache_manager, cache_key, input_files, output_files, parameters,
                        metadata={'stage': 'salmon_index', 'organism': config.get('organism', 'unknown')}
                    )
        else:
            # Run without caching
            shell("""
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
            """)


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
    run:
        # Check cache before running
        if cache_manager and CACHE_ENABLED:
            input_files = [Path(input.fastq1)]
            if input.fastq2 and input.fastq2 != "None":
                input_files.append(Path(input.fastq2))
            input_files.append(Path(input.index))

            output_files = [Path(output.quant), Path(output.libjson)]
            parameters = {
                'libtype': config["salmon"]["libtype"],
                'extra': config["salmon"]["extra"],
                'threads': threads,
                'sequencing_mode': SEQUENCING_MODE
            }

            should_skip, cache_key = should_skip_stage(
                cache_manager, f"salmon_quant_{wildcards.sample}", input_files, output_files, parameters
            )

            if should_skip:
                print(f"Using cached Salmon quantification for {wildcards.sample} (key: {cache_key})")
                shell("touch {output.quant} {output.libjson}")  # Ensure outputs exist
            else:
                print(f"Running Salmon quantification for {wildcards.sample} (key: {cache_key})")
                shell("""
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
                """)

                # Mark as complete in cache
                if cache_key:
                    mark_stage_complete(
                        cache_manager, cache_key, input_files, output_files, parameters,
                        metadata={'stage': 'salmon_quant', 'sample': wildcards.sample}
                    )
        else:
            # Run without caching
            shell("""
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
            """)
