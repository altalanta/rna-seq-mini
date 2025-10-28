"""
Single-cell RNA-seq analysis rules for RNASEQ-MINI.
Supports multiple quantification methods and analysis workflows.
"""

from pathlib import Path

# Single-cell configuration
SINGLECELL_ENABLED = config.get("singlecell", {}).get("enabled", False)
SINGLECELL_TECH = config.get("singlecell", {}).get("technology", "auto")
SINGLECELL_METHOD = config.get("singlecell", {}).get("quantification", {}).get("method", "auto")

# Output directories
SINGLECELL_DIR = Path(config.get("singlecell", {}).get("paths", {}).get("quantification", "results/singlecell/quantification"))
SINGLECELL_QC_DIR = Path(config.get("singlecell", {}).get("paths", {}).get("qc", "results/singlecell/qc"))
SINGLECELL_CLUSTER_DIR = Path(config.get("singlecell", {}).get("paths", {}).get("clustering", "results/singlecell/clustering"))

# Create directories
for directory in [SINGLECELL_DIR, SINGLECELL_QC_DIR, SINGLECELL_CLUSTER_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# Detect single-cell samples (samples with multiple FASTQ files in non-standard naming)
def detect_singlecell_samples():
    """Detect samples that appear to be single-cell based on file patterns."""
    singlecell_samples = []

    for sample in SAMPLES:
        fastq_info = FASTQ_LOOKUP[sample]
        # Look for patterns indicating single-cell data
        # This is a heuristic - could be improved with better metadata
        if any("10x" in str(fq).lower() or "chromium" in str(fq).lower()
               for fq in fastq_info.values() if fq):
            singlecell_samples.append(sample)

    return singlecell_samples

SINGLECELL_SAMPLES = detect_singlecell_samples() if SINGLECELL_ENABLED else []

# Rule to check if single-cell analysis should run
rule singlecell_check:
    """Check if single-cell analysis is enabled and samples are available."""
    output:
        flag = touch(SINGLECELL_DIR / ".singlecell_enabled")
    run:
        if not SINGLECELL_ENABLED:
            shell("echo 'Single-cell analysis disabled' > {output.flag}")
        elif not SINGLECELL_SAMPLES:
            shell("echo 'No single-cell samples detected' > {output.flag}")
        else:
            shell("echo 'Single-cell analysis enabled for samples: {SINGLECELL_SAMPLES}' > {output.flag}")

# Single-cell quantification rules
if SINGLECELL_ENABLED and SINGLECELL_SAMPLES:

    # CellRanger quantification (10x Genomics)
    if SINGLECELL_METHOD in ["auto", "cellranger"]:
        rule singlecell_cellranger_count:
            """Run CellRanger count for 10x Genomics data."""
            input:
                fastq_r1 = lambda wc: FASTQ_LOOKUP[wc.sample]["R1"] if FASTQ_LOOKUP[wc.sample].get("R1") else [],
                fastq_r2 = lambda wc: FASTQ_LOOKUP[wc.sample]["R2"] if FASTQ_LOOKUP[wc.sample].get("R2") else [],
                reference = config["reference"]["transcripts_fa"].replace("transcripts.fa.gz", "cellranger_reference")
            output:
                matrix = SINGLECELL_DIR / "{sample}" / "outs" / "filtered_feature_bc_matrix.h5",
                metrics = SINGLECELL_DIR / "{sample}" / "outs" / "metrics_summary.csv"
            params:
                sample_id = "{sample}",
                output_dir = str(SINGLECELL_DIR / "{sample}"),
                chemistry = config.get("singlecell", {}).get("chemistry", "auto"),
                threads = config.get("singlecell", {}).get("threads", 8),
                memory = config.get("singlecell", {}).get("memory_gb", 32)
            conda:
                "../envs/singlecell.yml"
            shell:
                """
                cellranger count \
                    --id={params.sample_id} \
                    --transcriptome={input.reference} \
                    --fastqs=$(dirname {input.fastq_r1}) \
                    --sample={params.sample_id} \
                    --chemistry={params.chemistry} \
                    --localcores={params.threads} \
                    --localmem={params.memory} \
                    --output-dir={params.output_dir}
                """

    # Kallisto|bustools quantification
    if SINGLECELL_METHOD in ["auto", "kallisto"]:
        rule singlecell_kallisto_bus:
            """Run kallisto bus for single-cell quantification."""
            input:
                fastq_r1 = lambda wc: FASTQ_LOOKUP[wc.sample]["R1"],
                fastq_r2 = lambda wc: FASTQ_LOOKUP[wc.sample]["R2"],
                index = config["reference"]["salmon_index"],
                t2g = config["reference"]["transcripts_fa"].replace("transcripts.fa.gz", "t2g.txt")
            output:
                busfile = SINGLECELL_DIR / "{sample}" / "output.bus"
            params:
                output_dir = str(SINGLECELL_DIR / "{sample}"),
                threads = config.get("singlecell", {}).get("threads", 8)
            conda:
                "../envs/singlecell.yml"
            shell:
                """
                kallisto bus \
                    -i {input.index} \
                    -o {params.output_dir} \
                    -x {SINGLECELL_TECH} \
                    -t {params.threads} \
                    {input.fastq_r1} {input.fastq_r2}
                """

        rule singlecell_bustools_count:
            """Convert kallisto bus output to count matrix."""
            input:
                busfile = SINGLECELL_DIR / "{sample}" / "output.bus",
                t2g = config["reference"]["transcripts_fa"].replace("transcripts.fa.gz", "t2g.txt")
            output:
                matrix = SINGLECELL_DIR / "{sample}" / "counts.mtx",
                barcodes = SINGLECELL_DIR / "{sample}" / "barcodes.txt",
                genes = SINGLECELL_DIR / "{sample}" / "genes.txt"
            params:
                output_dir = str(SINGLECELL_DIR / "{sample}")
            conda:
                "../envs/singlecell.yml"
            shell:
                """
                bustools count \
                    -o {params.output_dir}/counts \
                    -g {input.t2g} \
                    -e {params.output_dir}/matrix.ec \
                    -t {params.output_dir}/transcripts.txt \
                    {input.busfile}
                """

    # Quality control for single-cell data
    rule singlecell_qc:
        """Run quality control on single-cell data."""
        input:
            matrix = SINGLECELL_DIR / "{sample}" / "counts.mtx",
            barcodes = SINGLECELL_DIR / "{sample}" / "barcodes.txt"
        output:
            qc_report = SINGLECELL_QC_DIR / "{sample}" / "qc_report.html",
            filtered_matrix = SINGLECELL_QC_DIR / "{sample}" / "filtered_counts.mtx"
        params:
            output_dir = str(SINGLECELL_QC_DIR / "{sample}"),
            min_genes = config.get("singlecell", {}).get("qc", {}).get("min_genes_per_cell", 200),
            max_genes = config.get("singlecell", {}).get("qc", {}).get("max_genes_per_cell", 6000),
            mito_threshold = config.get("singlecell", {}).get("qc", {}).get("mito_percent_threshold", 10.0)
        conda:
            "../envs/singlecell.yml"
        script:
            "../../scripts/singlecell_qc.py"

    # Single-cell clustering and visualization
    rule singlecell_clustering:
        """Run clustering analysis on single-cell data."""
        input:
            matrix = SINGLECELL_QC_DIR / "{sample}" / "filtered_counts.mtx",
            barcodes = SINGLECELL_DIR / "{sample}" / "barcodes.txt",
            genes = SINGLECELL_DIR / "{sample}" / "genes.txt"
        output:
            clustering_results = SINGLECELL_CLUSTER_DIR / "{sample}" / "clustering_results.rds",
            umap_plot = SINGLECELL_CLUSTER_DIR / "{sample}" / "umap_clusters.png",
            tsne_plot = SINGLECELL_CLUSTER_DIR / "{sample}" / "tsne_clusters.png"
        params:
            output_dir = str(SINGLECELL_CLUSTER_DIR / "{sample}"),
            normalization = config.get("singlecell", {}).get("clustering", {}).get("normalization_method", "log_normalize"),
            n_pcs = config.get("singlecell", {}).get("clustering", {}).get("n_pcs", 50),
            resolution = config.get("singlecell", {}).get("clustering", {}).get("resolution", 1.0),
            method = config.get("singlecell", {}).get("clustering", {}).get("clustering_method", "leiden")
        conda:
            "../envs/singlecell.yml"
        script:
            "../../scripts/singlecell_clustering.py"

# Single-cell analysis completion marker
rule singlecell_complete:
    """Mark single-cell analysis as complete."""
    input:
        flag = SINGLECELL_DIR / ".singlecell_enabled",
        matrices = expand(SINGLECELL_QC_DIR / "{sample}" / "filtered_counts.mtx", sample=SINGLECELL_SAMPLES) if SINGLECELL_SAMPLES else [],
        clustering = expand(SINGLECELL_CLUSTER_DIR / "{sample}" / "clustering_results.rds", sample=SINGLECELL_SAMPLES) if SINGLECELL_SAMPLES else []
    output:
        complete = touch(SINGLECELL_DIR / ".singlecell_complete")
    run:
        if SINGLECELL_SAMPLES:
            print(f"Single-cell analysis completed for samples: {SINGLECELL_SAMPLES}")
        else:
            print("No single-cell samples processed")








