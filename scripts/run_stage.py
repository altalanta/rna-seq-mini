#!/usr/bin/env python3
"""
Centralized script runner for the RNASEQ-MINI pipeline.
This script provides a unified command-line interface for all core pipeline stages,
ensuring that the same logic is used by both Snakemake and Nextflow workflows.
"""

import subprocess
import sys
from pathlib import Path
import click
import yaml

# Add src to the Python path to allow importing modules
sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))
from rnaseq_mini.logger import get_logger
from scripts import cache_manager

log = get_logger(__name__)

# --- Utility Functions ---

def run_command(cmd: list, **kwargs):
    """
    Executes a command, logs it, captures its output, and checks for errors.
    """
    cmd_str = ' '.join(map(str, cmd))
    log.info("executing_command", command=cmd_str)
    
    try:
        # Use PIPE to capture stdout/stderr
        process = subprocess.run(
            cmd, check=True, capture_output=True, text=True, **kwargs
        )
        # Log stdout/stderr from the completed process
        if process.stdout:
            log.debug("command_stdout", command=cmd_str, output=process.stdout.strip())
        if process.stderr:
            log.debug("command_stderr", command=cmd_str, output=process.stderr.strip())

    except subprocess.CalledProcessError as e:
        log.error(
            "command_failed",
            command=cmd_str,
            return_code=e.returncode,
            stdout=e.stdout.strip(),
            stderr=e.stderr.strip(),
        )
        sys.exit(e.returncode)

# --- CLI Definition ---

@click.group()
@click.option("--use-cache/--no-cache", default=True, help="Enable or disable caching for this run.")
@click.pass_context
def cli(ctx, use_cache):
    """A central CLI for running rnaseq-mini pipeline stages."""
    ctx.ensure_object(dict)
    
    # Load the main config to check if caching is globally enabled
    with open("config/params.yaml", 'r') as f:
        config = yaml.safe_load(f)
    
    global_cache_enabled = config.get("cache", {}).get("enabled", True)
    ctx.obj['USE_CACHE'] = global_cache_enabled and use_cache

@cli.command()
@click.option("--fastq", "fastq_path", required=True, type=click.Path(exists=True, path_type=Path), help="Path to the input FASTQ file.")
@click.option("--outdir", "output_dir", required=True, type=click.Path(path_type=Path), help="Output directory for FastQC results.")
@click.option("--threads", default=1, type=int, help="Number of threads to use.")
@click.pass_context
def fastqc(ctx, fastq_path: Path, output_dir: Path, threads: int):
    """Run FastQC on a single FASTQ file."""
    log.info("fastqc_started", input_fastq=str(fastq_path), out_dir=str(output_dir))
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Caching Logic
    if ctx.obj['USE_CACHE']:
        input_files = [fastq_path]
        params = {"threads": threads, "command": "fastqc"}
        cache_hash = cache_manager.compute_hash(input_files, params)
        
        # NOTE: Caching FastQC is tricky because its output is two files (.html, .zip) inside a dir.
        # For simplicity, we will treat the whole directory as the output.
        # A more robust solution might handle individual files.
        output_target = output_dir / f"{fastq_path.stem}_fastqc.html" # Sentinel file
        if cache_manager.check_cache(cache_hash):
            log.info("cache_hit", stage="fastqc", hash=cache_hash)
            # FastQC outputs to a directory, so we can't easily symlink.
            # We'll skip for now, but a real implementation might zip/unzip the output dir.
            log.warn("fastqc_caching_not_fully_implemented", reason="Directory output is complex to cache.")

    cmd = [
        "fastqc",
        "--outdir", output_dir,
        "--threads", threads,
        fastq_path,
    ]
    run_command(cmd)
    
    # Caching Logic: Store result
    # if ctx.obj['USE_CACHE']:
    #     cache_manager.store_in_cache(cache_hash, output_target)

    log.info("fastqc_finished", input_fastq=str(fastq_path))

@cli.command()
@click.option("--analysis-dir", "analysis_dir", required=True, type=click.Path(exists=True, path_type=Path), help="Directory containing analysis files to scan (e.g., FastQC zips).")
@click.option("--outdir", "output_dir", required=True, type=click.Path(path_type=Path), help="Output directory for the MultiQC report.")
@click.option("--title", default="MultiQC Report", help="Title for the MultiQC report.")
@click.pass_context
def multiqc(ctx, analysis_dir: Path, output_dir: Path, title: str):
    """Run MultiQC on a directory of analysis files."""
    log.info("multiqc_started", analysis_dir=str(analysis_dir))
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Caching is complex for multiqc as inputs are a whole directory.
    # We will skip implementing it for this stage for now.

    cmd = [
        "multiqc",
        "--force",
        "--title", title,
        "--outdir", output_dir,
        analysis_dir,
    ]
    run_command(cmd)
    log.info("multiqc_finished", output_dir=str(output_dir))

@cli.command()
@click.option("--index", "index_path", required=True, type=click.Path(exists=True, path_type=Path), help="Path to the Salmon index.")
@click.option("--outdir", "output_dir", required=True, type=click.Path(path_type=Path), help="Output directory for quantification.")
@click.option("--libtype", default="A", help="Library type for Salmon.")
@click.option("--threads", default=4, type=int, help="Number of threads to use.")
@click.option("--fastq1", type=click.Path(exists=True, path_type=Path), help="Path to R1 FASTQ file.")
@click.option("--fastq2", type=click.Path(exists=True, path_type=Path), help="Path to R2 FASTQ file (for paired-end).")
@click.option("--extra-opts", default="", help="Extra options to pass to Salmon.")
@click.pass_context
def salmon_quant(ctx, index_path: Path, output_dir: Path, libtype: str, threads: int, fastq1: Path, fastq2: Path, extra_opts: str):
    """Run Salmon quantification."""
    log.info("salmon_quant_started", sample=fastq1.name, out_dir=str(output_dir))
    output_dir.mkdir(parents=True, exist_ok=True)

    # Caching Logic
    if ctx.obj['USE_CACHE']:
        input_files = [fastq1, fastq2] if fastq2 else [fastq1]
        input_files.append(index_path) # The index is also a key input
        params = {"libtype": libtype, "threads": threads, "extra": extra_opts, "command": "salmon_quant"}
        cache_hash = cache_manager.compute_hash(input_files, params)
        
        # Salmon outputs a directory. We will cache the key file `quant.sf`.
        output_target = output_dir / "quant.sf"
        if cache_manager.check_cache(cache_hash):
            log.info("cache_hit", stage="salmon_quant", hash=cache_hash)
            # This is a simplification. A real implementation would need to handle the whole directory.
            # For now, we assume if quant.sf is cached, the rest is fine.
            cache_manager.retrieve_from_cache(cache_hash, output_target)
            log.info("salmon_quant_finished_from_cache", sample=fastq1.name)
            return

    cmd = [
        "salmon", "quant",
        "--index", index_path,
        "--libType", libtype,
        "--threads", threads,
        "--output", output_dir,
        "--validateMappings",
        "--gcBias",
    ]
    if fastq1 and fastq2:
        cmd.extend(["-1", fastq1, "-2", fastq2])
    elif fastq1:
        cmd.extend(["-r", fastq1])
    
    if extra_opts:
        cmd.extend(extra_opts.split())
        
    run_command(cmd)

    # Caching Logic: Store result
    if ctx.obj['USE_CACHE']:
        cache_manager.store_in_cache(cache_hash, output_target)

    log.info("salmon_quant_finished", sample=fastq1.name)


@cli.command()
@click.option("--deseq2-file", "deseq2_file", required=True, type=click.Path(exists=True, path_type=Path), help="Path to DESeq2 result file.")
@click.option("--geneset-file", "geneset_file", required=True, type=click.Path(exists=True, path_type=Path), help="Path to gene set file (.gmt).")
@click.option("--outdir", "output_dir", required=True, type=click.Path(path_type=Path), help="Output directory for GSEA results.")
@click.pass_context
def pathway_analysis(ctx, deseq2_file: Path, geneset_file: Path, output_dir: Path):
    """Run Gene Set Enrichment Analysis (GSEA) using fgsea."""
    log.info("gsea_started", contrast=deseq2_file.stem)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Caching for GSEA is complex due to multiple outputs and plot generation.
    # We will skip implementing it for this stage for now.

    cmd = [
        "Rscript",
        str(Path(__file__).parent / "run_fgsea.R"),
        "--deseq2_file", deseq2_file,
        "--geneset_file", geneset_file,
        "--outdir", output_dir,
    ]
    run_command(cmd)

    log.info("gsea_finished", contrast=deseq2_file.stem)


@cli.command()
@click.option("--salmon-dir", required=True, type=click.Path(exists=True, path_type=Path), help="Path to the root directory of Salmon outputs.")
@click.option("--tx2gene", required=True, type=click.Path(exists=True, path_type=Path), help="Path to the transcript-to-gene mapping file.")
@click.option("--samples", required=True, type=click.Path(exists=True, path_type=Path), help="Path to the samples.tsv file.")
@click.option("--contrasts", required=True, type=click.Path(exists=True, path_type=Path), help="Path to the contrasts.tsv file.")
@click.option("--outdir", required=True, type=click.Path(path_type=Path), help="Directory to save DESeq2 results.")
@click.option("--design", default="~ condition", help="Design formula for DESeq2.")
@click.option("--contrast-variable", default="condition", help="Variable to use for contrasts.")
def tximport(salmon_dir: Path, tx2gene: Path, samples: Path, contrasts: Path, outdir: Path, design: str, contrast_variable: str):
    """Run tximport and DESeq2 R script."""
    log.info("tximport_deseq2_started", out_dir=str(outdir))
    outdir.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "Rscript", "scripts/tximport_deseq2.R",
        "--salmon_dir", salmon_dir,
        "--tx2gene", tx2gene,
        "--samples", samples,
        "--contrasts", contrasts,
        "--outdir", outdir,
        "--design", f'"{design}"', # Pass formula with quotes
        "--variable", contrast_variable,
    ]
    run_command(cmd)
    log.info("tximport_deseq2_finished", out_dir=str(outdir))

@cli.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True, path_type=Path), help="Path to the main params.yaml config file.")
@click.option("--samples", required=True, type=click.Path(exists=True, path_type=Path), help="Path to the samples.tsv file.")
@click.option("--r-env", required=True, type=click.Path(exists=True, path_type=Path), help="Path to a file capturing the R session info.")
def render_report(config_path: Path, samples: Path, r_env: Path):
    """Render the final HTML report using R Markdown."""
    log.info("render_report_started")
    
    cmd = [
        "Rscript", "scripts/render_report.R",
        "--config", config_path,
        "--samples", samples,
        "--r_env", r_env,
    ]
    run_command(cmd)
    log.info("render_report_finished")

if __name__ == "__main__":
    cli()
