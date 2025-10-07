# RNA-seq Mini Pipeline

Dual-engine (Snakemake & Nextflow) RNA-seq workflow that performs QC, transcript quantification, differential expression, pathway analysis, and report generation from FASTQ files to a final HTML report. The repository ships with a tiny yeast dataset for smoke tests and CI, plus reference scaffolding for human-scale runs.

## Features
- Paired/single-end Illumina support with decoy-aware Salmon indexing
- FastQC + MultiQC, Salmon quantification, DESeq2 + tximport, fgsea pathway analysis
- Deterministic R Markdown report summarising QC, DE, and pathway results
- Reproducible environments via Conda/Mamba and optional Docker image
- Turnkey execution on local machines or Slurm clusters with shared params file
- CI smoke test exercises both engines using tiny synthetic reads
- **Enhanced validation**: Cross-engine reproducibility checking with numerical validation
- **Intelligent reference management**: Auto-download references for any species
- **Advanced configuration**: Interactive wizards and comprehensive validation
- **Resource-aware execution**: Automatic resource estimation and optimization
- **Real-time monitoring**: Live progress tracking with quality metrics dashboard

## Enhanced Features

### 🔍 Cross-Engine Validation
Ensure both Snakemake and Nextflow produce identical results:
```bash
make validate                    # Compare outputs between engines
python scripts/validate_determinism.py  # Detailed validation report
```

### 📥 Intelligent Reference Management
Auto-download references for any species:
```bash
python scripts/download_references.py --list    # List available species
python scripts/download_references.py human     # Download human references
make download-refs SPECIES=mouse               # Download with Makefile
```

### ⚙️ Advanced Configuration
Interactive setup and comprehensive validation:
```bash
make wizard                    # Interactive configuration wizard
make validate-full             # Comprehensive validation checks
python scripts/validate_config.py --comprehensive  # Detailed validation
```

### 📊 Resource-Aware Execution
Optimize resource allocation automatically:
```bash
make estimate                  # Show resource recommendations
make optimize                  # Generate optimized configuration
python scripts/estimate_resources.py --report-only  # Resource analysis only
```

### 📈 Real-Time Monitoring
Track pipeline progress with live dashboard:
```bash
make monitor                   # Start real-time monitoring
make monitor-once             # Show current status once
python scripts/monitor_progress.py  # Custom monitoring options
```

## Quickstart

### 1. Clone and set up environments
```bash
mamba env create -f envs/base.yml
mamba env create -f envs/qc.yml
mamba env create -f envs/salmon.yml
mamba env create -f envs/r.yml
```
(Optional) build the Docker image:
```bash
docker build -t rnaseq-mini -f containers/Dockerfile .
```

### 2. Configure inputs
Edit `config/params.yaml`, `config/samples.tsv`, and `config/contrasts.tsv`:
- **samples.tsv** must contain `sample`, `condition`, `fastq_1`, optional `fastq_2`, and optional `batch`.
- **contrasts.tsv** lists two columns (`groupA`, `groupB`) describing pairwise comparisons for DESeq2/fgsea.
- **params.yaml** controls execution (threads, organism presets, Salmon settings, R design formula, report metadata). Organism presets live in `config/genome.yaml`.

### 3. Run the pipeline

#### Snakemake (local)
```bash
snakemake -s pipeline/snakemake/Snakefile --configfile config/params.yaml --use-conda --cores 8
```

#### Snakemake (Slurm)
```bash
snakemake -s pipeline/snakemake/Snakefile \
  --configfile config/params.yaml \
  --profile config/profiles/slurm.smk.yaml
```

#### Nextflow (local)
```bash
nextflow run pipeline/nextflow/main.nf -params-file config/params.yaml -with-conda -profile local
```

#### Nextflow (Slurm)
```bash
nextflow run pipeline/nextflow/main.nf \
  -params-file config/params.yaml \
  -profile slurm \
  -with-conda
```

### Outputs
All results land under `results/`:
- `qc/fastqc/` individual FastQC reports
- `qc/multiqc/multiqc_report.html`
- `salmon/<sample>/quant.sf` expression estimates
- `counts/` gene-level counts (`counts.tsv`, `tpm.tsv`, `txi.rds`)
- `de/` DESeq2 outputs (per-contrast tables, `de_summary.tsv`, plots)
- `fgsea/` pathway tables and enrichment plots
- `report.html` final analysis summary

### Validation & Quality Control
The enhanced validation system provides:
- **Cross-engine reproducibility**: Ensures identical results between Snakemake and Nextflow
- **Numerical validation**: Validates that statistical results are within tolerance
- **Configuration validation**: Comprehensive checks for common setup errors
- **Resource validation**: Ensures resource requirements match available hardware
- **Progress monitoring**: Real-time tracking of pipeline execution and quality metrics

## Methods Summary
1. **QC:** FastQC per FASTQ, aggregated by MultiQC.
2. **Quantification:** Optional Salmon decoy index (`scripts/build_salmon_index.sh`), `salmon quant` with bias correction, stored per sample.
3. **Gene counts + DE:** `scripts/tximport_deseq2.R` performs tximport aggregation, DESeq2 using formula from params, generates per-contrast statistics and plots.
4. **Pathways:** `scripts/fgsea_pathways.R` runs fgsea on log2FC rankings for each contrast.
5. **Reporting:** `report/rnaseq_report.Rmd` composes HTML narrative with QC links, DE tables, pathway hits, and session info.

## Extending
- Update `config/genome.yaml` with new organism presets (FASTA, GTF, decoy, pre-built Salmon index).
- Modify `params.yaml` to change DE design (e.g., include covariates in the formula) or adjust fgsea gene sets and thresholds.
- Add user-specific gene set files (GMT/TSV) and reference in `params.fgsea.genesets`.
- Swap to containerized execution via `-with-docker rnaseq-mini` (Nextflow) or `--use-singularity` / `--use-apptainer` in Snakemake.
- Add new species to `scripts/download_references.py` for automatic reference management.
- Customize resource estimation in `scripts/estimate_resources.py` for specialized hardware.

## Troubleshooting

### Common Issues
- **Configuration validation errors**: Run `make validate-full` to identify and fix issues
- **Resource allocation problems**: Use `make estimate` to optimize resource settings
- **Cross-engine differences**: Run `make validate` to ensure both engines produce identical results
- **Missing references**: Use `make download-refs SPECIES=<organism>` to auto-download
- **Pipeline monitoring**: Run `make monitor` to track progress in real-time

### Performance Tips
- Use `make optimize` to generate resource-optimized configuration
- Monitor resource usage with `make monitor` during execution
- For large datasets (>50 samples), consider HPC clusters with Slurm profiles
- Enable `auto_download: true` in `config/params.yaml` for automatic reference management

## Continuous Integration
- **Main CI** (`.github/workflows/ci.yml`): Runs smoke tests for both engines and enhanced validation
- **Enhanced Validation** (`.github/workflows/enhanced_validation.yml`): Comprehensive testing of all new features

Both workflows ensure cross-engine reproducibility, validate configurations, and test all enhanced features including reference management, resource estimation, and progress monitoring.
