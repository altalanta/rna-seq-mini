# RNA-seq Mini Pipeline

Dual-engine (Snakemake & Nextflow) RNA-seq workflow that performs QC, transcript quantification, differential expression, pathway analysis, and report generation from FASTQ files to a final HTML report. The repository ships with a tiny yeast dataset for smoke tests and CI, plus reference scaffolding for human-scale runs.

## Features
- Paired/single-end Illumina support with decoy-aware Salmon indexing
- FastQC + MultiQC, Salmon quantification, DESeq2 + tximport, fgsea pathway analysis
- Deterministic R Markdown report summarising QC, DE, and pathway results
- Reproducible environments via Conda/Mamba and optional Docker image
- Turnkey execution on local machines or Slurm clusters with shared params file
- CI smoke test exercises both engines using tiny synthetic reads

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

## Continuous Integration
GitHub Actions (`.github/workflows/ci.yml`) runs the `tests/run_smoke.sh` script in a matrix over Snakemake and Nextflow, ensuring both engines process the bundled yeast dataset end-to-end and emit the expected artifacts.
