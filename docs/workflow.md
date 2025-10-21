# 🔬 Standard Workflow Guide

Comprehensive RNA-seq analysis workflow with detailed explanations and best practices.

## Overview

This guide covers the complete analysis pipeline from raw FASTQ files to publication-ready results. Each step includes quality checkpoints and optimization recommendations.

## 1. Data Preparation & Quality Control

### Input Requirements

**File Format**: FASTQ files (`.fastq.gz` compressed)
**Naming Convention**: `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz` for paired-end data

**Sample Metadata**: Ensure your `samples.tsv` includes:
- `sample`: Unique sample identifier
- `condition`: Experimental condition (e.g., "treated", "control")
- `fastq_1`: Path to forward reads
- `fastq_2`: Path to reverse reads (omit for single-end)
- `batch`: Optional batch identifier for batch correction

### Quality Control Pipeline

The pipeline automatically runs:

1. **FastQC** - Per-sample quality assessment
2. **MultiQC** - Aggregated quality report across all samples

**Quality Gates**: Automatic checks for:
- Minimum read count per sample (>100K reads recommended)
- Per-base sequence quality (Phred score >28 recommended)
- Adapter contamination (<5% recommended)

### Interpreting QC Results

```bash
# View aggregated quality report
open results/qc/multiqc_report.html

# Check for problematic samples
make assess-quality
```

**Warning Signs**:
- ⚠️ Low read counts (< 1M reads)
- ⚠️ High adapter content (>10%)
- ⚠️ Poor per-base quality (median Phred < 25)

## 2. Expression Quantification

### Salmon Quantification

**Key Parameters**:
```yaml
salmon:
  libtype: "A"  # Auto-detect for most Illumina data
  extra: "--validateMappings --gcBias"  # Recommended for accuracy
  threads: 4  # Adjust based on available cores
```

**What Salmon Does**:
- Pseudo-alignment to transcriptome
- Bias correction for GC content and positional bias
- Quantification uncertainty estimation
- Gene-level and transcript-level estimates

### Output Files
- `results/salmon/<sample>/quant.sf` - Transcript-level quantification
- `results/counts/counts.tsv` - Gene-level counts matrix
- `results/counts/tpm.tsv` - TPM normalized expression

## 3. Differential Expression Analysis

### DESeq2 Analysis

**Design Formula** (configure in `params.yaml`):
```yaml
r:
  design: "~ condition + batch"  # Includes batch correction
  # Alternative: "~ condition" for simple two-group comparison
```

**Key Parameters**:
```yaml
r:
  alpha: 0.05  # Significance threshold
  lfc_shrink: true  # Shrink log fold changes for better ranking
  pvalue_adjust: "BH"  # Benjamini-Hochberg multiple testing correction
```

### Output Interpretation

**Essential Results Files**:
- `results/de/deseq2_results.tsv` - Complete DE results table
- `results/de/volcano_plot.png` - Visual summary of DE results
- `results/de/pca_plot.png` - Sample clustering visualization

**Key Columns to Examine**:
- `padj` - Adjusted p-value (use <0.05 as threshold)
- `log2FoldChange` - Effect size (typically |LFC| > 1 for biological relevance)
- `baseMean` - Expression level indicator

### Advanced DESeq2 Features

```bash
# Enable batch correction optimization
make optimize-batch

# Perform power analysis
make power-analysis

# Cross-validation for result stability
make cross-validation
```

## 4. Pathway Enrichment Analysis

### fgsea Configuration

**Gene Set Collections**: Configure in `params.yaml`:
```yaml
fgsea:
  genesets: "config/custom_genesets.gmt"  # Custom gene sets
  # Built-in options: "hallmark", "kegg", "reactome", "go_bp"
  min_size: 15  # Minimum genes per pathway
  max_size: 500  # Maximum genes per pathway
  nperm: 1000  # Permutations for significance testing
```

### Interpreting Pathway Results

**Key Metrics**:
- `NES` (Normalized Enrichment Score) - Effect size
- `padj` - Multiple testing corrected p-value
- `leadingEdge` - Genes driving the enrichment

**Output Files**:
- `results/fgsea/fgsea_results.tsv` - Complete enrichment results
- `results/fgsea/enrichment_plots/` - Visual pathway representations

## 5. Quality Assessment & Validation

### Comprehensive Quality Check

```bash
# Run full quality assessment
make assess-quality

# Benchmark against reference datasets
make benchmark-analysis

# Evaluate publication readiness
make quality-gate
```

### Quality Metrics

**QC Assessment**:
- Read quality distribution
- Mapping rates and coverage
- Batch effect detection

**DE Quality**:
- Effect size distribution
- Significance rate validation
- Replicate correlation

**Count Quality**:
- Expression range assessment
- Missing data analysis
- Library size normalization

## 6. Reporting & Visualization

### HTML Report

The final report (`results/report.html`) includes:
- **QC Dashboard** - Links to FastQC and MultiQC reports
- **DE Results** - Interactive tables and volcano plots
- **Pathway Analysis** - Enrichment results and visualizations
- **Methods Summary** - Complete parameter documentation
- **Session Info** - Reproducible environment details

### Interactive Web Dashboard

```bash
# Launch interactive exploration
make serve-results
```

**Features**:
- Real-time filtering of DE results
- Interactive volcano plots with dynamic thresholds
- Export capabilities for custom result sets
- Responsive design for mobile access

## 7. Best Practices & Optimization

### Performance Optimization

```bash
# Analyze your dataset characteristics
make estimate

# Generate optimized configuration
make optimize

# Enable intelligent caching
cache:
  enabled: true
  max_age_days: 30
```

### Reproducibility

- **Version Control**: Track all configuration files
- **Environment Isolation**: Use conda environments
- **Seed Setting**: Set random seeds for reproducible results
- **Documentation**: Record all analysis decisions

### Troubleshooting Common Issues

See [Troubleshooting Guide](troubleshooting.md) for solutions to:
- Memory allocation errors
- Reference genome issues
- Statistical power problems
- Visualization rendering issues

## 📚 Related Documentation

- **[Quick Start](quickstart.md)** - If you're new to the pipeline
- **[Advanced Features](advanced.md)** - Enterprise features and multi-omics
- **[Configuration Reference](configuration.md)** - Detailed parameter documentation
- **[API Guide](api.md)** - Programmatic access and integration





