# 🚀 Quick Start Guide

Get your first RNA-seq analysis running in **5 minutes** with this streamlined setup.

## Prerequisites

- **Conda/Mamba** installed (for environment management)
- **Basic familiarity** with command line

## 1. Clone and Setup (2 minutes)

```bash
# Clone the repository
git clone <repository-url>
cd rnaseq-mini

# Create analysis environments (downloads ~2GB of software)
make setup

# Verify installation
make smoke
```

**What just happened?** Created isolated environments with all required bioinformatics tools (Salmon, DESeq2, etc.).

## 2. Configure Your Data (2 minutes)

Edit three key files in the `config/` directory:

### `config/samples.tsv` - Tell us about your samples
```tsv
sample	condition	fastq_1	fastq_2
sample1	treated	data/sample1_R1.fastq.gz	data/sample1_R2.fastq.gz
sample2	treated	data/sample2_R1.fastq.gz	data/sample2_R2.fastq.gz
sample3	control	data/sample3_R1.fastq.gz	data/sample3_R2.fastq.gz
sample4	control	data/sample4_R2.fastq.gz	data/sample4_R2.fastq.gz
```

### `config/contrasts.tsv` - Define comparisons
```tsv
groupA	groupB
treated	control
```

### `config/params.yaml` - Set basic parameters
```yaml
project: "my-first-analysis"
organism: "human"  # or "mouse", "yeast", etc.
threads: 8  # adjust based on your computer
```

## 3. Run Analysis (1 minute)

```bash
# Start the analysis (takes 10-30 minutes depending on data size)
make run
```

**What happens next?**
- Quality control with FastQC
- Expression quantification with Salmon
- Differential expression analysis with DESeq2
- Pathway enrichment analysis
- HTML report generation

## 4. Explore Results (immediate)

```bash
# Launch interactive web dashboard
make serve-results

# Or view the static HTML report
open results/report.html
```

## 🎯 Expected Results

Your `results/` directory will contain:
- `qc/` - Quality control reports
- `salmon/` - Expression quantification results
- `de/` - Differential expression tables and plots
- `fgsea/` - Pathway enrichment results
- `report.html` - Complete analysis summary

## 🚨 Common First-Time Issues

**Problem**: "Command not found" errors
**Solution**: Run `conda activate rnaseq-mini-base` first

**Problem**: Analysis fails with "reference not found"
**Solution**: Check that your organism is supported in `config/genome.yaml`

**Problem**: Slow performance
**Solution**: Reduce `threads` in `config/params.yaml` to match your CPU cores

## 📚 Next Steps

- **[Standard Workflow Guide](workflow.md)** - For comprehensive analysis options
- **[Advanced Features](advanced.md)** - Multi-omics, single-cell, cloud deployment
- **[Troubleshooting Guide](troubleshooting.md)** - Solutions to common problems
- **[Configuration Reference](configuration.md)** - Detailed parameter explanations

## 💡 Pro Tips

- **Start small**: Use the included test data first (`tests/data/fastq/`)
- **Monitor progress**: Run `make monitor` for real-time updates
- **Iterate faster**: Enable caching with `cache: enabled: true` in params.yaml
- **Get help**: Check the troubleshooting guide for common solutions

**Ready to dive deeper?** → [Standard Workflow Guide](workflow.md)





















