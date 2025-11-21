# 🔧 Troubleshooting Guide

Common issues and solutions for RNASEQ-MINI.

## Installation Issues

### Conda Environment Problems

**Problem**: `make setup` fails with conda errors
```bash
# Solution 1: Use mamba (faster conda replacement)
conda install mamba -n base -c conda-forge
make setup

# Solution 2: Clean conda state
conda clean --all
make setup

# Solution 3: Create environments manually
mamba env create -f envs/base.yml --force
mamba env create -f envs/qc.yml --force
mamba env create -f envs/salmon.yml --force
mamba env create -f envs/r.yml --force
```

**Problem**: "Conda environment not found" when running commands
```bash
# Activate the base environment first
conda activate rnaseq-mini-base

# Then run commands
make run
```

### Missing Dependencies

**Problem**: Import errors for bioinformatics tools
```bash
# Reinstall environments
make setup

# Or install specific packages
conda install -n rnaseq-mini-base fastqc multiqc salmon

# Check environment contents
conda list -n rnaseq-mini-base
```

## Configuration Issues

### Sample File Problems

**Problem**: "Sample file not found" or format errors
```bash
# Validate sample file format
python scripts/validate_config.py

# Check file exists and has correct headers
head -1 config/samples.tsv
# Expected: sample<tab>condition<tab>fastq_1<tab>fastq_2<tab>batch

# Use setup wizard for guided configuration
python scripts/setup_wizard.py
```

**Problem**: FASTQ files not found
```bash
# Check file paths exist
ls -la data/*.fastq.gz

# Update paths in samples.tsv
# Use absolute paths if relative paths don't work
```

### Reference Genome Issues

**Problem**: "Reference not found" errors
```bash
# Check organism is supported
cat config/genome.yaml

# Download missing references
make download-refs SPECIES=human

# Use auto-download feature
# In params.yaml:
reference:
  auto_download: true
```

**Problem**: Salmon index not found
```bash
# Build custom Salmon index
bash scripts/build_salmon_index.sh \
  -t references/human/transcripts.fa.gz \
  -a references/human/annotation.gtf.gz \
  -d references/human/decoys.fa.gz \
  -o references/human/salmon_index
```

## Runtime Issues

### Memory Problems

**Problem**: Out of memory errors
```bash
# Reduce memory allocation
# In params.yaml:
memory_gb: 16  # Down from 32
threads: 4     # Down from 8

# Use swap file
sudo fallocate -l 16G /swapfile
sudo chmod 600 /swapfile
sudo swapon /swapfile
```

**Problem**: Large dataset memory issues
```bash
# Process samples in batches
# Modify samples.tsv to process smaller groups

# Use streaming mode for large files
salmon:
  extra: "--validateMappings --gcBias --streamMode"
```

### Disk Space Issues

**Problem**: "No space left on device"
```bash
# Clean old results
make clean

# Clear cache
make cache-clear

# Use external storage
# In params.yaml:
paths:
  outdir: "/external/drive/results"
```

**Problem**: Cache taking too much space
```bash
# Check cache size
make cache-stats

# Reduce cache retention
# In params.yaml:
cache:
  max_age_days: 7  # Down from 30

# Clear old cache entries
make cache-cleanup-force
```

### Performance Issues

**Problem**: Analysis running very slowly
```bash
# Increase thread count
# In params.yaml:
threads: 16  # Up from 8

# Enable caching
cache:
  enabled: true

# Use faster storage for intermediate files
cache:
  dir: "/tmp/rnaseq_cache"
```

**Problem**: Salmon quantification slow
```bash
# Optimize Salmon parameters
salmon:
  threads: 8  # Increase from 4
  extra: "--validateMappings --gcBias"

# Use pre-built index if available
reference:
  salmon_index: "path/to/prebuilt_index"
```

## Statistical Issues

### DESeq2 Problems

**Problem**: "Design matrix not full rank" error
```bash
# Fix design formula in params.yaml
r:
  design: "~ condition"  # Remove problematic covariates

# Check for collinear factors
# Remove redundant batch variables

# Use setup wizard for guidance
python scripts/setup_wizard.py
```

**Problem**: Low number of differentially expressed genes
```bash
# Adjust significance threshold
r:
  alpha: 0.1  # Up from 0.05

# Check for batch effects
r:
  design: "~ condition + batch"

# Increase sample size or biological replicates
```

**Problem**: "Dispersion outliers" warnings
```bash
# Use alternative dispersion estimation
# This is usually informational, not an error

# Check for sample quality issues
make assess-quality
```

### Pathway Analysis Issues

**Problem**: fgsea not finding enrichments
```bash
# Check gene set file format
# GMT format: pathway_name<tab>description<tab>gene1<tab>gene2...

# Adjust pathway size limits
fgsea:
  min_size: 10   # Down from 15
  max_size: 1000 # Up from 500

# Use different gene ranking
fgsea:
  score_column: "stat"  # Instead of log2FoldChange
```

## Output Issues

### Report Generation Problems

**Problem**: HTML report not generating
```bash
# Check R environment
conda activate rnaseq-mini-r
R --version

# Install missing R packages
R -e "install.packages(c('rmarkdown', 'ggplot2', 'pheatmap'))"

# Run report generation manually
Rscript scripts/render_report.R
```

**Problem**: Web dashboard not loading
```bash
# Install web dependencies
make web-requirements

# Check port availability
lsof -i :8000

# Try different port
python web_app/app.py --port 8001
```

### File Permission Issues

**Problem**: "Permission denied" errors
```bash
# Fix file permissions
chmod +x scripts/*.sh
chmod +x scripts/*.py

# Run with sudo if needed (not recommended)
# Better to fix permissions properly

# Check user owns files
chown -R $USER:$USER rnaseq-mini/
```

## Network Issues

### Reference Download Problems

**Problem**: "Network timeout" when downloading references
```bash
# Use local references instead
reference:
  auto_download: false

# Download manually with wget
wget -O references/human/transcripts.fa.gz \
  https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
```

**Problem**: Proxy server issues
```bash
# Configure proxy in conda
conda config --set proxy_servers.http http://proxy.company.com:8080
conda config --set proxy_servers.https https://proxy.company.com:8080
```

## Engine-Specific Issues

### Snakemake Problems

**Problem**: "Missing rule" errors
```bash
# Update Snakemake
conda update snakemake

# Check syntax
snakemake --lint -s pipeline/snakemake/Snakefile

# Dry run to check dependencies
snakemake -n -s pipeline/snakemake/Snakefile
```

**Problem**: Cluster submission issues
```bash
# Test cluster configuration
snakemake --profile config/profiles/slurm.smk.yaml --dry-run

# Check cluster status
squeue -u $USER

# Fix cluster config file permissions
chmod 644 config/profiles/slurm.smk.yaml
```

### Nextflow Issues

**Problem**: "Process failed" errors
```bash
# Check Nextflow version
nextflow -version

# Update Nextflow
conda update nextflow

# Test with different profile
nextflow run pipeline/nextflow/main.nf -profile local -with-conda
```

## Getting Help

### Diagnostic Commands

```bash
# Full system diagnostics
python scripts/validate_config.py --comprehensive

# Resource estimation
make estimate

# Quality assessment
make assess-quality

# Cache statistics
make cache-stats
```

### Log Analysis

```bash
# Check main log files
tail -f logs/snakemake.log
tail -f logs/nextflow.log

# Check for error patterns
grep -i error logs/*.log

# Monitor resource usage
htop
```

### Community Support

- **[Discussion Forum](https://discourse.rnaseq-mini.org)** - Ask questions
- **[GitHub Issues](https://github.com/rnaseq-mini/rnaseq-mini/issues)** - Bug reports
- **[Email Support](support@rnaseq-mini.org)** - Direct assistance

### Professional Support

For enterprise users:
- **Priority Support**: Fast-tracked issue resolution
- **Custom Troubleshooting**: Personalized debugging sessions
- **Training**: Hands-on configuration and optimization

## Prevention

### Best Practices

1. **Validate before running**: `python scripts/validate_config.py`
2. **Start small**: Test with subset of data first
3. **Monitor resources**: Use `make monitor` during runs
4. **Enable caching**: Prevents redundant computations
5. **Keep environments updated**: Regular `conda update --all`

### Regular Maintenance

```bash
# Weekly maintenance
make cache-cleanup
conda clean --all

# Monthly maintenance
conda update --all -n rnaseq-mini-base
conda update --all -n rnaseq-mini-qc
conda update --all -n rnaseq-mini-salmon
conda update --all -n rnaseq-mini-r

# Before major analyses
make validate-full
```

## Emergency Recovery

### Complete Reset

```bash
# If everything else fails
make clean
rm -rf .cache .snakemake .nextflow*
conda env remove -n rnaseq-mini-base --all
conda env remove -n rnaseq-mini-qc --all
conda env remove -n rnaseq-mini-salmon --all
conda env remove -n rnaseq-mini-r --all

# Start fresh
make setup
python scripts/setup_wizard.py
```

## 📚 Related Documentation

- **[Quick Start](quickstart.md)** - Basic setup troubleshooting
- **[Standard Workflow](workflow.md)** - Analysis-specific issues
- **[Advanced Features](advanced.md)** - Enterprise troubleshooting
- **[Configuration](configuration.md)** - Parameter-specific problems












