# ⚙️ Configuration Reference

Complete parameter reference for customizing RNASEQ-MINI analyses.

## Project Configuration

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `project` | Project name for reports | `"my-rnaseq-analysis"` | `"cancer-study-2024"` |
| `engine` | Workflow engine | `"snakemake"` | `"nextflow"` |
| `threads` | Number of threads | `8` | `16` |
| `memory_gb` | Memory allocation | `32` | `64` |
| `organism` | Reference organism | `"human"` | `"mouse"`, `"yeast"` |

## Reference Configuration

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `transcripts_fa` | Transcript FASTA file | `"references/human/transcripts.fa.gz"` | `"data/custom_transcripts.fa.gz"` |
| `annotation_gtf` | Gene annotation GTF | `"references/human/annotation.gtf.gz"` | `"data/custom_annotation.gtf.gz"` |
| `decoy_fasta` | Decoy sequences for Salmon | `"references/human/decoys.fa.gz"` | `"data/custom_decoys.fa.gz"` |
| `salmon_index` | Pre-built Salmon index | `"references/human/salmon_index"` | `"data/custom_index"` |
| `auto_download` | Auto-download missing references | `true` | `false` |

## Quality Control Parameters

### FastQC Configuration
```yaml
fastqc:
  extra: ""  # Additional FastQC parameters
```

### MultiQC Configuration
```yaml
multiqc:
  title: "RNA-seq QC"  # Report title
```

## Quantification Parameters

### Salmon Configuration
```yaml
salmon:
  libtype: "A"  # Library type (A, ISR, ISF, SR, SF)
  extra: "--validateMappings --gcBias"  # Additional parameters
  threads: 4  # Threads for Salmon
```

**Library Type Options:**
- `A`: Auto-detect (recommended for most cases)
- `ISR`: Inward, Stranded, Reverse
- `ISF`: Inward, Stranded, Forward
- `SR`: Stranded, Reverse
- `SF`: Stranded, Forward

## Statistical Analysis Parameters

### DESeq2 Configuration
```yaml
r:
  design: "~ condition"  # Statistical design formula
  contrast_variable: "condition"  # Variable for contrasts
  alpha: 0.05  # Significance threshold
  lfc_shrink: true  # Shrink log fold changes
  gene_id_column: "gene_id"  # Gene ID column name
  pvalue_adjust: "BH"  # Multiple testing correction
```

**Design Formula Examples:**
- `"~ condition"` - Simple two-group comparison
- `"~ condition + batch"` - With batch correction
- `"~ condition + batch + covariate"` - Multiple factors

### Contrast Configuration
```yaml
# config/contrasts.tsv
groupA	groupB
treated	control
disease	healthy
timepoint1	timepoint2
```

## Pathway Analysis Parameters

### fgsea Configuration
```yaml
fgsea:
  genesets: null  # Path to gene set file (GMT/TSV)
  min_size: 15  # Minimum genes per pathway
  max_size: 500  # Maximum genes per pathway
  nperm: 1000  # Permutations for significance
  padj_cutoff: 0.05  # Adjusted p-value threshold
  score_column: "log2FoldChange"  # Column for ranking
```

**Gene Set Sources:**
- Built-in: `"hallmark"`, `"kegg"`, `"reactome"`, `"go_bp"`
- Custom: Path to GMT or TSV file

## Caching Configuration

```yaml
cache:
  enabled: true  # Enable intelligent caching
  dir: ".cache"  # Cache directory
  max_age_days: 30  # Auto-clean old entries
```

## Output Paths

```yaml
paths:
  samples: "config/samples.tsv"  # Sample metadata file
  outdir: "results"  # Main output directory
  logs: "logs"  # Log files directory
  qc: "results/qc"  # Quality control outputs
  salmon: "results/salmon"  # Salmon quantification
  counts: "results/counts"  # Gene counts
  de: "results/de"  # Differential expression
  fgsea: "results/fgsea"  # Pathway analysis
  report_dir: "results"  # Report output directory
```

## Container Configuration

```yaml
containers:
  enabled: false  # Use containerized execution
  image: "ghcr.io/example/rnaseq-mini:latest"  # Container image
```

## Profile Configuration

```yaml
profiles:
  snakemake: "config/profiles/local.smk.yaml"  # Snakemake profile
  nextflow: "config/profiles/local.nf.config"  # Nextflow profile
```

## Advanced Parameters

### Single-End Reads
```yaml
se: false  # Set to true for single-end data
```

### Custom Scripts
```yaml
# Add custom processing steps
custom_scripts:
  pre_processing: "scripts/custom_preprocessing.R"
  post_processing: "scripts/custom_postprocessing.py"
```

### Parallel Processing
```yaml
parallel:
  max_jobs: 100  # Maximum concurrent jobs
  job_script: "scripts/parallel_job.sh"  # Custom job script
```

## Environment Variables

| Variable | Description | Example |
|----------|-------------|---------|
| `RNASEQ_MINI_THREADS` | Override thread count | `export RNASEQ_MINI_THREADS=16` |
| `RNASEQ_MINI_MEMORY` | Override memory allocation | `export RNASEQ_MINI_MEMORY=64` |
| `RNASEQ_MINI_CACHE_DIR` | Custom cache directory | `export RNASEQ_MINI_CACHE_DIR=/tmp/cache` |

## Validation & Diagnostics

### Configuration Validation
```bash
# Validate configuration syntax
python scripts/validate_config.py

# Interactive configuration wizard
python scripts/setup_wizard.py

# Comprehensive validation
python scripts/validate_config.py --comprehensive
```

### Resource Estimation
```bash
# Estimate optimal resource allocation
make estimate

# Generate optimized configuration
make optimize
```

## Troubleshooting Configuration

### Common Issues

**Memory Errors:**
```yaml
# Reduce memory allocation
memory_gb: 16  # Down from 32
threads: 4     # Down from 8
```

**Slow Performance:**
```yaml
# Optimize for speed
salmon:
  threads: 8  # Increase Salmon threads

cache:
  enabled: true  # Enable caching
```

**Disk Space Issues:**
```yaml
# Use faster storage
cache:
  dir: "/tmp/cache"  # Use RAM disk if available

# Reduce cache retention
cache:
  max_age_days: 7  # Shorter retention
```

## Examples

### Minimal Configuration
```yaml
project: "quick-analysis"
organism: "human"
threads: 4
r:
  design: "~ condition"
```

### Production Configuration
```yaml
project: "production-study"
engine: "snakemake"
threads: 16
memory_gb: 64
organism: "human"
cache:
  enabled: true
  max_age_days: 30
r:
  design: "~ condition + batch"
  alpha: 0.01
salmon:
  libtype: "ISR"
  extra: "--validateMappings --gcBias"
```

### Multi-Omics Configuration
```yaml
project: "multiomics-study"
multiomics:
  enabled: true
  data_types: ["rnaseq", "atacseq", "proteomics"]
  normalization: "combat"
cache:
  enabled: true
quality:
  assessment_enabled: true
  benchmarking_enabled: true
```

## Best Practices

1. **Start Simple**: Begin with default parameters, then optimize
2. **Version Control**: Track configuration changes with git
3. **Documentation**: Comment complex parameter choices
4. **Validation**: Always validate configuration before running
5. **Backups**: Keep backups of working configurations

## Support

For configuration help:
- **[Interactive Wizard](scripts/setup_wizard.py)** - Guided setup
- **[Validation Tool](scripts/validate_config.py)** - Configuration checking
- **[Troubleshooting Guide](troubleshooting.md)** - Common solutions










