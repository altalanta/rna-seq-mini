#!/usr/bin/env python3
"""
Config validation script for RNASEQ-MINI.
Validates config/params.yaml and related files against a schema.
Usage: python scripts/validate_config.py
"""

import yaml
import sys
from schema import Schema, And, Or, Use, Optional
from pathlib import Path

# Define schema for params.yaml
PARAMS_SCHEMA = Schema({
    'project': str,
    'engine': And(str, Or('snakemake', 'nextflow')),
    'se': bool,
    'threads': And(int, lambda x: x > 0),
    'memory_gb': And(int, lambda x: x > 0),
    'organism': str,
    'reference': {
        'transcripts_fa': str,
        'annotation_gtf': str,
        'decoy_fasta': str,
        'salmon_index': str
    },
    'fastqc': {
        Optional('extra'): str
    },
    'multiqc': {
        'title': str
    },
    'salmon': {
        'libtype': str,
        Optional('extra'): str,
        'threads': And(int, lambda x: x > 0)
    },
    'r': {
        'design': str,
        'contrast_variable': str,
        Optional('batch_column'): str,
        'alpha': And(float, lambda x: 0 < x < 1),
        'lfc_shrink': bool,
        'contrasts_file': str,
        'gene_id_column': str,
        'pvalue_adjust': str
    },
    'fgsea': {
        Optional('genesets'): Or(str, None),
        'min_size': And(int, lambda x: x > 0),
        'max_size': And(int, lambda x: x > 0),
        'nperm': And(int, lambda x: x > 0),
        'padj_cutoff': And(float, lambda x: 0 < x <= 1),
        'score_column': str
    },
    'report': {
        'author': str,
        'title': str,
        'output_html': str
    },
    'paths': {
        'samples': str,
        'outdir': str,
        'logs': str,
        'qc': str,
        'salmon': str,
        'counts': str,
        'de': str,
        'fgsea': str,
        'report_dir': str
    },
    'containers': {
        'enabled': bool,
        Optional('image'): str
    },
    'profiles': {
        'snakemake': str,
        'nextflow': str
    }
})

def validate_config(config_path, schema):
    """Validate a config file against a schema."""
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: Config file '{config_path}' not found.")
        return False
    except yaml.YAMLError as e:
        print(f"Error parsing YAML in '{config_path}': {e}")
        return False

    try:
        schema.validate(config)
        print(f"✓ '{config_path}' is valid.")
        return True
    except Exception as e:
        print(f"✗ Validation failed for '{config_path}': {e}")
        return False

def check_file_exists(file_path, description):
    """Check if a referenced file exists."""
    path = Path(file_path)
    if not path.exists():
        print(f"Warning: {description} '{file_path}' does not exist.")
        return False
    return True

def main():
    config_file = "config/params.yaml"
    contrasts_file = "config/contrasts.tsv"
    samples_file = "config/samples.tsv"

    all_valid = True

    # Validate main config
    all_valid &= validate_config(config_file, PARAMS_SCHEMA)

    # Check for required files
    if all_valid:
        config = yaml.safe_load(open(config_file, 'r'))
        all_valid &= check_file_exists(config['r']['contrasts_file'], "Contrasts file")
        all_valid &= check_file_exists(config['paths']['samples'], "Samples file")

    if all_valid:
        print("✓ All configurations are valid and files exist.")
        sys.exit(0)
    else:
        print("✗ Configuration validation failed.")
        sys.exit(1)

if __name__ == "__main__":
    main()
