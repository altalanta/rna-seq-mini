#!/usr/bin/env python3
"""
Enhanced config validation script for RNASEQ-MINI.
Validates config/params.yaml and related files against a comprehensive schema.
Provides interactive configuration wizard for common setups.
Usage: python scripts/validate_config.py [--wizard] [--fix] [config_file]
"""

import yaml
import sys
import json
import argparse
import pandas as pd
from schema import Schema, And, Or, Use, Optional, SchemaError
from pathlib import Path
import subprocess
import os

# Define schema for params.yaml
PARAMS_SCHEMA = Schema({
    'project': str,
    'engine': And(str, Or('snakemake', 'nextflow')),
    'se': bool,
    'threads': And(int, lambda x: x > 0),
    'memory_gb': And(int, lambda x: x > 0),
    'organism': str,
    'reference': {
        'transcripts_fa': And(str, lambda x: Path(x).suffix in ['.fa', '.fasta', '.fa.gz', '.fasta.gz']),
        'annotation_gtf': And(str, lambda x: Path(x).suffix in ['.gtf', '.gtf.gz']),
        'decoy_fasta': And(str, lambda x: Path(x).suffix in ['.fa', '.fasta', '.fa.gz', '.fasta.gz']),
        'salmon_index': And(str, lambda x: Path(x).is_dir() or not Path(x).exists()),  # Can be non-existent for auto-creation
        Optional('auto_download'): bool
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


def validate_file_extensions(config):
    """Validate that all file paths have correct extensions."""
    errors = []

    # Reference files
    ref_files = [
        (config['reference']['transcripts_fa'], 'transcripts_fa', ['.fa', '.fasta', '.fa.gz', '.fasta.gz']),
        (config['reference']['annotation_gtf'], 'annotation_gtf', ['.gtf', '.gtf.gz']),
        (config['reference']['decoy_fasta'], 'decoy_fasta', ['.fa', '.fasta', '.fa.gz', '.fasta.gz'])
    ]

    for file_path, field_name, valid_exts in ref_files:
        path = Path(file_path)
        if path.exists() and path.suffix not in valid_exts:
            errors.append(f"File {field_name} has invalid extension '{path.suffix}'. Expected: {valid_exts}")

    return errors


def validate_resource_requirements(config):
    """Validate resource requirements against system capabilities."""
    errors = []
    warnings = []

    # Check threads vs available CPUs
    available_cpus = os.cpu_count() or 4
    requested_threads = config['threads']
    salmon_threads = config['salmon']['threads']

    if requested_threads > available_cpus:
        errors.append(f"Requested {requested_threads} threads but system has only {available_cpus} CPUs")

    if salmon_threads > requested_threads:
        errors.append(f"Salmon threads ({salmon_threads}) cannot exceed main threads ({requested_threads})")

    # Check memory requirements
    requested_memory = config['memory_gb']
    try:
        # Get available memory (this is approximate)
        with open('/proc/meminfo', 'r') as f:
            for line in f:
                if line.startswith('MemTotal:'):
                    available_memory_gb = int(line.split()[1]) / (1024 * 1024)
                    break
            else:
                available_memory_gb = 8  # Default assumption

        if requested_memory > available_memory_gb * 0.9:  # Leave 10% buffer
            warnings.append(f"Requested {requested_memory}GB memory, system has ~{available_memory_gb:.1f}GB available")
    except:
        pass  # Skip memory check if we can't determine available memory

    return errors, warnings


def validate_sample_consistency(config):
    """Validate that samples and contrasts are consistent."""
    errors = []
    warnings = []

    samples_file = config['paths']['samples']
    contrasts_file = config['r']['contrasts_file']

    try:
        # Read samples
        samples_df = pd.read_csv(samples_file, sep='\t')
        required_cols = ['sample', 'condition']
        missing_cols = [col for col in required_cols if col not in samples_df.columns]
        if missing_cols:
            errors.append(f"Samples file missing required columns: {missing_cols}")
            return errors, warnings

        samples = set(samples_df['sample'].unique())
        conditions = set(samples_df['condition'].unique())

        # Read contrasts
        contrasts_df = pd.read_csv(contrasts_file, sep='\t')
        if 'groupA' not in contrasts_df.columns or 'groupB' not in contrasts_df.columns:
            errors.append("Contrasts file must have 'groupA' and 'groupB' columns")
            return errors, warnings

        contrast_conditions = set(contrasts_df['groupA']).union(set(contrasts_df['groupB']))

        # Check if all contrast conditions exist in samples
        missing_conditions = contrast_conditions - conditions
        if missing_conditions:
            errors.append(f"Contrasts reference conditions not found in samples: {missing_conditions}")

        # Check for single-sample conditions
        condition_counts = samples_df['condition'].value_counts()
        single_sample_conditions = condition_counts[condition_counts == 1].index.tolist()
        if single_sample_conditions:
            warnings.append(f"Conditions with only one sample: {single_sample_conditions}")

    except FileNotFoundError as e:
        errors.append(f"Could not read required file: {e.filename}")
    except Exception as e:
        errors.append(f"Error validating sample consistency: {e}")

    return errors, warnings


def run_comprehensive_validation(config_path):
    """Run all validation checks on a config file."""
    print("🔍 Running comprehensive configuration validation...")

    all_errors = []
    all_warnings = []

    # Basic schema validation
    config = None
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        PARAMS_SCHEMA.validate(config)
        print("✅ Schema validation passed")
    except SchemaError as e:
        print(f"❌ Schema validation failed: {e}")
        return False, ["Schema validation failed"]
    except FileNotFoundError:
        print(f"❌ Config file '{config_path}' not found")
        return False, [f"Config file '{config_path}' not found"]
    except Exception as e:
        print(f"❌ Error loading config: {e}")
        return False, [f"Error loading config: {e}"]

    # File extension validation
    file_errors = validate_file_extensions(config)
    all_errors.extend(file_errors)
    if file_errors:
        print(f"❌ Found {len(file_errors)} file extension errors")
        for error in file_errors:
            print(f"   {error}")

    # Resource validation
    resource_errors, resource_warnings = validate_resource_requirements(config)
    all_errors.extend(resource_errors)
    all_warnings.extend(resource_warnings)
    if resource_errors:
        print(f"❌ Found {len(resource_errors)} resource errors")
        for error in resource_errors:
            print(f"   {error}")
    if resource_warnings:
        print(f"⚠️  Found {len(resource_warnings)} resource warnings")
        for warning in resource_warnings:
            print(f"   {warning}")

    # Sample consistency validation
    sample_errors, sample_warnings = validate_sample_consistency(config)
    all_errors.extend(sample_errors)
    all_warnings.extend(sample_warnings)
    if sample_errors:
        print(f"❌ Found {len(sample_errors)} sample consistency errors")
        for error in sample_errors:
            print(f"   {error}")
    if sample_warnings:
        print(f"⚠️  Found {len(sample_warnings)} sample consistency warnings")
        for warning in sample_warnings:
            print(f"   {warning}")

    # Summary
    if not all_errors and not all_warnings:
        print("✅ All validations passed!")
        return True, []
    elif not all_errors:
        print(f"✅ Validation passed with {len(all_warnings)} warnings")
        return True, all_warnings
    else:
        print(f"❌ Validation failed with {len(all_errors)} errors")
        return False, all_errors + all_warnings


def interactive_config_wizard():
    """Interactive configuration wizard for common experimental designs."""
    print("🧙 Interactive RNA-seq Configuration Wizard")
    print("=" * 50)

    # Get project info
    project_name = input("Project name: ").strip() or "rnaseq-analysis"
    organism = input("Organism (yeast, human, mouse, zebrafish, fruitfly): ").strip().lower() or "yeast"

    # Sequencing parameters
    print("\nSequencing parameters:")
    se = input("Single-end reads? (y/N): ").strip().lower() == 'y'
    threads = int(input("Number of threads (default: 8): ").strip() or "8")
    memory = int(input("Memory in GB (default: 32): ").strip() or "32")

    # Experimental design
    print("\nExperimental design:")
    conditions_input = input("Condition names (comma-separated, e.g., 'control,treatment'): ").strip()
    conditions = [c.strip() for c in conditions_input.split(',') if c.strip()] or ['control', 'treatment']

    replicates = int(input("Replicates per condition (default: 2): ").strip() or "2")

    # Create sample sheet
    print("\n📋 Generating sample sheet...")
    samples_data = []
    for i, condition in enumerate(conditions, 1):
        for rep in range(1, replicates + 1):
            sample_name = f"{condition}_{rep}"
            samples_data.append({
                'sample': sample_name,
                'condition': condition,
                'fastq_1': f"data/{sample_name}_R1.fastq.gz",
                'fastq_2': f"data/{sample_name}_R2.fastq.gz" if not se else None
            })

    # Write sample sheet
    samples_file = Path("config/samples.tsv")
    samples_file.parent.mkdir(exist_ok=True)
    with open(samples_file, 'w') as f:
        if samples_data and samples_data[0].get('fastq_2'):
            f.write("sample\tcondition\tfastq_1\tfastq_2\n")
            for sample in samples_data:
                f.write(f"{sample['sample']}\t{sample['condition']}\t{sample['fastq_1']}\t{sample['fastq_2']}\n")
        else:
            f.write("sample\tcondition\tfastq_1\n")
            for sample in samples_data:
                f.write(f"{sample['sample']}\t{sample['condition']}\t{sample['fastq_1']}\n")

    # Create contrasts
    print("🔗 Generating contrasts...")
    contrasts_data = []
    for i, cond1 in enumerate(conditions):
        for cond2 in conditions[i+1:]:
            contrasts_data.append({'groupA': cond1, 'groupB': cond2})

    contrasts_file = Path("config/contrasts.tsv")
    with open(contrasts_file, 'w') as f:
        f.write("groupA\tgroupB\n")
        for contrast in contrasts_data:
            f.write(f"{contrast['groupA']}\t{contrast['groupB']}\n")

    # Generate config
    print("⚙️  Generating configuration...")
    config = {
        'project': project_name,
        'engine': 'snakemake',
        'se': se,
        'threads': threads,
        'memory_gb': memory,
        'organism': organism,
        'reference': {
            'transcripts_fa': f"references/{organism}/transcripts.fa.gz",
            'annotation_gtf': f"references/{organism}/annotation.gtf.gz",
            'decoy_fasta': f"references/{organism}/decoys.fa.gz",
            'salmon_index': f"references/{organism}/salmon_index",
            'auto_download': True
        },
        'fastqc': {'extra': ''},
        'multiqc': {'title': f"{project_name} - RNA-seq QC"},
        'salmon': {
            'libtype': 'A',
            'extra': '--validateMappings --gcBias',
            'threads': min(threads // 2, 4)  # Use half threads for salmon, max 4
        },
        'r': {
            'design': '~ condition',
            'contrast_variable': 'condition',
            'alpha': 0.05,
            'lfc_shrink': True,
            'contrasts_file': 'config/contrasts.tsv',
            'gene_id_column': 'gene_id',
            'pvalue_adjust': 'BH'
        },
        'fgsea': {
            'min_size': 15,
            'max_size': 500,
            'nperm': 1000,
            'padj_cutoff': 0.25,
            'score_column': 'log2FoldChange'
        },
        'report': {
            'author': 'Your Name',
            'title': f'{project_name} - RNA-seq Analysis Report',
            'output_html': 'results/report.html'
        },
        'paths': {
            'samples': 'config/samples.tsv',
            'outdir': 'results',
            'logs': 'logs',
            'qc': 'results/qc',
            'salmon': 'results/salmon',
            'counts': 'results/counts',
            'de': 'results/de',
            'fgsea': 'results/fgsea',
            'report_dir': 'results'
        },
        'containers': {
            'enabled': False,
            'image': 'ghcr.io/example/rnaseq-mini:latest'
        },
        'profiles': {
            'snakemake': 'config/profiles/local.smk.yaml',
            'nextflow': 'config/profiles/local.nf.config'
        }
    }

    # Write config
    config_file = Path("config/params.yaml")
    config_file.parent.mkdir(exist_ok=True)
    with open(config_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    print("\n✅ Configuration files generated!")
    print(f"   Config: {config_file}")
    print(f"   Samples: {samples_file}")
    print(f"   Contrasts: {contrasts_file}")
    print("\n🚀 You can now run:")
    print(f"   make run")
    print(f"   # or download references: make download-refs SPECIES={organism}")

    return True

def main():
    parser = argparse.ArgumentParser(description="Enhanced configuration validation for RNASEQ-MINI")
    parser.add_argument('--wizard', action='store_true', help='Run interactive configuration wizard')
    parser.add_argument('--comprehensive', action='store_true', help='Run comprehensive validation checks')
    parser.add_argument('config_file', nargs='?', default='config/params.yaml', help='Configuration file to validate')

    args = parser.parse_args()

    if args.wizard:
        # Run interactive wizard
        success = interactive_config_wizard()
        sys.exit(0 if success else 1)

    elif args.comprehensive:
        # Run comprehensive validation
        success, issues = run_comprehensive_validation(args.config_file)
        if not success:
            print(f"\n❌ Validation failed with {len(issues)} issues:")
            for issue in issues:
                print(f"   {issue}")
            sys.exit(1)
        else:
            print(f"\n✅ Validation successful{' with warnings' if issues else ''}!")
            if issues:
                print(f"⚠️  {len(issues)} warnings:")
                for issue in issues:
                    print(f"   {issue}")
            sys.exit(0)

    else:
        # Run basic validation (legacy mode)
        config_file = args.config_file
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

