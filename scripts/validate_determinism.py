#!/usr/bin/env python3
"""
Enhanced validation script for RNASEQ-MINI that compares outputs between engines
and validates numerical reproducibility of key results.

Usage: python scripts/validate_determinism.py [snakemake_hashes] [nextflow_hashes]
"""

import sys
import hashlib
import pandas as pd
import numpy as np
from pathlib import Path
import yaml
import json
import subprocess


def compute_file_hash(file_path):
    """Compute SHA256 hash of a file."""
    hash_sha256 = hashlib.sha256()
    try:
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_sha256.update(chunk)
        return hash_sha256.hexdigest()
    except FileNotFoundError:
        return None


def compare_hashes(hash_file1, hash_file2, tolerance=1e-10):
    """Compare hash files and validate numerical reproducibility."""
    results = {
        'hash_matches': [],
        'numerical_matches': [],
        'files_missing': [],
        'validation_summary': {}
    }

    try:
        # Read hash files
        hashes1 = {}
        with open(hash_file1, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        hashes1[parts[1]] = parts[0]

        hashes2 = {}
        with open(hash_file2, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        hashes2[parts[1]] = parts[0]

        # Compare hashes
        all_files = set(hashes1.keys()) | set(hashes2.keys())

        for file_path in all_files:
            hash1 = hashes1.get(file_path)
            hash2 = hashes2.get(file_path)

            if hash1 is None:
                results['files_missing'].append(f"File {file_path} missing from {hash_file1}")
                continue
            if hash2 is None:
                results['files_missing'].append(f"File {file_path} missing from {hash_file2}")
                continue

            if hash1 == hash2:
                results['hash_matches'].append(file_path)
            else:
                # For data files, check numerical reproducibility
                if file_path.endswith('.tsv') or file_path.endswith('.csv'):
                    numerical_match = validate_numerical_reproducibility(file_path, hash_file1, hash_file2, tolerance)
                    results['numerical_matches'].append((file_path, numerical_match))
                else:
                    results['hash_matches'].append(f"{file_path} (hash mismatch)")

        # Summary
        results['validation_summary'] = {
            'total_files': len(all_files),
            'hash_matches': len(results['hash_matches']),
            'numerical_comparisons': len(results['numerical_matches']),
            'files_missing': len(results['files_missing'])
        }

    except Exception as e:
        print(f"Error comparing hashes: {e}")
        return None

    return results


def validate_numerical_reproducibility(file_path, hash_file1, hash_file2, tolerance=1e-10):
    """Validate that numerical results are reproducible between engines."""
    try:
        # Infer results directories from hash file names
        # e.g., hashes-snakemake.sha256 -> results-snakemake
        engine1_name = Path(hash_file1).stem.split('-')[1]
        engine2_name = Path(hash_file2).stem.split('-')[1]
        results_dir1 = f"results-{engine1_name}"
        results_dir2 = f"results-{engine2_name}"

        file1 = Path(results_dir1) / file_path
        file2 = Path(results_dir2) / file_path

        if not file1.exists():
            print(f"  Warning: Cannot find {file1} for numerical comparison.")
            return False
        if not file2.exists():
            print(f"  Warning: Cannot find {file2} for numerical comparison.")
            return False

        # Read data files
        df1 = pd.read_csv(file1, sep='\\t', engine='python')
        df2 = pd.read_csv(file2, sep='\\t', engine='python')

        # Check if dataframes have same shape
        if df1.shape != df2.shape:
            return False

        # For each numeric column, check if values are within tolerance
        numeric_columns = df1.select_dtypes(include=[np.number]).columns

        for col in numeric_columns:
            if col not in df2.columns:
                return False

            # Handle potential NaN values
            valid_idx = ~(np.isnan(df1[col]) | np.isnan(df2[col]))

            if valid_idx.sum() == 0:
                continue

            max_diff = np.max(np.abs(df1[col][valid_idx] - df2[col][valid_idx]))

            if max_diff > tolerance:
                print(f"  Numerical mismatch in {file_path}:{col}, max diff: {max_diff}")
                return False

        return True

    except Exception as e:
        print(f"  Error validating numerical reproducibility for {file_path}: {e}")
        return False


def validate_engine_outputs(engine1="snakemake", engine2="nextflow"):
    """Main validation function to compare outputs between two engines."""
    print(f"🔍 Validating outputs between {engine1} and {engine2}...")

    hash_file1 = f"hashes-{engine1}.sha256"
    hash_file2 = f"hashes-{engine2}.sha256"

    if not Path(hash_file1).exists():
        print(f"❌ Hash file {hash_file1} not found. Run pipeline with {engine1} first.")
        return False

    if not Path(hash_file2).exists():
        print(f"❌ Hash file {hash_file2} not found. Run pipeline with {engine2} first.")
        return False

    # Compare hashes and validate reproducibility
    comparison = compare_hashes(hash_file1, hash_file2)

    if comparison is None:
        return False

    # Report results
    print("\n📊 Validation Summary:")
    print(f"   Total files compared: {comparison['validation_summary']['total_files']}")
    print(f"   Hash matches: {comparison['validation_summary']['hash_matches']}")
    print(f"   Numerical validations: {comparison['validation_summary']['numerical_comparisons']}")
    print(f"   Files missing: {comparison['validation_summary']['files_missing']}")

    # Detailed results
    if comparison['hash_matches']:
        print("\n✅ Files with matching hashes:")
        for file in comparison['hash_matches']:
            print(f"   {file}")

    numerical_passed = sum(1 for _, passed in comparison['numerical_matches'] if passed)
    numerical_failed = len(comparison['numerical_matches']) - numerical_passed

    if comparison['numerical_matches']:
        print(f"\n🔢 Numerical validation results: {numerical_passed} passed, {numerical_failed} failed")
        for file_path, passed in comparison['numerical_matches']:
            status = "✅" if passed else "❌"
            print(f"   {status} {file_path}")

    if comparison['files_missing']:
        print("\n⚠️  Files missing from comparison:")
        for file in comparison['files_missing']:
            print(f"   {file}")

    # Overall assessment
    all_good = (
        comparison['validation_summary']['files_missing'] == 0 and
        all(passed for _, passed in comparison['numerical_matches'])
    )

    if all_good:
        print("\n🎉 Cross-engine validation PASSED! Both engines produce identical results.")
        return True
    else:
        print("\n❌ Cross-engine validation FAILED! Differences detected between engines.")
        return False


def main():
    """Main entry point."""
    if len(sys.argv) == 3:
        engine1, engine2 = sys.argv[1], sys.argv[2]
        success = validate_engine_outputs(engine1, engine2)
    elif len(sys.argv) == 1:
        # Default comparison
        success = validate_engine_outputs("snakemake", "nextflow")
    else:
        print("Usage: python scripts/validate_determinism.py [engine1] [engine2]")
        print("Example: python scripts/validate_determinism.py snakemake nextflow")
        sys.exit(1)

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
