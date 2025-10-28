#!/usr/bin/env bash
set -euo pipefail

# Script to compute deterministic hashes for key output files
# Usage: bash scripts/compute_hashes.sh <results-directory> <output-hash-file>

if [[ $# -ne 2 ]]; then
    echo "Usage: bash scripts/compute_hashes.sh <results-directory> <output-hash-file>"
    exit 1
fi

RESULTS_DIR=$1
HASH_FILE=$2

if [[ ! -d "$RESULTS_DIR" ]]; then
    echo "Error: Results directory '$RESULTS_DIR' not found."
    exit 1
fi

# Key files to hash (order matters for consistency)
FILES_TO_HASH=(
    "report.html"
    "qc/multiqc/multiqc_report.html"
    "counts/counts.tsv"
    "counts/tpm.tsv"
    "de/de_summary.tsv"
    "fgsea/fgsea_summary.tsv"
)

# Additional files for numerical validation
NUMERIC_FILES=(
    "counts/counts.tsv"
    "counts/tpm.tsv"
    "de/de_summary.tsv"
    "fgsea/fgsea_summary.tsv"
)

cd "$RESULTS_DIR" || { echo "Error: Cannot enter results directory '$RESULTS_DIR'"; exit 1; }

echo "[hash] Computing hashes for key outputs in '$RESULTS_DIR'..."

# Ensure hash file is empty before starting
rm -f "../$HASH_FILE"
touch "../$HASH_FILE"

# Compute SHA256 hashes
for file in "${FILES_TO_HASH[@]}"; do
    if [[ -f "$file" ]]; then
        # Using relative paths in the hash file for consistency
        sha256sum "$file" >> "../$HASH_FILE"
    else
        echo "Warning: $file not found in '$RESULTS_DIR', skipping"
    fi
done

echo "[hash] Hashes written to $HASH_FILE"

