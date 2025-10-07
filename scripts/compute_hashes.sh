#!/usr/bin/env bash
set -euo pipefail

# Script to compute deterministic hashes for key output files
# Usage: bash scripts/compute_hashes.sh [output_file]

OUTPUT_DIR="results"
HASH_FILE="hashes.sha256"

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

cd "$OUTPUT_DIR" || { echo "Error: Cannot find results directory"; exit 1; }

echo "[hash] Computing hashes for key outputs..."

# Compute SHA256 hashes
for file in "${FILES_TO_HASH[@]}"; do
    if [[ -f "$file" ]]; then
        sha256sum "$file" >> "../$HASH_FILE"
    else
        echo "Warning: $file not found, skipping"
    fi
done

echo "[hash] Hashes written to $HASH_FILE"

# If an output file is specified, save hashes there
if [[ $# -eq 1 ]]; then
    cp "../$HASH_FILE" "$1"
    echo "[hash] Hashes also saved to $1"
fi

