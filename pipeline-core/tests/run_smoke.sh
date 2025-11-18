#!/usr/bin/env bash
set -euo pipefail

# This script should be run from the pipeline-core directory
ROOT_DIR=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT_DIR"

RESULTS_DIR=${1:-results}
echo "[smoke] Using results directory: $RESULTS_DIR"

# Ensure results directory exists and is empty
rm -rf "$RESULTS_DIR"
mkdir -p "$RESULTS_DIR"

TRANSCRIPTS="references/yeast/transcripts.fa.gz"
ANNOTATION="references/yeast/annotation.gtf.gz"
DECOY="references/yeast/decoys.fa.gz"
INDEX_DIR="references/yeast/salmon_index"

if [[ ! -d "$INDEX_DIR" ]]; then
  echo "[smoke] Building Salmon index..."
  mkdir -p "$INDEX_DIR"
  scripts/build_salmon_index.sh \
    -t "$TRANSCRIPTS" \
    -a "$ANNOTATION" \
    -d "$DECOY" \
    -o "$INDEX_DIR" \
    -p 2
fi

echo "[smoke] Running Snakemake..."
snakemake \
  -s workflow/Snakefile \
  --configfile config/params.yaml \
  --use-conda \
  --cores 2 \
  --config results_dir="$RESULTS_DIR"

[[ -f "$RESULTS_DIR/report.html" ]] || { echo "Snakemake report missing"; exit 1; }
[[ -f "$RESULTS_DIR/counts/counts.tsv" ]] || { echo "Snakemake counts missing"; exit 1; }
ls "$RESULTS_DIR/de/DE_"*vs_*.tsv >/dev/null 2>&1 || { echo "Snakemake DE results missing"; exit 1; }

echo "[smoke] Smoke test finished successfully."
