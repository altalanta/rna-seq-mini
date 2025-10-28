#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT_DIR"

ENGINE=${PIPELINE_ENGINE:-both}
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

if [[ "$ENGINE" == "both" || "$ENGINE" == "snakemake" ]]; then
  echo "[smoke] Running Snakemake..."
  snakemake \
    -s pipeline/snakemake/Snakefile \
    --configfile config/params.yaml \
    --use-conda \
    --cores 2 \
    --config results_dir="$RESULTS_DIR"

  [[ -f "$RESULTS_DIR/report.html" ]] || { echo "Snakemake report missing"; exit 1; }
  [[ -f "$RESULTS_DIR/counts/counts.tsv" ]] || { echo "Snakemake counts missing"; exit 1; }
  ls "$RESULTS_DIR/de/DE_"*vs_*.tsv >/dev/null 2>&1 || { echo "Snakemake DE results missing"; exit 1; }

fi

if [[ "$ENGINE" == "both" || "$ENGINE" == "nextflow" ]]; then
  echo "[smoke] Running Nextflow..."
  nextflow run pipeline/nextflow/main.nf \
    -params-file config/params.yaml \
    --results_dir "$RESULTS_DIR" \
    -with-conda \
    -profile local

  [[ -f "$RESULTS_DIR/report.html" ]] || { echo "Nextflow report missing"; exit 1; }
  [[ -f "$RESULTS_DIR/counts/counts.tsv" ]] || { echo "Nextflow counts missing"; exit 1; }
  ls "$RESULTS_DIR/de/DE_"*vs_*.tsv >/dev/null 2>&1 || { echo "Nextflow DE results missing"; exit 1; }
fi

echo "[smoke] Smoke test finished successfully."
