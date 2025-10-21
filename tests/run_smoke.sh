#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT_DIR"

ENGINE=${PIPELINE_ENGINE:-both}

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
    --cores 2

  [[ -f "results/report.html" ]] || { echo "Snakemake report missing"; exit 1; }
  [[ -f "results/counts/counts.tsv" ]] || { echo "Snakemake counts missing"; exit 1; }
  ls results/de/DE_*_vs_*.tsv >/dev/null 2>&1 || { echo "Snakemake DE results missing"; exit 1; }

  rm -rf work .nextflow* results
fi

if [[ "$ENGINE" == "both" || "$ENGINE" == "nextflow" ]]; then
  rm -rf results work .nextflow* || true
  echo "[smoke] Running Nextflow..."
  nextflow run pipeline/nextflow/main.nf \
    -params-file config/params.yaml \
    -with-conda \
    -profile local

  [[ -f "results/report.html" ]] || { echo "Nextflow report missing"; exit 1; }
  [[ -f "results/counts/counts.tsv" ]] || { echo "Nextflow counts missing"; exit 1; }
  ls results/de/DE_*_vs_*.tsv >/dev/null 2>&1 || { echo "Nextflow DE results missing"; exit 1; }
fi

echo "[smoke] Smoke test finished successfully."
