#!/usr/bin/env bash
set -euo pipefail

# Script to auto-download human or mouse references if missing and auto_download is enabled
# Usage: bash scripts/download_references.sh [organism]

ORGANISM=${1:-$(python -c 'import yaml; print(yaml.safe_load(open("config/params.yaml"))["organism"])')}
AUTO_DOWNLOAD=$(python -c 'import yaml; print(yaml.safe_load(open("config/params.yaml"))["reference"]["auto_download"])')

if [[ "$AUTO_DOWNLOAD" != "true" ]]; then
    echo "[download] Auto-download disabled in config. Skipping."
    exit 0
fi

REF_DIR="references/$ORGANISM"
mkdir -p "$REF_DIR"

# Define download URLs (example for human and mouse; add more as needed)
declare -A TRANSCRIPTS=(
    ["human"]="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.transcripts.fa.gz"
    ["mouse"]="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.transcripts.fa.gz"
)

declare -A ANNOTATION=(
    ["human"]="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz"
    ["mouse"]="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz"
)

declare -A DECOY=(
    ["human"]="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.transcripts.fa.gz"
    ["mouse"]="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.transcripts.fa.gz"
)

if [[ -z "${TRANSCRIPTS[$ORGANISM]}" ]]; then
    echo "[download] Unsupported organism: $ORGANISM. Supported: human, mouse."
    exit 1
fi

echo "[download] Downloading references for $ORGANISM..."

# Download transcripts
TRANSCRIPTS_FILE="$REF_DIR/transcripts.fa.gz"
if [[ ! -f "$TRANSCRIPTS_FILE" ]]; then
    wget -O "$TRANSCRIPTS_FILE" "${TRANSCRIPTS[$ORGANISM]}"
    echo "[download] Downloaded transcripts to $TRANSCRIPTS_FILE"
fi

# Download annotation
ANNOTATION_FILE="$REF_DIR/annotation.gtf.gz"
if [[ ! -f "$ANNOTATION_FILE" ]]; then
    wget -O "$ANNOTATION_FILE" "${ANNOTATION[$ORGANISM]}"
    echo "[download] Downloaded annotation to $ANNOTATION_FILE"
fi

# Download decoy (use transcripts as decoy for simplicity)
DECOY_FILE="$REF_DIR/decoys.fa.gz"
if [[ ! -f "$DECOY_FILE" ]]; then
    cp "$TRANSCRIPTS_FILE" "$DECOY_FILE"
    echo "[download] Created decoy file from transcripts: $DECOY_FILE"
fi

echo "[download] References ready in $REF_DIR"
