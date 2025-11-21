#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage: $0 -t transcripts.fa[.gz] -a annotation.gtf[.gz] -o index_dir [-d decoy.fa[.gz]] [-p threads]
Build a decoy-aware Salmon index.
USAGE
}

THREADS=4
TRANSCRIPTS=""
ANNOTATION=""
DECOY=""
OUTDIR=""

while getopts ":t:a:d:o:p:h" opt; do
  case ${opt} in
    t) TRANSCRIPTS=${OPTARG} ;;
    a) ANNOTATION=${OPTARG} ;;
    d) DECOY=${OPTARG} ;;
    o) OUTDIR=${OPTARG} ;;
    p) THREADS=${OPTARG} ;;
    h) usage; exit 0 ;;
    :) echo "Missing argument for -${OPTARG}" >&2; usage; exit 1 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "${TRANSCRIPTS}" || -z "${ANNOTATION}" || -z "${OUTDIR}" ]]; then
  echo "Error: transcripts (-t), annotation (-a), and output (-o) are required." >&2
  usage
  exit 1
fi

mkdir -p "${OUTDIR}"
TMPDIR=$(mktemp -d)
trap 'rm -rf "${TMPDIR}"' EXIT

extract() {
  local input=$1
  local output=$2
  case "${input}" in
    *.gz) gzip -cd "${input}" > "${output}" ;;
    *) cp "${input}" "${output}" ;;
  esac
}

TRANSCRIPTS_FA="${TMPDIR}/transcripts.fa"
ANNOTATION_GTF="${TMPDIR}/annotation.gtf"
extract "${TRANSCRIPTS}" "${TRANSCRIPTS_FA}"
extract "${ANNOTATION}" "${ANNOTATION_GTF}"

GENTROME="${TMPDIR}/gentrome.fa"
DECOYS_TXT="${TMPDIR}/decoys.txt"
cp "${TRANSCRIPTS_FA}" "${GENTROME}"

touch "${DECOYS_TXT}"
if [[ -n "${DECOY}" ]]; then
  DECOY_FA="${TMPDIR}/decoys.fa"
  extract "${DECOY}" "${DECOY_FA}"
  awk '/^>/{print substr($0,2)}' "${DECOY_FA}" > "${DECOYS_TXT}"
  cat "${DECOY_FA}" >> "${GENTROME}"
else
  echo "Warning: no decoy FASTA provided; building transcript-only index." >&2
fi

awk 'BEGIN{FS="\t";OFS="\t"}
     $3=="transcript" {
       match($0, /transcript_id "([^"]+)"/, tx)
       match($0, /gene_id "([^"]+)"/, gene)
       if (tx[1] != "" && gene[1] != "") {
         print tx[1], gene[1]
       }
     }' "${ANNOTATION_GTF}" | sort -u > "${OUTDIR}/tx2gene.tsv"

if [[ -s "${DECOYS_TXT}" ]]; then
  salmon index --transcripts "${GENTROME}" --decoys "${DECOYS_TXT}" --index "${OUTDIR}" --threads "${THREADS}" --kmerLen 31
else
  salmon index --transcripts "${GENTROME}" --index "${OUTDIR}" --threads "${THREADS}" --kmerLen 31
fi

echo "Salmon index built at ${OUTDIR}" >&2
