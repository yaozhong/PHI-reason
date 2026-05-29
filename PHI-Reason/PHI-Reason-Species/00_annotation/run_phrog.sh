#!/usr/bin/env bash
# run_phrog.sh — PHROG annotation via DIAMOND blastp
# Annotates predicted proteins against the PHROG database.
# PHROG DB download: https://phrogs.lmge.uca.fr/
# Build DB first: diamond makedb --in phrogs_rep_seq.fasta -d phrogs_db
set -euo pipefail

PROT_DIR=""
PHROG_DB=""
OUT_DIR=""
THREADS=8

usage() {
    cat <<EOF
Usage: bash $0 --prot-dir <dir> --phrog-db <path> --out-dir <dir> [--threads 8]

  --prot-dir   Directory with protein FASTAs ({acc}/proteins.faa)
  --phrog-db   Path to DIAMOND-formatted PHROG database (.dmnd)
  --out-dir    Output directory for PHROG hit tables
  --threads    Number of threads (default: 8)
  -h, --help   Show this help message

Note: Download PHROG sequences from https://phrogs.lmge.uca.fr/
      Then build: diamond makedb --in phrogs_rep_seq.fasta -d phrogs_db
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --prot-dir)  PROT_DIR="$2"; shift 2 ;;
        --phrog-db)  PHROG_DB="$2"; shift 2 ;;
        --out-dir)   OUT_DIR="$2";  shift 2 ;;
        --threads)   THREADS="$2";  shift 2 ;;
        -h|--help)   usage 0 ;;
        *)           echo "Unknown option: $1"; usage 1 ;;
    esac
done

if [[ -z "$PROT_DIR" || -z "$PHROG_DB" || -z "$OUT_DIR" ]]; then
    echo "Error: --prot-dir, --phrog-db, and --out-dir are required."
    usage 1
fi

command -v diamond >/dev/null 2>&1 || { echo "Error: diamond not found in PATH"; exit 1; }

mkdir -p "$OUT_DIR"

mapfile -t FAA_FILES < <(find "$PROT_DIR" -name "proteins.faa" -type f | sort)
total=${#FAA_FILES[@]}
echo "Found $total protein file(s) under $PROT_DIR"

count=0
for faa in "${FAA_FILES[@]}"; do
    # Extract accession from parent directory name
    acc=$(basename "$(dirname "$faa")")
    mkdir -p "$OUT_DIR/$acc"

    diamond blastp \
        --db "$PHROG_DB" \
        --query "$faa" \
        --out "$OUT_DIR/$acc/phrog_hits.tsv" \
        --outfmt 6 qseqid sseqid pident qcovhsp bitscore \
        --threads "$THREADS" \
        --quiet

    count=$((count + 1))
    if (( count % 50 == 0 )); then
        echo "  Progress: $count / $total files processed"
    fi
done

echo "Done. Processed $count file(s). Output in $OUT_DIR"
