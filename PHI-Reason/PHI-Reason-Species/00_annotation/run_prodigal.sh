#!/usr/bin/env bash
# run_prodigal.sh — Gene prediction with Prodigal
# Runs Prodigal on each genome FASTA in a directory.
set -euo pipefail

GENOME_DIR=""
OUT_DIR=""
MODE="single"

usage() {
    cat <<EOF
Usage: bash $0 --genome-dir <dir> --out-dir <dir> [--mode single|meta]

  --genome-dir  Directory containing .fna/.fa/.fasta genome files
  --out-dir     Output directory (proteins.faa per genome)
  --mode        Prodigal mode: single (default) or meta
  -h, --help    Show this help message
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome-dir) GENOME_DIR="$2"; shift 2 ;;
        --out-dir)    OUT_DIR="$2";    shift 2 ;;
        --mode)       MODE="$2";       shift 2 ;;
        -h|--help)    usage 0 ;;
        *)            echo "Unknown option: $1"; usage 1 ;;
    esac
done

if [[ -z "$GENOME_DIR" || -z "$OUT_DIR" ]]; then
    echo "Error: --genome-dir and --out-dir are required."
    usage 1
fi

if [[ "$MODE" != "single" && "$MODE" != "meta" ]]; then
    echo "Error: --mode must be 'single' or 'meta'."
    exit 1
fi

command -v prodigal >/dev/null 2>&1 || { echo "Error: prodigal not found in PATH"; exit 1; }

mkdir -p "$OUT_DIR"

count=0
total=$(find "$GENOME_DIR" -maxdepth 1 -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) | wc -l)
echo "Found $total genome(s) in $GENOME_DIR"

for genome in "$GENOME_DIR"/*.fna "$GENOME_DIR"/*.fa "$GENOME_DIR"/*.fasta; do
    [[ -f "$genome" ]] || continue
    acc=$(basename "$genome" | sed 's/\.\(fna\|fa\|fasta\)$//')
    mkdir -p "$OUT_DIR/$acc"

    prodigal -i "$genome" -a "$OUT_DIR/$acc/proteins.faa" -p "$MODE" -o /dev/null 2>/dev/null

    count=$((count + 1))
    if (( count % 50 == 0 )); then
        echo "  Progress: $count / $total genomes processed"
    fi
done

echo "Done. Processed $count genome(s). Output in $OUT_DIR"
