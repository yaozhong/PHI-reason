#!/usr/bin/env bash
# run_eggnog.sh — Functional annotation with eggNOG-mapper
# Requires: eggNOG-mapper installed (pip install eggnog-mapper)
#           eggNOG databases downloaded (download_eggnog_data.py)
set -euo pipefail

PROT_DIR=""
EGGNOG_DATA=""
OUT_DIR=""
THREADS=8

usage() {
    cat <<EOF
Usage: bash $0 --prot-dir <dir> --eggnog-data <dir> --out-dir <dir> [--threads 8]

  --prot-dir     Directory with protein FASTAs ({acc}/proteins.faa)
  --eggnog-data  Path to eggNOG database directory
  --out-dir      Output directory for eggNOG annotations
  --threads      Number of CPU threads (default: 8)
  -h, --help     Show this help message

Note: Requires eggNOG-mapper (pip install eggnog-mapper) and databases
      downloaded via: download_eggnog_data.py -y --data_dir <dir>
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --prot-dir)    PROT_DIR="$2";    shift 2 ;;
        --eggnog-data) EGGNOG_DATA="$2"; shift 2 ;;
        --out-dir)     OUT_DIR="$2";     shift 2 ;;
        --threads)     THREADS="$2";     shift 2 ;;
        -h|--help)     usage 0 ;;
        *)             echo "Unknown option: $1"; usage 1 ;;
    esac
done

if [[ -z "$PROT_DIR" || -z "$EGGNOG_DATA" || -z "$OUT_DIR" ]]; then
    echo "Error: --prot-dir, --eggnog-data, and --out-dir are required."
    usage 1
fi

command -v emapper.py >/dev/null 2>&1 || { echo "Error: emapper.py not found in PATH"; exit 1; }

mkdir -p "$OUT_DIR"

mapfile -t FAA_FILES < <(find "$PROT_DIR" -name "proteins.faa" -type f | sort)
total=${#FAA_FILES[@]}
echo "Found $total protein file(s) under $PROT_DIR"

count=0
for faa in "${FAA_FILES[@]}"; do
    acc=$(basename "$(dirname "$faa")")
    mkdir -p "$OUT_DIR/$acc"

    emapper.py \
        -i "$faa" \
        -o "$OUT_DIR/$acc/eggnog" \
        --data_dir "$EGGNOG_DATA" \
        --cpu "$THREADS" \
        --override

    count=$((count + 1))
    if (( count % 50 == 0 )); then
        echo "  Progress: $count / $total files processed"
    fi
done

echo "Done. Processed $count file(s). Output in $OUT_DIR"
