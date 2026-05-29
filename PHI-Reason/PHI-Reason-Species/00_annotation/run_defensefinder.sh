#!/usr/bin/env bash
# run_defensefinder.sh — Defense system detection with DefenseFinder
# Requires: DefenseFinder installed (pip install mdmparis-defense-finder)
#           Models updated via: defense-finder update
set -euo pipefail

PROT_DIR=""
OUT_DIR=""

usage() {
    cat <<EOF
Usage: bash $0 --prot-dir <dir> --out-dir <dir>

  --prot-dir   Directory with host protein FASTAs ({host}/proteins.faa)
  --out-dir    Output directory for DefenseFinder results
  -h, --help   Show this help message

Note: Requires DefenseFinder (pip install mdmparis-defense-finder)
      Update models before first run: defense-finder update
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --prot-dir) PROT_DIR="$2"; shift 2 ;;
        --out-dir)  OUT_DIR="$2";  shift 2 ;;
        -h|--help)  usage 0 ;;
        *)          echo "Unknown option: $1"; usage 1 ;;
    esac
done

if [[ -z "$PROT_DIR" || -z "$OUT_DIR" ]]; then
    echo "Error: --prot-dir and --out-dir are required."
    usage 1
fi

command -v defense-finder >/dev/null 2>&1 || { echo "Error: defense-finder not found in PATH"; exit 1; }

mkdir -p "$OUT_DIR"

mapfile -t FAA_FILES < <(find "$PROT_DIR" -name "proteins.faa" -type f | sort)
total=${#FAA_FILES[@]}
echo "Found $total protein file(s) under $PROT_DIR"

count=0
for faa in "${FAA_FILES[@]}"; do
    host=$(basename "$(dirname "$faa")")
    mkdir -p "$OUT_DIR/$host"

    defense-finder run -o "$OUT_DIR/$host/" "$faa"

    count=$((count + 1))
    if (( count % 50 == 0 )); then
        echo "  Progress: $count / $total hosts processed"
    fi
done

echo "Done. Processed $count host(s). Output in $OUT_DIR"
