#!/usr/bin/env bash
# run_blastn.sh — Whole-genome BLASTN for phage-host prediction pipeline
# Requires: BLAST+ (makeblastdb, blastn) available in PATH
set -euo pipefail

TRAIN_DIR=""; QUERY_DIR=""; QUERY_LIST=""; OUT_DIR=""; THREADS=8

usage() {
    cat <<EOF
Usage: bash $0 --train-genomes <dir> --query-genomes <dir> --query-list <txt> --out-dir <dir> [--threads 8]
  --train-genomes  Directory of training phage genomes (one FASTA per phage)
  --query-genomes  Directory of test phage genomes
  --query-list     Text file listing test phage IDs (one per line)
  --out-dir        Output directory (produces blastn_hits.tsv)
  --threads        CPU threads (default: 8)
  -h, --help       Show this help
Genome lookup: <dir>/<acc>/genome.fa > <acc>.fasta > <acc>.fna > <acc>.fa
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --train-genomes) TRAIN_DIR="$2";  shift 2 ;;
        --query-genomes) QUERY_DIR="$2";  shift 2 ;;
        --query-list)    QUERY_LIST="$2"; shift 2 ;;
        --out-dir)       OUT_DIR="$2";    shift 2 ;;
        --threads)       THREADS="$2";    shift 2 ;;
        -h|--help)       usage 0 ;;
        *)               echo "Unknown option: $1"; usage 1 ;;
    esac
done

if [[ -z "$TRAIN_DIR" || -z "$QUERY_DIR" || -z "$QUERY_LIST" || -z "$OUT_DIR" ]]; then
    echo "Error: --train-genomes, --query-genomes, --query-list, and --out-dir are required."
    usage 1
fi
for cmd in makeblastdb blastn; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "Error: $cmd not found in PATH"; exit 1; }
done

mkdir -p "$OUT_DIR"

# ── 1. Build BLAST database from training genomes ──────────────────────────
COMBINED="$OUT_DIR/train_combined.fna"
: > "$COMBINED"
n_train=0
for f in $(find -L "$TRAIN_DIR" -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) | sort); do
    cat "$f" >> "$COMBINED"
    n_train=$((n_train + 1))
done
echo "Concatenated $n_train training genome file(s) into $COMBINED"

makeblastdb -in "$COMBINED" -dbtype nucl -out "$OUT_DIR/train_db"
echo "BLAST database built: $OUT_DIR/train_db"

# ── 2. Find genome file for a given accession ─────────────────────────────
find_genome() {
    local dir="$1" acc="$2"
    for candidate in "$dir/$acc/genome.fa" "$dir/$acc.fasta" "$dir/$acc.fna" "$dir/$acc.fa"; do
        [[ -f "$candidate" ]] && echo "$candidate" && return 0
    done
    return 1
}

# ── 3. Run BLASTN for each query phage ─────────────────────────────────────
OUTFILE="$OUT_DIR/blastn_hits.tsv"
: > "$OUTFILE"
n_query=0 n_skip=0
while IFS= read -r acc || [[ -n "$acc" ]]; do
    acc="${acc%%[[:space:]]*}"
    [[ -z "$acc" || "$acc" == \#* ]] && continue
    genome=$(find_genome "$QUERY_DIR" "$acc") || { echo "  SKIP $acc: genome not found"; n_skip=$((n_skip+1)); continue; }
    blastn -query "$genome" -db "$OUT_DIR/train_db" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
        -max_target_seqs 10 -evalue 1e-5 -num_threads "$THREADS" \
        >> "$OUTFILE"
    n_query=$((n_query + 1))
done < "$QUERY_LIST"
echo "Done. Queried $n_query phage(s), skipped $n_skip. Output: $OUTFILE"
