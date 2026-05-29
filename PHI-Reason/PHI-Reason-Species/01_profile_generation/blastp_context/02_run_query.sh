#!/bin/bash
# 02_run_query.sh — Run diamond blastp: query test-set proteins vs train-set DB.
#
# Usage:
#   bash 02_run_query.sh \
#     --db       <path/to/train_phage_prot.dmnd> \
#     --phage-dir <dir_with_per-phage_subdirs>   \
#     --phage-list <test_phage_ids.txt>          \
#     --out-dir  <output_directory>              \
#     [--threads 8] [--use-prodigal]
#
# Each phage subdir should contain proteins.faa. If --use-prodigal is set,
# phage subdirs should contain genome.fa and prodigal will be called to
# generate proteins on the fly (for datasets without pre-computed proteins).

set -euo pipefail

# ─── Argument parsing ─────────────────────────────────────────────────────
DB="" ; PHAGE_DIR="" ; PHAGE_LIST="" ; OUT_DIR="" ; THREADS=8 ; USE_PRODIGAL=false

usage() {
  echo "Usage: $0 --db <dmnd> --phage-dir <dir> --phage-list <txt> --out-dir <dir> [--threads N] [--use-prodigal]"
  exit 2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --db)          DB="$2"; shift 2 ;;
    --phage-dir)   PHAGE_DIR="$2"; shift 2 ;;
    --phage-list)  PHAGE_LIST="$2"; shift 2 ;;
    --out-dir)     OUT_DIR="$2"; shift 2 ;;
    --threads)     THREADS="$2"; shift 2 ;;
    --use-prodigal) USE_PRODIGAL=true; shift ;;
    -h|--help)     usage ;;
    *)             echo "[ERR] unknown arg: $1"; usage ;;
  esac
done

[[ -z "$DB" || -z "$PHAGE_DIR" || -z "$PHAGE_LIST" || -z "$OUT_DIR" ]] && usage

DIAMOND="${DIAMOND_BIN:-$(command -v diamond 2>/dev/null || echo diamond)}"
PRODIGAL="${PRODIGAL_BIN:-$(command -v prodigal 2>/dev/null || echo prodigal)}"

[[ -f "$DB" ]] || { echo "[ERR] DB missing: $DB — run 01_build_db.py first"; exit 1; }
[[ -f "$PHAGE_LIST" ]] || { echo "[ERR] phage list missing: $PHAGE_LIST"; exit 1; }

mkdir -p "$OUT_DIR"
QFAA="$OUT_DIR/query_proteins.faa"
TSV="$OUT_DIR/blastp_hits.tsv"
LOG="$OUT_DIR/02_run_query.log"
exec > >(tee -a "$LOG") 2>&1
echo "[$(date -Is)] BEGIN"

# ─── Step A: assemble query proteins ──────────────────────────────────────
: > "$QFAA"
n_phages=0; n_missing=0

while read -r acc; do
  [[ -z "$acc" ]] && continue

  if [[ "$USE_PRODIGAL" == true ]]; then
    genome="$PHAGE_DIR/$acc/genome.fa"
    if [[ ! -f "$genome" ]]; then
      # try alternative extensions
      for ext in fna fasta; do
        [[ -f "$PHAGE_DIR/$acc/genome.$ext" ]] && genome="$PHAGE_DIR/$acc/genome.$ext" && break
      done
    fi
    if [[ -f "$genome" ]]; then
      tmpfaa=$(mktemp)
      "$PRODIGAL" -i "$genome" -a "$tmpfaa" -p meta -q -o /dev/null 2>>"$LOG"
      awk -v acc="$acc" 'BEGIN{n=0} /^>/ {n++; print ">" acc "__gene" n; next} {print}' "$tmpfaa" >> "$QFAA"
      rm -f "$tmpfaa"
      n_phages=$((n_phages+1))
    else
      n_missing=$((n_missing+1))
    fi
  else
    src="$PHAGE_DIR/$acc/proteins.faa"
    if [[ -f "$src" ]]; then
      awk -v acc="$acc" 'BEGIN{n=0} /^>/ {n++; print ">" acc "__gene" n; next} {print}' "$src" >> "$QFAA"
      n_phages=$((n_phages+1))
    else
      n_missing=$((n_missing+1))
    fi
  fi
done < "$PHAGE_LIST"

echo "  assembled: $n_phages phages, $n_missing missing"
n_seq=$(grep -c "^>" "$QFAA" || true)
echo "  query FAA: $QFAA ($n_seq sequences)"

# ─── Step B: diamond blastp ───────────────────────────────────────────────
echo "[$(date -Is)] diamond blastp"
"$DIAMOND" blastp \
  --query    "$QFAA" \
  --db       "$DB" \
  --outfmt   6 qseqid sseqid pident length evalue bitscore \
  --evalue   1e-5 \
  --max-target-seqs 50 \
  --threads  "$THREADS" \
  --block-size 2.0 \
  --index-chunks 4 \
  --out      "$TSV"
echo "[$(date -Is)] DONE — output: $TSV ($(wc -l < "$TSV") rows)"
