#!/bin/bash
# 02_blastp_query.sh — Run diamond blastp: query test set vs cherry1306 protein DB.
#
# Test sets supported (TESTSET arg):
#   cherry634   — uses RefSeq-634 test phages, proteins already in data/Cherry_data/phage/<ACC>/proteins.faa
#   vhdb        — uses VHDB test phages, proteins to be merged from per-phage .faa
#
# Usage:
#   bash 02_blastp_query.sh <TESTSET>
#
# Inputs:
#   db/cherry1306_phage_prot.dmnd  (built by 01_build_db.py)
# Outputs (under runs/, queries/, logs/):
#   queries/<TESTSET>/proteins.faa
#   runs/<TESTSET>_vs_cherry1306.tsv
#   logs/02_blastp_query_<TESTSET>.log

set -euo pipefail

TESTSET="${1:-}"
if [[ -z "$TESTSET" ]]; then
  echo "usage: $0 <cherry634|vhdb>"; exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PIPELINE_DIR}/config.sh"

EXP_DIR="${PHI_PROJECT_ROOT}/experiments/blastp_cherry1306"
DIAMOND="${DIAMOND_BIN:-$(command -v diamond 2>/dev/null || echo diamond)}"
PRODIGAL="${PRODIGAL_BIN:-$(command -v prodigal 2>/dev/null || echo prodigal)}"

DB="$EXP_DIR/db/cherry1306_phage_prot.dmnd"
[[ -f "$DB" ]] || { echo "[ERR] DB missing: $DB — run 01_build_db.py first"; exit 1; }

QDIR="$EXP_DIR/queries/$TESTSET"
RDIR="$EXP_DIR/runs"
LDIR="$EXP_DIR/logs"
mkdir -p "$QDIR" "$RDIR" "$LDIR"

QFAA="$QDIR/proteins.faa"
LOG="$LDIR/02_blastp_query_${TESTSET}.log"
exec > >(tee -a "$LOG") 2>&1
echo "[$(date -Is)] BEGIN $TESTSET"

# ─── Step A: assemble query proteins ───────────────────────────────────────
case "$TESTSET" in
  cherry634)
    # RefSeq-634 test list and proteins (already prodigal'd)
    CHERRY_PHAGE_DIR="${PHI_PROJECT_ROOT}/data/Cherry_data/phage"
    TEST_LIST="${PHI_PROJECT_ROOT}/data/Cherry_data/cherry_test634_list.txt"
    [[ -f "$TEST_LIST" ]] || { echo "[ERR] test list missing: $TEST_LIST"; exit 1; }
    : > "$QFAA"
    n_phages=0; n_missing=0
    while read -r acc; do
      [[ -z "$acc" ]] && continue
      src="$CHERRY_PHAGE_DIR/$acc/proteins.faa"
      if [[ -f "$src" ]]; then
        awk -v acc="$acc" 'BEGIN{n=0}
          /^>/ {n++; print ">" acc "__gene" n; next}
          {print}' "$src" >> "$QFAA"
        n_phages=$((n_phages+1))
      else
        n_missing=$((n_missing+1))
      fi
    done < "$TEST_LIST"
    echo "  cherry634: $n_phages phages assembled, $n_missing missing"
    ;;
  vhdb)
    # VHDB test phages — call prodigal on each test phage genome
    VHDB_DIR="${PHI_PROJECT_ROOT}/experiments/baseline_inputs_vhdb/split_virus"
    [[ -d "$VHDB_DIR" ]] || { echo "[ERR] VHDB dir missing: $VHDB_DIR"; exit 1; }
    : > "$QFAA"
    tmpdir=$(mktemp -d)
    n_phages=0
    for fa in "$VHDB_DIR"/*.fa "$VHDB_DIR"/*.fasta; do
      [[ -e "$fa" ]] || continue
      acc=$(basename "$fa")
      acc="${acc%.fa}"; acc="${acc%.fasta}"
      tmp_faa="$tmpdir/${acc}.faa"
      "$PRODIGAL" -i "$fa" -a "$tmp_faa" -p meta -q -o /dev/null 2>>"$LOG"
      awk -v acc="$acc" 'BEGIN{n=0}
        /^>/ {n++; print ">" acc "__gene" n; next}
        {print}' "$tmp_faa" >> "$QFAA"
      n_phages=$((n_phages+1))
    done
    rm -rf "$tmpdir"
    echo "  vhdb: $n_phages phages assembled (prodigal -p meta)"
    ;;
  *)
    echo "[ERR] unknown TESTSET: $TESTSET"; exit 2;;
esac

n_seq=$(grep -c "^>" "$QFAA" || true)
echo "  query FAA: $QFAA  ($n_seq sequences)"

# ─── Step B: diamond blastp ────────────────────────────────────────────────
TSV="$RDIR/${TESTSET}_vs_cherry1306.tsv"
echo "[$(date -Is)] diamond blastp"
"$DIAMOND" blastp \
  --query    "$QFAA" \
  --db       "$DB" \
  --outfmt   6 qseqid sseqid pident length evalue bitscore \
  --evalue   1e-5 \
  --max-target-seqs 50 \
  --threads  32 \
  --block-size 2.0 \
  --index-chunks 4 \
  --out      "$TSV"
echo "[$(date -Is)] DONE — output: $TSV ($(wc -l < "$TSV") rows)"
