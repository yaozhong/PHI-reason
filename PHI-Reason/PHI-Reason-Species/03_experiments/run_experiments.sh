#!/bin/bash
# PHI-Reason-Species — Example experiment run queue
#
# Four experiments from the paper:
#   exp2c  RefSeq-634   RBP + BLASTN + BLASTP  (primary result)
#   exp5a  HiC Sp-1     (46 phages)
#   exp5b  HiC Sp-2     (82 phages)
#   exp8c  VHDB-3150    RBP + BLASTN + BLASTP
#
# Prerequisites:
#   1. Ollama running with the model loaded
#   2. Phage profiles built (01_profile_generation/)
#   3. Host list generated (01_profile_generation/host_profiles/)
#   4. Phage-host pair CSV prepared
#
# Usage:
#   export PHI_PROJECT_ROOT=/path/to/project
#   bash run_experiments.sh
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
source "${PIPELINE_DIR}/config.sh"

cd "${PHI_PROJECT_ROOT}"

EVAL="${PIPELINE_DIR}/02_inference/run_inference.py"
LOG="experiments/logs"
mkdir -p "$LOG"

run_exp() {
    local exp_id=$1; shift
    echo "========================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] START: $exp_id"
    echo "========================================"
    python3 -u "$EVAL" --exp-id "$exp_id" "$@" \
        >> "${LOG}/${exp_id}.log" 2>&1
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE:  $exp_id"
}

# ── exp2c: RefSeq-634, RBP + BLASTN + BLASTP ────────────────────────────────
run_exp exp2c_cherry_rbp_blastn_blastp \
    --host-list    ws/Cherry/host_list_v4G.txt \
    --pair-csv     data/Cherry_data/phage1940_host_pair.csv \
    --phage-list   ws/Cherry/phage_list_634.txt \
    --phage-prof-dir ws/Cherry/phage_profiles_rbp_blastn_blastp \
    --concurrency 4 --num-ctx 40960 --num-predict 4096

# ── exp5a: HiC Sp-1 (46 phages) ─────────────────────────────────────────────
run_exp exp5a_hic_sp1 \
    --host-list    ws/HiC/host_list_hic_paper52.txt \
    --pair-csv     ws/HiC/pair.csv \
    --phage-list   ws/HiC/phage_list_paper_sp1_46.txt \
    --phage-prof-dir ws/HiC/phage_profiles_rbp_blastn_sp1 \
    --concurrency 4 --num-ctx 40960 --num-predict 4096

# ── exp5b: HiC Sp-2 (82 phages) ─────────────────────────────────────────────
run_exp exp5b_hic_sp2 \
    --host-list    ws/HiC/host_list_hic_paper52.txt \
    --pair-csv     ws/HiC/pair.csv \
    --phage-list   ws/HiC/phage_list_paper_sp2_82.txt \
    --phage-prof-dir ws/HiC/phage_profiles_rbp_blastn_sp2 \
    --concurrency 4 --num-ctx 40960 --num-predict 4096

# ── exp8c: VHDB-3150, RBP + BLASTN + BLASTP ─────────────────────────────────
run_exp exp8c_vhdb_rbp_blastn_blastp \
    --host-list    ws/VHDB/host_list_v4G.txt \
    --pair-csv     ws/VHDB/pair.csv \
    --phage-list   ws/VHDB/phage_list.txt \
    --phage-prof-dir ws/VHDB/phage_profiles_rbp_blastn_blastp \
    --concurrency 4 --num-ctx 40960 --num-predict 4096

echo "[$(date '+%Y-%m-%d %H:%M:%S')] All experiments done."
