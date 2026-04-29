#!/bin/bash
# PHIreason core experiment run queue
#
# Four main experiments:
#   exp2c  Cherry-634    RBP + BLASTN + BLASTP  (primary result)
#   exp5a  HiC Sp-1 (46 phages)
#   exp5b  HiC Sp-2 (82 phages)
#   exp8c  VHDB          RBP + BLASTN + BLASTP
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

# ── exp2c: Cherry-634, RBP + BLASTN + BLASTP ─────────────────────────────────
run_exp exp2c_cherry_rbp_blastn_blastp \
    --host-list  ws/Cherry/host_list_v4G.txt \
    --pair-csv   "${PHI_PROJECT_ROOT}/data/Cherry_data/phage1940_host_pair.csv" \
    --phage-list ws/Cherry/phage_list_634.txt \
    --phage-prof-dir ws/Cherry/textGeneProfile_v2cm_R1_rbp_blastn_id_masked_blastp \
    --concurrency 12 --num-ctx 40960 --num-predict 4096

# ── exp5a: HiC Sp-1 (46 phages) ──────────────────────────────────────────────
run_exp exp5a_hic_paper_sp1 \
    --host-list  ws/HiC/host_list_hic_paper52.txt \
    --pair-csv   ws/HiC/pair.csv \
    --phage-list ws/HiC/phage_list_paper_sp1_46.txt \
    --phage-prof-dir ws/HiC/phage_profiles_rbp_blastn_sp1 \
    --concurrency 12 --num-ctx 40960 --num-predict 4096

# ── exp5b: HiC Sp-2 (82 phages) ──────────────────────────────────────────────
run_exp exp5b_hic_paper_sp2 \
    --host-list  ws/HiC/host_list_hic_paper52.txt \
    --pair-csv   ws/HiC/pair.csv \
    --phage-list ws/HiC/phage_list_paper_sp2_82.txt \
    --phage-prof-dir ws/HiC/phage_profiles_rbp_blastn_sp2 \
    --concurrency 12 --num-ctx 40960 --num-predict 4096

# ── exp8c: VHDB, RBP + BLASTN + BLASTP ───────────────────────────────────────
run_exp exp8c_vhdb_rbp_blastn_blastp \
    --host-list  ws/VHDB/host_list_v4G.txt \
    --pair-csv   ws/VHDB/pair.csv \
    --phage-list ws/VHDB/phage_list.txt \
    --phage-prof-dir ws/VHDB/phage_profiles_rbp_blastn_blastp_leak_refdb \
    --concurrency 12 --num-ctx 40960 --num-predict 4096 --think auto

echo "[$(date '+%Y-%m-%d %H:%M:%S')] All 4 experiments done."
