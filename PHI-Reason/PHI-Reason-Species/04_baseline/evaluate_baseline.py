#!/usr/bin/env python3
"""
Unified evaluator for baseline re-runs on RefSeq-634 under PHIreason_final/experiments/.

Imports parsers + compute_metrics from the original evaluate_full.py so metric
definitions stay identical. Output path for each tool is passed via CLI.

Usage:
  python evaluate_baseline.py --tool php      --outputs <exp_dir>/outputs
  python evaluate_baseline.py --tool wish     --outputs <exp_dir>/outputs
  python evaluate_baseline.py --tool phist    --outputs <exp_dir>/outputs
  python evaluate_baseline.py --tool deephost --outputs <exp_dir>/outputs
  python evaluate_baseline.py --tool phabox2  --outputs <exp_dir>/outputs  --mode db|mag
  python evaluate_baseline.py --tool vhmnet   --outputs <exp_dir>/outputs

Emits <exp_dir>/results/METRICS.json and prints a compact summary.
"""
from __future__ import annotations
import argparse, json, os, sys, importlib.util, csv, re

_DEFAULT_REVIEW_ROOT = os.environ.get("HOST_PREDICTION_REVIEW_DIR", "")

def _load_eval_module(eval_full_path: str):
    spec = importlib.util.spec_from_file_location("eval_full", eval_full_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def parse_phist(predictions_csv, ef):
    """
    PHIST predictions.csv columns:
      phage, host, common kmers, <others...>
    Rank hosts by 'common kmers' descending per phage. Top-1 = host with most shared kmers.
    Returns score_matrix[vid] = {host_norm: score}.
    """
    scores = {}
    if not os.path.exists(predictions_csv):
        return scores
    with open(predictions_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            phage = row.get("phage") or row.get("Phage") or ""
            host  = row.get("host")  or row.get("Host")  or ""
            kmers = (row.get("common kmers") or row.get("#common_kmer")
                     or row.get("common_kmers") or row.get("#common-kmers") or "0")
            vid = ef.normalize_vid(phage.strip().split()[0]) if phage else ""
            if not vid or not host:
                continue
            try:
                k = float(kmers)
            except ValueError:
                continue
            host_stripped = re.sub(r"\.fa(sta)?$", "", host.strip(), flags=re.IGNORECASE)
            host_norm = ef.normalize(host_stripped)
            scores.setdefault(vid, {})
            if host_norm not in scores[vid] or scores[vid][host_norm] < k:
                scores[vid][host_norm] = k
    return scores

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True,
                    choices=["php","wish","phist","deephost","phabox2","vhmnet"])
    ap.add_argument("--outputs", required=True, help="<exp_dir>/outputs directory")
    ap.add_argument("--mode", choices=["db","mag"], default="db")
    ap.add_argument("--results", default=None, help="<exp_dir>/results output dir (default: ../results)")
    ap.add_argument("--review-dir", default=_DEFAULT_REVIEW_ROOT,
                    help="Path to HostPredictionReview-main directory "
                         "(or set HOST_PREDICTION_REVIEW_DIR env var)")
    args = ap.parse_args()

    eval_full = os.path.join(args.review_dir,
                             "CHERRY_benchmark_results", "evaluate_full.py")
    datasplit = os.path.join(args.review_dir, "CHERRY_benchmark_datasplit")

    if not args.review_dir:
        sys.exit("[ERR] Set --review-dir or HOST_PREDICTION_REVIEW_DIR env var to the "
                 "HostPredictionReview directory.\n"
                 "  git clone https://github.com/KennthShang/HostPredictionReview.git")
    if not os.path.isfile(eval_full):
        sys.exit(f"[ERR] evaluate_full.py not found: {eval_full}\n"
                 f"  Set --review-dir or HOST_PREDICTION_REVIEW_DIR to the "
                 f"HostPredictionReview directory.")

    ef = _load_eval_module(eval_full)

    virus_ids    = ef.load_test_fasta_ids(os.path.join(datasplit, "CHERRY_test.fasta"))
    ground_truth = ef.load_ground_truth(os.path.join(datasplit, "CHERRY_y_test.csv"), virus_ids)

    outputs = os.path.abspath(args.outputs)
    results_dir = args.results or os.path.join(os.path.dirname(outputs), "results")
    os.makedirs(results_dir, exist_ok=True)

    tool = args.tool
    if tool == "php":
        p = os.path.join(outputs, "cherry_host_kmer4_Prediction_Allhost.csv")
        sm = ef.parse_PHP_matrix(p)
        r  = ef.compute_metrics("PHP", sm, ground_truth, virus_ids)

    elif tool == "wish":
        p = os.path.join(outputs, "predictions", "llikelihood.matrix")
        sm = ef.parse_WIsH_matrix(p)
        r  = ef.compute_metrics("WIsH", sm, ground_truth, virus_ids)

    elif tool == "phist":
        p = os.path.join(outputs, "predictions.csv")
        sm = parse_phist(p, ef)
        split_host_dir = os.path.join(os.path.dirname(outputs), "inputs", "split_host")
        host_universe = set()
        if os.path.isdir(split_host_dir):
            for fn in os.listdir(split_host_dir):
                name = re.sub(r"\.fa(sta)?$", "", fn, flags=re.IGNORECASE)
                host_universe.add(ef.normalize(name))
        r  = ef.compute_metrics("PHIST", sm, ground_truth, virus_ids,
                                n_labels=len(host_universe) if host_universe else len(ground_truth),
                                gt_labels_set=host_universe or None)

    elif tool == "deephost":
        p_multi = os.path.join(outputs, "DeepHost_multihost.tsv")
        sm = ef.parse_DeepHost_multi(p_multi)
        r  = ef.compute_metrics("DeepHost", sm, ground_truth, virus_ids)

    elif tool == "phabox2":
        tsv = os.path.join(outputs, "results", "final_prediction", "cherry_prediction.tsv")
        if not os.path.exists(tsv):
            tsv = os.path.join(outputs, "final_prediction", "cherry_prediction.tsv")
        top1 = {}
        with open(tsv) as f:
            next(f)
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 3: continue
                vid = ef.normalize_vid(parts[0].strip())
                raw = parts[2].strip()
                if raw and raw != "-":
                    host = ef.normalize(raw.split(":",1)[1].strip() if ":" in raw else raw)
                    if host: top1[vid] = host
        label = f"phabox2 ({args.mode.upper()} mode)"
        r = ef.compute_metrics(label, top1, ground_truth, virus_ids,
                               top1_only=True, n_labels=len(ground_truth))

    elif tool == "vhmnet":
        pred_dir = os.path.join(outputs, "output", "predictions")
        if not os.path.isdir(pred_dir):
            pred_dir = os.path.join(outputs, "predictions")
        sm = ef.parse_VHMnet_top30(pred_dir)
        r  = ef.compute_metrics("VHMnet (top-30)", sm, ground_truth, virus_ids,
                                n_labels=len(ground_truth))
        r["AUC"]  = "N/A"
        r["AUPR"] = "N/A"

    metrics_json = os.path.join(results_dir, "METRICS.json")
    with open(metrics_json, "w") as f:
        json.dump(r, f, indent=2, ensure_ascii=False)

    print(json.dumps(r, indent=2, ensure_ascii=False))
    print(f"\n[OK] Wrote {metrics_json}")

if __name__ == "__main__":
    main()
