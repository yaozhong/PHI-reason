#!/usr/bin/env python3
"""
eval_multi.py — parameterized driver for Cherry / HiC / VHDB evaluations.
=========================================================================
Differences from eval_v4G.py:
  - prompt imported from prompt/prompt_template.py (not inlined)
  - --host-list and --pair-csv are CLI args (not hardcoded)
  - NUM_CTX defaults to 40960
  - pair CSV supports optional truth_level column ("species" | "genus")
    → species metrics computed only over species-truth rows
    → genus metrics computed over all rows

Usage (run from project root, or set PHI_PROJECT_ROOT):
  export PHI_PROJECT_ROOT=/path/to/your/project
  python3 PHIreason_pipeline/02_inference/run_inference.py \\
      --exp-id       exp1_cherry_base \\
      --host-list    ws/Cherry/host_list_v4G.txt \\
      --pair-csv     data/Cherry_data/phage1940_host_pair.csv \\
      --phage-list   data/Cherry_data/cherry_test634_list.txt \\
      --phage-prof-dir ws/Cherry/textGeneProfile_v2cm_R1
"""
from __future__ import annotations

import argparse
import asyncio
import csv
import json
import os
import re
import sys
import time
from collections import defaultdict
from datetime import datetime
from json import JSONDecoder
from pathlib import Path

import aiohttp
import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

_AUTO_ROOT = Path(__file__).resolve().parent.parent.parent
FINAL_DIR = Path(os.environ.get("PHI_PROJECT_ROOT", str(_AUTO_ROOT)))
WS_DIR    = FINAL_DIR / "ws"
EXP_ROOT  = FINAL_DIR / "experiments"

sys.path.insert(0, str(Path(__file__).resolve().parent / "prompt"))
from prompt_template import SYSTEM_PROMPT, USER_PREFIX_TEMPLATE, USER_VARIABLE_TEMPLATE, TOP_K

DEFAULT_MODEL = os.environ.get("PHI_MODEL",       "qwen3-coder-next:q4_K_M")
DEFAULT_URLS  = os.environ.get("PHI_OLLAMA_URLS", "http://127.0.0.1:11435,http://127.0.0.1:11436")
DEFAULT_CONC  = int(os.environ.get("PHI_CONCURRENCY", "12"))
NUM_CTX       = 40960
NUM_PREDICT   = 4096
SELFCHECK_EVERY = 50


def load_pair_csv(path: Path) -> dict[str, tuple[str, str]]:
    """Return {phage_id: (true_host, truth_level)}.

    Supports two CSV flavors:
      A) phage,host                              → truth_level = "species"
      B) phage_id,true_host,truth_level          → per-row
    """
    out: dict[str, tuple[str, str]] = {}
    with path.open() as f:
        first = f.readline()
        cols = [c.strip() for c in first.split(",")]
        has_level = len(cols) >= 3 and cols[2].lower() in ("truth_level", "level")
        header_row = any(c.lower() in ("phage_id", "phage", "accession") for c in cols)
        if not header_row:
            f.seek(0)
        for line in f:
            parts = [p.strip() for p in line.rstrip("\n").split(",")]
            if len(parts) < 2 or not parts[0]:
                continue
            pid = parts[0]
            host = parts[1]
            level = parts[2] if has_level and len(parts) >= 3 else "species"
            out[pid] = (host, level)
    return out


def load_host_names(host_list_path: Path) -> list[str]:
    names = []
    for line in host_list_path.read_text().splitlines():
        m = re.match(r"^\d+\.\s+(\S+)", line.strip())
        if m:
            names.append(m.group(1))
    return sorted(names)


def load_phage_profile(prof_dir: Path, phage_id: str) -> str | None:
    p = prof_dir / f"{phage_id}.md"
    return p.read_text(errors="replace") if p.exists() else None


def parse_response(raw: str, phage_id: str, host_names: list[str],
                   num_predict: int, gen_tokens: int) -> dict:
    text = raw.strip()
    text = re.sub(r"<think>.*?</think>", "", text, flags=re.DOTALL).strip()
    for tok in ("<|endoftext|>", "<|im_start|>", "<|im_end|>"):
        if tok in text:
            text = text.split(tok, 1)[0].strip()

    scores = {h: 0.0 for h in host_names}
    reasoning = gram_decision = ""
    cutoff = gen_tokens >= num_predict - 50

    json_start = text.find("{")
    payload = None
    parse_err = False
    if json_start != -1:
        try:
            payload, _ = JSONDecoder().raw_decode(text[json_start:])
        except Exception:
            parse_err = True
            cutoff = True
            truncated = text[json_start:]
            m = re.search(r'"reasoning"\s*:\s*"([^"]*)"', truncated)
            if m: reasoning = m.group(1)
            m = re.search(r'"gram_type_decision"\s*:\s*"([^"]*)"', truncated)
            if m: gram_decision = m.group(1)
            for m in re.finditer(
                r'"host"\s*:\s*"([^"]+)"[^}]*?"score"\s*:\s*([\d.]+)', truncated
            ):
                h = m.group(1).strip().replace(" ", "_")
                try:
                    s = float(m.group(2))
                    if h in scores:
                        scores[h] = max(scores[h], min(1.0, max(0.0, s)))
                except Exception:
                    pass
    else:
        parse_err = True

    if isinstance(payload, dict):
        reasoning = payload.get("reasoning", "")
        gram_decision = payload.get("gram_type_decision", "")
        for row in payload.get("predictions", []):
            try:
                h = str(row.get("host", "")).strip().replace(" ", "_")
                s = float(row.get("score", 0.0))
                if h in scores:
                    scores[h] = max(0.0, min(1.0, s))
            except Exception:
                pass

    no_pred = not any(v > 0 for v in scores.values())
    err = None
    if no_pred and parse_err:
        err = "parse_fail_no_predictions"
    elif no_pred:
        err = "no_predictions"
    elif cutoff:
        err = "cutoff_partial"

    return {
        "phage_id": phage_id, "gram_decision": gram_decision,
        "reasoning": reasoning, "scores": scores, "raw_response": raw,
        "error": err, "_cutoff": cutoff, "_parse_err": parse_err,
    }


async def call_ollama(session, url, system, user, model, num_predict, num_ctx, think=None):
    prompt = user
    if think is True:
        prompt = "/think\n" + user
    elif think is False:
        prompt = "/no_think\n" + user
    payload = {
        "model": model, "prompt": prompt, "system": system,
        "stream": False, "keep_alive": "120m",
        "options": {
            "temperature": 0.1, "num_predict": num_predict, "num_ctx": num_ctx,
            "num_gpu": 99,
            "stop": ["<|endoftext|>", "<|im_start|>", "<|im_end|>"],
        },
    }
    async with session.post(f"{url}/api/generate", json=payload,
                            timeout=aiohttp.ClientTimeout(total=1800)) as resp:
        return await resp.json()


async def infer_one(session, url_queue, phage_id, prof_dir, user_prefix,
                    host_names, model, num_predict, num_ctx,
                    cache_dir, input_dir, counter, total, t_start, think=None):
    safe = phage_id.replace("/", "_")
    cache_file = cache_dir / f"{safe}.json"
    input_file = input_dir / f"{safe}.txt"

    if cache_file.exists():
        try:
            r = json.loads(cache_file.read_text())
            if not input_file.exists():
                pt = load_phage_profile(prof_dir, phage_id)
                if pt:
                    um = user_prefix + USER_VARIABLE_TEMPLATE.format(phage_profile=pt)
                    input_file.write_text(
                        f"=== SYSTEM PROMPT ===\n{SYSTEM_PROMPT}\n=== USER MESSAGE ===\n{um}",
                        encoding="utf-8")
            counter[0] += 1
            return r
        except Exception:
            pass

    t0 = time.perf_counter()
    pt = load_phage_profile(prof_dir, phage_id)
    if pt is None:
        r = {"phage_id": phage_id, "gram_decision": "", "reasoning": "",
             "scores": {h: 0.0 for h in host_names}, "raw_response": "",
             "error": "missing_profile", "_cutoff": False, "_parse_err": False,
             "_prompt_tokens": 0, "_gen_tokens": 0,
             "_llm_process_s": 0.0, "_elapsed_s": 0.0}
        cache_file.write_text(json.dumps(r, ensure_ascii=False, indent=2))
        counter[0] += 1
        return r

    user_msg = user_prefix + USER_VARIABLE_TEMPLATE.format(phage_profile=pt)
    if not input_file.exists():
        input_file.write_text(
            f"=== SYSTEM PROMPT ===\n{SYSTEM_PROMPT}\n=== USER MESSAGE ===\n{user_msg}",
            encoding="utf-8")

    url = await url_queue.get()
    try:
        data = await call_ollama(session, url, SYSTEM_PROMPT, user_msg,
                                 model, num_predict, num_ctx, think=think)
        raw = data.get("response", "")
        gen_tokens = data.get("eval_count", 0)
        r = parse_response(raw, phage_id, host_names, num_predict, gen_tokens)
        r["_llm_prefill_s"] = round(data.get("prompt_eval_duration", 0) / 1e9, 2)
        r["_llm_gen_s"]     = round(data.get("eval_duration", 0) / 1e9, 2)
        r["_prompt_tokens"] = data.get("prompt_eval_count", 0)
        r["_gen_tokens"]    = gen_tokens
        r["_llm_process_s"] = round(r["_llm_prefill_s"] + r["_llm_gen_s"], 2)
    except Exception as exc:
        r = {"phage_id": phage_id, "gram_decision": "", "reasoning": "",
             "scores": {h: 0.0 for h in host_names}, "raw_response": "",
             "error": str(exc), "_cutoff": False, "_parse_err": False,
             "_prompt_tokens": 0, "_gen_tokens": 0, "_llm_process_s": 0.0}
    finally:
        await url_queue.put(url)

    r["_elapsed_s"] = round(time.perf_counter() - t0, 1)
    cache_file.write_text(json.dumps(r, ensure_ascii=False, indent=2))
    counter[0] += 1

    done = counter[0]
    if done % SELFCHECK_EVERY == 0 or done == total:
        elapsed = time.perf_counter() - t_start
        rate = done / elapsed if elapsed > 0 else 0
        eta = (total - done) / rate if rate > 0 else 0
        sample = list(cache_dir.glob("*.json"))[-min(200, done):]
        nc = npe = ne = gs = 0
        for cf in sample:
            try:
                d = json.loads(cf.read_text())
                if d.get("_cutoff"):    nc  += 1
                if d.get("_parse_err"): npe += 1
                if d.get("error"):      ne  += 1
                gs += d.get("_gen_tokens", 0)
            except Exception: pass
        avg_gen = gs // len(sample) if sample else 0
        print(f"  [自查 {done:4d}/{total}] 速率={rate:.2f}/s ETA={eta/60:.1f}min "
              f"| sample{len(sample)}: err={ne} parse_err={npe} cutoff={nc} avg_gen={avg_gen}tok",
              flush=True)
    return r


async def infer_batch(phage_ids, prof_dir, cache_dir, input_dir, host_names,
                      user_prefix, urls, model, num_predict, num_ctx,
                      concurrency, t_start, think=None):
    slots = max(1, concurrency // len(urls))
    q: asyncio.Queue = asyncio.Queue()
    for u in urls:
        for _ in range(slots):
            await q.put(u)
    counter = [0]
    async with aiohttp.ClientSession() as s:
        tasks = [infer_one(s, q, pid, prof_dir, user_prefix, host_names,
                           model, num_predict, num_ctx, cache_dir, input_dir,
                           counter, len(phage_ids), t_start, think=think)
                 for pid in phage_ids]
        return await asyncio.gather(*tasks)


def evaluate(results, true_pairs: dict[str, tuple[str, str]]) -> dict:
    """Compute species + genus metrics separately.

    species metrics: only over phages with truth_level='species'.
    genus metrics:   over all phages with a true_host in scores.
    """
    sp_top1, sp_top5, sp_top10, sp_top30 = [], [], [], []
    gn_top1, gn_top5, gn_top10, gn_top30 = [], [], [], []
    sp_ranks, sp_rr, gn_rr = [], [], []
    labels, flat = [], []

    for row in results:
        pid = row["phage_id"]
        pair = true_pairs.get(pid)
        if pair is None:
            continue
        true_host, level = pair
        if true_host not in row["scores"]:
            continue
        ranked = sorted(row["scores"].items(), key=lambda kv: (-kv[1], kv[0]))
        order  = [h for h, _ in ranked]
        rank   = order.index(true_host) + 1

        true_genus = true_host.split("_")[0]
        genus_rank = min(
            (i+1 for i, h in enumerate(order) if h.split("_")[0] == true_genus),
            default=len(order)+1)

        if level == "species":
            sp_ranks.append(rank); sp_rr.append(1.0/rank)
            sp_top1.append(int(rank<=1));  sp_top5.append(int(rank<=5))
            sp_top10.append(int(rank<=10)); sp_top30.append(int(rank<=30))
            for h, s in ranked:
                labels.append(int(h == true_host)); flat.append(s)

        gn_rr.append(1.0/genus_rank)
        gn_top1.append(int(genus_rank<=1));  gn_top5.append(int(genus_rank<=5))
        gn_top10.append(int(genus_rank<=10)); gn_top30.append(int(genus_rank<=30))

    def avg(lst): return float(np.mean(lst)) if lst else 0.0
    roc = float(roc_auc_score(labels, flat)) if labels and len(set(labels)) > 1 else 0.0
    prc = float(average_precision_score(labels, flat)) if labels and len(set(labels)) > 1 else 0.0

    return {
        "n_species_truth": len(sp_ranks),
        "n_genus_truth":   len(gn_rr),
        "species_top1":  round(avg(sp_top1), 4),  "species_top5":  round(avg(sp_top5), 4),
        "species_top10": round(avg(sp_top10), 4), "species_top30": round(avg(sp_top30), 4),
        "species_mrr":   round(avg(sp_rr), 4),
        "species_mean_rank":   round(float(np.mean(sp_ranks)) if sp_ranks else 0.0, 4),
        "species_median_rank": round(float(np.median(sp_ranks)) if sp_ranks else 0.0, 4),
        "genus_top1":  round(avg(gn_top1), 4),  "genus_top5":  round(avg(gn_top5), 4),
        "genus_top10": round(avg(gn_top10), 4), "genus_top30": round(avg(gn_top30), 4),
        "genus_mrr":   round(avg(gn_rr), 4),
        "roc_auc": round(roc, 4), "pr_auc": round(prc, 4),
    }


def write_outputs(results_dir, results, true_pairs):
    rows = []
    for r in results:
        pair = true_pairs.get(r["phage_id"])
        true_host = pair[0] if pair else ""
        for h, s in r["scores"].items():
            rows.append((r["phage_id"], h, int(h == true_host), s))
    with (results_dir / "all_scores.csv").open("w") as f:
        f.write("phage_id,host_id,label,score\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]:.6f}\n")
    pid_rows: dict[str, list] = defaultdict(list)
    for r in rows: pid_rows[r[0]].append(r)
    with (results_dir / "top_predictions.tsv").open("w") as f:
        f.write("phage_id\thost_id\tlabel\tscore\trank\n")
        for pid, pr in sorted(pid_rows.items()):
            for rank, (phage, host, label, score) in enumerate(
                    sorted(pr, key=lambda x: -x[3]), 1):
                if rank <= TOP_K:
                    f.write(f"{phage}\t{host}\t{label}\t{score:.6f}\t{rank}\n")


def write_summary(results_dir, exp_id, m, model, n_phages,
                  elapsed_h, n_err, n_cutoff):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M")
    md = f"""# 实验结果汇总 — {exp_id}

**时间**: {ts}
**模型**: {model}
**噬菌体数**: {n_phages}
**耗时**: {elapsed_h:.2f}h
**错误/截断**: errors={n_err}  cutoffs={n_cutoff}
**样本数**: species-truth={m['n_species_truth']}  genus-truth={m['n_genus_truth']}

## 评价指标

| 指标 | 值 |
|------|----|
| species_top1 | {m['species_top1']} |
| species_top5 | {m['species_top5']} |
| species_top10 | {m['species_top10']} |
| species_top30 | {m['species_top30']} |
| species_mrr | {m['species_mrr']} |
| species_mean_rank | {m['species_mean_rank']} |
| species_median_rank | {m['species_median_rank']} |
| genus_top1 | {m['genus_top1']} |
| genus_top5 | {m['genus_top5']} |
| genus_top10 | {m['genus_top10']} |
| genus_top30 | {m['genus_top30']} |
| genus_mrr | {m['genus_mrr']} |
| roc_auc (species) | {m['roc_auc']} |
| pr_auc (species) | {m['pr_auc']} |
"""
    (results_dir / "RESULT_SUMMARY.md").write_text(md)
    (results_dir / "metrics.json").write_text(
        json.dumps({"exp_id": exp_id, "model": model, "n_phages": n_phages,
                    "elapsed_h": round(elapsed_h, 3), **m}, indent=2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-id",          required=True)
    ap.add_argument("--host-list",       required=True, type=Path)
    ap.add_argument("--pair-csv",        required=True, type=Path)
    ap.add_argument("--phage-list",      required=True, type=Path)
    ap.add_argument("--phage-prof-dir",  required=True, type=Path)
    ap.add_argument("--model",           default=DEFAULT_MODEL)
    ap.add_argument("--ollama-urls",     default=DEFAULT_URLS)
    ap.add_argument("--concurrency",     type=int, default=DEFAULT_CONC)
    ap.add_argument("--num-predict",     type=int, default=NUM_PREDICT)
    ap.add_argument("--num-ctx",         type=int, default=NUM_CTX)
    ap.add_argument("--pilot",           type=int, default=0)
    ap.add_argument("--overwrite",       action="store_true")
    ap.add_argument("--think",           choices=["yes", "no", "auto"], default="auto")
    args = ap.parse_args()

    exp_dir     = EXP_ROOT / args.exp_id
    results_dir = exp_dir / "results"
    cache_dir   = results_dir / "cache"
    input_dir   = exp_dir / "data" / "inputs"
    results_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)

    for p, label in [(args.host_list, "host-list"), (args.pair_csv, "pair-csv"),
                     (args.phage_list, "phage-list"), (args.phage_prof_dir, "prof-dir")]:
        if not p.exists():
            raise SystemExit(f"ERROR: {label} not found: {p}")

    urls = [u.strip() for u in args.ollama_urls.split(",") if u.strip()]
    true_pairs = load_pair_csv(args.pair_csv)
    host_names = load_host_names(args.host_list)
    host_list_text = args.host_list.read_text()
    user_prefix = USER_PREFIX_TEMPLATE.format(n_hosts=len(host_names),
                                               host_list=host_list_text)

    phage_ids = [l.strip() for l in args.phage_list.read_text().splitlines() if l.strip()]
    if args.pilot > 0:
        phage_ids = phage_ids[:args.pilot]

    missing = [pid for pid in phage_ids if not (args.phage_prof_dir / f"{pid}.md").exists()]
    if missing:
        print(f"  ⚠️  {len(missing)} 个 phage 缺 profile: {missing[:5]}...")
    phage_ids = [pid for pid in phage_ids if (args.phage_prof_dir / f"{pid}.md").exists()]

    if args.overwrite:
        for pid in phage_ids:
            (cache_dir / f"{pid.replace('/','_')}.json").unlink(missing_ok=True)

    print("=" * 72)
    print(f"  PHIreason_final · {args.exp_id}")
    print(f"  时间   : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  模型   : {args.model}")
    print(f"  噬菌体 : {len(phage_ids)}  宿主: {len(host_names)}")
    print(f"  并发   : {args.concurrency}  URLs: {urls}")
    print(f"  num_ctx={args.num_ctx}  num_predict={args.num_predict}  think={args.think}")
    print(f"  prof dir : {args.phage_prof_dir}")
    print(f"  host list: {args.host_list} ({len(host_list_text)/1024:.1f} KB)")
    print(f"  pair csv : {args.pair_csv}  ({len(true_pairs)} entries)")
    print("=" * 72)

    think_flag = {"yes": True, "no": False, "auto": None}[args.think]
    t0 = time.perf_counter()
    results = asyncio.run(infer_batch(
        phage_ids=phage_ids, prof_dir=args.phage_prof_dir,
        cache_dir=cache_dir, input_dir=input_dir,
        host_names=host_names, user_prefix=user_prefix,
        urls=urls, model=args.model,
        num_predict=args.num_predict, num_ctx=args.num_ctx,
        concurrency=args.concurrency, t_start=t0, think=think_flag))
    elapsed_h = (time.perf_counter() - t0) / 3600

    n_err    = sum(1 for r in results if r.get("error") and r["error"] != "cutoff_partial")
    n_cutoff = sum(1 for r in results if r.get("_cutoff"))

    m = evaluate(results, true_pairs)

    print("\n" + "=" * 72)
    print(f"  完成: {len(results)}  耗时: {elapsed_h:.2f}h")
    print(f"  errors={n_err}  cutoffs={n_cutoff}")
    print(f"  样本: species={m['n_species_truth']}  genus={m['n_genus_truth']}")
    for k in ["species_top1","species_top5","species_top10","species_mrr",
              "genus_top1","genus_top5","genus_top10","genus_mrr",
              "roc_auc"]:
        print(f"  {k:20s} = {m[k]}")
    print("=" * 72)

    write_outputs(results_dir, results, true_pairs)
    write_summary(results_dir, args.exp_id, m, args.model,
                  len(results), elapsed_h, n_err, n_cutoff)
    print(f"\n结果已保存至: {results_dir}")


if __name__ == "__main__":
    main()
