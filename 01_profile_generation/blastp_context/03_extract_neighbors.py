#!/usr/bin/env python3
"""
03_neighbors_top5.py — Aggregate diamond blastp TSV → top-5 neighbor JSON.

Input:
  --tsv       runs/<TESTSET>_vs_cherry1306.tsv
              outfmt 6 qseqid sseqid pident length evalue bitscore
              IDs are formatted as ACC__geneN (both query and subject).
  --pair_csv  db/cherry1306_phage_host_pair.csv   (acc,host — DB-side only)
  --out_json  outputs/neighbors_top5_<TESTSET>.json

Aggregation:
  For each (query_phage, subject_phage) pair in the DB:
    summed_bitscore = Σ bitscore of all significant protein hits
  Rank subject phages by summed bitscore; keep top-K (default 5).
  Annotate each neighbor with its known host label (from DB-side pair csv).

No LOO filter: query is an independent test set, DB is cherry1306 train.
If a test-phage acc happens to coincide with a DB acc (should not, if PLAN is
respected), the script warns and still keeps the hit — leakage is a pipeline
problem, not something to hide in aggregation.

Output schema:
  {
    "<query_acc>": [
      {"acc": "<subject_acc>", "score": <float>, "host": "<host_or_unknown>"},
      ...  (up to top_k)
    ],
    ...
  }
"""
from __future__ import annotations
import argparse, csv, json, sys
from collections import defaultdict
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv",      type=Path, required=True)
    ap.add_argument("--pair_csv", type=Path, required=True)
    ap.add_argument("--out_json", type=Path, required=True)
    ap.add_argument("--top_k",    type=int, default=5)
    args = ap.parse_args()

    # 1. Load DB-side host labels (cherry1306 only)
    host_of: dict[str, str] = {}
    with args.pair_csv.open() as f:
        for line in f:
            line = line.strip()
            if not line or "," not in line:
                continue
            acc, host = line.split(",", 1)
            host_of[acc] = host
    print(f"[1] loaded {len(host_of)} DB-side host labels from {args.pair_csv}")

    # 2. Parse diamond TSV, aggregate per (q_phage, s_phage)
    phage_sim: dict[tuple, float] = defaultdict(float)
    n_rows = 0
    q_accs: set[str] = set()
    leak_accs: set[str] = set()
    with args.tsv.open() as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 6:
                continue
            q_prot, s_prot = row[0], row[1]
            bitscore = float(row[5])
            q_phage = q_prot.split("__")[0]
            s_phage = s_prot.split("__")[0]
            q_accs.add(q_phage)
            if q_phage == s_phage:
                leak_accs.add(q_phage)
                continue  # self-acc must not be treated as neighbor
            phage_sim[(q_phage, s_phage)] += bitscore
            n_rows += 1
    print(f"[2] parsed {n_rows:,} protein hits across "
          f"{len(phage_sim):,} (query,subject) phage pairs")
    print(f"    queries with ≥1 hit: {len(q_accs)}")
    if leak_accs:
        print(f"    [WARN] {len(leak_accs)} query accs coincide with DB accs "
              f"(expected 0 under PLAN): first 5 → {sorted(leak_accs)[:5]}")

    # 3. Rank subject phages per query, keep top_k
    neighbors_by_query: dict[str, list] = defaultdict(list)
    for (q, s), score in phage_sim.items():
        neighbors_by_query[q].append({
            "acc":   s,
            "score": round(score, 1),
            "host":  host_of.get(s, "unknown"),
        })
    result: dict[str, list] = {}
    for q, hits in neighbors_by_query.items():
        top = sorted(hits, key=lambda x: -x["score"])[:args.top_k]
        result[q] = top

    # 4. Ensure every query acc that appeared anywhere gets a key (empty list if no neighbors)
    for q in q_accs:
        result.setdefault(q, [])

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    with args.out_json.open("w") as f:
        json.dump(result, f, indent=2)

    covered = sum(1 for v in result.values() if v)
    print(f"[3] wrote {args.out_json}")
    print(f"    queries         : {len(result)}")
    print(f"    with neighbors  : {covered}  ({covered/max(1,len(result)):.1%})")

    unknown_ratio = sum(
        1 for hits in result.values() for h in hits if h["host"] == "unknown"
    ) / max(1, sum(len(v) for v in result.values()))
    print(f"    unknown-host hits ratio: {unknown_ratio:.1%}")

    # 5. Sanity print (first query with neighbors)
    first = next((q for q in sorted(result) if result[q]), None)
    if first:
        print(f"\nSample — {first}:")
        for n in result[first]:
            print(f"  {n['acc']}  score={n['score']}  host={n['host']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
