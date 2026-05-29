#!/usr/bin/env python3
"""
build_blastn_neighbors_json.py
==============================
Convert whole-genome BLASTN output (outfmt 6 + qcovs) and a phage-host
metadata CSV into a JSON file consumed by build_phage_profiles.py.

Output schema
-------------
{
  "<query_phage_acc>": [
    {"sid": "<subject_phage_acc>", "pident": 85.3, "qcovs": 72, "host": "Escherichia_coli"},
    ...
  ]
}

BLASTN TSV columns (outfmt 6 + qcovs):
  qseqid sseqid pident length mismatch gapopen qstart qend sstart send
  evalue bitscore qcovs

Host CSV (no header):  accession,host_species
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build BLASTN neighbors JSON from TSV + host CSV.",
    )
    p.add_argument("--blastn-tsv", required=True, type=Path,
                    help="BLASTN output TSV (outfmt 6 + qcovs, 13 columns)")
    p.add_argument("--host-csv", required=True, type=Path,
                    help="Phage-to-host CSV: accession,host_species (no header)")
    p.add_argument("--out-json", required=True, type=Path,
                    help="Output JSON path")
    p.add_argument("--top-k", type=int, default=5,
                    help="Number of top neighbors per query (default: 5)")
    p.add_argument("--min-identity", type=float, default=70.0,
                    help="Minimum percent identity threshold (default: 70.0)")
    p.add_argument("--min-qcov", type=int, default=50,
                    help="Minimum query coverage threshold (default: 50)")
    p.add_argument("--exclude-list", type=Path, default=None,
                    help="Text file of accessions to exclude (one per line)")
    return p.parse_args()


def load_host_map(csv_path: Path) -> dict[str, str]:
    """Load accession -> host_species mapping from a headerless CSV."""
    host_map: dict[str, str] = {}
    with open(csv_path, newline="") as fh:
        reader = csv.reader(fh)
        for row in reader:
            if len(row) >= 2:
                host_map[row[0].strip()] = row[1].strip()
    return host_map


def load_exclude_set(path: Path | None) -> set[str]:
    """Load a set of accessions to exclude (one per line)."""
    if path is None:
        return set()
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def main() -> None:
    args = parse_args()

    host_map = load_host_map(args.host_csv)
    exclude = load_exclude_set(args.exclude_list)

    # ── Parse BLASTN TSV ─────────────────────────────────────────────
    # Collect all passing hits per query, keyed by (query, subject).
    # Keep the best bitscore per query-subject pair (in case of dupes).
    hits: dict[str, list[dict]] = defaultdict(list)

    with open(args.blastn_tsv) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 13:
                print(f"WARNING: skipping line {lineno} ({len(cols)} columns)", file=sys.stderr)
                continue

            qseqid = cols[0]
            sseqid = cols[1]
            pident = float(cols[2])
            bitscore = float(cols[11])
            qcovs = float(cols[12])

            # Filters
            if qseqid == sseqid:
                continue
            if sseqid in exclude:
                continue
            if pident < args.min_identity:
                continue
            if qcovs < args.min_qcov:
                continue
            if sseqid not in host_map:
                continue

            hits[qseqid].append({
                "sid": sseqid,
                "pident": round(pident, 2),
                "qcovs": int(qcovs),
                "host": host_map[sseqid],
                "_bitscore": bitscore,
            })

    # ── Rank & select top-k per query ────────────────────────────────
    result: dict[str, list[dict]] = {}
    for qid, hit_list in hits.items():
        # Sort by bitscore descending; keep top-k unique subjects
        hit_list.sort(key=lambda h: h["_bitscore"], reverse=True)
        seen: set[str] = set()
        top: list[dict] = []
        for h in hit_list:
            if h["sid"] in seen:
                continue
            seen.add(h["sid"])
            top.append({k: v for k, v in h.items() if k != "_bitscore"})
            if len(top) >= args.top_k:
                break
        result[qid] = top

    # ── Write output ─────────────────────────────────────────────────
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_json, "w") as fh:
        json.dump(result, fh, indent=2)

    # ── Summary ──────────────────────────────────────────────────────
    n_queries = len(result)
    n_with_neighbors = sum(1 for v in result.values() if v)
    neighbor_counts = [len(v) for v in result.values()]
    avg_neighbors = sum(neighbor_counts) / max(n_queries, 1)

    print(f"Queries with ≥1 neighbor : {n_with_neighbors} / {n_queries}")
    print(f"Avg neighbors per query  : {avg_neighbors:.1f}")
    print(f"Host CSV entries loaded  : {len(host_map)}")
    print(f"Excluded accessions      : {len(exclude)}")
    print(f"Output written to        : {args.out_json}")


if __name__ == "__main__":
    main()
