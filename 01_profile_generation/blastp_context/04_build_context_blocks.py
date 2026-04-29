#!/usr/bin/env python3
"""
04_context_blocks.py — Convert top-5 neighbor JSON → LLM-ready context blocks.

Input:
  --neighbors_json  outputs/neighbors_top5_<TESTSET>.json   (from 03_neighbors_top5.py)
  --pair_csv        db/cherry1306_phage_host_pair.csv       (for gram lookup)
  --out_json        outputs/blastp_context_<TESTSET>.json

Output schema (one entry per query acc that had ≥1 neighbor):
  {
    "<query_acc>": {
      "context_block": "## Phylogenetic Cluster Context (...)\\n...",
      "gram_hint":     "Gram-positive" | "Gram-negative" | ""
    },
    ...
  }

Queries with 0 neighbors are OMITTED (caller should fall back to BLASTN / RBP
profile). This matches the behavior of the prior blastp_context_v6 pipeline.

Format mirrors /mnt/hdd/remote_phi_reason/cherry1940/blast/blastp_context_v6.json
but the header is updated to reflect the clean cherry1306 reference scope:

  ## Phylogenetic Cluster Context (proteome similarity via BLASTP, CHERRY-1306 reference)
  Top proteomically similar train-set phages (BLASTP, query not in DB):
    - <acc> (score=<float>): infects <host>
    ...
  Host signal: STRONG | MODERATE | MIXED | No close neighbors
  Gram hint:   Gram-(positive|negative) host strongly suggested | Gram type unclear ...

Host signal rule:
  genera = [host.split("_")[0] for host in neighbors if host != "unknown"]
  top_frac = count(top_genus) / len(genera)
    1.0   → STRONG
    ≥0.6  → MODERATE
    else  → MIXED

Gram hint rule:
  vote over neighbor hosts by curated gram dict; if ≥60% agree → emit
  "Gram-positive" / "Gram-negative"; else "".
"""
from __future__ import annotations
import argparse, hashlib, json, random, re, sys
from collections import Counter
from pathlib import Path

# Mask NCBI accession IDs (e.g. AB002632.1) → PHAGE_<MD5[:8]> so LLM cannot
# look up the training-set phage by accession.  Host names are kept intact.
_ACC_RE = re.compile(r'\b([A-Z]{1,2}_?\d{5,9}\.\d+)\b')

def _mask_acc(acc: str) -> str:
    return "PHAGE_" + hashlib.md5(acc.encode()).hexdigest()[:8].upper()


GRAM_NEG = {
    "escherichia", "salmonella", "klebsiella", "pseudomonas", "vibrio",
    "acinetobacter", "shigella", "yersinia", "serratia", "enterobacter",
    "citrobacter", "proteus", "morganella", "providencia", "hafnia",
    "pectobacterium", "dickeya", "erwinia", "xanthomonas", "caulobacter",
    "agrobacterium", "brucella", "rhizobium", "burkholderia", "ralstonia",
    "delftia", "comamonas", "stenotrophomonas", "acidovorax", "campylobacter",
    "helicobacter", "neisseria", "moraxella", "haemophilus", "pasteurella",
    "mannheimia", "actinobacillus", "aggregatibacter", "myxococcus",
    "bdellovibrio", "cyanobacteria", "synechococcus", "prochlorococcus",
    "cronobacter", "sodalis", "pantoea", "rahnella", "cedecea",
    "aliivibrio", "alteromonas", "aeromonas", "sinorhizobium",
    "mesorhizobium", "azospirillum",
}
GRAM_POS = {
    "bacillus", "staphylococcus", "streptococcus", "lactococcus",
    "lactobacillus", "enterococcus", "listeria", "clostridium",
    "mycobacterium", "mycolicibacterium", "corynebacterium", "arthrobacter",
    "brevibacterium", "micrococcus", "paenibacillus", "geobacillus",
    "brevibacillus", "streptomyces", "propionibacterium", "bifidobacterium",
    "cutibacterium", "actinomyces", "rhodococcus", "nocardia",
}


def build_host_gram(pair_csv: Path) -> dict[str, str]:
    """host_string → 'Gram-negative' | 'Gram-positive'."""
    host_gram: dict[str, str] = {}
    for line in pair_csv.read_text().splitlines():
        line = line.strip()
        if not line or "," not in line:
            continue
        _, host = line.split(",", 1)
        genus = host.split("_")[0].lower()
        if genus in GRAM_NEG:
            host_gram[host] = "Gram-negative"
        elif genus in GRAM_POS:
            host_gram[host] = "Gram-positive"
    return host_gram


def signal_quality(neighbors: list[dict]) -> str:
    if not neighbors:
        return "No close neighbors"
    genera = [n["host"].split("_")[0]
              for n in neighbors if n.get("host", "unknown") != "unknown"]
    if not genera:
        return "No close neighbors"
    top_genus, top_count = Counter(genera).most_common(1)[0]
    frac = top_count / len(genera)
    if frac == 1.0:
        return f"STRONG (all {top_count}/{len(genera)} neighbors → {top_genus})"
    if frac >= 0.6:
        return f"MODERATE ({top_count}/{len(genera)} neighbors agree → {top_genus})"
    return f"MIXED (top genus {top_genus}: only {top_count}/{len(genera)} neighbors agree)"


def gram_hint(neighbors: list[dict], host_gram: dict[str, str]) -> str:
    grams = [host_gram.get(n["host"], "") for n in neighbors]
    gc = Counter(g for g in grams if g)
    if not gc:
        return ""
    top_g, cnt = gc.most_common(1)[0]
    if cnt / len(neighbors) >= 0.6:
        return "Gram-positive" if "positive" in top_g else "Gram-negative"
    return ""


def build_block(neighbors: list[dict], host_gram: dict[str, str]) -> dict | None:
    if not neighbors:
        return None
    lines = [
        "## Phylogenetic Cluster Context (proteome similarity via BLASTP, CHERRY-1306 reference)",
        "Top proteomically similar train-set phages (BLASTP, query not in DB):",
    ]
    for n in neighbors:
        masked = _ACC_RE.sub(lambda m: _mask_acc(m.group(1)), n["acc"])
        lines.append(f"  - {masked} (score={n['score']}): infects {n['host']}")
    lines.append("")
    lines.append(f"Host signal: {signal_quality(neighbors)}")
    g = gram_hint(neighbors, host_gram)
    if g:
        lines.append(f"Gram hint:   {g} host strongly suggested")
    else:
        lines.append("Gram hint:   Gram type unclear from neighbor hosts")
    lines.append("")
    return {"context_block": "\n".join(lines), "gram_hint": g}


def build_block_no_host_label(neighbors: list[dict]) -> dict | None:
    """Build BLASTP context block with host labels removed (similarity scores only).

    Retains masked accession IDs and BLASTP scores; omits host species/genus,
    Host signal, and Gram hint lines. Used for the host-label-free ablation.
    """
    if not neighbors:
        return None
    lines = [
        "## Phylogenetic Cluster Context (proteome similarity via BLASTP, CHERRY-1306 reference)",
        "Top proteomically similar train-set phages (BLASTP, query not in DB):",
    ]
    for n in neighbors:
        masked = _ACC_RE.sub(lambda m: _mask_acc(m.group(1)), n["acc"])
        lines.append(f"  - {masked} (score={n['score']})")
    lines.append("")
    return {"context_block": "\n".join(lines), "gram_hint": ""}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--neighbors_json", type=Path, required=True)
    ap.add_argument("--pair_csv",       type=Path, required=True)
    ap.add_argument("--out_json",       type=Path, required=True)
    ap.add_argument("--blast_context_mode", default="original",
                    choices=["original", "no_host_label", "scrambled"],
                    help="'original': include host labels (default); "
                         "'no_host_label': strip host/Gram from output; "
                         "'scrambled': globally permute host labels across all queries (seed=42)")
    args = ap.parse_args()

    host_gram = build_host_gram(args.pair_csv)
    print(f"[1] host_gram: {len(host_gram)} hosts classified "
          f"(Gram-neg={sum(1 for v in host_gram.values() if v=='Gram-negative')}, "
          f"Gram-pos={sum(1 for v in host_gram.values() if v=='Gram-positive')})")

    neighbors_raw: dict[str, list[dict]] = json.loads(args.neighbors_json.read_text())
    print(f"[2] neighbors_json: {len(neighbors_raw)} queries  mode={args.blast_context_mode}")

    # ── Scrambled mode: globally permute host labels across all queries ─────────
    if args.blast_context_mode == "scrambled":
        all_slots = [(acc, i) for acc, hits in neighbors_raw.items() for i in range(len(hits))]
        hosts_pool = [neighbors_raw[acc][i]["host"] for acc, i in all_slots]
        random.Random(42).shuffle(hosts_pool)
        neighbors_raw = {acc: [dict(h) for h in hits] for acc, hits in neighbors_raw.items()}
        for (acc, i), host in zip(all_slots, hosts_pool):
            neighbors_raw[acc][i]["host"] = host

    neighbors = neighbors_raw

    out: dict[str, dict] = {}
    n_strong = n_moderate = n_mixed = n_none = 0
    for acc, hits in neighbors.items():
        if args.blast_context_mode == "no_host_label":
            block = build_block_no_host_label(hits)
        else:
            block = build_block(hits, host_gram)
        if block is None:
            n_none += 1
            continue
        out[acc] = block
        sq = block["context_block"]
        if "STRONG" in sq:     n_strong   += 1
        elif "MODERATE" in sq: n_moderate += 1
        elif "MIXED" in sq:    n_mixed    += 1

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(out, indent=2))

    total_with = len(out)
    print(f"[3] wrote {args.out_json}")
    print(f"    queries with context : {total_with}")
    print(f"    queries w/o neighbors: {n_none}")
    print(f"    signal STRONG   : {n_strong}   ({n_strong/max(1,total_with):.1%})")
    print(f"    signal MODERATE : {n_moderate} ({n_moderate/max(1,total_with):.1%})")
    print(f"    signal MIXED    : {n_mixed}    ({n_mixed/max(1,total_with):.1%})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
