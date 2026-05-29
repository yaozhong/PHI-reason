#!/usr/bin/env python3
"""
build_phage_profiles_blastn.py
==============================
Extends phage genome text profiles with whole-genome BLASTN phylogenetic
neighbor context, producing the BLASTN-augmented variant for LLM-based
phage-host prediction.

This script takes the output of build_phage_profiles.py (no-BLASTN profiles)
as input and inserts a "Phylogenetic Cluster Context" section summarizing
genomically similar phages with experimentally confirmed host genera.

BLASTN is run against a reference set of phages with known hosts (e.g.,
the training/dev split of the evaluation dataset). Test-set phage accessions
must NOT be included in the reference set to avoid data leakage.

Phylogenetic section format:
  ## Phylogenetic Cluster Context
    (whole-genome BLASTN vs experimentally confirmed phage-host pairs)
  Genomically similar phages with known hosts from training/dev reference set:
    - ACC (identity=99.7%, coverage=13%): infects Lactobacillus paracasei
    - ...
  Host signal: STRONG (all 5/5 neighbors → Lactobacillus)
  Gram hint:   Gram-positive host strongly suggested

Host signal classification:
  STRONG   — all neighbors agree on the same genus
  MODERATE — ≥60% of neighbors (and ≥2) agree on one genus
  MIXED    — no clear consensus

Gram hint is inferred from the known Gram-positive / Gram-negative
status of the neighbor host genera (hard-coded lookup table).

BLASTN neighbors JSON format (--blastn-json):
  {
    "phage_id": [
      {"sid": "KC171647.1", "pident": 99.7, "qcovs": 13.0, "host": "Lactobacillus paracasei"},
      ...
    ],
    ...
  }

Usage:
    # Run BLASTN first (example using BLAST+ makeblastdb / blastn):
    makeblastdb -in train_dev_phages.fna -dbtype nucl -out train_dev_db
    blastn -query test_phages.fna -db train_dev_db \\
           -outfmt "6 qseqid sseqid pident qcovs length" \\
           -perc_identity 70 -qcov_hsp_perc 3 -max_hsps 1 \\
           -out blastn_hits.tsv

    # Build the neighbors JSON from BLASTN output + host metadata:
    # (see companion script: build_blastn_neighbors_json.py)

    # Build BLASTN-augmented profiles:
    python build_phage_profiles_blastn.py \\
        --base-profiles /path/to/phage_profiles_no_blastn \\
        --blastn-json   /path/to/blastn_neighbors.json \\
        --out-dir       /path/to/phage_profiles_blastn \\
        [--phage-list   phage_ids.json] \\
        [--dry-run]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import re
from collections import Counter
from pathlib import Path

# ── Neighbour phage accession masking ──────────────────────────────���──────────
# Replace NCBI accession IDs in BLASTN neighbour lines with deterministic hashes
# (e.g. JQ182730.1 → PHAGE_187861D0) to prevent LLM from querying the accession.
# Host names are kept intact.
_ACC_RE = re.compile(r'\b([A-Z]{1,2}_?\d{5,9}\.\d+)\b')

def _mask_accession(acc: str) -> str:
    return "PHAGE_" + hashlib.md5(acc.encode()).hexdigest()[:8].upper()

# ── Gram classification of host genera ────────────────────────────────────────
# Used to infer Gram type from BLASTN neighbor host genera.
# Extend as needed for your dataset.
GRAM_POS: frozenset[str] = frozenset({
    "listeria", "bacillus", "staphylococcus", "streptococcus",
    "lactococcus", "enterococcus", "mycolicibacterium", "mycobacterium",
    "streptomyces", "corynebacterium", "rhodococcus", "nocardia",
    "clostridium", "clostridioides", "lactobacillus", "leuconostoc",
    "bifidobacterium", "arthrobacter", "microbacterium", "propionibacterium",
    "brevibacillus", "paenibacillus", "haloarcula", "halorubrum",
    "haloferax", "natrialba", "gordonia", "tsukamurella",
    "lacticaseibacillus", "ligilactobacillus", "limosilactobacillus",
})

GRAM_NEG: frozenset[str] = frozenset({
    "escherichia", "salmonella", "klebsiella", "vibrio", "pseudomonas",
    "acinetobacter", "aeromonas", "cronobacter", "citrobacter", "enterobacter",
    "serratia", "shigella", "yersinia", "burkholderia", "stenotrophomonas",
    "xanthomonas", "xylella", "helicobacter", "campylobacter", "bacteroides",
    "flavobacterium", "cellulophaga", "synechococcus", "prochlorococcus",
    "caulobacter", "rhizobium", "agrobacterium", "sinorhizobium",
    "bradyrhizobium", "sulfitobacter", "achromobacter", "acidianus",
    "alteromonas", "pseudoalteromonas", "shewanella", "photobacterium",
    "erwinia", "pantoea", "pectobacterium", "dickeya", "hafnia",
    "providencia", "proteus", "morganella", "edwardsiella",
})


def _genus_of(host: str) -> str:
    """Extract genus (first word) from a host string like 'Lactobacillus casei'."""
    return host.replace("_", " ").split()[0].lower()


def _classify_signal(genera: list[str]) -> str:
    """Return 'STRONG', 'MODERATE', or 'MIXED' based on genus consensus."""
    g_count = Counter(genera)
    top_g, top_n = g_count.most_common(1)[0]
    n_total = len(genera)

    if top_n == n_total:
        return f"STRONG (all {n_total}/{n_total} neighbors → {top_g.capitalize()})"
    elif top_n >= max(2, n_total * 0.6):
        return f"MODERATE ({top_n}/{n_total} neighbors → {top_g.capitalize()})"
    else:
        top3 = ", ".join(
            f"{g.capitalize()}×{c}" for g, c in g_count.most_common(3)
        )
        return f"MIXED ({top3})"


def _gram_hint(genera: list[str]) -> str:
    """Return Gram type hint based on neighbor host genera."""
    n_pos = sum(1 for g in genera if g in GRAM_POS)
    n_neg = sum(1 for g in genera if g in GRAM_NEG)
    if n_pos > n_neg:
        return "Gram-positive host strongly suggested"
    if n_neg > n_pos:
        return "Gram-negative host strongly suggested"
    return "Gram type unclear from neighbor hosts"


def build_phylo_section(neighbors: list[dict]) -> str:
    """Build the '## Phylogenetic Cluster Context' Markdown block.

    Args:
        neighbors: list of dicts with keys 'sid', 'pident', 'qcovs', 'host'.
                   Empty list → 'no neighbors found' message.

    Returns:
        Multi-line string (with trailing newline).
    """
    header = (
        "## Phylogenetic Cluster Context"
        " (whole-genome BLASTN vs experimentally confirmed phage-host pairs)"
    )

    if not neighbors:
        return (
            f"{header}\n"
            "No close genomic neighbors found in train/dev reference set "
            "(BLASTN, 70% identity, 3% coverage threshold).\n"
        )

    genera = [_genus_of(nb["host"]) for nb in neighbors]

    lines = [
        header,
        "Genomically similar phages with known hosts from training/dev reference set:",
    ]
    for nb in neighbors:
        sid_masked = _ACC_RE.sub(lambda m: _mask_accession(m.group(1)), nb['sid'])
        lines.append(
            f"  - {sid_masked}"
            f" (identity={nb['pident']:.1f}%, coverage={nb['qcovs']:.0f}%):"
            f" infects {nb['host']}"
        )

    lines.append(f"\nHost signal: {_classify_signal(genera)}")
    lines.append(f"Gram hint:   {_gram_hint(genera)}")
    return "\n".join(lines) + "\n"


def build_phylo_section_no_host_label(neighbors: list[dict]) -> str:
    """Build BLASTN context block with host labels removed (similarity scores only).

    Retains masked neighbor IDs, nucleotide identity, and query coverage; omits
    host species/genus, Host signal, and Gram hint lines to test whether BLASTN
    context provides useful signal beyond explicit host-label retrieval.

    Args:
        neighbors: list of dicts with keys 'sid', 'pident', 'qcovs'.
                   Empty list → 'no neighbors found' message.

    Returns:
        Multi-line string (with trailing newline).
    """
    header = (
        "## Phylogenetic Cluster Context"
        " (whole-genome BLASTN vs experimentally confirmed phage-host pairs)"
    )

    if not neighbors:
        return (
            f"{header}\n"
            "No close genomic neighbors found in train/dev reference set "
            "(BLASTN, 70% identity, 3% coverage threshold).\n"
        )

    lines = [
        header,
        "Genomically similar phages from training/dev reference set:",
    ]
    for nb in neighbors:
        sid_masked = _ACC_RE.sub(lambda m: _mask_accession(m.group(1)), nb['sid'])
        lines.append(
            f"  - {sid_masked}"
            f" (identity={nb['pident']:.1f}%, coverage={nb['qcovs']:.0f}%)"
        )
    return "\n".join(lines) + "\n"


def add_phylo_to_profile(src_text: str, neighbors: list[dict], mode: str = "original") -> str:
    """Insert the phylogenetic context section into an existing profile.

    The section is inserted immediately after the header block (the first
    non-empty block of lines, ending at the first blank line).

    Args:
        src_text:  Full text of the base profile.
        neighbors: BLASTN neighbor list (may be empty).
        mode: 'original' (include host labels) or 'no_host_label' (scores only).
              For 'scrambled', pass pre-scrambled neighbors and use 'original'.

    Returns:
        Modified profile text with phylogenetic section inserted.
    """
    if mode == "no_host_label":
        phylo_sec = build_phylo_section_no_host_label(neighbors)
    else:
        phylo_sec = build_phylo_section(neighbors)
    text_lines = src_text.rstrip("\n").split("\n")

    # Find the end of the header block (first blank line after line 0)
    insert_at = 1
    while insert_at < len(text_lines) and text_lines[insert_at].strip():
        insert_at += 1

    new_lines = (
        text_lines[:insert_at]
        + ["", phylo_sec.rstrip("\n")]
        + text_lines[insert_at:]
    )
    return "\n".join(new_lines) + "\n"


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build BLASTN-augmented phage profiles for LLM host prediction. "
            "Takes no-BLASTN profiles (from build_phage_profiles.py) and adds "
            "a whole-genome phylogenetic neighbor section."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--base-profiles", required=True,
                        help="Input profile directory (output of build_phage_profiles.py)")
    parser.add_argument("--blastn-json",   required=True,
                        help="JSON: {phage_id: [{sid, pident, qcovs, host}, ...]}")
    parser.add_argument("--out-dir",       required=True,
                        help="Output directory for BLASTN-augmented profiles")
    parser.add_argument("--phage-list",    default=None,
                        help="JSON list of phage IDs to process (default: all in base-profiles)")
    parser.add_argument("--dry-run",       action="store_true",
                        help="Print preview of first 2 profiles without writing files")
    parser.add_argument("--blast_context_mode", default="original",
                        choices=["original", "no_host_label", "scrambled"],
                        help="'original': include host labels (default); "
                             "'no_host_label': strip host species/genus/Gram from context output; "
                             "'scrambled': globally permute host labels across all phages (seed=42)")
    args = parser.parse_args()

    base_dir  = Path(args.base_profiles)
    out_dir   = Path(args.out_dir)
    blastn_nb = json.loads(Path(args.blastn_json).read_text())

    if args.phage_list:
        phage_ids: list[str] = json.loads(Path(args.phage_list).read_text())
    else:
        phage_ids = [p.stem for p in base_dir.glob("*.md")]

    # ── Scrambled mode: globally permute host labels across all phages ─────────
    if args.blast_context_mode == "scrambled":
        all_slots: list[tuple[str, int]] = [
            (pid, i)
            for pid in phage_ids
            for i in range(len(blastn_nb.get(pid, [])))
        ]
        hosts_pool = [blastn_nb[pid][i]["host"] for pid, i in all_slots]
        random.Random(42).shuffle(hosts_pool)
        blastn_nb = {pid: [dict(nb) for nb in blastn_nb.get(pid, [])] for pid in phage_ids}
        for (pid, i), host in zip(all_slots, hosts_pool):
            blastn_nb[pid][i]["host"] = host

    n_with_nb = sum(1 for pid in phage_ids if pid in blastn_nb)
    print(f"Base profiles : {base_dir}")
    print(f"Output dir    : {out_dir}")
    print(f"Phages        : {len(phage_ids)}")
    print(f"With neighbors: {n_with_nb} / {len(phage_ids)}")
    print(f"Context mode  : {args.blast_context_mode}")
    print(f"Dry run       : {args.dry_run}")
    print()

    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)

    built = skipped_missing = preview_shown = 0

    for phage_id in phage_ids:
        src_path = base_dir / f"{phage_id}.md"
        if not src_path.exists():
            print(f"  MISSING base profile: {phage_id}")
            skipped_missing += 1
            continue

        neighbors = blastn_nb.get(phage_id, [])
        new_text  = add_phylo_to_profile(src_path.read_text(), neighbors, mode=args.blast_context_mode)

        if args.dry_run:
            if preview_shown < 2:
                print(f"  [DRY] {phage_id} ({len(neighbors)} neighbors)")
                print(new_text[:600])
                print("  ...")
                preview_shown += 1
        else:
            (out_dir / f"{phage_id}.md").write_text(new_text)

        built += 1

    print()
    if not args.dry_run:
        print(f"Done.")
        print(f"  Profiles written: {built}")
        print(f"  Missing base:     {skipped_missing}")

        # Signal distribution summary
        strong = moderate = mixed_ = no_nb = 0
        for pid in phage_ids:
            nbs = blastn_nb.get(pid, [])
            if not nbs:
                no_nb += 1
                continue
            genera = [_genus_of(nb["host"]) for nb in nbs]
            top_n  = Counter(genera).most_common(1)[0][1]
            n_tot  = len(genera)
            if top_n == n_tot:
                strong += 1
            elif top_n >= max(2, n_tot * 0.6):
                moderate += 1
            else:
                mixed_ += 1

        print(f"\nPhylogenetic signal distribution (n={len(phage_ids)}):")
        print(f"  STRONG   (all neighbors agree)      : {strong}")
        print(f"  MODERATE (≥60% neighbors agree)     : {moderate}")
        print(f"  MIXED    (no clear consensus)        : {mixed_}")
        print(f"  NO NEIGHBORS                         : {no_nb}")
        print(f"\nOutput: {out_dir}")
    else:
        print(f"[DRY RUN] Would write {built} profiles.")


if __name__ == "__main__":
    main()
