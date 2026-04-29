#!/usr/bin/env python3
"""
05_inject_context.py
========================
Inject BLASTP phylogenetic context blocks into existing phage profiles
(rbp+blastn versions), producing a new rbp_blastn_blastp version.

Source profiles already contain:
  ## Phylogenetic Cluster Context (whole-genome BLASTN ...)
  ...
  ## Tail & Host Adsorption [KEY ...]
  ...

This script inserts a SECOND context block (BLASTP) immediately after
the BLASTN block:

  ## Phylogenetic Cluster Context (BLASTN ...)
  ...
                                  ← inserted here
  ## Phylogenetic Cluster Context (proteome similarity via BLASTP, CHERRY-1306 reference)
  ...

  ## Tail & Host Adsorption [KEY ...]
  ...

Usage:
  python 05_inject_context.py \
    --src_dir    ws/Cherry/textGeneProfile_v2cm_R1_rbp_blastn_masked \
    --context_json experiments/blastp_cherry1306/outputs/blastp_context_cherry634.json \
    --out_dir    ws/Cherry/textGeneProfile_v2cm_R1_rbp_blastn_blastp

  python 05_inject_context.py \
    --src_dir    ws/VHDB/phage_profiles_rbp_blastn_noleak \
    --context_json experiments/blastp_cherry1306/outputs/blastp_context_vhdb.json \
    --out_dir    ws/VHDB/phage_profiles_rbp_blastn_blastp_noleak
"""
from __future__ import annotations
import argparse, json, re, sys
from pathlib import Path

# Header regex matching BLASTN section start
_BLASTN_HDR = re.compile(
    r"^##\s+Phylogenetic Cluster Context\s*\((?:whole-genome\s+)?BLASTN",
    re.MULTILINE | re.IGNORECASE,
)
# Header regex matching ANY ## section start
_SEC_HDR = re.compile(r"^##\s", re.MULTILINE)

NO_BLASTP_BLOCK = (
    "## Phylogenetic Cluster Context (proteome similarity via BLASTP, CHERRY-1306 reference)\n"
    "No proteomically similar train-set phage found in CHERRY-1306 reference "
    "(no significant BLASTP hits).\n"
    "\n"
    "Host signal: No close neighbors\n"
    "Gram hint:   Gram type unclear from neighbor hosts\n"
)

NO_BLASTP_BLOCK_NO_HOST_LABEL = (
    "## Phylogenetic Cluster Context (proteome similarity via BLASTP, CHERRY-1306 reference)\n"
    "No proteomically similar train-set phage found in CHERRY-1306 reference "
    "(no significant BLASTP hits).\n"
)


def find_blastn_block_end(text: str, blastn_start: int) -> int:
    """Return the position just after the BLASTN section (start of next ## or EOF)."""
    rest = text[blastn_start:]
    # skip the first ## (the BLASTN header itself)
    m = _SEC_HDR.search(rest, 1)
    if m:
        return blastn_start + m.start()
    return len(text)


def inject(profile_text: str, blastp_block: str) -> str:
    """Insert BLASTP context block right after the BLASTN block."""
    m = _BLASTN_HDR.search(profile_text)
    if m is None:
        # No BLASTN section at all — insert before the first ## section
        m2 = _SEC_HDR.search(profile_text)
        insert_at = m2.start() if m2 else len(profile_text)
        head = profile_text[:insert_at].rstrip("\n") + "\n\n"
        tail = profile_text[insert_at:]
        return head + blastp_block.rstrip("\n") + "\n\n" + tail

    end = find_blastn_block_end(profile_text, m.start())
    before = profile_text[:end].rstrip("\n") + "\n\n"
    after  = profile_text[end:]
    return before + blastp_block.rstrip("\n") + "\n\n" + after


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--src_dir",      type=Path, required=True)
    ap.add_argument("--context_json", type=Path, required=True)
    ap.add_argument("--out_dir",      type=Path, required=True)
    ap.add_argument("--no_placeholder", action="store_true",
                    help="Skip files without a BLASTP hit (do not write placeholder)")
    ap.add_argument("--blast_context_mode", default="original",
                    choices=["original", "no_host_label"],
                    help="'original' (default): preserve host labels in placeholder block; "
                         "'no_host_label': use host-label-free placeholder for zero-hit phages")
    args = ap.parse_args()

    if not args.src_dir.exists():
        print(f"[ERR] src_dir not found: {args.src_dir}"); return 1
    if not args.context_json.exists():
        print(f"[ERR] context_json not found: {args.context_json}"); return 1

    context: dict[str, dict] = json.loads(args.context_json.read_text())
    print(f"[1] context_json: {len(context)} entries  ({args.context_json})")

    fallback_block = (
        NO_BLASTP_BLOCK_NO_HOST_LABEL
        if args.blast_context_mode == "no_host_label"
        else NO_BLASTP_BLOCK
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    md_files = sorted(args.src_dir.glob("*.md"))
    n_total = len(md_files)
    print(f"[2] src profiles: {n_total}  →  {args.out_dir}")

    n_injected = n_placeholder = n_skipped = 0
    for src in md_files:
        acc = src.stem
        entry = context.get(acc)
        if entry:
            blastp_block = entry["context_block"].rstrip("\n") + "\n"
            n_injected += 1
        elif args.no_placeholder:
            n_skipped += 1
            import shutil; shutil.copy2(src, args.out_dir / src.name)
            continue
        else:
            blastp_block = fallback_block
            n_placeholder += 1

        merged = inject(src.read_text(), blastp_block)
        (args.out_dir / src.name).write_text(merged)

    print(f"[3] done")
    print(f"    injected BLASTP context : {n_injected}")
    print(f"    placeholder (no hit)    : {n_placeholder}")
    print(f"    copied (no-placeholder) : {n_skipped}")
    print(f"    total written           : {n_injected + n_placeholder + n_skipped} / {n_total}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
