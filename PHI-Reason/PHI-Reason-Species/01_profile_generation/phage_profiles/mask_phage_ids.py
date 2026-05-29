#!/usr/bin/env python3
"""
mask_phage_ids.py — Mask test phage accession IDs in profile files to prevent
LLM recognition.  Replaces header accession with PHAGE_<MD5> pseudonym and
scrubs stray NCBI accession patterns in the BLASTN context block.
"""
from __future__ import annotations
import argparse, hashlib, re
from pathlib import Path

_ACC_RE = re.compile(r"\b([A-Z]{1,2}_?\d{5,9}(?:\.\d+)?)\b")
_HEADER_RE = re.compile(r"^(# Phage Genome Profile:\s+)(\S+)")


def phage_hash(pid: str) -> str:
    return "PHAGE_" + hashlib.md5(pid.encode()).hexdigest()[:8].upper()


def mask_profile(text: str) -> tuple[str, bool]:
    """Mask header accession and stray accessions in BLASTN context block."""
    modified = False
    m = _HEADER_RE.match(text)
    if m:
        text = m.group(1) + phage_hash(m.group(2)) + text[m.end():]
        modified = True
    cs = text.find("## Phylogenetic Cluster Context")
    if cs != -1:
        ns = text.find("\n## ", cs + 1)
        ns = ns if ns != -1 else len(text)
        block = text[cs:ns]
        masked = _ACC_RE.sub(lambda m: phage_hash(m.group(1)), block)
        if masked != block:
            text = text[:cs] + masked + text[ns:]
            modified = True
    return text, modified


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Mask test phage accession IDs in profile files.",
    )
    parser.add_argument("--src-dir", required=True, type=Path,
                        help="Input profiles directory")
    parser.add_argument("--out-dir", required=True, type=Path,
                        help="Output directory for masked profiles")
    parser.add_argument("--phage-list", type=Path, default=None,
                        help="Text file of phage IDs to process (one per line; default: all .md)")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.phage_list:
        ids = [l.strip() for l in args.phage_list.read_text().splitlines() if l.strip()]
        md_files = [args.src_dir / f"{pid}.md" for pid in ids]
    else:
        md_files = sorted(args.src_dir.glob("*.md"))

    n_masked = n_unchanged = n_missing = 0
    for md in md_files:
        if not md.exists():
            print(f"  MISSING: {md.name}")
            n_missing += 1
            continue
        masked_text, was_modified = mask_profile(md.read_text())
        (args.out_dir / md.name).write_text(masked_text)
        if was_modified:
            n_masked += 1
        else:
            n_unchanged += 1

    total = n_masked + n_unchanged
    print(f"Done: {total} profile(s) processed, {n_missing} missing")
    print(f"  Masked:    {n_masked}")
    print(f"  Unchanged: {n_unchanged}")
    print(f"  Output:    {args.out_dir}")


if __name__ == "__main__":
    main()
