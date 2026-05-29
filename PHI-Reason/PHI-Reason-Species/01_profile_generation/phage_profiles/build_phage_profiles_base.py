#!/usr/bin/env python3
"""
Stage 0 — Build per-phage base profile (.md) using PHROG annotations only.

This script generates a structured Markdown genome profile for each phage,
summarising its ORFs by functional category (PHROG) with biophysical
properties computed from the amino-acid sequence.

Inputs
------
--phage-dir   Directory containing per-phage subdirectories.
              Each subdirectory is named by phage ID and must contain:
                - proteins.faa        (Prodigal amino-acid FASTA)
                - phrog_hits.tsv      (optional; DIAMOND vs PHROG:
                  columns gene_id, phrog_id [phrog_NNN], pident, qcov, bitscore)

--phrog-lookup  TSV file with columns: phrog_id (int), annotation, category.
                Used to map PHROG hit IDs to human-readable annotations and
                functional categories.

--out-dir     Output directory. One file per phage: {out_dir}/{phage_id}.md

--phage-list  (optional) Text file listing phage IDs to process (one per line).
              If omitted, every subdirectory of --phage-dir is processed.

Output profile section order
----------------------------
  1. Tail & Host Adsorption          (PHROG category 'tail')
  2. Lysis Cassette                  (PHROG category 'lysis')
  3. Head & DNA Packaging            (PHROG categories 'head and packaging', 'connector')
  4. DNA Replication & Metabolism     (PHROG category 'DNA, RNA and nucleotide metabolism')
  5. Gene Expression & Regulation    (PHROG category 'transcription regulation')
  6. Lysogeny / Integration          (PHROG category 'integration and excision')
  7. Moron / Host-Derived Factors    (PHROG category 'moron, auxiliary metabolic gene and host takeover')
  8. Other Annotated Genes           (PHROG categories 'other', 'unknown function')
  9. Hypothetical Proteins: N (no annotation)

Gene line format (no eggNOG, no phage-family):
  - gene{N} [LOW]: unknown [PHROG: <annot>]. (Naa, pI=X.X, X.XkDa[, hydrophobic(GRAVY=...)][, ~NTM-helix])

Dependencies
------------
  - BioPython (Bio.SeqUtils.ProtParam.ProteinAnalysis)
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ---------------------------------------------------------------------------
# Section ordering and PHROG-category mapping
# ---------------------------------------------------------------------------
SECTION_ORDER = [
    ("tail", "## Tail & Host Adsorption [KEY for host specificity]"),
    ("lysis", "## Lysis Cassette [Gram-type indicator]"),
    ("head_pack", "## Head & DNA Packaging"),
    ("dna_meta", "## DNA Replication & Metabolism"),
    ("transcription", "## Gene Expression & Regulation"),
    ("integration", "## Lysogeny / Integration"),
    ("moron", "## Moron / Host-Derived Factors"),
    ("other", "## Other Annotated Genes"),
]
CATEGORY_TO_SECTION = {
    "tail": "tail",
    "lysis": "lysis",
    "head and packaging": "head_pack",
    "connector": "head_pack",
    "DNA, RNA and nucleotide metabolism": "dna_meta",
    "transcription regulation": "transcription",
    "integration and excision": "integration",
    "moron, auxiliary metabolic gene and host takeover": "moron",
    "other": "other",
    "unknown function": "other",
}

# ---------------------------------------------------------------------------
# Biophysics constants
# ---------------------------------------------------------------------------
# Kyte-Doolittle hydropathy scale
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
TM_WINDOW = 19
TM_MIN_MEAN = 1.6
HYDROPHOBIC_GRAVY = 0.4
HYDROPHILIC_GRAVY = -0.8

PHROG_ID_RE = re.compile(r"phrog_(\d+)")


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def load_phrog_lookup(lookup_path: Path) -> dict[int, tuple[str, str]]:
    """Load phrog_lookup.tsv -> {phrog_id: (annotation, category)}."""
    m: dict[int, tuple[str, str]] = {}
    with lookup_path.open() as fh:
        next(fh)  # skip header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                pid = int(f[0])
            except ValueError:
                continue
            m[pid] = (f[1], f[2])
    return m


def parse_faa(path: Path) -> list[tuple[str, str]]:
    """Parse a Prodigal FASTA and return [(gene_id, sequence), ...] in order."""
    out: list[tuple[str, str]] = []
    cur_seq: list[str] = []
    idx = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_seq:
                    out.append((f"gene{idx}", "".join(cur_seq).rstrip("*")))
                idx += 1
                cur_seq = []
            else:
                cur_seq.append(line.strip())
    if cur_seq:
        out.append((f"gene{idx}", "".join(cur_seq).rstrip("*")))
    return out


def parse_phrog_hits(path: Path) -> dict[str, int]:
    """Parse DIAMOND phrog_hits.tsv -> {gene_id: phrog_id}.

    Best hit per gene (DIAMOND already returns top-1).
    """
    m: dict[str, int] = {}
    if not path.exists():
        return m
    with path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            gid, sid = f[0], f[1]
            mt = PHROG_ID_RE.search(sid)
            if mt:
                m[gid] = int(mt.group(1))
    return m


def biophysics(seq: str) -> str:
    """Format the trailing biophysics string for one gene.

    Returns a parenthesised string like:
      (Naa, pI=X.X, X.XkDa[, hydrophobic(GRAVY=...)][, ~NTM-helix])
    """
    seq_clean = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())
    n = len(seq_clean)
    if n < 5:
        return f"({len(seq)}aa)"
    pa = ProteinAnalysis(seq_clean)
    pi = pa.isoelectric_point()
    mw_kda = pa.molecular_weight() / 1000.0
    gravy = pa.gravy()
    parts = [f"{n}aa", f"pI={pi:.1f}", f"{mw_kda:.1f}kDa"]
    if gravy >= HYDROPHOBIC_GRAVY:
        parts.append(f"hydrophobic(GRAVY={gravy:.2f})")
    elif gravy <= HYDROPHILIC_GRAVY:
        parts.append(f"hydrophilic(GRAVY={gravy:.2f})")
    # TM-helix estimation: non-overlapping windows with mean KD > threshold
    n_tm = 0
    i = 0
    while i + TM_WINDOW <= n:
        window = seq_clean[i:i + TM_WINDOW]
        mean = sum(KD.get(c, 0.0) for c in window) / TM_WINDOW
        if mean >= TM_MIN_MEAN:
            n_tm += 1
            i += TM_WINDOW
        else:
            i += 1
    if n_tm >= 1:
        parts.append(f"~{n_tm}TM-helix")
    return "(" + ", ".join(parts) + ")"


def build_profile(
    phage_id: str,
    phage_dir: Path,
    lookup: dict[int, tuple[str, str]],
) -> str:
    """Build the full Markdown profile for one phage."""
    pdir = phage_dir / phage_id
    genes = parse_faa(pdir / "proteins.faa")
    hits = parse_phrog_hits(pdir / "phrog_hits.tsv")

    total = len(genes)
    n_phrog = sum(1 for g, _ in genes if g in hits)
    # 'Annotated' = PHROG hits with category != 'unknown function'
    n_annotated = 0
    sections: dict[str, list[str]] = {k: [] for k, _ in SECTION_ORDER}
    n_hypothetical = 0

    for gid, seq in genes:
        bio = biophysics(seq)
        pid = hits.get(gid)
        if pid is None:
            n_hypothetical += 1
            continue
        annot, cat = lookup.get(pid, ("NA", "unknown function"))
        if cat != "unknown function":
            n_annotated += 1
        section = CATEGORY_TO_SECTION.get(cat, "other")
        if annot in ("NA", "", None):
            phrog_part = f"[PHROG: phrog_{pid}]"
        else:
            phrog_part = f"[PHROG: {annot}]"
        line = f"- {gid} [LOW]: unknown {phrog_part}. {bio}"
        sections[section].append(line)

    # Assemble the profile text
    out: list[str] = []
    out.append(f"# Phage Genome Profile: {phage_id}")
    out.append(f"Total ORFs: {total} | Annotated: {n_annotated}/{total} | PHROG hits: {n_phrog}/{total}")
    out.append("Phage family (eggNOG): (eggnog unavailable)")
    out.append("Taxonomic context: (eggnog unavailable)")
    out.append("")

    for key, header in SECTION_ORDER:
        if not sections[key]:
            continue
        out.append(header)
        out.extend(sections[key])
        out.append("")

    out.append(f"## Hypothetical Proteins: {n_hypothetical} proteins (no annotation)")
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build per-phage base Markdown profiles from PHROG annotations.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--phage-dir",
        required=True,
        type=Path,
        help="Directory containing per-phage subdirectories, each with proteins.faa "
             "and optionally phrog_hits.tsv.",
    )
    parser.add_argument(
        "--phrog-lookup",
        required=True,
        type=Path,
        help="Path to phrog_lookup.tsv (columns: phrog_id, annotation, category).",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        type=Path,
        help="Output directory for generated .md profiles (flat: {out_dir}/{phage_id}.md).",
    )
    parser.add_argument(
        "--phage-list",
        default=None,
        type=Path,
        help="Optional text file listing phage IDs to process (one per line). "
             "If omitted, all subdirectories of --phage-dir are processed.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    phage_dir: Path = args.phage_dir
    phrog_lookup_path: Path = args.phrog_lookup
    out_dir: Path = args.out_dir
    phage_list_path: Path | None = args.phage_list

    # Validate inputs
    if not phage_dir.is_dir():
        print(f"ERROR: --phage-dir does not exist or is not a directory: {phage_dir}", file=sys.stderr)
        sys.exit(1)
    if not phrog_lookup_path.exists():
        print(f"ERROR: --phrog-lookup file not found: {phrog_lookup_path}", file=sys.stderr)
        sys.exit(1)

    # Determine which phage IDs to process
    if phage_list_path is not None:
        if not phage_list_path.exists():
            print(f"ERROR: --phage-list file not found: {phage_list_path}", file=sys.stderr)
            sys.exit(1)
        ids = [
            line.strip()
            for line in phage_list_path.read_text().splitlines()
            if line.strip()
        ]
    else:
        ids = sorted(p.name for p in phage_dir.iterdir() if p.is_dir())

    if not ids:
        print("WARNING: no phage IDs to process.", file=sys.stderr)
        sys.exit(0)

    # Create output directory
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load PHROG lookup table
    print(f"Loading PHROG lookup ({phrog_lookup_path})...")
    lookup = load_phrog_lookup(phrog_lookup_path)
    print(f"  {len(lookup)} phrog_id entries")

    # Build profiles
    print(f"Building base profiles for {len(ids)} phages...")
    n_ok = 0
    n_skip = 0

    for i, phage_id in enumerate(ids, 1):
        pdir = phage_dir / phage_id
        faa = pdir / "proteins.faa"
        if not faa.exists() or faa.stat().st_size == 0:
            n_skip += 1
            print(f"  SKIP {phage_id}: no proteins.faa")
            continue
        text = build_profile(phage_id, phage_dir, lookup)
        out_path = out_dir / f"{phage_id}.md"
        out_path.write_text(text, encoding="utf-8")
        n_ok += 1
        if i % 50 == 0 or i == len(ids):
            print(f"  [{i}/{len(ids)}] ok={n_ok} skip={n_skip}")

    print(f"Done. ok={n_ok} skip={n_skip}")
    print(f"  Output directory: {out_dir}")


if __name__ == "__main__":
    main()
