#!/usr/bin/env python3
"""
build_host_profiles_base.py
============================
Stage 1 of the PHI-Reason host profile pipeline.

Generates base host text profiles from:
  1. eggNOG-mapper functional annotations (.emapper.annotations)
  2. DefenseFinder anti-phage defense system predictions

Output: one Markdown file per host with two sections:
  ## Defense Systems [Phage resistance — DefenseFinder]
      Lists detected defense systems grouped by type (RM, CRISPR-Cas, Thoeris, etc.)
  ## Surface Receptors & Cell Envelope [Phage attachment targets]
      Lists surface-relevant proteins with [HIGH] confidence tags
      (outer membrane proteins, LPS biosynthesis, pili, flagella, capsule, WTA, etc.)

These base profiles are the input to build_host_profiles.py (Stage 2), which
enriches them with NCBI taxonomy, curated receptor knowledge, and phage biology
context for LLM-based phage-host prediction.

Prerequisites (external tools):
  1. Gene prediction: Prodigal or Prokka
       prodigal -i {host}.fna -a {host}.faa -p meta
  2. Functional annotation: eggNOG-mapper v2
       emapper.py -i {host}.faa --output {host} -m diamond --cpu 8
  3. Defense system prediction: DefenseFinder v1.2+
       defense-finder run -i {host}.faa -o {host_dir}/

File naming conventions (adjust --anno-suffix and --defense-suffix as needed):
  eggNOG:       {anno-dir}/{host_id}{anno-suffix}.emapper.annotations
  DefenseFinder:{defense-dir}/{host_id}/{host_id}{defense-suffix}_defense_finder_systems.tsv

Usage:
    python build_host_profiles_base.py \\
        --anno-dir     /path/to/emapper_results \\
        --defense-dir  /path/to/defensefinder_results \\
        --out-dir      /path/to/host_base_profiles \\
        [--host-list   hosts.json] \\
        [--anno-suffix ""] \\
        [--defense-suffix ""] \\
        [--dry-run]
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

# ── eggNOG-mapper output column names (v2, default 21-column format) ─────────
ANNO_COLS = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
    "max_annot_lvl", "COG_category", "Description", "Preferred_name",
    "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
    "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC",
    "CAZy", "BiGG_Reaction", "PFAMs",
]

# ── Keyword sets for gene classification ──────────────────────────────────────

# Surface / cell-envelope genes relevant to phage adsorption
_SURFACE_KW = [
    # Outer membrane proteins (Gram-negative)
    "outer membrane protein", "outer membrane porin", "outer membrane receptor",
    "outer membrane channel", "outer membrane assembly", "tonb-dependent",
    "ompa", "ompb", "ompc", "ompd", "ompf", "ompt", "omps", "ompn",
    "fhua", "fepa", "cira", "bfeb", "bfec", "bfed", "fpva", "fpvb",
    "ferrichrome outer membrane", "iron-siderophore outer membrane",
    # Porins
    "porin", "general diffusion porin",
    # LPS / surface polymers
    "lipopolysaccharide biosynthesis", "lps biosynthesis",
    "o-antigen", "lipid a biosynthesis", "kdo",
    "capsular polysaccharide", "capsule biosynthesis",
    "peptidoglycan", "murein", "n-acetylglucosamine",
    # Flagella (structural)
    "flagellin", "flagellar hook", "flagellar basal", "flagellar filament",
    "flagellar motor", "flagellar export",
    # Pili / fimbriae (structural)
    "type iv pilus", "type iv pilin", "pilus assembly",
    "fimbrial adhesin", "chaperone-usher fimbrial",
    "type i fimbriae", "p fimbriae", "curli", "conjugative pili",
    "type iii secretion", "needle complex",
    # Cell wall components (Gram-positive)
    "teichoic acid", "wall teichoic", "lipoteichoic",
    "peptidoglycan hydrolase", "murein hydrolase",
    "autolysin", "N-acetylmuramidase",
    # Gram-positive surface proteins
    "sortase", "lpxtg", "MSCRAMM",
    # Glycan / capsule (general)
    "exopolysaccharide", "glucan synthase", "cellulose synthase",
    "alginate", "colanic acid",
    # S-layer (archaea + some bacteria)
    "s-layer", "surface array protein",
    # TonB complex (siderophore-receptor link)
    "tonb", "exbb", "exbd",
    # Lipoprotein / outer membrane biogenesis
    "lipoprotein", "bam complex", "beta-barrel assembly",
    "signal peptidase", "lgt", "lsp",
]

# Defense system genes (eggNOG keyword-based; DefenseFinder TSV takes priority)
_DEFENSE_KW = [
    "restriction", "modification enzyme", "methyltransferase type i",
    "methyltransferase type ii", "methyltransferase type iii",
    "anti-crispr", "crispr", "cas1", "cas2", "cas3", "cas9", "cas12", "cas13",
    "abortive infection", "abi",
    "thoeris", "pycsar", "cbass", "retron",
    "defense island", "anti-phage",
    "bacteriophage exclusion", "hsd", "hsdr", "hsdm", "hsds",
    "dpn", "mbol", "ecori",
    "phosphorothioation", "dnd",
    "zorya", "druantia", "gabija", "shedu", "lamassu",
    "wadjet", "doron", "kiwa", "avs",
    "bacteriocin",
    "dna modification", "dam methylase", "dcm methylase",
]

# ── Confidence tiers ──────────────────────────────────────────────────────────
# [HIGH]: gene is unambiguously surface-exposed and phage-relevant
_HIGH_KW = [
    "outer membrane protein", "outer membrane porin", "tonb-dependent",
    "ompa", "ompb", "ompc", "ompd", "ompf", "ompt", "ompn", "fhua", "fepa",
    "lipopolysaccharide", "o-antigen", "lipid a", "kdo",
    "capsular polysaccharide", "capsule biosynthesis",
    "flagellin", "flagellar hook", "flagellar filament",
    "type iv pilus", "type iv pilin",
    "teichoic acid", "wall teichoic", "lipoteichoic",
    "porin", "bam complex", "outer membrane assembly",
    "s-layer", "sortase",
]


def _parse_annotations(anno_path: Path) -> dict[str, dict]:
    """Parse eggNOG-mapper .emapper.annotations → {gene_id: {col: val}}."""
    genes: dict[str, dict] = {}
    if not anno_path.exists():
        return genes
    with open(anno_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            rec = dict(zip(ANNO_COLS, parts + [""] * max(0, len(ANNO_COLS) - len(parts))))
            genes[rec["query"]] = rec
    return genes


def _gene_label(rec: dict) -> str | None:
    """Pick best human-readable label from an eggNOG annotation record."""
    name = rec.get("Preferred_name", "-").strip()
    desc = rec.get("Description", "-").strip()
    pfam = rec.get("PFAMs", "-").strip()
    if name and name != "-":
        return f"{name} ({desc})" if (desc and desc != "-") else name
    if desc and desc != "-":
        return desc
    if pfam and pfam != "-":
        return f"domain: {pfam}"
    return None


def _classify(rec: dict) -> str | None:
    """Return 'surface', 'defense', or None."""
    label = (_gene_label(rec) or "").lower()
    pfam  = rec.get("PFAMs", "").lower()
    desc  = rec.get("Description", "").lower()
    combined = f"{label} {pfam} {desc}"

    for kw in _DEFENSE_KW:
        if kw in combined:
            return "defense"
    for kw in _SURFACE_KW:
        if kw in combined:
            return "surface"
    return None


def _is_high_confidence(label: str) -> bool:
    ll = label.lower()
    return any(kw in ll for kw in _HIGH_KW)


def _load_defensefinder(tsv_path: Path) -> list[dict]:
    """Parse a DefenseFinder *_defense_finder_systems.tsv into a list of records."""
    if not tsv_path.exists():
        return []
    systems = []
    header = None
    with open(tsv_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if header is None:
                header = cols
                continue
            systems.append(dict(zip(header, cols)))
    return systems


def _format_defense_section(systems: list[dict]) -> str:
    """Format DefenseFinder results as the Defense Systems profile section."""
    if not systems:
        return "## Defense Systems [Phage resistance — DefenseFinder]\n- (none detected)\n"

    by_type: dict[str, list] = defaultdict(list)
    for s in systems:
        by_type[s.get("type", "?")].append(s)

    lines = ["## Defense Systems [Phage resistance — DefenseFinder]"]
    for stype, entries in sorted(by_type.items()):
        subtypes = sorted(set(e.get("subtype", stype) for e in entries))
        n = len(entries)
        if len(subtypes) <= 3:
            sub_str = ", ".join(subtypes)
        else:
            sub_str = f"{subtypes[0]}...+{len(subtypes) - 1}"
        lines.append(f"- {stype} ({n} system{'s' if n > 1 else ''}): {sub_str}")
    lines.append("")
    return "\n".join(lines)


def build_base_profile(
    host_id: str,
    anno_path: Path,
    defense_tsv: Path | None,
) -> str:
    """Build the base host profile for one host species.

    Args:
        host_id:     Species identifier, e.g. 'Escherichia_coli'
        anno_path:   Path to eggNOG .emapper.annotations file
        defense_tsv: Path to DefenseFinder *_defense_finder_systems.tsv (may be None)

    Returns:
        Markdown string with Defense Systems + Surface Receptors sections.
    """
    anno = _parse_annotations(anno_path)

    # ── Surface receptor genes ─────────────────────────────────────────────────
    surface: list[tuple[str, str, bool]] = []  # (gene_id, label, is_high)
    n_total = 0

    for gene_id, rec in anno.items():
        n_total += 1
        cat = _classify(rec)
        if cat == "surface":
            lbl = _gene_label(rec)
            if lbl:
                surface.append((gene_id, lbl, _is_high_confidence(lbl)))

    # ── Defense systems (prefer DefenseFinder over eggNOG keywords) ───────────
    defense_systems = _load_defensefinder(defense_tsv) if defense_tsv else []

    # ── Assemble profile ───────────────────────────────────────────────────────
    n_surface  = len(surface)
    n_defense  = sum(1 for rec in anno.values() if _classify(rec) == "defense")

    lines = [
        f"# Host Interaction Profile: {host_id}",
        f"Total ORFs: {n_total} | Surface-relevant: {n_surface} | Defense-relevant: {n_defense}",
        "",
    ]

    # Defense Systems section
    lines.append(_format_defense_section(defense_systems))

    # Surface Receptors section
    lines.append("## Surface Receptors & Cell Envelope [Phage attachment targets]")
    if surface:
        seen: set[str] = set()
        shown = 0
        for gene_id, lbl, high in surface:
            key = lbl[:50].lower()
            if key in seen:
                continue
            seen.add(key)
            tag = "[HIGH] " if high else ""
            lines.append(f"- {tag}{gene_id} ({lbl[:120]})")
            shown += 1
            if shown >= 30:
                remaining = n_surface - shown
                if remaining > 0:
                    lines.append(f"- ... and {remaining} more surface/envelope-related proteins")
                break
    else:
        lines.append("- (No surface receptor annotations detected)")
    lines.append("")

    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build base host text profiles for LLM phage-host prediction "
            "(Stage 1: eggNOG annotations + DefenseFinder → base profiles)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--anno-dir",      required=True,
                        help="Directory of eggNOG-mapper output (.emapper.annotations files)")
    parser.add_argument("--defense-dir",   default=None,
                        help="Directory of DefenseFinder results "
                             "({defense-dir}/{host_id}/*_defense_finder_systems.tsv). "
                             "If omitted, defense sections are left empty.")
    parser.add_argument("--out-dir",       required=True,
                        help="Output directory for base profile .md files")
    parser.add_argument("--host-list",     default=None,
                        help="JSON file listing host IDs to process "
                             "(default: all hosts with annotations in --anno-dir)")
    parser.add_argument("--anno-suffix",   default="",
                        help="Suffix between host_id and .emapper.annotations, "
                             "e.g. '_emapper' → {host_id}_emapper.emapper.annotations")
    parser.add_argument("--defense-suffix", default="",
                        help="Suffix inserted before _defense_finder_systems.tsv, "
                             "e.g. '_df' → {host_id}_df_defense_finder_systems.tsv")
    parser.add_argument("--dry-run",       action="store_true",
                        help="Preview first 2 profiles without writing files")
    args = parser.parse_args()

    anno_dir    = Path(args.anno_dir)
    defense_dir = Path(args.defense_dir) if args.defense_dir else None
    out_dir     = Path(args.out_dir)

    # Build host list
    if args.host_list:
        host_ids: list[str] = json.loads(Path(args.host_list).read_text())
        print(f"Host list: {len(host_ids)} hosts from {args.host_list}")
    else:
        host_ids = sorted(
            p.stem.removesuffix(args.anno_suffix)
            for p in anno_dir.glob(f"*{args.anno_suffix}.emapper.annotations")
        )
        print(f"Host list: {len(host_ids)} hosts discovered in {anno_dir}")

    print(f"Anno dir    : {anno_dir}")
    print(f"Defense dir : {defense_dir or '(none)'}")
    print(f"Output dir  : {out_dir}")
    print()

    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)

    written = missing_anno = preview_shown = 0

    for host_id in host_ids:
        anno_path = anno_dir / f"{host_id}{args.anno_suffix}.emapper.annotations"
        if not anno_path.exists():
            print(f"  MISSING annotation: {host_id}")
            missing_anno += 1
            continue

        # DefenseFinder TSV: {defense_dir}/{host_id}/{host_id}{suffix}_defense_finder_systems.tsv
        defense_tsv = None
        if defense_dir:
            defense_tsv = defense_dir / host_id / f"{host_id}{args.defense_suffix}_defense_finder_systems.tsv"

        text = build_base_profile(host_id, anno_path, defense_tsv)

        if args.dry_run:
            if preview_shown < 2:
                print(f"  [DRY] {host_id}")
                print(text[:600])
                print("  ...")
                preview_shown += 1
        else:
            (out_dir / f"{host_id}.md").write_text(text)

        written += 1

    print()
    if not args.dry_run:
        print("Done.")
        print(f"  Profiles written    : {written}")
        print(f"  Missing annotations : {missing_anno}")
        print(f"  Output              : {out_dir}")
    else:
        print(f"[DRY RUN] Would write {written} profiles.")


if __name__ == "__main__":
    main()
