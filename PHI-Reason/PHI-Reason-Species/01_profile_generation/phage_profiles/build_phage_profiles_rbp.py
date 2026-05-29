#!/usr/bin/env python3
"""
build_phage_profiles.py
=======================
Generates phage genome text profiles for LLM-based phage-host prediction.

This script is the primary profile builder (no-BLASTN variant). It produces
one Markdown profile per phage containing:
  - Genome overview (ORF count, taxonomic context from eggNOG)
  - Functional sections (Tail & Host Adsorption, Lysis Cassette, DNA Packaging, etc.)
  - Inline RBP homology evidence on tail/adsorption genes, from BLASTP against
    a curated phage RBP reference database with experimentally confirmed host genera
  - Hypothetical protein count summary

RBP evidence format:
  ← RBP matches (BLASTP): Genus (identity=X%, qcov=Y%, N hits)
  where "N hits" is the number of reference RBPs from phages infecting that genus
  that matched the query gene by sequence similarity.

Data leakage policy:
  The RBP reference database (rbp_phage_host_table.csv) must have zero phage
  accession overlap with the evaluation dataset. Verify before use.

Pipeline phases:
  1. Extract RBP sequences from ColabFold/AlphaFold PDB structures → FASTA
  2. Build DIAMOND protein database from RBP FASTA
  3. Extract tail/adsorption proteins per phage (PHROG-guided)
  4. Run DIAMOND BLASTP (records identity, query coverage, e-value)
  5. Build final profiles: annotate tail genes with RBP BLASTP evidence,
     fix annotation quality tags, collapse hypothetical proteins section

Usage:
    python build_phage_profiles.py \\
        --base-profiles  /path/to/textGeneProfile_v2_correct_mask \\
        --prot-dir       /path/to/phage_idr/prot \\
        --phrog-json     /path/to/phrog_by_phage.json \\
        --rbp-table      /path/to/rbp_phage_host_table.csv \\
        --rbp-structs    /path/to/rbp_structures \\
        --work-dir       /path/to/work_dir \\
        --out-dir        /path/to/output_profiles \\
        [--phage-list    phage_ids.json] \\
        [--skip-db] [--skip-blast] [--threads 16]

    # Skip DB and BLAST if already computed:
    python build_phage_profiles.py ... --skip-db --skip-blast

Output:
    One .md file per phage in --out-dir, ready for LLM host prediction.
    Intermediate files (DIAMOND db, BLASTP results) are written to --work-dir.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

DIAMOND = Path(os.environ.get("DIAMOND_BIN",
               shutil.which("diamond") or "diamond"))

# ── BLASTP filtering thresholds ───────────────────────────────────────────────
MIN_PIDENT   = 40.0   # minimum % identity to include a hit (aligns with PHIstruct's lowest similarity bin)
MIN_QCOV     = 5.0    # minimum query coverage (%) — loose floor only; real ranking uses rank_score
MAX_EVALUE   = 1e-5   # maximum e-value
TOP_HOSTS    = 5      # max distinct host genera shown per gene

# ── PHROG product keywords identifying tail / host-adsorption genes ───────────
RBP_PRODUCTS = [
    "tail fiber", "tailspike", "tail spike", "receptor binding",
    "distal tail protein", "baseplate spike", "host range",
    "long tail fiber", "short tail fiber", "host specificity",
    "spike protein", "needle protein", "tail needle",
    "tip protein", "adsorption",
]

# ── Three-letter to one-letter amino acid code (for PDB parsing) ──────────────
_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O", "UNK": "X",
}

# ── [LOW]-tag logic: PHROG canonical short product names ─────────────────────
_PHROG_SHORT_FUNCS = {
    "holin", "endolysin", "lysin", "spanin", "lysozyme", "muramidase",
    "tail fiber", "tail fiber protein", "tail spike", "tailspike",
    "tail spike protein", "baseplate spike", "receptor binding protein",
    "major capsid protein", "minor capsid protein", "major head protein",
    "minor head protein", "minor tail protein", "tape measure protein",
    "portal protein", "terminase large subunit", "terminase small subunit",
    "head completion protein", "neck protein", "collar protein",
    "lysis inhibition protein", "rz-like spanin", "spanin outer",
    "spanin inner", "holin-like", "anti-holin",
}


# ─────────────────────────────────────────────────────────────────────────────
# Phase 1: Extract RBP sequences from PDB structures
# ─────────────────────────────────────────────────────────────────────────────

def _extract_seq_from_pdb(pdb_path: Path) -> str:
    """Read CA atoms from a PDB file and return the amino acid sequence."""
    residues: dict[tuple, str] = {}
    try:
        with open(pdb_path) as fh:
            for line in fh:
                if not line.startswith("ATOM "):
                    continue
                if line[12:16].strip() != "CA":
                    continue
                try:
                    res_num = int(line[22:26].strip())
                except ValueError:
                    continue
                key = (line[21], res_num)
                if key not in residues:
                    residues[key] = line[17:20].strip()
    except Exception:
        return ""
    if not residues:
        return ""
    return "".join(_AA3TO1.get(r, "X") for _, r in sorted(residues.items()))


def build_rbp_fasta(table_path: Path, structs_dir: Path, out_fasta: Path) -> int:
    """Extract RBP sequences from PDB structures and write a FASTA file.

    FASTA header format: >protein_id|phage_accession|host_genus
    Returns the number of sequences written.
    """
    print("Phase 1: Extracting RBP sequences from PDB structures...")
    rows: list[tuple[str, str, str]] = []
    with open(table_path) as fh:
        next(fh)  # skip header line
        for line in fh:
            parts = line.strip().split(",")
            if len(parts) >= 3:
                rows.append((parts[0].strip(), parts[1].strip(), parts[2].strip().lower()))

    written = failed = 0
    with open(out_fasta, "w") as fout:
        for i, (pid, pacc, host) in enumerate(rows):
            if i % 2000 == 0:
                print(f"  {i}/{len(rows)} processed...", flush=True)
            pdb_file = structs_dir / f"{pid}.pdb"
            if not pdb_file.exists():
                failed += 1
                continue
            seq = _extract_seq_from_pdb(pdb_file)
            if len(seq) < 30:
                failed += 1
                continue
            fout.write(f">{pid}|{pacc}|{host}\n{seq}\n")
            written += 1

    print(f"  Written: {written}  |  Failed/short: {failed}")
    return written


# ─────────────────────────────────────────────────────────────────────────────
# Phase 2: Build DIAMOND protein database
# ─────────────────────────────────────────────────────────────────────────────

def build_diamond_db(fasta_path: Path, db_path: Path, threads: int = 8) -> None:
    """Build a DIAMOND protein database from the RBP FASTA."""
    print(f"Phase 2: Building DIAMOND database → {db_path}")
    subprocess.run(
        [str(DIAMOND), "makedb",
         "--in", str(fasta_path),
         "--db", str(db_path),
         "--threads", str(threads),
         "--quiet"],
        check=True,
    )
    print("  Done.")


# ─────────────────────────────────────────────────────────────────────────────
# Phase 3: Extract tail/adsorption proteins per phage
# ─────────────────────────────────────────────────────────────────────────────

def load_phrog_tail_genes(phrog_json: Path) -> dict[str, dict[str, str]]:
    """Return {phage_id: {gene_id: product_str}} for all PHROG-annotated tail genes."""
    print("Phase 3a: Loading PHROG tail annotations...")
    with open(phrog_json) as fh:
        phrogs: dict = json.load(fh)

    tail_genes: dict[str, dict[str, str]] = {}
    for phage_id, genes in phrogs.items():
        tg = {
            gid: ann.get("product", "tail protein")
            for gid, ann in genes.items()
            if ann.get("function", "").lower() == "tail"
            or any(kw in ann.get("product", "").lower() for kw in RBP_PRODUCTS)
        }
        if tg:
            tail_genes[phage_id] = tg

    print(f"  Phages with ≥1 PHROG tail gene: {len(tail_genes)}/{len(phrogs)}")
    return tail_genes


def _read_fasta(path: Path) -> dict[str, str]:
    """Parse a FASTA file into {seq_id: sequence}."""
    seqs: dict[str, str] = {}
    cur_id, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(buf)
                cur_id = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur_id:
        seqs[cur_id] = "".join(buf)
    return seqs


def extract_tail_proteins(
    phage_ids: list[str],
    tail_genes: dict[str, dict[str, str]],
    prot_dir: Path,
    out_fasta: Path,
) -> dict[str, dict[str, int]]:
    """Extract tail-gene protein sequences and write to a combined FASTA.

    Header format: >phage_id::gene_id
    Returns {phage_id: {gene_id: aa_length}} for sequences written.
    """
    print(f"Phase 3b: Extracting tail proteins for {len(phage_ids)} phages...")
    written_map: dict[str, dict[str, int]] = {}
    no_tail: list[str] = []

    with open(out_fasta, "w") as fout:
        for phage_id in phage_ids:
            faa_path = prot_dir / phage_id / "proteins.faa"
            if not faa_path.exists():
                continue
            seqs = _read_fasta(faa_path)
            target = {g: p for g, p in tail_genes.get(phage_id, {}).items() if g in seqs}
            if not target:
                no_tail.append(phage_id)
                continue
            written_map[phage_id] = {}
            for gene_id in target:
                seq = seqs[gene_id]
                if len(seq) < 20:
                    continue
                fout.write(f">{phage_id}::{gene_id}\n{seq}\n")
                written_map[phage_id][gene_id] = len(seq)

    total_genes = sum(len(v) for v in written_map.values())
    print(f"  Phages with tail genes extracted : {len(written_map)}")
    print(f"  Phages with no PHROG tail genes  : {len(no_tail)}")
    print(f"  Total tail gene sequences        : {total_genes}")
    return written_map


def _reconstruct_written_genes(tail_fasta: Path) -> dict[str, dict[str, int]]:
    """Rebuild written_map from an existing tail protein FASTA (for --skip-blast)."""
    written: dict[str, dict[str, int]] = defaultdict(dict)
    cur_phage = cur_gene = None
    cur_len = 0
    with open(tail_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_phage and cur_gene:
                    written[cur_phage][cur_gene] = cur_len
                parts = line[1:].strip().split("::")
                cur_phage, cur_gene = (parts[0], parts[1]) if len(parts) >= 2 else (None, None)
                cur_len = 0
            else:
                cur_len += len(line.strip())
    if cur_phage and cur_gene:
        written[cur_phage][cur_gene] = cur_len
    return dict(written)


# ─────────────────────────────────────────────────────────────────────────────
# Phase 4: DIAMOND BLASTP
# ─────────────────────────────────────────────────────────────────────────────

def run_diamond_blastp(
    query_fasta: Path,
    db_path: Path,
    out_tsv: Path,
    threads: int = 16,
) -> None:
    """Run DIAMOND BLASTP: query phage tail proteins vs RBP reference database.

    Output columns: qseqid sseqid pident qcovhsp evalue length
    """
    print("Phase 4: Running DIAMOND BLASTP...")
    subprocess.run(
        [str(DIAMOND), "blastp",
         "--db",      str(db_path),
         "--query",   str(query_fasta),
         "--out",     str(out_tsv),
         "--outfmt",  "6",
         "qseqid", "sseqid", "pident", "qcovhsp", "evalue", "length",
         "--evalue",  str(MAX_EVALUE),
         "--id",      str(int(MIN_PIDENT)),
         "-k",        "20",
         "--threads", str(threads),
         "--quiet"],
        check=True,
    )
    print(f"  Done → {out_tsv}")


def parse_blast_results(
    blast_tsv: Path,
    rbp_table: Path,
) -> dict[str, dict[str, list[tuple[str, float, float, int]]]]:
    """Parse DIAMOND output and aggregate hits by phage/gene/host genus.

    Returns:
        {phage_id: {gene_id: [(host_genus, best_pident, best_qcov, hit_count), ...]}}
        sorted by rank_score DESC, where rank_score = pident^0.6 * qcov^0.4.
    """
    print("Phase 4b: Parsing BLAST results...")

    # Fallback protein_id → host_genus lookup (in case header parsing fails)
    pid2host: dict[str, str] = {}
    with open(rbp_table) as fh:
        next(fh)
        for line in fh:
            parts = line.strip().split(",")
            if len(parts) >= 3:
                pid2host[parts[0].strip()] = parts[2].strip().lower()

    # Columns: qseqid sseqid pident qcovhsp evalue length
    gene_hits: dict = defaultdict(lambda: defaultdict(list))
    with open(blast_tsv) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            qseqid = parts[0]
            sseqid = parts[1]
            pident = float(parts[2])
            qcov   = float(parts[3])
            # evalue = float(parts[4])  # already filtered by DIAMOND

            if "::" not in qseqid:
                continue
            phage_id, gene_id = qseqid.split("::", 1)

            if qcov < MIN_QCOV:
                continue

            sparts = sseqid.split("|")
            host_genus = sparts[2] if len(sparts) >= 3 else pid2host.get(sparts[0], "unknown")

            gene_hits[phage_id][gene_id].append((host_genus, pident, qcov))

    # Aggregate per host genus: best_pident, best_qcov, hit_count
    result: dict[str, dict[str, list]] = {}
    for phage_id, genes in gene_hits.items():
        result[phage_id] = {}
        for gene_id, hits in genes.items():
            agg: dict[str, list] = defaultdict(lambda: [0.0, 0.0, 0])
            for host, pident, qcov in hits:
                a = agg[host]
                if pident > a[0]:
                    a[0] = pident
                if qcov > a[1]:
                    a[1] = qcov
                a[2] += 1
            sorted_hosts = sorted(
                [(h, v[0], v[1], v[2]) for h, v in agg.items()],
                key=lambda x: -(x[1] ** 0.6 * x[2] ** 0.4),
            )
            result[phage_id][gene_id] = sorted_hosts[:TOP_HOSTS]

    n_genes = sum(len(v) for v in result.values())
    print(f"  Phages with BLAST hits: {len(result)}")
    print(f"  Genes  with BLAST hits: {n_genes}")
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Phase 5: Build annotated profiles
# ─────────────────────────────────────────────────────────────────────────────

def _format_genus_hit(genus: str, pident: float, qcov: float, count: int) -> str:
    """Format a single genus entry in the BLASTP annotation.

    Example: 'Lactobacillus (identity=100%, qcov=87%, 5 hits)'
             'Streptomyces (identity=73%, qcov=62%)'   ← single hit, count omitted
    """
    parts = [f"identity={pident:.0f}%", f"qcov={qcov:.0f}%"]
    if count > 1:
        parts.append(f"{count} hits")
    return f"{genus.capitalize()} ({', '.join(parts)})"


def _build_blastp_annotation(gene_hits: list[tuple[str, float, float, int]]) -> str:
    """Return the inline BLASTP annotation string for one gene."""
    items = [_format_genus_hit(h, p, c, n) for h, p, c, n in gene_hits]
    return " ← RBP matches (BLASTP): " + ", ".join(items)


def _is_phrog_short_name(label: str) -> bool:
    """Return True if label looks like a bare PHROG product name (no 'putative' yet)."""
    ll = label.lower().strip().rstrip(".")
    if "putative" in ll or "[phrog:" in ll:
        return False
    if ll.startswith("domain:") or ll.startswith("pfam:"):
        return False
    if any(ll == fn or ll.startswith(fn + " ") for fn in _PHROG_SHORT_FUNCS):
        return True
    return len(ll) <= 30 and not any(
        c in ll for c in ["(", ")", "family", "domain", "protein i", "protein iv"]
    )


def _fix_low_tag(line: str) -> str:
    """Convert '- geneN [LOW]: label ...' to '- geneN : [putative] label ...'

    [LOW] signals that the eggNOG annotation has low confidence. We either
    prefix 'putative' (when PHROG short name) or simply drop the tag.
    """
    m = re.match(r"^(- \S+) \[LOW\]: (.+)$", line)
    if not m:
        return line
    prefix, rest = m.group(1), m.group(2)

    # Split label from trailing metadata (" [hint] (features) ← ...")
    label_end = re.search(r"\.\s+(?=[\[\(]|←)", rest)
    if label_end:
        label = rest[: label_end.start()]
        tail  = rest[label_end.start():]
    else:
        label = rest.rstrip(".")
        tail  = "."

    ll = label.lower()
    if "putative" in ll or "[phrog:" in ll:
        fixed_label = label
    elif _is_phrog_short_name(label):
        fixed_label = f"putative {label}"
    else:
        fixed_label = label

    return f"{prefix} : {fixed_label}{tail}"


def process_profile(
    base_text: str,
    blastp_hits: dict[str, list],   # {gene_id: [(host, pident, qcov, count), ...]}
) -> str:
    """Apply all transformations to a base v2cm profile text.

    Transformations:
      1. Insert inline BLASTP RBP annotations on tail/adsorption gene lines
      2. Fix [LOW] quality tags
      3. Collapse '## Hypothetical Proteins' section to a count summary
    """
    lines = base_text.split("\n")
    out: list[str] = []
    in_hypo = False
    hypo_count = 0

    for line in lines:
        # ── Hypothetical Proteins section ─────────────────────────────────────
        if re.match(r"^## Hypothetical Proteins", line):
            in_hypo = True
            hypo_count = 0
            out.append("__HYPO_PLACEHOLDER__")
            continue

        if in_hypo:
            if line.startswith("## "):
                # End of hypothetical block — replace placeholder with count
                for i in range(len(out) - 1, -1, -1):
                    if out[i] == "__HYPO_PLACEHOLDER__":
                        out[i] = f"## Hypothetical Proteins: {hypo_count} proteins (no annotation)"
                        break
                in_hypo = False
                # Fall through to process current section header normally
            elif line.startswith("- "):
                hypo_count += 1
                continue
            else:
                continue  # skip blank lines within hypo section

        # ── Fix [LOW] confidence tag ───────────────────────────────────────────
        if " [LOW]: " in line:
            line = _fix_low_tag(line)

        # ── Append BLASTP annotation to tail gene lines ────────────────────────
        m = re.match(r"^(- (gene\d+) )", line)
        if m:
            gene_id = m.group(2)
            if gene_id in blastp_hits:
                line = line.rstrip() + _build_blastp_annotation(blastp_hits[gene_id])

        out.append(line)

    # End of file still inside hypothetical section
    if in_hypo:
        for i in range(len(out) - 1, -1, -1):
            if out[i] == "__HYPO_PLACEHOLDER__":
                out[i] = f"## Hypothetical Proteins: {hypo_count} proteins (no annotation)"
                break

    return "\n".join(out)


def build_profiles(
    phage_ids: list[str],
    blast_hits: dict[str, dict[str, list]],
    base_profiles_dir: Path,
    out_dir: Path,
) -> dict[str, int]:
    """Read base v2cm profiles, apply BLASTP annotations, write final profiles."""
    print(f"Phase 5: Building profiles → {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    stats = {"augmented": 0, "no_blast_hit": 0, "missing_base": 0}
    for phage_id in phage_ids:
        base_file = base_profiles_dir / f"{phage_id}.md"
        if not base_file.exists():
            stats["missing_base"] += 1
            continue

        gene_hits = blast_hits.get(phage_id, {})
        text = process_profile(base_file.read_text(), gene_hits)
        (out_dir / f"{phage_id}.md").write_text(text)

        if gene_hits:
            stats["augmented"] += 1
        else:
            stats["no_blast_hit"] += 1

    return stats


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build phage genome text profiles for LLM host prediction (no-BLASTN variant).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--base-profiles", required=True,
                        help="Directory of base v2cm phage profiles (.md files)")
    parser.add_argument("--prot-dir",      required=True,
                        help="Directory of per-phage protein FASTA files")
    parser.add_argument("--phrog-json",    required=True,
                        help="JSON file: {phage_id: {gene_id: {function, product, ...}}}")
    parser.add_argument("--rbp-table",     required=True,
                        help="CSV: protein_id, phage_accession, host_genus  (RBP reference)")
    parser.add_argument("--rbp-structs",   required=True,
                        help="Directory of RBP PDB structure files")
    parser.add_argument("--work-dir",      required=True,
                        help="Working directory for intermediate files (DIAMOND db, BLASTP TSV)")
    parser.add_argument("--out-dir",       required=True,
                        help="Output directory for final .md profiles")
    parser.add_argument("--phage-list",    default=None,
                        help="JSON file listing phage IDs to process (default: all in base-profiles)")
    parser.add_argument("--skip-db",       action="store_true",
                        help="Skip phases 1-2 (DIAMOND DB already built in work-dir)")
    parser.add_argument("--skip-blast",    action="store_true",
                        help="Skip phases 3-4 (BLASTP results already in work-dir)")
    parser.add_argument("--threads",       type=int, default=16)
    args = parser.parse_args()

    base_profiles_dir = Path(args.base_profiles)
    prot_dir          = Path(args.prot_dir)
    phrog_json        = Path(args.phrog_json)
    rbp_table         = Path(args.rbp_table)
    rbp_structs       = Path(args.rbp_structs)
    work_dir          = Path(args.work_dir)
    out_dir           = Path(args.out_dir)

    work_dir.mkdir(parents=True, exist_ok=True)

    rbp_fasta  = work_dir / "rbp_reference.fasta"
    rbp_dmnd   = work_dir / "rbp_reference"        # DIAMOND appends .dmnd
    tail_fasta = work_dir / "tail_proteins.faa"
    blast_out  = work_dir / "blastp_results.tsv"

    # Determine phage list
    if args.phage_list:
        phage_ids: list[str] = json.loads(Path(args.phage_list).read_text())
        print(f"Phage list: {len(phage_ids)} phages from {args.phage_list}")
    else:
        phage_ids = [p.stem for p in base_profiles_dir.glob("*.md")]
        print(f"Phage list: all {len(phage_ids)} phages in {base_profiles_dir}")

    # ── Phases 1-2: Build RBP reference database ──────────────────────────────
    if not args.skip_db:
        n = build_rbp_fasta(rbp_table, rbp_structs, rbp_fasta)
        if n < 100:
            print(f"ERROR: Only {n} sequences extracted — check PDB files.", file=sys.stderr)
            sys.exit(1)
        build_diamond_db(rbp_fasta, rbp_dmnd, threads=args.threads)
    else:
        print("Phases 1-2 skipped (--skip-db)")

    # ── Phase 3: Extract tail proteins ────────────────────────────────────────
    phrog_tail = load_phrog_tail_genes(phrog_json)

    if not args.skip_blast:
        written_genes = extract_tail_proteins(phage_ids, phrog_tail, prot_dir, tail_fasta)

        # ── Phase 4: DIAMOND BLASTP ───────────────────────────────────────────
        if tail_fasta.stat().st_size > 0:
            run_diamond_blastp(tail_fasta, rbp_dmnd, blast_out, threads=args.threads)
        else:
            print("WARNING: No tail proteins extracted — BLASTP skipped.")
            blast_out.write_text("")
    else:
        print("Phases 3-4 skipped (--skip-blast)")
        written_genes = _reconstruct_written_genes(tail_fasta) if tail_fasta.exists() else {}

    # ── Phase 4b + 5: Parse results and build profiles ────────────────────────
    blast_hits: dict = {}
    if blast_out.exists() and blast_out.stat().st_size > 0:
        blast_hits = parse_blast_results(blast_out, rbp_table)

    stats = build_profiles(phage_ids, blast_hits, base_profiles_dir, out_dir)

    total = len(phage_ids)
    print("\nDone.")
    print(f"  Augmented with BLASTP hits : {stats['augmented']}/{total} ({stats['augmented']/total:.1%})")
    print(f"  No BLASTP hits             : {stats['no_blast_hit']}/{total}")
    print(f"  Missing base profile       : {stats['missing_base']}")
    print(f"  Output                     : {out_dir}")


if __name__ == "__main__":
    main()
