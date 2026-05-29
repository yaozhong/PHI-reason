#!/usr/bin/env python3
"""
01_build_db.py — Build a protein DIAMOND database from a training set of phages.

Inputs (read-only):
  --train-dir    Directory holding train-set phage subdirectories (each with proteins.faa)
  --train-list   (optional) Text file listing train accessions (one per line).
                 If omitted, all subdirectories in --train-dir are used.
  --pair-csv     Phage-host pair CSV (accession,host_species)
  --out-dir      Output directory for the database files

Outputs (written under OUT_DIR):
  train_acc.txt                  - sorted list of accessions included
  train_phage_proteins.faa       - merged FAA, headers rewritten as >ACC__geneN
  train_phage_prot.dmnd          - DIAMOND DB
  train_phage_host_pair.csv      - subset of pair CSV (train only)
  build_db.log                   - summary
"""
from __future__ import annotations
import argparse, csv, os, shutil, subprocess, sys
from pathlib import Path

_DEFAULT_BASE = Path(os.environ.get("PHI_PROJECT_ROOT",
                     str(Path(__file__).resolve().parents[3])))

DIAMOND = Path(os.environ.get("DIAMOND_BIN",
               shutil.which("diamond") or "diamond"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--train-dir", type=Path, required=True,
                    help="Directory of train-set phage subdirs (each with proteins.faa)")
    ap.add_argument("--pair-csv",  type=Path, required=True,
                    help="Phage-host pair CSV (accession,host)")
    ap.add_argument("--out-dir",   type=Path, required=True,
                    help="Output directory for database files")
    ap.add_argument("--train-list", type=Path, default=None,
                    help="Text file listing train accessions (one per line). "
                         "If omitted, all subdirectories in --train-dir are used.")
    ap.add_argument("--faa-pattern", default="{acc}/proteins.faa",
                    help="FAA file path pattern relative to train-dir. "
                         "Use {acc} as placeholder. Default: '{acc}/proteins.faa'")
    ap.add_argument("--diamond",   type=Path, default=DIAMOND)
    ap.add_argument("--threads",   type=int,  default=8)
    ap.add_argument("--dry-run",   action="store_true",
                    help="Validate inputs and report counts; skip large file writes")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    log = open(args.out_dir / "build_db.log", "w") if not args.dry_run else sys.stdout
    def say(msg: str) -> None:
        print(msg, file=log, flush=True)
        if not args.dry_run:
            print(msg, flush=True)

    # 1. Enumerate train accessions
    if args.train_list:
        train_accs = sorted(
            l.strip() for l in args.train_list.read_text().splitlines() if l.strip()
        )
        say(f"[1] train-list: {args.train_list}")
    else:
        train_accs = sorted(p.name for p in args.train_dir.iterdir() if p.is_dir())
        say(f"[1] train-dir: {args.train_dir}")
    say(f"    acc directories: {len(train_accs)}")
    if not train_accs:
        say("[ERR] no accessions found")
        return 1

    # 2. Verify each acc has a FAA file
    faa_present, faa_missing = [], []
    for acc in train_accs:
        faa = args.train_dir / args.faa_pattern.format(acc=acc)
        (faa_present if faa.exists() else faa_missing).append(acc)
    say(f"[2] FAA check:")
    say(f"    present: {len(faa_present)} / {len(train_accs)}")
    if faa_missing:
        say(f"    missing ({len(faa_missing)}): first 5 → {faa_missing[:5]}")

    (args.out_dir / "train_acc.txt").write_text("\n".join(faa_present) + "\n")

    # 3. Subset host pair CSV
    host_of: dict[str, str] = {}
    with args.pair_csv.open() as f:
        for line in f:
            line = line.strip()
            if not line or "," not in line:
                continue
            acc, host = line.split(",", 1)
            host_of[acc] = host
    pair_subset = [(acc, host_of[acc]) for acc in faa_present if acc in host_of]
    pair_missing = [acc for acc in faa_present if acc not in host_of]
    say(f"[3] pair-csv: {args.pair_csv} ({len(host_of)} total rows)")
    say(f"    matched: {len(pair_subset)} / {len(faa_present)}; missing host: {len(pair_missing)}")

    pair_out = args.out_dir / "train_phage_host_pair.csv"
    with pair_out.open("w", newline="") as f:
        w = csv.writer(f)
        for acc, host in pair_subset:
            w.writerow([acc, host])
    say(f"    wrote {pair_out}")

    if args.dry_run:
        say("[dry_run] stopping before FAA merge / diamond makedb.")
        return 0

    # 4. Merge FAA (rewrite headers as >ACC__geneN)
    out_faa = args.out_dir / "train_phage_proteins.faa"
    n_phages, n_proteins = 0, 0
    with out_faa.open("w") as out:
        for acc in faa_present:
            faa = args.train_dir / args.faa_pattern.format(acc=acc)
            n_phages += 1
            gene_idx = 0
            for line in faa.read_text().splitlines():
                if line.startswith(">"):
                    gene_idx += 1
                    out.write(f">{acc}__gene{gene_idx}\n")
                    n_proteins += 1
                elif line.strip():
                    out.write(line + "\n")
    say(f"[4] merged FAA: {out_faa}")
    say(f"    phages   : {n_phages}")
    say(f"    proteins : {n_proteins}")

    # 5. diamond makedb
    out_dmnd = args.out_dir / "train_phage_prot.dmnd"
    cmd = [
        str(args.diamond), "makedb",
        "--in",  str(out_faa),
        "--db",  str(out_dmnd.with_suffix("")),
        "--threads", str(args.threads),
    ]
    say(f"[5] diamond makedb")
    say(f"    cmd: {' '.join(cmd)}")
    rc = subprocess.run(cmd, check=False).returncode
    if rc != 0:
        say(f"    [ERR] diamond exit {rc}")
        return rc
    say(f"    wrote {out_dmnd}")

    say("[OK] done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
