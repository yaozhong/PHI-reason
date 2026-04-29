#!/usr/bin/env python3
"""
01_build_db.py — Build CHERRY-1306 protein DIAMOND database (clean, no test leakage).

Inputs (read-only):
  TRAIN_DIR  - dir holding 1306 acc subdirectories (PHP retrain trainData)
  PROT_POOL  - dir holding all 1940 phages' .faa files (CHERRY_phage_1940.part_<ACC>.faa)
  PAIR_CSV   - phage1940_host_pair.csv (acc,host)

Outputs (written under OUT_DIR):
  trainData_cherry1306_acc.txt        - sorted list of acc actually included
  cherry1306_phage_proteins.faa       - merged FAA, headers rewritten as >ACC__geneN
  cherry1306_phage_prot.dmnd          - DIAMOND DB
  cherry1306_phage_host_pair.csv      - subset of phage1940_host_pair.csv (1306 only)
  build_db.log                        - summary
"""
from __future__ import annotations
import argparse, csv, os, shutil, subprocess, sys
from pathlib import Path

_DEFAULT_BASE = Path(os.environ.get("PHI_PROJECT_ROOT",
                     str(Path(__file__).resolve().parents[3])))

DEFAULT_TRAIN = Path(os.environ.get("CHERRY_TRAIN_DIR",
                     str(_DEFAULT_BASE / "data/Cherry_data/trainData_cherry1306")))
DEFAULT_POOL  = Path(os.environ.get("CHERRY_PROT_POOL",
                     str(_DEFAULT_BASE / "data/Cherry_data/phage/prot")))
DEFAULT_PAIR  = Path(os.environ.get("CHERRY_PAIR_CSV",
                     str(_DEFAULT_BASE / "data/Cherry_data/phage1940_host_pair.csv")))
DEFAULT_OUT   = _DEFAULT_BASE / "experiments/blastp_cherry1306/db"
DIAMOND       = Path(os.environ.get("DIAMOND_BIN",
                     shutil.which("diamond") or "diamond"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--train_dir", type=Path, default=DEFAULT_TRAIN)
    ap.add_argument("--prot_pool", type=Path, default=DEFAULT_POOL)
    ap.add_argument("--pair_csv",  type=Path, default=DEFAULT_PAIR)
    ap.add_argument("--out_dir",   type=Path, default=DEFAULT_OUT)
    ap.add_argument("--diamond",   type=Path, default=DIAMOND)
    ap.add_argument("--threads",   type=int,  default=8)
    ap.add_argument("--dry_run",   action="store_true",
                    help="Validate inputs and report counts; do not write large files")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    log = open(args.out_dir / "build_db.log", "w") if not args.dry_run else sys.stdout
    def say(msg: str) -> None:
        print(msg, file=log, flush=True)
        if not args.dry_run:
            print(msg, flush=True)

    # 1. Acc list from train_dir (1306 nominal, 1298 actual)
    train_accs = sorted(p.name for p in args.train_dir.iterdir() if p.is_dir())
    say(f"[1] train_dir: {args.train_dir}")
    say(f"    acc directories: {len(train_accs)}")
    if not train_accs:
        say("[ERR] no acc subdirectories found")
        return 1

    # 2. Verify each acc has a .faa in prot_pool
    faa_present, faa_missing = [], []
    for acc in train_accs:
        faa = args.prot_pool / f"CHERRY_phage_1940.part_{acc}.faa"
        (faa_present if faa.exists() else faa_missing).append(acc)
    say(f"[2] prot_pool: {args.prot_pool}")
    say(f"    .faa present: {len(faa_present)} / {len(train_accs)}")
    if faa_missing:
        say(f"    .faa missing ({len(faa_missing)}): first 5 → {faa_missing[:5]}")

    (args.out_dir / "trainData_cherry1306_acc.txt").write_text(
        "\n".join(faa_present) + "\n"
    )

    # 3. Subset host pair csv
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
    say(f"[3] pair_csv: {args.pair_csv} ({len(host_of)} total rows)")
    say(f"    matched: {len(pair_subset)} / {len(faa_present)}; missing host: {len(pair_missing)}")
    if pair_missing:
        say(f"    missing first 5 → {pair_missing[:5]}")

    pair_out = args.out_dir / "cherry1306_phage_host_pair.csv"
    with pair_out.open("w", newline="") as f:
        w = csv.writer(f)
        for acc, host in pair_subset:
            w.writerow([acc, host])
    say(f"    wrote {pair_out}")

    if args.dry_run:
        say("[dry_run] stopping before FAA merge / diamond makedb.")
        return 0

    # 4. Merge FAA (rewrite headers >ACC__<orig>)
    out_faa = args.out_dir / "cherry1306_phage_proteins.faa"
    n_phages, n_proteins = 0, 0
    with out_faa.open("w") as out:
        for acc in faa_present:
            faa = args.prot_pool / f"CHERRY_phage_1940.part_{acc}.faa"
            n_phages += 1
            for line in faa.read_text().splitlines():
                if line.startswith(">"):
                    out.write(f">{acc}__{line[1:]}\n")
                    n_proteins += 1
                elif line.strip():
                    out.write(line + "\n")
    say(f"[4] merged FAA: {out_faa}")
    say(f"    phages   : {n_phages}")
    say(f"    proteins : {n_proteins}")

    # 5. diamond makedb
    out_dmnd = args.out_dir / "cherry1306_phage_prot.dmnd"
    cmd = [
        str(args.diamond), "makedb",
        "--in",  str(out_faa),
        "--db",  str(out_dmnd.with_suffix("")),  # diamond appends .dmnd
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
