#!/usr/bin/env python3
"""
gen_hostlist_v4G.py
===================
基于 Stage-2 host profiles (默认 v3_ncbi_R1 格式) 生成 v4G host list。

输入字段名（v3_ncbi_R1，docs/host_textProfile_gen/build_host_profiles.py 输出）：
  GRAM:, TAXONOMY:, ECOLOGY:, PRIMARY_RECEPTORS:, DEFENSE:, KNOWN_PHAGE_FAMILIES:

输出格式：
  === Gram-negative|Gram-positive|Archaea: Family (N) ===
  N. HostName | habitat | receptor | defense | phage_families

用法：
  python gen_hostlist_v4G.py [--host-prof DIR] [--out FILE]

历史备注：
  先前版本读取 ws/textGeneProfile_v4（字段 RECEPTOR / PHAGE_FAM），
  那是 v3_ncbi_R1 的 rename 版本，现已 abandon。此脚本直接读 v3_ncbi_R1
  字段，同时向后兼容旧 v4 字段。
"""
import argparse
import os
import re
from pathlib import Path
from collections import Counter

_DEFAULT_BASE = Path(os.environ.get("PHI_PROJECT_ROOT",
                     str(Path(__file__).resolve().parents[3])))
DEFAULT_HOST_PROF = _DEFAULT_BASE / "ws/textGeneProfile_v3_ncbi_R1"
DEFAULT_OUT_FILE  = Path(__file__).resolve().parent / "host_list_v4G.txt"

RECEP_EMPTY_MARKERS = ("see surface genes below", "surface genes below", "see below")

_GRAM_SHORT = {
    "GRAM-NEGATIVE": "Gram-negative",
    "GRAM-POSITIVE": "Gram-positive",
    "ARCHAEA":       "Archaea",
}

# 字段别名：(canonical_key, [可能的 prefix 列表，按优先级])
# v3_ncbi_R1 使用 PRIMARY_RECEPTORS / KNOWN_PHAGE_FAMILIES，
# archived v4 使用 RECEPTOR / PHAGE_FAM —— 两者都支持。
_FIELD_ALIASES = [
    ("gram",      ["GRAM:"]),
    ("taxonomy",  ["TAXONOMY:"]),
    ("ecology",   ["ECOLOGY:"]),
    ("receptor",  ["PRIMARY_RECEPTORS:", "RECEPTOR:"]),
    ("defense",   ["DEFENSE:"]),
    ("phage_fam", ["KNOWN_PHAGE_FAMILIES:", "PHAGE_FAM:"]),
]

def _parse_quick_profile_v4(md_text: str) -> dict:
    out = {
        "gram": "?", "taxonomy": "", "ecology": "",
        "receptor": "-", "defense": "-", "phage_fam": "-",
    }
    in_block = False
    for line in md_text.splitlines():
        s = line.strip()
        if s.startswith("## QUICK_PROFILE"):
            in_block = True; continue
        if in_block and s.startswith("##"):
            break
        if not in_block or not s:
            continue
        for key, prefixes in _FIELD_ALIASES:
            for prefix in prefixes:
                if s.startswith(prefix):
                    out[key] = s.split(":", 1)[1].strip()
                    break
    return out

def _infer_gram_group(gram_str: str) -> str:
    g = gram_str.lower()
    if "gram-positive" in g or "gram+" in g: return "GRAM-POSITIVE"
    if "gram-negative" in g or "gram-" in g: return "GRAM-NEGATIVE"
    if "archaea" in g:                        return "ARCHAEA"
    return "GRAM-POSITIVE"

def _last_family(taxonomy: str) -> str:
    parts = [p.strip() for p in taxonomy.split(">") if p.strip()]
    if not parts: return ""
    for p in reversed(parts):
        if p.endswith("aceae"): return p
    return parts[-1]

def _clean_receptor(r: str) -> str:
    s = r.strip()
    if not s: return "-"
    for m in RECEP_EMPTY_MARKERS:
        if m in s.lower(): return "-"
    return s

def _clean_defense(d: str) -> str:
    s = d.strip()
    if not s or s == "-": return "-"
    s = re.sub(r'^\d+\s+systems?\s*[—\-]\s*', '', s, flags=re.IGNORECASE)
    s = re.sub(r'\(\d+\)', '', s)
    parts = [p.strip() for p in re.split(r'[;,]', s) if p.strip()]
    return "; ".join(parts) if parts else "-"

def _clean_phagefam(p: str) -> str:
    s = p.strip()
    if not s: return "-"
    return s

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host-prof", type=Path, default=DEFAULT_HOST_PROF,
                    help=f"host profile 目录（默认 {DEFAULT_HOST_PROF}）")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT_FILE,
                    help=f"输出文件（默认 {DEFAULT_OUT_FILE}）")
    args = ap.parse_args()

    if not args.host_prof.exists():
        raise SystemExit(f"ERROR: host profile 目录不存在: {args.host_prof}")

    host_data = []

    for md_path in sorted(args.host_prof.glob("*.md")):
        name = md_path.stem
        text = md_path.read_text(encoding="utf-8", errors="replace")
        f = _parse_quick_profile_v4(text)

        habitat  = f["ecology"].strip() if f["ecology"] else "-"
        receptor = _clean_receptor(f["receptor"])
        defense  = _clean_defense(f["defense"])
        phagefam = _clean_phagefam(f["phage_fam"])
        gram_raw = f["gram"]
        family   = _last_family(f["taxonomy"]) or "?"

        host_data.append((name, habitat, receptor, defense, phagefam, gram_raw, family))

    bucket_order = ["GRAM-NEGATIVE", "GRAM-POSITIVE", "ARCHAEA"]
    entries = []
    for i, (name, habitat, receptor, defense, phagefam, gram_raw, family) in \
            enumerate(sorted(host_data, key=lambda x: x[0]), 1):
        gb = _infer_gram_group(gram_raw)
        entries.append((gb, family, i, name, habitat, receptor, defense, phagefam))

    entries.sort(key=lambda x: (
        bucket_order.index(x[0]) if x[0] in bucket_order else 99,
        x[1]
    ))

    fam_count = Counter((e[0], e[1]) for e in entries)

    lines = []
    cur_header = ""
    for (gb, family, idx, name, habitat, receptor, defense, phagefam) in entries:
        header = f"{_GRAM_SHORT.get(gb, gb)}: {family}"
        n_sp = fam_count[(gb, family)]

        if header != cur_header:
            if cur_header: lines.append("")
            lines.append(f"=== {header} ({n_sp}) ===")
            cur_header = header

        row = f"{idx}. {name} | {habitat} | {receptor} | {defense} | {phagefam}"
        lines.append(row)

    lines.append("")
    text = "\n".join(lines)
    args.out.write_text(text, encoding="utf-8")

    print(f"生成完成: {args.out}")
    print(f"  宿主数  : {len(host_data)}")
    print(f"  文件大小: {len(text)/1024:.1f} KB  (~{len(text)//4} tokens)")
    for line in text.splitlines()[:12]:
        print(" ", line)

if __name__ == "__main__":
    main()
