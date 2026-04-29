"""
prompt_v22G.py
==============
PHI-reason v22G prompt 模板定义（三段式结构）。

结构说明：
  ① SYSTEM_PROMPT         — 固定，每次请求都发（Ollama system 字段）
  ② USER_PREFIX_TEMPLATE  — 固定前缀，含宿主列表；跨 phage 共享 → KV cache 复用
  ③ USER_VARIABLE_TEMPLATE — 每个噬菌体不同，含 phage profile + 推断指令

用法示例:
  from prompt_v22G import SYSTEM_PROMPT, USER_PREFIX_TEMPLATE, USER_VARIABLE_TEMPLATE

  user_prefix = USER_PREFIX_TEMPLATE.format(n_hosts=223, host_list=host_list_text)
  user_msg    = user_prefix + USER_VARIABLE_TEMPLATE.format(phage_profile=profile_text)
"""

# ─────────────────────────────────────────────────────────────────────────────
# ① SYSTEM_PROMPT
#    固定，每次请求都发（Ollama /api/generate 的 system 字段）
# ─────────────────────────────────────────────────────────────────────────────
SYSTEM_PROMPT = """\
/no_think
You are a bacteriophage-host interaction prediction system.

Given a phage genome profile and a list of candidate bacterial hosts, predict the infection
probability for each host based on the provided data.

Respond as fast as possible. Output ONLY the JSON object below.
"""

# ─────────────────────────────────────────────────────────────────────────────
# ② USER_PREFIX_TEMPLATE
#    固定前缀，含宿主列表（Format C genus-grouped）
#    跨 phage 完全相同 → Ollama/llama.cpp 自动复用 KV cache
#    占位符: {n_hosts}, {host_list}
# ─────────────────────────────────────────────────────────────────────────────
USER_PREFIX_TEMPLATE = """\
=== {n_hosts} CANDIDATE HOSTS (family-grouped, Gram annotated) ===
Layout:
  === Gram-[negative|positive]: Family (N species) ===
    N. HostName | habitat | receptor | phage_families

Field meanings:
  • habitat         — ecological niche / host of isolation / notable lifestyle
  • receptor        — known surface features phages adsorb to (LPS O-antigen, capsule, pili, teichoic acid, etc.)
  • phage_families  — representative known phages + receptor summary for this host

'-' means field absent or unknown.

{host_list}

---

"""

# ─────────────────────────────────────────────────────────────────────────────
# ③ USER_VARIABLE_TEMPLATE
#    每个噬菌体不同，含 phage profile + 推断指令 + JSON 格式要求
#    占位符: {phage_profile}
# ─────────────────────────────────────────────────────────────────────────────
USER_VARIABLE_TEMPLATE = """\
=== PHAGE GENOME PROFILE ===
{phage_profile}

---

INSTRUCTIONS:
  ⚠️  ANTI-BIAS RULE: Do not favor common hosts (Escherichia/Salmonella/Klebsiella) without evidence.

  For EACH host, evaluate:
    1. TAIL/RBP match — strongest signal (check `← RBP matches (BLASTP): Genus (identity=X%, qcov=Y%)` lines)
    2. LYSIS compatibility — Gram type match (Gram- endolysin cannot lyse Gram+ cell walls)
    3. TAXONOMIC signals — BLASTN neighbor hosts, eggNOG `[Gram-/Gram+/Genus]` tags on individual genes

  Scoring guide:
    0.0      = incompatible cell envelope
    0.01-0.2 = compatible envelope, no receptor evidence
    0.2-0.5  = weak / indirect evidence
    0.5-0.7  = receptor class matches
    0.7-0.9  = specific receptor gene or genus-level match
    0.9-1.0  = near-certain (direct annotation or literature-grade)

Output ONLY valid JSON (no text before or after):
{{
  "gram_type_decision": "Gram- | Gram+ | Acid-fast | Archaea | unknown",
  "reasoning": "<3-4 sentence mechanistic summary: RBP name, inferred receptor, top genus rationale>",
  "predictions": [
    {{"rank": 1, "host": "<exact name from candidate list>", "score": <0.0-1.0>}},
    ...
  ]
}}

Rules:
- Host names must match exactly as shown in the candidate list above.
- Include exactly 30 entries, ordered by score descending.
- Hosts not listed implicitly score 0.0.
- After the closing `}}` stop immediately.
"""

TOP_K = 30  # JSON 中要求输出的 host 数量
