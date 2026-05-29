"""
prompt_template.py
==================
PHI-reason v22G prompt template definitions (three-part structure).

Structure:
  1. SYSTEM_PROMPT         — fixed, sent with every request (Ollama system field)
  2. USER_PREFIX_TEMPLATE  — fixed prefix containing host list; shared across phages for KV cache reuse
  3. USER_VARIABLE_TEMPLATE — varies per phage, contains phage profile + inference instructions

Usage example:
  from prompt_template import SYSTEM_PROMPT, USER_PREFIX_TEMPLATE, USER_VARIABLE_TEMPLATE

  user_prefix = USER_PREFIX_TEMPLATE.format(n_hosts=223, host_list=host_list_text)
  user_msg    = user_prefix + USER_VARIABLE_TEMPLATE.format(phage_profile=profile_text)
"""

# ─────────────────────────────────────────────────────────────────────────────
# 1. SYSTEM_PROMPT
#    Fixed; sent with every request (Ollama /api/generate system field)
# ─────────────────────────────────────────────────────────────────────────────
SYSTEM_PROMPT = """\
/no_think
You are a bacteriophage-host interaction prediction system.

Given a phage genome profile and a list of candidate bacterial hosts, predict the infection
probability for each host based on the provided data.

Respond as fast as possible. Output ONLY the JSON object below.
"""

# ─────────────────────────────────────────────────────────────────────────────
# 2. USER_PREFIX_TEMPLATE
#    Fixed prefix containing host list (Format C genus-grouped)
#    Identical across phages -> Ollama/llama.cpp automatically reuses KV cache
#    Placeholders: {n_hosts}, {host_list}
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
# 3. USER_VARIABLE_TEMPLATE
#    Varies per phage; contains phage profile + inference instructions + JSON format
#    Placeholder: {phage_profile}
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

TOP_K = 30  # Number of hosts required in JSON output
