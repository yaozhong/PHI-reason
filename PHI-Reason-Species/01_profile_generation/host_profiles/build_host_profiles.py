#!/usr/bin/env python3
"""
build_host_profiles.py
=======================
Stage 2 of the PHI-Reason host profile pipeline.

Takes base host profiles (from build_host_profiles_base.py) and enriches them
with NCBI taxonomy, curated phage-host biology knowledge, and optionally PubMed
literature snippets to produce final text profiles for LLM-based phage-host
prediction.

Output format — two-layer progressive disclosure:
  ## QUICK_PROFILE   (~200 tokens)
      GRAM, TAXONOMY, ECOLOGY, PRIMARY_RECEPTORS, DEFENSE, KNOWN_PHAGE_FAMILIES,
      SPECIES_SPECIFICITY, ANNOTATION_WARNING, HOST_RANGE_BOUNDARY
  ## PHAGE_INTERACTION_DETAILS  (~1500 tokens)
      Receptor Biochemistry, Defense Systems, Known Infecting Phage Families,
      Species Specificity vs Closest Relatives, Cross-Infection Boundaries,
      Phage Recognition Signals

Post-processing applied (v3_ncbi_R1 equivalent):
  - NCBI TaxID strings removed from output (not useful for LLM reasoning)
  - Key Literature section omitted (titles too generic; adds noise)

Curated knowledge embedded in this script (no external files required):
  GRAM_MAP            — Gram type for 114 genera
  HABITAT_MAP         — Ecological context for 206 species
  KNOWN_RECEPTORS     — Experimentally confirmed receptors for 78 species
  PHAGE_FAMILIES      — Known infecting phage groups for 70 genera/species
  INTRASPECIES_DISTINCTION — Discriminating features for 50+ closely-related species
  PHAGE_RECEPTOR_NOTES     — Annotation warnings to prevent LLM misinterpretation
  CROSS_INFECTION_ALERT    — Host-range boundary notes
  PHAGE_RECOGNITION_GUIDE  — If-phage-has-X → indicates/excludes-host signals

Prerequisites:
  - Stage 1 complete: base profiles in --base-profiles dir
    (run build_host_profiles_base.py first)
  - Internet access for NCBI eutils (taxonomy fetch)
  - Optional: NCBI API key for higher rate limits (--ncbi-api-key)

Usage:
    python build_host_profiles.py \\
        --base-profiles /path/to/host_base_profiles \\
        --out-dir       /path/to/host_profiles_v3_ncbi_R1 \\
        [--host-list    hosts.json] \\
        [--ncbi-api-key YOUR_KEY] \\
        [--threads      2] \\
        [--skip-pubmed] \\
        [--force] \\
        [--dry-run]
"""

from __future__ import annotations

import argparse
import json
import re
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# ── NCBI eutils settings ───────────────────────────────────────────────────────
NCBI_EMAIL = "phi-reason-build@local"
NCBI_TOOL  = "phi_reason_host_profile"
RATE_DELAY = 0.35   # seconds between requests without API key (≤3/s)
RATE_DELAY_KEY = 0.11  # seconds with API key (≤10/s)

# ── Intra-genus / intra-family discriminating features ────────────────────────
INTRASPECIES_DISTINCTION = {
    "Pseudomonas_aeruginosa":
        "DISTINCT from P. fluorescens/syringae/putida: (1) type IV pili (PilA/PilQ/PilB) — primary receptor "
        "for ΦKZ, LUZ19, PA phages; (2) PAO1/PA14-specific LPS O-antigen (O5 most common) targeted by LPS-specific "
        "phages; (3) alginate capsule in mucoid strains; (4) clinical lung/wound isolate context. "
        "P. fluorescens lacks type IV pili as phage receptor and has different O-antigen chemistry.",
    "Pseudomonas_fluorescens":
        "DISTINCT from P. aeruginosa: no type IV pili as major phage receptor; soil/food spoilage environment "
        "(not clinical); different LPS O-antigen structure. Phages specific to P. fluorescens "
        "often use flagella or outer membrane proteins as receptors rather than pili.",
    "Pseudomonas_putida":
        "DISTINCT from P. aeruginosa/fluorescens: soil/rhizosphere specialist; no type IV pili phage receptor; "
        "rarely causes human infection; flagella-mediated phage adsorption predominant.",
    "Pseudomonas_syringae":
        "DISTINCT from P. aeruginosa: plant pathogen (phyllosphere), not animal/clinical host; "
        "type III secretion system for plant cell entry; LPS with plant host-adapted O-antigen.",
    "Klebsiella_pneumoniae":
        "DISTINCT from K. oxytoca/aerogenes: (1) highly diverse K-antigen capsule (K1-K147) — phage "
        "tail spikes with PECTATE LYASE / polysaccharide depolymerase domains target K-antigen "
        "(NOT plant pectin); (2) hypervirulent strains with hypermucoviscous capsule; (3) primary "
        "hospital pathogen. Phage 'pectate lyase' annotation = K-antigen capsule depolymerase.",
    "Klebsiella_oxytoca":
        "DISTINCT from K. pneumoniae: produces oxazoline (causes antibiotic-associated colitis); "
        "K-antigen capsule but different serotype spectrum from K. pneumoniae; OmpK35/OmpK36 porins "
        "similar but different capsule K-types targeted by specific phage tail spikes.",
    "Klebsiella_aerogenes":
        "DISTINCT from K. pneumoniae/oxytoca: formerly Enterobacter aerogenes; less K-antigen diversity; "
        "more environmental/opportunistic; shares Enterobacteriaceae OmpC/OmpF porins.",
    "Escherichia_coli":
        "DISTINCT from Shigella/Citrobacter/Salmonella: (1) F-pilus present (receptor for ssRNA/ssDNA phages: "
        "MS2, Qβ, M13); (2) diverse O-antigen (O1–O187); (3) smooth vs. rough LPS affects T4 adsorption; "
        "(4) K-antigen capsule in some strains. Tail fibers targeting OmpC/OmpF/LamB are E. coli-specific.",
    "Salmonella_enterica":
        "DISTINCT from E. coli/Citrobacter: (1) P22-like phages use O-antigen LPS (Salmonella-specific "
        "O-antigen serotypes); (2) Vi-antigen capsule in typhi serovar; (3) no F-pilus (no ssRNA phage "
        "susceptibility); (4) host range often determined by specific LPS O-antigen chemistry.",
    "Citrobacter_freundii":
        "DISTINCT from E. coli/Salmonella: Citrobacter-specific O-antigen LPS; intimin-like outer membrane "
        "adhesins absent; phage tail fibers often annotated as 'Citrobacter-specific' in PHROG due to "
        "divergent receptor-binding domain. Shares Enterobacteriaceae porins (OmpC/OmpF) with E. coli.",
    "Citrobacter_rodentium":
        "DISTINCT from C. freundii: mouse colonic pathogen; intimin (eae) outer membrane adhesin; "
        "rodent-adapted; smaller host range for phages than C. freundii.",
    "Shigella_flexneri":
        "DISTINCT from E. coli/S. sonnei: (1) unique O-antigen (type 1-6 LPS) absent in E. coli; "
        "(2) no flagella (non-motile); (3) IcsA/VirG surface protein; (4) intracellular pathogen. "
        "Phages specific to S. flexneri target serotype-specific O-antigen.",
    "Shigella_sonnei":
        "DISTINCT from S. flexneri: single LPS O-antigen form I/II; less O-antigen diversity; "
        "frequently converts between rough and smooth phenotype affecting phage susceptibility.",
    "Shigella_boydii":
        "DISTINCT from S. flexneri/sonnei: O-antigen serotype C-group; different phage susceptibility "
        "profile; Shigella-specific but distinct LPS from other Shigella species.",
    "Staphylococcus_aureus":
        "DISTINCT from S. epidermidis/haemolyticus: (1) WTA with GlcNAc modification on ribitol phosphate "
        "backbone — primary Kayvirus receptor; (2) SpA (protein A) on surface; (3) clinical pathogen "
        "(MRSA/MSSA); (4) Protein A blocks IgG, relevant to phage-opsonization interactions. "
        "DISTINCT from Geobacillus kaustophilus: Staphylococcus phages are MESOPHILIC (37°C); Geobacillus "
        "requires thermophilic phages active at 55-70°C. Non-thermophilic Bacillales phages → Staphylococcus, "
        "NOT Geobacillus. Podovirus, Siphovirus, and Myovirus types all infect Staphylococcus.",
    "Staphylococcus_epidermidis":
        "DISTINCT from S. aureus: (1) WTA with GlcNAc on glycerol phosphate (not ribitol) backbone — "
        "different receptor chemistry from S. aureus; (2) biofilm-former (PNAG/PIA); (3) commensal/device "
        "pathogen; Kayvirus phages from S. aureus typically do NOT infect S. epidermidis.",
    "Mycolicibacterium_smegmatis":
        "DISTINCT from Mycobacterium avium/abscessus: (1) rapid grower (3-5 day colonies); (2) non-pathogenic "
        "model organism; (3) glycopeptidolipids (GPL) absent in M. avium — GPL mediate phage adsorption in "
        "smegmatis; (4) ~7,000+ isolated phages spanning Clusters A-Z use mycolic acid cell wall binding. "
        "Phage tropism in mycobacteria is driven by cell-wall glycolipid chemistry, NOT protein receptors.",
    "Mycobacterium_avium":
        "DISTINCT from Mycolicibacterium smegmatis: (1) slow grower (weeks); (2) pathogen (MAC complex); "
        "(3) glycopeptidolipid (GPL) surface — GPL serotype determines phage host range within MAC; "
        "(4) much fewer phages isolated compared to M. smegmatis model system.",
    "Mycobacteroides_abscessus":
        "DISTINCT from Mycolicibacterium smegmatis: (1) rapid grower but pathogen; (2) rough vs. smooth "
        "morphotype — smooth has GPL (susceptible to more phages), rough lacks GPL (restricted host range); "
        "(3) CF lung pathogen; (4) acquired resistance through GPL loss is key phage evasion mechanism.",
    "Lactococcus_lactis":
        "DISTINCT from Streptococcus thermophilus/Lactobacillus: (1) phage infection protein (Pip) on "
        "cell surface — primary receptor for 936/c2/P335 phage groups; (2) CWPS (cell-wall polysaccharide) "
        "is secondary receptor determining phage species specificity; (3) dairy fermentation starter.",
    "Streptococcus_thermophilus":
        "DISTINCT from L. lactis/S. mutans: (1) yogurt/dairy starter; (2) phage receptor = rhamnosyl "
        "cell-wall polysaccharide (not Pip); (3) Sfi phages (Siphoviridae) dominant; (4) thermophile "
        "(42°C optimum) vs. mesophile L. lactis.",
    "Bacillus_subtilis":
        "DISTINCT from B. anthracis/cereus/thuringiensis: (1) non-pathogenic model organism; (2) WTA "
        "(glucosylated poly-ribitol phosphate) as receptor; (3) SPO1/phi29-type phages; (4) sporulation "
        "model — phages can specifically target vegetative vs. sporulating cells.",
    "Bacillus_thuringiensis":
        "DISTINCT from B. subtilis/cereus: (1) produces Cry toxin insecticidal crystal proteins; "
        "(2) biocontrol agent; (3) shares WTA receptor with B. cereus group; (4) phages often "
        "cross-infect between B. thuringiensis and B. cereus due to similar surface chemistry.",
    "Vibrio_cholerae":
        "DISTINCT from V. parahaemolyticus/harveyi: (1) ICP1/ICP2/ICP3 phages are O1 El Tor-specific "
        "using O1 O-antigen as receptor; (2) MSHA pili (mannose-sensitive hemagglutinin pili) = receptor "
        "for some phages; (3) O139 strains use capsule; (4) cholera toxin phage (CTXφ) uses TCP pili. "
        "ICP phages do NOT infect other Vibrio species.",
    "Streptomyces_coelicolor":
        "DISTINCT from S. griseus/lividans: (1) spore coat-specific phage adsorption; (2) ΦC31 integrase "
        "model (integrates into attB site); (3) produces actinorhodin/undecylprodigiosin pigments as markers; "
        "(4) SapB aerial hyphae protein involved in phage access to spore surface.",
    "Gordonia_terrae":
        "DISTINCT from G. rubripertincta/sputi: (1) soil hydrocarbon degrader; (2) mycolic acid cell "
        "wall shared with all Gordonia — phage tropism determined by precise cell-wall glycolipid chemistry; "
        "(3) Cluster G phages from SEA-PHAGES often species-specific within Gordonia. "
        "KEY SIGNAL INTERPRETATION: 'Nocardiaceae', 'Nocardia', 'Rhodococcus', 'Corynebacterium (Actinobacteria)' "
        "PHROG annotations in phage = Mycobacteriales (Gordonia/Rhodococcus/Nocardia branch) indicator — these "
        "ARE the mycolic acid cell wall signals. You do NOT need an explicit 'mycolic acid binding' annotation; "
        "Nocardiaceae family membership IS the mycolic acid signal. DISTINCT from Cutibacterium: "
        "Propionibacteriales (skin) have NO mycolic acid cell wall. If phage has Nocardiaceae/Actinobacteria "
        "PHROG hits → Gordonia_terrae, NOT Cutibacterium.",
    "Pectobacterium_carotovorum":
        "DISTINCT from Salmonella/E. coli: (1) plant soft-rot pathogen (NOT mammalian gut pathogen); "
        "(2) LPS O-antigen with plant-adapted structure; (3) pectate lyase IS a virulence factor here "
        "(degrades plant pectin) — phage tail spikes with pectate lyase target Pectobacterium O-antigen, "
        "NOT mammalian Enterobacteriaceae. Phages from Salmonella/E. coli do NOT infect Pectobacterium.",
    "Pectobacterium_atrosepticum":
        "DISTINCT from Salmonella/Citrobacter: plant pathogen (potato blackleg); pectate lyase is "
        "a Pectobacterium virulence factor (pectin degradation in plant tissue). Phage tail spikes "
        "with pectate lyase domain = targeting Pectobacterium O-antigen polysaccharide. "
        "NOT cross-infectious with mammalian Enterobacteriaceae like Salmonella or Citrobacter.",
    "Acinetobacter_baumannii":
        "DISTINCT from A. johnsonii/soli: (1) primary nosocomial pathogen (ESKAPE); (2) OmpA is major "
        "phage receptor; (3) K-locus capsular polysaccharide (CPS) — phage tail fibers target "
        "K-type-specific CPS depolymerase; (4) carbapenem resistance common; (5) dry-environment survivor.",
    "Burkholderia_cenocepacia":
        "DISTINCT from B. pseudomallei/thailandensis: (1) CF lung pathogen (Burkholderia cepacia complex); "
        "(2) cable pili (unique to epidemic ET12 lineage) — receptor for some phages; (3) LPS with "
        "Burkholderia-specific O-antigen; (4) NOT the same as B. pseudomallei (BSL-3).",
    "Synechococcus_sp._WH_8102":
        "DISTINCT from Klebsiella/Enterobacteriaceae: CYANOBACTERIUM (marine photosynthetic); "
        "NO peptidoglycan-based surface like Enterobacteriaceae; phage receptors = cyanobacterial outer "
        "membrane proteins (OmpA-like); ONLY infected by cyanophages (S-PM2, Syn5 type). "
        "Completely different phage host range from Enterobacteriaceae.",
    "Microcystis_aeruginosa":
        "DISTINCT from Escherichia/Enterobacteriaceae: CYANOBACTERIUM (freshwater bloom-forming); "
        "gas vesicle proteins on surface; cyanophage-specific receptors; produces microcystin toxin; "
        "COMPLETELY different phage host range from E. coli/Shigella — no cross-infection possible.",
    "Flavobacterium_psychrophilum":
        "DISTINCT from Morganella/Enterobacteriaceae: (1) Bacteroidota phylum (NOT Proteobacteria); "
        "(2) psychrophile (optimal 15-18°C); (3) gliding motility; (4) freshwater fish pathogen; "
        "(5) cell envelope: Bacteroidota outer membrane with SusC/SusD nutrient transporters, NOT "
        "Enterobacteriaceae OmpC/OmpF. Phages infecting Flavobacterium do NOT infect Morganella.",
    "Cutibacterium_acnes":
        "DISTINCT from Gordonia/Mycolicibacterium/Rhodococcus: (1) Propionibacteriales order — NO mycolic acid "
        "cell wall (unlike Gordonia/Mycobacteriales); (2) skin commensal/pathogen (sebaceous follicles), NOT "
        "soil/environmental; (3) ferments propionate (anaerobic/microaerophilic); (4) phages: PAD20/PHL112M "
        "(Siphoviridae) — Propionibacteriales-specific, do NOT cross-infect Gordonia or Mycolicibacterium. "
        "Actinobacteria taxonomic markers in phage genome do NOT imply Gordonia/Mycobacteriales host.",
    "Lactobacillus_fermentum":
        "DISTINCT from Gordonia/Actinomycetota: (1) Bacillota phylum (Firmicutes), NOT Actinomycetota; "
        "(2) lactic acid fermentation (dairy/gut), NOT mycolic acid metabolism; (3) no mycolic acid cell wall; "
        "Gordonia phages (Cluster G, mycolic acid-binding) do NOT infect Lactobacillus. "
        "DISTINCT from Streptococcus: no Pip phage receptor; fermentum = heterofermentative.",
    "Pseudoalteromonas_marina":
        "DISTINCT from Proteus/Enterobacteriaceae: (1) marine Gammaproteobacteria (Alteromonadales, NOT "
        "Enterobacteriales); (2) Pseudoalteromonadaceae (completely different family from Morganellaceae); "
        "(3) marine/cold-adapted organism; phages are marine-adapted (PM2-like, cold-active Podoviridae) — "
        "NOT clinical enteric phages. Proteus phages target flagella-specific swarming LPS; Pseudoalteromonas "
        "phages use marine OMP and LPS in cold-water ecosystems.",
    "Brevibacillus_laterosporus":
        "DISTINCT from Staphylococcus aureus: (1) Paenibacillaceae (Bacillales order), NOT Staphylococcaceae; "
        "(2) spore-forming; insect pathogen; soil organism — NOT skin commensal; (3) phage PHROG annotations "
        "referencing Paenibacillaceae or Bacillus = Brevibacillus indicator, NOT Staphylococcus indicator; "
        "(4) WTA structure differs from Staphylococcus aureus GlcNAc-ribitol WTA. Bacillus/Paenibacillus "
        "PHROG hits → Brevibacillus, NOT Staphylococcus.",
    "Brucella_inopinata":
        "DISTINCT from Sinorhizobium/Rhizobiaceae: (1) Brucellaceae (intracellular animal pathogen); "
        "(2) Brucella phages (BTP1, Tb) recognize smooth LPS (unique S-LPS); (3) Brucellaceae PHROG hits "
        "are a STRONGER indicator than Rhizobiaceae hits for choosing Brucella vs Sinorhizobium; "
        "(4) no nitrogen-fixing capability (unlike Sinorhizobium); zoonotic pathogen, NOT soil nitrogen-fixer.",
}

# ── Phage receptor annotation notes ───────────────────────────────────────────
PHAGE_RECEPTOR_NOTES = {
    "Klebsiella_pneumoniae":
        "NOTE: 'Pectate lyase' or 'polysaccharide lyase' domains in phage TAIL SPIKES targeting "
        "Klebsiella = CAPSULE (K-antigen) DEPOLYMERASES, not plant pectin enzymes. Klebsiella "
        "K-antigen is a polysaccharide structurally similar to pectin — phage enzymes evolved to "
        "degrade it for adsorption. This is the primary infection mechanism for KP phages.",
    "Klebsiella_oxytoca":
        "NOTE: Phage 'pectate lyase' / polysaccharide lyase tail spikes target Klebsiella K-antigen "
        "capsular polysaccharide. This enzymatic tail spike degrades the capsule to allow "
        "phage adsorption to the underlying outer membrane.",
    "Pectobacterium_carotovorum":
        "NOTE: Pectobacterium IS a plant pathogen and pectate lyase is its OWN virulence factor. "
        "However, phage tail spikes with pectate lyase domains targeting Pectobacterium attack its "
        "LPS O-antigen polysaccharide (not plant pectin). Cross-infection from mammalian "
        "Enterobacteriaceae phages to Pectobacterium is extremely rare. "
        "CRITICAL: 'Peptidoglycan catabolic process' in phage tail fiber ≠ pectate lyase. "
        "Peptidoglycan hydrolysis (muramic acid) targets cell wall of ALL bacteria — it does NOT "
        "indicate Pectobacterium specifically. Pectobacterium requires pectate LYASE/polysaccharide lyase "
        "(plant pectin-related enzyme). Do NOT use 'peptidoglycan catabolic' as Pectobacterium signal.",
    "Pectobacterium_atrosepticum":
        "NOTE: Like P. carotovorum — pectate lyase in phage context = O-antigen depolymerase targeting "
        "Pectobacterium specifically. Salmonella/E. coli phages do NOT cross-infect Pectobacterium "
        "despite sharing Enterobacteriaceae family membership. "
        "CRITICAL: 'Peptidoglycan catabolic process' in tail fiber ≠ pectate lyase; do NOT use this as "
        "Pectobacterium indicator. Only pectate LYASE/polysaccharide lyase → Pectobacterium. "
        "CRITICAL: Pectobacterium annotations in HEAD proteins (terminase, portal, capsid, regulatory genes) "
        "are NOT host-specific — head proteins conserved across Enterobacteriaceae. If TAIL proteins annotate "
        "as Rahnella/Enterobacter/Salmonella, the host is mammalian Enterobacteriaceae, NOT Pectobacterium.",
    "Klebsiella_aerogenes":
        "NOTE: 'NosD' (copper-binding nosD gene) annotated on phage tail fiber = bacterial gene "
        "HGT artifact in phage genome — NOT evidence that the phage targets NosD-expressing hosts "
        "(like Pasteurella). This annotation does NOT indicate Pasteurella/Mannheimia tropism. "
        "Ignore NosD domain when assessing host range; use K-antigen/LPS O-antigen as the relevant receptor.",
    "Salmonella_enterica":
        "NOTE: 'N-acetylmuramoyl-L-alanine amidase' (muramidase) in phage tail = PEPTIDOGLYCAN hydrolase "
        "used to digest cell wall during infection — NOT a K-antigen capsule depolymerase. Phages with "
        "amidase tails target LPS/peptidoglycan of Enterobacteriaceae (Salmonella, E. coli), NOT "
        "Klebsiella K-antigen capsule (which requires polysaccharide depolymerase/lyase). "
        "'Concanavalin A-like lectin/glucanase' domain in phage gene = carbohydrate-binding for LPS O-antigen "
        "targeting — present in broad-host Myoviridae (Felix-O1 type) infecting Salmonella. "
        "Felix-O1 (Myoviridae) is the classic broad-host Salmonella phage; SP6-like infect Salmonella specifically.",
    "Escherichia_coli":
        "NOTE: T4-like Myoviridae have an INTERNAL lysozyme (gpe5/e-type) and tail-associated "
        "muramidase domains — these are Gram-NEGATIVE host markers (E. coli, Shigella, Salmonella), "
        "NOT Gram-positive indicators. Muramoyl/N-acetylmuramoyl domains in phage tail proteins targeting "
        "E. coli are used to digest the thin Gram-negative peptidoglycan layer during injection.",
}

# ── Cross-infection alert ──────────────────────────────────────────────────────
CROSS_INFECTION_ALERT = {
    "Salmonella_enterica":
        "CROSS-INFECTION BOUNDARY: Salmonella phages (esp. P22-like) use Salmonella-specific LPS O-antigen. "
        "Cross-infection to E. coli possible for some broad-host phages (LPS-independent entry via BtuB/FhuA). "
        "Cross-infection to Pectobacterium/plant pathogens: essentially zero despite Enterobacteriaceae membership.",
    "Escherichia_coli":
        "CROSS-INFECTION BOUNDARY: E. coli phages may cross-infect Shigella (nearly identical LPS/porins) "
        "and sometimes Citrobacter (similar OmpC/OmpF). Cross-infection to Pseudomonas/Vibrio/plant pathogens: "
        "extremely rare. F-pilus phages (MS2, M13) are E. coli specific — Salmonella lacks F-pilus.",
    "Pseudomonas_aeruginosa":
        "CROSS-INFECTION BOUNDARY: Pa-specific phages (ΦKZ, PaP1) use Pa-specific LPS O5 or type IV pili. "
        "Cross-infection to P. fluorescens: possible for some broad-host Myoviridae but NOT for pili-specific "
        "phages. Cross-infection to Vibrio/Enterobacteriaceae: essentially zero.",
    "Klebsiella_pneumoniae":
        "CROSS-INFECTION BOUNDARY: KP phages use K-antigen capsule (highly specific to K-type serotype). "
        "Cross-infection within Klebsiella genus (Kp/Ko/Ka) depends on capsule K-type similarity. "
        "Cross-infection to E. coli/Salmonella via OmpC/OmpF is possible for broad-host phages but "
        "capsule depolymerase-based phages are Klebsiella-specific.",
    "Flavobacterium_psychrophilum":
        "CROSS-INFECTION BOUNDARY: Flavobacterium is Bacteroidota — completely different outer membrane "
        "architecture from Proteobacteria (Morganella, Enterobacteriaceae). No cross-infection possible "
        "between Flavobacterium phages and Enterobacteriaceae. Psychrophile (15°C) vs. mesophile context "
        "also eliminates cross-infection likelihood.",
    "Synechococcus_sp._WH_8102":
        "CROSS-INFECTION BOUNDARY: Cyanobacterium — entirely different cell envelope from Enterobacteriaceae. "
        "Only cyanophages infect Synechococcus. No cross-infection with any Proteobacteria phage.",
    "Microcystis_aeruginosa":
        "CROSS-INFECTION BOUNDARY: Cyanobacterium — no cross-infection with Proteobacteria/Firmicutes. "
        "Freshwater bloom organism; cyanophage-specific receptors only.",
    "Mycolicibacterium_smegmatis":
        "CROSS-INFECTION BOUNDARY: Actinomycetota with mycolic acid cell wall — completely different "
        "from Firmicutes (Clostridium, Lactobacillus). No cross-infection between mycobacterium phages "
        "and Gram-positive Firmicutes. Cluster assignment determines cross-infection within mycobacteria.",
    "Pseudoalteromonas_marina":
        "CROSS-INFECTION BOUNDARY: Marine Alteromonadales — completely different ecology and receptor "
        "biochemistry from clinical Enterobacteriaceae (Proteus, E. coli, Salmonella). Marine Podoviridae "
        "and PM2-like phages do NOT infect Proteus or other enteric bacteria. Proteus swarming flagella "
        "receptors are NOT present in Pseudoalteromonas.",
    "Cutibacterium_acnes":
        "CROSS-INFECTION BOUNDARY: Propionibacteriales (skin) — distinct from Gordonia/Mycobacteriales. "
        "To distinguish: if phage PHROG annotations reference Propionibacterium/Cutibacterium-specific genes → Cutibacterium. "
        "If TAIL proteins reference Rhodococcus/Nocardia/Gordonia/Mycolicibacterium → Gordonia/Mycobacteriales, NOT Cutibacterium. "
        "CRITICAL: Rhodococcus/Nocardia annotation in HOLIN or ENDOLYSIN (lysis cassette [LOW] tag) = cross-annotation artifact; does NOT indicate Gordonia. Only TAIL protein genus annotations determine host range. "
        "Mycolic acid cell wall binding is Gordonia/Mycolicibacterium-specific; Cutibacterium lacks mycolic acids.",
    "Gordonia_terrae":
        "CROSS-INFECTION BOUNDARY: Mycobacteriales with mycolic acid cell wall — phages binding mycolic acid "
        "cell wall target Gordonia/Mycolicibacterium/Rhodococcus, NOT Cutibacterium (Propionibacteriales). "
        "Cluster G phages are Gordonia-specific. If phage encodes mycolic acid-binding proteins → Gordonia, not skin pathogens.",
    "Proteus_mirabilis":
        "CROSS-INFECTION BOUNDARY: Clinical Enterobacteriaceae — phages targeting Proteus use swarming-specific "
        "flagellar LPS (unique to Proteus). Cross-infection to marine Gammaproteobacteria (Pseudoalteromonas, "
        "Alteromonas) is essentially zero; marine organisms lack Proteus-specific flagellar antigens.",
    "Geobacillus_kaustophilus":
        "CROSS-INFECTION BOUNDARY: Thermophile (optimal 55-70°C) — ONLY thermostable phages infect Geobacillus. "
        "Mesophilic Bacillales phages (Staphylococcus, Listeria, Bacillus subtilis, Paenibacillus) do NOT "
        "infect Geobacillus. Known Geobacillus phages: ΦGP1 (Siphoviridae, thermophile-adapted). "
        "If phage has no thermophile adaptation → NOT Geobacillus; prefer mesophilic Gram-positive host.",
    "Vibrio_splendidus":
        "CROSS-INFECTION BOUNDARY: Marine Vibrionales — V. splendidus phages are typically small Myoviridae "
        "or Podoviridae (ICP-type, VP882-like), NOT giant Myoviridae (>100 ORF). Giant Myoviruses with Rz-like "
        "spanin (ΦKZ/LMA2-like, >150 ORF) specifically target P. aeruginosa, NOT Vibrio. Cross-infection "
        "from Pseudomonas-adapted phages to Vibrio: essentially zero (different LPS O-antigen, no type IV pili in Vibrio).",
    "Xanthomonas_vesicatoria":
        "CROSS-INFECTION BOUNDARY: Plant-pathogenic Xanthomonadaceae — completely different ecology and "
        "surface chemistry from mammalian Enterobacteriaceae. T4-like giant Myoviridae (PinA, gp19, gp15-like) "
        "infect Enterobacteriaceae (E. coli, Citrobacter, Salmonella), NOT Xanthomonas. No cross-infection "
        "between T4/T7-type Enterobacteriaceae phages and Xanthomonas plant pathogens.",
    "Morganella_morganii":
        "CROSS-INFECTION BOUNDARY: Clinical Enterobacteriaceae (Proteobacteria). Vibrio annotation in phage "
        "genes indicates marine/aquatic ecology, NOT clinical Morganella. Bacteroidota (Flavobacterium) phages "
        "CANNOT infect Morganella (completely different phylum). If phage has Vibrio/marine signals or "
        "Verrucomicrobiae/Bacteroidota annotations → NOT Morganella (clinical Enterobacteriaceae).",
    "Klebsiella_oxytoca":
        "CROSS-INFECTION BOUNDARY: Klebsiella phages require K-antigen capsule-targeting enzyme "
        "(polysaccharide depolymerase/lyase) in tail to penetrate the capsule. Podoviridae/Siphoviridae "
        "WITHOUT a capsule depolymerase/pectate lyase cannot efficiently infect encapsulated Klebsiella. "
        "If phage lacks capsule depolymerase → prefer non-encapsulated Enterobacteriaceae (Salmonella, "
        "Citrobacter, E. coli) over Klebsiella.",
}

# ── Phage-to-host recognition signal guide ────────────────────────────────────
PHAGE_RECOGNITION_GUIDE = {
    "Escherichia_coli": [
        "✓ INDICATES: T4 INTERNAL LYSOZYME (gpe5/e-type) or baseplate muramidase genes → E. coli/Shigella target",
        "✓ INDICATES: Lambda (λ) phage-like genes (cI repressor, N antitermination) → E. coli/Shigella",
        "✓ INDICATES: F-pilus receptor (ssRNA phages MS2/Qβ, ssDNA M13/fd) → E. coli ONLY",
        "✓ INDICATES: LamB/OmpC/OmpF tail fiber → Enterobacteriaceae, primary E. coli",
        "✓ INDICATES: muramoyl/N-acetylmuramoyl domain in tail → Gram-negative thin PG (E. coli/Salmonella/Shigella)",
        "✗ EXCLUDES: pectate lyase / depolymerase tail spike → not E. coli primary (Klebsiella capsule target)",
        "✗ EXCLUDES: Giant Myoviridae (>300 ORF) with Alphaproteobacteria PHROG hits → NOT E. coli (may infect Sinorhizobium/Rhizobiales instead)",
    ],
    "Salmonella_enterica": [
        "✓ INDICATES: P22-like genes (P22 tailspike, eject protein) → Salmonella O-antigen specific",
        "✓ INDICATES: SP6-like (LPS O-antigen), epsilon15-like → Salmonella",
        "✓ INDICATES: Felix-O1-like (Myoviridae, broad Salmonella host) → Salmonella",
        "✓ INDICATES: Concanavalin A-like lectin on phage tail → carbohydrate LPS binding (Felix-O1 type, Salmonella Myoviridae)",
        "✗ EXCLUDES: F-pilus targeting → NOT Salmonella (Salmonella has no F-pilus)",
        "✗ EXCLUDES: pectate lyase tail spike → NOT Salmonella (that's Klebsiella capsule / Pectobacterium)",
    ],
    "Klebsiella_pneumoniae": [
        "✓ INDICATES: pectate lyase / polysaccharide depolymerase / capsule depolymerase tail spike → K-antigen targeting (KP primary)",
        "✓ INDICATES: Drulisvirus / Sugarlandvirus gene signatures → Klebsiella",
        "✗ EXCLUDES: NosD (copper-binding) on tail fiber → ARTIFACT, not receptor signal",
    ],
    "Klebsiella_aerogenes": [
        "✓ INDICATES: KP-like phage genes with K-antigen targeting → Klebsiella aerogenes",
        "✗ EXCLUDES: NosD / periplasmic copper-binding on tail fiber → ARTIFACT, not Pasteurella indicator",
    ],
    "Pseudomonas_aeruginosa": [
        "✓ INDICATES: ΦKZ-like giant phage genes → P. aeruginosa",
        "✓ INDICATES: LUZ19/LKD16-like (Autographiviridae) → P. aeruginosa",
        "✓ INDICATES: type IV pili binding domain in tail fiber → P. aeruginosa (P. fluorescens lacks pili receptor)",
        "✓ INDICATES: PAO1/PA14-specific LPS O5 antigen → P. aeruginosa",
        "✗ EXCLUDES: Propionibacterium/Cutibacterium PHROG → NOT P. aeruginosa",
    ],
    "Staphylococcus_aureus": [
        "✓ INDICATES: WTA-binding domain (GlcNAc-ribitol phosphate) → Staphylococcus aureus specific",
        "✓ INDICATES: Kayvirus/K-like phage genes → S. aureus",
        "✗ EXCLUDES: no WTA binding → less likely S. aureus (WTA is primary receptor)",
    ],
    "Mycolicibacterium_smegmatis": [
        "✓ INDICATES: Cluster A (D29-like) to Cluster Z mycobacteriophage genes → Mycolicibacterium/Mycobacterium specifically",
        "✓ INDICATES: mycolic acid / arabinogalactan binding proteins → Mycobacteriales",
        "NOTE: Nocardia/Rhodococcus PHROG = Mycobacteriales signal shared with Gordonia. D29-like / actinobacteriophage cluster (A-Z) genes → Mycolicibacterium. 'Cluster G' → Gordonia. Without cluster markers, both Gordonia and Mycolicibacterium are valid; use lab/soil ecology.",
        "✗ EXCLUDES: Propionibacterium/Cutibacterium PHROG → NOT Mycolicibacterium",
        "✗ EXCLUDES: Firmicutes WTA-binding → NOT Mycolicibacterium",
        "✗ EXCLUDES: 'Cluster G' Gordonia-specific genes → NOT Mycolicibacterium",
    ],
    "Gordonia_terrae": [
        "✓ INDICATES: 'Cluster G' phage genes / Gordonia-specific PHROG → Gordonia specifically within Mycobacteriales",
        "✓ INDICATES: Rhodococcus/Nocardia/Nocardiaceae PHROG in TAIL → Mycobacteriales signal (Gordonia + Mycolicibacterium + Rhodococcus all share mycolic acid cell wall; nocardia IS a valid genus indicator, not a false positive)",
        "✓ INDICATES: Corynebacterium/Nocardiaceae tail protein annotation → Mycobacteriales branch, NOT Cutibacterium (Propionibacteriales)",
        "NOTE: Nocardia/Rhodococcus PHROG alone = Mycobacteriales signal shared by Gordonia AND Mycolicibacterium. 'Cluster G' marker or environmental/industrial ecology → Gordonia. Lab-context / D29-like → Mycolicibacterium.",
        "✗ EXCLUDES: Propionibacterium/Cutibacterium PHROG → NOT Gordonia (different Actinobacteria order)",
    ],
    "Sinorhizobium_meliloti": [
        "✓ INDICATES: Alphaproteobacteria PHROG hits in phage → Sinorhizobium candidate",
        "✓ INDICATES: Giant Myoviridae (>300 ORF) with Alphaproteobacteria (n>1) taxonomic context → can infect Sinorhizobium/Rhizobiales, NOT exclusively E. coli",
        "✓ INDICATES: T4-like Myovirus + Rhizobiaceae/Rhizobiales PHROG hits → Sinorhizobium host",
        "✗ EXCLUDES: Brucellaceae PHROG hits > Rhizobiaceae hits → likely Brucella, not Sinorhizobium",
    ],
    "Cutibacterium_acnes": [
        "✓ INDICATES: Propionibacterium/Cutibacterium-specific PHROG → Cutibacterium",
        "✓ INDICATES: PAD20/PHL112M-like (skin-adapted Siphoviridae) genes → Cutibacterium",
        "✓ INDICATES: Microbacteriaceae/Actinomycetes tail protein hits CAN indicate Cutibacterium — Propionibacteriales phages cross-annotate to Microbacteriaceae in tail proteins. Ecology distinguishes: skin/sebaceous → Cutibacterium; soil/industrial → Microbacterium/Arthrobacter.",
        "✓ INDICATES: Carbohydrate-active enzyme (glucanase, polysaccharide-binding domain) in small Actinobacteria Siphoviridae + Microbacteriaceae tail cross-annotation → may indicate Cutibacterium phage targeting CWPS cell wall polysaccharide receptor.",
        "✗ EXCLUDES: Rhodococcus/Gordonia/Mycolicibacterium annotation in TAIL proteins → NOT Cutibacterium (tail proteins are host-specific; mycolic acid cell wall binding = Gordonia/Mycobacteriales)",
        "NOTE: Rhodococcus/Nocardia annotation in HOLIN or ENDOLYSIN (lysis cassette) → LOW WEIGHT cross-annotation; does NOT exclude Cutibacterium. Only TAIL protein Rhodococcus/Gordonia annotations matter.",
        "✗ EXCLUDES: mycolic acid binding proteins in tail → NOT Cutibacterium (lacks mycolic acids)",
    ],
    "Flavobacterium_psychrophilum": [
        "✓ INDICATES: Bacteroidota-adapted phage genes → Flavobacterium",
        "✗ EXCLUDES: OmpC/OmpF/LPS Enterobacteriaceae signals → NOT Flavobacterium (Bacteroidota outer membrane differs)",
    ],
    "Pseudoalteromonas_marina": [
        "✓ INDICATES: PM2-like (membrane-containing) / marine Podoviridae genes → Pseudoalteromonas",
        "✓ INDICATES: cold-active / marine-adapted phage annotations → marine Alteromonadales",
        "✗ EXCLUDES: Proteus/swarming-flagella targeting → NOT Pseudoalteromonas (marine, lacks swarming LPS)",
        "✗ EXCLUDES: clinical Enterobacteriaceae phage genes → NOT Pseudoalteromonas",
    ],
    "Vibrio_cholerae": [
        "✓ INDICATES: ICP1/ICP2/ICP3 phage genes → Vibrio cholerae O1 El Tor specific",
        "✓ INDICATES: O1 O-antigen targeting → V. cholerae",
        "✗ EXCLUDES: Agrobacterium/Rhizobium (soil symbiont) context → NOT Vibrio",
    ],
    "Pectobacterium_carotovorum": [
        "✓ INDICATES: pectate lyase / polysaccharide depolymerase TAIL SPIKE in PLANT pathogen context → Pectobacterium",
        "✗ EXCLUDES: P22-like (Salmonella-specific) genes → NOT Pectobacterium",
        "✗ EXCLUDES: Salmonella/E. coli OmpC/OmpF signals → NOT Pectobacterium",
        "✗ EXCLUDES: 'peptidoglycan catabolic process' in tail fiber → NOT Pectobacterium indicator",
        "✗ EXCLUDES: Rahnella/Enterobacter/Salmonella TAIL protein annotations → NOT Pectobacterium (mammalian Enterobacteriaceae)",
        "NOTE: Pectobacterium annotations in HEAD genes (terminase, capsid, regulators) are NOT host-specific; only TAIL spike pectate lyase → Pectobacterium",
    ],
    "Brucella_inopinata": [
        "✓ INDICATES: 'brucella' in Gram-evidence section + Brucellaceae TAIL proteins → Brucella host",
        "✓ INDICATES: Brucellaceae PHROG tail protein hits > Rhizobiaceae head protein hits → Brucella wins",
        "✓ INDICATES: GTA (gene transfer agent) phage-like genes [common in Alphaproteobacteria including Brucella] → Brucella/Rhodobacter branch",
        "✗ EXCLUDES: nitrogen-fixing symbiosis genes → NOT Brucella (Sinorhizobium-specific)",
        "NOTE: HEAD protein annotations (Rhizobiaceae) are less reliable for host prediction than TAIL protein annotations (Brucellaceae). Trust tail section over head section for host specificity.",
    ],
    "Brevibacillus_laterosporus": [
        "✓ INDICATES: Paenibacillaceae PHROG annotations → Brevibacillus (Paenibacillaceae)",
        "✓ INDICATES: Bacillus phage gene hits + spore-forming context → Brevibacillus/Bacillus, NOT Staphylococcus",
        "✗ EXCLUDES: GlcNAc-ribitol WTA-binding → NOT Brevibacillus (WTA receptor = Staphylococcus/Listeria specific)",
    ],
}

# ── Curated known phage receptors by species ───────────────────────────────────
KNOWN_RECEPTORS = {
    "Escherichia_coli":
        "LPS (O-antigen/core), OmpC, OmpF, LamB (maltoporin), BtuB, FhuA, FepA, TolC, flagella (FliC), F-pilus",
    "Salmonella_enterica":
        "LPS (O-antigen/core), OmpC, OmpF, BtuB, FhuA, flagella, Vi-antigen (typhi)",
    "Klebsiella_pneumoniae":
        "K-antigen (capsular polysaccharide — primary), LPS O-antigen, OmpK35/OmpK36, type 1 fimbriae",
    "Klebsiella_aerogenes":
        "K-antigen capsule (less diverse than K. pneumoniae), LPS O-antigen, OmpC/OmpF porins; flagella (fENko phages use flagellum as receptor)",
    "Klebsiella_oxytoca":
        "K-antigen capsule (different serotype spectrum than K. pneumoniae), LPS O-antigen, OmpC/OmpF porins",
    "Pseudomonas_aeruginosa":
        "LPS O-antigen (O-specific chain), flagella, type IV pili, outer-membrane proteins (OprD/OprF)",
    "Pseudomonas_fluorescens": "LPS O-antigen, outer membrane proteins, flagella",
    "Pseudomonas_putida":      "LPS O-antigen, outer membrane proteins",
    "Pseudomonas_syringae":    "LPS O-antigen, type IV pili, flagella",
    "Staphylococcus_aureus":
        "WTA (wall teichoic acid, GlcNAc-modified), LTA, SpA (protein A), peptidoglycan",
    "Staphylococcus_epidermidis": "WTA (ribitol-phosphate backbone), LTA, PNAG biofilm polysaccharide",
    "Staphylococcus_haemolyticus": "WTA, LTA, surface proteins",
    "Streptococcus_thermophilus": "rhamnosyl polysaccharide (cell-wall), Pip (phage infection protein)",
    "Streptococcus_pneumoniae":   "choline-containing WTA, capsular polysaccharide, PspC",
    "Streptococcus_mutans":   "cell-wall polysaccharide, surface proteins",
    "Streptococcus_pyogenes": "hyaluronic acid capsule, M-protein, cell-wall polysaccharide",
    "Lactococcus_lactis":
        "polysaccharide pellicle (pip), cell-wall polysaccharide (CWPS), membrane lipids",
    "Lactobacillus_plantarum":    "S-layer protein, cell-wall polysaccharide, teichoic acids",
    "Listeria_monocytogenes":
        "WTA (rhamnose moiety of GlcNAc-rhamnosyl polymer), N-acetylglucosamine moieties",
    "Bacillus_subtilis":
        "WTA, glucosylated poly-ribitol phosphate, flagella, peptidoglycan",
    "Bacillus_anthracis":    "WTA, poly-γ-D-glutamic acid capsule, S-layer",
    "Mycolicibacterium_smegmatis":
        "mycolic acid cell wall, arabinogalactan-peptidoglycan, glycopeptidolipids (GPL)",
    "Mycobacterium_avium":
        "glycopeptidolipids (GPL), mycolic acids, arabinogalactan-peptidoglycan",
    "Mycobacteroides_abscessus":
        "glycopeptidolipids (GPL — rough/smooth morphotype), mycolic acids",
    "Gordonia_terrae":       "mycolic acid cell wall, arabinogalactan, surface glycolipids",
    "Gordonia_rubripertincta": "mycolic acid cell wall, surface glycolipids",
    "Gordonia_sputi":        "mycolic acid cell wall, surface glycolipids",
    "Gordonia_alkanivorans": "mycolic acid cell wall, surface glycolipids",
    "Gordonia_malaquae":     "mycolic acid cell wall, surface glycolipids",
    "Gordonia_neofelifaecis": "mycolic acid cell wall, surface glycolipids",
    "Vibrio_cholerae":
        "O-antigen (O1/O139 predominant), mannose-sensitive hemagglutinin pili (MSHA), toxin-coregulated pili (TCP)",
    "Vibrio_parahaemolyticus": "O-antigen, flagella, type IV pili, outer membrane proteins",
    "Vibrio_harveyi":          "O-antigen, outer membrane proteins, pili",
    "Vibrio_natriegens":       "O-antigen, outer membrane proteins",
    "Vibrio_alginolyticus":    "O-antigen, flagella, outer membrane proteins",
    "Vibrio_splendidus":       "O-antigen, outer membrane proteins",
    "Vibrio_vulnificus":       "O-antigen, capsular polysaccharide, outer membrane proteins",
    "Helicobacter_pylori":
        "LPS (Lewis antigen-decorated), outer membrane proteins (BabA, HopQ, OipA), flagella",
    "Campylobacter_jejuni":
        "flagella (FlaA/FlaB), capsular polysaccharide (CPS), LPS core",
    "Acinetobacter_baumannii":
        "OmpA (major outer membrane protein), capsular polysaccharide (K-locus), LPS (lipooligosaccharide)",
    "Burkholderia_cenocepacia": "LPS O-antigen, type IV pili, cable pili, outer membrane proteins",
    "Burkholderia_pseudomallei": "LPS O-antigen, type IV pili, capsular polysaccharide",
    "Streptomyces_coelicolor":  "spore surface glycan (SapB), rodlet layer, BldN",
    "Streptomyces_griseus":     "spore surface glycan, spore coat proteins",
    "Mycoplasma_pulmonis":      "spike proteins (P97-like adhesins), membrane proteins",
    "Caulobacter_vibrioides":   "holdfast polysaccharide, flagellum, pili (CbpA)",
    "Synechococcus_sp._WH_8102": "outer membrane proteins (OmpA-like), carbohydrate-binding",
    "Rhizobium_leguminosarum":  "LPS O-antigen, K-antigen EPS, flagella, type IV pili",
    "Xanthomonas_campestris":   "LPS O-antigen, xanthan EPS, type IV pili",
    "Erwinia_amylovora":        "LPS O-antigen, amylovoran EPS, type IV pili",
    "Pectobacterium_carotovorum": "LPS O-antigen, outer membrane proteins",
    "Cutibacterium_acnes":
        "cell wall polysaccharide (CWPS), surface proteins (CAMP factor, lipase); Propionibacteriales-specific phage receptors; skin-adapted",
    "Lactobacillus_fermentum":
        "cell wall polysaccharide, WTA (wall teichoic acids), S-layer proteins (when present)",
    "Pseudoalteromonas_marina":
        "LPS O-antigen, outer membrane proteins (OMP); PM2-like phages use membrane lipid bilayer; Podoviridae and Siphoviridae in marine environments",
    "Citrobacter_freundii":     "LPS O-antigen, OmpC/OmpF porins, flagella",
    "Ralstonia_solanacearum":
        "LPS O-antigen, type IV pili (pilA), outer membrane proteins; plant pathogen — phage receptors distinct from mammalian pathogens",
    "Stenotrophomonas_maltophilia":
        "LPS O-antigen, OMPs (StmOmpA), flagella; Xanthomonadaceae (distinct from Pseudomonas)",
    "Sinorhizobium_meliloti":  "LPS O-antigen, K-antigen EPS, flagella, pili; nitrogen-fixing symbiont",
    "Morganella_morganii":     "LPS O-antigen, outer membrane proteins, flagella",
    "Proteus_mirabilis":       "LPS O-antigen, flagella (FliC), outer membrane proteins",
    "Citrobacter_rodentium":    "LPS O-antigen, intimin (outer membrane adhesin), flagella",
    "Shigella_flexneri":        "LPS O-antigen (type-specific), IcsA surface protein",
    "Shigella_boydii":          "LPS O-antigen (type-specific)",
    "Shigella_sonnei":          "LPS O-antigen (form I/II), outer membrane proteins",
    "Yersinia_enterocolitica":  "LPS O-antigen, Ail outer membrane protein, flagella",
    "Yersinia_pestis":          "LPS (defective O-chain), Psa fimbriae, Caf1 capsule",
    "Bacteroides_fragilis":     "surface polysaccharides (PSA–PSH), outer membrane proteins",
    "Clostridioides_difficile": "S-layer (SlpA), cell-wall polysaccharide (PS-II), flagella",
    "Mycoplasma_arthritidis":   "membrane lipoproteins, surface-exposed proteins (MAA1/2)",
    "Haloarcula_hispanica":  "S-layer glycoprotein, halocell envelope, archaeal lipids",
    "Haloarcula_californiae": "S-layer glycoprotein, halocell envelope",
    "Haloarcula_vallismortis": "S-layer glycoprotein, halocell envelope",
    "Haloarcula_sinaiiensis": "S-layer glycoprotein, halocell envelope",
    "Halorubrum_coriense":   "S-layer glycoprotein, ether-linked archaeal lipids",
    "Sulfolobus_islandicus": "S-layer (SlaA/SlaB glycoprotein), archaeal glycolipids",
    "Saccharolobus_solfataricus": "S-layer glycoprotein, archaeal tetraether lipids",
    "Acidianus_hospitalis":  "S-layer, archaeal glycolipids, extreme-acidophile cell envelope",
    "Aeropyrum_pernix":      "S-layer (hyperthermophile-adapted), archaeal ether lipids",
}

# ── Gram type by genus ─────────────────────────────────────────────────────────
GRAM_MAP = {
    # Gram-negative
    "Achromobacter":"-","Acinetobacter":"-","Aeromonas":"-","Agrobacterium":"-",
    "Aliivibrio":"-","Alteromonas":"-","Azospirillum":"-","Bacteroides":"-",
    "Bdellovibrio":"-","Brucella":"-","Burkholderia":"-","Campylobacter":"-",
    "Candidatus_Hamiltonella":"-","Candidatus_Liberibacter":"-",
    "Candidatus_Pelagibacter":"-","Candidatus_Puniceispirillum":"-",
    "Caulobacter":"-","Cellulophaga":"-","Chlamydia":"-","Citrobacter":"-",
    "Clavibacter":"+",
    "Colwellia":"-","Croceibacter":"-","Cronobacter":"-","Delftia":"-",
    "Dickeya":"-","Dinoroseobacter":"-","Edwardsiella":"-","Enterobacter":"-",
    "Erwinia":"-","Escherichia":"-","Flavobacterium":"-","Glaesserella":"-",
    "Helicobacter":"-","Klebsiella":"-","Mannheimia":"-","Mesorhizobium":"-",
    "Microcystis":"-","Morganella":"-","Myxococcus":"-","Pantoea":"-",
    "Parabacteroides":"-","Pasteurella":"-","Pectobacterium":"-",
    "Prochlorococcus":"-","Proteus":"-","Providencia":"-",
    "Pseudoalteromonas":"-","Pseudomonas":"-","Ralstonia":"-","Rhizobium":"-",
    "Rhodobacter":"-","Rhodovulum":"-","Roseobacter":"-","Ruegeria":"-",
    "Salinivibrio":"-","Salmonella":"-","Serratia":"-","Shewanella":"-",
    "Shigella":"-","Sinorhizobium":"-","Sodalis":"-","Stenotrophomonas":"-",
    "Sulfitobacter":"-","Synechococcus":"-","Vibrio":"-","Xanthomonas":"-",
    "Xylella":"-","Yersinia":"-","Planktothrix":"-",
    # Gram-positive
    "Actinomyces":"+","Arthrobacter":"+","Bacillus":"+","Brevibacillus":"+",
    "Brochothrix":"+","Clostridium":"+","Clostridioides":"+",
    "Cutibacterium":"+","Enterococcus":"+","Erysipelothrix":"+",
    "Geobacillus":"+","Gordonia":"+","Kitasatospora":"+","Lactobacillus":"+",
    "Lactococcus":"+","Leuconostoc":"+","Listeria":"+","Microbacterium":"+",
    "Mycobacterium":"+","Mycobacteroides":"+","Mycolicibacterium":"+",
    "Mycoplasma":"+","Nocardia":"+","Oenococcus":"+","Paenibacillus":"+",
    "Propionibacterium":"+","Rhodococcus":"+","Staphylococcus":"+",
    "Streptococcus":"+","Streptomyces":"+","Thermoanaerobacterium":"+",
    "Thermus":"+","Trichormus":"+","Tsukamurella":"+","Weissella":"+",
    # Archaea
    "Acidianus":"archaea","Aeropyrum":"archaea","Haloarcula":"archaea",
    "Halorubrum":"archaea","Methanothermobacter":"archaea",
    "Pyrobaculum":"archaea","Pyrococcus":"archaea","Saccharolobus":"archaea",
    "Sulfolobus":"archaea",
}

# ── Habitat / ecological context ──────────────────────────────────────────────
HABITAT_MAP = {
    "Escherichia_coli": "mammalian gut (commensal/pathogen); soil; water; food",
    "Salmonella_enterica": "mammalian gut (enteric pathogen); food; environment",
    "Klebsiella_pneumoniae": "mammalian gut; respiratory tract; hospital environment; soil",
    "Klebsiella_aerogenes": "mammalian gut; soil; water; nosocomial",
    "Klebsiella_oxytoca": "mammalian gut; soil; water; nosocomial",
    "Pseudomonas_aeruginosa": "soil; water; plant rhizosphere; hospital (opportunistic pathogen)",
    "Pseudomonas_fluorescens": "soil; plant rhizosphere; freshwater; food spoilage",
    "Pseudomonas_putida": "soil; plant rhizosphere; freshwater",
    "Pseudomonas_syringae": "plant pathogen; phyllosphere; soil",
    "Pseudomonas_tolaasii": "soil; mushroom pathogen",
    "Staphylococcus_aureus": "human skin/nares (commensal/pathogen); hospital; food",
    "Staphylococcus_epidermidis": "human skin (commensal); hospital; biofilm-former",
    "Staphylococcus_haemolyticus": "human skin; opportunistic nosocomial pathogen",
    "Staphylococcus_saprophyticus": "human urinary tract; skin; food",
    "Staphylococcus_capitis": "human skin/scalp; neonatal nosocomial pathogen",
    "Staphylococcus_hominis": "human skin; nosocomial",
    "Staphylococcus_pasteuri": "soil; water; clinical isolates",
    "Staphylococcus_xylosus": "food fermentation; animal skin",
    "Streptococcus_pyogenes": "human throat/skin (Group A Strep pathogen)",
    "Streptococcus_pneumoniae": "human nasopharynx (pathogen); community/hospital",
    "Streptococcus_mutans": "human oral cavity; dental caries",
    "Streptococcus_thermophilus": "dairy fermentation (yogurt/cheese starter)",
    "Streptococcus_suis": "pig respiratory tract/tonsils; zoonotic pathogen",
    "Streptococcus_dysgalactiae": "animal/human (Group C/G Strep); bovine mastitis",
    "Streptococcus_gordonii": "human oral cavity; dental plaque biofilm",
    "Streptococcus_mitis": "human oral cavity; opportunistic",
    "Streptococcus_oralis": "human oral cavity; dental plaque",
    "Streptococcus_salivarius": "human oral cavity/gut; probiotic",
    "Streptococcus_parauberis": "fish/marine aquaculture pathogen",
    "Lactococcus_lactis": "dairy fermentation (cheese/butter starter); plant surfaces",
    "Lactococcus_garvieae": "fish/bovine mastitis pathogen; food",
    "Lactobacillus_plantarum": "fermented foods; plant material; gut",
    "Lactobacillus_casei": "dairy; gut probiotic; fermented foods",
    "Lactobacillus_rhamnosus": "gut (probiotic); dairy; vaginal",
    "Lactobacillus_delbrueckii": "dairy fermentation (yogurt)",
    "Lactobacillus_fermentum": "fermented foods; gut",
    "Lactobacillus_gasseri": "human gut; vaginal",
    "Lactobacillus_jensenii": "vaginal; gut",
    "Lactobacillus_johnsonii": "gut; dairy",
    "Lactobacillus_paracasei": "dairy; gut; plant material",
    "Bacillus_subtilis": "soil; plant rhizosphere; food spoilage; model organism",
    "Bacillus_anthracis": "soil (spores); anthrax pathogen",
    "Bacillus_cereus": "soil; food (food poisoning pathogen)",
    "Bacillus_thuringiensis": "soil; insect pathogen (Bt toxin); biocontrol",
    "Bacillus_megaterium": "soil; plant rhizosphere; industrial applications",
    "Bacillus_pumilus": "soil; marine; plant surfaces",
    "Bacillus_alcalophilus": "alkaline soil; soda lakes",
    "Bacillus_halmapalus": "soil; alkaliphile",
    "Bacillus_mycoides": "soil; plant rhizosphere",
    "Listeria_monocytogenes": "soil; food (listeriosis pathogen); hospital",
    "Acinetobacter_baumannii": "hospital environment (nosocomial pathogen); soil; water",
    "Acinetobacter_johnsonii": "soil; water; human skin",
    "Acinetobacter_soli": "soil",
    "Vibrio_cholerae": "aquatic (brackish/estuarine); cholera pathogen",
    "Vibrio_parahaemolyticus": "marine/estuarine; seafood (gastroenteritis)",
    "Vibrio_alginolyticus": "marine; fish/shellfish pathogen",
    "Vibrio_harveyi": "marine; shrimp/fish aquaculture pathogen; bioluminescent",
    "Vibrio_natriegens": "salt marsh; estuarine; fastest-growing bacterium",
    "Vibrio_splendidus": "marine; bivalve pathogen",
    "Vibrio_vulnificus": "marine/estuarine; shellfish; wound/septicemia pathogen",
    "Helicobacter_pylori": "human gastric mucosa (peptic ulcer/gastric cancer)",
    "Campylobacter_jejuni": "poultry/cattle gut (commensal); foodborne pathogen",
    "Aeromonas_hydrophila": "freshwater; fish pathogen; opportunistic human pathogen",
    "Aeromonas_salmonicida": "freshwater fish (furunculosis pathogen)",
    "Aeromonas_media": "freshwater; environmental",
    "Burkholderia_cenocepacia": "soil; lung (CF patient opportunistic pathogen)",
    "Burkholderia_cepacia": "soil; plant rhizosphere; CF lung pathogen",
    "Burkholderia_pseudomallei": "soil/water SE Asia (melioidosis pathogen)",
    "Burkholderia_ambifaria": "soil; plant rhizosphere; CF lung",
    "Burkholderia_pyrrocinia": "soil; plant rhizosphere",
    "Burkholderia_thailandensis": "soil SE Asia; low virulence BSL-2 surrogate",
    "Mycolicibacterium_smegmatis": "soil; water; model mycobacterium (non-pathogenic)",
    "Mycolicibacterium_phlei": "soil; plant material; non-pathogenic",
    "Mycobacterium_avium": "soil; water; NTM pathogen (immunocompromised/MAC)",
    "Mycobacteroides_abscessus": "soil; water; NTM lung pathogen (CF patients)",
    "Gordonia_terrae": "soil; hydrocarbon degradation",
    "Gordonia_rubripertincta": "soil; hydrocarbon degradation",
    "Gordonia_sputi": "soil; rare human respiratory infections",
    "Gordonia_alkanivorans": "soil; oil bioremediation",
    "Gordonia_malaquae": "water treatment; nosocomial rare",
    "Gordonia_neofelifaecis": "animal gut; feces",
    "Streptomyces_coelicolor": "soil; antibiotic production model organism",
    "Streptomyces_griseus": "soil; streptomycin producer",
    "Streptomyces_avermitilis": "soil; avermectin producer",
    "Streptomyces_lividans": "soil; industrial expression host",
    "Streptomyces_venezuelae": "soil; chloramphenicol producer",
    "Streptomyces_flavovirens": "soil",
    "Mycoplasma_pulmonis": "rat/mouse respiratory tract (pathogen); cell-wall-less",
    "Mycoplasma_arthritidis": "rat/mouse joints/respiratory (pathogen); cell-wall-less",
    "Xanthomonas_campestris": "plant pathogen (brassicas); soil",
    "Xanthomonas_citri": "citrus canker pathogen; phyllosphere",
    "Xanthomonas_oryzae": "rice pathogen (bacterial blight); phyllosphere",
    "Xanthomonas_vesicatoria": "solanaceous plant pathogen",
    "Xylella_fastidiosa": "xylem-limited plant pathogen; insect-transmitted",
    "Erwinia_amylovora": "plant pathogen (fire blight of rosaceous plants)",
    "Erwinia_pyrifoliae": "Asian pear pathogen",
    "Pectobacterium_carotovorum": "plant pathogen (soft rot); soil",
    "Pectobacterium_atrosepticum": "potato pathogen (blackleg); soil",
    "Dickeya_solani": "potato/ornamental plant pathogen; soil",
    "Ralstonia_solanacearum": "soil; vascular plant pathogen (bacterial wilt)",
    "Ralstonia_pickettii": "soil; water; nosocomial opportunistic",
    "Agrobacterium_tumefaciens": "soil; plant pathogen (crown gall); rhizosphere",
    "Rhizobium_leguminosarum": "soil; plant root nodules (N-fixation with legumes)",
    "Rhizobium_etli": "soil; bean root nodules (N-fixation)",
    "Rhizobium_gallicum": "soil; legume root nodules",
    "Sinorhizobium_meliloti": "soil; alfalfa root nodules (N-fixation)",
    "Sinorhizobium_sp._LM21": "soil; legume rhizosphere",
    "Mesorhizobium_loti": "soil; lotus root nodules (N-fixation)",
    "Rhodobacter_capsulatus": "freshwater sediment; anoxygenic phototroph; model organism",
    "Dinoroseobacter_shibae": "marine; dinoflagellate symbiont",
    "Roseobacter_denitrificans": "marine; anoxygenic phototroph",
    "Ruegeria_pomeroyi": "marine; DMSP degrader; Roseobacter clade",
    "Sulfitobacter_sp._CB2047": "marine; Roseobacter clade",
    "Sulfitobacter_sp._EE-36": "marine; Roseobacter clade",
    "Rhodovulum_sp._P5": "marine/saline; phototrophic",
    "Alteromonas_macleodii": "marine; heterotroph; surface ocean",
    "Pseudoalteromonas_atlantica": "marine; biofilm-forming; EPS producer",
    "Pseudoalteromonas_marina": "marine; heterotroph",
    "Colwellia_psychrerythraea": "polar marine/sea-ice; psychrophile",
    "Shewanella_putrefaciens": "freshwater/marine; iron-reducing; psychrotolerant",
    "Cellulophaga_baltica": "marine; cellulose/algae degrader; Baltic Sea",
    "Flavobacterium_columnare": "freshwater fish pathogen (columnaris disease)",
    "Flavobacterium_psychrophilum": "freshwater fish pathogen (cold-water disease); psychrophile",
    "Synechococcus_sp._CB0101": "marine; cyanobacterium; photosynthetic",
    "Synechococcus_sp._WH_7803": "marine; open ocean cyanobacterium",
    "Synechococcus_sp._WH_7805": "marine; open ocean cyanobacterium",
    "Synechococcus_sp._WH_8102": "marine; open ocean cyanobacterium; model picocyanobacterium",
    "Synechococcus_sp._WH_8109": "marine; open ocean cyanobacterium",
    "Prochlorococcus_marinus": "marine; most abundant photosynthetic organism on Earth",
    "Microcystis_aeruginosa": "freshwater; bloom-forming cyanobacterium; microcystin producer",
    "Planktothrix_agardhii": "freshwater; bloom-forming cyanobacterium",
    "Trichormus_variabilis": "freshwater/soil; N-fixing cyanobacterium; akinete-forming",
    "Haloarcula_californiae": "saline lake (California); extreme halophile archaea",
    "Haloarcula_hispanica": "solar saltern (Spain); extreme halophile archaea",
    "Haloarcula_sinaiiensis": "hypersaline lake (Sinai); extreme halophile archaea",
    "Haloarcula_vallismortis": "Death Valley salt flat; extreme halophile archaea",
    "Halorubrum_coriense": "solar saltern (Australia); extreme halophile archaea",
    "Acidianus_hospitalis": "acidic hot spring (sulfuric); hyperthermoacidophile archaea",
    "Aeropyrum_pernix": "marine hydrothermal vent; hyperthermophile archaea; aerobic",
    "Sulfolobus_islandicus": "acidic hot spring; thermoacidophile archaea; CRISPR model",
    "Saccharolobus_solfataricus": "acidic hot spring; thermoacidophile archaea",
    "Pyrobaculum_arsenaticum": "anaerobic hot spring; hyperthermophile archaea; arsenate reducer",
    "Pyrococcus_abyssi": "deep-sea hydrothermal vent; hyperthermophile archaea; anaerobic",
    "Methanothermobacter_marburgensis": "anaerobic sludge; thermophilic methanogen archaea",
    "Thermus_thermophilus": "hot spring; thermophile bacterium; biotechnology source (Taq pol family)",
    "Geobacillus_kaustophilus": "compost/hot spring; thermophile; spore-forming",
    "Thermoanaerobacterium_saccharolyticum": "hot spring; thermophilic anaerobe; cellulose/xylan fermentation",
    "Enterococcus_faecalis": "mammalian gut (commensal/pathogen); nosocomial; dairy",
    "Enterococcus_faecium": "mammalian gut (commensal/pathogen); nosocomial; VRE",
    "Clostridioides_difficile": "mammalian gut (CDI pathogen); hospital; soil spores",
    "Clostridium_perfringens": "soil; gut; food poisoning/gas gangrene pathogen",
    "Clostridium_botulinum": "soil; food (botulism pathogen); spore-forming",
    "Clostridium_sporogenes": "soil; gut; non-toxigenic C. botulinum relative",
    "Clostridium_tetani": "soil; wound (tetanus pathogen)",
    "Bacteroides_fragilis": "mammalian gut (commensal dominant); opportunistic pathogen",
    "Parabacteroides_distasonis": "mammalian gut (commensal); gut microbiome",
    "Parabacteroides_merdae": "mammalian gut (commensal); gut microbiome",
    "Myxococcus_xanthus": "soil; predatory; fruiting body-forming model organism",
    "Bdellovibrio_bacteriovorus": "soil; freshwater; predatory bacterium (preys on Gram-)",
    "Brucella_abortus": "cattle/livestock (brucellosis zoonosis); intracellular pathogen",
    "Brucella_melitensis": "goat/sheep (brucellosis); most virulent Brucella",
    "Brucella_canis": "dog (brucellosis); zoonotic",
    "Brucella_suis": "pig (brucellosis); zoonotic",
    "Brucella_inopinata": "rare human isolate; atypical Brucella",
    "Pasteurella_multocida": "animal respiratory tract; zoonotic bite pathogen",
    "Mannheimia_haemolytica": "bovine respiratory tract (BRD pathogen); ruminant",
    "Glaesserella_parasuis": "pig respiratory tract (Glässer's disease); farm pathogen",
    "Aggregatibacter_actinomycetemcomitans": "human oral cavity; periodontal pathogen",
    "Actinomyces_naeslundii": "human oral cavity; dental plaque; root caries",
    "Cutibacterium_acnes": "human skin (sebaceous follicles); acne; low-oxygen",
    "Propionibacterium_freudenreichii": "dairy (Swiss cheese eye-forming); gut",
    "Nocardia_brasiliensis": "soil; cutaneous nocardiosis pathogen",
    "Kitasatospora_aureofaciens": "soil; antibiotic-producing actinomycete",
    "Tsukamurella_paurometabola": "soil; water; rare human opportunistic pathogen",
    "Rhodococcus_erythropolis": "soil; hydrocarbon degradation; cold-tolerant",
    "Rhodococcus_hoagii": "soil; animal infections (rare)",
    "Rhodococcus_rhodochrous": "soil; hydrocarbon degradation; industrial biocatalysis",
    "Microbacterium_oxydans": "soil; plant material; food; rare clinical",
    "Arthrobacter_sp._ATCC_21022": "soil; cold-tolerant; hydrocarbon degradation",
    "Azospirillum_brasilense": "soil; plant rhizosphere; N-fixing; plant growth-promoting",
    "Croceibacter_atlanticus": "marine; Atlantic Ocean; Flavobacteriia",
    "Candidatus_Hamiltonella_defensa": "insect gut (aphid endosymbiont); obligate intracellular",
    "Candidatus_Liberibacter_asiaticus": "citrus phloem (HLB citrus greening pathogen); insect-vectored",
    "Candidatus_Pelagibacter_ubique": "marine open ocean; most abundant heterotroph on Earth; SAR11",
    "Candidatus_Puniceispirillum_marinum": "marine; SAR116 clade; phototrophic",
    "Delftia_acidovorans": "soil; freshwater; gold-precipitating; rare clinical",
    "Delftia_sp._670": "environmental; soil/water",
    "Delftia_tsuruhatensis": "soil; wastewater; rare clinical",
    "Leuconostoc_mesenteroides": "plant surfaces; fermented vegetables/dairy",
    "Leuconostoc_pseudomesenteroides": "fermented foods; dairy",
    "Oenococcus_oeni": "wine (malolactic fermentation); acid-tolerant",
    "Weissella_cibaria": "fermented food; plant material; gut",
    "Brochothrix_thermosphacta": "meat spoilage; food refrigeration environments",
    "Erysipelothrix_rhusiopathiae": "animal tissues (swine erysipelas); zoonotic",
    "Paenibacillus_larvae": "honeybee larvae (American foulbrood pathogen)",
    "Brevibacillus_laterosporus": "soil; insect pathogen; spore-forming",
    "Stenotrophomonas_maltophilia": "soil; water; plant rhizosphere; nosocomial MDR pathogen",
    "Caulobacter_vibrioides": "oligotrophic freshwater; dimorphic (stalked/swarmer); model organism",
    "Morganella_morganii": "mammalian gut; opportunistic nosocomial pathogen",
    "Proteus_mirabilis": "mammalian gut; soil; UTI pathogen; swarming",
    "Providencia_stuartii": "mammalian gut; nosocomial UTI pathogen",
    "Salinivibrio_costicola": "saline/hypersaline environments; moderate halophile",
    "Sodalis_glossinidius": "tsetse fly gut (secondary endosymbiont)",
    "Chlamydia_abortus": "ruminant placenta (ovine enzootic abortion); obligate intracellular",
    "Chlamydia_pecorum": "ruminant gut/respiratory; obligate intracellular",
    "Chlamydia_pneumoniae": "human respiratory tract (pneumonia); obligate intracellular",
}

# ── Known phage families / groups by genus ─────────────────────────────────────
PHAGE_FAMILIES = {
    "Escherichia": "T4-like (Myoviridae), T7/T3-like (Autographiviridae), λ-like (Siphoviridae), P1-like, P2-like, Mu-like, Ff-like (M13/fd ssDNA), Qβ/MS2 (ssRNA), ΦX174 (ssDNA)",
    "Shigella": "λ-like (Siphoviridae), T4-like (Myoviridae), Sf6/P22-like (Podoviridae); Shigella-specific O-antigen LPS receptors; related to E. coli but serotype-specific",
    "Salmonella": "P22-like (Podoviridae), SP6-like, ε15-like, Myoviridae (Felix-O1), Siphoviridae (HK97-like), P2-like Myoviridae",
    "Klebsiella": "KP series (Myoviridae/Siphoviridae/Podoviridae), Drulisvirus, Sugarlandvirus; capsule (K-antigen) is primary determinant",
    "Pseudomonas_aeruginosa": "ΦKZ/LMA2 giant Myoviridae (>150 ORF, DEFINITIVE P. aeruginosa signal), LUZ19/LKD16 (Autographiviridae), PaP1 (Myoviridae), Siphoviridae (LUZ24-like); receptors: LPS O-antigen (O5 most common), type IV pili, LPS core",
    "Pseudomonas_fluorescens": "ΦS1 (Myoviridae), PfCl01/LUZ24-like (Podoviridae); soil/food spoilage; LPS O-antigen (distinct from P. aeruginosa)",
    "Pseudomonas_syringae": "ΦPsyM2 (Myoviridae), plant-pathogen-adapted phages; LPS O-antigen",
    "Pseudomonas_putida": "diverse Siphoviridae, Podoviridae; soil saprophyte; LPS O-antigen",
    "Pseudomonas": "Pseudomonas phages (Myoviridae/Siphoviridae/Autographiviridae); LPS O-antigen and pili receptors",
    "Staphylococcus": "Kayvirus/K-like (Myoviridae), Twort-like, phiMR11-like (Siphoviridae); WTA (wall teichoic acid) is primary receptor",
    "Streptococcus": "Skunavirus (P335-like), Siphoviridae (λSa series), C1 phage; polysaccharide capsule and cell-wall carbohydrates as receptors",
    "Lactococcus": "936/P335/c2 groups (Siphoviridae); phage-host coevolution classic model; receptor: polysaccharide pellicle (pip)",
    "Lactobacillus": "A2-like, ΦadE, Myoviridae and Siphoviridae groups; polysaccharide and WTA receptors",
    "Bacillus": "SPO1/SPβ (Myoviridae), PBS series, PMBT28, φ29 (Podoviridae, transducing); WTA, polysaccharide, flagella receptors",
    "Listeria": "A511 (Myoviridae), P100-like, PSA (Siphoviridae), φA118; WTA (rhamnose moiety) is primary receptor",
    "Mycolicibacterium": "Cluster A (D29-like) to Cluster Z; remarkable diversity; ~7,000+ isolated phages; cell-wall (mycolic acid) binding",
    "Mycobacterium": "D29-like, Bxz1, Trixie; mycobacterial cell envelope; cluster classification system (phagesDB)",
    "Mycobacteroides": "similar to Mycolicibacterium phage clusters; rough vs smooth morphotype susceptibility",
    "Gordonia": "Cluster G Gordonia phages; Siphoviridae; wide diversity; mycolic-acid cell wall",
    "Vibrio": "ICP1/ICP2/ICP3 (O1 El Tor-specific), K139, VP882-like; O-antigen (O1/O139) and pili receptors",
    "Campylobacter": "NCTC 12673, CP220, CPS5; flagella and capsular polysaccharide receptors",
    "Acinetobacter": "AB1 (Myoviridae), Fri1, SH-Ab15519; OmpA and capsule receptors",
    "Burkholderia": "ΦE125 (Myoviridae), BcepMu, KL3; LPS and pili receptors",
    "Aeromonas": "Aeromonas phages (Myoviridae/Siphoviridae); O-antigen and pili",
    "Enterococcus": "IMEEF1 (Myoviridae), phiEF24C, EFRM31; polysaccharide and WTA",
    "Clostridium": "ΦCD27, phi3626 (Siphoviridae); cell-wall polysaccharide",
    "Clostridioides": "CDKM9, phi-CD119, phiC2; S-layer and cell-wall polysaccharide",
    "Helicobacter": "KHP30, HP1, ΦHP33; Helicobacter-specific lipopolysaccharide",
    "Streptomyces": "ΦC31 (integrating), SV1 (Myoviridae), RP2/RP3; spore coat receptors",
    "Synechococcus": "Syn5, S-PM2 (Myoviridae), S-RIM (Siphoviridae); OMP receptors; cyanophage",
    "Prochlorococcus": "P-SSM2/P-SSP7; cyanophages; OMP receptors",
    "Mycoplasma": "P1 phages; spike proteins; limited due to parasitic lifestyle",
    "Haloarcula": "SH1 (icosahedral), ΦH (tailless); archaeal lipid membrane",
    "Halorubrum": "His1 (lemon-shaped), HF1/HF2; archaeal; lipid and s-layer receptors",
    "Sulfolobus": "SIRV1/2 (rod-shaped), STIV, SSV1 (spindle); crenarchaeal S-layer",
    "Acidianus": "ATV (spindle/extracellular tail), SIRV-like; S-layer; extreme conditions",
    "Xanthomonas": "ΦXv, ΦXca, XF phages; LPS and pili receptors",
    "Erwinia": "Era103 (Myoviridae), φEa1 (Podoviridae); O-antigen receptors",
    "Pectobacterium": "ΦTE, PP90 (Myoviridae); O-antigen; similar to Erwinia phages",
    "Rhizobium": "ΦM12, RL38 (Myoviridae/Siphoviridae); LPS O-antigen; flagella",
    "Agrobacterium": "PW8 (Myoviridae); LPS receptors",
    "Caulobacter": "ΦCbK (giant Myoviridae), Cr30; holdfast and pili receptors",
    "Rhodobacter": "ΦRS1, RcapMu; LPS and pili receptors",
    "Flavobacterium": "FCL-2 (Myoviridae), FpsP1 (Podoviridae); cell-wall glycan",
    "Bacteroides": "ΦB40-8, ΦNT1 (Myoviridae); O-antigen and surface polysaccharide",
    "Clavibacter": "CMP1, CM1/55A (Myoviridae/Siphoviridae); cell-wall polysaccharide",
    "Myxococcus": "Mx8 (integrating Siphoviridae), Mx4; LPS O-antigen",
    "Brucella": "BTP1 (Podoviridae), Tb (BTP-like), Wb; LPS; intracellular pathogen reduces phage access",
    "Pasteurella": "ΦPM2, Φ14 (Myoviridae/Siphoviridae); LPS O-antigen",
    "Citrobacter": "Cf1/CfP1 (Siphoviridae/Podoviridae), T4-like giant Myoviridae; O-antigen LPS receptors; Enterobacteriaceae (mammalian gut)",
    "Ralstonia": "ΦRSS1/ΦRSLseries (Podoviridae/Myoviridae); LPS O-antigen; plant pathogen",
    "Cutibacterium": "PAD20/PHL112M-like (Siphoviridae), pac1/phi11-like; skin-adapted Actinobacteria phages; Propionibacteriales host",
    "Stenotrophomonas": "S1 (Siphoviridae), IME-SM1, Smp14 (Myoviridae); LPS O-antigen; Xanthomonadaceae",
    "Sinorhizobium": "ΦM9 (Siphoviridae), ΦN3; LPS O-antigen; nitrogen-fixing symbiont",
    "Morganella": "ΦMP-1 (Podoviridae); O-antigen LPS; Enterobacteriaceae",
    "Proteus": "ΦPV22 (Myoviridae), ΦKO2; LPS O-antigen; flagella receptor",
    "Achromobacter": "ΦAX (Myoviridae); diverse Betaproteobacteria phages; LPS receptors",
    "Pantoea": "LiMac (Siphoviridae); O-antigen LPS; Enterobacteriaceae",
    "Paenibacillus": "ΦPaX (Myoviridae/Siphoviridae); WTA and polysaccharide receptors; Gram-positive",
    "Geobacillus": "ΦGP1 (Siphoviridae); thermophile (55-70°C optimal); WTA receptor",
    "Brevibacillus": "ΦNY36 (Myoviridae); WTA and polysaccharide; insecticidal strains common",
    "Rhodococcus": "diverse Actinobacteria phages (Siphoviridae); mycolic acid cell wall; soil degrader",
    "Microcystis": "cyanophages (Ma-LMM01); cyanobacterial outer membrane; NO cross-infection with Enterobacteriaceae",
    "Pseudoalteromonas": "PM2 (membrane-containing Corticoviridae), H105/1 (Myoviridae), diverse Siphoviridae/Podoviridae; marine-adapted Alteromonadales",
    "Shewanella": "Sfn1 (Podoviridae); outer membrane receptors; iron-reducing; marine/aquatic",
}


# ── Parse base profile ─────────────────────────────────────────────────────────
def parse_base_profile(profile_path: Path) -> dict:
    """Extract defense systems and HIGH-confidence surface receptor genes.

    Parses a base profile produced by build_host_profiles_base.py.

    Returns:
        dict with keys:
          'defense'    : {system_type: count}
          'receptors'  : list of "gene: description" strings (top 12 [HIGH] genes)
          'total_orfs' : int
    """
    if not profile_path.exists():
        return {"defense": {}, "receptors": [], "total_orfs": 0}

    text = profile_path.read_text()

    m = re.search(r"Total ORFs: (\d+)", text)
    total_orfs = int(m.group(1)) if m else 0

    defense: dict[str, int] = {}
    in_defense = False
    for line in text.splitlines():
        if "## Defense Systems" in line:
            in_defense = True
            continue
        if in_defense and line.startswith("## "):
            break
        if in_defense and line.startswith("- ") and "(" in line:
            m2 = re.match(r"- (\S+) \((\d+) system", line)
            if m2:
                defense[m2.group(1)] = int(m2.group(2))

    receptors: list[str] = []
    in_surface = False
    for line in text.splitlines():
        if "## Surface Receptors" in line:
            in_surface = True
            continue
        if in_surface and line.startswith("## "):
            break
        if in_surface and "[HIGH]" in line:
            m3 = re.match(r"- \[HIGH\] (\w+) \(([^)]{10,80})", line)
            if m3:
                name, desc = m3.group(1), m3.group(2)
                receptors.append(f"{name}: {desc[:60].rstrip()}")
        if len(receptors) >= 12:
            break

    return {"defense": defense, "total_orfs": total_orfs, "receptors": receptors}


# ── NCBI fetchers ──────────────────────────────────────────────────────────────
def _ncbi_get(url: str, api_key: str | None, retries: int = 3) -> bytes | None:
    if api_key:
        url += f"&api_key={api_key}"
    delay = RATE_DELAY_KEY if api_key else RATE_DELAY
    for attempt in range(retries):
        try:
            time.sleep(delay)
            with urllib.request.urlopen(url, timeout=20) as r:
                return r.read()
        except Exception as e:
            if attempt == retries - 1:
                print(f"  NCBI fetch failed: {url[:80]}... → {e}")
            time.sleep(1.0 * (attempt + 1))
    return None


def fetch_taxonomy(species_name: str, api_key: str | None = None) -> dict:
    """Fetch full taxonomy lineage from NCBI eutils for a species."""
    query = urllib.parse.quote(species_name.replace("_", " "))
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
           f"?db=taxonomy&term={query}[Scientific+Name]"
           f"&retmode=json&tool={NCBI_TOOL}&email={NCBI_EMAIL}")
    data = _ncbi_get(url, api_key)
    if not data:
        return {}
    try:
        d = json.loads(data)
        ids = d["esearchresult"]["idlist"]
        if not ids:
            url2 = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                    f"?db=taxonomy&term={query}"
                    f"&retmode=json&tool={NCBI_TOOL}&email={NCBI_EMAIL}")
            data = _ncbi_get(url2, api_key)
            if not data:
                return {}
            d = json.loads(data)
            ids = d["esearchresult"]["idlist"]
        if not ids:
            return {}
        taxid = ids[0]
    except Exception:
        return {}

    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=taxonomy&id={taxid}&retmode=xml"
           f"&tool={NCBI_TOOL}&email={NCBI_EMAIL}")
    data = _ncbi_get(url, api_key)
    if not data:
        return {}
    try:
        root = ET.fromstring(data)
        taxon = root.find("Taxon")
        if taxon is None:
            return {}

        lineage_ranks: dict[str, str] = {}
        lineage_ex = taxon.find("LineageEx")
        if lineage_ex:
            for t in lineage_ex.findall("Taxon"):
                r = t.findtext("Rank", "")
                n = t.findtext("ScientificName", "")
                if r in ("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"):
                    lineage_ranks[r] = n

        return {
            "taxid": taxid,
            "sci_name": taxon.findtext("ScientificName", ""),
            "rank": taxon.findtext("Rank", ""),
            "division": taxon.findtext("Division", ""),
            "lineage": lineage_ranks,
            "lineage_str": taxon.findtext("Lineage", ""),
        }
    except Exception as e:
        print(f"  XML parse error for {species_name}: {e}")
        return {}


def fetch_pubmed_snippets(
    species_name: str,
    api_key: str | None = None,
    max_papers: int = 5,
) -> list[str]:
    """Fetch phage receptor/susceptibility title snippets from PubMed."""
    genus = species_name.split("_")[0]
    sp = species_name.replace("_", " ")
    ids: list[str] = []
    for query_str in [f'"{sp}" phage receptor', f'"{genus}" bacteriophage receptor']:
        query = urllib.parse.quote(query_str)
        url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
               f"?db=pubmed&term={query}&retmax={max_papers}&sort=relevance"
               f"&retmode=json&tool={NCBI_TOOL}&email={NCBI_EMAIL}")
        data = _ncbi_get(url, api_key)
        if not data:
            continue
        try:
            d = json.loads(data)
            ids = d["esearchresult"]["idlist"]
            if ids:
                break
        except Exception:
            continue

    if not ids:
        return []

    url2 = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=pubmed&id={','.join(ids)}&rettype=abstract&retmode=xml"
            f"&tool={NCBI_TOOL}&email={NCBI_EMAIL}")
    data = _ncbi_get(url2, api_key)
    if not data:
        return []

    snippets: list[str] = []
    try:
        root = ET.fromstring(data)
        for article in root.findall(".//PubmedArticle"):
            title_el = article.find(".//ArticleTitle")
            if title_el is not None and title_el.text:
                snippets.append(f"• {title_el.text.strip()[:120]}")
    except Exception:
        pass
    return snippets[:max_papers]


# ── Infer Gram type ────────────────────────────────────────────────────────────
def infer_gram(host_id: str, tax: dict) -> str:
    genus = host_id.split("_")[0]
    if genus in GRAM_MAP:
        return {"+": "Gram-positive (+)", "-": "Gram-negative (-)",
                "archaea": "Archaea (no peptidoglycan)"}[GRAM_MAP[genus]]
    division = tax.get("division", "")
    if "Archaea" in division:
        return "Archaea (no peptidoglycan)"
    phylum = tax.get("lineage", {}).get("phylum", "")
    for p in ("Pseudomonadota", "Proteobacteria", "Bacteroidota",
              "Spirochaetota", "Fusobacteriota", "Campylobacterota",
              "Cyanobacteriota", "Chlorobiota"):
        if p in phylum:
            return "Gram-negative (-)"
    for p in ("Bacillota", "Actinomycetota", "Tenericutes",
              "Firmicutes", "Actinobacteria"):
        if p in phylum:
            return "Gram-positive (+)"
    return "Unknown"


# ── Build host profile ─────────────────────────────────────────────────────────
def build_profile(host_id: str, tax: dict, base: dict, pubmed: list[str]) -> str:
    """Assemble two-layer host profile (QUICK_PROFILE + PHAGE_INTERACTION_DETAILS).

    R1 cleanup is applied inline:
      - NCBI TaxID strings are not included in output
      - Key Literature section is not written
    """
    genus = host_id.split("_")[0]
    gram = infer_gram(host_id, tax)
    habitat = HABITAT_MAP.get(host_id, "")
    phage_fams = PHAGE_FAMILIES.get(host_id, PHAGE_FAMILIES.get(genus, ""))

    lin = tax.get("lineage", {})
    tax_parts = [lin[r] for r in ("phylum", "class", "order", "family") if r in lin]
    taxonomy_compact = " > ".join(tax_parts) if tax_parts else tax.get("lineage_str", "")[:80]

    defense = base.get("defense", {})
    defense_count = sum(defense.values())
    defense_compact = ("; ".join(f"{k}({v})" for k, v in sorted(defense.items()))
                       if defense else "none detected")

    curated_rec = KNOWN_RECEPTORS.get(host_id, "")
    receptor_str = curated_rec if curated_rec else "see surface gene annotations in base profile"

    intra        = INTRASPECIES_DISTINCTION.get(host_id, "")
    rec_note     = PHAGE_RECEPTOR_NOTES.get(host_id, "")
    cross        = CROSS_INFECTION_ALERT.get(host_id, "")
    phage_signals = PHAGE_RECOGNITION_GUIDE.get(host_id, [])

    lines: list[str] = []

    # ── LAYER 1: QUICK PROFILE ─────────────────────────────────────────────────
    lines.append(f"# HOST: {host_id}")
    lines.append("")
    lines.append("## QUICK_PROFILE")
    lines.append(f"GRAM: {gram}")
    if taxonomy_compact:
        lines.append(f"TAXONOMY: {taxonomy_compact}")
    if habitat:
        lines.append(f"ECOLOGY: {habitat}")
    lines.append(f"PRIMARY_RECEPTORS: {receptor_str}")
    lines.append(f"DEFENSE: {defense_count} systems — {defense_compact}")
    if phage_fams:
        lines.append(f"KNOWN_PHAGE_FAMILIES: {phage_fams}")
    if intra:
        lines.append(f"SPECIES_SPECIFICITY: {intra}")
    if rec_note:
        lines.append(f"ANNOTATION_WARNING: {rec_note}")
    if cross:
        lines.append(f"HOST_RANGE_BOUNDARY: {cross}")
    lines.append("")

    # ── LAYER 2: PHAGE INTERACTION DETAILS ────────────────────────────────────
    lines.append("## PHAGE_INTERACTION_DETAILS")
    lines.append("")
    lines.append(f"**Species:** *{host_id.replace('_', ' ')}*")
    if taxonomy_compact:
        lines.append(f"**Full taxonomy:** {taxonomy_compact}")
    lines.append("")

    lines.append("### Receptor Biochemistry")
    if curated_rec:
        lines.append(f"Confirmed receptors: {curated_rec}")
    else:
        lines.append("No curated receptor data; predicted from surface gene annotations.")
    if rec_note:
        lines.append(f"**Annotation note:** {rec_note}")
    lines.append("")

    if defense:
        lines.append("### Defense Systems (DefenseFinder)")
        lines.append(f"Total: {defense_count} — {defense_compact}")
        lines.append("")

    if phage_fams:
        lines.append("### Known Infecting Phage Families")
        lines.append(phage_fams)
        lines.append("")

    if intra:
        lines.append("### Species Specificity vs Closest Relatives")
        lines.append(intra)
        lines.append("")

    if cross:
        lines.append("### Cross-Infection Boundaries")
        lines.append(cross)
        lines.append("")

    if phage_signals:
        lines.append("### Phage Recognition Signals")
        for sig in phage_signals:
            lines.append(f"  {sig}")
        lines.append("")

    lines.append("---")
    return "\n".join(lines)


# ── Per-host pipeline ──────────────────────────────────────────────────────────
def process_host(
    host_id: str,
    base_dir: Path,
    out_dir: Path,
    api_key: str | None,
    skip_pubmed: bool,
    force: bool,
) -> str:
    out_path = out_dir / f"{host_id}.md"
    if out_path.exists() and not force:
        return f"[skip] {host_id}"

    base_path = base_dir / f"{host_id}.md"
    base = parse_base_profile(base_path)
    tax  = fetch_taxonomy(host_id, api_key)
    pubmed: list[str] = []
    if not skip_pubmed:
        pubmed = fetch_pubmed_snippets(host_id, api_key)

    profile = build_profile(host_id, tax, base, pubmed)
    out_path.write_text(profile)

    return (f"[done] {host_id} "
            f"(taxid={tax.get('taxid', '?')}, "
            f"defense={sum(base['defense'].values())}, "
            f"pubmed={len(pubmed)})")


# ── Main ───────────────────────────────────────────────────────────────────────
def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build NCBI-enriched host text profiles for LLM phage-host prediction "
            "(Stage 2: base profiles + NCBI taxonomy + curated knowledge → final profiles)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--base-profiles", required=True,
                        help="Directory of base host profiles from build_host_profiles_base.py")
    parser.add_argument("--out-dir",       required=True,
                        help="Output directory for enriched host profiles")
    parser.add_argument("--host-list",     default=None,
                        help="JSON file listing host IDs to process "
                             "(default: all hosts with base profiles in --base-profiles)")
    parser.add_argument("--ncbi-api-key",  default=None,
                        help="NCBI API key for higher rate limits (10/s vs 3/s). "
                             "Register at https://www.ncbi.nlm.nih.gov/account/")
    parser.add_argument("--threads",       type=int, default=1,
                        help="Parallel threads for NCBI fetching "
                             "(max 2 recommended without API key; max 5 with key)")
    parser.add_argument("--skip-pubmed",   action="store_true",
                        help="Skip PubMed literature fetch (faster; for offline use)")
    parser.add_argument("--force",         action="store_true",
                        help="Overwrite existing output profiles")
    parser.add_argument("--dry-run",       action="store_true",
                        help="Print preview of first 2 profiles without writing files")
    args = parser.parse_args()

    base_dir = Path(args.base_profiles)
    out_dir  = Path(args.out_dir)

    if args.host_list:
        host_ids: list[str] = json.loads(Path(args.host_list).read_text())
        print(f"Host list : {len(host_ids)} hosts from {args.host_list}")
    else:
        host_ids = sorted(p.stem for p in base_dir.glob("*.md"))
        print(f"Host list : {len(host_ids)} hosts discovered in {base_dir}")

    print(f"Base dir  : {base_dir}")
    print(f"Output    : {out_dir}")
    print(f"API key   : {'yes' if args.ncbi_api_key else 'no (≤3 req/s)'}")
    print(f"PubMed    : {'skip' if args.skip_pubmed else 'yes'}")
    print(f"Threads   : {args.threads}")
    print(f"Force     : {args.force}")
    print(f"Dry run   : {args.dry_run}")
    print()

    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        preview = 0
        for host_id in host_ids:
            if preview >= 2:
                break
            base = parse_base_profile(base_dir / f"{host_id}.md")
            tax  = fetch_taxonomy(host_id, args.ncbi_api_key)
            profile = build_profile(host_id, tax, base, [])
            print(f"  [DRY] {host_id}")
            print(profile[:600])
            print("  ...")
            preview += 1
        print(f"\n[DRY RUN] Would write {len(host_ids)} profiles.")
        return

    if args.threads > 1:
        max_workers = min(args.threads, 5 if args.ncbi_api_key else 2)
        with ThreadPoolExecutor(max_workers=max_workers) as exe:
            futures = {
                exe.submit(
                    process_host,
                    h, base_dir, out_dir,
                    args.ncbi_api_key, args.skip_pubmed, args.force,
                ): h
                for h in host_ids
            }
            for fut in as_completed(futures):
                print(" ", fut.result())
    else:
        for host_id in host_ids:
            result = process_host(
                host_id, base_dir, out_dir,
                args.ncbi_api_key, args.skip_pubmed, args.force,
            )
            print(" ", result)

    total = len(list(out_dir.glob("*.md")))
    print(f"\nDone. {total} host profiles written to: {out_dir}")


if __name__ == "__main__":
    main()
