# PHI-Reason-Species

![PHI-Reason-Species Poster](poster.png)

Code for **PHI-Reason-Species**: an LLM-based phage–host interaction (PHI) prediction framework that integrates multi-channel genomic evidence (RBP annotations, BLASTN neighbours, BLASTP neighbours) into structured text profiles and queries a reasoning LLM to predict phage host range.

---

## Repository structure

```
PHIreason_pipeline/
├── config.sh                          # Path configuration — edit before first run
├── requirements.txt                   # Python dependencies
│
├── 01_profile_generation/             # Build phage and host profiles
│   ├── phage_profiles/
│   │   └── build_phage_profiles.py    # Generate phage text profiles (RBP + BLASTN)
│   ├── host_profiles/
│   │   ├── build_host_profiles.py     # Generate host text profiles from NCBI annotations
│   │   └── build_host_list.py         # Compile candidate host list for the prompt
│   └── blastp_context/                # BLASTP phylogenetic context (run in order)
│       ├── 01_build_db.py             # Build DIAMOND protein database (CHERRY-1306)
│       ├── 02_run_query.sh            # Run diamond blastp: test phages vs. reference DB
│       ├── 03_extract_neighbors.py    # Extract top-5 protein-similarity neighbours
│       ├── 04_build_context_blocks.py # Format neighbours into LLM-ready context blocks
│       └── 05_inject_context.py       # Inject BLASTP blocks into phage profiles
│
├── 02_inference/                      # LLM inference
│   ├── run_inference.py               # Main inference driver (async, concurrent, resumable)
│   └── prompt/
│       └── prompt_template.py         # System prompt and user message templates
│
├── 03_experiments/
│   └── run_experiments.sh             # Run the four core experiments sequentially
│
└── 04_baseline/
    └── evaluate_baseline.py           # Evaluate baseline tools (PHP, WIsH, PHIST, etc.)
```

---

## Requirements

### System dependencies

| Tool | Purpose | Install |
|------|---------|---------|
| [Ollama](https://ollama.com) | Serve local LLM | See ollama.com |
| [DIAMOND](https://github.com/bbuchfink/diamond) | Protein similarity search | `conda install -c bioconda diamond` |
| [Prodigal](https://github.com/hyattpd/Prodigal) | Gene prediction (VHDB only) | `conda install -c bioconda prodigal` |

### Python dependencies

```bash
pip install -r requirements.txt
```

Requires Python 3.10+.

---

## Setup

**1. Configure project root**

Edit `config.sh` or set the environment variable before running any script:

```bash
export PHI_PROJECT_ROOT=/path/to/your/project
```

`PHI_PROJECT_ROOT` should be the directory containing `ws/`, `data/`, and `experiments/`. If `PHIreason_pipeline/` is cloned directly inside that directory, the path is auto-detected.

**2. (Optional) Override tool paths**

```bash
export DIAMOND_BIN=/path/to/diamond    # default: looked up from PATH
export PRODIGAL_BIN=/path/to/prodigal  # default: looked up from PATH
```

---

## Pipeline walkthrough

### Step 1 — Build host profiles and host list

```bash
# Generate host text profiles from NCBI genome annotations
python 01_profile_generation/host_profiles/build_host_profiles.py \
    --genome-dir  ${PHI_PROJECT_ROOT}/data/host_genomes \
    --out-dir     ${PHI_PROJECT_ROOT}/ws/textGeneProfile_v3_ncbi_R1

# Compile the candidate host list used in the LLM prompt
python 01_profile_generation/host_profiles/build_host_list.py \
    --host-prof   ${PHI_PROJECT_ROOT}/ws/textGeneProfile_v3_ncbi_R1 \
    --out          ws/Cherry/host_list_v4G.txt
```

### Step 2 — Build phage profiles

```bash
python 01_profile_generation/phage_profiles/build_phage_profiles.py \
    --phage-dir  ${PHI_PROJECT_ROOT}/data/Cherry_data/phage \
    --host-list  ws/Cherry/host_list_v4G.txt \
    --out-dir    ws/Cherry/phage_profiles_rbp_blastn
```

### Step 3 — Add BLASTP context (optional, improves accuracy)

Run the five scripts in order from `01_profile_generation/blastp_context/`:

```bash
# 1. Build protein reference database
python 01_build_db.py --out_dir ${PHI_PROJECT_ROOT}/experiments/blastp_cherry1306/db

# 2. Query test phage proteins against the reference DB
bash 02_run_query.sh cherry634

# 3–5. Extract neighbours, build context blocks, inject into profiles
python 03_extract_neighbors.py
python 04_build_context_blocks.py
python 05_inject_context.py \
    --src_dir     ws/Cherry/phage_profiles_rbp_blastn \
    --context_json experiments/blastp_cherry1306/outputs/blastp_context_cherry634.json \
    --out_dir     ws/Cherry/phage_profiles_rbp_blastn_blastp
```

### Step 4 — Run inference

```bash
export PHI_PROJECT_ROOT=/path/to/project

python 02_inference/run_inference.py \
    --exp-id       my_experiment \
    --host-list    ws/Cherry/host_list_v4G.txt \
    --pair-csv     data/Cherry_data/phage1940_host_pair.csv \
    --phage-list   data/Cherry_data/cherry_test634_list.txt \
    --phage-prof-dir ws/Cherry/phage_profiles_rbp_blastn_blastp \
    --model        qwen3-coder-next:q4_K_M \
    --concurrency  12 \
    --num-ctx      40960 \
    --num-predict  4096
```

Results are written to `experiments/my_experiment/results/metrics.json`. Inference is resumable — re-running the same `--exp-id` skips already-completed phages.

**Key inference options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--model` | `qwen3-coder-next:q4_K_M` | Ollama model tag |
| `--ollama-urls` | `http://127.0.0.1:11435,...` | Comma-separated Ollama endpoints (for multi-GPU) |
| `--concurrency` | `12` | Number of parallel requests |
| `--think` | — | `yes` / `no` / `auto` — enable chain-of-thought (if model supports it) |
| `--pilot` | — | Run on first N phages only (for testing) |
| `--overwrite` | — | Re-run even if cached result exists |

---


## Baseline evaluation

See [`04_baseline/README.md`](04_baseline/README.md) for step-by-step instructions on installing each baseline tool, running it on the Cherry RefSeq-634 test set, and computing metrics with `evaluate_baseline.py`.

Quick example:

```bash
# 1. Clone the shared benchmark repository (provides metric code and test split)
git clone https://github.com/KennthShang/HostPredictionReview.git
export HOST_PREDICTION_REVIEW_DIR=/path/to/HostPredictionReview

# 2. After running a baseline tool, evaluate its output
python 04_baseline/evaluate_baseline.py \
    --tool       wish \
    --outputs    experiments/baseline_wish/outputs \
    --review-dir ${HOST_PREDICTION_REVIEW_DIR}
```

Supported tools: `php`, `wish`, `phist`, `deephost`, `phabox2`, `vhmnet`.

---

## Citation

> (manuscript in preparation)
