# Baseline Evaluation

This directory contains `evaluate_baseline.py`, a unified evaluator that applies consistent metric definitions across six baseline phage–host prediction tools.

---

## Step 1 — Clone HostPredictionReview

All metric computations (`sp@1`, `g@1`, `sp@5`, `MRR`, `AUC`, `AUPR`) are delegated to `evaluate_full.py` from the [HostPredictionReview](https://github.com/KennthShang/HostPredictionReview) repository, which provides the shared benchmark split and metric code.

```bash
git clone https://github.com/KennthShang/HostPredictionReview.git
export HOST_PREDICTION_REVIEW_DIR=/path/to/HostPredictionReview
```

The repository provides two items used by `evaluate_baseline.py`:
- `CHERRY_benchmark_datasplit/CHERRY_test.fasta` — test phage IDs
- `CHERRY_benchmark_datasplit/CHERRY_y_test.csv` — ground-truth host labels
- `CHERRY_benchmark_results/evaluate_full.py` — metric computation code

---

## Step 2 — Install and run each baseline tool

Run each tool on the **Cherry RefSeq-634** test set. Then call `evaluate_baseline.py` on the output directory.

---

### PHP

**Install:** https://github.com/dengzq1234/PHP

```bash
# Run PHP on Cherry test phages
python php.py \
    --phage  data/Cherry_data/phage_test.fasta \
    --host   data/Cherry_data/host.fasta \
    --kmer   4 \
    --out    experiments/baseline_php/outputs/

# Evaluate
python 04_baseline/evaluate_baseline.py \
    --tool    php \
    --outputs experiments/baseline_php/outputs \
    --review-dir ${HOST_PREDICTION_REVIEW_DIR}
```

Expected output file: `outputs/cherry_host_kmer4_Prediction_Allhost.csv`

---

### WIsH

**Install:** https://github.com/soedinglab/WIsH

```bash
# Build host models
WIsH -c build -g data/Cherry_data/host/ -m wish_models/

# Predict
WIsH -c predict -g data/Cherry_data/phage_test/ -m wish_models/ \
    -r experiments/baseline_wish/outputs/predictions/

# Evaluate
python 04_baseline/evaluate_baseline.py \
    --tool    wish \
    --outputs experiments/baseline_wish/outputs \
    --review-dir ${HOST_PREDICTION_REVIEW_DIR}
```

Expected output file: `outputs/predictions/llikelihood.matrix`

---

### PHIST

**Install:** https://github.com/refresh-bio/PHIST

```bash
# Run PHIST
phist.py data/Cherry_data/phage_test/ data/Cherry_data/host/ \
    --out experiments/baseline_phist/outputs/predictions.csv

# Evaluate
python 04_baseline/evaluate_baseline.py \
    --tool    phist \
    --outputs experiments/baseline_phist/outputs \
    --review-dir ${HOST_PREDICTION_REVIEW_DIR}
```

Expected output file: `outputs/predictions.csv` (columns: `phage, host, common kmers`)

---

### DeepHost

**Install:** https://github.com/deepomicslab/DeepHost

```bash
# Run DeepHost
python predict.py \
    --phage data/Cherry_data/phage_test.fasta \
    --host  data/Cherry_data/host.fasta \
    --out   experiments/baseline_deephost/outputs/

# Evaluate
python 04_baseline/evaluate_baseline.py \
    --tool    deephost \
    --outputs experiments/baseline_deephost/outputs \
    --review-dir ${HOST_PREDICTION_REVIEW_DIR}
```

Expected output file: `outputs/DeepHost_multihost.tsv`

---

### PhaBox2

**Install:** https://github.com/KennthShang/PhaBOX

```bash
# Run PhaBox2 (host prediction module)
python main.py \
    --contigs  data/Cherry_data/phage_test.fasta \
    --dbdir    /path/to/phabox2_db \
    --outpth   experiments/baseline_phabox2/outputs/ \
    --task     cherry

# Evaluate (db mode or mag mode)
python 04_baseline/evaluate_baseline.py \
    --tool    phabox2 \
    --mode    db \
    --outputs experiments/baseline_phabox2/outputs \
    --review-dir ${HOST_PREDICTION_REVIEW_DIR}
```

Expected output file: `outputs/results/final_prediction/cherry_prediction.tsv`

---

### VHMnet

**Install:** https://github.com/yolandalalala/VHMnet

```bash
# Run VHMnet
python predict.py \
    --phage data/Cherry_data/phage_test.fasta \
    --host  data/Cherry_data/host.fasta \
    --out   experiments/baseline_vhmnet/outputs/

# Evaluate
python 04_baseline/evaluate_baseline.py \
    --tool    vhmnet \
    --outputs experiments/baseline_vhmnet/outputs \
    --review-dir ${HOST_PREDICTION_REVIEW_DIR}
```

Expected output directory: `outputs/output/predictions/` or `outputs/predictions/`

---

## Output

Each evaluation call writes `results/METRICS.json` alongside the `outputs/` directory and prints a summary:

```json
{
  "tool": "WIsH",
  "sp@1": 0.213,
  "g@1":  0.341,
  "sp@5": 0.398,
  "MRR":  0.287,
  "AUC":  0.671,
  "AUPR": 0.244
}
```

## Environment variable reference

| Variable | Description |
|----------|-------------|
| `HOST_PREDICTION_REVIEW_DIR` | Path to cloned HostPredictionReview repository |
| `PHI_PROJECT_ROOT` | PHIreason project root (see `../config.sh`) |
