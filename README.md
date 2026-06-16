# AgAbGated — Gated Tri-Stream Antibody–Antigen Binding-Affinity Prediction

A sequence-only model that predicts antibody–antigen binding affinity (pK_d) directly from heavy,
light and antigen sequences. Frozen **ESM-2 (650M)** embeddings are pooled (heavy chain over its
**CDR-H1/H2/H3** loops, light and antigen by mean), kept as three separate streams, coupled by a
heavy↔light **fusion gate** and a **gated cross-attention** interaction, and read out by a bounded
**cosine** head. One model spans conventional paired antibodies and single-domain nanobodies.

## Headline results

| Dataset | Setting | Seeds | Pearson r | Spearman ρ | RMSE (pK_d) |
|---|---|---|---|---|---|
| **SAaIntDB** (2,575) | random 10-fold CV | 3 | **0.858 ± 0.006** | 0.844 ± 0.009 | 0.694 ± 0.014 |
| **SAaIntDB** (2,575) | cold (no-PDB-overlap) 10-fold | 3 | **0.568 ± 0.020** | 0.551 ± 0.016 | 1.102 ± 0.015 |
| **SAbDab** | 10-fold CV | 5 | **0.593 ± 0.015** | 0.582 ± 0.014 | 1.326 ± 0.018 |
| **AB-Bind** | 10-fold CV | 5 | **0.781 ± 0.006** | 0.753 ± 0.008 | 1.288 ± 0.009 |
| **SKEMPI 2.0** | 10-fold CV | 5 | **0.722 ± 0.010** | 0.656 ± 0.015 | 1.021 ± 0.015 |
| **Benchmark** | held-out (train SAbDab) | 5 | **0.551 ± 0.056** | 0.547 ± 0.062 | 1.339 ± 0.075 |

**Seeds.** SAaIntDB is summarised over **3 random seeds** (`42 114 144`); the external MVSF-AB
datasets over **5 seeds** (`42 114 144 314 777`) — these are the seed counts used in the paper, and
they are the defaults of the respective runners (see below).

## Repository layout

```
AgAbGated/
├── agabgated/                # the importable Python package
│   ├── models/               # MutualTriStreamStrong (final) + variants
│   ├── utils/                # data loading, metrics, training loops
│   ├── training/             # CV / benchmark / multi-seed runners, ablations
│   ├── explain/              # integrated-gradients attribution + gate importance
│   ├── data_prep/            # ESM-2 embedding precompute / All-CDR caching
│   └── figures/              # notebook / figure builders
├── configs/
│   ├── saaintdb/             # SAaIntDB configs (All-CDR / mean-pool × random / cold)
│   └── benchmark/            # external-dataset configs (CV + SAbDab→benchmark)
├── data/                     # pair tables, SAaIntDB rows, complexes.json, AND prebuilt caches:
│   ├── allcdr_natural_650M.pkl          # SAbDab + benchmark, All-CDR pooled
│   ├── allcdr_mutation_650M.pkl         # AB-Bind + SKEMPI, All-CDR pooled
│   ├── esm2_embeddings_saaintdb_650M.pkl
│   └── saaintdb_heavy_cdr_embeddings.pkl
├── model_weights/
│   └── saaintdb_allcdr_random_bestfold.pt   # final model, SAaIntDB random best fold
├── notebooks/                # custom run + publication panels + explainability
├── environment.yml / requirements.txt / pyproject.toml
```

The prebuilt embedding caches and the final checkpoint **are committed**, so both the notebook and
the cross-validation runners work without any ESM-2 recompute. (Large per-residue caches and fold
outputs are not committed — see `.gitignore`.)

## Installation

```bash
conda env create -f environment.yml      # exact pins used for every reported number
conda activate agabgated
pip install -e .                          # makes `import agabgated` available
```

---

## Quick start — run the model on **your own** data

Open **`notebooks/AgAbGated_Custom_ZeroShot_FewShot_IG.ipynb`**. It loads
`model_weights/saaintdb_allcdr_random_bestfold.pt` and, from a CSV of your complexes
(`light, heavy, antigen, Y`), runs end to end:

1. **Zero-shot** affinity prediction (frozen model) → Pearson / Spearman / RMSE,
2. **Few-shot** fine-tuning on a fraction of your labels (validation early-stopping, held-out test),
3. **Integrated Gradients** per-residue attribution → heatmap,
4. **Structure mapping** — top-attribution residues highlighted on the 3D structure, rendered inline.

This is the recommended entry point and needs only the bundled checkpoint.

---

## Reproduce the paper results

The committed caches in `data/` mean these run directly — **no precompute step required**.

### SAaIntDB — headline (r = 0.858 random / 0.568 cold), 3 seeds

```bash
python -m agabgated.training.run_saaintdb_multiseed --config configs/saaintdb/sa_ours_allcdr_random.yaml
python -m agabgated.training.run_saaintdb_multiseed --config configs/saaintdb/sa_ours_allcdr_cold.yaml
```
Uses `data/saaintdb_with_antigen_names.csv` + the SAaIntDB caches; default seeds `42 114 144`.
Results → `results/results_saaintdb/<name>/aggregated_summary.csv`.

### External MVSF-AB datasets — SAbDab / AB-Bind / SKEMPI / benchmark, 5 seeds

The same final model (MutualTriStreamStrong, All-CDR pooling), driven by `configs/benchmark/`:

```bash
# 10-fold CV on SAbDab, AB-Bind and SKEMPI:
python -m agabgated.training.run_multiseed --config configs/benchmark/exp02_ours_allcdr_cv.yaml

# Train on 100% SAbDab, evaluate on the held-out benchmark:
python -m agabgated.training.run_multiseed --config configs/benchmark/exp04_ours_allcdr_benchmark.yaml

# Or every benchmark config at once:
python -m agabgated.training.run_multiseed --all
```
Uses `data/allcdr_natural_650M.pkl` and `data/allcdr_mutation_650M.pkl`; default seeds
`42 114 144 314 777`. Results → `results/<name>/aggregated_summary.csv` + `results/MASTER_SUMMARY.csv`.

> To rebuild the caches from scratch (e.g. on new data), run
> `python -m agabgated.data_prep.precompute_embeddings_saaintdb` and
> `python -m agabgated.data_prep.precompute_allcdr_cache` first; this requires the large per-residue
> ESM-2 caches.

---

## Notebooks

- **`AgAbGated_Custom_ZeroShot_FewShot_IG.ipynb`** — run the model on a custom CSV: zero-shot,
  few-shot, IG heatmap, inline structure mapping (above).
- **`BALM_AbAg_Panels.ipynb`** — the publication result panels (self-contained).
- **`BALM_AbAg_Explainability_Analysis.ipynb`** — integrated-gradients attribution on five canonical
  complexes (1VFB, 3HFM, 5GRJ, 5Y9J, 4ETQ).

## Data format

All pair tables and the custom-notebook CSV use four columns: `light`, `heavy`, `antigen`, `Y`
(affinity in pK_d). For nanobodies, place the heavy sequence in the `light` column too.



