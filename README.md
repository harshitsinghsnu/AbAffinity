# AgAbGated — Gated Tri-Stream Antibody–Antigen Binding-Affinity Prediction

A sequence-only model that predicts antibody–antigen binding affinity (pK_d) directly from the heavy,
light and antigen sequences. Frozen **ESM-2 (650M)** embeddings are pooled (heavy chain by its
**CDR-H1/H2/H3 loops**, light and antigen by mean), kept as three separate streams, coupled by a
heavy↔light fusion gate and a **gated cross-attention** interaction, and read out by a bounded cosine
head. A single model spans conventional paired antibodies and single-domain nanobodies.

## Headline results (SAaIntDB, 2,575 complexes)

| Setting | Pearson r | Spearman ρ | RMSE (pK_d) |
|---|---|---|---|
| Random 10-fold CV | **0.858 ± 0.006** | 0.844 | 0.694 |
| Cold (no-PDB-overlap) 10-fold | **0.568 ± 0.020** | 0.551 | 1.102 |

On identical data, splits and ESM-2 features, the model leads the code-matched baselines
(two-stream, concat+MLP), the re-implemented recent SOTA (DuaDeep-SeqAffinity), and a structure-based
interface baseline computed on the *experimental* crystal structures — see
[`sota_comparison/SOTA_comparison.md`](sota_comparison/SOTA_comparison.md) and
[`docs/ESMFOLD2_positioning.md`](docs/ESMFOLD2_positioning.md). It also exceeds MVSF-AB on the external
SAbDab / AB-Bind / SKEMPI / benchmark sets (`docs/RESULTS.md`, Table 4).

## Repository layout

```
AgAbGated/
├── agabgated/              # the Python package
│   ├── models/             # gated tri-stream model + two-stream / mutual-strong variants
│   ├── utils/              # data loading, splits, metrics, cached-embedding dataset
│   ├── training/           # CV drivers, ablations, baselines, multi-seed runners
│   ├── explain/            # integrated-gradients attribution + gate-importance
│   ├── data_prep/          # ESM-2 embedding precompute / caching
│   └── figures/            # notebook + figure builders
├── configs/                # YAML experiment configs (benchmark/, saaintdb/)
├── notebooks/              # ★ the two publication notebooks (run from these)
├── docs/                   # METHODS, RESULTS (+ tables), FIGURE_CAPTIONS, ESMFOLD2 positioning
├── data/                   # small CSVs (SAaIntDB table, pairs_*), complexes.json
├── experiments/            # result artifacts the notebooks read (explainability npz, stats)
├── results_saaintdb_allcdr/, advanced_results/   # more notebook artifacts
├── figures/                # generated publication figures (panels, explainability)
├── sota_comparison/        # recent-SOTA + structure-based baseline pipeline
├── pymol_groundtruth/      # structural interface figures (PyMOL) + composites
├── requirements.txt / environment.yml / pyproject.toml
```

## Installation

```bash
conda create -n agabgated python=3.10 && conda activate agabgated
pip install -r requirements.txt
pip install -e .            # makes `import agabgated` available
```
(Install the `torch` build matching your CUDA. The notebooks need only the first five packages in
`requirements.txt`; `torch/transformers/captum/biopython/ANARCI` are for full reproduction.)

## Quick start — reproduce the publication figures (no GPU, minutes)

Both key notebooks are **self-contained**: they read the small result artifacts shipped in this repo
and need only numpy/pandas/scipy/matplotlib. Launch Jupyter anywhere in the repo:

```bash
jupyter lab     # then open notebooks/
```

- **`notebooks/BALM_AbAg_Panels.ipynb`** — the seven Nature-style result panels (architecture,
  antibody-vs-nanobody, MVSF/overall-vs-intra, transfer, gating + stream interventions + per-CDR +
  top-K, per-residue mutagenesis heatmaps, best-fold + benchmark scatter). Numbers are embedded in the
  cells; only the scatter/heatmap panels read CSV/NPZ from `results_saaintdb_allcdr/`,
  `experiments/results_allcdr_stats/` and `advanced_results/`.
- **`notebooks/BALM_AbAg_Explainability_Analysis.ipynb`** — integrated-gradients attribution comparing
  the final All-CDR model vs the two-stream baseline on five canonical complexes
  (1VFB, 3HFM, 5GRJ, 5Y9J, 4ETQ): CDR-enrichment, heavy-chain attribution profiles, and per-residue
  importance heatmaps. Reads `experiments/results_explainability/*.npz` + `cdr_enrichment.csv` and the
  bundled `data/complexes.json`.

Every cell runs top-to-bottom and regenerates the figures into `figures/`.

## Full reproduction from scratch (GPU recommended)

```bash
# 1. cache ESM-2 650M embeddings for SAaIntDB
python -m agabgated.data_prep.precompute_embeddings_saaintdb
python -m agabgated.data_prep.precompute_allcdr_cache        # All-CDR pooled features

# 2. cross-validation of the model + baselines (random + cold)
python -m agabgated.training.run_saaintdb_multiseed --config configs/saaintdb/sa_ours_allcdr_random.yaml
python -m agabgated.training.run_saaintdb_multiseed --config configs/saaintdb/sa_ours_allcdr_cold.yaml

# 3. integrated-gradients attribution for the explainability notebook
python -m agabgated.explain.export_ig_attributions

# 4. rebuild the notebooks / figures
python -m agabgated.figures.build_panels_notebook
python -m agabgated.figures.build_explainability_notebook
```

Configs in `configs/` define every ablation (pooling, architecture, PLM backbone, gating, antigen
dependence). Full methodology and the complete results tables are in `docs/METHODS.md` and
`docs/RESULTS.md`.

## Recent-SOTA & structure-based comparison

`sota_comparison/` retrains recent sequence-based SOTA (DuaDeep-SeqAffinity, DG-Affinity — no public
code, re-implemented from the papers) on the **same SAaIntDB rows, same seed-42 10-fold splits and
same frozen ESM-2 650M features**, and computes a structure-based interface baseline on the deposited
experimental structures (the ESMFold2 angle — an upper bound on any predicted-structure pipeline). See
`sota_comparison/README.md`.

## Data availability

Small tables (the SAaIntDB affinity table, pair lists, the five explainability sequences) ship in
`data/`. Large ESM-2 embedding caches and fetched PDB structures are **not** committed (see
`.gitignore`); regenerate them with `agabgated/data_prep/` and `sota_comparison/precompute_*` /
`structure_baseline.py`. SAaIntDB and the external datasets are derived from their original public
releases.

## Citation

If you use this code, please cite the accompanying manuscript (in preparation). Model and analysis by
**Harshit Singh**. Built on ESM-2 (Lin et al., 2023) and captum integrated gradients.
