# Recent-SOTA comparison (same dataset, same splits, same features)

A fair, apples-to-apples comparison of our model against **recent sequence-based SOTA** on
**SAaIntDB**, under the **same seed-42 10-fold random and cold splits** and the **same frozen
ESM-2 650M residue features** used by our model. The only thing that varies is the prediction
architecture.

## Which models, and why

| Model | Year | Included here | Reason |
|---|---|---|---|
| **MVSF-AB** | 2024 | no | already compared in RESULTS **Table 4** (SAbDab/AB-Bind/SKEMPI/Benchmark) |
| **DuaDeep-SeqAffinity** | 2025 | **yes** | sequence-only pK_d regression on ESM-2 650M (arXiv:2512.22007) |
| **DG-Affinity** | 2023/24 | **yes** | sequence-only pK_d regression, PLM + ConvNeXt family |
| **LlamaAffinity** | 2025 | no | binary binder **classifier** (accuracy/F1/AUC), not pK_d regression — not comparable |

No public code is available for DuaDeep-SeqAffinity / DG-Affinity, so they are **re-implemented
faithfully from the papers** (`sota_models.py`). To guarantee fairness, every model receives the
**identical frozen ESM-2 650M features** (`experiments/cache/perres_saaintdb_650M.pkl`) and the
**identical fold indices** (`get_fold_splits`, seed 42), and is trained with the **same budget**
(AdamW, MSE, ≤40 epochs, early stopping on validation Pearson).

## Files
- `precompute_perres_saaintdb.py` — per-residue ESM-2 650M for all unique SAaIntDB chains (one-time).
- `sota_models.py` — `DuaDeepSeqAffinity`, `DGAffinityConvNeXt` (paper-faithful; `reduced=True`
  projects 1280→256 with a 1-layer Transformer so it is trainable **on CPU**; `--full` uses the
  paper dimensions and is recommended on GPU).
- `run_sota_cv.py` — trains/evaluates a model under our exact folds; saves per-fold CSV incrementally.
- `make_comparison_table.py` — assembles `SOTA_comparison.md` from the CSVs + our Table-1 numbers.

## Run
```
# one-time features (CPU ~30 min, GPU ~2 min)
python sota_comparison/precompute_perres_saaintdb.py

# full matrix (GPU strongly recommended; CPU is hours per fold)
python sota_comparison/run_sota_cv.py --model DuaDeep-SeqAffinity --split random --full
python sota_comparison/run_sota_cv.py --model DuaDeep-SeqAffinity --split cold   --full
python sota_comparison/run_sota_cv.py --model DG-Affinity         --split random --full
python sota_comparison/run_sota_cv.py --model DG-Affinity         --split cold   --full

python sota_comparison/make_comparison_table.py   # -> SOTA_comparison.md
```

## Compute note
This box has **no GPU**; ESM-2 feature precompute (~30 min) is done and cached, but full-CV training
of the per-residue Transformer/ConvNeXt baselines is **hours per fold on CPU**. The matrix is launched
with the CPU-trainable `reduced` config and writes results incrementally; on a GPU the same commands
finish in minutes per fold with the `--full` (paper-dimension) config. `make_comparison_table.py` can
be re-run at any time to refresh `SOTA_comparison.md` from whatever folds have completed.
