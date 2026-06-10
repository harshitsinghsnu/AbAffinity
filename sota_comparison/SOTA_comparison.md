# Recent-SOTA comparison on SAaIntDB (identical 10-fold splits, identical ESM-2 650M features)

All models are trained and evaluated on the **same SAaIntDB rows**, the **same seed-42 10-fold random and cold (no-PDB-overlap) splits**, and the **same frozen ESM-2 650M** residue features — the only difference is the prediction architecture. MVSF-AB is already compared in RESULTS Table 4; LlamaAffinity is a binary binder *classifier* (accuracy/AUC) and is not directly comparable to pK_d regression, so it is omitted here.

## Random split

| Model | Pearson r | Spearman ρ | RMSE | folds |
|---|---|---|---|---|
| **Ours (All-CDR)** | **0.858 ± 0.006** | 0.844 ± 0.009 | 0.694 ± 0.014 | 10 |
| Two-stream | 0.815 ± 0.027 | 0.804 ± 0.028 | 0.783 ± 0.055 | 10 |
| Concat+MLP | 0.817 ± 0.031 | 0.787 ± 0.046 | 0.750 ± 0.063 | 10 |
| DuaDeep-SeqAffinity (reimpl.) | 0.764 ± 0.040 | 0.741 ± 0.047 | 0.851 ± 0.094 | 6 |
| Structure interface (experimental, RF) | 0.489 ± 0.076 | 0.442 ± 0.073 | 1.094 ± 0.067 | 10 |

## Cold split

| Model | Pearson r | Spearman ρ | RMSE | folds |
|---|---|---|---|---|
| **Ours (All-CDR)** | **0.568 ± 0.020** | 0.551 ± 0.016 | 1.102 ± 0.015 | 10 |
| Two-stream | 0.550 ± 0.107 | 0.528 ± 0.077 | 1.101 ± 0.071 | 10 |
| Concat+MLP | 0.549 ± 0.084 | 0.520 ± 0.088 | 1.445 ± 1.116 | 10 |
| Structure interface (experimental, RF) | 0.267 ± 0.091 | 0.241 ± 0.064 | 1.222 ± 0.157 | 10 |

