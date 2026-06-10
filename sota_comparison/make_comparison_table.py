"""
make_comparison_table.py
========================
Assemble the recent-SOTA comparison on SAaIntDB under identical 10-fold splits (seed 42) and
identical frozen ESM-2 650M features. Reads the per-fold CSVs written by run_sota_cv.py and
combines them with our model / baselines (from RESULTS Table 1) into a markdown table.

Out -> sota_comparison/SOTA_comparison.md
"""
import os, glob, sys
import numpy as np, pandas as pd

HD = os.path.dirname(os.path.abspath(__file__))

# Ours + code-matched baselines (paper_reproducibility/RESULTS_updated.md, Table 1)
OURS = {
    ('random', 'Ours (All-CDR)'): (0.858, 0.006, 0.844, 0.009, 0.694, 0.014),
    ('random', 'Two-stream'):     (0.815, 0.027, 0.804, 0.028, 0.783, 0.055),
    ('random', 'Concat+MLP'):     (0.817, 0.031, 0.787, 0.046, 0.750, 0.063),
    ('cold',   'Ours (All-CDR)'): (0.568, 0.020, 0.551, 0.016, 1.102, 0.015),
    ('cold',   'Two-stream'):     (0.550, 0.107, 0.528, 0.077, 1.101, 0.071),
    ('cold',   'Concat+MLP'):     (0.549, 0.084, 0.520, 0.088, 1.445, 1.116),
}

def fmt(m, s):
    return f"{m:.3f} ± {s:.3f}"

def _agg(df):
    return (df.pearson.mean(), df.pearson.std(), df.spearman.mean(), df.spearman.std(),
            df.rmse.mean(), df.rmse.std(), len(df))

def load_sota():
    rows = {}
    for f in glob.glob(os.path.join(HD, 'results_*_*.csv')):
        name = os.path.basename(f)[len('results_'):-4]
        model, split = name.rsplit('_', 1)
        df = pd.read_csv(f)
        if len(df):
            rows[(split, model)] = _agg(df)
    # structure-based baseline (one file, 'split' column) — the ESMFold2-angle upper bound
    sb = os.path.join(HD, 'results_StructureBaseline.csv')
    if os.path.exists(sb):
        d = pd.read_csv(sb)
        for split, g in d.groupby('split'):
            rows[(split, 'Structure interface (experimental, RF)')] = _agg(g)
    return rows

def main():
    sota = load_sota()
    lines = ["# Recent-SOTA comparison on SAaIntDB (identical 10-fold splits, identical ESM-2 650M features)",
             "",
             "All models are trained and evaluated on the **same SAaIntDB rows**, the **same seed-42 "
             "10-fold random and cold (no-PDB-overlap) splits**, and the **same frozen ESM-2 650M** "
             "residue features — the only difference is the prediction architecture. MVSF-AB is already "
             "compared in RESULTS Table 4; LlamaAffinity is a binary binder *classifier* (accuracy/AUC) "
             "and is not directly comparable to pK_d regression, so it is omitted here.",
             ""]
    for split in ('random', 'cold'):
        lines += [f"## {split.capitalize()} split", "",
                  "| Model | Pearson r | Spearman ρ | RMSE | folds |", "|---|---|---|---|---|"]
        # ours/baselines first
        for (sp, model), v in OURS.items():
            if sp != split:
                continue
            star = "**" if "Ours" in model else ""
            lines.append(f"| {star}{model}{star} | {star}{fmt(v[0],v[1])}{star} | {fmt(v[2],v[3])} | {fmt(v[4],v[5])} | 10 |")
        # reimplemented SOTA + structure baseline
        for (sp, model), v in sorted(sota.items()):
            if sp != split:
                continue
            tag = "" if model.startswith("Structure") else " (reimpl.)"
            lines.append(f"| {model}{tag} | {fmt(v[0],v[1])} | {fmt(v[2],v[3])} | {fmt(v[4],v[5])} | {v[6]} |")
        lines.append("")
    out = os.path.join(HD, 'SOTA_comparison.md')
    open(out, 'w', encoding='utf-8').write("\n".join(lines) + "\n")
    print("wrote", out)
    sys.stdout.buffer.write(("\n".join(lines) + "\n").encode('utf-8'))

if __name__ == '__main__':
    main()
