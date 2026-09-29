"""
benchmark_fixed_epoch_matched.py
==================================
Protocol-matched comparison with MVSF-AB (Li et al., 2024): their train.py
trains for a fixed --num-epochs (default 30) with NO early stopping and NO
checkpoint selection on any held-out signal (verified directly against
https://github.com/TAI-Medical-Lab/MVSF-AB/blob/main/train.py -- the loop is
`for epoch in range(num_epochs)` with metrics logged/reported at epoch 30,
no `best`/`patience` tracking, model-saving code commented out).

This trains AbAffinity (All-CDR, cosine head) on the full SAbDab set for
exactly 30 epochs with no early stopping and no validation-based checkpoint
selection, then evaluates once on the held-out Benchmark set at that final
epoch -- the same protocol MVSF-AB uses, so training-procedure asymmetry is
removed from the comparison. Run across 3 seeds (9999, 114, 144) for
consistency with the rest of the revision.

This is the headline "MVSF-AB-matched protocol" Benchmark result reported in
the manuscript and rebuttal: Pearson r = 0.564 +- 0.026, Spearman
rho = 0.570 +- 0.033, RMSE = 1.364 +- 0.116 pKd.
"""
import os, sys, pickle, gc
import numpy as np, pandas as pd, torch
from scipy.stats import pearsonr, spearmanr
import warnings; warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from AbAffinity.models.mutual_strong import MutualTriStreamStrong
from AbAffinity.utils.main_symmetric_mean import load_data, CachedEmbeddingDataset, collate_fn, DEFAULT_CONFIG
from torch.utils.data import DataLoader
from Bio import SeqIO

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
OUT = os.path.join(HERE, 'results_fixed_epoch_matched')
os.makedirs(OUT, exist_ok=True)
print(f"Device: {DEVICE} | Output: {OUT}")

IMGT_H1 = set(range(27, 39)); IMGT_H2 = set(range(56, 66)); IMGT_H3 = set(range(105, 118))

def cdr_positions(seq):
    try:
        import anarci
        res = anarci.anarci([('q', seq)], scheme='imgt', assign_germline=False,
                             allow=set('ACDEFGHIKLMNPQRSTVWY'))
        if res and res[0] and res[0][0] and res[0][0][0] is not None:
            numbered = res[0][0][0]
            idxs, sp = [], 0
            for (pos, _ins), aa in numbered:
                if aa == '-':
                    continue
                if pos in IMGT_H1 or pos in IMGT_H2 or pos in IMGT_H3:
                    idxs.append(sp)
                sp += 1
            if idxs:
                return sorted(set(idxs))
    except Exception:
        pass
    return _regex_cdr(seq)

def _regex_cdr(seq):
    import re
    out = []
    m1 = re.search(r'C[A-Z]{2,5}[SAGTV]([A-Z]{8,15})WVRQ', seq)
    out += list(range(*m1.span(1))) if m1 else list(range(26, 36))
    m2 = re.search(r'W[VI]RQ[A-Z]{6,14}W[VL][AS]([A-Z]{10,20})VKGRF', seq)
    out += list(range(*m2.span(1))) if m2 else list(range(50, 64))
    m3 = re.search(r'WYYCA([A-Z]+)WGQGT', seq) or re.search(r'WYYC[A-Z]([A-Z]+)WGQG', seq)
    out += list(range(*m3.span(1))) if m3 else list(range(95, 110))
    return sorted(set(out))

def build_combined_emb(seqs, per_res, mean_emb, heavy_ids):
    combined = dict(mean_emb)
    for hid in heavy_ids:
        if hid not in per_res or hid not in seqs:
            continue
        pr = per_res[hid]
        seq = seqs[hid]
        idxs = [i for i in cdr_positions(seq) if 0 <= i < pr.shape[0]]
        combined[hid] = (pr[idxs, :].mean(0) if idxs else pr.mean(0)).astype(np.float32)
    return combined

class DictLoader:
    def __init__(self, d):
        self.embeddings = d
        self.embedding_dim = next(iter(d.values())).shape[0]
    def get_embedding(self, k):
        return self.embeddings[k]


print("[1] Loading sequences + natural embeddings...")
seqs = {}
for fn in ['data/seq_natural.fasta', 'data/seq.fasta']:
    fp = os.path.join(HERE, fn)
    if os.path.exists(fp):
        for r in SeqIO.parse(fp, 'fasta'):
            seqs.setdefault(r.id, str(r.seq))

with open(os.path.join(HERE, 'data/esm2_per_residue_embeddings_natural_650M.pkl'), 'rb') as f:
    pr_nat = pickle.load(f)
with open(os.path.join(HERE, 'data/esm2_embeddings_natural_650M.pkl'), 'rb') as f:
    mean_nat = pickle.load(f)

df_sab = load_data(os.path.join(HERE, 'data/pairs_sabdab.csv'))
df_bench = load_data(os.path.join(HERE, 'data/pairs_benchmark.csv'))
heavy_ids_nat = set(df_sab['heavy_id']) | set(df_bench['heavy_id'])
comb_nat = build_combined_emb(seqs, pr_nat, mean_nat, heavy_ids_nat)
loader = DictLoader(comb_nat)
del pr_nat, mean_nat; gc.collect()

idc = ['heavy_id', 'light_id', 'antigen_id']
keep = df_sab[idc].apply(lambda c: c.isin(loader.embeddings)).all(axis=1)
df_tr_full = df_sab[keep].reset_index(drop=True)
keepb = df_bench[idc].apply(lambda c: c.isin(loader.embeddings)).all(axis=1)
df_b = df_bench[keepb].reset_index(drop=True)
print(f"  SAbDab train n={len(df_tr_full)}  Benchmark test n={len(df_b)}")


def train_fixed_epochs_no_earlystop(df_train, num_epochs, seed, batch_size=16, lr=1e-4, weight_decay=1e-5):
    """MVSF-AB-matched: fixed epoch count, NO early stopping, NO validation-based
    checkpoint selection -- train on ALL of df_train, return the model exactly
    as it stands after the final epoch."""
    torch.manual_seed(seed); np.random.seed(seed)
    lo, hi = df_train['binding_affinity'].min(), df_train['binding_affinity'].max()
    model = MutualTriStreamStrong(esm_dim=loader.embedding_dim, projected_size=256,
                                   num_heads=8, dropout=0.1, n_layers=2, device=DEVICE).to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    dl = DataLoader(CachedEmbeddingDataset(df_train, loader), batch_size=batch_size,
                     shuffle=True, collate_fn=collate_fn)
    import torch.nn.functional as F
    model.train()
    for ep in range(num_epochs):
        for b in dl:
            h, l, a = b['heavy_emb'].to(DEVICE), b['light_emb'].to(DEVICE), b['antigen_emb'].to(DEVICE)
            y = b['affinity'].to(DEVICE)
            tgt = 2 * (y - lo) / (hi - lo) - 1
            opt.zero_grad()
            pred = model(l, h, a)['cosine_similarity']
            loss = F.mse_loss(pred, tgt)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
    return model, (lo, hi)


def evaluate(model, df_test, bounds):
    lo, hi = bounds
    model.eval()
    dl = DataLoader(CachedEmbeddingDataset(df_test, loader, return_ids=False),
                     batch_size=64, collate_fn=collate_fn)
    p, t = [], []
    with torch.no_grad():
        for b in dl:
            cos = model(b['light_emb'].to(DEVICE), b['heavy_emb'].to(DEVICE),
                        b['antigen_emb'].to(DEVICE))['cosine_similarity'].cpu().numpy()
            p.extend(((cos + 1) / 2 * (hi - lo) + lo).tolist())
            t.extend(b['affinity'].tolist())
    p, t = np.array(p), np.array(t)
    r, _ = pearsonr(t, p); rho, _ = spearmanr(t, p)
    rmse = float(np.sqrt(np.mean((t - p) ** 2)))
    return r, rho, rmse


NUM_EPOCHS = 30  # matches MVSF-AB's --num-epochs default exactly
SEEDS = [9999, 114, 144]
rows = []
for seed in SEEDS:
    print(f"\n=== fixed {NUM_EPOCHS}-epoch, no early stopping, seed={seed} ===")
    model, bounds = train_fixed_epochs_no_earlystop(df_tr_full, NUM_EPOCHS, seed)
    r, rho, rmse = evaluate(model, df_b, bounds)
    print(f"  epoch-{NUM_EPOCHS} result: r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}")
    rows.append({'seed': seed, 'pearson': r, 'spearman': rho, 'rmse': rmse})
    del model; gc.collect()
    if DEVICE == 'cuda':
        torch.cuda.empty_cache()

df_res = pd.DataFrame(rows)
df_res.to_csv(os.path.join(OUT, 'benchmark_fixed_epoch_raw.csv'), index=False)
print("\n" + "=" * 70)
print(f"MVSF-AB-PROTOCOL-MATCHED (fixed {NUM_EPOCHS} epochs, no early stopping): "
      f"SAbDab -> Benchmark, {len(SEEDS)} seeds")
print("=" * 70)
print(f"Pearson r = {df_res.pearson.mean():.4f} +- {df_res.pearson.std():.4f}")
print(f"Spearman rho = {df_res.spearman.mean():.4f} +- {df_res.spearman.std():.4f}")
print(f"RMSE = {df_res.rmse.mean():.4f} +- {df_res.rmse.std():.4f}")
pd.DataFrame([{
    'pearson_mean': df_res.pearson.mean(), 'pearson_std': df_res.pearson.std(),
    'spearman_mean': df_res.spearman.mean(), 'spearman_std': df_res.spearman.std(),
    'rmse_mean': df_res.rmse.mean(), 'rmse_std': df_res.rmse.std(),
}]).to_csv(os.path.join(OUT, 'benchmark_fixed_epoch_summary.csv'), index=False)
print(f"\nExpected (from the revision): Pearson 0.564+-0.026, Spearman 0.570+-0.033, RMSE 1.364+-0.116")
