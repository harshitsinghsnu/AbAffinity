"""
post_hoc_corrected.py
=====================
CORRECTED post-hoc gating + stream-intervention ablation, producing
Supplementary Table 15 (gate ablation) and Table 16 (stream intervention).

A previous version of this analysis re-created the outer folds with
get_fold_splits(df, 10, 9999, 'random'), which does not reproduce the exact
row-level fold assignment the checkpoints were trained on (KFold's shuffle
is sensitive to df's exact row order and any upstream filtering) -- the
baseline was evaluated partly in-sample, inflating results (observed
inflation from a true ~0.84 to ~0.91-0.94). Here, each fold checkpoint is
instead evaluated on its TRUE held-out rows, defined by the `fold` column
of results_saaintdb_allcdr/random/all_preds.csv (the exact out-of-fold
assignment recorded during training). The 'learned'/'full' baseline
therefore reproduces the corrected 10-fold CV headline (~0.838).

Outputs -> results_allcdr_stats/
  allcdr_gating_saaintdb.csv           (Table 15)
  allcdr_stream_interventions.csv      (Table 16)
  allcdr_stream_interventions_folds.csv
"""
import os, sys, glob, pickle
import numpy as np, pandas as pd, torch
from scipy.stats import pearsonr, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from AbAffinity.models.mutual_strong import MutualTriStreamStrong, GatedCrossAttention

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
OUT = os.path.join(HERE, 'results_allcdr_stats')
os.makedirs(OUT, exist_ok=True)

with open(os.path.join(HERE, 'data/esm2_embeddings_saaintdb_650M.pkl'), 'rb') as f:
    emb = dict(pickle.load(f))
with open(os.path.join(HERE, 'results_saaintdb_allcdr/saaintdb_heavy_cdr_embeddings.pkl'), 'rb') as f:
    emb.update(pickle.load(f))

ap = pd.read_csv(os.path.join(HERE, 'results_saaintdb_allcdr/random/all_preds.csv'))
idc = ['heavy_id', 'light_id', 'antigen_id']
ap = ap[ap[idc].apply(lambda c: c.isin(emb)).all(axis=1)].reset_index(drop=True)
ag_mean = np.stack([emb[a] for a in ap['antigen_id'].unique()]).mean(0).astype(np.float32)


def set_gate(model, mode):
    for m in model.modules():
        if isinstance(m, GatedCrossAttention):
            if mode is None:
                if hasattr(m, '_orig'):
                    m.forward = m._orig
                continue
            if not hasattr(m, '_orig'):
                m._orig = m.forward

            def mk(mod, mode):
                def fwd(q, kv):
                    res = q; Q = mod.W_q(q); K = mod.W_k(kv); V = mod.W_v(kv)
                    sd = torch.tanh(Q * K) * V
                    if mode == 'open':
                        g = torch.ones_like(sd)
                    elif mode == 'closed':
                        g = torch.zeros_like(sd)
                    elif mode == 'fixed':
                        g = 0.5 * torch.ones_like(sd)
                    elif mode == 'random':
                        g = torch.rand_like(sd)
                    else:
                        g = torch.sigmoid(mod.W_gate(q))
                    return mod.layer_norm(res + mod.dropout(mod.W_o(sd * g)))
                return fwd
            m.forward = mk(m, mode)


def load_fold(fp):
    ck = torch.load(fp, map_location=DEVICE, weights_only=False); cfg = ck.get('config', {})
    m = MutualTriStreamStrong(esm_dim=1280, projected_size=cfg.get('projected_size', 256),
        num_heads=cfg.get('num_heads', 8), dropout=cfg.get('dropout', 0.1),
        n_layers=cfg.get('n_layers', 2), device=DEVICE).to(DEVICE)
    m.load_state_dict(ck['model_state_dict']); m.eval(); return m, ck['pkd_bounds']


@torch.no_grad()
def evalc(m, dv, bounds, cond, seed):
    lo, hi = bounds
    H = np.stack([emb[i] for i in dv['heavy_id']]).astype(np.float32)
    L = np.stack([emb[i] for i in dv['light_id']]).astype(np.float32)
    A = np.stack([emb[i] for i in dv['antigen_id']]).astype(np.float32)
    if cond == 'zero_light':
        L = np.zeros_like(L)
    elif cond == 'zero_heavy':
        H = np.zeros_like(H)
    elif cond == 'zero_antigen':
        A = np.zeros_like(A)
    elif cond == 'mean_antigen':
        A = np.tile(ag_mean, (len(A), 1))
    elif cond == 'shuffled_antigen':
        A = A[np.random.default_rng(seed).permutation(len(A))]
    P = []
    for i in range(0, len(dv), 256):
        cos = m(torch.tensor(L[i:i+256]).to(DEVICE), torch.tensor(H[i:i+256]).to(DEVICE),
                torch.tensor(A[i:i+256]).to(DEVICE))['cosine_similarity'].cpu().numpy()
        P.extend(((cos + 1) / 2 * (hi - lo) + lo).tolist())
    t = dv['binding_affinity'].values; p = np.array(P)
    return float(pearsonr(t, p)[0]), float(spearmanr(t, p)[0]), float(np.sqrt(np.mean((t - p) ** 2)))


folds = sorted(glob.glob(os.path.join(HERE, 'results_saaintdb_allcdr/random/fold_*/model.pt')))
gate_modes = ['learned', 'fixed', 'open', 'random', 'closed']
conds = ['full', 'zero_light', 'zero_heavy', 'zero_antigen', 'mean_antigen', 'shuffled_antigen']
grows, srows = [], []
for fi, fp in enumerate(folds, 1):
    dv = ap[ap['fold'] == fi].reset_index(drop=True)
    if len(dv) < 2:
        continue
    m, b = load_fold(fp)
    for gm in gate_modes:
        set_gate(m, None if gm == 'learned' else gm)
        r, rho, rmse = evalc(m, dv, b, 'full', 9999 + fi)
        grows.append({'fold': fi, 'gate_mode': gm, 'pearson': r, 'spearman': rho, 'rmse': rmse})
    set_gate(m, None)
    for c in conds:
        r, rho, rmse = evalc(m, dv, b, c, 9999 + fi)
        srows.append({'fold': fi, 'condition': c, 'pearson': r, 'spearman': rho, 'rmse': rmse})
    base = [x for x in grows if x['fold'] == fi and x['gate_mode'] == 'learned'][0]['pearson']
    clo = [x for x in grows if x['fold'] == fi and x['gate_mode'] == 'closed'][0]['pearson']
    print(f"fold {fi}: learned={base:.3f} closed={clo:.3f}", flush=True)

gd, sd = pd.DataFrame(grows), pd.DataFrame(srows)
agg = lambda d, k, order: (d.groupby(k).agg(pearson_mean=('pearson', 'mean'), pearson_std=('pearson', 'std'),
        spearman_mean=('spearman', 'mean'), spearman_std=('spearman', 'std'),
        rmse_mean=('rmse', 'mean'), rmse_std=('rmse', 'std')).reindex(order).reset_index())
gagg = agg(gd, 'gate_mode', gate_modes); gagg['dataset'] = 'saaintdb'
sagg = agg(sd, 'condition', conds)
gagg.to_csv(os.path.join(OUT, 'allcdr_gating_saaintdb.csv'), index=False)
sagg.to_csv(os.path.join(OUT, 'allcdr_stream_interventions.csv'), index=False)
sd.to_csv(os.path.join(OUT, 'allcdr_stream_interventions_folds.csv'), index=False)
print("\n=== CORRECTED GATING (true held-out) -- Table 15 ===")
print(gagg[['gate_mode', 'pearson_mean', 'pearson_std', 'spearman_mean', 'rmse_mean']].to_string(index=False))
print("\n=== CORRECTED STREAM (true held-out) -- Table 16 ===")
print(sagg[['condition', 'pearson_mean', 'pearson_std', 'spearman_mean', 'rmse_mean']].to_string(index=False))
print("saved corrected CSVs")
print("\nExpected (from the revision): Learned gate r=0.838+-0.033; Shuffled-antigen r=0.365+-0.056")
