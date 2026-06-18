"""
saaintdb_allcdr.py
==================
Re-runs SAaIntDB 10-fold CV (random + cold splits) using CDR-aware
ESM-2 embeddings: pool only over CDR-H1, CDR-H2, CDR-H3 token positions
extracted via ANARCI IMGT numbering.

Compares directly to mean-pool baseline (results_mutual_strong_saaintdb/).
"""

import os, sys, glob, pickle, json, tempfile, shutil
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import KFold
from scipy.stats import pearsonr, spearmanr
from transformers import EsmTokenizer, EsmModel
import warnings; warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from AbAffinity.models.mutual_strong import MutualTriStreamStrong
from AbAffinity.models.mutual_strong_saaintdb import (
    load_saaintdb, _make_heavy_id, _make_light_id, _make_antigen_id,
    _train_fold, get_fold_splits, filter_missing_embeddings
)
from AbAffinity.utils.main_symmetric_mean import CachedEmbeddingLoader, CachedEmbeddingDataset, collate_fn, DEFAULT_CONFIG

OUT = os.path.join(HERE, 'results_saaintdb_allcdr')
os.makedirs(OUT, exist_ok=True)
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f"Device: {DEVICE}  |  Output: {OUT}")

# ─────────────────────────────────────────────────────────────────────────────
# CDR detection via ANARCI (IMGT numbering)
# ─────────────────────────────────────────────────────────────────────────────
# IMGT CDR definitions for VH:
# CDR-H1: positions 27-38
# CDR-H2: positions 56-65
# CDR-H3: positions 105-117 (variable end)
IMGT_CDR_H1 = set(range(27, 39))
IMGT_CDR_H2 = set(range(56, 66))
IMGT_CDR_H3 = set(range(105, 138))  # generous upper bound

def get_cdr_seq_indices(seq, scheme='imgt'):
    """
    Returns 0-indexed sequence positions for CDR-H1, CDR-H2, CDR-H3
    using ANARCI IMGT numbering.
    Falls back to regex pattern if ANARCI fails.
    """
    try:
        import anarci
        results = anarci.anarci([('q', seq)], scheme=scheme, assign_germline=False,
                                  allow=set('ACDEFGHIKLMNPQRSTVWY'))
        if (results and results[0] and results[0][0] and
                results[0][0][0] is not None):
            numbered = results[0][0][0]  # list of ((pos, ins_code), aa)
            h1_idx, h2_idx, h3_idx = [], [], []
            seq_pos = 0
            for (pos, _ins), aa in numbered:
                if aa == '-':
                    continue
                if pos in IMGT_CDR_H1:
                    h1_idx.append(seq_pos)
                elif pos in IMGT_CDR_H2:
                    h2_idx.append(seq_pos)
                elif pos in IMGT_CDR_H3:
                    h3_idx.append(seq_pos)
                seq_pos += 1
            all_cdr = sorted(set(h1_idx + h2_idx + h3_idx))
            if all_cdr:
                return all_cdr, h1_idx, h2_idx, h3_idx
    except Exception:
        pass
    return _regex_cdr_fallback(seq)

def _regex_cdr_fallback(seq):
    """Regex-based CDR detection for sequences ANARCI cannot number."""
    import re
    all_cdr = []
    # CDR-H1
    m1 = re.search(r'C[A-Z]{2,5}[SAGTV]([A-Z]{8,15})WVRQ', seq)
    if m1:
        s, e = m1.start(1), m1.end(1)
        h1 = list(range(s, e))
        all_cdr.extend(h1)
    else:
        h1 = list(range(26, 36))
        all_cdr.extend(h1)
    # CDR-H2
    m2 = re.search(r'W[VI]RQ[A-Z]{6,14}W[VL][AS]([A-Z]{10,20})VKGRF', seq)
    if m2:
        s, e = m2.start(1), m2.end(1)
        h2 = list(range(s, e))
        all_cdr.extend(h2)
    else:
        h2 = list(range(50, 64))
        all_cdr.extend(h2)
    # CDR-H3
    m3 = re.search(r'WYYCA([A-Z]+)WGQGT', seq)
    if not m3:
        m3 = re.search(r'WYYC[A-Z]([A-Z]+)(?:WGQGT|WGQG)', seq)
    if m3:
        s, e = m3.start(1), m3.end(1)
        h3 = list(range(s, e))
    else:
        h3 = list(range(95, 110))
    all_cdr.extend(h3)
    return sorted(set(all_cdr)), h1, h2, h3

# ─────────────────────────────────────────────────────────────────────────────
# Compute CDR-aware embeddings for SAaIntDB
# ─────────────────────────────────────────────────────────────────────────────
def compute_cdr_embeddings_saaintdb(id_to_seq, device, batch_size=8,
                                     model_name='facebook/esm2_t33_650M_UR50D',
                                     cache_path=None):
    if cache_path and os.path.isfile(cache_path):
        print(f"  Loading cached CDR embeddings: {cache_path}")
        with open(cache_path, 'rb') as f:
            return pickle.load(f)

    print(f"  Computing CDR-aware ESM-2 embeddings for {len(id_to_seq)} sequences...")
    tokenizer = EsmTokenizer.from_pretrained(model_name)
    esm = EsmModel.from_pretrained(model_name).to(device)
    esm.eval()

    unique_seqs = list(dict.fromkeys(id_to_seq.values()))
    seq_to_emb  = {}
    n_fallback  = 0

    for i in range(0, len(unique_seqs), batch_size):
        batch = unique_seqs[i:i+batch_size]
        inputs = tokenizer(batch, return_tensors='pt', padding=True,
                           truncation=True, max_length=512).to(device)
        with torch.no_grad():
            tok_embs = esm(**inputs).last_hidden_state  # [B, L, 1280]
        attn_mask = inputs['attention_mask']

        for j, seq in enumerate(batch):
            seq_len   = int(attn_mask[j].sum().item()) - 2
            cdr_idxs, *_ = get_cdr_seq_indices(seq)
            # Convert seq positions → token positions (+1 for CLS), clamp
            tok_idxs = [p + 1 for p in cdr_idxs if 0 <= p < seq_len]
            if tok_idxs:
                seq_to_emb[seq] = tok_embs[j, tok_idxs, :].mean(0).cpu().numpy()
            else:
                # Fallback: full mean-pool
                seq_to_emb[seq] = tok_embs[j, 1:seq_len+1, :].mean(0).cpu().numpy()
                n_fallback += 1

        if (i // batch_size) % 20 == 0:
            print(f"    {min(i+batch_size, len(unique_seqs))}/{len(unique_seqs)}", end='\r')

    del esm; torch.cuda.empty_cache()
    print(f"\n  Done.  Fallback (full mean-pool): {n_fallback}/{len(unique_seqs)} sequences")

    emb_dict = {sid: seq_to_emb[seq] for sid, seq in id_to_seq.items()}
    if cache_path:
        with open(cache_path, 'wb') as f:
            pickle.dump(emb_dict, f)
        print(f"  Saved: {cache_path}")
    return emb_dict

# ─────────────────────────────────────────────────────────────────────────────
# Load SAaIntDB
# ─────────────────────────────────────────────────────────────────────────────
print("\n[1] Loading SAaIntDB...")
csv_path = os.path.join(HERE, 'data/saaintdb_with_antigen_names.csv')
df_full  = load_saaintdb(csv_path)

# Build id→seq mapping for heavy chains only (CDR-aware)
# Light chains and antigens keep mean-pool (CDR-L1/L2/L3 are in light chain,
# but the key signal for binding comes from heavy chain CDRs)
print("\n[2] Computing CDR-aware embeddings for heavy chains...")
heavy_id_to_seq = {row['heavy_id']: row['H_seq'] for _, row in df_full.iterrows()
                   if isinstance(row.get('H_seq'), str)}

cdr_emb_cache = os.path.join(OUT, 'saaintdb_heavy_cdr_embeddings.pkl')
heavy_cdr_emb = compute_cdr_embeddings_saaintdb(
    heavy_id_to_seq, DEVICE, batch_size=8, cache_path=cdr_emb_cache)

# Load existing mean-pool embeddings for light and antigen
print("\n[3] Loading mean-pool embeddings for light/antigen...")
with open(os.path.join(HERE, 'data/esm2_embeddings_saaintdb_650M.pkl'), 'rb') as f:
    mean_emb = pickle.load(f)

# Build combined embedding dict: CDR-aware heavy + mean-pool light/antigen
combined_emb = {}
combined_emb.update(heavy_cdr_emb)
for sid, emb in mean_emb.items():
    if '_L_' in sid or '_Ag_' in sid:
        combined_emb[sid] = emb
    elif sid not in combined_emb:
        combined_emb[sid] = emb  # fallback for any missing

print(f"  Combined embedding dict: {len(combined_emb)} entries")
emb_dim = next(iter(combined_emb.values())).shape[0]
print(f"  Embedding dim: {emb_dim}")

# ─────────────────────────────────────────────────────────────────────────────
# Cross-validation helper
# ─────────────────────────────────────────────────────────────────────────────
from AbAffinity.utils.main_symmetric_mean import CachedEmbeddingLoader, CachedEmbeddingDataset, collate_fn

class DictEmbeddingLoader:
    """Wraps a plain dict to work as a CachedEmbeddingLoader."""
    def __init__(self, emb_dict):
        self.embeddings    = emb_dict
        self.embedding_dim = next(iter(emb_dict.values())).shape[0]

    def get_embedding(self, key):
        return self.embeddings[key]

def run_cv(df_full, emb_loader, config, n_folds, split_type, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    splits   = get_fold_splits(df_full, n_folds, config['seed'], split_type)
    id_cols  = ['heavy_id', 'light_id', 'antigen_id']
    fold_res = []
    all_preds_dfs = []

    print(f"\n{'='*65}")
    print(f"SAaIntDB All-CDR — {n_folds}-fold {split_type} CV")
    print(f"{'='*65}")

    for fold, (tr_idx, va_idx) in enumerate(splits, 1):
        df_tr = df_full.iloc[tr_idx].reset_index(drop=True)
        df_va = df_full.iloc[va_idx].reset_index(drop=True)

        # Filter missing embeddings
        keep_tr = df_tr[id_cols].apply(lambda col: col.isin(emb_loader.embeddings)).all(axis=1)
        keep_va = df_va[id_cols].apply(lambda col: col.isin(emb_loader.embeddings)).all(axis=1)
        df_tr = df_tr[keep_tr].reset_index(drop=True)
        df_va = df_va[keep_va].reset_index(drop=True)

        print(f"\n  FOLD {fold}/{n_folds}  train={len(df_tr)} val={len(df_va)}")
        model, pkd_bounds = _train_fold(df_tr, df_va, emb_loader, config, DEVICE)

        # Evaluate
        val_loader = DataLoader(CachedEmbeddingDataset(df_va, emb_loader, return_ids=False),
                                batch_size=config['batch_size'], collate_fn=collate_fn)
        lo, hi = pkd_bounds
        preds, trues = [], []
        model.eval()
        with torch.no_grad():
            for batch in val_loader:
                le = batch['light_emb'].to(DEVICE)
                he = batch['heavy_emb'].to(DEVICE)
                ae = batch['antigen_emb'].to(DEVICE)
                af = batch['affinity']
                cos = model(le, he, ae)['cosine_similarity'].cpu().numpy()
                pkd = (cos + 1.0) / 2.0 * (hi - lo) + lo
                preds.extend(pkd.tolist()); trues.extend(af.tolist())

        preds = np.array(preds); trues = np.array(trues)
        r,   _ = pearsonr(trues, preds)
        rho, _ = spearmanr(trues, preds)
        rmse   = float(np.sqrt(np.mean((trues - preds)**2)))
        print(f"  Fold {fold}: r={r:.4f}  rho={rho:.4f}  RMSE={rmse:.4f}")
        fold_res.append({'fold': fold, 'pearson': r, 'spearman': rho, 'rmse': rmse})

        # Save fold checkpoint
        fold_dir = os.path.join(output_dir, f'fold_{fold:02d}')
        os.makedirs(fold_dir, exist_ok=True)
        torch.save({'model_state_dict': {k: v.cpu() for k,v in model.state_dict().items()},
                    'pkd_bounds': pkd_bounds, 'config': config,
                    'val_metrics': {'pearson': r, 'spearman': rho, 'rmse': rmse}},
                   os.path.join(fold_dir, 'model.pt'))

        df_pred = df_va.copy()
        df_pred['predicted_affinity'] = preds
        df_pred['error'] = trues - preds
        df_pred['fold'] = fold
        all_preds_dfs.append(df_pred)

    fold_df = pd.DataFrame(fold_res)
    fold_df.to_csv(os.path.join(output_dir, 'cv_summary.csv'), index=False)
    pd.concat(all_preds_dfs).to_csv(os.path.join(output_dir, 'all_preds.csv'), index=False)

    mean_r   = fold_df['pearson'].mean()
    std_r    = fold_df['pearson'].std()
    mean_rho = fold_df['spearman'].mean()
    mean_rmse= fold_df['rmse'].mean()

    print(f"\n  RESULT ({split_type}): Pearson r={mean_r:.4f}±{std_r:.4f}  "
          f"Spearman={mean_rho:.4f}  RMSE={mean_rmse:.4f}")

    # Compare to mean-pool baseline
    baseline = pd.read_csv(os.path.join(HERE, 'results_mutual_strong_saaintdb/cv_summary.csv'))
    def _base(m):
        row = baseline[baseline.metric == m]
        return float(row['mean'].values[0]) if len(row) else np.nan

    print(f"\n  Baseline (mean-pool ESM-2):  r={_base('pearson'):.4f}  "
          f"rho={_base('spearman'):.4f}  RMSE={_base('rmse'):.4f}")
    print(f"  All-CDR gain:  Δr={mean_r - _base('pearson'):+.4f}  "
          f"Δrho={mean_rho - _base('spearman'):+.4f}")
    return fold_df

# ─────────────────────────────────────────────────────────────────────────────
# Run CV
# ─────────────────────────────────────────────────────────────────────────────
emb_loader = DictEmbeddingLoader(combined_emb)
config = {
    **DEFAULT_CONFIG,
    'device': DEVICE,
    'epochs': 50, 'patience': 10,
    'batch_size': 32, 'learning_rate': 1e-4, 'weight_decay': 0.01,
    'projected_size': 256, 'num_heads': 8, 'dropout': 0.1,
    'n_layers': 2, 'seed': 9999,
}

# Random split
rand_res = run_cv(df_full, emb_loader, config, n_folds=10,
                  split_type='random', output_dir=os.path.join(OUT, 'random'))

# Cold split
cold_res = run_cv(df_full, emb_loader, config, n_folds=10,
                  split_type='cold', output_dir=os.path.join(OUT, 'cold'))

# ─────────────────────────────────────────────────────────────────────────────
# Final comparison table
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'='*65}")
print("FINAL COMPARISON: Mean-pool vs CDR-aware (All-CDR H1+H2+H3)")
print(f"{'='*65}")

baseline  = pd.read_csv(os.path.join(HERE, 'results_mutual_strong_saaintdb/cv_summary.csv'))
def _b(m): return float(baseline[baseline.metric==m]['mean'].values[0]) if len(baseline[baseline.metric==m]) else np.nan

print(f"\n  {'Method':<30} {'Random r':>10} {'Random rho':>12} {'Cold r':>10}")
print(f"  {'-'*65}")

cold_baseline = pd.read_csv(os.path.join(HERE,
    'figures_out_saaintdb/additional_experiments/fold_level_errorbars_all_models.csv'))
cold_r = float(cold_baseline[(cold_baseline.model=='BALM cold')&(cold_baseline.metric=='pearson')]['mean'].values[0]) \
         if len(cold_baseline[(cold_baseline.model=='BALM cold')&(cold_baseline.metric=='pearson')]) else np.nan

print(f"  {'ESM-2 mean-pool (baseline)':<30} {_b('pearson'):>10.4f} {_b('spearman'):>12.4f} {cold_r:>10.4f}")
print(f"  {'All-CDR H1+H2+H3 (ours)':<30} {rand_res.pearson.mean():>10.4f} {rand_res.spearman.mean():>12.4f} {cold_res.pearson.mean():>10.4f}")
delta_rand = rand_res.pearson.mean() - _b('pearson')
delta_cold = cold_res.pearson.mean() - cold_r
print(f"  {'Δ (CDR − baseline)':<30} {delta_rand:>+10.4f} {'':>12} {delta_cold:>+10.4f}")

comparison = pd.DataFrame([
    {'method': 'ESM-2 mean-pool', 'split': 'random', 'pearson': _b('pearson'),
     'spearman': _b('spearman'), 'rmse': _b('rmse')},
    {'method': 'All-CDR H1+H2+H3', 'split': 'random', 'pearson': rand_res.pearson.mean(),
     'spearman': rand_res.spearman.mean(), 'rmse': rand_res.rmse.mean()},
    {'method': 'ESM-2 mean-pool', 'split': 'cold', 'pearson': cold_r,
     'spearman': np.nan, 'rmse': np.nan},
    {'method': 'All-CDR H1+H2+H3', 'split': 'cold', 'pearson': cold_res.pearson.mean(),
     'spearman': cold_res.spearman.mean(), 'rmse': cold_res.rmse.mean()},
])
comparison.to_csv(os.path.join(OUT, 'saaintdb_allcdr_comparison.csv'), index=False)
print(f"\nSaved: {OUT}/saaintdb_allcdr_comparison.csv")
