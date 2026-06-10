#!/usr/bin/env python3
"""
Two‑Stream Symmetric – Ablation of chain separation.

Runs:
  1. 10‑fold CV on SAbDab, AbBind, SKEMPI
  2. Benchmark: train on 100% SAbDab → test on Benchmark (3 seeds)

All configuration matches the three‑stream symmetric_mean.
"""

import os
import random
import pickle
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
from copy import deepcopy

# ==============================================================================
# UTILITIES (identical to three‑stream)
# ==============================================================================

def setup_reproducibility(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
    os.environ['PYTHONHASHSEED'] = str(seed)

def pearson_correlation(y_true, y_pred):
    y_true, y_pred = np.array(y_true, dtype=float), np.array(y_pred, dtype=float)
    y_true_c = y_true - y_true.mean()
    y_pred_c = y_pred - y_pred.mean()
    num = np.sum(y_true_c * y_pred_c)
    den = np.sqrt(np.sum(y_true_c**2)) * np.sqrt(np.sum(y_pred_c**2))
    return float(num / den) if den > 1e-12 else 0.0

def spearman_correlation(y_true, y_pred):
    def ranks(data):
        order = np.argsort(data)
        r = np.empty_like(order, dtype=float)
        r[order] = np.arange(len(data))
        for val in np.unique(data):
            mask = data == val
            r[mask] = r[mask].mean()
        return r + 1
    return pearson_correlation(ranks(np.array(y_true)), ranks(np.array(y_pred)))

def rmse(y_true, y_pred):
    return float(np.sqrt(mean_squared_error(y_true, y_pred)))

def r_squared(y_true, y_pred):
    y_true, y_pred = np.array(y_true, dtype=float), np.array(y_pred, dtype=float)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - y_true.mean()) ** 2)
    return float(1.0 - ss_res / ss_tot) if ss_tot > 1e-12 else 0.0

def affinity_to_cosine(pkd, pkd_lower, pkd_upper):
    return 2.0 * (pkd - pkd_lower) / (pkd_upper - pkd_lower) - 1.0

def cosine_to_affinity(cosine_sim, pkd_lower, pkd_upper):
    return (cosine_sim + 1.0) / 2.0 * (pkd_upper - pkd_lower) + pkd_lower

# ==============================================================================
# DATA LOADING (copied from three‑stream)
# ==============================================================================

def load_data(pairs_csv):
    print(f"  Loading {pairs_csv}")
    df = pd.read_csv(pairs_csv)

    light_col = heavy_col = antigen_col = affinity_col = None
    for col in df.columns:
        c = col.lower().strip()
        if 'light' in c and light_col is None:
            light_col = col
        elif 'heavy' in c and heavy_col is None:
            heavy_col = col
        elif 'antigen' in c and antigen_col is None:
            antigen_col = col
        elif c in ('y', 'affinity', 'binding_affinity', 'pkd', 'ddg', 'dg', 'kd', 'ic50') and affinity_col is None:
            affinity_col = col

    missing = []
    if light_col is None: missing.append("light chain ID (col with 'light')")
    if heavy_col is None: missing.append("heavy chain ID (col with 'heavy')")
    if antigen_col is None: missing.append("antigen ID (col with 'antigen')")
    if affinity_col is None: missing.append("affinity column")
    if missing:
        raise ValueError(f"Missing columns in {pairs_csv}:\n  {missing}\nAvailable: {list(df.columns)}")

    df = df.rename(columns={
        light_col: 'light_id',
        heavy_col: 'heavy_id',
        antigen_col: 'antigen_id',
        affinity_col: 'binding_affinity'
    })
    df = df.dropna(subset=['light_id', 'heavy_id', 'antigen_id', 'binding_affinity'])
    print(f"  Loaded {len(df)} pairs | affinity [{df['binding_affinity'].min():.3f}, {df['binding_affinity'].max():.3f}]")
    return df

# ==============================================================================
# GATED CROSS-ATTENTION (identical to three‑stream)
# ==============================================================================

class GatedCrossAttention(nn.Module):
    def __init__(self, query_dim, kv_dim, num_heads=8, dropout=0.1):
        super().__init__()
        assert query_dim % num_heads == 0
        self.num_heads = num_heads
        self.query_dim = query_dim

        self.W_q = nn.Linear(query_dim, query_dim, bias=False)
        self.W_k = nn.Linear(kv_dim, query_dim, bias=False)
        self.W_v = nn.Linear(kv_dim, query_dim, bias=False)
        self.W_o = nn.Linear(query_dim, query_dim, bias=False)
        self.W_gate = nn.Linear(query_dim, query_dim, bias=True)

        self.dropout = nn.Dropout(dropout)
        self.layer_norm = nn.LayerNorm(query_dim)
        self._init_weights()

    def _init_weights(self):
        for m in [self.W_q, self.W_k, self.W_v, self.W_o]:
            nn.init.xavier_uniform_(m.weight)
        nn.init.normal_(self.W_gate.weight, mean=0.0, std=0.01)
        nn.init.zeros_(self.W_gate.bias)

    def forward(self, query, key_value):
        residual = query
        q_proj = self.W_q(query)
        k_proj = self.W_k(key_value)
        v_proj = self.W_v(key_value)
        interaction = torch.tanh(q_proj * k_proj)
        sdpa_out = interaction * v_proj
        gate = torch.sigmoid(self.W_gate(query))
        gated_out = sdpa_out * gate
        output = self.W_o(gated_out)
        output = self.dropout(output)
        output = self.layer_norm(residual + output)
        return output

# ==============================================================================
# TWO-STREAM EMBEDDING LOADER & DATASET (unchanged)
# ==============================================================================
class TwoStreamEmbeddingLoader:
    def __init__(self, cache_file):
        print(f"  Loading two-stream embeddings from {cache_file}")
        with open(cache_file, 'rb') as f:
            package = pickle.load(f)
        self.antibody_embeddings = package['antibody']
        self.antigen_embeddings = package['antigen']
        self.embedding_dim = next(iter(self.antibody_embeddings.values())).shape[0]
        print(f"  Loaded {len(self.antibody_embeddings)} antibody, {len(self.antigen_embeddings)} antigen embeddings")

    def get_antibody_embedding(self, antibody_id):
        return self.antibody_embeddings[antibody_id]

    def get_antigen_embedding(self, antigen_id):
        return self.antigen_embeddings[antigen_id]

    def validate_ids(self, df, antibody_col='antibody_id', antigen_col='antigen_id'):
        missing = []
        for aid in df[antibody_col].unique():
            if aid not in self.antibody_embeddings:
                missing.append(('antibody', aid))
        for agid in df[antigen_col].unique():
            if agid not in self.antigen_embeddings:
                missing.append(('antigen', agid))
        if missing:
            raise ValueError(f"Missing embeddings: {missing[:5]}")
        print(f"  Cache validation OK.")

class TwoStreamDataset(Dataset):
    def __init__(self, df, embedding_loader, return_ids=False):
        self.df = df.reset_index(drop=True)
        self.loader = embedding_loader
        self.return_ids = return_ids

    def __len__(self): return len(self.df)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        antibody_id = f"{row['heavy_id']}_{row['light_id']}"
        ab_emb = self.loader.get_antibody_embedding(antibody_id)
        ag_emb = self.loader.get_antigen_embedding(row['antigen_id'])

        item = {
            'antibody_emb': torch.tensor(ab_emb, dtype=torch.float32),
            'antigen_emb': torch.tensor(ag_emb, dtype=torch.float32),
            'affinity': torch.tensor(row['binding_affinity'], dtype=torch.float32)
        }
        if self.return_ids:
            item['antibody_id'] = antibody_id
            item['antigen_id'] = row['antigen_id']
        return item

def collate_fn(batch):
    out = {
        'antibody_emb': torch.stack([b['antibody_emb'] for b in batch]),
        'antigen_emb': torch.stack([b['antigen_emb'] for b in batch]),
        'affinity': torch.tensor([b['affinity'] for b in batch], dtype=torch.float32),
    }
    if 'antibody_id' in batch[0]:
        out['antibody_id'] = [b['antibody_id'] for b in batch]
        out['antigen_id'] = [b['antigen_id'] for b in batch]
    return out

# ==============================================================================
# TWO-STREAM MODEL (Main Change)
# ==============================================================================
class TwoStreamSymmetric(nn.Module):
    """
    Fair Two-Stream model for ablation:
    - Antibody = Heavy + CLS + Light (pre-concatenated)
    - Bidirectional gated cross-attention (same as three-stream)
    - Small learnable fusion gate after attention (makes it fairer)
    - Dual heads for cosine (same anti-collapse mechanism)
    """
    def __init__(self, esm_dim=1280, projected_size=256, num_heads=8, dropout=0.1, device='cuda'):
        super().__init__()
        self.projected_size = projected_size

        self.antibody_projection = nn.Sequential(
            nn.Linear(esm_dim, projected_size), nn.LayerNorm(projected_size),
            nn.GELU(), nn.Dropout(dropout)
        )
        self.antigen_projection = nn.Sequential(
            nn.Linear(esm_dim, projected_size), nn.LayerNorm(projected_size),
            nn.GELU(), nn.Dropout(dropout)
        )

        # Bidirectional gated cross-attention
        self.ab_to_ag_attn = GatedCrossAttention(projected_size, projected_size, num_heads, dropout)
        self.ag_to_ab_attn = GatedCrossAttention(projected_size, projected_size, num_heads, dropout)

        # Small learnable fusion gate (makes comparison fairer)
        self.fusion_gate = nn.Linear(projected_size, 1)   # learns how much to keep from original vs attended

        # Dual heads (same as three-stream)
        self.antibody_head = nn.Sequential(
            nn.Linear(projected_size, projected_size),
            nn.LayerNorm(projected_size), nn.GELU(), nn.Dropout(dropout)
        )
        self.antigen_head = nn.Sequential(
            nn.Linear(projected_size, projected_size),
            nn.LayerNorm(projected_size), nn.GELU(), nn.Dropout(dropout)
        )

        self._init_weights()

    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)

    def forward(self, antibody_emb, antigen_emb):
        antibody_emb = F.layer_norm(antibody_emb, antibody_emb.shape[-1:])
        antigen_emb  = F.layer_norm(antigen_emb,  antigen_emb.shape[-1:])

        ab_proj = self.antibody_projection(antibody_emb)
        ag_proj = self.antigen_projection(antigen_emb)

        # Bidirectional gated cross-attention
        ab_to_ag_ctx = self.ab_to_ag_attn(ab_proj, ag_proj)      # Ab → Ag
        ag_to_ab_ctx = self.ag_to_ab_attn(ag_proj, ab_proj)      # Ag → Ab

        # Use conditioned representations
        ab_ctx = ag_to_ab_ctx
        ag_ctx = ab_to_ag_ctx

        # === Small Learnable Fusion Gate (Fairness addition) ===
        gate = torch.sigmoid(self.fusion_gate(ab_ctx))           # [B, 1]
        ab_ctx = gate * ab_ctx + (1 - gate) * ab_proj            # residual-style fusion

        # Dual heads + Cosine
        antibody_ctx = self.antibody_head(ab_ctx)
        antigen_ctx  = self.antigen_head(ag_ctx)

        ab_norm = F.normalize(antibody_ctx, p=2, dim=-1)
        ag_norm = F.normalize(antigen_ctx,  p=2, dim=-1)
        cosine_sim = F.cosine_similarity(ab_norm, ag_norm, dim=-1)

        return {
            'cosine_similarity': cosine_sim,
            'antibody_context': ab_norm,
            'antigen_context': ag_norm
        }

# ==============================================================================
# CONCAT TWO-STREAM MODEL
# Replaces heavy+light self-attention with direct concatenation of CLS tokens.
# Architecture after concatenation mirrors MutualTriStreamStrong exactly.
# ==============================================================================

class ConcatTwoStream(nn.Module):
    """
    Two-stream antibody–antigen affinity model using pre-combined antibody embeddings.

    The antibody embedding is precomputed by running ESM-2 over the full
    "heavy + <cls><cls> + light" sequence and mean-pooling, so it arrives as
    a single esm_dim-dimensional vector — no internal concatenation needed.

    Bidirectional gated cross-attention + dual heads + cosine are identical
    to MutualTriStreamStrong (mutual_strong.py).
    """
    def __init__(self, esm_dim=1280, projected_size=256, num_heads=8,
                 dropout=0.1, n_layers=2, device='cuda'):
        super().__init__()
        self.projected_size = projected_size
        self.n_layers = n_layers

        self.antibody_projection = nn.Sequential(
            nn.Linear(esm_dim, projected_size), nn.LayerNorm(projected_size),
            nn.GELU(), nn.Dropout(dropout)
        )
        self.antigen_projection = nn.Sequential(
            nn.Linear(esm_dim, projected_size), nn.LayerNorm(projected_size),
            nn.GELU(), nn.Dropout(dropout)
        )

        self.ab_to_ag_layers = nn.ModuleList([
            GatedCrossAttention(projected_size, projected_size, num_heads, dropout)
            for _ in range(n_layers)
        ])
        self.ag_to_ab_layers = nn.ModuleList([
            GatedCrossAttention(projected_size, projected_size, num_heads, dropout)
            for _ in range(n_layers)
        ])

        self.antibody_head = nn.Sequential(
            nn.Linear(projected_size, projected_size),
            nn.LayerNorm(projected_size), nn.GELU(), nn.Dropout(dropout)
        )
        self.antigen_head = nn.Sequential(
            nn.Linear(projected_size, projected_size),
            nn.LayerNorm(projected_size), nn.GELU(), nn.Dropout(dropout)
        )

        self._init_weights()

    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)

    def forward(self, antibody_emb, antigen_emb):
        antibody_emb = F.layer_norm(antibody_emb, antibody_emb.shape[-1:])
        antigen_emb  = F.layer_norm(antigen_emb,  antigen_emb.shape[-1:])

        ab_proj = self.antibody_projection(antibody_emb)  # [B, projected_size]
        ag_proj = self.antigen_projection(antigen_emb)    # [B, projected_size]

        ag_ctx = ag_proj
        for i in range(self.n_layers):
            # Antibody queries antigen (ab_proj conditioned on ag_ctx)
            new_ag_ctx = self.ab_to_ag_layers[i](ab_proj, ag_ctx)
            # Antigen queries antibody — mirrors mutual_strong.py (output unused for ab update)
            self.ag_to_ab_layers[i](ag_ctx, ab_proj)
            ag_ctx = new_ag_ctx

        ab_out = self.antibody_head(ab_proj)
        ag_out = self.antigen_head(ag_ctx)

        ab_norm = F.normalize(ab_out, p=2, dim=-1)
        ag_norm = F.normalize(ag_out, p=2, dim=-1)
        cosine_sim = F.cosine_similarity(ab_norm, ag_norm, dim=-1)

        return {
            'cosine_similarity': cosine_sim,
            'antibody_context': ab_norm,
            'antigen_context': ag_norm,
        }


def train_epoch_concat(model, loader, optimizer, device, pkd_bounds):
    """train_epoch for ConcatTwoStream: batch keys antibody_emb, antigen_emb."""
    model.train()
    total_loss = 0.0
    pkd_lower, pkd_upper = pkd_bounds
    for batch in loader:
        ab  = batch['antibody_emb'].to(device)
        ae  = batch['antigen_emb'].to(device)
        aff = batch['affinity'].to(device)

        target = affinity_to_cosine(aff, pkd_lower, pkd_upper)
        out = model(ab, ae)
        loss = F.mse_loss(out['cosine_similarity'], target)

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        total_loss += loss.item()
    return total_loss / len(loader)


def evaluate_concat(model, loader, device, pkd_bounds):
    """evaluate for ConcatTwoStream: batch keys antibody_emb, antigen_emb."""
    model.eval()
    preds, labels = [], []
    total_loss = 0.0
    pkd_lower, pkd_upper = pkd_bounds
    with torch.no_grad():
        for batch in loader:
            ab  = batch['antibody_emb'].to(device)
            ae  = batch['antigen_emb'].to(device)
            aff = batch['affinity'].to(device)

            target = affinity_to_cosine(aff, pkd_lower, pkd_upper)
            out = model(ab, ae)
            loss = F.mse_loss(out['cosine_similarity'], target)
            total_loss += loss.item()

            pred_pkd = cosine_to_affinity(out['cosine_similarity'], pkd_lower, pkd_upper)
            preds.extend(pred_pkd.cpu().numpy().tolist())
            labels.extend(aff.cpu().numpy().tolist())

    metrics = {
        'loss': total_loss / len(loader),
        'r2': r_squared(labels, preds),
        'pearson': pearson_correlation(labels, preds),
        'spearman': spearman_correlation(labels, preds),
        'rmse': rmse(labels, preds),
    }
    return metrics, preds, labels


# ==============================================================================
# TRAINING & EVALUATION (identical to three‑stream)
# ==============================================================================

def train_epoch(model, loader, optimizer, device, pkd_bounds):
    model.train()
    total_loss = 0.0
    pkd_lower, pkd_upper = pkd_bounds
    for batch in loader:
        ab_emb = batch['antibody_emb'].to(device)
        ag_emb = batch['antigen_emb'].to(device)
        aff = batch['affinity'].to(device)

        target = affinity_to_cosine(aff, pkd_lower, pkd_upper)
        out = model(ab_emb, ag_emb)
        loss = F.mse_loss(out['cosine_similarity'], target)

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        total_loss += loss.item()
    return total_loss / len(loader)

def evaluate(model, loader, device, pkd_bounds, return_ids=False):
    model.eval()
    preds, labels = [], []
    ab_ids, ag_ids = [], []
    total_loss = 0.0
    pkd_lower, pkd_upper = pkd_bounds

    with torch.no_grad():
        for batch in loader:
            ab_emb = batch['antibody_emb'].to(device)
            ag_emb = batch['antigen_emb'].to(device)
            aff = batch['affinity'].to(device)

            target = affinity_to_cosine(aff, pkd_lower, pkd_upper)
            out = model(ab_emb, ag_emb)
            loss = F.mse_loss(out['cosine_similarity'], target)
            total_loss += loss.item()

            pred_pkd = cosine_to_affinity(out['cosine_similarity'], pkd_lower, pkd_upper)
            preds.extend(pred_pkd.cpu().numpy().tolist())
            labels.extend(aff.cpu().numpy().tolist())

            if return_ids and 'antibody_id' in batch:
                ab_ids.extend(batch['antibody_id'])
                ag_ids.extend(batch['antigen_id'])

    metrics = {
        'loss': total_loss / len(loader),
        'r2': r_squared(labels, preds),
        'pearson': pearson_correlation(labels, preds),
        'spearman': spearman_correlation(labels, preds),
        'rmse': rmse(labels, preds),
    }
    if return_ids:
        return metrics, preds, labels, ab_ids, ag_ids
    return metrics, preds, labels

def train_model_twostream(df_train, df_val, embedding_loader, config, device):
    pkd_lower = df_train['binding_affinity'].min()
    pkd_upper = df_train['binding_affinity'].max()
    pkd_bounds = (pkd_lower, pkd_upper)
    print(f"  Affinity range (train): [{pkd_lower:.3f}, {pkd_upper:.3f}]")

    train_loader = DataLoader(
        TwoStreamDataset(df_train, embedding_loader),
        batch_size=config['batch_size'], shuffle=True, collate_fn=collate_fn,
        num_workers=0, pin_memory=(config['device'] == 'cuda')
    )
    val_loader = DataLoader(
        TwoStreamDataset(df_val, embedding_loader),
        batch_size=config['batch_size'], collate_fn=collate_fn,
        num_workers=0, pin_memory=(config['device'] == 'cuda')
    )

    model = TwoStreamSymmetric(
        esm_dim=embedding_loader.embedding_dim,
        projected_size=config['projected_size'],
        num_heads=config['num_heads'],
        dropout=config['dropout'],
        device=config['device']
    ).to(config['device'])

    n = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"  Trainable parameters: {n:,}")

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config['learning_rate'], weight_decay=config['weight_decay']
    )

    best_pearson, best_state, patience = -float('inf'), None, 0
    for epoch in range(config['epochs']):
        train_loss = train_epoch(model, train_loader, optimizer, config['device'], pkd_bounds)
        val_metrics, _, _ = evaluate(model, val_loader, config['device'], pkd_bounds)

        print(f"  Epoch {epoch+1:3d}/{config['epochs']} | "
              f"Loss: {train_loss:.4f} | "
              f"Val Pearson: {val_metrics['pearson']:.4f} | "
              f"Val RMSE: {val_metrics['rmse']:.4f}")

        if val_metrics['pearson'] > best_pearson:
            best_pearson = val_metrics['pearson']
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            patience = 0
        else:
            patience += 1
            if patience >= config['patience']:
                print(f"  Early stopping (best val Pearson: {best_pearson:.4f})")
                break

    if best_state:
        model.load_state_dict({k: v.to(config['device']) for k, v in best_state.items()})
    return model, pkd_bounds

# ==============================================================================
# CROSS‑VALIDATION (10-fold)
# ==============================================================================

def cross_validate_twostream(df, embedding_loader, config, n_folds=10, output_dir='results_twostream', dataset_name=''):
    kfold = KFold(n_splits=n_folds, shuffle=True, random_state=config['seed'])
    fold_results = []
    all_preds = []
    device = config['device']
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{n_folds}-fold CV | dataset={dataset_name} | n={len(df)}")

    for fold, (tr_idx, va_idx) in enumerate(kfold.split(df), 1):
        print(f"\n{'='*70}\nFOLD {fold}/{n_folds}\n{'='*70}")
        df_tr = df.iloc[tr_idx].reset_index(drop=True)
        df_va = df.iloc[va_idx].reset_index(drop=True)
        print(f"  Train: {len(df_tr)} | Test: {len(df_va)}")

        model, pkd_bounds = train_model_twostream(df_tr, df_va, embedding_loader, config, device)

        val_loader = DataLoader(
            TwoStreamDataset(df_va, embedding_loader, return_ids=True),
            batch_size=config['batch_size'], collate_fn=collate_fn
        )
        metrics, preds, labels, ab_ids, ag_ids = evaluate(
            model, val_loader, device, pkd_bounds, return_ids=True
        )

        fold_dir = os.path.join(output_dir, f'fold_{fold:02d}')
        os.makedirs(fold_dir, exist_ok=True)

        fold_df = pd.DataFrame({
            'antibody_id': ab_ids, 'antigen_id': ag_ids,
            'true_affinity': labels, 'predicted_affinity': preds,
            'error': np.array(labels) - np.array(preds)
        })
        fold_df.to_csv(os.path.join(fold_dir, 'predictions.csv'), index=False)
        pd.DataFrame([metrics]).to_csv(os.path.join(fold_dir, 'metrics.csv'), index=False)
        pd.DataFrame([{'pkd_lower': pkd_bounds[0], 'pkd_upper': pkd_bounds[1]}]).to_csv(
            os.path.join(fold_dir, 'pkd_bounds.csv'), index=False)

        all_preds.append(fold_df)
        fold_results.append({'fold': fold, 'metrics': metrics, 'pkd_bounds': pkd_bounds})

        print(f"\n  Fold {fold}:")
        for k, v in metrics.items():
            print(f"    {k.upper():10s}: {v:.4f}")

    pd.concat(all_preds, ignore_index=True).to_csv(
        os.path.join(output_dir, 'all_folds_predictions.csv'), index=False)

    summary = {}
    for metric in ['r2', 'pearson', 'spearman', 'rmse', 'loss']:
        vals = [f['metrics'][metric] for f in fold_results]
        summary[metric] = {'mean': float(np.mean(vals)), 'std': float(np.std(vals))}

    print(f"\n{'='*70}\nCV SUMMARY — {dataset_name}\n{'='*70}")
    for metric, s in summary.items():
        print(f"  {metric.upper():10s}: {s['mean']:.4f} ± {s['std']:.4f}")

    pd.DataFrame([{'metric': k, 'mean': v['mean'], 'std': v['std']}
                  for k, v in summary.items()]).to_csv(
        os.path.join(output_dir, 'cv_summary.csv'), index=False)
    return fold_results, summary

# ==============================================================================
# BENCHMARK: SINGLE SEED (used internally)
# ==============================================================================

def train_on_sabdab_test_on_benchmark_twostream_single_seed(
    sabdab_csv, benchmark_csv,
    train_cache_file, test_cache_file,
    config, seed, save_model=True,
    model_path=None, output_dir=None
):
    """Run one seed of benchmark."""
    setup_reproducibility(seed)
    config['seed'] = seed
    device = config['device']

    df_train = load_data(sabdab_csv)
    df_test = load_data(benchmark_csv)

    df_train['antibody_id'] = df_train['heavy_id'] + '_' + df_train['light_id']
    df_test['antibody_id'] = df_test['heavy_id'] + '_' + df_test['light_id']

    train_loader_emb = TwoStreamEmbeddingLoader(train_cache_file)
    test_loader_emb = TwoStreamEmbeddingLoader(test_cache_file)

    # Filter missing
    train_antibody_ids = set(df_train['antibody_id'].unique())
    train_antigen_ids = set(df_train['antigen_id'].unique())
    train_valid_antibody = [aid for aid in train_antibody_ids if aid in train_loader_emb.antibody_embeddings]
    train_valid_antigen = [agid for agid in train_antigen_ids if agid in train_loader_emb.antigen_embeddings]
    df_train = df_train[df_train['antibody_id'].isin(train_valid_antibody) & df_train['antigen_id'].isin(train_valid_antigen)]

    test_antibody_ids = set(df_test['antibody_id'].unique())
    test_antigen_ids = set(df_test['antigen_id'].unique())
    test_valid_antibody = [aid for aid in test_antibody_ids if aid in test_loader_emb.antibody_embeddings]
    test_valid_antigen = [agid for agid in test_antigen_ids if agid in test_loader_emb.antigen_embeddings]
    df_test = df_test[df_test['antibody_id'].isin(test_valid_antibody) & df_test['antigen_id'].isin(test_valid_antigen)]

    if len(df_train) == 0 or len(df_test) == 0:
        print(f"  Seed {seed}: No valid samples, skipping.")
        return None

    train_lower = df_train['binding_affinity'].min()
    train_upper = df_train['binding_affinity'].max()
    train_bounds = (train_lower, train_upper)

    test_lower = df_test['binding_affinity'].min()
    test_upper = df_test['binding_affinity'].max()
    test_bounds = (test_lower, test_upper)

    train_loader = DataLoader(
        TwoStreamDataset(df_train, train_loader_emb),
        batch_size=config['batch_size'], shuffle=True,
        collate_fn=collate_fn, num_workers=0, pin_memory=(device == 'cuda')
    )
    test_loader = DataLoader(
        TwoStreamDataset(df_test, test_loader_emb, return_ids=True),
        batch_size=config['batch_size'], shuffle=False,
        collate_fn=collate_fn, num_workers=0, pin_memory=(device == 'cuda')
    )

    model = TwoStreamSymmetric(
        esm_dim=train_loader_emb.embedding_dim,
        projected_size=config['projected_size'],
        num_heads=config['num_heads'],
        dropout=config['dropout'],
        device=device
    ).to(device)

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config['learning_rate'], weight_decay=config['weight_decay']
    )

    print(f"\n  Training seed {seed} on full SAbDab for {config['epochs']} epochs...")
    for epoch in range(config['epochs']):
        train_loss = train_epoch(model, train_loader, optimizer, device, train_bounds)
        print(f"  Seed {seed} Epoch {epoch+1:3d}/{config['epochs']} | Train Loss: {train_loss:.4f}")

    train_metrics, _, _ = evaluate(model, train_loader, device, train_bounds)
    test_metrics, test_preds, test_labels, ab_ids, ag_ids = evaluate(
        model, test_loader, device, test_bounds, return_ids=True
    )

    if save_model and model_path:
        os.makedirs(os.path.dirname(model_path), exist_ok=True)
        torch.save({
            'model_state_dict': model.state_dict(),
            'config': config,
            'train_bounds': train_bounds,
            'test_bounds': test_bounds,
            'train_metrics': train_metrics,
            'test_metrics': test_metrics,
        }, model_path)
        print(f"  Model saved to {model_path}")

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        pd.DataFrame({
            'antibody_id': ab_ids, 'antigen_id': ag_ids,
            'true_affinity': test_labels, 'predicted_affinity': test_preds,
            'error': np.array(test_labels) - np.array(test_preds)
        }).to_csv(os.path.join(output_dir, f'predictions_seed{seed}.csv'), index=False)

    return test_metrics

# ==============================================================================
# BENCHMARK: MULTIPLE SEEDS (matching three‑stream)
# ==============================================================================

def run_benchmark_multiple_seeds_twostream(
    sabdab_csv, benchmark_csv,
    train_cache_file, test_cache_file,
    config, seeds=[314, 114, 144],
    save_models=True, base_model_path='models/twostream_balm_seed',
    output_dir='results_twostream/benchmark'
):
    print("=" * 80)
    print(f"100% SABDAB → BENCHMARK | {len(seeds)} SEEDS | Two‑Stream")
    print("=" * 80)

    all_test_metrics = []
    for run_idx, seed in enumerate(seeds, 1):
        print(f"\n{'#'*90}")
        print(f"RUN {run_idx}/{len(seeds)} — SEED = {seed}")
        print(f"{'#'*90}\n")

        model_path = f"{base_model_path}{seed}.pt"
        seed_output_dir = os.path.join(output_dir, f'seed{seed}')
        os.makedirs(seed_output_dir, exist_ok=True)

        metrics = train_on_sabdab_test_on_benchmark_twostream_single_seed(
            sabdab_csv=sabdab_csv,
            benchmark_csv=benchmark_csv,
            train_cache_file=train_cache_file,
            test_cache_file=test_cache_file,
            config=config.copy(),
            seed=seed,
            save_model=save_models,
            model_path=model_path,
            output_dir=seed_output_dir
        )
        if metrics is not None:
            all_test_metrics.append(metrics)
            print(f"✓ Seed {seed} done — Benchmark Pearson: {metrics['pearson']:.4f}")

    if not all_test_metrics:
        print("No benchmark seeds completed.")
        return None

    # Aggregate
    print(f"\n{'='*90}\nMULTI-SEED BENCHMARK SUMMARY (Two‑Stream)\n{'='*90}")
    summary = {}
    for metric in ['pearson', 'spearman', 'rmse', 'r2', 'loss']:
        vals = [m[metric] for m in all_test_metrics]
        summary[metric] = {'mean': float(np.mean(vals)), 'std': float(np.std(vals))}
        print(f"  {metric.upper():12s}: {summary[metric]['mean']:.4f} ± {summary[metric]['std']:.4f}")

    os.makedirs(output_dir, exist_ok=True)
    pd.DataFrame([{'metric': m, 'mean': v['mean'], 'std': v['std']}
                  for m, v in summary.items()]).to_csv(
        os.path.join(output_dir, 'benchmark_multi_seed_summary.csv'), index=False)
    return summary

# ==============================================================================
# MAIN
# ==============================================================================

def main():
    config = {
        'device': 'cuda' if torch.cuda.is_available() else 'cpu',
        'seed': 9999,
        'projected_size': 256,
        'num_heads': 8,
        'dropout': 0.1,
        'batch_size': 32,
        'epochs': 50,
        'learning_rate': 1e-4,
        'weight_decay': 0.01,
        'patience': 10,
    }

    # Dataset configurations – adjust paths to your two‑stream cache files
    dataset_configs = {
        'sabdab': {
            'csv': 'datasets/pairs_sabdab.csv',
            'cache': 'datasets/esm2_embeddings_twostream_natural_650M.pkl',
            'n_folds': 10,
        },
        'abbind': {
            'csv': 'datasets/pairs_abbind.csv',
            'cache': 'datasets/esm2_embeddings_twostream_abbind_650M.pkl',
            'n_folds': 10,
        },
        'skempi': {
            'csv': 'datasets/pairs_skempi.csv',
            'cache': 'datasets/esm2_embeddings_twostream_skempi_650M.pkl',
            'n_folds': 10,
        },
    }

    print("=" * 80)
    print("TWO-STREAM SYMMETRIC – COMPREHENSIVE EVALUATION")
    print("=" * 80)
    print("\nConfiguration (identical to three‑stream):")
    for k, v in config.items():
        print(f"  {k}: {v}")

    all_summaries = {}

    # --- 10‑fold cross-validation ---
    for name, ds in dataset_configs.items():
        if not os.path.exists(ds['cache']):
            print(f"\n⚠️ Skipping {name}: cache file {ds['cache']} not found.")
            continue
        print(f"\n{'#'*70}\n# CV: {name.upper()} (Two‑Stream)\n{'#'*70}")
        try:
            loader = TwoStreamEmbeddingLoader(ds['cache'])
            df = load_data(ds['csv'])
            df['antibody_id'] = df['heavy_id'] + '_' + df['light_id']
            loader.validate_ids(df, antibody_col='antibody_id', antigen_col='antigen_id')

            output_dir = f"results_twostream/cv_{name}"
            _, summary = cross_validate_twostream(
                df, loader, config,
                n_folds=ds['n_folds'],
                output_dir=output_dir,
                dataset_name=name
            )
            all_summaries[name] = summary
        except Exception as e:
            print(f"  Error on {name}: {e}")
            import traceback; traceback.print_exc()

    # --- Benchmark with multiple seeds ---
    print("\n" + "#"*70)
    print("# BENCHMARK: Train on SAbDab → Test on Benchmark (Two‑Stream, 3 seeds)")
    print("#"*70)

    # You must generate these two caches separately:
    #   python precompute_embeddings_twostream.py --csv datasets/pairs_sabdab.csv --fasta datasets/seq_natural.fasta --output datasets/esm2_embeddings_twostream_sabdab_650M.pkl
    #   python precompute_embeddings_twostream.py --csv datasets/pairs_benchmark.csv --fasta datasets/seq_natural.fasta --output datasets/esm2_embeddings_twostream_benchmark_650M.pkl
    train_cache = 'datasets/esm2_embeddings_twostream_sabdab_650M.pkl'
    test_cache = 'datasets/esm2_embeddings_twostream_benchmark_650M.pkl'

    if os.path.exists(train_cache) and os.path.exists(test_cache):
        try:
            bench_summary = run_benchmark_multiple_seeds_twostream(
                sabdab_csv='datasets/pairs_sabdab.csv',
                benchmark_csv='datasets/pairs_benchmark.csv',
                train_cache_file=train_cache,
                test_cache_file=test_cache,
                config=config,
                seeds=[314, 114, 144],
                save_models=True,
                base_model_path='models/twostream_balm_seed',
                output_dir='results_twostream/benchmark'
            )
            if bench_summary:
                all_summaries['benchmark'] = bench_summary
        except Exception as e:
            print(f"  Benchmark error: {e}")
            import traceback; traceback.print_exc()
    else:
        print(f"⚠️ Benchmark cache files not found: {train_cache} or {test_cache}")
        print("   Please generate them using the precomputation script with the respective CSVs.")

    # --- Final summary ---
    print("\n" + "="*80)
    print("FINAL SUMMARY – TWO-STREAM SYMMETRIC")
    print("="*80)
    for name, summary in all_summaries.items():
        if name == 'benchmark':
            print(f"\nBenchmark (test, 3 seeds):")
            for k, v in summary.items():
                if k in ['pearson', 'spearman', 'rmse', 'r2']:
                    print(f"  {k.upper()}: {v['mean']:.4f} ± {v['std']:.4f}")
        else:
            print(f"\n{name.upper()} CV (10‑fold):")
            for metric in ['pearson', 'spearman', 'rmse', 'r2']:
                if metric in summary:
                    print(f"  {metric.upper()}: {summary[metric]['mean']:.4f} ± {summary[metric]['std']:.4f}")

if __name__ == "__main__":
    main()