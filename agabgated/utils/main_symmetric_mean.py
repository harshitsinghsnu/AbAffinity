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


# ==============================================================================
# UTILITIES
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


# ==============================================================================
# AFFINITY CONVERSION
# ==============================================================================

def affinity_to_cosine(pkd, pkd_lower, pkd_upper):
    """Map affinity linearly into [-1, 1] cosine space."""
    return 2.0 * (pkd - pkd_lower) / (pkd_upper - pkd_lower) - 1.0


def cosine_to_affinity(cosine_sim, pkd_lower, pkd_upper):
    """Inverse: cosine similarity → affinity units."""
    return (cosine_sim + 1.0) / 2.0 * (pkd_upper - pkd_lower) + pkd_lower


# ==============================================================================
# CACHED EMBEDDING LOADER
# ==============================================================================

class CachedEmbeddingLoader:
    def __init__(self, cache_file):
        print(f"  Loading embeddings from {cache_file}")
        with open(cache_file, 'rb') as f:
            self.embeddings = pickle.load(f)
        print(f"  Loaded {len(self.embeddings)} embeddings")
        self.embedding_dim = next(iter(self.embeddings.values())).shape[0]

    def get_embedding(self, seq_id):
        if seq_id not in self.embeddings:
            raise KeyError(
                f"Sequence ID '{seq_id}' not found in embedding cache. "
                f"Re-run precompute_embeddings.py."
            )
        return self.embeddings[seq_id]

    def validate_ids(self, df, id_cols):
        missing = []
        for col in id_cols:
            for sid in df[col].unique():
                if sid not in self.embeddings:
                    missing.append((col, sid))
        if missing:
            raise ValueError(
                f"{len(missing)} sequence IDs missing from cache.\n"
                f"First 5: {missing[:5]}\nRe-run precompute_embeddings.py."
            )
        print(f"  Cache validation OK for {len(id_cols)} ID columns.")


# ==============================================================================
# GATED CROSS-ATTENTION  (paper-correct G1, adapted for mean-pooled embeddings)
# ==============================================================================
#
# Reference: Qiu et al., NeurIPS 2025, Eq. 5
#   Y' = Y ⊙ σ(X_query · W_gate)
#
# Since ESM-2 embeddings are mean-pooled (seq_len=1), standard SDPA is trivial:
#   softmax([B,H,1,1]) = 1.0 always → output = V, W_q/W_k get zero gradients.
#
# Fix: bilinear interaction replaces Q@K^T:
#   Y = tanh(W_q(chain) ⊙ W_k(antigen)) * W_v(antigen)
# Then G1 gate:
#   Y' = Y ⊙ σ(W_gate(chain))
#
# This gives W_q and W_k real gradient signal and captures co-activation
# between chain and antigen embedding dimensions.

class GatedCrossAttention(nn.Module):
    """
    Gated cross-attention for mean-pooled embeddings.
    Bilinear interaction + G1 gate (Qiu et al., NeurIPS 2025).
    """

    def __init__(self, query_dim, kv_dim, num_heads=4, dropout=0.1):
        super().__init__()
        assert query_dim % num_heads == 0
        self.num_heads = num_heads
        self.head_dim  = query_dim // num_heads
        self.query_dim = query_dim

        # Bilinear interaction projections (both → query_dim)
        self.W_q = nn.Linear(query_dim, query_dim, bias=False)
        self.W_k = nn.Linear(kv_dim,   query_dim, bias=False)
        self.W_v = nn.Linear(kv_dim,   query_dim, bias=False)
        self.W_o = nn.Linear(query_dim, query_dim, bias=False)

        # G1 gate: query-driven, elementwise, sigmoid
        self.W_gate = nn.Linear(query_dim, query_dim, bias=True)

        self.dropout    = nn.Dropout(dropout)
        self.layer_norm = nn.LayerNorm(query_dim)

        self._init_weights()

    def _init_weights(self):
        for m in [self.W_q, self.W_k, self.W_v, self.W_o]:
            nn.init.xavier_uniform_(m.weight)
        # Small normal: gates start near 0.5 with slight query-dependence
        nn.init.normal_(self.W_gate.weight, mean=0.0, std=0.01)
        nn.init.zeros_(self.W_gate.bias)

    def forward(self, query, key_value):
        """
        query     : chain embedding  [B, query_dim]
        key_value : antigen embedding [B, kv_dim]
        returns   : gated cross-attended chain repr [B, query_dim]
        """
        residual = query

        # Bilinear interaction: Y = tanh(W_q(chain) ⊙ W_k(antigen)) * W_v(antigen)
        q_proj      = self.W_q(query)              # [B, D]
        k_proj      = self.W_k(key_value)          # [B, D]
        v_proj      = self.W_v(key_value)          # [B, D]
        interaction = torch.tanh(q_proj * k_proj)  # [B, D]  Hadamard
        sdpa_out    = interaction * v_proj          # [B, D]

        # G1 gate: Y' = Y ⊙ σ(X_query · W_gate)
        gate      = torch.sigmoid(self.W_gate(query))  # [B, D]
        gated_out = sdpa_out * gate                     # [B, D]

        # Output projection + residual + LayerNorm
        output = self.W_o(gated_out)
        output = self.dropout(output)
        output = self.layer_norm(residual + output)
        return output


# ==============================================================================
# BASE MODEL
# ==============================================================================

class BaseBALM(nn.Module):
    @staticmethod
    def cosine_to_affinity(c, lo, hi): return cosine_to_affinity(c, lo, hi)
    @staticmethod
    def affinity_to_cosine(p, lo, hi): return affinity_to_cosine(p, lo, hi)


# ==============================================================================
# ASYMMETRIC BALM
# ==============================================================================

class AsymmetricBALM(BaseBALM):
    """
    Biologically motivated: separate weights for heavy and light chains.
    Heavy chain (VH/CDR-H3) dominates antigen contact; light chain is secondary.

    Architecture:
      1. Separate projections: heavy_projection ≠ light_projection
      2. Bidirectional self-attention: heavy ↔ light context exchange
      3. Separate gated cross-attention: heavy→antigen, light→antigen (different weights)
      4. Antigen context:  antigen_fusion(cat[heavy_antigen_ctx, light_antigen_ctx])
      5. Antibody context: antibody_proj(cat[heavy_antigen_ctx, light_antigen_ctx])
         ← DIFFERENT projection from antigen_fusion, so antibody ≠ antigen vector
      6. cosine_similarity(normalize(antibody_ctx), normalize(antigen_ctx))
    """

    def __init__(self, esm_dim=1280, projected_size=256, num_heads=4,
                 dropout=0.1, device='cuda'):
        super().__init__()
        self.projected_size = projected_size

        self.heavy_projection = nn.Sequential(
            nn.Linear(esm_dim, projected_size), nn.LayerNorm(projected_size),
            nn.GELU(), nn.Dropout(dropout)
        )
        self.light_projection = nn.Sequential(
            nn.Linear(esm_dim, projected_size), nn.LayerNorm(projected_size),
            nn.GELU(), nn.Dropout(dropout)
        )
        self.antigen_projection = nn.Sequential(
            nn.Linear(esm_dim, projected_size), nn.LayerNorm(projected_size),
            nn.GELU(), nn.Dropout(dropout)
        )

        self.antibody_self_attention = nn.MultiheadAttention(
            projected_size, num_heads, dropout=dropout, batch_first=True
        )

        # Separate cross-attention modules (asymmetric)
        self.heavy_antigen_attention = GatedCrossAttention(
            projected_size, projected_size, num_heads, dropout)
        self.light_antigen_attention = GatedCrossAttention(
            projected_size, projected_size, num_heads, dropout)

        # Two DIFFERENT projections from the same cross-attended features:
        # antigen_fusion  → represents the antigen side
        # antibody_proj   → represents the antibody side
        # They are DIFFERENT linear maps → different vectors → cosine ≠ 1
        self.antigen_fusion = nn.Sequential(
            nn.Linear(projected_size * 2, projected_size),
            nn.LayerNorm(projected_size), nn.GELU(), nn.Dropout(dropout)
        )
        self.antibody_proj = nn.Sequential(
            nn.Linear(projected_size * 2, projected_size),
            nn.LayerNorm(projected_size), nn.GELU(), nn.Dropout(dropout)
        )

        self._init_weights()

    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)

    def forward(self, light_emb, heavy_emb, antigen_emb):
        light_emb   = F.layer_norm(light_emb,   light_emb.shape[-1:])
        heavy_emb   = F.layer_norm(heavy_emb,   heavy_emb.shape[-1:])
        antigen_emb = F.layer_norm(antigen_emb, antigen_emb.shape[-1:])

        light_proj   = self.light_projection(light_emb)
        heavy_proj   = self.heavy_projection(heavy_emb)
        antigen_proj = self.antigen_projection(antigen_emb)

        # 1. Bidirectional antibody self-attention
        stacked = torch.stack([heavy_proj, light_proj], dim=1)  # [B, 2, D]
        ab_ctx, _ = self.antibody_self_attention(stacked, stacked, stacked)
        heavy_ctx = ab_ctx[:, 0, :]   # [B, D]
        light_ctx = ab_ctx[:, 1, :]   # [B, D]

        # 2. Chain-specific gated cross-attention with antigen
        heavy_ag_ctx = self.heavy_antigen_attention(heavy_ctx, antigen_proj)
        light_ag_ctx = self.light_antigen_attention(light_ctx, antigen_proj)

        # 3. Two DIFFERENT projections of the same cross-attended features
        combined = torch.cat([heavy_ag_ctx, light_ag_ctx], dim=-1)  # [B, 2D]
        antigen_ctx  = self.antigen_fusion(combined)   # [B, D]  — antigen side
        antibody_ctx = self.antibody_proj(combined)    # [B, D]  — antibody side
        # antigen_ctx ≠ antibody_ctx because antigen_fusion ≠ antibody_proj

        # 4. Cosine similarity
        ab_norm  = F.normalize(antibody_ctx, p=2, dim=-1)
        ag_norm  = F.normalize(antigen_ctx,  p=2, dim=-1)
        cosine_sim = F.cosine_similarity(ab_norm, ag_norm, dim=-1)

        return {'cosine_similarity': cosine_sim,
                'antibody_context': ab_norm, 'antigen_context': ag_norm}


# ==============================================================================
# SYMMETRIC BALM
# ==============================================================================

class Symmetric(BaseBALM):
    """
    Shared weights for heavy and light chains (symmetric treatment).
    Swapping light/heavy gives identical output for mean/weighted fusion.

    CRITICAL FIX: antibody_ctx and antigen_ctx are computed via TWO DIFFERENT
    projection heads from the same cross-attended features. This ensures
    cosine_similarity ≠ 1.0 always.

    Fusion methods:
      'mean'     : (light_ag + heavy_ag) / 2  → then two different projections
      'weighted' : learned normalized weighted sum → then two different projections
      'concat'   : cat → two different linear projections
    """

    def __init__(self, esm_dim=1280, projected_size=256, num_heads=4,
                 dropout=0.1, fusion_method='mean', device='cuda'):
        super().__init__()
        self.projected_size = projected_size
        self.fusion_method  = fusion_method

        self.chain_projection = nn.Sequential(
            nn.Linear(esm_dim, projected_size), nn.LayerNorm(projected_size),
            nn.GELU(), nn.Dropout(dropout)
        )
        self.antigen_projection = nn.Sequential(
            nn.Linear(esm_dim, projected_size), nn.LayerNorm(projected_size),
            nn.GELU(), nn.Dropout(dropout)
        )

        self.antibody_self_attention = nn.MultiheadAttention(
            projected_size, num_heads, dropout=dropout, batch_first=True
        )

        # Shared cross-attention (symmetric)
        self.chain_antigen_attention = GatedCrossAttention(
            projected_size, projected_size, num_heads, dropout)

        # Fusion-specific layers
        if fusion_method == 'weighted':
            self.fusion_gate = nn.Linear(projected_size * 2, projected_size * 2)
            in_dim = projected_size
        elif fusion_method == 'concat':
            in_dim = projected_size * 2
        else:  # mean
            in_dim = projected_size


        self.antigen_head = nn.Sequential(
            nn.Linear(in_dim, projected_size),
            nn.LayerNorm(projected_size), nn.GELU(), nn.Dropout(dropout)
        )
        self.antibody_head = nn.Sequential(
            nn.Linear(in_dim, projected_size),
            nn.LayerNorm(projected_size), nn.GELU(), nn.Dropout(dropout)
        )

        self._init_weights()

    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)

    def forward(self, light_emb, heavy_emb, antigen_emb):
        light_emb   = F.layer_norm(light_emb,   light_emb.shape[-1:])
        heavy_emb   = F.layer_norm(heavy_emb,   heavy_emb.shape[-1:])
        antigen_emb = F.layer_norm(antigen_emb, antigen_emb.shape[-1:])

        light_proj   = self.chain_projection(light_emb)
        heavy_proj   = self.chain_projection(heavy_emb)
        antigen_proj = self.antigen_projection(antigen_emb)

        # 1. Bidirectional antibody self-attention
        stacked = torch.stack([light_proj, heavy_proj], dim=1)
        ab_ctx, _ = self.antibody_self_attention(stacked, stacked, stacked)
        light_ctx = ab_ctx[:, 0, :]
        heavy_ctx = ab_ctx[:, 1, :]

        # 2. Shared gated cross-attention
        light_ag_ctx = self.chain_antigen_attention(light_ctx, antigen_proj)
        heavy_ag_ctx = self.chain_antigen_attention(heavy_ctx, antigen_proj)

        # 3. Fuse cross-attended features
        if self.fusion_method == 'mean':
            fused = (light_ag_ctx + heavy_ag_ctx) / 2.0          # [B, D]

        elif self.fusion_method == 'weighted':
            combined = torch.cat([light_ag_ctx, heavy_ag_ctx], dim=-1)
            gw = torch.sigmoid(self.fusion_gate(combined))
            lw = gw[:, :self.projected_size]
            hw = gw[:, self.projected_size:]
            total = lw + hw + 1e-8
            fused = (lw / total) * light_ag_ctx + (hw / total) * heavy_ag_ctx

        elif self.fusion_method == 'concat':
            fused = torch.cat([light_ag_ctx, heavy_ag_ctx], dim=-1)  # [B, 2D]

        else:
            raise ValueError(f"Unknown fusion_method: '{self.fusion_method}'")

        # 4. TWO DIFFERENT heads → antibody_ctx ≠ antigen_ctx
        antigen_ctx  = self.antigen_head(fused)   # [B, D]
        antibody_ctx = self.antibody_head(fused)  # [B, D]  ← different weights!

        # 5. Cosine similarity
        ab_norm  = F.normalize(antibody_ctx, p=2, dim=-1)
        ag_norm  = F.normalize(antigen_ctx,  p=2, dim=-1)
        cosine_sim = F.cosine_similarity(ab_norm, ag_norm, dim=-1)

        return {'cosine_similarity': cosine_sim,
                'antibody_context': ab_norm, 'antigen_context': ag_norm}


# ==============================================================================
# MODEL FACTORY
# ==============================================================================

def get_model(model_name, **kwargs):
    """
    Available: 'asymmetric', 'symmetric_mean', 'symmetric_weighted', 'symmetric_concat'
    """
    kwargs.pop('gate_type', None)
    kwargs.pop('use_gated_attention', None)

    registry = {
        'asymmetric':         lambda **kw: AsymmetricBALM(**kw),
        'symmetric_mean':     lambda **kw: Symmetric(fusion_method='mean',     **kw),
        'symmetric_weighted': lambda **kw: Symmetric(fusion_method='weighted', **kw),
        'symmetric_concat':   lambda **kw: Symmetric(fusion_method='concat',   **kw),
    }
    if model_name not in registry:
        raise ValueError(f"Unknown model '{model_name}'. Available: {list(registry)}")
    return registry[model_name](**kwargs)


# ==============================================================================
# DATASET
# ==============================================================================

class CachedEmbeddingDataset(Dataset):
    def __init__(self, df, embedding_loader, return_ids=False):
        self.df         = df.reset_index(drop=True)
        self.loader     = embedding_loader
        self.return_ids = return_ids

    def __len__(self): return len(self.df)

    def __getitem__(self, idx):
        row  = self.df.iloc[idx]
        item = {
            'light_emb':   self.loader.get_embedding(row['light_id']),
            'heavy_emb':   self.loader.get_embedding(row['heavy_id']),
            'antigen_emb': self.loader.get_embedding(row['antigen_id']),
            'affinity':    row['binding_affinity']
        }
        if self.return_ids:
            item['light_id']   = row['light_id']
            item['heavy_id']   = row['heavy_id']
            item['antigen_id'] = row['antigen_id']
        return item


def collate_fn(batch):
    out = {
        'light_emb':   torch.tensor(np.array([b['light_emb']   for b in batch]), dtype=torch.float32),
        'heavy_emb':   torch.tensor(np.array([b['heavy_emb']   for b in batch]), dtype=torch.float32),
        'antigen_emb': torch.tensor(np.array([b['antigen_emb'] for b in batch]), dtype=torch.float32),
        'affinity':    torch.tensor([b['affinity'] for b in batch], dtype=torch.float32),
    }
    if 'light_id' in batch[0]:
        out['light_id']   = [b['light_id']   for b in batch]
        out['heavy_id']   = [b['heavy_id']   for b in batch]
        out['antigen_id'] = [b['antigen_id'] for b in batch]
    return out


# ==============================================================================
# TRAINING
# ==============================================================================

def train_epoch(model, loader, optimizer, device, pkd_bounds):
    model.train()
    total_loss = 0.0
    pkd_lower, pkd_upper = pkd_bounds
    for batch in loader:
        le = batch['light_emb'].to(device)
        he = batch['heavy_emb'].to(device)
        ae = batch['antigen_emb'].to(device)
        af = batch['affinity'].to(device)

        target = affinity_to_cosine(af, pkd_lower, pkd_upper)
        out    = model(le, he, ae)
        loss   = F.mse_loss(out['cosine_similarity'], target)

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        total_loss += loss.item()
    return total_loss / len(loader)


def evaluate(model, loader, device, pkd_bounds, return_ids=False):
    """
    Evaluate model.

    pkd_bounds: used for cosine→affinity conversion of PREDICTIONS.
    The labels are returned in their original affinity units (not converted).
    This means pkd_bounds must match the affinity scale of the dataset being evaluated.
    """
    model.eval()
    preds, labels = [], []
    l_ids, h_ids, ag_ids = [], [], []
    total_loss = 0.0
    pkd_lower, pkd_upper = pkd_bounds

    with torch.no_grad():
        for batch in loader:
            le = batch['light_emb'].to(device)
            he = batch['heavy_emb'].to(device)
            ae = batch['antigen_emb'].to(device)
            af = batch['affinity'].to(device)

            target = affinity_to_cosine(af, pkd_lower, pkd_upper)
            out    = model(le, he, ae)
            loss   = F.mse_loss(out['cosine_similarity'], target)
            total_loss += loss.item()

            pred_pkd = cosine_to_affinity(out['cosine_similarity'], pkd_lower, pkd_upper)
            preds.extend(pred_pkd.cpu().numpy().tolist())
            labels.extend(af.cpu().numpy().tolist())

            if return_ids and 'light_id' in batch:
                l_ids.extend(batch['light_id'])
                h_ids.extend(batch['heavy_id'])
                ag_ids.extend(batch['antigen_id'])

    metrics = {
        'loss':     total_loss / len(loader),
        'r2':       r_squared(labels, preds),
        'pearson':  pearson_correlation(labels, preds),
        'spearman': spearman_correlation(labels, preds),
        'rmse':     rmse(labels, preds),
    }
    if return_ids:
        return metrics, preds, labels, l_ids, h_ids, ag_ids
    return metrics, preds, labels


def train_model(df_train, df_val, embedding_loader, config, device, model_name, save_path=None):
    """Train with validation-Pearson early stopping. pkd_bounds from train only.
    If save_path is provided, save the best model state dict to that path.
    """
    pkd_lower  = df_train['binding_affinity'].min()
    pkd_upper  = df_train['binding_affinity'].max()
    pkd_bounds = (pkd_lower, pkd_upper)
    print(f"  Affinity range (train): [{pkd_lower:.3f}, {pkd_upper:.3f}]")

    train_loader = DataLoader(
        CachedEmbeddingDataset(df_train, embedding_loader),
        batch_size=config['batch_size'], shuffle=True, collate_fn=collate_fn,
        num_workers=0, pin_memory=(device == 'cuda')
    )
    val_loader = DataLoader(
        CachedEmbeddingDataset(df_val, embedding_loader),
        batch_size=config['batch_size'], collate_fn=collate_fn,
        num_workers=0, pin_memory=(device == 'cuda')
    )

    model = get_model(
        model_name, esm_dim=embedding_loader.embedding_dim,
        projected_size=config['projected_size'], num_heads=config['num_heads'],
        dropout=config['dropout'], device=device
    ).to(device)

    n = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"  Trainable parameters: {n:,}")

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config['learning_rate'],
        weight_decay=config['weight_decay']
    )

    best_pearson, best_state, patience = -float('inf'), None, 0

    for epoch in range(config['epochs']):
        train_loss = train_epoch(model, train_loader, optimizer, device, pkd_bounds)
        val_metrics, _, _ = evaluate(model, val_loader, device, pkd_bounds)

        print(f"  Epoch {epoch+1:3d}/{config['epochs']} | "
              f"Loss: {train_loss:.4f} | "
              f"Val Pearson: {val_metrics['pearson']:.4f} | "
              f"Val R²: {val_metrics['r2']:.4f} | "
              f"Val RMSE: {val_metrics['rmse']:.4f}")

        if val_metrics['pearson'] > best_pearson:
            best_pearson = val_metrics['pearson']
            best_state   = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            patience     = 0
        else:
            patience += 1
            if patience >= config['patience']:
                print(f"  Early stopping (best val Pearson: {best_pearson:.4f})")
                break

    if best_state:
        model.load_state_dict({k: v.to(device) for k, v in best_state.items()})
        # Save the best model if a path is provided
        if save_path:
            os.makedirs(os.path.dirname(save_path), exist_ok=True)
            torch.save({
                'model_state_dict': best_state,  # already on CPU
                'config': config,
                'model_name': model_name,
                'pkd_bounds': pkd_bounds,
                'val_metrics': val_metrics,
            }, save_path)
            print(f"  Best model saved to {save_path}")

    return model, pkd_bounds


# ==============================================================================
# CROSS-VALIDATION
# ==============================================================================

def cross_validate(df, embedding_loader, config, n_folds=10,
                   output_dir='results', model_name='symmetric_mean',
                   save_models=False):   # new parameter
    kfold = KFold(n_splits=n_folds, shuffle=True, random_state=config['seed'])
    fold_results = []
    all_preds    = []
    device       = config['device']
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{n_folds}-fold CV | model={model_name} | n={len(df)}")

    for fold, (tr_idx, va_idx) in enumerate(kfold.split(df), 1):
        print(f"\n{'='*70}\nFOLD {fold}/{n_folds}\n{'='*70}")
        df_tr = df.iloc[tr_idx].reset_index(drop=True)
        df_va = df.iloc[va_idx].reset_index(drop=True)
        print(f"  Train: {len(df_tr)} | Test (held-out fold): {len(df_va)}")

        # Prepare model save path if requested
        model_save_path = None
        if save_models:
            model_save_path = os.path.join(output_dir, f'fold_{fold:02d}', 'model.pt')

        model, pkd_bounds = train_model(
            df_tr, df_va, embedding_loader, config, device, model_name,
            save_path=model_save_path   # pass the path
        )

        val_loader = DataLoader(
            CachedEmbeddingDataset(df_va, embedding_loader, return_ids=True),
            batch_size=config['batch_size'], collate_fn=collate_fn
        )
        metrics, preds, labels, l_ids, h_ids, ag_ids = evaluate(
            model, val_loader, device, pkd_bounds, return_ids=True)

        fold_dir = os.path.join(output_dir, f'fold_{fold:02d}')
        os.makedirs(fold_dir, exist_ok=True)

        fold_df = pd.DataFrame({
            'light_id': l_ids, 'heavy_id': h_ids, 'antigen_id': ag_ids,
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

    print(f"\n{'='*70}\nCV SUMMARY — {model_name}\n{'='*70}")
    for metric, s in summary.items():
        print(f"  {metric.upper():10s}: {s['mean']:.4f} ± {s['std']:.4f}")

    pd.DataFrame([{'metric': k, 'mean': v['mean'], 'std': v['std']}
                  for k, v in summary.items()]).to_csv(
        os.path.join(output_dir, 'cv_summary.csv'), index=False)

    return fold_results, summary


# ==============================================================================
# FULL TRAIN ON SABDAB + TEST ON BENCHMARK
# ==============================================================================

def train_on_sabdab_test_on_benchmark(
    sabdab_csv    = 'datasets/pairs_sabdab.csv',
    benchmark_csv = 'datasets/pairs_benchmark.csv',
    cache_file    = 'datasets/esm2_embeddings_natural_650M.pkl',
    config        = None,
    model_name    = 'symmetric_mean',
    save_model    = True,
    model_path    = 'models/balm_sabdab_final.pt',
    output_dir    = 'results/benchmark',
    
):
    embedding_loader = CachedEmbeddingLoader(cache_file)
    """
    100% TRAIN ON SABDAB + TEST ON BENCHMARK (as requested).
    No validation split / no early stopping — full epochs on entire SAbDab.
    """
    if config is None:
        config = DEFAULT_CONFIG.copy()

    print("=" * 80)
    print(f"100% TRAIN ON SABDAB + BENCHMARK TEST | Model: {model_name}")
    print("=" * 80)
    for k, v in config.items():
        print(f"  {k}: {v}")

    setup_reproducibility(config['seed'])
    device = config['device']

    embedding_loader = CachedEmbeddingLoader(cache_file)

    df_sabdab    = load_data(sabdab_csv)
    df_benchmark = load_data(benchmark_csv)

    print(f"\n  SAbDab:    {len(df_sabdab)} samples")
    print(f"  Benchmark: {len(df_benchmark)} samples")

    embedding_loader.validate_ids(df_sabdab,    ['light_id', 'heavy_id', 'antigen_id'])
    embedding_loader.validate_ids(df_benchmark, ['light_id', 'heavy_id', 'antigen_id'])

    # === 100% SAbDab for training ===
    df_train = df_sabdab.reset_index(drop=True)
    sabdab_lower = df_train['binding_affinity'].min()
    sabdab_upper = df_train['binding_affinity'].max()
    sabdab_bounds = (sabdab_lower, sabdab_upper)

    bench_lower = df_benchmark['binding_affinity'].min()
    bench_upper = df_benchmark['binding_affinity'].max()
    bench_bounds = (bench_lower, bench_upper)

    print(f"  SAbDab affinity bounds:   [{sabdab_lower:.3f}, {sabdab_upper:.3f}]")
    print(f"  Benchmark affinity bounds: [{bench_lower:.3f}, {bench_upper:.3f}]")

    train_loader = DataLoader(
        CachedEmbeddingDataset(df_train, embedding_loader),
        batch_size=config['batch_size'], shuffle=True, collate_fn=collate_fn,
        num_workers=0, pin_memory=(device == 'cuda')
    )

    benchmark_loader = DataLoader(
        CachedEmbeddingDataset(df_benchmark, embedding_loader, return_ids=True),
        batch_size=config['batch_size'], collate_fn=collate_fn
    )

    model = get_model(
        model_name, esm_dim=embedding_loader.embedding_dim,
        projected_size=config['projected_size'], num_heads=config['num_heads'],
        dropout=config['dropout'], device=device
    ).to(device)

    print(f"\n  Trainable parameters: {sum(p.numel() for p in model.parameters() if p.requires_grad):,}")

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config['learning_rate'], weight_decay=config['weight_decay']
    )

    print(f"\n  Training for {config['epochs']} epochs on 100% SAbDab (no early stopping)...")

    for epoch in range(config['epochs']):
        train_loss = train_epoch(model, train_loader, optimizer, device, sabdab_bounds)
        print(f"  Epoch {epoch+1:3d}/{config['epochs']} | Train Loss: {train_loss:.4f}")

    # Final evaluation on full SAbDab train set
    print(f"\n{'='*60}\n  SABDAB FULL TRAINING SET EVALUATION\n{'='*60}")
    train_eval_loader = DataLoader(
        CachedEmbeddingDataset(df_train, embedding_loader),
        batch_size=config['batch_size'], collate_fn=collate_fn
    )
    train_metrics, _, _ = evaluate(model, train_eval_loader, device, sabdab_bounds)
    for k, v in train_metrics.items():
        print(f"    {k.upper():12s}: {v:.4f}")

    # Benchmark evaluation (uses benchmark bounds)
    print(f"\n{'='*60}\n  BENCHMARK TEST SET EVALUATION\n{'='*60}")
    print(f"  Using benchmark bounds [{bench_lower:.3f}, {bench_upper:.3f}]")

    test_metrics, test_preds, test_labels, t_l, t_h, t_ag = evaluate(
        model, benchmark_loader, device, bench_bounds, return_ids=True
    )
    for k, v in test_metrics.items():
        print(f"    {k.upper():12s}: {v:.4f}")

    # Save model
    if save_model:
        os.makedirs(os.path.dirname(model_path) or '.', exist_ok=True)
        torch.save({
            'model_state_dict': model.state_dict(),
            'config': config,
            'model_name': model_name,
            'sabdab_bounds': sabdab_bounds,
            'bench_bounds': bench_bounds,
            'train_metrics': train_metrics,
            'test_metrics': test_metrics,
        }, model_path)
        print(f"\n  Model saved to {model_path}")

    # Save results
    test_dir = os.path.join(output_dir, model_name)
    os.makedirs(test_dir, exist_ok=True)

    pd.DataFrame({
        'light_id': t_l, 'heavy_id': t_h, 'antigen_id': t_ag,
        'true_affinity': test_labels, 'predicted_affinity': test_preds,
        'error': np.array(test_labels) - np.array(test_preds)
    }).to_csv(os.path.join(test_dir, 'predictions.csv'), index=False)

    pd.DataFrame([{
        'dataset': 'sabdab_full', 'n_samples': len(df_train), **train_metrics
    }, {
        'dataset': 'benchmark', 'n_samples': len(df_benchmark), **test_metrics
    }]).to_csv(os.path.join(test_dir, 'metrics.csv'), index=False)

    return {
        'model': model,
        'train_metrics': train_metrics,
        'test_metrics': test_metrics,
        'test_predictions': test_preds,
        'test_labels': test_labels,
        'sabdab_bounds': sabdab_bounds,
        'bench_bounds': bench_bounds,
        'config': config,
    }


# ==============================================================================
# DATA LOADING
# ==============================================================================

def load_data(pairs_csv):
    print(f"  Loading {pairs_csv}")
    df = pd.read_csv(pairs_csv)

    light_col = heavy_col = antigen_col = affinity_col = None
    for col in df.columns:
        c = col.lower().strip()
        if   'light'   in c and light_col    is None: light_col    = col
        elif 'heavy'   in c and heavy_col    is None: heavy_col    = col
        elif 'antigen' in c and antigen_col  is None: antigen_col  = col
        elif c in ('y', 'affinity', 'binding_affinity', 'pkd', 
                   'ddg', 'dg', 'kd', 'ic50') and affinity_col is None:
            affinity_col = col

    missing = []
    if light_col    is None: missing.append("light chain ID (col with 'light')")
    if heavy_col    is None: missing.append("heavy chain ID (col with 'heavy')")
    if antigen_col  is None: missing.append("antigen ID (col with 'antigen')")
    if affinity_col is None: missing.append(
        "affinity (col named: y, affinity, binding_affinity, pkd,  ddg, dg, kd, ic50)")
    if missing:
        raise ValueError(
            f"Missing columns in {pairs_csv}:\n  {missing}\n"
            f"Available: {list(df.columns)}"
        )

    df = df.rename(columns={
        light_col: 'light_id', heavy_col: 'heavy_id',
        antigen_col: 'antigen_id', affinity_col: 'binding_affinity'
    })
    n_before = len(df)
    df = df.dropna(subset=['light_id', 'heavy_id', 'antigen_id', 'binding_affinity'])
    if len(df) < n_before:
        print(f"  Dropped {n_before - len(df)} rows with NaN.")

    print(f"  Loaded {len(df)} pairs | "
          f"affinity [{df['binding_affinity'].min():.3f}, {df['binding_affinity'].max():.3f}]")
    return df


# ==============================================================================
# DEFAULT CONFIG & DATASET CONFIGS
# ==============================================================================

DEFAULT_CONFIG = {
    'device':         'cuda' if torch.cuda.is_available() else 'cpu',
    'seed':           9999,
    'projected_size': 256,
    'num_heads':      8,
    'dropout':        0.1,
    'batch_size':     32,
    'epochs':         50,
    'learning_rate':  1e-4,
    'weight_decay':   0.01,
    'patience':       10,
}

DATASET_CONFIGS = {
    'sabdab': {
        'pairs_csv':  'datasets/pairs_sabdab.csv',
        'cache_file': 'datasets/esm2_embeddings_natural_650M.pkl',
        'n_folds':    10,
    },
    'abbind': {
        'pairs_csv':  'datasets/pairs_abbind.csv',
        'cache_file': 'datasets/esm2_embeddings_mutation_650M.pkl',
        'n_folds':    10,
    },
    'skempi': {
        'pairs_csv':  'datasets/pairs_skempi.csv',
        'cache_file': 'datasets/esm2_embeddings_mutation_650M.pkl',
        'n_folds':    10,
    },
}


# ==============================================================================
# CONVENIENCE ENTRY POINTS
# ==============================================================================

def run_cross_validation(dataset_name='sabdab', model_name='symmetric_mean', config=None, save_models=False):
    if config is None: config = DEFAULT_CONFIG.copy()
    ds = DATASET_CONFIGS[dataset_name]
    setup_reproducibility(config['seed'])
    loader = CachedEmbeddingLoader(ds['cache_file'])
    df     = load_data(ds['pairs_csv'])
    loader.validate_ids(df, ['light_id', 'heavy_id', 'antigen_id'])
    return cross_validate(
        df, loader, config, n_folds=ds['n_folds'],
        output_dir=f"results/{dataset_name}_cv_{model_name}",
        model_name=model_name,
        save_models=save_models   # pass the flag
    )


def run_all_cv(model_name='symmetric_mean', config=None, save_models=False):
    if config is None: config = DEFAULT_CONFIG.copy()
    summaries = {}
    for ds_name in ['sabdab', 'abbind', 'skempi']:
        print(f"\n{'#'*70}\n# CV: {ds_name.upper()} | {model_name}\n{'#'*70}")
        try:
            _, summary = run_cross_validation(ds_name, model_name, config, save_models=save_models)
            summaries[ds_name] = summary
        except FileNotFoundError as e:
            print(f"  Skipping {ds_name}: {e}")
        except Exception:
            import traceback; traceback.print_exc()

    print(f"\n{'='*70}\nALL DATASETS CV SUMMARY\n{'='*70}")
    rows = []
    for ds, s in summaries.items():
        row = {'dataset': ds}
        for m, v in s.items():
            row[f'{m}_mean'] = v['mean']; row[f'{m}_std'] = v['std']
        rows.append(row)
        print(f"\n  {ds.upper()}:")
        for m, v in s.items():
            print(f"    {m.upper():10s}: {v['mean']:.4f} ± {v['std']:.4f}")

    if rows:
        os.makedirs('results', exist_ok=True)
        pd.DataFrame(rows).to_csv(f'results/all_cv_summary_{model_name}.csv', index=False)
    return summaries


def run_benchmark(model_name='symmetric_mean', config=None):
    if config is None: config = DEFAULT_CONFIG.copy()
    return train_on_sabdab_test_on_benchmark(
        sabdab_csv='datasets/pairs_sabdab.csv',
        benchmark_csv='datasets/pairs_benchmark.csv',
        cache_file='datasets/esm2_embeddings_natural_650M.pkl',
        config=config, model_name=model_name, save_model=True,
        model_path=f'models/balm_{model_name}_final.pt',
        output_dir='results/benchmark'
    )

# ==============================================================================
# ADD THIS NEW FUNCTION TO YOUR CODE
# (Insert it right after the `run_benchmark` function definition)
# ==============================================================================

def run_benchmark_multiple_seeds(
    model_name='symmetric_mean',
    seeds=None,
    config=None,
    save_all_models=True
):
    if seeds is None:
        seeds = [314, 114, 144]     # Change seeds here if needed

    if config is None:
        config = DEFAULT_CONFIG.copy()

    print("=" * 90)
    print(f"100% SABDAB → BENCHMARK | {len(seeds)} SEEDS | Model: {model_name}")
    print("=" * 90)

    all_test_metrics = []

    for run_idx, seed in enumerate(seeds, 1):
        print(f"\n{'#'*90}")
        print(f"RUN {run_idx}/{len(seeds)} — SEED = {seed}")
        print(f"{'#'*90}\n")

        config_seed = config.copy()
        config_seed['seed'] = seed

        model_path = f"models/balm_{model_name}_seed{seed}.pt"
        output_dir = f"results/benchmark_seed{seed}"

        result = train_on_sabdab_test_on_benchmark(
            config=config_seed,
            model_name=model_name,
            save_model=save_all_models,
            model_path=model_path,
            output_dir=output_dir
        )

        all_test_metrics.append(result['test_metrics'])
        print(f"✓ Seed {seed} done — Benchmark Pearson: {result['test_metrics']['pearson']:.4f}")

    # Aggregate
    print(f"\n{'='*90}\nMULTI-SEED BENCHMARK SUMMARY (100% SAbDab training)\n{'='*90}")
    metrics_to_report = ['pearson', 'spearman', 'r2', 'rmse', 'loss']
    summary = {}

    for metric in metrics_to_report:
        vals = [r[metric] for r in all_test_metrics]
        mean_val = float(np.mean(vals))
        std_val  = float(np.std(vals, ddof=0))
        summary[metric] = {'mean': mean_val, 'std': std_val}
        print(f"  {metric.upper():12s}: {mean_val:.4f} ± {std_val:.4f}")

    os.makedirs('results', exist_ok=True)
    pd.DataFrame([{'metric': m, 'mean': v['mean'], 'std': v['std']} 
                  for m, v in summary.items()]).to_csv(
        'results/benchmark_multi_seed_summary.csv', index=False)

    print(f"\nSummary saved to: results/benchmark_multi_seed_summary.csv")
    return summary, all_test_metrics


# ==============================================================================
# UPDATE THE CLI SECTION (replace the entire if __name__ == "__main__" block)
# ==============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='AgAbGated — Binding Affinity model',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

        """
    )
    parser.add_argument('--save_cv_models', action='store_true',
                    help='Save best model weights from each CV fold')
    parser.add_argument('--mode',    choices=['cv','cv_all','benchmark','all'], default='cv')
    parser.add_argument('--dataset', choices=['sabdab','abbind','skempi'],      default='sabdab')
    parser.add_argument('--model',   choices=['asymmetric','symmetric_mean',
                                              'symmetric_weighted','symmetric_concat'],
                        default='symmetric_mean')
    parser.add_argument('--multi-seed', action='store_true',
                        help='Run benchmark mode for 3 different seeds (recommended)')
    parser.add_argument('--projected_size', type=int,   default=256)
    parser.add_argument('--num_heads',      type=int,   default=8)
    parser.add_argument('--dropout',        type=float, default=0.1)
    parser.add_argument('--batch_size',     type=int,   default=32)
    parser.add_argument('--epochs',         type=int,   default=50)
    parser.add_argument('--lr',             type=float, default=1e-4)
    parser.add_argument('--weight_decay',   type=float, default=0.01)
    parser.add_argument('--patience',       type=int,   default=10)
    parser.add_argument('--seed',           type=int,   default=9999)

    args   = parser.parse_args()
    
    config = {
        'device':         'cuda' if torch.cuda.is_available() else 'cpu',
        'seed':           args.seed,
        'projected_size': args.projected_size,
        'num_heads':      args.num_heads,
        'dropout':        args.dropout,
        'batch_size':     args.batch_size,
        'epochs':         args.epochs,
        'learning_rate':  args.lr,
        'weight_decay':   args.weight_decay,
        'patience':       args.patience,
    }

    if args.mode == 'cv':
        run_cross_validation(args.dataset, args.model, config, save_models=args.save_cv_models)
    elif args.mode == 'cv_all':
        run_all_cv(args.model, config, save_models=args.save_cv_models)
    elif args.mode == 'benchmark':
        if args.multi_seed:
            run_benchmark_multiple_seeds(model_name=args.model, config=config)
        else:
            run_benchmark(args.model, config)
    elif args.mode == 'all':
        run_all_cv(args.model, config, save_models=args.save_cv_models)
        run_benchmark(args.model, config)


# ==============================================================================
# CLI
# ==============================================================================

# if __name__ == "__main__":
#     import argparse

#     parser = argparse.ArgumentParser(
#         description='AgAbGated — Binding Affinity model',
#         formatter_class=argparse.RawDescriptionHelpFormatter,
#         epilog="""
# Examples:

#         """
#     )
#     parser.add_argument('--mode',    choices=['cv','cv_all','benchmark','all'], default='cv')
#     parser.add_argument('--dataset', choices=['sabdab','abbind','skempi'],      default='sabdab')
#     parser.add_argument('--model',   choices=['asymmetric','symmetric_mean',
#                                               'symmetric_weighted','symmetric_concat'],
#                         default='symmetric_mean')
#     parser.add_argument('--projected_size', type=int,   default=256)
#     parser.add_argument('--num_heads',      type=int,   default=8)
#     parser.add_argument('--dropout',        type=float, default=0.1)
#     parser.add_argument('--batch_size',     type=int,   default=32)
#     parser.add_argument('--epochs',         type=int,   default=50)
#     parser.add_argument('--lr',             type=float, default=1e-4)
#     parser.add_argument('--weight_decay',   type=float, default=0.01)
#     parser.add_argument('--patience',       type=int,   default=10)
#     parser.add_argument('--seed',           type=int,   default=9999)

#     args   = parser.parse_args()
#     config = {
#         'device':         'cuda' if torch.cuda.is_available() else 'cpu',
#         'seed':           args.seed,
#         'projected_size': args.projected_size,
#         'num_heads':      args.num_heads,
#         'dropout':        args.dropout,
#         'batch_size':     args.batch_size,
#         'epochs':         args.epochs,
#         'learning_rate':  args.lr,
#         'weight_decay':   args.weight_decay,
#         'patience':       args.patience,
#     }

#     if args.mode == 'cv':
#         run_cross_validation(args.dataset, args.model, config, save_models=args.save_cv_models)
#     elif args.mode == 'cv_all':
#         run_all_cv(args.model, config, save_models=args.save_cv_models)
#     elif args.mode == 'benchmark':
#         if args.multi_seed:
#             run_benchmark_multiple_seeds(model_name=args.model, config=config)
#         else:
#             run_benchmark(args.model, config)
#     elif args.mode == 'all':
#         run_all_cv(args.model, config, save_models=args.save_cv_models)
#         run_benchmark(args.model, config)
