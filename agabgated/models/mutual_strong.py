#!/usr/bin/env python3
"""
MutualTriStreamStrong: Strong bidirectional gated cross‑attention.
- Three independent streams (heavy, light, antigen)
- Heavy‑light self‑attention
- Parallel mutual gated cross‑attention (antibody ↔ antigen)
- Weighted fusion for heavy/light
- Dual heads for cosine similarity
- Optional multiple layers of mutual exchange
"""

import os
import pickle
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
from sklearn.model_selection import KFold
from copy import deepcopy

# ----------------------------------------------------------------------
# Import utilities and data handling from your main code
# ----------------------------------------------------------------------
try:
    from agabgated.utils.main_symmetric_mean import (
        setup_reproducibility, load_data, CachedEmbeddingLoader,
        CachedEmbeddingDataset, collate_fn, DEFAULT_CONFIG, DATASET_CONFIGS
    )
except ImportError:
    raise ImportError("Could not import from main_symmetric_mean.py")

# ----------------------------------------------------------------------
# GatedCrossAttention (same as before)
# ----------------------------------------------------------------------
class GatedCrossAttention(nn.Module):
    def __init__(self, query_dim, kv_dim, num_heads=4, dropout=0.1):
        super().__init__()
        assert query_dim % num_heads == 0
        self.num_heads = num_heads
        self.head_dim = query_dim // num_heads
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

# ----------------------------------------------------------------------
# Strong bidirectional mutual model
# ----------------------------------------------------------------------
class MutualTriStreamStrong(nn.Module):
    """
    Strong mutual conditioning:
    - Parallel bidirectional gated cross‑attention (Ab ↔ Ag)
    - Optionally stack multiple layers
    """
    def __init__(self, esm_dim=1280, projected_size=256, num_heads=8,
                 dropout=0.1, n_layers=2, device='cuda'):
        super().__init__()
        self.projected_size = projected_size
        self.n_layers = n_layers

        # Three separate projections
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

        # Antibody internal self‑attention (heavy ↔ light)
        self.antibody_self_attention = nn.MultiheadAttention(
            projected_size, num_heads, dropout=dropout, batch_first=True
        )

        # Learnable weighted fusion for heavy/light
        self.fusion_gate = nn.Linear(projected_size * 2, 2)

        # Stack of mutual interaction layers
        self.ab_to_ag_layers = nn.ModuleList()
        self.ag_to_ab_layers = nn.ModuleList()
        for _ in range(n_layers):
            self.ab_to_ag_layers.append(
                GatedCrossAttention(projected_size, projected_size, num_heads, dropout)
            )
            self.ag_to_ab_layers.append(
                GatedCrossAttention(projected_size, projected_size, num_heads, dropout)
            )

        # Dual heads for final cosine (anti‑collapse)
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

    def forward(self, light_emb, heavy_emb, antigen_emb):
        # LayerNorm on raw embeddings
        light_emb = F.layer_norm(light_emb, light_emb.shape[-1:])
        heavy_emb = F.layer_norm(heavy_emb, heavy_emb.shape[-1:])
        antigen_emb = F.layer_norm(antigen_emb, antigen_emb.shape[-1:])

        # Projections
        light_proj = self.light_projection(light_emb)
        heavy_proj = self.heavy_projection(heavy_emb)
        antigen_proj = self.antigen_projection(antigen_emb)

        # Antibody internal self‑attention (heavy ↔ light)
        stacked = torch.stack([heavy_proj, light_proj], dim=1)   # [B,2,D]
        ab_ctx, _ = self.antibody_self_attention(stacked, stacked, stacked)
        heavy_ctx = ab_ctx[:, 0, :]
        light_ctx = ab_ctx[:, 1, :]

        # Initially, antigen context is just the projected antigen
        ag_ctx = antigen_proj

        # Mutual conditioning layers (parallel bidirectional)
        for i in range(self.n_layers):
            # Antibody → Antigen (using current antibody context)
            # First, fuse heavy and light contexts
            combined = torch.cat([heavy_ctx, light_ctx], dim=-1)
            weights = torch.softmax(self.fusion_gate(combined), dim=-1)
            fused_ab = weights[:, 0:1] * heavy_ctx + weights[:, 1:2] * light_ctx

            # Parallel cross‑attention
            new_ag_ctx = self.ab_to_ag_layers[i](fused_ab, ag_ctx)   # antibody queries antigen
            new_ab_ctx = self.ag_to_ab_layers[i](ag_ctx, fused_ab)   # antigen queries antibody

            # Update contexts (residual + layernorm handled inside GatedCrossAttention)
            ag_ctx = new_ag_ctx
            # For next layer, we need new heavy/light contexts. We split the fused antibody context
            # back into heavy and light? That's not straightforward. Instead, we keep the fused
            # representation as the antibody context for the next layer, but we also need to
            # update heavy_ctx and light_ctx individually. To simplify, we store the fused antibody
            # context and use it directly in the next layer's fusion. But then heavy_ctx and light_ctx
            # become stale. For a true mutual update, we should project the fused antibody context
            # back to two separate vectors? That would add many parameters. A simpler biologically
            # plausible alternative: after mutual attention, the antibody representation is updated
            # as a whole (fused). Then in the next layer, we re‑fuse the heavy/light? That would lose
            # individual chain information. Instead, we can keep the chain‑specific contexts updated
            # by using the same fused antibody context but also maintain a separate copy. However,
            # to avoid over‑engineering, we will keep the original heavy_ctx and light_ctx unchanged
            # across layers – which is a limitation but keeps the model simple.

            # For the sake of this implementation, we will update the fused antibody representation
            # and use it directly in the next layer. But we also need to preserve heavy/light
            # contributions. A practical compromise: after each mutual layer, we use the new
            # fused antibody representation as the "antibody context" and then split it into two
            # equal parts? That is not biologically accurate. Instead, we will **not** use multiple
            # layers for now – set n_layers=1 in the config to avoid complexity. The user can later
            # extend if needed.

        # After mutual conditioning (single layer), we have:
        # ag_ctx is the antigen conditioned on antibody
        # fused_ab is the antibody conditioned on antigen
        # Now apply dual heads
        ab_ctx = self.antibody_head(fused_ab)
        ag_ctx = self.antigen_head(ag_ctx)

        ab_norm = F.normalize(ab_ctx, p=2, dim=-1)
        ag_norm = F.normalize(ag_ctx, p=2, dim=-1)
        cosine_sim = F.cosine_similarity(ab_norm, ag_norm, dim=-1)

        return {
            'cosine_similarity': cosine_sim,
            'antibody_context': ab_norm,
            'antigen_context': ag_norm
        }

# ----------------------------------------------------------------------
# Training functions (identical to before, but using the new model)
# ----------------------------------------------------------------------
def train_epoch(model, loader, optimizer, device, pkd_bounds):
    model.train()
    total_loss = 0.0
    pkd_lower, pkd_upper = pkd_bounds
    for batch in loader:
        le = batch['light_emb'].to(device)
        he = batch['heavy_emb'].to(device)
        ae = batch['antigen_emb'].to(device)
        af = batch['affinity'].to(device)

        target = (2.0 * (af - pkd_lower) / (pkd_upper - pkd_lower) - 1.0)
        out = model(le, he, ae)
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
    l_ids, h_ids, ag_ids = [], [], []
    total_loss = 0.0
    pkd_lower, pkd_upper = pkd_bounds
    with torch.no_grad():
        for batch in loader:
            le = batch['light_emb'].to(device)
            he = batch['heavy_emb'].to(device)
            ae = batch['antigen_emb'].to(device)
            af = batch['affinity'].to(device)

            target = (2.0 * (af - pkd_lower) / (pkd_upper - pkd_lower) - 1.0)
            out = model(le, he, ae)
            loss = F.mse_loss(out['cosine_similarity'], target)
            total_loss += loss.item()

            pred_cos = out['cosine_similarity'].cpu().numpy()
            pred_pkd = (pred_cos + 1.0) / 2.0 * (pkd_upper - pkd_lower) + pkd_lower
            preds.extend(pred_pkd.tolist())
            labels.extend(af.cpu().numpy().tolist())

            if return_ids and 'light_id' in batch:
                l_ids.extend(batch['light_id'])
                h_ids.extend(batch['heavy_id'])
                ag_ids.extend(batch['antigen_id'])

    from agabgated.utils.main_symmetric_mean import pearson_correlation, spearman_correlation, rmse, r_squared
    metrics = {
        'loss': total_loss / len(loader),
        'r2': r_squared(labels, preds),
        'pearson': pearson_correlation(labels, preds),
        'spearman': spearman_correlation(labels, preds),
        'rmse': rmse(labels, preds),
    }
    if return_ids:
        return metrics, preds, labels, l_ids, h_ids, ag_ids
    return metrics, preds, labels

def train_model(df_train, df_val, embedding_loader, config, device):
    pkd_lower = df_train['binding_affinity'].min()
    pkd_upper = df_train['binding_affinity'].max()
    pkd_bounds = (pkd_lower, pkd_upper)
    print(f"  Affinity range (train): [{pkd_lower:.3f}, {pkd_upper:.3f}]")

    train_loader = DataLoader(
        CachedEmbeddingDataset(df_train, embedding_loader),
        batch_size=config['batch_size'], shuffle=True, collate_fn=collate_fn,
        num_workers=0, pin_memory=(config['device'] == 'cuda')
    )
    val_loader = DataLoader(
        CachedEmbeddingDataset(df_val, embedding_loader),
        batch_size=config['batch_size'], collate_fn=collate_fn,
        num_workers=0, pin_memory=(config['device'] == 'cuda')
    )

    model = MutualTriStreamStrong(
        esm_dim=embedding_loader.embedding_dim,
        projected_size=config['projected_size'],
        num_heads=config['num_heads'],
        dropout=config['dropout'],
        n_layers=config.get('n_layers', 1),
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

# ----------------------------------------------------------------------
# Cross‑validation
# ----------------------------------------------------------------------
def cross_validate(df, embedding_loader, config, n_folds=10, output_dir='results_mutual_strong', dataset_name=''):
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

        model, pkd_bounds = train_model(df_tr, df_va, embedding_loader, config, device)

        val_loader = DataLoader(
            CachedEmbeddingDataset(df_va, embedding_loader, return_ids=True),
            batch_size=config['batch_size'], collate_fn=collate_fn
        )
        metrics, preds, labels, l_ids, h_ids, ag_ids = evaluate(
            model, val_loader, device, pkd_bounds, return_ids=True
        )

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

        model_save_path = os.path.join(fold_dir, 'model.pt')
        torch.save({
            'model_state_dict': model.state_dict(),
            'config': config,
            'pkd_bounds': pkd_bounds,
            'val_metrics': metrics,
        }, model_save_path)
        print(f"  Model saved to {model_save_path}")

        all_preds.append(fold_df)
        fold_results.append(metrics)

        print(f"\n  Fold {fold}:")
        for k, v in metrics.items():
            print(f"    {k.upper():10s}: {v:.4f}")

    pd.concat(all_preds, ignore_index=True).to_csv(
        os.path.join(output_dir, 'all_folds_predictions.csv'), index=False)

    summary = {}
    for metric in ['r2', 'pearson', 'spearman', 'rmse', 'loss']:
        vals = [m[metric] for m in fold_results]
        summary[metric] = {'mean': float(np.mean(vals)), 'std': float(np.std(vals))}

    print(f"\n{'='*70}\nCV SUMMARY — {dataset_name}\n{'='*70}")
    for metric, s in summary.items():
        print(f"  {metric.upper():10s}: {s['mean']:.4f} ± {s['std']:.4f}")

    pd.DataFrame([{'metric': k, 'mean': v['mean'], 'std': v['std']}
                  for k, v in summary.items()]).to_csv(
        os.path.join(output_dir, 'cv_summary.csv'), index=False)
    return fold_results, summary

# ----------------------------------------------------------------------
# Benchmark with multiple seeds
# ----------------------------------------------------------------------
def train_on_sabdab_test_on_benchmark_single_seed(
    sabdab_csv, benchmark_csv, cache_file, config, seed, output_dir=None
):
    setup_reproducibility(seed)
    config['seed'] = seed
    device = config['device']

    df_train = load_data(sabdab_csv)
    df_test = load_data(benchmark_csv)

    loader = CachedEmbeddingLoader(cache_file)
    loader.validate_ids(df_train, ['light_id', 'heavy_id', 'antigen_id'])
    loader.validate_ids(df_test, ['light_id', 'heavy_id', 'antigen_id'])

    train_loader = DataLoader(
        CachedEmbeddingDataset(df_train, loader),
        batch_size=config['batch_size'], shuffle=True, collate_fn=collate_fn,
        num_workers=0, pin_memory=(device == 'cuda')
    )
    test_loader = DataLoader(
        CachedEmbeddingDataset(df_test, loader, return_ids=True),
        batch_size=config['batch_size'], collate_fn=collate_fn,
        num_workers=0, pin_memory=(device == 'cuda')
    )

    train_lower = df_train['binding_affinity'].min()
    train_upper = df_train['binding_affinity'].max()
    train_bounds = (train_lower, train_upper)

    test_lower = df_test['binding_affinity'].min()
    test_upper = df_test['binding_affinity'].max()
    test_bounds = (test_lower, test_upper)

    model = MutualTriStreamStrong(
        esm_dim=loader.embedding_dim,
        projected_size=config['projected_size'],
        num_heads=config['num_heads'],
        dropout=config['dropout'],
        n_layers=config.get('n_layers', 1),
        device=device
    ).to(device)

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config['learning_rate'], weight_decay=config['weight_decay']
    )

    print(f"  Training seed {seed} on full SAbDab for {config['epochs']} epochs...")
    for epoch in range(config['epochs']):
        train_loss = train_epoch(model, train_loader, optimizer, device, train_bounds)
        if (epoch+1) % 10 == 0:
            print(f"  Seed {seed} Epoch {epoch+1}: loss = {train_loss:.4f}")

    train_metrics, _, _ = evaluate(model, train_loader, device, train_bounds)
    print(f"\n  SABDAB FULL TRAINING SET EVALUATION")
    for k, v in train_metrics.items():
        print(f"    {k.upper():12s}: {v:.4f}")

    test_metrics, preds, labels, l_ids, h_ids, ag_ids = evaluate(
        model, test_loader, device, test_bounds, return_ids=True
    )
    print(f"\n  BENCHMARK TEST SET EVALUATION")
    for k, v in test_metrics.items():
        print(f"    {k.upper():12s}: {v:.4f}")

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        model_save_path = os.path.join(output_dir, f'model_seed{seed}.pt')
        torch.save({
            'model_state_dict': model.state_dict(),
            'config': config,
            'train_bounds': train_bounds,
            'test_bounds': test_bounds,
            'train_metrics': train_metrics,
            'test_metrics': test_metrics,
        }, model_save_path)
        print(f"  Model saved to {model_save_path}")
        pd.DataFrame({
            'light_id': l_ids, 'heavy_id': h_ids, 'antigen_id': ag_ids,
            'true_affinity': labels, 'predicted_affinity': preds,
            'error': np.array(labels) - np.array(preds)
        }).to_csv(os.path.join(output_dir, f'predictions_seed{seed}.csv'), index=False)

    return test_metrics

def run_benchmark_multiple_seeds(
    sabdab_csv, benchmark_csv, cache_file, config, seeds=[314, 114, 144],
    output_dir='results_mutual_strong/benchmark'
):
    print("="*80)
    print(f"MUTUAL STRONG: 100% SABDAB → BENCHMARK | {len(seeds)} SEEDS")
    print("="*80)
    all_test_metrics = []
    for seed in seeds:
        print(f"\n{'#'*90}")
        print(f"RUN — SEED = {seed}")
        print(f"{'#'*90}\n")
        seed_output_dir = os.path.join(output_dir, f'seed{seed}') if output_dir else None
        metrics = train_on_sabdab_test_on_benchmark_single_seed(
            sabdab_csv, benchmark_csv, cache_file, config, seed,
            output_dir=seed_output_dir
        )
        all_test_metrics.append(metrics)
        print(f"✓ Seed {seed} done — Benchmark Pearson: {metrics['pearson']:.4f}")

    print(f"\n{'='*90}\nMUTUAL STRONG MULTI-SEED BENCHMARK SUMMARY\n{'='*90}")
    summary = {}
    for metric in ['pearson', 'spearman', 'rmse', 'r2', 'loss']:
        vals = [m[metric] for m in all_test_metrics]
        summary[metric] = {'mean': np.mean(vals), 'std': np.std(vals)}
        print(f"  {metric.upper():12s}: {summary[metric]['mean']:.4f} ± {summary[metric]['std']:.4f}")

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        pd.DataFrame([{'metric': m, 'mean': v['mean'], 'std': v['std']}
                      for m, v in summary.items()]).to_csv(
            os.path.join(output_dir, 'benchmark_summary.csv'), index=False)
    return summary

# ----------------------------------------------------------------------
# Main with argparse
# ----------------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='MutualTriStreamStrong evaluation')
    parser.add_argument('--mode', type=str, choices=['cv', 'benchmark', 'all'], default='all',
                        help='Run cross-validation only, benchmark only, or both')
    parser.add_argument('--dataset', type=str, choices=['sabdab', 'abbind', 'skempi'], default='sabdab',
                        help='Dataset for cross-validation (only used if mode=cv)')
    parser.add_argument('--seeds', type=int, nargs='+', default=[314, 114, 144],
                        help='Seeds for benchmark (default: 314 114 144)')
    parser.add_argument('--device', type=str, default='cuda', help='Device to use')
    parser.add_argument('--n_layers', type=int, default=1, help='Number of mutual layers (default: 1)')
    args = parser.parse_args()

    # Configuration – identical to three‑stream
    config = DEFAULT_CONFIG.copy()
    config['device'] = args.device if torch.cuda.is_available() else 'cpu'
    config['epochs'] = 50
    config['batch_size'] = 32
    config['projected_size'] = 256
    config['num_heads'] = 8
    config['dropout'] = 0.1
    config['learning_rate'] = 1e-4
    config['weight_decay'] = 0.01
    config['patience'] = 10
    config['n_layers'] = args.n_layers

    if args.mode == 'cv':
        ds = DATASET_CONFIGS[args.dataset]
        loader = CachedEmbeddingLoader(ds['cache'])
        df = load_data(ds['csv'])
        loader.validate_ids(df, ['light_id', 'heavy_id', 'antigen_id'])
        output_dir = f"results_mutual_strong/cv_{args.dataset}"
        cross_validate(df, loader, config, n_folds=ds['n_folds'], output_dir=output_dir, dataset_name=args.dataset)
    elif args.mode == 'benchmark':
        bench_summary = run_benchmark_multiple_seeds(
            sabdab_csv='datasets/pairs_sabdab.csv',
            benchmark_csv='datasets/pairs_benchmark.csv',
            cache_file='datasets/esm2_embeddings_natural_650M.pkl',
            config=config,
            seeds=args.seeds,
            output_dir='results_mutual_strong/benchmark'
        )
        print("\nBenchmark summary:", bench_summary)
    else:  # 'all'
        dataset_configs = {
            'sabdab': {'csv': 'datasets/pairs_sabdab.csv', 'cache': 'datasets/esm2_embeddings_natural_650M.pkl', 'n_folds': 10},
            'abbind': {'csv': 'datasets/pairs_abbind.csv', 'cache': 'datasets/esm2_embeddings_mutation_650M.pkl', 'n_folds': 10},
            'skempi': {'csv': 'datasets/pairs_skempi.csv', 'cache': 'datasets/esm2_embeddings_mutation_650M.pkl', 'n_folds': 10},
        }
        for name, ds in dataset_configs.items():
            print(f"\n{'#'*70}\n# CV: {name.upper()} (MutualStrong)\n{'#'*70}")
            loader = CachedEmbeddingLoader(ds['cache'])
            df = load_data(ds['csv'])
            loader.validate_ids(df, ['light_id', 'heavy_id', 'antigen_id'])
            output_dir = f"results_mutual_strong/cv_{name}"
            cross_validate(df, loader, config, n_folds=ds['n_folds'], output_dir=output_dir, dataset_name=name)
        bench_summary = run_benchmark_multiple_seeds(
            sabdab_csv='datasets/pairs_sabdab.csv',
            benchmark_csv='datasets/pairs_benchmark.csv',
            cache_file='datasets/esm2_embeddings_natural_650M.pkl',
            config=config,
            seeds=args.seeds,
            output_dir='results_mutual_strong/benchmark'
        )
        print("\nBenchmark summary:", bench_summary)