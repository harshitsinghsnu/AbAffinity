#!/usr/bin/env python3
"""
Train MutualTriStreamStrong on SAaIntDB with 10-fold cross-validation.

Identical model architecture to mutual_strong.py.
Results go to results_mutual_strong_saaintdb/ (separate from other datasets).
predictions.csv includes ALL original CSV columns + predicted_affinity + error.

Dataset columns used
--------------------
  H_seq   → heavy chain sequence  → heavy_id  = <PDB_ID>_H_<H_chain_ID>
  L_seq   → light chain sequence  → light_id  = <PDB_ID>_L_<L_chain_ID>
  Ag_seq  → antigen sequence      → antigen_id = <PDB_ID>_Ag_<Ag_chain_ID(s)>
  pKD     → binding affinity

Usage
-----
  # Step 1 – precompute embeddings (run once)
  python precompute_embeddings_saaintdb.py

  # Step 2 – train with 10-fold CV
  python mutual_strong_saaintdb.py
  python mutual_strong_saaintdb.py --device cpu --n_layers 2
"""

import os
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader
from sklearn.model_selection import KFold

try:
    from agabgated.models.mutual_strong import MutualTriStreamStrong, train_epoch, evaluate
    from agabgated.utils.main_symmetric_mean import (
        setup_reproducibility,
        CachedEmbeddingLoader,
        CachedEmbeddingDataset,
        collate_fn,
        DEFAULT_CONFIG,
    )
except ImportError as e:
    raise ImportError(
        f"Could not import dependencies: {e}\n"
        "Make sure mutual_strong.py and main_symmetric_mean.py are in the same directory."
    )


# ---------------------------------------------------------------------------
# ID helpers  (must match precompute_embeddings_saaintdb.py)
# ---------------------------------------------------------------------------

def _make_heavy_id(pdb_id, h_chain_id):
    return f"{pdb_id}_H_{h_chain_id}"

def _make_light_id(pdb_id, l_chain_id):
    return f"{pdb_id}_L_{l_chain_id}"

def _make_antigen_id(pdb_id, ag_chain_ids):
    chains = str(ag_chain_ids).replace(',', '_').replace(' ', '')
    return f"{pdb_id}_Ag_{chains}"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _is_valid_seq(val) -> bool:
    """Return True only for non-null, non-placeholder sequence strings."""
    if pd.isna(val):
        return False
    s = str(val).strip()
    return s not in ('', 'N.A.', 'NA', 'N/A', 'nan', 'None')


def load_saaintdb(csv_path: str) -> pd.DataFrame:
    """
    Load SAaIntDB CSV and add the ID / binding_affinity columns needed by
    CachedEmbeddingDataset, while preserving every original column.

    H_seq and Ag_seq are complete for all rows (no filtering needed there).
    L_seq may be absent for nanobodies — those rows are kept and the heavy
    chain embedding is reused in the light chain slot.
    """
    print(f"  Loading {csv_path}")
    df = pd.read_csv(csv_path)
    print(f"  Columns : {list(df.columns)}")
    print(f"  Rows    : {len(df)}")

    # Heavy and antigen IDs are always valid
    df['heavy_id']   = df.apply(lambda r: _make_heavy_id(r['PDB_ID'], r['H_chain_ID']),       axis=1)
    df['antigen_id'] = df.apply(lambda r: _make_antigen_id(r['PDB_ID'], r['Ag_chain_ID(s)']), axis=1)

    # Light chain: reuse heavy embedding for nanobodies that have no light chain
    n_no_light = (~df['L_seq'].apply(_is_valid_seq)).sum()
    if n_no_light:
        print(f"  {n_no_light} nanobody rows have no light chain — "
              f"heavy chain embedding will be reused in the light chain slot.")

    def _light_id(row):
        if _is_valid_seq(row['L_seq']):
            return _make_light_id(row['PDB_ID'], row['L_chain_ID'])
        return _make_heavy_id(row['PDB_ID'], row['H_chain_ID'])   # fallback → heavy

    df['light_id']         = df.apply(_light_id, axis=1)
    df['binding_affinity'] = df['pKD']

    n_before = len(df)
    df = df.dropna(subset=['binding_affinity'])
    if len(df) < n_before:
        print(f"  Dropped {n_before - len(df)} rows with NaN pKD.")

    aff = df['binding_affinity']
    print(f"  Valid pairs : {len(df)} | affinity [{aff.min():.3f}, {aff.max():.3f}]")
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Cache validation (soft — drops unembedded rows instead of raising)
# ---------------------------------------------------------------------------

def filter_missing_embeddings(df: pd.DataFrame, loader: 'CachedEmbeddingLoader',
                              id_cols: list) -> pd.DataFrame:
    """
    Remove rows whose sequence IDs are absent from the embedding cache.
    Prints a warning for every missing ID instead of raising an error.
    """
    keep = pd.Series(True, index=df.index)
    for col in id_cols:
        missing_mask = ~df[col].apply(lambda sid: sid in loader.embeddings)
        n_miss = missing_mask.sum()
        if n_miss:
            missing_ids = df.loc[missing_mask, col].unique().tolist()
            print(f"  Warning: {n_miss} rows have no cached embedding for "
                  f"'{col}' — skipping. IDs: {missing_ids}")
        keep &= ~missing_mask

    n_dropped = (~keep).sum()
    if n_dropped:
        print(f"  Total rows skipped due to missing embeddings: {n_dropped} "
              f"(out of {len(df)})")
    return df[keep].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Split strategy
# ---------------------------------------------------------------------------

def get_fold_splits(df: pd.DataFrame, n_folds: int, seed: int,
                    split_type: str = 'random') -> list:
    """
    Return a list of (tr_idx, va_idx) integer-array pairs for n_folds folds.

    'random' — standard KFold on rows (rows from the same PDB may appear in
               both train and val).
    'cold'   — KFold is performed on unique PDB_IDs; every row belonging to a
               PDB_ID is kept together, so no PDB overlaps between train/val.
    """
    kfold = KFold(n_splits=n_folds, shuffle=True, random_state=seed)

    if split_type == 'cold':
        unique_pdbs = df['PDB_ID'].unique()
        splits = []
        for tr_pdb_idx, va_pdb_idx in kfold.split(unique_pdbs):
            tr_pdbs = set(unique_pdbs[tr_pdb_idx])
            va_pdbs = set(unique_pdbs[va_pdb_idx])
            tr_row_idx = np.where(df['PDB_ID'].isin(tr_pdbs))[0]
            va_row_idx = np.where(df['PDB_ID'].isin(va_pdbs))[0]
            splits.append((tr_row_idx, va_row_idx))
        return splits

    return list(kfold.split(df))   # random


# ---------------------------------------------------------------------------
# Training (mirrors mutual_strong.py but without re-importing train_model)
# ---------------------------------------------------------------------------

def _train_fold(df_train, df_val, embedding_loader, config, device):
    """Train one fold; return (model, pkd_bounds)."""
    pkd_lower = df_train['binding_affinity'].min()
    pkd_upper = df_train['binding_affinity'].max()
    pkd_bounds = (pkd_lower, pkd_upper)
    print(f"  Affinity range (train): [{pkd_lower:.3f}, {pkd_upper:.3f}]")

    pin = (device == 'cuda')
    train_loader = DataLoader(
        CachedEmbeddingDataset(df_train, embedding_loader),
        batch_size=config['batch_size'], shuffle=True,
        collate_fn=collate_fn, num_workers=0, pin_memory=pin,
    )
    val_loader = DataLoader(
        CachedEmbeddingDataset(df_val, embedding_loader),
        batch_size=config['batch_size'],
        collate_fn=collate_fn, num_workers=0, pin_memory=pin,
    )

    model = MutualTriStreamStrong(
        esm_dim        = embedding_loader.embedding_dim,
        projected_size = config['projected_size'],
        num_heads      = config['num_heads'],
        dropout        = config['dropout'],
        n_layers       = config.get('n_layers', 1),
        device         = device,
    ).to(device)

    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"  Trainable parameters: {n_params:,}")

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=config['learning_rate'],
        weight_decay=config['weight_decay'],
    )

    best_pearson, best_state, patience = -float('inf'), None, 0

    for epoch in range(config['epochs']):
        train_loss = train_epoch(model, train_loader, optimizer, device, pkd_bounds)
        val_metrics, _, _ = evaluate(model, val_loader, device, pkd_bounds)

        print(f"  Epoch {epoch+1:3d}/{config['epochs']} | "
              f"Loss: {train_loss:.4f} | "
              f"Val Pearson: {val_metrics['pearson']:.4f} | "
              f"Val RMSE:    {val_metrics['rmse']:.4f}")

        if val_metrics['pearson'] > best_pearson:
            best_pearson = val_metrics['pearson']
            best_state   = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            patience     = 0
        else:
            patience += 1
            if patience >= config['patience']:
                print(f"  Early stopping at epoch {epoch+1} "
                      f"(best val Pearson: {best_pearson:.4f})")
                break

    if best_state:
        model.load_state_dict({k: v.to(device) for k, v in best_state.items()})
    return model, pkd_bounds


# ---------------------------------------------------------------------------
# Cross-validation
# ---------------------------------------------------------------------------

def cross_validate_saaintdb(
    df_full:          pd.DataFrame,
    embedding_loader: 'CachedEmbeddingLoader',
    config:           dict,
    n_folds:          int = 10,
    output_dir:       str = 'results_mutual_strong_saaintdb',
    split_type:       str = 'random',
):
    """
    10-fold CV on SAaIntDB.
    split_type='random' : standard row-level KFold.
    split_type='cold'   : KFold on unique PDB_IDs; no PDB overlap between folds.
    predictions.csv per fold contains ALL original CSV columns
    plus predicted_affinity and error.
    """
    splits = get_fold_splits(df_full, n_folds, config['seed'], split_type)
    device = config['device']
    os.makedirs(output_dir, exist_ok=True)

    fold_results = []
    all_fold_dfs = []

    print(f"\n{'='*70}")
    print(f"MutualTriStreamStrong — SAaIntDB — {n_folds}-fold CV ({split_type} split)")
    if split_type == 'cold':
        print(f"  Unique PDB_IDs : {df_full['PDB_ID'].nunique()}")
    print(f"  Samples : {len(df_full)}")
    print(f"  Device  : {device}")
    print(f"  Output  : {output_dir}")
    print(f"{'='*70}")

    for fold, (tr_idx, va_idx) in enumerate(splits, 1):
        print(f"\n{'='*70}\nFOLD {fold}/{n_folds}\n{'='*70}")

        df_tr = df_full.iloc[tr_idx].reset_index(drop=True)
        df_va = df_full.iloc[va_idx].reset_index(drop=True)
        print(f"  Train : {len(df_tr)} | Val : {len(df_va)}")

        model, pkd_bounds = _train_fold(df_tr, df_va, embedding_loader, config, device)

        # Evaluate with IDs returned (val loader must NOT shuffle)
        val_loader_ids = DataLoader(
            CachedEmbeddingDataset(df_va, embedding_loader, return_ids=True),
            batch_size=config['batch_size'],
            collate_fn=collate_fn,
            shuffle=False,
        )
        metrics, preds, labels, l_ids, h_ids, ag_ids = evaluate(
            model, val_loader_ids, device, pkd_bounds, return_ids=True
        )

        # ---------------------------------------------------------------
        # Build predictions DataFrame with ALL original columns
        # ---------------------------------------------------------------
        fold_df = df_va.copy().reset_index(drop=True)
        fold_df['predicted_affinity'] = preds
        fold_df['error']              = np.array(labels) - np.array(preds)

        # ---------------------------------------------------------------
        # Save fold outputs
        # ---------------------------------------------------------------
        fold_dir = os.path.join(output_dir, f'fold_{fold:02d}')
        os.makedirs(fold_dir, exist_ok=True)

        fold_df.to_csv(os.path.join(fold_dir, 'predictions.csv'), index=False)
        pd.DataFrame([metrics]).to_csv(os.path.join(fold_dir, 'metrics.csv'), index=False)
        pd.DataFrame([{'pkd_lower': pkd_bounds[0], 'pkd_upper': pkd_bounds[1]}]).to_csv(
            os.path.join(fold_dir, 'pkd_bounds.csv'), index=False)

        model_path = os.path.join(fold_dir, 'model.pt')
        torch.save({
            'model_state_dict': model.state_dict(),
            'config':           config,
            'pkd_bounds':       pkd_bounds,
            'val_metrics':      metrics,
        }, model_path)
        print(f"  Saved model → {model_path}")

        all_fold_dfs.append(fold_df)
        fold_results.append(metrics)

        print(f"\n  Fold {fold} results:")
        for k, v in metrics.items():
            print(f"    {k.upper():10s}: {v:.4f}")

    # -------------------------------------------------------------------
    # Aggregate across folds
    # -------------------------------------------------------------------
    pd.concat(all_fold_dfs, ignore_index=True).to_csv(
        os.path.join(output_dir, 'all_folds_predictions.csv'), index=False)

    summary = {}
    for metric in ['r2', 'pearson', 'spearman', 'rmse', 'loss']:
        vals = [m[metric] for m in fold_results]
        summary[metric] = {'mean': float(np.mean(vals)), 'std': float(np.std(vals))}

    print(f"\n{'='*70}\nCV SUMMARY — SAaIntDB (MutualTriStreamStrong)\n{'='*70}")
    for metric, s in summary.items():
        print(f"  {metric.upper():10s}: {s['mean']:.4f} ± {s['std']:.4f}")

    pd.DataFrame([
        {'metric': k, 'mean': v['mean'], 'std': v['std']}
        for k, v in summary.items()
    ]).to_csv(os.path.join(output_dir, 'cv_summary.csv'), index=False)

    return fold_results, summary


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='MutualTriStreamStrong — SAaIntDB 10-fold CV')
    parser.add_argument('--csv',        default='data/saaintdb_with_antigen_names.csv',
                        help='Path to SAaIntDB CSV')
    parser.add_argument('--cache',      default='data/esm2_embeddings_saaintdb_650M.pkl',
                        help='Pre-computed ESM-2 embedding cache')
    parser.add_argument('--output_dir', default=None,
                        help='Output directory (default: results_mutual_strong_saaintdb[_cold])')
    parser.add_argument('--split_type', choices=['random', 'cold'], default='random',
                        help='random: row-level KFold | cold: PDB-ID-level, no PDB overlap')
    parser.add_argument('--device',     default='cuda')
    parser.add_argument('--n_folds',    type=int, default=10)
    parser.add_argument('--n_layers',   type=int, default=2,
                        help='Number of mutual interaction layers (default: 1)')
    parser.add_argument('--epochs',         type=int,   default=50)
    parser.add_argument('--batch_size',     type=int,   default=32)
    parser.add_argument('--projected_size', type=int,   default=256)
    parser.add_argument('--num_heads',      type=int,   default=8)
    parser.add_argument('--dropout',        type=float, default=0.1)
    parser.add_argument('--lr',             type=float, default=1e-4)
    parser.add_argument('--weight_decay',   type=float, default=0.01)
    parser.add_argument('--patience',       type=int,   default=10)
    parser.add_argument('--seed',           type=int,   default=9999)
    args = parser.parse_args()

    config = DEFAULT_CONFIG.copy()
    config.update({
        'device':         args.device if torch.cuda.is_available() else 'cpu',
        'seed':           args.seed,
        'epochs':         args.epochs,
        'batch_size':     args.batch_size,
        'projected_size': args.projected_size,
        'num_heads':      args.num_heads,
        'dropout':        args.dropout,
        'learning_rate':  args.lr,
        'weight_decay':   args.weight_decay,
        'patience':       args.patience,
        'n_layers':       args.n_layers,
    })

    setup_reproducibility(config['seed'])

    print(f"\nConfig:")
    for k, v in config.items():
        print(f"  {k:20s}: {v}")

    # Load and prepare dataset
    df = load_saaintdb(args.csv)

    # Load embedding cache and drop any rows whose IDs are not in the cache
    loader = CachedEmbeddingLoader(args.cache)
    df = filter_missing_embeddings(df, loader, ['light_id', 'heavy_id', 'antigen_id'])

    # Resolve output directory (auto-suffix _cold when using cold split)
    output_dir = args.output_dir or (
        'results_mutual_strong_saaintdb_cold' if args.split_type == 'cold'
        else 'results_mutual_strong_saaintdb'
    )

    # Run cross-validation
    cross_validate_saaintdb(
        df_full          = df,
        embedding_loader = loader,
        config           = config,
        n_folds          = args.n_folds,
        output_dir       = output_dir,
        split_type       = args.split_type,
    )
