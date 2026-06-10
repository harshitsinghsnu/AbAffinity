#!/usr/bin/env python3
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
from scipy.stats import pearsonr, spearmanr
from copy import deepcopy

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

def get_metrics(y_true, y_pred):
    y_true, y_pred = np.array(y_true), np.array(y_pred)
    # Basic protection against constant predictions
    if np.std(y_true) < 1e-6 or np.std(y_pred) < 1e-6:
        p_corr, s_corr = 0.0, 0.0
    else:
        p_corr, _ = pearsonr(y_true, y_pred)
        s_corr, _ = spearmanr(y_true, y_pred)
    
    mse = mean_squared_error(y_true, y_pred)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot > 1e-6 else 0.0
    
    return {
        'pearson': float(p_corr),
        'spearman': float(s_corr),
        'rmse': float(np.sqrt(mse)),
        'r2': float(r2)
    }

# ==============================================================================
# DATA LOADING
# ==============================================================================

def load_data(pairs_csv):
    df = pd.read_csv(pairs_csv)
    col_map = {}
    for c in df.columns:
        low = c.lower().strip()
        if 'light' in low: col_map[c] = 'light_id'
        elif 'heavy' in low: col_map[c] = 'heavy_id'
        elif 'antigen' in low: col_map[c] = 'antigen_id'
        elif low in ('y', 'affinity', 'pkd', 'binding_affinity'): col_map[c] = 'affinity'
    df = df.rename(columns=col_map)
    return df[['light_id', 'heavy_id', 'antigen_id', 'affinity']].dropna()

class EmbeddingDataset(Dataset):
    def __init__(self, df, cache_file):
        self.df = df.reset_index(drop=True)
        with open(cache_file, 'rb') as f:
            self.embeddings = pickle.load(f)
        self.dim = next(iter(self.embeddings.values())).shape[0]

    def __len__(self): return len(self.df)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        # Direct concatenation and conversion to float32
        x = np.concatenate([
            self.embeddings[row['light_id']],
            self.embeddings[row['heavy_id']],
            self.embeddings[row['antigen_id']]
        ]).astype(np.float32)
        
        return {
            'x': torch.from_numpy(x),
            'y': torch.tensor(row['affinity'], dtype=torch.float32)
        }

# ==============================================================================
# BASELINE MODEL: CONCAT + MLP
# ==============================================================================

class ConcatMLP(nn.Module):
    def __init__(self, input_dim, hidden_dims=[1280, 256], dropout=0.1):
        super().__init__()
        layers = []
        curr_dim = input_dim
        for h in hidden_dims:
            layers.append(nn.Linear(curr_dim, h))
            layers.append(nn.BatchNorm1d(h))
            layers.append(nn.ReLU()) # Standard ReLU for MLP
            layers.append(nn.Dropout(dropout))
            curr_dim = h
        layers.append(nn.Linear(curr_dim, 1)) # No activation at output (Linear Regression)
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x).squeeze(-1)

# ==============================================================================
# ENGINE
# ==============================================================================

def run_train_eval(model, train_loader, val_loader, config):
    device = config['device']
    optimizer = torch.optim.AdamW(model.parameters(), lr=config['lr'], weight_decay=config['wd'])
    
    best_pearson = -1.0
    best_weights = None
    patience_counter = 0

    for epoch in range(config['epochs']):
        # Train
        model.train()
        for batch in train_loader:
            x, y = batch['x'].to(device), batch['y'].to(device)
            optimizer.zero_grad()
            loss = F.mse_loss(model(x), y)
            loss.backward()
            optimizer.step()

        # Val
        model.eval()
        preds, truths = [], []
        with torch.no_grad():
            for batch in val_loader:
                x, y = batch['x'].to(device), batch['y'].to(device)
                preds.extend(model(x).cpu().numpy())
                truths.extend(y.cpu().numpy())
        
        m = get_metrics(truths, preds)
        if m['pearson'] > best_pearson:
            best_pearson = m['pearson']
            best_weights = deepcopy(model.state_dict())
            patience_counter = 0
        else:
            patience_counter += 1
        
        if patience_counter >= config['patience']: break

    if best_weights:
        model.load_state_dict(best_weights)
    return model

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

def main():
    config = {
        'device': 'cuda' if torch.cuda.is_available() else 'cpu',
        'lr': 1e-4,
        'wd': 0.01,
        'epochs': 50,
        'batch_size': 32,
        'patience': 10,
        'hidden_dims': [1280, 256],
        'dropout': 0.1
    }

    datasets = {
        'sabdab':  ('datasets/pairs_sabdab.csv',  'datasets/esm2_embeddings_natural_650M.pkl'),
        'abbind':  ('datasets/pairs_abbind.csv',  'datasets/esm2_embeddings_mutation_650M.pkl'),
        'skempi':  ('datasets/pairs_skempi.csv',  'datasets/esm2_embeddings_mutation_650M.pkl'),
        'bench':   ('datasets/pairs_benchmark.csv','datasets/esm2_embeddings_natural_650M.pkl')
    }

    # 1. 10-Fold CV
    for name in ['sabdab', 'abbind', 'skempi']:
        if not os.path.exists(datasets[name][0]): continue
        print(f"\n>>> 10-Fold CV on {name.upper()}")
        df = load_data(datasets[name][0])
        ds = EmbeddingDataset(df, datasets[name][1])
        kf = KFold(n_splits=10, shuffle=True, random_state=42)
        
        cv_results = []
        for fold, (t_idx, v_idx) in enumerate(kf.split(df)):
            setup_reproducibility(42 + fold)
            t_loader = DataLoader(torch.utils.data.Subset(ds, t_idx), batch_size=config['batch_size'], shuffle=True)
            v_loader = DataLoader(torch.utils.data.Subset(ds, v_idx), batch_size=config['batch_size'])
            
            model = ConcatMLP(ds.dim * 3, config['hidden_dims'], config['dropout']).to(config['device'])
            model = run_train_eval(model, t_loader, v_loader, config)
            
            # Eval Fold
            preds, truths = [], []
            model.eval()
            with torch.no_grad():
                for batch in v_loader:
                    preds.extend(model(batch['x'].to(config['device'])).cpu().numpy())
                    truths.extend(batch['y'].numpy())
            
            fold_m = get_metrics(truths, preds)
            cv_results.append(fold_m)
            print(f" Fold {fold+1:02d} | Pearson: {fold_m['pearson']:.4f}")

        # Dataset Summary
        print(f"\n--- {name.upper()} FINAL SUMMARY ---")
        for m_name in ['pearson', 'spearman', 'rmse', 'r2']:
            vals = [r[m_name] for r in cv_results]
            print(f" {m_name.upper():10s}: {np.mean(vals):.4f} ± {np.std(vals):.4f}")

    # 2. Benchmark (SAbDab -> Benchmark)
    if os.path.exists(datasets['bench'][0]):
        print("\n" + "="*40 + "\n>>> BENCHMARK (SAbDab -> Benchmark)\n" + "="*40)
        df_tr = load_data(datasets['sabdab'][0])
        ds_tr = EmbeddingDataset(df_tr, datasets['sabdab'][1])
        df_te = load_data(datasets['bench'][0])
        ds_te = EmbeddingDataset(df_te, datasets['bench'][1])
        
        bench_scores = []
        for seed in [314, 114, 144]:
            setup_reproducibility(seed)
            t_ldr = DataLoader(ds_tr, batch_size=config['batch_size'], shuffle=True)
            v_ldr = DataLoader(ds_te, batch_size=config['batch_size'])
            
            model = ConcatMLP(ds_tr.dim * 3, config['hidden_dims'], config['dropout']).to(config['device'])
            model = run_train_eval(model, t_ldr, v_ldr, config)
            
            preds, truths = [], []
            with torch.no_grad():
                for batch in v_ldr:
                    preds.extend(model(batch['x'].to(config['device'])).cpu().numpy())
                    truths.extend(batch['y'].numpy())
            
            m = get_metrics(truths, preds)
            bench_scores.append(m)
            print(f" Seed {seed} | Pearson: {m['pearson']:.4f} | RMSE: {m['rmse']:.4f}")

        print(f"\n--- BENCHMARK FINAL SUMMARY (3 SEEDS) ---")
        for m_name in ['pearson', 'spearman', 'rmse', 'r2']:
            vals = [r[m_name] for r in bench_scores]
            print(f" {m_name.upper():10s}: {np.mean(vals):.4f} ± {np.std(vals):.4f}")

if __name__ == "__main__":
    main()