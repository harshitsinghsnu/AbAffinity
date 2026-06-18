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
        'pairs_csv':  'data/pairs_sabdab.csv',
        'cache_file': 'data/esm2_embeddings_natural_650M.pkl',
        'n_folds':    10,
    },
    'abbind': {
        'pairs_csv':  'data/pairs_abbind.csv',
        'cache_file': 'data/esm2_embeddings_mutation_650M.pkl',
        'n_folds':    10,
    },
    'skempi': {
        'pairs_csv':  'data/pairs_skempi.csv',
        'cache_file': 'data/esm2_embeddings_mutation_650M.pkl',
        'n_folds':    10,
    },
}


