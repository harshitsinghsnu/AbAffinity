"""
run_multiseed.py
================
YAML-driven multi-seed experiment runner for the MutualTriStreamStrong (Ours)
model on the natural antibody-antigen datasets.

Usage
-----
  python experiments/run_multiseed.py --config experiments/configs/exp01_ours_meanpool_cv.yaml
  python experiments/run_multiseed.py --all          # run every exp*.yaml
  python experiments/run_multiseed.py --all --seeds 42 114 144   # override seeds

Outputs (organized)
-------------------
  experiments/results/<name>/
      seed_<s>/<dataset>_cv.csv        per-seed, per-dataset fold metrics  (cv mode)
      seed_<s>/benchmark.csv           per-seed benchmark metrics          (benchmark mode)
      aggregated_summary.csv           mean ± std across seeds (and folds)
      config_resolved.yaml             the fully-resolved config used
"""
import os, sys, glob, json, pickle, argparse, copy
import numpy as np, pandas as pd, torch
import yaml
from torch.utils.data import DataLoader
from scipy.stats import pearsonr, spearmanr

HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # 3_stream/
sys.path.insert(0, HERE)
from AbAffinity.models.mutual_strong import MutualTriStreamStrong, train_epoch, evaluate
from AbAffinity.models.mutual_strong_saaintdb import _train_fold, get_fold_splits
from AbAffinity.utils.main_symmetric_mean import (load_data, CachedEmbeddingDataset, collate_fn,
                                 DEFAULT_CONFIG, setup_reproducibility)

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
CFG_DIR = os.path.join(HERE, 'configs', 'benchmark')
RES_DIR = os.path.join(HERE, 'results')

# ── config loading ───────────────────────────────────────────────────────────
def load_base():
    with open(os.path.join(CFG_DIR, 'base.yaml')) as f:
        return yaml.safe_load(f)

def resolve(cfg_path, base, seeds_override=None):
    with open(cfg_path) as f:
        exp = yaml.safe_load(f)
    hp = dict(base['defaults']['hyperparams'])
    hp.update(exp.get('hyperparams', {}))
    seeds = seeds_override or exp.get('seeds', base['defaults']['seeds'])
    return {
        'name': exp['name'], 'description': exp.get('description', ''),
        'mode': exp['mode'], 'pooling': exp['pooling'],
        'ablation': exp.get('ablation', None),
        'datasets': exp.get('datasets', []),
        'train_dataset': exp.get('train_dataset'),
        'test_dataset':  exp.get('test_dataset'),
        'n_folds': exp.get('n_folds', 10),
        'seeds': seeds, 'hp': hp,
        'ds_registry': base['datasets'], 'emb_registry': base['embeddings'],
    }

# ── embedding loaders (cached per family/pooling) ────────────────────────────
_EMB_CACHE = {}
def get_loader(pooling, family, emb_registry):
    key = (pooling, family)
    if key in _EMB_CACHE:
        return _EMB_CACHE[key]
    path = os.path.join(HERE, emb_registry[pooling][family])
    with open(path, 'rb') as f:
        d = pickle.load(f)
    class _L:
        def __init__(self, dd): self.embeddings = dd; self.embedding_dim = next(iter(dd.values())).shape[0]
        def get_embedding(self, k): return self.embeddings[k]
    loader = _L(d)
    _EMB_CACHE[key] = loader
    print(f"    loaded embeddings [{pooling}/{family}] dim={loader.embedding_dim} n={len(d)}")
    return loader

def make_cfg(hp, seed):
    return {**DEFAULT_CONFIG, **hp, 'device': DEVICE, 'seed': seed}

# ── antigen-shuffle wrapper ──────────────────────────────────────────────────
def shuffled_loader(loader, df, seed):
    """Return a loader whose antigen_ids are randomly permuted across rows."""
    rng = np.random.default_rng(seed)
    ag_ids = df['antigen_id'].unique().tolist()
    perm = rng.permutation(ag_ids)
    remap = dict(zip(ag_ids, perm))
    df2 = df.copy()
    df2['antigen_id'] = df2['antigen_id'].map(remap)
    return df2

# ── CV runner ─────────────────────────────────────────────────────────────────
def run_cv_one(df, loader, hp, seed, n_folds):
    cfg = make_cfg(hp, seed)
    setup_reproducibility(seed)
    splits = get_fold_splits(df, n_folds, seed, 'random')
    idc = ['heavy_id', 'light_id', 'antigen_id']
    rows = []
    for fold, (tr, va) in enumerate(splits, 1):
        d_tr = df.iloc[tr].reset_index(drop=True); d_va = df.iloc[va].reset_index(drop=True)
        d_tr = d_tr[d_tr[idc].apply(lambda c: c.isin(loader.embeddings)).all(axis=1)].reset_index(drop=True)
        d_va = d_va[d_va[idc].apply(lambda c: c.isin(loader.embeddings)).all(axis=1)].reset_index(drop=True)
        if len(d_tr) < 5 or len(d_va) < 2: continue
        model, bounds = _train_fold(d_tr, d_va, loader, cfg, DEVICE)
        lo, hi = bounds
        vl = DataLoader(CachedEmbeddingDataset(d_va, loader), batch_size=cfg['batch_size'], collate_fn=collate_fn)
        p, t = [], []
        model.eval()
        with torch.no_grad():
            for b in vl:
                cos = model(b['light_emb'].to(DEVICE), b['heavy_emb'].to(DEVICE),
                            b['antigen_emb'].to(DEVICE))['cosine_similarity'].cpu().numpy()
                p.extend(((cos+1)/2*(hi-lo)+lo).tolist()); t.extend(b['affinity'].tolist())
        p, t = np.array(p), np.array(t)
        rows.append({'fold': fold, 'pearson': pearsonr(t,p)[0],
                     'spearman': spearmanr(t,p)[0], 'rmse': float(np.sqrt(np.mean((t-p)**2)))})
    return pd.DataFrame(rows)

# ── Benchmark runner (full epochs, no early stop) ────────────────────────────
def run_benchmark_one(df_tr, df_b, loader, hp, seed):
    cfg = make_cfg(hp, seed)
    setup_reproducibility(seed)
    idc = ['heavy_id','light_id','antigen_id']
    df_tr = df_tr[df_tr[idc].apply(lambda c: c.isin(loader.embeddings)).all(axis=1)].reset_index(drop=True)
    df_b  = df_b[df_b[idc].apply(lambda c: c.isin(loader.embeddings)).all(axis=1)].reset_index(drop=True)
    sb = (df_tr['binding_affinity'].min(), df_tr['binding_affinity'].max())
    bb = (df_b['binding_affinity'].min(),  df_b['binding_affinity'].max())
    tl = DataLoader(CachedEmbeddingDataset(df_tr, loader), batch_size=cfg['batch_size'],
                    shuffle=True, collate_fn=collate_fn)
    bl = DataLoader(CachedEmbeddingDataset(df_b, loader), batch_size=cfg['batch_size'], collate_fn=collate_fn)
    model = MutualTriStreamStrong(esm_dim=loader.embedding_dim, projected_size=cfg['projected_size'],
            num_heads=cfg['num_heads'], dropout=cfg['dropout'], n_layers=cfg['n_layers'], device=DEVICE).to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=cfg['learning_rate'], weight_decay=cfg['weight_decay'])
    for _ in range(cfg['epochs']):
        train_epoch(model, tl, opt, DEVICE, sb)
    met, preds, labels = evaluate(model, bl, DEVICE, bb)
    return {'pearson': met['pearson'], 'spearman': met['spearman'], 'rmse': met['rmse'], 'n': len(df_b)}, labels, preds

# ── per-experiment driver ─────────────────────────────────────────────────────
def run_experiment(cfg):
    name = cfg['name']
    out = os.path.join(RES_DIR, name); os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'config_resolved.yaml'), 'w') as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    print(f"\n{'='*70}\nEXPERIMENT: {name}\n  {cfg['description']}\n  seeds={cfg['seeds']}\n{'='*70}")

    agg_rows = []
    if cfg['mode'] == 'cv':
        for ds in cfg['datasets']:
            reg = cfg['ds_registry'][ds]
            loader = get_loader(cfg['pooling'], reg['family'], cfg['emb_registry'])
            df = load_data(os.path.join(HERE, reg['pairs']))
            per_seed_means = []
            for s in cfg['seeds']:
                df_use = shuffled_loader(loader, df, s) if cfg['ablation'] == 'shuffle_antigen' else df
                fold_df = run_cv_one(df_use, loader, cfg['hp'], s, cfg['n_folds'])
                sd = os.path.join(out, f'seed_{s}'); os.makedirs(sd, exist_ok=True)
                fold_df.to_csv(os.path.join(sd, f'{ds}_cv.csv'), index=False)
                m = {k: fold_df[k].mean() for k in ['pearson','spearman','rmse']}
                per_seed_means.append(m)
                print(f"  [{ds}] seed {s}: r={m['pearson']:.4f} rho={m['spearman']:.4f}")
            pm = pd.DataFrame(per_seed_means)
            agg_rows.append({'dataset': ds,
                'pearson_mean': pm.pearson.mean(), 'pearson_std': pm.pearson.std(),
                'spearman_mean': pm.spearman.mean(), 'spearman_std': pm.spearman.std(),
                'rmse_mean': pm.rmse.mean(), 'rmse_std': pm.rmse.std(),
                'n_seeds': len(cfg['seeds'])})
    elif cfg['mode'] == 'benchmark':
        reg_tr = cfg['ds_registry'][cfg['train_dataset']]
        reg_b  = cfg['ds_registry'][cfg['test_dataset']]
        loader = get_loader(cfg['pooling'], reg_tr['family'], cfg['emb_registry'])
        df_tr = load_data(os.path.join(HERE, reg_tr['pairs']))
        df_b  = load_data(os.path.join(HERE, reg_b['pairs']))
        per_seed = []
        for s in cfg['seeds']:
            res, labels, preds = run_benchmark_one(df_tr, df_b, loader, cfg['hp'], s)
            sd = os.path.join(out, f'seed_{s}'); os.makedirs(sd, exist_ok=True)
            pd.DataFrame([res]).to_csv(os.path.join(sd, 'benchmark.csv'), index=False)
            pd.DataFrame({'true': labels, 'pred': preds}).to_csv(os.path.join(sd, 'benchmark_preds.csv'), index=False)
            per_seed.append(res)
            print(f"  benchmark seed {s}: r={res['pearson']:.4f} rho={res['spearman']:.4f}")
        pm = pd.DataFrame(per_seed)
        agg_rows.append({'dataset': cfg['test_dataset'],
            'pearson_mean': pm.pearson.mean(), 'pearson_std': pm.pearson.std(),
            'spearman_mean': pm.spearman.mean(), 'spearman_std': pm.spearman.std(),
            'rmse_mean': pm.rmse.mean(), 'rmse_std': pm.rmse.std(), 'n_seeds': len(cfg['seeds'])})

    agg = pd.DataFrame(agg_rows)
    agg.insert(0, 'experiment', name)
    agg.to_csv(os.path.join(out, 'aggregated_summary.csv'), index=False)
    print(f"\n  AGGREGATED ({name}):")
    print(agg[['dataset','pearson_mean','pearson_std','spearman_mean','rmse_mean']].to_string(index=False))
    return agg

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config'); ap.add_argument('--all', action='store_true')
    ap.add_argument('--seeds', type=int, nargs='+', default=None)
    args = ap.parse_args()
    base = load_base()
    cfgs = (sorted(glob.glob(os.path.join(CFG_DIR, 'exp*.yaml')))
            if args.all else [args.config])
    all_agg = []
    for cp in cfgs:
        cfg = resolve(cp, base, args.seeds)
        all_agg.append(run_experiment(cfg))
    # master summary across all experiments
    master = pd.concat(all_agg, ignore_index=True)
    master.to_csv(os.path.join(RES_DIR, 'MASTER_SUMMARY.csv'), index=False)
    print(f"\n{'#'*70}\nMASTER SUMMARY -> experiments/results/MASTER_SUMMARY.csv\n{'#'*70}")
    print(master.to_string(index=False))

if __name__ == '__main__':
    main()
