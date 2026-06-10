"""
run_saaintdb_multiseed.py
=========================
YAML-driven multi-seed re-runs of the SAaIntDB experiments (MutualTriStreamStrong),
mirroring run_multiseed.py for the natural datasets.

Experiments (configs in experiments/configs_saaintdb/):
  - Ours mean-pool ESM-2  | random + cold split | 10-fold CV
  - Ours All-CDR (final)  | random + cold split | 10-fold CV

Usage
-----
  python experiments/run_saaintdb_multiseed.py --all --seeds 42 114 144
  python experiments/run_saaintdb_multiseed.py --config experiments/configs_saaintdb/sa_ours_meanpool_random.yaml

Outputs (organized, reproducible):
  experiments/results_saaintdb/<name>/
      seed_<s>/cv_summary.csv
      aggregated_summary.csv      (mean ± std over seeds)
      config_resolved.yaml
  experiments/results_saaintdb/MASTER_SUMMARY_SAAINTDB.csv
"""
import os, sys, glob, pickle, argparse
import numpy as np, pandas as pd, torch, yaml
from torch.utils.data import DataLoader
from scipy.stats import pearsonr, spearmanr

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from agabgated.models.mutual_strong_saaintdb import load_saaintdb, get_fold_splits, _train_fold
from agabgated.utils.main_symmetric_mean import CachedEmbeddingDataset, collate_fn, DEFAULT_CONFIG, setup_reproducibility

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
EXP = os.path.join(HERE, 'experiments')
RES = os.path.join(EXP, 'results_saaintdb'); os.makedirs(RES, exist_ok=True)
CSV = os.path.join(HERE, 'datasets/saaintdb_with_antigen_names.csv')

class DictLoader:
    def __init__(s, d): s.embeddings=d; s.embedding_dim=next(iter(d.values())).shape[0]
    def get_embedding(s, k): return s.embeddings[k]

_CACHE = {}
def get_loader(pooling):
    if pooling in _CACHE: return _CACHE[pooling]
    with open(os.path.join(HERE,'datasets/esm2_embeddings_saaintdb_650M.pkl'),'rb') as f:
        mean = pickle.load(f)
    if pooling == 'meanpool':
        loader = DictLoader(dict(mean))
    else:  # allcdr: overwrite heavy ids with CDR-pooled heavy
        with open(os.path.join(HERE,'results_saaintdb_allcdr/saaintdb_heavy_cdr_embeddings.pkl'),'rb') as f:
            cdr = pickle.load(f)
        comb = dict(mean); comb.update(cdr)
        loader = DictLoader(comb)
    _CACHE[pooling] = loader
    print(f"  loaded embeddings [{pooling}] dim={loader.embedding_dim} n={len(loader.embeddings)}")
    return loader

def run_cv_one(df, loader, cfg, n_folds, split_type):
    splits = get_fold_splits(df, n_folds, cfg['seed'], split_type)
    idc = ['heavy_id','light_id','antigen_id']
    rows = []
    for fold,(tr,va) in enumerate(splits,1):
        d_tr=df.iloc[tr].reset_index(drop=True); d_va=df.iloc[va].reset_index(drop=True)
        d_tr=d_tr[d_tr[idc].apply(lambda c:c.isin(loader.embeddings)).all(axis=1)].reset_index(drop=True)
        d_va=d_va[d_va[idc].apply(lambda c:c.isin(loader.embeddings)).all(axis=1)].reset_index(drop=True)
        if len(d_tr)<5 or len(d_va)<2: continue
        model,bounds=_train_fold(d_tr,d_va,loader,cfg,DEVICE); lo,hi=bounds
        vl=DataLoader(CachedEmbeddingDataset(d_va,loader),batch_size=cfg['batch_size'],collate_fn=collate_fn)
        P,T=[],[]; model.eval()
        with torch.no_grad():
            for b in vl:
                cos=model(b['light_emb'].to(DEVICE),b['heavy_emb'].to(DEVICE),b['antigen_emb'].to(DEVICE))['cosine_similarity'].cpu().numpy()
                P.extend(((cos+1)/2*(hi-lo)+lo).tolist()); T.extend(b['affinity'].tolist())
        P,T=np.array(P),np.array(T)
        rows.append({'fold':fold,'pearson':pearsonr(T,P)[0],'spearman':spearmanr(T,P)[0],
                     'rmse':float(np.sqrt(np.mean((T-P)**2)))})
    return pd.DataFrame(rows)

def resolve(path, base):
    with open(path) as f: e=yaml.safe_load(f)
    hp=dict(base['hyperparams']); hp.update(e.get('hyperparams',{}))
    return {'name':e['name'],'pooling':e['pooling'],'split':e['split'],
            'n_folds':e.get('n_folds',10),'hp':hp}

def run_experiment(cfg, seeds):
    out=os.path.join(RES,cfg['name']); os.makedirs(out,exist_ok=True)
    with open(os.path.join(out,'config_resolved.yaml'),'w') as f: yaml.safe_dump({**cfg,'seeds':seeds},f,sort_keys=False)
    print(f"\n{'='*66}\n{cfg['name']}  | pooling={cfg['pooling']} split={cfg['split']} seeds={seeds}\n{'='*66}")
    loader=get_loader(cfg['pooling'])
    df=load_saaintdb(CSV)
    per_seed=[]
    for s in seeds:
        c={**DEFAULT_CONFIG,**cfg['hp'],'device':DEVICE,'seed':s}
        setup_reproducibility(s)
        fold_df=run_cv_one(df,loader,c,cfg['n_folds'],cfg['split'])
        sd=os.path.join(out,f'seed_{s}'); os.makedirs(sd,exist_ok=True)
        fold_df.to_csv(os.path.join(sd,'cv_summary.csv'),index=False)
        m={k:fold_df[k].mean() for k in ['pearson','spearman','rmse']}
        per_seed.append(m); print(f"  seed {s}: r={m['pearson']:.4f} rho={m['spearman']:.4f} rmse={m['rmse']:.4f}")
    pm=pd.DataFrame(per_seed)
    agg=pd.DataFrame([{'experiment':cfg['name'],'pooling':cfg['pooling'],'split':cfg['split'],
        'pearson_mean':pm.pearson.mean(),'pearson_std':pm.pearson.std(),
        'spearman_mean':pm.spearman.mean(),'spearman_std':pm.spearman.std(),
        'rmse_mean':pm.rmse.mean(),'rmse_std':pm.rmse.std(),'n_seeds':len(seeds)}])
    agg.to_csv(os.path.join(out,'aggregated_summary.csv'),index=False)
    print(agg.to_string(index=False))
    return agg

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--config'); ap.add_argument('--all',action='store_true')
    ap.add_argument('--seeds',type=int,nargs='+',default=[42,114,144]); a=ap.parse_args()
    base={'hyperparams':{'epochs':50,'patience':10,'batch_size':32,'learning_rate':1e-4,'weight_decay':0.01,
          'projected_size':256,'num_heads':8,'dropout':0.1,'n_layers':2}}
    cfgs=(sorted(glob.glob(os.path.join(EXP,'configs_saaintdb','*.yaml'))) if a.all else [a.config])
    aggs=[run_experiment(resolve(c,base),a.seeds) for c in cfgs]
    master=pd.concat(aggs,ignore_index=True)
    master.to_csv(os.path.join(RES,'MASTER_SUMMARY_SAAINTDB.csv'),index=False)
    print(f"\n{'#'*66}\nMASTER -> experiments/results_saaintdb/MASTER_SUMMARY_SAAINTDB.csv\n{'#'*66}")
    print(master.to_string(index=False))

if __name__=='__main__': main()
