"""
allcdr_natural_oof_preds.py
===========================
Generate All-CDR out-of-fold predictions for SAbDab / AbBind / SKEMPI (10-fold,
seed 42) for the publication scatter plots.
Outputs -> experiments/results_allcdr_stats/oof_preds_{ds}.csv  (true, pred)
"""
import os, sys, pickle
import numpy as np, pandas as pd, torch
from torch.utils.data import DataLoader
from scipy.stats import pearsonr, spearmanr
HERE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,HERE)
from agabgated.models.mutual_strong_saaintdb import get_fold_splits, _train_fold
from agabgated.utils.main_symmetric_mean import load_data, CachedEmbeddingDataset, collate_fn, DEFAULT_CONFIG, setup_reproducibility
DEVICE='cuda' if torch.cuda.is_available() else 'cpu'
OUT=os.path.join(HERE,'experiments','results_allcdr_stats'); os.makedirs(OUT,exist_ok=True)
DS={'sabdab':('datasets/pairs_sabdab.csv','natural'),
    'abbind':('datasets/pairs_abbind.csv','mutation'),
    'skempi':('datasets/pairs_skempi.csv','mutation')}
class DL:
    def __init__(s,d): s.embeddings=d; s.embedding_dim=next(iter(d.values())).shape[0]
    def get_embedding(s,k): return s.embeddings[k]
_cache={}
def loader(fam):
    if fam in _cache: return _cache[fam]
    with open(os.path.join(HERE,f'experiments/cache/allcdr_{fam}_650M.pkl'),'rb') as f: d=pickle.load(f)
    _cache[fam]=DL(d); return _cache[fam]

for ds,(csv,fam) in DS.items():
    ld=loader(fam); df=load_data(os.path.join(HERE,csv))
    idc=['heavy_id','light_id','antigen_id']
    df=df[df[idc].apply(lambda c:c.isin(ld.embeddings)).all(axis=1)].reset_index(drop=True)
    setup_reproducibility(42)
    cfg={**DEFAULT_CONFIG,'device':DEVICE,'seed':42,'epochs':50,'patience':10,'batch_size':32,
         'learning_rate':1e-4,'weight_decay':0.01,'projected_size':256,'num_heads':8,'dropout':0.1,'n_layers':2}
    T,P,FD=[],[],[]
    for fi,(tr,va) in enumerate(get_fold_splits(df,10,42,'random')):
        dtr=df.iloc[tr].reset_index(drop=True); dva=df.iloc[va].reset_index(drop=True)
        m,b=_train_fold(dtr,dva,ld,cfg,DEVICE); lo,hi=b
        vl=DataLoader(CachedEmbeddingDataset(dva,ld),batch_size=64,collate_fn=collate_fn)
        m.eval()
        with torch.no_grad():
            for bt in vl:
                cos=m(bt['light_emb'].to(DEVICE),bt['heavy_emb'].to(DEVICE),bt['antigen_emb'].to(DEVICE))['cosine_similarity'].cpu().numpy()
                pr=((cos+1)/2*(hi-lo)+lo).tolist()
                P.extend(pr); T.extend(bt['affinity'].tolist()); FD.extend([fi]*len(pr))
    T,P=np.array(T),np.array(P)
    pd.DataFrame({'true':T,'pred':P,'fold':FD}).to_csv(os.path.join(OUT,f'oof_preds_{ds}.csv'),index=False)
    print(f"{ds}: n={len(T)} r={pearsonr(T,P)[0]:.4f} rho={spearmanr(T,P)[0]:.4f} rmse={np.sqrt(np.mean((T-P)**2)):.4f}")
print("Saved OOF preds -> experiments/results_allcdr_stats/")
