"""
allcdr_stream_interventions.py
==============================
Post-hoc stream-intervention ablation on the FINAL All-CDR SAaIntDB model
(supplementary): full vs zero-light vs zero-heavy vs zero-antigen vs
mean-antigen vs shuffled-antigen. Uses trained All-CDR folds; no retraining.
Also per-CDR-region note. Outputs -> experiments/results_allcdr_stats/
"""
import os, sys, glob, pickle
import numpy as np, pandas as pd, torch
from scipy.stats import pearsonr, spearmanr
HERE=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,HERE)
from agabgated.models.mutual_strong import MutualTriStreamStrong
from agabgated.models.mutual_strong_saaintdb import load_saaintdb, get_fold_splits
DEVICE='cuda' if torch.cuda.is_available() else 'cpu'
OUT=os.path.join(HERE,'experiments','results_allcdr_stats'); os.makedirs(OUT,exist_ok=True)

# All-CDR combined embeddings (heavy CDR + mean light/antigen)
with open(os.path.join(HERE,'data/esm2_embeddings_saaintdb_650M.pkl'),'rb') as f: mean=pickle.load(f)
with open(os.path.join(HERE,'data/saaintdb_heavy_cdr_embeddings.pkl'),'rb') as f: cdr=pickle.load(f)
emb=dict(mean); emb.update(cdr)
df=load_saaintdb(os.path.join(HERE,'data/saaintdb_with_antigen_names.csv'))
idc=['heavy_id','light_id','antigen_id']
df=df[df[idc].apply(lambda c:c.isin(emb)).all(axis=1)].reset_index(drop=True)
# global mean antigen vector (for mean_antigen intervention)
ag_mat=np.stack([emb[a] for a in df['antigen_id'].unique()]); ag_mean=ag_mat.mean(0).astype(np.float32)

def load_fold(fp):
    ck=torch.load(fp,map_location=DEVICE); cfg=ck.get('config',{})
    m=MutualTriStreamStrong(esm_dim=1280,projected_size=cfg.get('projected_size',256),num_heads=cfg.get('num_heads',8),
        dropout=cfg.get('dropout',0.1),n_layers=cfg.get('n_layers',2),device=DEVICE).to(DEVICE)
    m.load_state_dict(ck['model_state_dict']); m.eval(); return m, ck['pkd_bounds']

@torch.no_grad()
def evalcond(m,dv,bounds,cond,seed):
    lo,hi=bounds
    H=np.stack([emb[i] for i in dv['heavy_id']]).astype(np.float32)
    L=np.stack([emb[i] for i in dv['light_id']]).astype(np.float32)
    A=np.stack([emb[i] for i in dv['antigen_id']]).astype(np.float32)
    if cond=='zero_light': L=np.zeros_like(L)
    elif cond=='zero_heavy': H=np.zeros_like(H)
    elif cond=='zero_antigen': A=np.zeros_like(A)
    elif cond=='mean_antigen': A=np.tile(ag_mean,(len(A),1))
    elif cond=='shuffled_antigen':
        rng=np.random.default_rng(seed); A=A[rng.permutation(len(A))]
    P=[]
    for i in range(0,len(dv),256):
        cos=m(torch.tensor(L[i:i+256]).to(DEVICE),torch.tensor(H[i:i+256]).to(DEVICE),
              torch.tensor(A[i:i+256]).to(DEVICE))['cosine_similarity'].cpu().numpy()
        P.extend(((cos+1)/2*(hi-lo)+lo).tolist())
    t=dv['binding_affinity'].values; p=np.array(P)
    return float(pearsonr(t,p)[0]), float(spearmanr(t,p)[0]), float(np.sqrt(np.mean((t-p)**2)))

conds=['full','zero_light','zero_heavy','zero_antigen','mean_antigen','shuffled_antigen']
folds=sorted(glob.glob(os.path.join(HERE,'results_saaintdb_allcdr/random/fold_*/model.pt')))
splits=get_fold_splits(df,10,9999,'random')   # matches All-CDR training seed
rows=[]
for fi,(fp,(tr,va)) in enumerate(zip(folds,splits),1):
    dv=df.iloc[va].reset_index(drop=True)
    dv=dv[dv[idc].apply(lambda c:c.isin(emb)).all(axis=1)].reset_index(drop=True)
    m,b=load_fold(fp)
    for c in conds:
        r,rho,rmse=evalcond(m,dv,b,c,9999+fi)
        rows.append({'fold':fi,'condition':c,'pearson':r,'spearman':rho,'rmse':rmse})
    print(f"fold {fi}: "+" ".join(f"{c}={[x for x in rows if x['fold']==fi and x['condition']==c][0]['pearson']:.3f}" for c in conds))
res=pd.DataFrame(rows)
agg=res.groupby('condition').agg(pearson_mean=('pearson','mean'),pearson_std=('pearson','std'),
    spearman_mean=('spearman','mean'),rmse_mean=('rmse','mean')).reindex(conds).reset_index()
res.to_csv(os.path.join(OUT,'allcdr_stream_interventions_folds.csv'),index=False)
agg.to_csv(os.path.join(OUT,'allcdr_stream_interventions.csv'),index=False)
print("\nAll-CDR stream interventions (mean over 10 folds):")
print(agg.to_string(index=False))
