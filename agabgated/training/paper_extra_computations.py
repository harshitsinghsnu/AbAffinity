"""
paper_extra_computations.py
===========================
(1) SAaIntDB gating ablation (post-hoc) with ALL metrics (Pearson/Spearman/RMSE).
(2) Regenerate per-antigen & per-PDB intra-target tables WITH RMSE for the
    "Overall vs Fisher-aggregated" bar plots.
Outputs -> experiments/results_allcdr_stats/
"""
import os, sys, glob, pickle
import numpy as np, pandas as pd, torch
from scipy.stats import pearsonr, spearmanr
HERE=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,HERE)
from agabgated.models.mutual_strong import MutualTriStreamStrong, GatedCrossAttention
from agabgated.models.mutual_strong_saaintdb import load_saaintdb, get_fold_splits
DEVICE='cuda' if torch.cuda.is_available() else 'cpu'
OUT=os.path.join(HERE,'experiments','results_allcdr_stats'); os.makedirs(OUT,exist_ok=True)

# All-CDR SAINTdb embeddings
with open(os.path.join(HERE,'data/esm2_embeddings_saaintdb_650M.pkl'),'rb') as f: mean=pickle.load(f)
with open(os.path.join(HERE,'data/saaintdb_heavy_cdr_embeddings.pkl'),'rb') as f: cdr=pickle.load(f)
emb=dict(mean); emb.update(cdr)
df=load_saaintdb(os.path.join(HERE,'data/saaintdb_with_antigen_names.csv'))
idc=['heavy_id','light_id','antigen_id']
df=df[df[idc].apply(lambda c:c.isin(emb)).all(axis=1)].reset_index(drop=True)

def set_gate(model, mode):
    for m in model.modules():
        if isinstance(m, GatedCrossAttention):
            if mode is None:
                if hasattr(m,'_orig'): m.forward=m._orig
                continue
            if not hasattr(m,'_orig'): m._orig=m.forward
            def mk(mod,mode):
                def fwd(q,kv):
                    res=q; Q=mod.W_q(q); K=mod.W_k(kv); V=mod.W_v(kv)
                    sd=torch.tanh(Q*K)*V
                    if mode=='open': g=torch.ones_like(sd)
                    elif mode=='closed': g=torch.zeros_like(sd)
                    elif mode=='fixed': g=0.5*torch.ones_like(sd)
                    elif mode=='random': g=torch.rand_like(sd)
                    else: g=torch.sigmoid(mod.W_gate(q))
                    return mod.layer_norm(res+mod.dropout(mod.W_o(sd*g)))
                return fwd
            m.forward=mk(m,mode)

def load_fold(fp):
    ck=torch.load(fp,map_location=DEVICE); cfg=ck.get('config',{})
    m=MutualTriStreamStrong(esm_dim=1280,projected_size=cfg.get('projected_size',256),num_heads=cfg.get('num_heads',8),
        dropout=cfg.get('dropout',0.1),n_layers=cfg.get('n_layers',2),device=DEVICE).to(DEVICE)
    m.load_state_dict(ck['model_state_dict']); m.eval(); return m, ck['pkd_bounds']

@torch.no_grad()
def evalmodel(m,dv,bounds,gate):
    lo,hi=bounds
    H=np.stack([emb[i] for i in dv['heavy_id']]).astype(np.float32)
    L=np.stack([emb[i] for i in dv['light_id']]).astype(np.float32)
    A=np.stack([emb[i] for i in dv['antigen_id']]).astype(np.float32)
    set_gate(m, None if gate=='learned' else gate); P=[]
    for i in range(0,len(dv),256):
        cos=m(torch.tensor(L[i:i+256]).to(DEVICE),torch.tensor(H[i:i+256]).to(DEVICE),
              torch.tensor(A[i:i+256]).to(DEVICE))['cosine_similarity'].cpu().numpy()
        P.extend(((cos+1)/2*(hi-lo)+lo).tolist())
    set_gate(m,None)
    t=dv['binding_affinity'].values; p=np.array(P)
    return float(pearsonr(t,p)[0]),float(spearmanr(t,p)[0]),float(np.sqrt(np.mean((t-p)**2)))

# (1) SAINTdb gating, all metrics
print("[1] SAaIntDB gating ablation (all metrics)...")
gate_modes=['learned','fixed','open','closed','random']
folds=sorted(glob.glob(os.path.join(HERE,'results_saaintdb_allcdr/random/fold_*/model.pt')))
splits=get_fold_splits(df,10,9999,'random'); rows=[]
for fi,(fp,(tr,va)) in enumerate(zip(folds,splits),1):
    dv=df.iloc[va].reset_index(drop=True)
    dv=dv[dv[idc].apply(lambda c:c.isin(emb)).all(axis=1)].reset_index(drop=True)
    m,b=load_fold(fp)
    for gm in gate_modes:
        r,rho,rmse=evalmodel(m,dv,b,gm)
        rows.append({'fold':fi,'gate_mode':gm,'pearson':r,'spearman':rho,'rmse':rmse})
res=pd.DataFrame(rows)
agg=res.groupby('gate_mode').agg(pearson_mean=('pearson','mean'),pearson_std=('pearson','std'),
    spearman_mean=('spearman','mean'),spearman_std=('spearman','std'),
    rmse_mean=('rmse','mean'),rmse_std=('rmse','std')).reindex(gate_modes).reset_index()
agg['dataset']='saaintdb'
agg.to_csv(os.path.join(OUT,'allcdr_gating_saaintdb.csv'),index=False)
print(agg.to_string(index=False))

# (2) per-antigen / per-PDB with RMSE + overall-vs-Fisher summary
print("\n[2] per-antigen & per-PDB intra-target (with RMSE) + Fisher aggregation...")
d=pd.read_csv(os.path.join(HERE,'results_saaintdb_allcdr/random/all_preds.csv'))
t_all=d['binding_affinity'].values; p_all=d['predicted_affinity'].values
overall=dict(pearson=float(pearsonr(t_all,p_all)[0]),spearman=float(spearmanr(t_all,p_all)[0]),
             rmse=float(np.sqrt(np.mean((t_all-p_all)**2))))
def group_metrics(col):
    out=[]
    for k,g in d.groupby(col):
        if len(g)>=3 and g['binding_affinity'].std()>1e-6:
            gt,gp=g['binding_affinity'].values,g['predicted_affinity'].values
            try:
                out.append({col:str(k)[:40],'n':len(g),'pearson':pearsonr(gt,gp)[0],
                            'spearman':spearmanr(gt,gp)[0],'rmse':float(np.sqrt(np.mean((gt-gp)**2)))})
            except: pass
    return pd.DataFrame(out)
ag=group_metrics('Antigen_Name'); pdb=group_metrics('PDB_ID')
ag.to_csv(os.path.join(OUT,'allcdr_per_antigen_intratarget.csv'),index=False)
pdb.to_csv(os.path.join(OUT,'allcdr_per_pdb.csv'),index=False)
def fisher_mean(vals):
    z=np.arctanh(np.clip(np.asarray(vals,float),-0.999,0.999))
    return float(np.tanh(np.nanmean(z))), float(np.nanstd(np.tanh(z)))
summ=[]
for grp,tab in [('per_antigen',ag),('per_pdb',pdb)]:
    for met in ['pearson','spearman']:
        fm,fs=fisher_mean(tab[met].dropna())
        summ.append({'grouping':grp,'metric':met,'overall':overall[met],
                     'fisher_mean':fm,'fisher_std':fs,'n_groups':int(tab[met].notna().sum())})
    summ.append({'grouping':grp,'metric':'rmse','overall':overall['rmse'],
                 'fisher_mean':float(tab['rmse'].mean()),'fisher_std':float(tab['rmse'].std()),
                 'n_groups':int(tab['rmse'].notna().sum())})
pd.DataFrame(summ).to_csv(os.path.join(OUT,'overall_vs_fisher.csv'),index=False)
print(pd.DataFrame(summ).to_string(index=False))
print("\nSaved -> experiments/results_allcdr_stats/")
