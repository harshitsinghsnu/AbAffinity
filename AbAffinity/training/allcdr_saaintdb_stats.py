"""
allcdr_saaintdb_stats.py
========================
Statistical analyses for the FINAL All-CDR SAaIntDB model from existing
out-of-fold predictions (results_saaintdb_allcdr/random/all_preds.csv):
  1. Bootstrap 95% CI + permutation p-value (pearson/spearman/rmse)
  2. Per-PDB metrics
  3. Per-antigen intra-target metrics
  4. AA 20x20 mutation-impact (true ΔpKd, predicted ΔpKd, error) from AAYL51 All-CDR

Outputs -> experiments/results_allcdr_stats/
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr
HERE=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT=os.path.join(HERE,'experiments','results_allcdr_stats'); os.makedirs(OUT,exist_ok=True)
rng=np.random.default_rng(0)

# ── 1+2+3 from SAINTdb All-CDR OOF predictions ──────────────────────────────
d=pd.read_csv(os.path.join(HERE,'results_saaintdb_allcdr/random/all_preds.csv'))
t=d['binding_affinity'].values; p=d['predicted_affinity'].values
def metrics(t,p):
    return dict(pearson=float(pearsonr(t,p)[0]),spearman=float(spearmanr(t,p)[0]),
                rmse=float(np.sqrt(np.mean((t-p)**2))),mae=float(np.mean(np.abs(t-p))))
pt=metrics(t,p)
# bootstrap CI
B=1000; boot={k:[] for k in pt}
n=len(t)
for _ in range(B):
    idx=rng.integers(0,n,n); bt,bp=t[idx],p[idx]
    boot['pearson'].append(pearsonr(bt,bp)[0]); boot['spearman'].append(spearmanr(bt,bp)[0])
    boot['rmse'].append(np.sqrt(np.mean((bt-bp)**2))); boot['mae'].append(np.mean(np.abs(bt-bp)))
# permutation p-value (label shuffle) for pearson
perm=[]
for _ in range(B):
    perm.append(abs(pearsonr(t,rng.permutation(p))[0]))
pval=(np.sum(np.array(perm)>=abs(pt['pearson']))+1)/(B+1)
rows=[]
for k in ['pearson','spearman','rmse','mae']:
    arr=np.array(boot[k])
    rows.append({'model':'All-CDR_random','metric':k,'point':pt[k],
                 'ci95_low':float(np.percentile(arr,2.5)),'ci95_high':float(np.percentile(arr,97.5)),
                 'permutation_pvalue':float(pval) if k=='pearson' else np.nan,'n':n})
pd.DataFrame(rows).to_csv(os.path.join(OUT,'allcdr_bootstrap_ci.csv'),index=False)
print("Bootstrap CI (All-CDR random):"); print(pd.DataFrame(rows).to_string(index=False))

# per-PDB
pdb_rows=[]
for pdb,g in d.groupby('PDB_ID'):
    if len(g)>=3:
        gt,gp=g['binding_affinity'].values,g['predicted_affinity'].values
        try: r=pearsonr(gt,gp)[0]; rho=spearmanr(gt,gp)[0]
        except: r=rho=np.nan
        pdb_rows.append({'pdb':pdb,'n':len(g),'pearson':r,'spearman':rho,
                         'rmse':float(np.sqrt(np.mean((gt-gp)**2)))})
pd.DataFrame(pdb_rows).to_csv(os.path.join(OUT,'allcdr_per_pdb.csv'),index=False)
print(f"\nPer-PDB groups (n>=3): {len(pdb_rows)}")

# per-antigen intra-target
ag_rows=[]
for ag,g in d.groupby('Antigen_Name'):
    if len(g)>=3 and g['binding_affinity'].std()>1e-6:
        gt,gp=g['binding_affinity'].values,g['predicted_affinity'].values
        try: r=pearsonr(gt,gp)[0]; rho=spearmanr(gt,gp)[0]
        except: r=rho=np.nan
        ag_rows.append({'antigen':str(ag)[:40],'n':len(g),'pearson':r,'spearman':rho})
pd.DataFrame(ag_rows).to_csv(os.path.join(OUT,'allcdr_per_antigen_intratarget.csv'),index=False)
print(f"Per-antigen intra-target groups (n>=3): {len(ag_rows)}")

# ── 4. AA 20x20 mutation impact from AAYL51 All-CDR ──────────────────────────
import pickle, torch
sys.path.insert(0,HERE)
from AbAffinity.models.mutual_strong import MutualTriStreamStrong
DEVICE='cuda' if torch.cuda.is_available() else 'cpu'
ADV=os.path.join(HERE,'advanced_results')
df=pd.read_csv(os.path.join(HERE,'data/formatted_aayl51_benchmarking_data.csv')); df['pKD']=df['Y'].astype(float)
df['heavy_id']=[f'aayl51_H_{i:05d}' for i in range(len(df))]; df['light_id']='aayl51_L_00000'; df['antigen_id']='aayl51_Ag_00000'
with open(os.path.join(ADV,'aayl51_allcdr_embeddings.pkl'),'rb') as f: cdr=pickle.load(f)
with open(os.path.join(HERE,'data/zf_embeddings/aayl51_embeddings.pkl'),'rb') as f: raw=pickle.load(f)
emb={f'aayl51_H_{i:05d}':cdr[df.iloc[i]['heavy']] for i in range(len(df))}; emb['aayl51_L_00000']=raw['light'][0]; emb['aayl51_Ag_00000']=raw['antigen'][0]
ck=torch.load(os.path.join(ADV,'allcdr_indomain_model.pt'),map_location=DEVICE); cfg=ck.get('config',{}); lo,hi=ck['pkd_bounds']
m=MutualTriStreamStrong(esm_dim=1280,projected_size=cfg.get('projected_size',256),num_heads=cfg.get('num_heads',8),dropout=cfg.get('dropout',0.1),n_layers=cfg.get('n_layers',2),device=DEVICE).to(DEVICE)
m.load_state_dict(ck['model_state_dict']); m.eval()
H=np.stack([emb[f'aayl51_H_{i:05d}'] for i in range(len(df))]).astype(np.float32)
L=np.tile(emb['aayl51_L_00000'],(len(df),1)).astype(np.float32); A=np.tile(emb['aayl51_Ag_00000'],(len(df),1)).astype(np.float32)
P=[]
with torch.no_grad():
    for i in range(0,len(df),256):
        cos=m(torch.tensor(L[i:i+256]).to(DEVICE),torch.tensor(H[i:i+256]).to(DEVICE),torch.tensor(A[i:i+256]).to(DEVICE))['cosine_similarity'].cpu().numpy()
        P.extend(((cos+1)/2*(hi-lo)+lo).tolist())
df['pred']=np.array(P)
wt=df['heavy'].iloc[0]; wt_pkd=float(df['pKD'].iloc[0]); wt_pred=float(df['pred'].iloc[0])
recs=[]
for _,r in df.iterrows():
    for i in range(min(len(r['heavy']),len(wt))):
        if r['heavy'][i]!=wt[i]:
            recs.append({'wt_aa':wt[i],'mut_aa':r['heavy'][i],
                         'dtrue':r['pKD']-wt_pkd,'dpred':r['pred']-wt_pred,
                         'err':(r['pKD']-wt_pkd)-(r['pred']-wt_pred)})
ex=pd.DataFrame(recs)
AAs=list('RHKDESTNQAVILMFYWCGP')  # ordered by category (charged/polar/hydrophobic/special)
def mat(col):
    M=np.full((20,20),np.nan)
    g=ex.groupby(['wt_aa','mut_aa'])[col].mean()
    for i,w in enumerate(AAs):
        for j,mm in enumerate(AAs):
            if (w,mm) in g.index: M[i,j]=g[(w,mm)]
    return M
np.savez(os.path.join(OUT,'aa_mutation_impact.npz'),
         true=mat('dtrue'),pred=mat('dpred'),err=mat('err'),aas=np.array(AAs))
ex.to_csv(os.path.join(OUT,'aa_mutation_events.csv'),index=False)
print(f"\nAA mutation impact: {len(ex)} events, 20x20 matrices saved (true/pred/err).")
print("Saved -> experiments/results_allcdr_stats/")
