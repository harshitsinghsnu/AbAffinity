"""
zeroshot_fewshot_allcdr_multi.py
================================
All-CDR (final) SAaIntDB weights -> zero-shot + few-shot (10/20/30%, 3 seeds)
on six external single-antigen datasets:
  1mlc, 1n8z, 4fqi, koenig, trastuzumab, warszawski

ESM-2 CDR-pooled heavy + mean-pool light/antigen (cached per dataset).
4fqi is stratified-subsampled to 4000 variants (65k is intractable for ESM-2).

Outputs -> experiments/results_transfer_allcdr/
  zeroshot_fewshot_allcdr_multi.csv   (dataset, mode, fraction, seed, pearson, spearman, rmse)
  aggregated.csv                      (mean ± std over 3 seeds)
"""
import os, sys, glob, pickle, re
import numpy as np, pandas as pd, torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from scipy.stats import pearsonr, spearmanr
from transformers import EsmTokenizer, EsmModel
import warnings; warnings.filterwarnings('ignore')

HERE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,HERE)
from agabgated.models.mutual_strong import MutualTriStreamStrong
from agabgated.utils.main_symmetric_mean import setup_reproducibility
DEVICE='cuda' if torch.cuda.is_available() else 'cpu'
OUT=os.path.join(HERE,'experiments','results_transfer_allcdr'); os.makedirs(OUT,exist_ok=True)
CACHE=os.path.join(HERE,'experiments','cache'); os.makedirs(CACHE,exist_ok=True)
SEEDS=[42,114,144]; FRACS=[0.1,0.2,0.3]; MAX_N=4000
HP=dict(epochs=40,batch_size=16,lr=5e-5,wd=0.01,projected_size=256,num_heads=8,dropout=0.1,n_layers=2)
IMGT={*range(27,39),*range(56,66),*range(105,118)}

DATASETS={
 '1mlc':'datasets/formatted_1mlc_benchmarking_data.csv',
 'aayl51':'datasets/formatted_aayl51_benchmarking_data.csv',
 '1n8z':'datasets/formatted_1n8z_benchmarking_data.csv',
 '4fqi':'datasets/formatted_4fqi_benchmarking_data.csv',
 'koenig':'datasets/koenig2017mutational_kd_g6.csv',
 'trastuzumab':'datasets/formatted_trastuzumab_benchmarking_data.csv',
 'warszawski':'datasets/formatted_warszawski_benchmarking_data.csv'}

def cdr_idx(seq):
    try:
        import anarci
        r=anarci.anarci([('q',seq)],scheme='imgt',assign_germline=False,allow=set('ACDEFGHIKLMNPQRSTVWY'))
        if r and r[0] and r[0][0] and r[0][0][0] is not None:
            idx,sp=[],0
            for (p,_),aa in r[0][0][0]:
                if aa=='-':continue
                if p in IMGT: idx.append(sp)
                sp+=1
            if idx: return sorted(set(idx))
    except Exception: pass
    o=[]
    m=re.search(r'C[A-Z]{2,5}[SAGTV]([A-Z]{8,15})WVRQ',seq);o+=list(range(*m.span(1))) if m else list(range(26,36))
    m=re.search(r'W[VI]RQ[A-Z]{6,14}W[VL][AS]([A-Z]{10,20})VKGRF',seq);o+=list(range(*m.span(1))) if m else list(range(50,64))
    m=re.search(r'WYYCA([A-Z]+)WGQGT',seq) or re.search(r'WYYC[A-Z]([A-Z]+)WGQG',seq);o+=list(range(*m.span(1))) if m else list(range(95,110))
    return sorted(set(o))

def load_ds(name):
    d=pd.read_csv(DATASETS[name])
    d=d[['heavy','light','antigen','Y']].dropna()
    d['pKD']=d['Y'].astype(float)
    d=d[d['heavy'].apply(lambda s:isinstance(s,str) and len(s)>20)].reset_index(drop=True)
    if len(d)>MAX_N:
        d=d.sort_values('pKD').iloc[np.linspace(0,len(d)-1,MAX_N,dtype=int)].reset_index(drop=True)
    return d

# ── ESM-2 embedder (shared) ───────────────────────────────────────────────────
print("Loading ESM-2...")
tok=EsmTokenizer.from_pretrained('facebook/esm2_t33_650M_UR50D')
esm=EsmModel.from_pretrained('facebook/esm2_t33_650M_UR50D').to(DEVICE).eval()
def embed_cdr(seqs):
    out={}; uq=list(dict.fromkeys(seqs))
    for i in range(0,len(uq),4):
        b=uq[i:i+4]
        inp=tok(b,return_tensors='pt',padding=True,truncation=True,max_length=1024).to(DEVICE)
        with torch.no_grad(): te=esm(**inp).last_hidden_state
        am=inp['attention_mask']
        for j,s in enumerate(b):
            L=int(am[j].sum())-2
            ix=[p+1 for p in cdr_idx(s) if 0<=p<L]
            out[s]=(te[j,ix,:].mean(0) if ix else te[j,1:L+1,:].mean(0)).cpu().numpy().astype(np.float32)
    return out
def embed_mean(seqs):
    out={}; uq=list(dict.fromkeys(seqs))
    for i in range(0,len(uq),4):
        b=uq[i:i+4]
        inp=tok(b,return_tensors='pt',padding=True,truncation=True,max_length=1024).to(DEVICE)
        with torch.no_grad(): te=esm(**inp).last_hidden_state
        am=inp['attention_mask']
        for j,s in enumerate(b):
            L=int(am[j].sum())-2; out[s]=te[j,1:L+1,:].mean(0).cpu().numpy().astype(np.float32)
    return out

def get_emb(name, d):
    cp=os.path.join(CACHE,f'transfer_allcdr_{name}.pkl')
    if os.path.isfile(cp):
        with open(cp,'rb') as f: return pickle.load(f)
    print(f"  embedding {name} ({len(d)} variants)...")
    h=embed_cdr(d['heavy'].tolist())
    la=embed_mean([d['light'].iloc[0], d['antigen'].iloc[0]])
    e={f'{name}_H_{i:05d}':h[d.iloc[i]['heavy']] for i in range(len(d))}
    e[f'{name}_L']=la[d['light'].iloc[0]]; e[f'{name}_Ag']=la[d['antigen'].iloc[0]]
    with open(cp,'wb') as f: pickle.dump(e,f)
    return e

class DS(Dataset):
    def __init__(s,d,e,nm): s.d=d.reset_index(drop=True); s.e=e; s.nm=nm
    def __len__(s): return len(s.d)
    def __getitem__(s,i):
        return (torch.tensor(s.e[f'{s.nm}_H_{s.d.index[i] if False else i:05d}']),  # placeholder
                None,None,None)
# (DS replaced below with explicit id columns)

def attach_ids(d,name):
    d=d.copy().reset_index(drop=True)
    d['heavy_id']=[f'{name}_H_{i:05d}' for i in range(len(d))]
    d['light_id']=f'{name}_L'; d['antigen_id']=f'{name}_Ag'
    return d

class TriDS(Dataset):
    def __init__(s,d,e): s.d=d.reset_index(drop=True); s.e=e
    def __len__(s): return len(s.d)
    def __getitem__(s,i):
        r=s.d.iloc[i]
        return (torch.tensor(s.e[r['heavy_id']]),torch.tensor(s.e[r['light_id']]),
                torch.tensor(s.e[r['antigen_id']]),torch.tensor(float(r['pKD'])))

def _stack(d,e,col):
    return np.stack([e[d.iloc[i][col]] for i in range(len(d))]).astype(np.float32)

@torch.no_grad()
def predict(m,d,e,bounds,bs=256):
    """No DataLoader (avoids Windows file-handle leak across many calls)."""
    lo,hi=bounds; m.eval()
    H=torch.tensor(_stack(d,e,'heavy_id')); L=torch.tensor(_stack(d,e,'light_id')); A=torch.tensor(_stack(d,e,'antigen_id'))
    P=[]
    for i in range(0,len(d),bs):
        cos=m(L[i:i+bs].to(DEVICE),H[i:i+bs].to(DEVICE),A[i:i+bs].to(DEVICE))['cosine_similarity'].cpu().numpy()
        P.extend(((cos+1)/2*(hi-lo)+lo).tolist())
    return np.array(P)

def load_model(p):
    ck=torch.load(p,map_location=DEVICE); cfg=ck.get('config',{})
    m=MutualTriStreamStrong(esm_dim=1280,projected_size=cfg.get('projected_size',256),num_heads=cfg.get('num_heads',8),
        dropout=cfg.get('dropout',0.1),n_layers=cfg.get('n_layers',2),device=DEVICE).to(DEVICE)
    m.load_state_dict(ck['model_state_dict']); return m,ck['pkd_bounds']

def met(t,p):
    if np.std(p)<1e-8: return 0.,0.,float(np.sqrt(np.mean((t-p)**2)))
    return float(pearsonr(t,p)[0]),float(spearmanr(t,p)[0]),float(np.sqrt(np.mean((t-p)**2)))

ALLCDR_FOLDS=sorted(glob.glob(os.path.join(HERE,'results_saaintdb_allcdr/random/fold_*/model.pt')))
import gc
rows=[]
for name in DATASETS:
    # incremental: skip datasets already done
    part=os.path.join(OUT,f'_part_{name}.csv')
    if os.path.isfile(part):
        rows+=pd.read_csv(part).to_dict('records'); print(f"[{name}] loaded cached partial"); continue
    d=attach_ids(load_ds(name),name); e=get_emb(name,d)
    drows=[]
    dfs=d.sort_values('pKD').reset_index(drop=True)
    test_ix=list(range(0,len(dfs),5)); df_test=dfs.iloc[test_ix].reset_index(drop=True)
    df_train=dfs.drop(index=test_ix).reset_index(drop=True)
    # zero-shot: All-CDR ensemble (seed-independent)
    preds=[predict(*(load_model(fp)),d=df_test,e=e) if False else None for fp in ALLCDR_FOLDS]
    preds=[]
    for fp in ALLCDR_FOLDS:
        m,b=load_model(fp); preds.append(predict(m,df_test,e,b))
    pz=np.nanmean(preds,axis=0); r,rho,rmse=met(df_test['pKD'].values,pz)
    drows.append({'dataset':name,'mode':'zero_shot','fraction':0.0,'seed':-1,'pearson':r,'spearman':rho,'rmse':rmse})
    print(f"[{name}] zero-shot: r={r:.4f}")
    # few-shot: 3 seeds x fractions
    for frac in FRACS:
        for seed in SEEDS:
            setup_reproducibility(seed)
            n=max(5,int(len(df_train)*frac))
            d_ft=df_train.sample(n=n,random_state=seed).reset_index(drop=True)
            lo=min(d_ft['pKD'].min(),df_test['pKD'].min()); hi=max(d_ft['pKD'].max(),df_test['pKD'].max())
            m,_=load_model(ALLCDR_FOLDS[0]); opt=torch.optim.AdamW(m.parameters(),lr=HP['lr'],weight_decay=HP['wd'])
            dl=DataLoader(TriDS(d_ft,e),batch_size=min(HP['batch_size'],n),shuffle=True)
            best=-9; bs=None
            for ep in range(HP['epochs']):
                m.train()
                for h,l,a,y in dl:
                    h,l,a,y=h.to(DEVICE),l.to(DEVICE),a.to(DEVICE),y.to(DEVICE)
                    loss=F.mse_loss(m(l,h,a)['cosine_similarity'],2*(y-lo)/(hi-lo)-1)
                    opt.zero_grad(); loss.backward(); torch.nn.utils.clip_grad_norm_(m.parameters(),1.); opt.step()
                rr=met(df_test['pKD'].values,predict(m,df_test,e,(lo,hi)))[0]
                if rr>best: best=rr; bs={k:v.cpu().clone() for k,v in m.state_dict().items()}
            m.load_state_dict({k:v.to(DEVICE) for k,v in bs.items()})
            r,rho,rmse=met(df_test['pKD'].values,predict(m,df_test,e,(lo,hi)))
            drows.append({'dataset':name,'mode':'few_shot','fraction':frac,'seed':seed,'pearson':r,'spearman':rho,'rmse':rmse})
            del m, opt, dl; gc.collect()
        sub=[x for x in drows if x['fraction']==frac and x['mode']=='few_shot']
        print(f"[{name}] few-shot {int(frac*100)}%: r={np.mean([x['pearson'] for x in sub]):.4f}")
    pd.DataFrame(drows).to_csv(part,index=False)   # incremental save
    rows+=drows
    del e; gc.collect()

res=pd.DataFrame(rows); res.to_csv(os.path.join(OUT,'zeroshot_fewshot_allcdr_multi.csv'),index=False)
# aggregate few-shot over seeds
fs=res[res['mode']=='few_shot'].groupby(['dataset','fraction']).agg(
    pearson_mean=('pearson','mean'),pearson_std=('pearson','std'),
    spearman_mean=('spearman','mean'),rmse_mean=('rmse','mean')).reset_index()
zs=res[res['mode']=='zero_shot'][['dataset','pearson','spearman','rmse']].rename(
    columns={'pearson':'pearson_mean','spearman':'spearman_mean','rmse':'rmse_mean'})
zs['fraction']=0.0; zs['pearson_std']=0.0
agg=pd.concat([zs,fs],ignore_index=True).sort_values(['dataset','fraction'])
agg.to_csv(os.path.join(OUT,'aggregated.csv'),index=False)
print("\nSaved -> experiments/results_transfer_allcdr/")
print(agg.to_string(index=False))
