"""
allcdr_final_ablations.py
=========================
All-CDR adopted as the FINAL model. Produces, multi-seed (3 seeds), 10-fold CV:

  1. Architecture comparison  (Ours / Two-stream / Concat+MLP) x (All-CDR / mean-pool)
     -> experiments/results/allcdr_architecture_comparison.csv
  2. Post-hoc gating ablation on the All-CDR Ours model
     (learned / open=1 / closed=0 / fixed=0.5 / random gate)
     -> experiments/results/allcdr_gating_ablation.csv
  3. All-CDR out-of-fold mutational-impact predictions (SKEMPI, AbBind)
     -> experiments/results/allcdr_mutational_preds_{skempi,abbind}.csv

Uses cached embeddings (experiments/cache/allcdr_*.pkl + datasets/esm2_embeddings_*).
"""
import os, sys, pickle, types
import numpy as np, pandas as pd, torch
import torch.nn as nn, torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from scipy.stats import pearsonr, spearmanr
import warnings; warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from agabgated.models.mutual_strong import MutualTriStreamStrong, GatedCrossAttention
from agabgated.models.mutual_strong_saaintdb import get_fold_splits
from agabgated.utils.main_symmetric_mean import load_data, setup_reproducibility

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
RES = os.path.join(HERE, 'experiments', 'results'); os.makedirs(RES, exist_ok=True)
CACHE = os.path.join(HERE, 'experiments', 'cache')
SEEDS = [42, 114, 144]
HP = dict(epochs=50, patience=10, batch_size=32, lr=1e-4, wd=0.01,
          projected_size=256, num_heads=8, dropout=0.1, n_layers=2)

DS = {'sabdab': ('datasets/pairs_sabdab.csv', 'natural'),
      'abbind': ('datasets/pairs_abbind.csv', 'mutation'),
      'skempi': ('datasets/pairs_skempi.csv', 'mutation')}

# ── embedding caches ─────────────────────────────────────────────────────────
_EMB = {}
def emb(pooling, family):
    key = (pooling, family)
    if key in _EMB: return _EMB[key]
    path = (os.path.join(CACHE, f'allcdr_{family}_650M.pkl') if pooling == 'allcdr'
            else os.path.join(HERE, f'datasets/esm2_embeddings_{family}_650M.pkl'))
    with open(path, 'rb') as f: _EMB[key] = pickle.load(f)
    return _EMB[key]

def vec(d, k):
    return np.asarray(d[k], dtype=np.float32)

# ── datasets ──────────────────────────────────────────────────────────────────
class TriDS(Dataset):
    def __init__(s, df, e): s.df=df.reset_index(drop=True); s.e=e
    def __len__(s): return len(s.df)
    def __getitem__(s, i):
        r=s.df.iloc[i]
        return (torch.tensor(vec(s.e,r['heavy_id'])), torch.tensor(vec(s.e,r['light_id'])),
                torch.tensor(vec(s.e,r['antigen_id'])), torch.tensor(float(r['binding_affinity'])))

# ── models ────────────────────────────────────────────────────────────────────
class ConcatMLP(nn.Module):
    def __init__(s, input_dim, hidden=(1280,256), dropout=0.1):
        super().__init__(); layers=[]; c=input_dim
        for h in hidden:
            layers += [nn.Linear(c,h), nn.BatchNorm1d(h), nn.ReLU(), nn.Dropout(dropout)]; c=h
        layers.append(nn.Linear(c,1)); s.net=nn.Sequential(*layers)
    def forward(s,x): return s.net(x).squeeze(-1)

class TwoStream(nn.Module):
    """Antibody = mean(heavy, light); bidirectional gated cross-attn; cosine head."""
    def __init__(s, esm_dim=1280, projected_size=256, num_heads=8, dropout=0.1):
        super().__init__()
        s.ab_proj=nn.Sequential(nn.Linear(esm_dim,projected_size),nn.LayerNorm(projected_size),nn.GELU(),nn.Dropout(dropout))
        s.ag_proj=nn.Sequential(nn.Linear(esm_dim,projected_size),nn.LayerNorm(projected_size),nn.GELU(),nn.Dropout(dropout))
        s.ab2ag=GatedCrossAttention(projected_size,projected_size,num_heads,dropout)
        s.ag2ab=GatedCrossAttention(projected_size,projected_size,num_heads,dropout)
        s.head_ab=nn.Sequential(nn.Linear(projected_size,projected_size),nn.LayerNorm(projected_size),nn.GELU(),nn.Dropout(dropout))
        s.head_ag=nn.Sequential(nn.Linear(projected_size,projected_size),nn.LayerNorm(projected_size),nn.GELU(),nn.Dropout(dropout))
    def forward(s, ab, ag):
        ab=F.layer_norm(ab,ab.shape[-1:]); ag=F.layer_norm(ag,ag.shape[-1:])
        ap=s.ab_proj(ab); gp=s.ag_proj(ag)
        ab_ctx=s.ag2ab(ap,gp); ag_ctx=s.ab2ag(gp,ap)
        a=F.normalize(s.head_ab(ab_ctx),dim=-1); g=F.normalize(s.head_ag(ag_ctx),dim=-1)
        return (a*g).sum(-1)

# ── training helpers ──────────────────────────────────────────────────────────
def metrics(t,p):
    if np.std(p)<1e-8: return 0.,0.,float(np.sqrt(np.mean((t-p)**2)))
    return float(pearsonr(t,p)[0]), float(spearmanr(t,p)[0]), float(np.sqrt(np.mean((t-p)**2)))

def train_cosine(model, df_tr, df_va, e, kind, seed):
    """kind in {'ours','two'} -> cosine model with bounds; returns (r,rho,rmse)."""
    setup_reproducibility(seed)
    lo,hi=df_tr['binding_affinity'].min(),df_tr['binding_affinity'].max()
    opt=torch.optim.AdamW(model.parameters(),lr=HP['lr'],weight_decay=HP['wd'])
    dl=DataLoader(TriDS(df_tr,e),batch_size=HP['batch_size'],shuffle=True)
    best=-np.inf; best_state=None; bad=0
    for ep in range(HP['epochs']):
        model.train()
        for h,l,a,y in dl:
            h,l,a,y=h.to(DEVICE),l.to(DEVICE),a.to(DEVICE),y.to(DEVICE)
            tgt=2*(y-lo)/(hi-lo)-1
            cos = model(l,h,a) if kind=='ours' else model((h+l)/2,a)
            cos = cos['cosine_similarity'] if isinstance(cos,dict) else cos
            loss=F.mse_loss(cos,tgt)
            opt.zero_grad(); loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(),1.0); opt.step()
        r,_,_=eval_cosine(model,df_va,e,kind,(lo,hi))
        if r>best: best=r; best_state={k:v.cpu().clone() for k,v in model.state_dict().items()}; bad=0
        else:
            bad+=1
            if bad>=HP['patience']: break
    model.load_state_dict({k:v.to(DEVICE) for k,v in best_state.items()})
    return (lo,hi)

@torch.no_grad()
def eval_cosine(model, df, e, kind, bounds, gate_mode=None):
    lo,hi=bounds; model.eval()
    dl=DataLoader(TriDS(df,e),batch_size=64,shuffle=False)
    P,T=[],[]
    for h,l,a,y in dl:
        h,l,a=h.to(DEVICE),l.to(DEVICE),a.to(DEVICE)
        if gate_mode is not None: set_gate(model, gate_mode)
        cos = model(l,h,a) if kind=='ours' else model((h+l)/2,a)
        cos = cos['cosine_similarity'] if isinstance(cos,dict) else cos
        cos=cos.cpu().numpy(); P.extend(((cos+1)/2*(hi-lo)+lo).tolist()); T.extend(y.tolist())
    if gate_mode is not None: set_gate(model, None)  # restore
    return metrics(np.array(T),np.array(P))[0], np.array(T), np.array(P)

def train_concat(df_tr, df_va, e, seed):
    setup_reproducibility(seed)
    dim=len(vec(e,df_tr.iloc[0]['heavy_id']))*3
    model=ConcatMLP(dim,dropout=HP['dropout']).to(DEVICE)
    opt=torch.optim.AdamW(model.parameters(),lr=HP['lr'],weight_decay=HP['wd'])
    def batches(df):
        idx=np.arange(len(df));
        for i in range(0,len(df),HP['batch_size']):
            b=df.iloc[idx[i:i+HP['batch_size']]]
            X=np.stack([np.concatenate([vec(e,r['heavy_id']),vec(e,r['light_id']),vec(e,r['antigen_id'])]) for _,r in b.iterrows()])
            yield torch.tensor(X,dtype=torch.float32).to(DEVICE), torch.tensor(b['binding_affinity'].values,dtype=torch.float32).to(DEVICE)
    best=-np.inf; best_state=None; bad=0
    for ep in range(HP['epochs']):
        model.train()
        for X,y in batches(df_tr.sample(frac=1,random_state=ep)):
            if len(X)<2: continue
            opt.zero_grad(); loss=F.mse_loss(model(X),y); loss.backward(); opt.step()
        # eval
        model.eval(); P,T=[],[]
        with torch.no_grad():
            for X,y in batches(df_va):
                P.extend(model(X).cpu().numpy().tolist()); T.extend(y.cpu().numpy().tolist())
        r=metrics(np.array(T),np.array(P))[0]
        if r>best: best=r; best_state={k:v.cpu().clone() for k,v in model.state_dict().items()}; bad=0
        else:
            bad+=1
            if bad>=HP['patience']: break
    model.load_state_dict({k:v.to(DEVICE) for k,v in best_state.items()})
    model.eval(); P,T=[],[]
    with torch.no_grad():
        for X,y in batches(df_va):
            P.extend(model(X).cpu().numpy().tolist()); T.extend(y.cpu().numpy().tolist())
    return metrics(np.array(T),np.array(P)), np.array(T), np.array(P)

# ── gate override (post-hoc) ─────────────────────────────────────────────────
def set_gate(model, mode):
    """Monkeypatch every GatedCrossAttention.forward to use a fixed gate value.
       mode: None=restore learned; 'open'=1; 'closed'=0; 'fixed'=0.5; 'random'."""
    for m in model.modules():
        if isinstance(m, GatedCrossAttention):
            if mode is None:
                if hasattr(m,'_orig_forward'): m.forward=m._orig_forward
                continue
            if not hasattr(m,'_orig_forward'): m._orig_forward=m.forward
            def make(mod, mode):
                def fwd(query, key_value):
                    residual=query
                    q=mod.W_q(query); k=mod.W_k(key_value); v=mod.W_v(key_value)
                    inter=torch.tanh(q*k); sdpa=inter*v
                    if mode=='open': gate=torch.ones_like(sdpa)
                    elif mode=='closed': gate=torch.zeros_like(sdpa)
                    elif mode=='fixed': gate=0.5*torch.ones_like(sdpa)
                    elif mode=='random': gate=torch.rand_like(sdpa)
                    else: gate=torch.sigmoid(mod.W_gate(query))
                    out=mod.W_o(sdpa*gate); out=mod.dropout(out)
                    return mod.layer_norm(residual+out)
                return fwd
            m.forward=make(m, mode)

# ─────────────────────────────────────────────────────────────────────────────
# 1. ARCHITECTURE COMPARISON
# ─────────────────────────────────────────────────────────────────────────────
print("="*70); print("1. ARCHITECTURE COMPARISON (All-CDR vs mean-pool)"); print("="*70)
arch_rows=[]
for ds,(csv,fam) in DS.items():
    df=load_data(os.path.join(HERE,csv))
    for pooling in ['allcdr','meanpool']:
        e=emb(pooling,fam)
        idc=['heavy_id','light_id','antigen_id']
        df_ok=df[df[idc].apply(lambda c:c.isin(e)).all(axis=1)].reset_index(drop=True)
        for model_name in ['Ours','Two-stream','Concat+MLP']:
            seed_r=[]; seed_rho=[]; seed_rmse=[]
            for seed in SEEDS:
                splits=get_fold_splits(df_ok,10,seed,'random')
                fr=[];frho=[];frmse=[]
                for tr,va in splits:
                    d_tr=df_ok.iloc[tr].reset_index(drop=True); d_va=df_ok.iloc[va].reset_index(drop=True)
                    if len(d_tr)<5 or len(d_va)<2: continue
                    if model_name=='Concat+MLP':
                        (r,rho,rmse),_,_=train_concat(d_tr,d_va,e,seed)
                    else:
                        kind='ours' if model_name=='Ours' else 'two'
                        if kind=='ours':
                            mdl=MutualTriStreamStrong(esm_dim=1280,projected_size=HP['projected_size'],
                                num_heads=HP['num_heads'],dropout=HP['dropout'],n_layers=HP['n_layers'],device=DEVICE).to(DEVICE)
                        else:
                            mdl=TwoStream().to(DEVICE)
                        bounds=train_cosine(mdl,d_tr,d_va,e,kind,seed)
                        r,T,P=eval_cosine(mdl,d_va,e,kind,bounds)
                        rho=metrics(T,P)[1]; rmse=metrics(T,P)[2]
                    fr.append(r);frho.append(rho);frmse.append(rmse)
                seed_r.append(np.mean(fr)); seed_rho.append(np.mean(frho)); seed_rmse.append(np.mean(frmse))
            arch_rows.append({'dataset':ds,'pooling':pooling,'model':model_name,
                'pearson_mean':np.mean(seed_r),'pearson_std':np.std(seed_r),
                'spearman_mean':np.mean(seed_rho),'spearman_std':np.std(seed_rho),
                'rmse_mean':np.mean(seed_rmse),'rmse_std':np.std(seed_rmse),'n_seeds':len(SEEDS)})
            print(f"  {ds:7s} {pooling:8s} {model_name:11s}: r={np.mean(seed_r):.4f}±{np.std(seed_r):.4f}")
pd.DataFrame(arch_rows).to_csv(os.path.join(RES,'allcdr_architecture_comparison.csv'),index=False)

# ─────────────────────────────────────────────────────────────────────────────
# 2. GATING ABLATION (post-hoc, All-CDR Ours)
# ─────────────────────────────────────────────────────────────────────────────
print("\n"+"="*70); print("2. GATING ABLATION (post-hoc, All-CDR Ours)"); print("="*70)
gate_modes=['learned','fixed','open','closed','random']
gate_rows=[]
for ds in ['sabdab','skempi']:
    csv,fam=DS[ds]; e=emb('allcdr',fam); df=load_data(os.path.join(HERE,csv))
    idc=['heavy_id','light_id','antigen_id']
    df_ok=df[df[idc].apply(lambda c:c.isin(e)).all(axis=1)].reset_index(drop=True)
    mode_seed={m:[] for m in gate_modes}
    for seed in SEEDS:
        splits=get_fold_splits(df_ok,10,seed,'random')
        mode_fold={m:[] for m in gate_modes}
        for tr,va in splits:
            d_tr=df_ok.iloc[tr].reset_index(drop=True); d_va=df_ok.iloc[va].reset_index(drop=True)
            if len(d_tr)<5 or len(d_va)<2: continue
            mdl=MutualTriStreamStrong(esm_dim=1280,projected_size=HP['projected_size'],num_heads=HP['num_heads'],
                dropout=HP['dropout'],n_layers=HP['n_layers'],device=DEVICE).to(DEVICE)
            bounds=train_cosine(mdl,d_tr,d_va,e,'ours',seed)
            for gm in gate_modes:
                r,_,_=eval_cosine(mdl,d_va,e,'ours',bounds,gate_mode=(None if gm=='learned' else gm))
                mode_fold[gm].append(r)
        for gm in gate_modes: mode_seed[gm].append(np.mean(mode_fold[gm]))
    for gm in gate_modes:
        gate_rows.append({'dataset':ds,'gate_mode':gm,
            'pearson_mean':np.mean(mode_seed[gm]),'pearson_std':np.std(mode_seed[gm]),'n_seeds':len(SEEDS)})
        print(f"  {ds:7s} gate={gm:8s}: r={np.mean(mode_seed[gm]):.4f}±{np.std(mode_seed[gm]):.4f}")
pd.DataFrame(gate_rows).to_csv(os.path.join(RES,'allcdr_gating_ablation.csv'),index=False)

# ─────────────────────────────────────────────────────────────────────────────
# 3. MUTATIONAL IMPACT (All-CDR out-of-fold preds, SKEMPI + AbBind)
# ─────────────────────────────────────────────────────────────────────────────
print("\n"+"="*70); print("3. MUTATIONAL IMPACT (All-CDR out-of-fold preds)"); print("="*70)
for ds in ['skempi','abbind']:
    csv,fam=DS[ds]; e=emb('allcdr',fam); df=load_data(os.path.join(HERE,csv))
    idc=['heavy_id','light_id','antigen_id']
    df_ok=df[df[idc].apply(lambda c:c.isin(e)).all(axis=1)].reset_index(drop=True)
    splits=get_fold_splits(df_ok,10,SEEDS[0],'random')
    recs=[]
    for tr,va in splits:
        d_tr=df_ok.iloc[tr].reset_index(drop=True); d_va=df_ok.iloc[va].reset_index(drop=True)
        if len(d_tr)<5 or len(d_va)<2: continue
        mdl=MutualTriStreamStrong(esm_dim=1280,projected_size=HP['projected_size'],num_heads=HP['num_heads'],
            dropout=HP['dropout'],n_layers=HP['n_layers'],device=DEVICE).to(DEVICE)
        bounds=train_cosine(mdl,d_tr,d_va,e,'ours',SEEDS[0])
        _,T,P=eval_cosine(mdl,d_va,e,'ours',bounds)
        for i,(t,p) in enumerate(zip(T,P)):
            recs.append({'heavy_id':d_va.iloc[i]['heavy_id'],'true':t,'pred':p,'error':t-p})
    out=pd.DataFrame(recs); out.to_csv(os.path.join(RES,f'allcdr_mutational_preds_{ds}.csv'),index=False)
    r,rho,rmse=metrics(out['true'].values,out['pred'].values)
    print(f"  {ds}: n={len(out)} r={r:.4f} MAE={out['error'].abs().mean():.4f}")

print("\nAll All-CDR final ablations saved to experiments/results/")
