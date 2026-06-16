"""
saaintdb_heavy_light_cdr.py
===========================
SAaIntDB All-CDR (heavy+light) CV (random, 3 seeds) for the pooling ablation:
All-CDR(heavy) vs All-CDR(heavy+light) vs full-chain mean-pool.
Light-chain CDRs pooled via ANARCI IMGT; heavy uses existing CDR cache.
Outputs -> experiments/results_allcdr_stats/saaintdb_heavy_light_cv.csv
"""
import os, sys, pickle
import numpy as np, pandas as pd, torch
from torch.utils.data import DataLoader
from scipy.stats import pearsonr, spearmanr
from transformers import EsmTokenizer, EsmModel
HERE=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,HERE)
from agabgated.models.mutual_strong_saaintdb import load_saaintdb, get_fold_splits, _train_fold
from agabgated.utils.main_symmetric_mean import CachedEmbeddingDataset, collate_fn, DEFAULT_CONFIG, setup_reproducibility
DEVICE='cuda' if torch.cuda.is_available() else 'cpu'
OUT=os.path.join(HERE,'experiments','results_allcdr_stats'); os.makedirs(OUT,exist_ok=True)
IMGT={*range(27,39),*range(56,66),*range(105,118)}

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
    import re;o=[]
    m=re.search(r'C[A-Z]{2,5}[SAGTV]([A-Z]{8,15})W[VFY]RQ',seq);o+=list(range(*m.span(1))) if m else list(range(26,36))
    return sorted(set(o)) if o else list(range(26,36))

df=load_saaintdb(os.path.join(HERE,'data/saaintdb_with_antigen_names.csv'))
# light_id -> sequence map
lid2seq={}
for _,r in df.iterrows():
    s=r.get('L_seq')
    if isinstance(s,str) and len(s)>5: lid2seq[r['light_id']]=s
print(f"{len(lid2seq)} light chains to CDR-pool")

cache=os.path.join(HERE,'experiments','cache','saaintdb_light_cdr.pkl')
if os.path.isfile(cache):
    light_cdr=pickle.load(open(cache,'rb'))
else:
    tok=EsmTokenizer.from_pretrained('facebook/esm2_t33_650M_UR50D')
    esm=EsmModel.from_pretrained('facebook/esm2_t33_650M_UR50D').to(DEVICE).eval()
    light_cdr={}; uq=list(dict.fromkeys(lid2seq.values()))
    seq2emb={}
    for i in range(0,len(uq),4):
        b=uq[i:i+4]
        inp=tok(b,return_tensors='pt',padding=True,truncation=True,max_length=512).to(DEVICE)
        with torch.no_grad(): te=esm(**inp).last_hidden_state
        am=inp['attention_mask']
        for j,s in enumerate(b):
            L=int(am[j].sum())-2; ix=[p+1 for p in cdr_idx(s) if 0<=p<L]
            seq2emb[s]=(te[j,ix,:].mean(0) if ix else te[j,1:L+1,:].mean(0)).cpu().numpy().astype(np.float32)
    for lid,s in lid2seq.items(): light_cdr[lid]=seq2emb[s]
    pickle.dump(light_cdr,open(cache,'wb')); del esm; torch.cuda.empty_cache()
    print(f"computed {len(light_cdr)} light-CDR embeddings")

# combined: mean baseline + heavy CDR + light CDR
with open(os.path.join(HERE,'data/esm2_embeddings_saaintdb_650M.pkl'),'rb') as f: mean=pickle.load(f)
with open(os.path.join(HERE,'data/saaintdb_heavy_cdr_embeddings.pkl'),'rb') as f: hcdr=pickle.load(f)
emb=dict(mean); emb.update(hcdr); emb.update(light_cdr)
class DL:
    def __init__(s,d): s.embeddings=d; s.embedding_dim=next(iter(d.values())).shape[0]
    def get_embedding(s,k): return s.embeddings[k]
loader=DL(emb)
idc=['heavy_id','light_id','antigen_id']
df=df[df[idc].apply(lambda c:c.isin(emb)).all(axis=1)].reset_index(drop=True)

rows=[]
for seed in [42,114,144]:
    cfg={**DEFAULT_CONFIG,'device':DEVICE,'seed':seed,'epochs':50,'patience':10,'batch_size':32,
         'learning_rate':1e-4,'weight_decay':0.01,'projected_size':256,'num_heads':8,'dropout':0.1,'n_layers':2}
    setup_reproducibility(seed)
    splits=get_fold_splits(df,10,seed,'random'); fr=[]
    for tr,va in splits:
        dtr=df.iloc[tr].reset_index(drop=True); dva=df.iloc[va].reset_index(drop=True)
        m,b=_train_fold(dtr,dva,loader,cfg,DEVICE); lo,hi=b
        vl=DataLoader(CachedEmbeddingDataset(dva,loader),batch_size=32,collate_fn=collate_fn)
        P,T=[],[]; m.eval()
        with torch.no_grad():
            for bt in vl:
                cos=m(bt['light_emb'].to(DEVICE),bt['heavy_emb'].to(DEVICE),bt['antigen_emb'].to(DEVICE))['cosine_similarity'].cpu().numpy()
                P.extend(((cos+1)/2*(hi-lo)+lo).tolist()); T.extend(bt['affinity'].tolist())
        fr.append((pearsonr(T,P)[0],spearmanr(T,P)[0],float(np.sqrt(np.mean((np.array(T)-np.array(P))**2)))))
    fr=np.array(fr); rows.append({'seed':seed,'pearson':float(fr[:,0].mean()),'spearman':float(fr[:,1].mean()),'rmse':float(fr[:,2].mean())})
    print(f"seed {seed}: r={np.mean(fr):.4f}")
res=pd.DataFrame(rows)
summary={'pooling':'allcdr_heavy_light','split':'random','pearson_mean':res.pearson.mean(),'pearson_std':res.pearson.std(),'spearman_mean':res.spearman.mean(),'spearman_std':res.spearman.std(),'rmse_mean':res.rmse.mean(),'rmse_std':res.rmse.std()}
pd.DataFrame([summary]).to_csv(os.path.join(OUT,'saaintdb_heavy_light_cv.csv'),index=False)
print(f"\nSAINTdb All-CDR(heavy+light) random: r={res.pearson.mean():.4f} ± {res.pearson.std():.4f}")
print("(compare: heavy-only=0.858, mean-pool=0.833)")
