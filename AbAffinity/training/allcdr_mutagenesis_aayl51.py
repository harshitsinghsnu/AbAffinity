"""
allcdr_mutagenesis_aayl51.py
============================
Regenerate the AAYL51 in-silico saturation-mutagenesis heatmap using the
FINAL All-CDR model weights + All-CDR embeddings (not the full-chain model).

Outputs -> advanced_results/
  cdr_heatmap_data_allcdr.csv
  cdr_region_performance_allcdr.csv
"""
import os, sys, pickle
import numpy as np, pandas as pd, torch
from torch.utils.data import Dataset, DataLoader
from scipy.stats import pearsonr, spearmanr
HERE=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,HERE)
from AbAffinity.models.mutual_strong import MutualTriStreamStrong
DEVICE='cuda' if torch.cuda.is_available() else 'cpu'
ADV=os.path.join(HERE,'advanced_results')

df=pd.read_csv(os.path.join(HERE,'data/formatted_aayl51_benchmarking_data.csv'))
df['pKD']=df['Y'].astype(float)
df['heavy_id']=[f'aayl51_H_{i:05d}' for i in range(len(df))]
df['light_id']='aayl51_L_00000'; df['antigen_id']='aayl51_Ag_00000'
with open(os.path.join(ADV,'aayl51_allcdr_embeddings.pkl'),'rb') as f: cdr=pickle.load(f)
with open(os.path.join(HERE,'data/zf_embeddings/aayl51_embeddings.pkl'),'rb') as f: raw=pickle.load(f)
emb={f'aayl51_H_{i:05d}':cdr[df.iloc[i]['heavy']] for i in range(len(df))}
emb['aayl51_L_00000']=raw['light'][0]; emb['aayl51_Ag_00000']=raw['antigen'][0]

class DS(Dataset):
    def __init__(s,d,e): s.d=d.reset_index(drop=True); s.e=e
    def __len__(s): return len(s.d)
    def __getitem__(s,i):
        r=s.d.iloc[i]
        return (torch.tensor(s.e[r['heavy_id']]),torch.tensor(s.e[r['light_id']]),
                torch.tensor(s.e[r['antigen_id']]),torch.tensor(float(r['pKD'])))

ck=torch.load(os.path.join(ADV,'allcdr_indomain_model.pt'),map_location=DEVICE)
cfg=ck.get('config',{}); lo,hi=ck['pkd_bounds']
m=MutualTriStreamStrong(esm_dim=1280,projected_size=cfg.get('projected_size',256),num_heads=cfg.get('num_heads',8),
    dropout=cfg.get('dropout',0.1),n_layers=cfg.get('n_layers',2),device=DEVICE).to(DEVICE)
m.load_state_dict(ck['model_state_dict']); m.eval()
P=[]
with torch.no_grad():
    for h,l,a,_ in DataLoader(DS(df,emb),batch_size=128):
        cos=m(l.to(DEVICE),h.to(DEVICE),a.to(DEVICE))['cosine_similarity'].cpu().numpy()
        P.extend(((cos+1)/2*(hi-lo)+lo).tolist())
df['balm_pred']=np.array(P)

wt=df['heavy'].iloc[0]
def all_muts(s):
    return [(i,wt[i],s[i]) for i in range(min(len(s),len(wt))) if s[i]!=wt[i]]
def region(p):
    if 25<=p<=34: return 'CDR-H1'
    if 49<=p<=62: return 'CDR-H2'
    if 95<=p<=112: return 'CDR-H3'
    return 'WT' if p==-1 else 'Framework'
rows=[]
for _,r in df.iterrows():
    for pos,wa,ma in all_muts(r['heavy']):
        rows.append({'pKD':r['pKD'],'balm_pred':r['balm_pred'],'mut_pos':pos,'wt_aa':wa,'mut_aa':ma,'cdr_region':region(pos)})
ex=pd.DataFrame(rows)
wt_pkd=float(df['pKD'].iloc[0]); wt_pred=float(df['balm_pred'].iloc[0])
ex['delta_pkd_true']=ex['pKD']-wt_pkd; ex['delta_pkd_pred']=ex['balm_pred']-wt_pred
cdr_mut=ex[ex.cdr_region.isin(['CDR-H1','CDR-H2','CDR-H3'])]
hm=cdr_mut.groupby(['mut_pos','mut_aa']).agg(
    delta_pkd_true=('delta_pkd_true','mean'),delta_pkd_pred=('delta_pkd_pred','mean'),
    n=('pKD','count'),wt_aa=('wt_aa','first'),cdr_region=('cdr_region','first')).reset_index()
hm.to_csv(os.path.join(ADV,'cdr_heatmap_data_allcdr.csv'),index=False)
perf=[]
for reg in ['CDR-H1','CDR-H2','CDR-H3']:
    sub=cdr_mut[cdr_mut.cdr_region==reg][['pKD','balm_pred','mut_pos']].drop_duplicates()
    if len(sub)>=10:
        r,_=pearsonr(sub['pKD'],sub['balm_pred']); rho,_=spearmanr(sub['pKD'],sub['balm_pred'])
        perf.append({'cdr':reg,'n_variants':len(sub),'n_positions':sub['mut_pos'].nunique(),
            'pearson_r':r,'spearman_rho':rho,'rmse':float(np.sqrt(np.mean((sub['pKD']-sub['balm_pred'])**2)))})
        print(f"  {reg}: n={len(sub)} r={r:.4f}")
pd.DataFrame(perf).to_csv(os.path.join(ADV,'cdr_region_performance_allcdr.csv'),index=False)
print("Saved All-CDR mutagenesis -> cdr_heatmap_data_allcdr.csv + cdr_region_performance_allcdr.csv")
