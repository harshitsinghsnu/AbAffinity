"""
benchmark_nlayers_check.py
==========================
Train All-CDR on full SAbDab, test on the independent Benchmark set, comparing
n_layers = 1 vs n_layers = 2 (3 seeds). Decides whether n_layers=1 transfers
better for the benchmark experiment.
Outputs -> experiments/results_allcdr_stats/benchmark_nlayers.csv
"""
import os, sys, pickle
import numpy as np, pandas as pd
HERE=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,HERE); sys.path.insert(0,os.path.join(HERE,'experiments'))
from AbAffinity.training import run_multiseed as R
from AbAffinity.utils.main_symmetric_mean import load_data
OUT=os.path.join(HERE,'experiments','results_allcdr_stats'); os.makedirs(OUT,exist_ok=True)

class DL:
    def __init__(s,d): s.embeddings=d; s.embedding_dim=next(iter(d.values())).shape[0]
    def get_embedding(s,k): return s.embeddings[k]
emb=pickle.load(open(os.path.join(HERE,'data/allcdr_natural_650M.pkl'),'rb'))
loader=DL(emb)
df_tr=load_data(os.path.join(HERE,'data/pairs_sabdab.csv'))
df_b =load_data(os.path.join(HERE,'data/pairs_benchmark.csv'))

# base hyperparams (same as benchmark experiment), only n_layers varies
base_hp=dict(projected_size=256,num_heads=8,dropout=0.1,batch_size=32,epochs=50,
             learning_rate=1e-4,weight_decay=0.01,patience=10)
SEEDS=[314,114,144]
rows=[]
for nl in [1,2]:
    hp={**base_hp,'n_layers':nl}
    per=[]
    for s in SEEDS:
        res,_,_=R.run_benchmark_one(df_tr,df_b,loader,hp,s)
        per.append(res); print(f"n_layers={nl} seed={s}: r={res['pearson']:.4f} rho={res['spearman']:.4f} rmse={res['rmse']:.4f}")
    pm=pd.DataFrame(per)
    rows.append({'n_layers':nl,'pearson_mean':pm.pearson.mean(),'pearson_std':pm.pearson.std(),
                 'spearman_mean':pm.spearman.mean(),'spearman_std':pm.spearman.std(),
                 'rmse_mean':pm.rmse.mean(),'rmse_std':pm.rmse.std()})
res=pd.DataFrame(rows); res.to_csv(os.path.join(OUT,'benchmark_nlayers.csv'),index=False)
print("\n=== SAbDab->Benchmark (All-CDR) n_layers comparison ==="); print(res.to_string(index=False))
