"""
gate_importance_saaintdb.py
===========================
Extract learned gate activations g = sigma(W_gate q) from the trained SAaIntDB
All-CDR model (the antibody->antigen gated cross-attention), averaged over all
complexes, reshaped to a 16x16 latent grid. Shows the gate's sparse, structured
amplify/suppress pattern (the denoising operator).
Outputs -> experiments/results_allcdr_stats/gate_importance_saaintdb.npz
"""
import os, sys, pickle
import numpy as np, torch
from torch.utils.data import DataLoader
HERE=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,HERE)
from agabgated.explain import integrated_gradients_explainability_saaintdb as ig3
from agabgated.models.mutual_strong_saaintdb import load_saaintdb
from agabgated.utils.main_symmetric_mean import CachedEmbeddingDataset, collate_fn
DEVICE=ig3.DEVICE
OUT=os.path.join(HERE,'experiments','results_allcdr_stats'); os.makedirs(OUT,exist_ok=True)

# final All-CDR model (same weights used for explainability)
model,bounds=ig3.load_model("model_weights/saaintdb_allcdr_random_bestfold.pt")
model=model.to(DEVICE).eval()

# All-CDR embeddings
with open(os.path.join(HERE,'data/esm2_embeddings_saaintdb_650M.pkl'),'rb') as f: mean=pickle.load(f)
with open(os.path.join(HERE,'data/saaintdb_heavy_cdr_embeddings.pkl'),'rb') as f: cdr=pickle.load(f)
emb=dict(mean); emb.update(cdr)
class DL:
    def __init__(s,d): s.embeddings=d; s.embedding_dim=next(iter(d.values())).shape[0]
    def get_embedding(s,k): return s.embeddings[k]
loader=DL(emb)
df=load_saaintdb(os.path.join(HERE,'data/saaintdb_with_antigen_names.csv'))
idc=['heavy_id','light_id','antigen_id']
df=df[df[idc].apply(lambda c:c.isin(emb)).all(axis=1)].reset_index(drop=True)

# hook gate of antibody->antigen GCA (primary prediction pathway), layer 0
gca=model.ab_to_ag_layers[0]
gates=[]
def hook(mod,inp,out): gates.append(torch.sigmoid(out).detach().cpu().numpy())
h=gca.W_gate.register_forward_hook(hook)

dl=DataLoader(CachedEmbeddingDataset(df,loader),batch_size=64,collate_fn=collate_fn)
with torch.no_grad():
    for b in dl:
        model(b['light_emb'].to(DEVICE),b['heavy_emb'].to(DEVICE),b['antigen_emb'].to(DEVICE))
h.remove()
G=np.concatenate(gates,0)              # (N, 256)
gmean=G.mean(0); gstd=G.std(0)
side=int(round(np.sqrt(gmean.shape[0])))   # 16
print(f"gate dims={gmean.shape[0]} -> {side}x{side} grid; "
      f"mean={gmean.mean():.3f}, frac suppressed(<0.4)={np.mean(gmean<0.4):.2f}, "
      f"frac amplified(>0.6)={np.mean(gmean>0.6):.2f}")
np.savez(os.path.join(OUT,'gate_importance_saaintdb.npz'),
         gate_mean=gmean,gate_std=gstd,grid=gmean.reshape(side,side),side=side,n=len(G))
print("saved gate_importance_saaintdb.npz")
