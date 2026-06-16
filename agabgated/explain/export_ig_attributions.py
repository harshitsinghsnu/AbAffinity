"""
export_ig_attributions.py
=========================
Compute & save raw IG attribution arrays + heavy-chain CDR spans for both
3-stream (final tri-stream) and two-stream models on the SAME 5 canonical
complexes, so the explainability notebook can build CDR-annotated comparison
figures and quantify CDR enrichment.

Outputs -> experiments/results_explainability/
  {complex}.npz  with: three_heavy, three_light, three_antigen,
                       two_heavy,  two_light,  two_antigen,
                       cdr_h1, cdr_h2, cdr_h3 (heavy index spans), heavy_len
  cdr_enrichment.csv  (per complex/model: fraction of |attr| in CDRs)
"""
import os, sys, re
import numpy as np, pandas as pd, torch
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from agabgated.explain import integrated_gradients_explainability_saaintdb as ig3
from agabgated.explain import integrated_gradient_twostream as ig2
from agabgated.models.two_stream_mutualstrong import ConcatTwoStream
import torch
DEVICE = ig3.DEVICE
OUT = os.path.join(HERE,'experiments','results_explainability'); os.makedirs(OUT, exist_ok=True)
COMPLEXES = ig2.COMPLEXES

# IMGT CDR spans on heavy chain (per-sequence via ANARCI)
def cdr_spans(seq):
    H1=set(range(27,39)); H2=set(range(56,66)); H3=set(range(105,118))
    out={'h1':[], 'h2':[], 'h3':[]}
    try:
        import anarci
        r=anarci.anarci([('q',seq)],scheme='imgt',assign_germline=False,allow=set('ACDEFGHIKLMNPQRSTVWY'))
        if r and r[0] and r[0][0] and r[0][0][0] is not None:
            sp=0
            for (p,_),aa in r[0][0][0]:
                if aa=='-':continue
                if p in H1: out['h1'].append(sp)
                elif p in H2: out['h2'].append(sp)
                elif p in H3: out['h3'].append(sp)
                sp+=1
            if out['h1'] or out['h2'] or out['h3']: return out
    except Exception: pass
    m=re.search(r'C[A-Z]{2,5}[SAGTV]([A-Z]{8,15})WVRQ',seq); out['h1']=list(range(*m.span(1))) if m else list(range(26,36))
    m=re.search(r'W[VI]RQ[A-Z]{6,14}W[VL][AS]([A-Z]{10,20})VKGRF',seq); out['h2']=list(range(*m.span(1))) if m else list(range(50,64))
    m=re.search(r'WYYCA([A-Z]+)WGQGT',seq) or re.search(r'WYYC[A-Z]([A-Z]+)WGQG',seq); out['h3']=list(range(*m.span(1))) if m else list(range(95,110))
    return out

def two_load(path):
    ck=torch.load(path,map_location=DEVICE); cfg=ck.get('config',{})
    m=ConcatTwoStream(esm_dim=1280,projected_size=cfg.get('projected_size',256),num_heads=cfg.get('num_heads',8),
        dropout=0.0,n_layers=cfg.get('n_layers',2),device=DEVICE).to(DEVICE)
    m.load_state_dict(ck['model_state_dict']); m.eval()
    return m, ck.get('pkd_bounds',(5.0,12.0))

print("Loading models + ESM-2...")
m3,b3=ig3.load_model("model_weights/saaintdb_allcdr_random_bestfold.pt")  # FINAL All-CDR model (random best fold)
w3=ig3.MutualTriStreamPerResidue(m3,b3).to(DEVICE).eval()
m2,b2=two_load("results_twostream_saaintdb/fold_01/model.pt")
w2=ig2.TwoStreamPerResidue(m2,b2).to(DEVICE).eval()
emb=ig3.ESM2PerResidueEmbedder(ig3.ESM2_MODEL,DEVICE)

def enrich(attr, spans, hlen):
    a=np.abs(attr[:hlen]); tot=a.sum()+1e-9
    cdr=sorted(set(spans['h1']+spans['h2']+spans['h3']))
    cdr=[i for i in cdr if 0<=i<hlen]
    cdr_mass=a[cdr].sum() if cdr else 0.0
    cdr_frac_len=len(cdr)/hlen if hlen else 0
    return float(cdr_mass/tot), float(cdr_frac_len)

rows=[]
for name,info in COMPLEXES.items():
    print(f"\n=== {name} ===")
    h,l,a=info['heavy'],info['light'],info['antigen']
    he=emb.embed(h); le=emb.embed(l); ae=emb.embed(a)
    spans=cdr_spans(h); hlen=he.shape[0]
    # 3-stream: returns (light, heavy, antigen)
    la3,ha3,aa3=ig3.compute_ig_attributions(w3,le,he,ae)
    # 2-stream: returns (heavy, light, antigen)
    ha2,la2,aa2=ig2.compute_ig_attributions(w2,he,le,ae)
    np.savez(os.path.join(OUT,f'{name}.npz'),
        three_heavy=ha3,three_light=la3,three_antigen=aa3,
        two_heavy=ha2,two_light=la2,two_antigen=aa2,
        cdr_h1=np.array(spans['h1']),cdr_h2=np.array(spans['h2']),cdr_h3=np.array(spans['h3']),
        heavy_len=hlen)
    e3,frac=enrich(ha3,spans,hlen); e2,_=enrich(ha2,spans,hlen)
    rows.append({'complex':name,'model':'three_stream','cdr_attr_frac':e3,'cdr_len_frac':frac,'enrichment':e3/(frac+1e-9)})
    rows.append({'complex':name,'model':'two_stream','cdr_attr_frac':e2,'cdr_len_frac':frac,'enrichment':e2/(frac+1e-9)})
    print(f"  CDR attr-mass: 3-stream={e3:.3f}  2-stream={e2:.3f}  (CDR len frac={frac:.3f})")
emb.cleanup()
pd.DataFrame(rows).to_csv(os.path.join(OUT,'cdr_enrichment.csv'),index=False)
print("\nSaved attributions + cdr_enrichment.csv ->",OUT)
