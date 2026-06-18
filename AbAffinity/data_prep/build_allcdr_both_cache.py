"""
build_allcdr_both_cache.py
==========================
Build 'All-CDR (heavy+light)' embedding caches: CDR-pool BOTH heavy and light
chains (IMGT L1/L2/L3 + H1/H2/H3 via ANARCI), antigen stays mean-pool.
Compare against heavy-only All-CDR and mean-pool.

Outputs:
  data/allcdrboth_natural_650M.pkl
  data/allcdrboth_mutation_650M.pkl
"""
import os, sys, pickle, gc
import numpy as np
from Bio import SeqIO
HERE=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,HERE)
from AbAffinity.utils.main_symmetric_mean import load_data
CACHE=os.path.join(HERE,'data'); os.makedirs(CACHE,exist_ok=True)

# IMGT CDR ranges are identical position windows for H and L chains
IMGT={*range(27,39),*range(56,66),*range(105,118)}
def cdr_positions(seq):
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
    m=re.search(r'W[VIFY]RQ[A-Z]{6,16}([A-Z]{7,20})[FW]G',seq);o+=list(range(*m.span(1))) if m else list(range(50,64))
    return sorted(set(o)) if o else list(range(26,36))

seqs={}
for fn in ['data/seq_natural.fasta','data/seq.fasta']:
    for r in SeqIO.parse(os.path.join(HERE,fn),'fasta'): seqs.setdefault(r.id,str(r.seq))
print(f"{len(seqs)} sequences")

def build(per_res_path, mean_path, hl_ids, out_name):
    out=os.path.join(CACHE,out_name)
    if os.path.isfile(out): print(f"  {out_name} exists, skip"); return
    print(f"  loading {os.path.basename(mean_path)}");
    with open(mean_path,'rb') as f: mean=pickle.load(f)
    print(f"  loading {os.path.basename(per_res_path)} (large)")
    with open(per_res_path,'rb') as f: pr=pickle.load(f)
    comb=dict(mean); n=0
    for cid in hl_ids:
        if cid in pr and cid in seqs:
            a=pr[cid]; ix=[i for i in cdr_positions(seqs[cid]) if 0<=i<a.shape[0]]
            comb[cid]=(a[ix,:].mean(0) if ix else a.mean(0)).astype(np.float32); n+=1
    with open(out,'wb') as f: pickle.dump(comb,f)
    print(f"  CDR-pooled {n} heavy+light chains -> {out_name} ({len(comb)} entries)")
    del pr,mean,comb; gc.collect()

# natural: sabdab heavy+light (+benchmark)
df_sab=load_data(os.path.join(HERE,'data/pairs_sabdab.csv'))
df_ben=load_data(os.path.join(HERE,'data/pairs_benchmark.csv'))
nat_ids=set(df_sab['heavy_id'])|set(df_sab['light_id'])|set(df_ben['heavy_id'])|set(df_ben['light_id'])
print("\n[natural] heavy+light")
build(os.path.join(HERE,'data/esm2_per_residue_embeddings_natural_650M.pkl'),
      os.path.join(HERE,'data/esm2_embeddings_natural_650M.pkl'),nat_ids,'allcdrboth_natural_650M.pkl')

# mutation: abbind+skempi heavy+light
df_abb=load_data(os.path.join(HERE,'data/pairs_abbind.csv'))
df_ske=load_data(os.path.join(HERE,'data/pairs_skempi.csv'))
mut_ids=set(df_abb['heavy_id'])|set(df_abb['light_id'])|set(df_ske['heavy_id'])|set(df_ske['light_id'])
print("\n[mutation] heavy+light")
build(os.path.join(HERE,'data/esm2_per_residue_embeddings_mutation_650M.pkl'),
      os.path.join(HERE,'data/esm2_embeddings_mutation_650M.pkl'),mut_ids,'allcdrboth_mutation_650M.pkl')
print("\nDone.")
