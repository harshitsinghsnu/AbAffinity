"""
precompute_allcdr_cache.py
==========================
One-time: build All-CDR (H1+H2+H3 heavy-pooled) embedding caches for the natural
and mutation datasets, so the multi-seed runner can load small caches instead of
the 5 GB per-residue files each time.

Outputs:
  experiments/cache/allcdr_natural_650M.pkl   (heavy=CDR-pool, light/antigen=mean-pool)
  experiments/cache/allcdr_mutation_650M.pkl
"""
import os, sys, pickle, gc
import numpy as np
from Bio import SeqIO

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # 3_stream/
sys.path.insert(0, HERE)
import pandas as pd
from agabgated.utils.main_symmetric_mean import load_data

CACHE = os.path.join(HERE, 'experiments', 'cache')
os.makedirs(CACHE, exist_ok=True)

IMGT = {*range(27, 39), *range(56, 66), *range(105, 118)}
def cdr_positions(seq):
    try:
        import anarci
        res = anarci.anarci([('q', seq)], scheme='imgt', assign_germline=False,
                            allow=set('ACDEFGHIKLMNPQRSTVWY'))
        if res and res[0] and res[0][0] and res[0][0][0] is not None:
            idx, sp = [], 0
            for (p, _), aa in res[0][0][0]:
                if aa == '-': continue
                if p in IMGT: idx.append(sp)
                sp += 1
            if idx: return sorted(set(idx))
    except Exception: pass
    import re; o = []
    m = re.search(r'C[A-Z]{2,5}[SAGTV]([A-Z]{8,15})WVRQ', seq); o += list(range(*m.span(1))) if m else list(range(26,36))
    m = re.search(r'W[VI]RQ[A-Z]{6,14}W[VL][AS]([A-Z]{10,20})VKGRF', seq); o += list(range(*m.span(1))) if m else list(range(50,64))
    m = re.search(r'WYYCA([A-Z]+)WGQGT', seq) or re.search(r'WYYC[A-Z]([A-Z]+)WGQG', seq); o += list(range(*m.span(1))) if m else list(range(95,110))
    return sorted(set(o))

# sequences (both fasta cover everything)
seqs = {}
for fn in ['datasets/seq_natural.fasta', 'datasets/seq.fasta']:
    for r in SeqIO.parse(os.path.join(HERE, fn), 'fasta'):
        seqs.setdefault(r.id, str(r.seq))
print(f"{len(seqs)} sequences loaded")

def build(per_res_path, mean_path, heavy_ids, out_name):
    out_path = os.path.join(CACHE, out_name)
    if os.path.isfile(out_path):
        print(f"  {out_name} already exists — skipping"); return
    print(f"  Loading {os.path.basename(mean_path)} ...")
    with open(mean_path, 'rb') as f: mean = pickle.load(f)
    print(f"  Loading {os.path.basename(per_res_path)} (large) ...")
    with open(per_res_path, 'rb') as f: pr = pickle.load(f)
    comb = dict(mean)
    n = 0
    for hid in heavy_ids:
        if hid in pr and hid in seqs:
            a = pr[hid]
            ix = [i for i in cdr_positions(seqs[hid]) if 0 <= i < a.shape[0]]
            comb[hid] = (a[ix, :].mean(0) if ix else a.mean(0)).astype(np.float32)
            n += 1
    with open(out_path, 'wb') as f: pickle.dump(comb, f)
    print(f"  CDR-pooled {n} heavy chains -> {out_name}  ({len(comb)} total entries)")
    del pr, mean, comb; gc.collect()

# Natural: sabdab + benchmark
df_sab = load_data(os.path.join(HERE, 'datasets/pairs_sabdab.csv'))
df_ben = load_data(os.path.join(HERE, 'datasets/pairs_benchmark.csv'))
print("\n[natural] sabdab + benchmark")
build(os.path.join(HERE, 'datasets/esm2_per_residue_embeddings_natural_650M.pkl'),
      os.path.join(HERE, 'datasets/esm2_embeddings_natural_650M.pkl'),
      set(df_sab['heavy_id']) | set(df_ben['heavy_id']),
      'allcdr_natural_650M.pkl')

# Mutation: abbind + skempi
df_abb = load_data(os.path.join(HERE, 'datasets/pairs_abbind.csv'))
df_ske = load_data(os.path.join(HERE, 'datasets/pairs_skempi.csv'))
print("\n[mutation] abbind + skempi")
build(os.path.join(HERE, 'datasets/esm2_per_residue_embeddings_mutation_650M.pkl'),
      os.path.join(HERE, 'datasets/esm2_embeddings_mutation_650M.pkl'),
      set(df_abb['heavy_id']) | set(df_ske['heavy_id']),
      'allcdr_mutation_650M.pkl')

print("\nDone. Caches in", CACHE)
