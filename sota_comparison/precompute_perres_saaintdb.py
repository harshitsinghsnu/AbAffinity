"""
precompute_perres_saaintdb.py
=============================
Per-residue ESM-2 650M embeddings for every unique chain sequence in SAaIntDB, so the recent
sequence-based SOTA baselines (DuaDeep-SeqAffinity, DG-Affinity) can be trained on the SAME data
and SAME features as our model. Stored as float16 to keep the cache small.

Out -> experiments/cache/perres_saaintdb_650M.pkl   {sequence_string: float16 (L,1280)}
"""
import os, sys, pickle, time
import numpy as np, pandas as pd, torch
from transformers import AutoModel, AutoTokenizer

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(HERE, 'datasets', 'saaintdb_with_antigen_names.csv')
OUT = os.path.join(HERE, 'experiments', 'cache', 'perres_saaintdb_650M.pkl')
MODEL = 'facebook/esm2_t33_650M_UR50D'
MAXLEN = 512
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'

def main():
    df = pd.read_csv(CSV)
    seqs = set()
    for col in ('H_seq', 'L_seq', 'Ag_seq'):
        seqs |= {s for s in df[col].dropna() if isinstance(s, str) and s}
    seqs = sorted(seqs)
    cache = {}
    if os.path.exists(OUT):
        cache = pickle.load(open(OUT, 'rb'))
    todo = [s for s in seqs if s not in cache]
    print(f"{len(seqs)} unique sequences, {len(todo)} to embed, device={DEVICE}")
    tok = AutoTokenizer.from_pretrained(MODEL)
    esm = AutoModel.from_pretrained(MODEL).to(DEVICE).eval()
    t0 = time.time()
    for i, s in enumerate(todo):
        with torch.no_grad():
            x = tok(s[:MAXLEN], return_tensors='pt').to(DEVICE)
            emb = esm(**x).last_hidden_state.squeeze(0)[1:-1].cpu().numpy().astype(np.float16)
        cache[s] = emb
        if (i + 1) % 200 == 0:
            pickle.dump(cache, open(OUT, 'wb'))
            print(f"  {i+1}/{len(todo)}  ({(time.time()-t0)/(i+1):.2f}s/seq)")
    pickle.dump(cache, open(OUT, 'wb'))
    print(f"Saved {len(cache)} -> {OUT}  ({time.time()-t0:.0f}s)")

if __name__ == '__main__':
    main()
