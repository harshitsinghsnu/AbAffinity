"""
build_antigen_clusters.py
==========================
Sequence-identity clustering of SAAINT-DB antigen sequences, to replace the
PDB-ID-only "cold" split with a split that also guarantees no near-identical
antigen crosses the train/test boundary (Reviewer 1 #1 / Reviewer 2 #2).

CD-HIT itself is not available in this environment (no win64 bioconda
package; WSL Ubuntu has no passwordless sudo to install it here), so this
implements a greedy incremental clustering in the same spirit as CD-HIT:
sequences are processed longest-first, and each sequence either joins the
first existing cluster whose representative is within the identity
threshold, or founds a new cluster. Identity is approximated with
difflib.SequenceMatcher.ratio() (Ratcliff/Obershelp), computed only after a
cheap length-ratio prefilter, which is fast enough for ~1,100 antigen
sequences without needing all-pairs alignment.

This is a documented substitute, not literal CD-HIT — flagged explicitly so
the numbers in the response letter aren't overclaimed. If CD-HIT becomes
available, rerun with it and diff the cluster assignments.

Outputs (per threshold):
  antigen_clusters_{pct}.csv   — Ag_seq -> cluster_id, cluster_size
"""
import os, sys, time
import pandas as pd
from difflib import SequenceMatcher

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, 'datasets/saaintdb_with_antigen_names.csv')
OUT = os.path.join(HERE, 'results_antigen_clusters')
os.makedirs(OUT, exist_ok=True)


def greedy_cluster(seqs, threshold):
    """seqs: list of unique sequences. Returns dict seq -> cluster_id (int)."""
    order = sorted(range(len(seqs)), key=lambda i: -len(seqs[i]))
    reps = []          # list of (rep_seq, cluster_id)
    assign = {}
    t0 = time.time()
    for n_done, i in enumerate(order):
        s = seqs[i]
        ls = len(s)
        best = None
        for rep_seq, cid in reps:
            lr = len(rep_seq)
            # length-ratio prefilter: identity cannot exceed min(ls,lr)/max(ls,lr)
            if min(ls, lr) / max(ls, lr) < threshold:
                continue
            sm = SequenceMatcher(None, s, rep_seq, autojunk=False)
            if sm.quick_ratio() < threshold:
                continue
            ident = sm.ratio()
            if ident >= threshold:
                best = cid
                break
        if best is None:
            best = len(reps)
            reps.append((s, best))
        assign[s] = best
        if (n_done + 1) % 100 == 0:
            print(f"    {n_done+1}/{len(seqs)} sequences, {len(reps)} clusters so far "
                  f"({time.time()-t0:.0f}s elapsed)")
    return assign, len(reps)


def main():
    df = pd.read_csv(CSV)
    seqs = sorted(df['Ag_seq'].dropna().unique().tolist())
    print(f"Unique antigen sequences: {len(seqs)}")

    for pct in (90, 70):
        threshold = pct / 100.0
        print(f"\n=== Clustering at {pct}% identity ===")
        assign, n_clusters = greedy_cluster(seqs, threshold)
        print(f"  -> {n_clusters} clusters at {pct}% identity threshold")

        rows = []
        sizes = {}
        for s, cid in assign.items():
            sizes[cid] = sizes.get(cid, 0) + 1
        for s, cid in assign.items():
            rows.append({'Ag_seq': s, 'cluster_id': cid, 'cluster_size': sizes[cid]})
        out_df = pd.DataFrame(rows)
        out_path = os.path.join(OUT, f'antigen_clusters_{pct}.csv')
        out_df.to_csv(out_path, index=False)
        print(f"  Saved: {out_path}")
        print(f"  Singleton clusters: {sum(1 for v in sizes.values() if v == 1)}/{n_clusters}")
        print(f"  Largest cluster size: {max(sizes.values())}")


if __name__ == '__main__':
    main()
