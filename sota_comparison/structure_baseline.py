"""
structure_baseline.py
=====================
Structure-based affinity baseline (the comparison ESMFold2 is meant to enable), computed on the
DEPOSITED EXPERIMENTAL complex structures of SAaIntDB. For every complex we measure geometric
antibody-antigen interface descriptors and ask how well they predict pK_d under our exact 10-fold
splits. Because these use the *true* crystal structures, the result is an OPTIMISTIC UPPER BOUND on
any predicted-structure pipeline (ESMFold2 / AlphaFold) which must additionally absorb modelling
error. We then compare against our sequence-only regressor (random 0.858 / cold 0.568).

Phase 1 (slow, network): fetch each unique PDB, compute per-row interface descriptors -> CSV cache.
Phase 2 (fast): RandomForest CV on the descriptors under random + cold splits -> results CSV.

Usage:
  python sota_comparison/structure_baseline.py --phase descriptors   # build the descriptor cache
  python sota_comparison/structure_baseline.py --phase cv            # cross-validate
  python sota_comparison/structure_baseline.py --phase all
"""
import os, sys, argparse, warnings
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(HERE, 'datasets', 'saaintdb_with_antigen_names.csv')
PDBDIR = os.path.join(HERE, 'sota_comparison', 'pdb_cache'); os.makedirs(PDBDIR, exist_ok=True)
DESC = os.path.join(HERE, 'sota_comparison', 'structure_descriptors.csv')

DESC_COLS = ['n_atom_contacts_4', 'n_atom_contacts_5', 'n_hbonds', 'n_iface_res_ab',
             'n_iface_res_ag', 'n_iface_res_tot', 'min_dist', 'mean_contact_dist']


def descriptors():
    from Bio.PDB import PDBList, PDBParser, NeighborSearch
    df = pd.read_csv(CSV)
    pl = PDBList(); parser = PDBParser(QUIET=True)
    cache = {}
    if os.path.exists(DESC):
        cache = {r['idx']: r for _, r in pd.read_csv(DESC).iterrows()}
    rows = []
    by_pdb = df.groupby('PDB_ID')
    done = 0
    for pdb, grp in by_pdb:
        fn = os.path.join(PDBDIR, f'pdb{pdb.lower()}.ent')
        if not os.path.exists(fn):
            try:
                pl.retrieve_pdb_file(pdb, pdir=PDBDIR, file_format='pdb')
            except Exception:
                pass
        st = None
        if os.path.exists(fn):
            try:
                st = parser.get_structure(pdb, fn)[0]
            except Exception:
                st = None
        for i, r in grp.iterrows():
            rec = {'idx': i, 'pKD': r['pKD']}
            if st is None:
                rows.append({**rec, **{c: np.nan for c in DESC_COLS}}); continue
            ab_ch = [c for c in [r['H_chain_ID'], r['L_chain_ID']] if isinstance(c, str)]
            ag_ch = [c for c in str(r['Ag_chain_ID(s)']).split(';') if c and c != 'nan']
            try:
                ab_at = [a for ch in ab_ch if ch in st for a in st[ch].get_atoms()]
                ag_at = [a for ch in ag_ch if ch in st for a in st[ch].get_atoms()]
            except Exception:
                ab_at = ag_at = []
            if not ab_at or not ag_at:
                rows.append({**rec, **{c: np.nan for c in DESC_COLS}}); continue
            ns = NeighborSearch(ag_at)
            c4 = c5 = hb = 0; dists = []; iab = set(); iag = set()
            ag_pol = lambda a: a.element in ('N', 'O')
            for a in ab_at:
                near = ns.search(a.coord, 5.0)
                for b in near:
                    d = float(np.linalg.norm(a.coord - b.coord)); dists.append(d)
                    if d < 5.0:
                        c5 += 1; iab.add(a.get_parent().id[1]); iag.add(b.get_parent().id[1])
                    if d < 4.0:
                        c4 += 1
                    if d < 3.5 and a.element in ('N', 'O') and ag_pol(b):
                        hb += 1
            rec.update({'n_atom_contacts_4': c4, 'n_atom_contacts_5': c5, 'n_hbonds': hb,
                        'n_iface_res_ab': len(iab), 'n_iface_res_ag': len(iag),
                        'n_iface_res_tot': len(iab) + len(iag),
                        'min_dist': min(dists) if dists else np.nan,
                        'mean_contact_dist': float(np.mean([d for d in dists if d < 5])) if any(d < 5 for d in dists) else np.nan})
            rows.append(rec)
        done += 1
        if done % 100 == 0:
            pd.DataFrame(rows).to_csv(DESC, index=False)
            print(f"  {done}/{by_pdb.ngroups} PDBs processed", flush=True)
    pd.DataFrame(rows).to_csv(DESC, index=False)
    print(f"descriptors -> {DESC}  ({len(rows)} rows, "
          f"{pd.DataFrame(rows)[DESC_COLS[0]].notna().sum()} with structure)")


def cv():
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import KFold
    from scipy.stats import pearsonr, spearmanr
    df = pd.read_csv(CSV)
    d = pd.read_csv(DESC).set_index('idx')
    df = df.join(d[DESC_COLS], how='left')
    df = df.dropna(subset=DESC_COLS + ['pKD']).reset_index(drop=True)
    print(f"complexes with usable experimental interface: {len(df)}")
    X = df[DESC_COLS].values; y = df['pKD'].values

    def splits(split_type):
        kf = KFold(10, shuffle=True, random_state=42)
        if split_type == 'cold':
            p = df['PDB_ID'].unique(); out = []
            for tr, va in kf.split(p):
                trp, vap = set(p[tr]), set(p[va])
                out.append((np.where(df['PDB_ID'].isin(trp))[0], np.where(df['PDB_ID'].isin(vap))[0]))
            return out
        return list(kf.split(df))

    allrows = []
    for st in ('random', 'cold'):
        res = []
        for fi, (tr, va) in enumerate(splits(st)):
            m = RandomForestRegressor(n_estimators=300, random_state=42, n_jobs=-1)
            m.fit(X[tr], y[tr]); p = m.predict(X[va])
            res.append((pearsonr(y[va], p)[0], spearmanr(y[va], p)[0],
                        float(np.sqrt(np.mean((y[va] - p) ** 2)))))
        res = np.array(res)
        print(f"[{st}] structure-RF  Pearson {res[:,0].mean():.3f}±{res[:,0].std():.3f} | "
              f"Spearman {res[:,1].mean():.3f} | RMSE {res[:,2].mean():.3f}")
        for fi, r in enumerate(res):
            allrows.append({'split': st, 'fold': fi, 'pearson': r[0], 'spearman': r[1], 'rmse': r[2]})
    # single-descriptor correlations (interpretability)
    print("\nSingle-descriptor |Pearson| vs pKd:")
    for c in DESC_COLS:
        print(f"  {c:20s} {abs(pearsonr(df[c], y)[0]):.3f}")
    out = os.path.join(HERE, 'sota_comparison', 'results_StructureBaseline.csv')
    pd.DataFrame(allrows).to_csv(out, index=False); print("saved", out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--phase', default='all',
                                                     choices=['descriptors', 'cv', 'all'])
    a = ap.parse_args()
    if a.phase in ('descriptors', 'all'):
        descriptors()
    if a.phase in ('cv', 'all'):
        cv()
