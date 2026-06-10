"""
run_sota_cv.py
==============
Train/evaluate the re-implemented recent SOTA baselines (DuaDeep-SeqAffinity, DG-Affinity) on
SAaIntDB under the SAME 10-fold splits (random + cold/no-PDB-overlap, seed 42) and the SAME frozen
ESM-2 650M features as our model, then report Pearson / Spearman / RMSE for an apples-to-apples
comparison against "Ours".

Usage:
  python sota_comparison/run_sota_cv.py --model DuaDeep-SeqAffinity --split random
  (flags: --folds N to run a subset of folds, --max_n N to subsample rows (smoke test),
   --epochs E, --reduced/--full, --out path)
"""
import os, sys, json, time, argparse, pickle
import numpy as np, pandas as pd, torch
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import KFold
from scipy.stats import pearsonr, spearmanr

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sota_models import MODELS

CSV = os.path.join(HERE, 'datasets', 'saaintdb_with_antigen_names.csv')
CACHE = os.path.join(HERE, 'experiments', 'cache', 'perres_saaintdb_650M.pkl')
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
MAXLEN = 512   # overridden by --maxlen


def get_fold_splits(df, n, seed, split_type):
    kf = KFold(n_splits=n, shuffle=True, random_state=seed)
    if split_type == 'cold':
        pdbs = df['PDB_ID'].unique()
        out = []
        for tr, va in kf.split(pdbs):
            trp, vap = set(pdbs[tr]), set(pdbs[va])
            out.append((np.where(df['PDB_ID'].isin(trp))[0], np.where(df['PDB_ID'].isin(vap))[0]))
        return out
    return list(kf.split(df))


class DS(Dataset):
    def __init__(self, df, cache):
        self.rows = df.reset_index(drop=True); self.cache = cache
    def __len__(self): return len(self.rows)
    def __getitem__(self, i):
        r = self.rows.iloc[i]
        h = self.cache[r['H_seq']]
        ab = h
        if isinstance(r['L_seq'], str) and r['L_seq'] in self.cache:
            ab = np.concatenate([h, self.cache[r['L_seq']]], axis=0)
        ag = self.cache[r['Ag_seq']]
        return (torch.from_numpy(ab[:MAXLEN].astype(np.float32)),
                torch.from_numpy(ag[:MAXLEN].astype(np.float32)),
                float(r['pKD']))


def pad(batch_seqs):
    L = max(s.shape[0] for s in batch_seqs)
    x = torch.zeros(len(batch_seqs), L, batch_seqs[0].shape[1])
    m = torch.zeros(len(batch_seqs), L)
    for i, s in enumerate(batch_seqs):
        x[i, :s.shape[0]] = s; m[i, :s.shape[0]] = 1
    return x, m


def collate(b):
    ab = [x[0] for x in b]; ag = [x[1] for x in b]; y = torch.tensor([x[2] for x in b])
    abx, abm = pad(ab); agx, agm = pad(ag)
    return abx, abm, agx, agm, y


def evaluate(model, dl):
    model.eval(); P, T = [], []
    with torch.no_grad():
        for abx, abm, agx, agm, y in dl:
            p = model(abx.to(DEVICE), abm.to(DEVICE), agx.to(DEVICE), agm.to(DEVICE))
            P += p.cpu().tolist(); T += y.tolist()
    P, T = np.array(P), np.array(T)
    return P, T


def metrics(P, T):
    return (float(pearsonr(T, P)[0]), float(spearmanr(T, P)[0]),
            float(np.sqrt(np.mean((T - P) ** 2))))


def train_fold(ModelCls, dtr, dva, cache, reduced, epochs, d_model=256, bs=32, patience=6, lr=3e-4):
    torch.manual_seed(42)
    model = ModelCls(reduced=reduced, d_model=d_model).to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    lossf = torch.nn.MSELoss()
    tl = DataLoader(DS(dtr, cache), batch_size=bs, shuffle=True, collate_fn=collate)
    vl = DataLoader(DS(dva, cache), batch_size=32, collate_fn=collate)
    best, best_state, bad = -1e9, None, 0
    for ep in range(epochs):
        model.train()
        for abx, abm, agx, agm, y in tl:
            opt.zero_grad()
            p = model(abx.to(DEVICE), abm.to(DEVICE), agx.to(DEVICE), agm.to(DEVICE))
            loss = lossf(p, y.to(DEVICE)); loss.backward(); opt.step()
        P, T = evaluate(model, vl); r = pearsonr(T, P)[0]
        if r > best:
            best, best_state, bad = r, {k: v.clone() for k, v in model.state_dict().items()}, 0
        else:
            bad += 1
            if bad >= patience:
                break
    if best_state:
        model.load_state_dict(best_state)
    return model


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default='DuaDeep-SeqAffinity', choices=list(MODELS))
    ap.add_argument('--split', default='random', choices=['random', 'cold'])
    ap.add_argument('--folds', type=int, default=10)
    ap.add_argument('--max_n', type=int, default=0)
    ap.add_argument('--epochs', type=int, default=40)
    ap.add_argument('--reduced', action='store_true', default=True)
    ap.add_argument('--full', dest='reduced', action='store_false')
    ap.add_argument('--dmodel', type=int, default=256)
    ap.add_argument('--maxlen', type=int, default=512)
    ap.add_argument('--bs', type=int, default=32)
    ap.add_argument('--out', default='')
    a = ap.parse_args()
    global MAXLEN; MAXLEN = a.maxlen

    df = pd.read_csv(CSV)
    cache = pickle.load(open(CACHE, 'rb'))
    df = df[df['H_seq'].isin(cache) & df['Ag_seq'].isin(cache)].reset_index(drop=True)
    df = df.dropna(subset=['pKD']).reset_index(drop=True)
    if a.max_n:
        df = df.sample(a.max_n, random_state=0).reset_index(drop=True)
    print(f"{a.model} | split={a.split} | n={len(df)} | device={DEVICE} | reduced={a.reduced}")

    splits = get_fold_splits(df, 10, 42, a.split)[:a.folds]
    ModelCls = MODELS[a.model]
    out = a.out or os.path.join(os.path.dirname(__file__), f"results_{a.model}_{a.split}.csv")
    rows = []; t0 = time.time()
    for fi, (tr, va) in enumerate(splits):
        ft = time.time()
        model = train_fold(ModelCls, df.iloc[tr], df.iloc[va], cache, a.reduced, a.epochs,
                           d_model=a.dmodel, bs=a.bs)
        P, T = evaluate(model, DataLoader(DS(df.iloc[va], cache), batch_size=32, collate_fn=collate))
        r, rho, rmse = metrics(P, T)
        rows.append({'fold': fi, 'pearson': r, 'spearman': rho, 'rmse': rmse, 'n': len(va)})
        print(f"  fold {fi}: r={r:.3f} rho={rho:.3f} rmse={rmse:.3f}  ({time.time()-ft:.0f}s)", flush=True)
        pd.DataFrame(rows).to_csv(out, index=False)          # incremental save
    res = pd.DataFrame(rows)
    print(f"\n{a.model} [{a.split}]  Pearson {res.pearson.mean():.3f}±{res.pearson.std():.3f} | "
          f"Spearman {res.spearman.mean():.3f} | RMSE {res.rmse.mean():.3f}  (total {time.time()-t0:.0f}s)")
    print("saved", out)


if __name__ == '__main__':
    main()
