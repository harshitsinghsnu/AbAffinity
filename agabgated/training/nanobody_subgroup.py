import os,pandas as pd,numpy as np
from scipy.stats import pearsonr,spearmanr
HERE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT=os.path.join(HERE,'experiments','results_allcdr_stats')
d=pd.read_csv(os.path.join(HERE,'results_saaintdb_allcdr/random/all_preds.csv'))
def valid(s): 
    return isinstance(s,str) and s.strip() not in ('','nan','N/A','NA','N.A.','None') and len(s.strip())>5
d['is_nano']=~d['L_seq'].apply(valid)
def met(g):
    t=g['binding_affinity'].values; p=g['predicted_affinity'].values
    return dict(n=len(g),pearson=float(pearsonr(t,p)[0]),spearman=float(spearmanr(t,p)[0]),rmse=float(np.sqrt(np.mean((t-p)**2))))
rows=[{'subgroup':'Overall',**met(d)},
      {'subgroup':'Antibody (paired H+L)',**met(d[~d.is_nano])},
      {'subgroup':'Nanobody (VHH)',**met(d[d.is_nano])}]
pd.DataFrame(rows).to_csv(os.path.join(OUT,'antibody_vs_nanobody.csv'),index=False)
print(pd.DataFrame(rows).to_string(index=False))
