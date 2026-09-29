import pandas as pd
p = r"c:\Users\hs494\OneDrive - Shiv Nadar Institution of Eminence\Desktop\saaintdb_with_antigen_names.csv"
d = pd.read_csv(p)
print("rows:", len(d))
print("cols:", list(d.columns))
print("unique PDB_ID:", d['PDB_ID'].nunique() if 'PDB_ID' in d.columns else 'NO PDB_ID')
for c in ['H_seq', 'L_seq', 'Ag_seq', 'pKD']:
    print(f"  has {c}:", c in d.columns)
print(d.head(2).to_string())
