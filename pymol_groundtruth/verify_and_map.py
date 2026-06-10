"""
verify_and_map.py
=================
1. Parse ground-truth residue lists (paper-reported binding residues).
2. Verify each against the deposited PDB author numbering (chain + resi -> resname).
3. Align the ESM-2 embedding sequence to the modelled structure sequence so each
   structure residue (author resnum) gets its IG-attribution value.
4. Emit pymol_groundtruth/mapping_{PDB}.json with, per chain/model, the GT residues
   and their signed-heatmap-normalised attribution (0..1).
"""
import os, re, json, warnings, sys
import numpy as np
warnings.filterwarnings('ignore')
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def _load_complexes():
    """Extract COMPLEXES dict literal from integrated_gradient_twostream.py without importing torch/captum."""
    src = open(os.path.join(HERE, 'integrated_gradient_twostream.py'), encoding='utf-8').read()
    i = src.index('COMPLEXES = {')
    j = src.index('\n}', i) + 2
    ns = {}
    exec(src[i:j], ns)
    return ns['COMPLEXES']
PDBDIR = os.path.join(HERE, 'pymol_groundtruth', 'pdb')
EXP = os.path.join(HERE, 'experiments', 'results_explainability')
SEQ = _load_complexes()

# chain roles (from structure inspection)
CHAINS = {
    '1VFB': {'heavy': 'B', 'light': 'A', 'antigen': 'C'},
    '4ETQ': {'heavy': 'H', 'light': 'L', 'antigen': 'C'},  # H/L bind antigen copy C (X binds the A/B copy)
    '5GRJ': {'heavy': 'H', 'light': 'L', 'antigen': 'A'},
    '5Y9J': {'heavy': 'H', 'light': 'L', 'antigen': 'A'},
}
PDBFILE = {'1VFB': 'pdb1vfb.ent', '4ETQ': 'pdb4etq.ent', '5GRJ': 'pdb5grj.ent', '5Y9J': 'pdb5y9j.ent'}

THREE2ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
NAME3 = {v:k.capitalize() for k,v in THREE2ONE.items()}  # 'T'->'Thr'
AA3SET = set(THREE2ONE)

# ---- ground-truth residues (verbatim from user) ----
GT_RAW = {
 '1VFB': {
   'heavy':'Thr30,Gly31,Tyr32,Trp52,Gly53,Asp54,Asn56,Asp58,Arg99,Asp100,Tyr101,Arg102',
   'light':'His30,Tyr32,Tyr49,Tyr50,Thr52,Tyr53,Phe91,Trp92,Ser93,Thr94',
   'antigen':'Asp18,Asp19,Gly22,Ser24,Gly26,Gly102,Lys116,Gln121,Leu129'},
 '4ETQ': {
   'heavy':'Ser29,Asn31,Phe32,Trp34,Met51,Asp53,Ser55,Glu56,Glu58,Arg60,Arg75,Tyr102,Arg103,Asp105',
   'light':'Ser31,Trp91,Thr92,Tyr94',
   'antigen':'Gln3,Leu5,Thr39,Gly40,Lys41,Arg44,Lys108,Asn145,Ile174,Asn175,His176,Ser177,Ser204,Leu205,Ile215,Glu217,Tyr219,Arg220,Asn221,Tyr223,Lys224'},
 '5GRJ': {
   'heavy':'F27,T28,S31,I33,Y52,P53,S54,I57,F59,L101,G102,T103,V104,T105,T106',
   'light':'Y32,Y34,Y93,S97,R99,S95',
   'antigen':'I54,E60,D61,K62,N63,Q66,Y56,E58,M59,H69,D73,K75,R113,A121,Y123,R125,M115,S117,H78'},
 '5Y9J': {
   'heavy':'N31,M54,F55,T57,D101,L102,L103,L104',
   'light':'N95,Y31,R28,K50',
   'antigen':'S162,Y163,Y206,R231,I233,L240,N242,R265,E266,L224,D222,S225'},
}

def parse_res(token):
    """'Thr30' or 'T28' -> (one_letter, resnum)."""
    t = token.strip()
    m = re.match(r'^([A-Za-z]{1,3})(\d+)$', t)
    aa, num = m.group(1), int(m.group(2))
    if len(aa) == 3:
        one = THREE2ONE[aa.upper()]
    else:
        one = aa.upper()
    return one, num

def chain_residues(chain):
    out = []  # (resnum, icode, one_letter)
    for r in chain:
        if r.id[0] != ' ':
            continue
        rn = r.resname.strip().upper()
        if rn not in THREE2ONE:
            continue
        out.append((r.id[1], r.id[2], THREE2ONE[rn]))
    return out

def build_map(emb_seq, struct_res):
    """Align embedding seq to structure seq; return dict author_resnum -> embedding_index."""
    sseq = ''.join(o for _,_,o in struct_res)
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    aln = aligner.align(emb_seq, sseq)[0]
    # aln.aligned -> blocks [[(e0,e1),...],[(s0,s1),...]]
    a_emb, a_str = aln.aligned
    res = {}
    for (e0, e1), (s0, s1) in zip(a_emb, a_str):
        for k in range(e1 - e0):
            ei = e0 + k
            si = s0 + k
            resnum = struct_res[si][0]
            res[resnum] = ei
    return res

def main():
    parser = PDBParser(QUIET=True)
    report = []
    for pdb, roles in CHAINS.items():
        st = parser.get_structure(pdb, os.path.join(PDBDIR, PDBFILE[pdb]))[0]
        d = np.load(os.path.join(EXP, f'{pdb}.npz'))
        out = {'pdb': pdb, 'chains': roles, 'residues': {}}
        for role in ('heavy', 'light', 'antigen'):
            ch_id = roles[role]
            chain = st[ch_id]
            sres = chain_residues(chain)
            sres_by_num = {n: o for n, _, o in sres}
            emb_seq = SEQ[pdb][role]
            amap = build_map(emb_seq, sres)
            norms = {}
            for mod in ('three', 'two'):
                a = np.abs(d[f'{mod}_{role}'].astype(float)) * 0  # placeholder
                raw = d[f'{mod}_{role}'].astype(float)
                lo, hi = raw.min(), raw.max()
                norms[mod] = (raw - lo) / (hi - lo + 1e-8)  # signed heatmap norm 0..1
            reslist = []
            for tok in GT_RAW[pdb][role].split(','):
                one, num = parse_res(tok)
                struct_one = sres_by_num.get(num, '?')
                ok = (struct_one == one)
                ei = amap.get(num, None)
                v3 = float(norms['three'][ei]) if (ei is not None and ei < len(norms['three'])) else None
                v2 = float(norms['two'][ei]) if (ei is not None and ei < len(norms['two'])) else None
                reslist.append({'token': tok, 'one': one, 'resnum': num,
                                'struct_resn': struct_one, 'match': ok,
                                'attr_three': v3, 'attr_two': v2})
            out['residues'][role] = reslist
            nbad = sum(1 for r in reslist if not r['match'])
            nunmap = sum(1 for r in reslist if r['attr_three'] is None)
            report.append(f'{pdb} {role:8s} chain {ch_id}: {len(reslist)} GT, mismatched={nbad}, unmapped={nunmap}')
            for r in reslist:
                if not r['match'] or r['attr_three'] is None:
                    report.append(f"    !! {r['token']} -> struct={r['struct_resn']}{r['resnum']} match={r['match']} attr3={r['attr_three']}")
        json.dump(out, open(os.path.join(HERE,'pymol_groundtruth',f'mapping_{pdb}.json'),'w'), indent=1)
    print('\n'.join(report))

if __name__ == '__main__':
    main()
