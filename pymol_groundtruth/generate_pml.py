"""
generate_pml.py
===============
Generate one PyMOL script (.pml) per (PDB, model) that:
  * loads the deposited structure, white background, publication ray-trace settings,
  * draws all three chains as pale cartoons (Heavy=green, Light=blue, Antigen=purple),
  * shows the paper ground-truth binding residues as sticks, each coloured by THAT
    model's IG attribution (signed heatmap normalisation) through the chain colormap
    (Greens / Blues / Purples) so residues with attribution > 0.4 appear saturated,
  * orients on the antibody-antigen interface and renders a 300-350 dpi PNG + .pse.

Reads pymol_groundtruth/mapping_{PDB}.json (from verify_and_map.py).
Writes pymol_groundtruth/pml/{PDB}_{ours|twostream}.pml
"""
import os, json, glob
import numpy as np
import matplotlib.cm as cm

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.join(HERE, 'pymol_groundtruth')
PMLDIR = os.path.join(ROOT, 'pml'); os.makedirs(PMLDIR, exist_ok=True)
IMGDIR = os.path.join(ROOT, 'figures'); os.makedirs(IMGDIR, exist_ok=True)
# paths inside the .pml are RELATIVE to the 3_stream dir (avoids spaces in the absolute path)
PDB_REL = 'pymol_groundtruth/pdb'
IMG_REL = 'pymol_groundtruth/figures'
PDBFILE = {'1VFB': 'pdb1vfb.ent', '4ETQ': 'pdb4etq.ent', '5GRJ': 'pdb5grj.ent', '5Y9J': 'pdb5y9j.ent'}

# chain colormap + pale cartoon base (matches explainability heatmaps: H=Greens,L=Blues,A=Purples)
# 'base' = light cartoon tint (kept BELOW the stick attribution range so sticks always stand out)
ROLE = {
    'heavy':   {'cmap': cm.get_cmap('Greens'),  'base': 0.26},
    'light':   {'cmap': cm.get_cmap('Blues'),   'base': 0.22},
    'antigen': {'cmap': cm.get_cmap('Purples'), 'base': 0.30},
}

def hexcol(cmap, t):
    c = cmap(float(np.clip(t, 0, 1)))
    return "0x{:02X}{:02X}{:02X}".format(int(c[0]*255), int(c[1]*255), int(c[2]*255))

def attr_to_t(v):
    """Map signed-heatmap-norm attribution (0..1) to a visible colormap position,
    kept clearly darker than the cartoon tint so highlighted residues pop; attr>0.4 saturates."""
    return 0.42 + 0.58 * float(np.clip(v, 0, 1))

def fwd(p):
    return p.replace('\\', '/')

def build(pdb, model, mapping):
    """model: 'ours' (-> attr_three) or 'twostream' (-> attr_two)"""
    attr_key = 'attr_three' if model == 'ours' else 'attr_two'
    roles = mapping['chains']
    obj = pdb
    L = []
    A = L.append
    A(f"# {pdb} -- {'Ours (All-CDR)' if model=='ours' else 'Two-stream'} -- ground-truth binding residues coloured by IG attribution")
    A("reinitialize")
    A(f'load {PDB_REL}/{PDBFILE[pdb]}, {obj}')
    A(f"remove not polymer")
    A(f"hide everything")
    # keep only the three relevant chains
    keep = "+".join(roles[r] for r in ('heavy', 'light', 'antigen'))
    A(f"remove {obj} and not chain {keep}")
    A(f"bg_color white")
    A(f"set ray_opaque_background, 1")
    A(f"set ray_shadows, 0")
    A(f"set ray_trace_mode, 1")          # black outlines (publication)
    A(f"set ray_trace_color, grey50")
    A(f"set antialias, 2")
    A(f"set cartoon_fancy_helices, 1")
    A(f"set cartoon_transparency, 0.30")
    A(f"set stick_radius, 0.30")
    A(f"set sphere_scale, 0.32")
    A(f"set valence, 0")
    # residue labels (e.g. ASP58): bold, black with white outline so they read on any chain colour
    A(f"set label_size, 13")
    A(f"set label_font_id, 7")
    A(f"set label_color, black")
    A(f"set label_outline_color, white")
    A(f"set label_bg_color, white")
    A(f"set label_bg_transparency, 0.30")
    A(f"set float_labels, 1")
    A(f"set label_position, (0, 0, 3)")
    A(f"show cartoon")
    # pale per-chain cartoon base colours (light but clearly distinguishable)
    for r in ('heavy', 'light', 'antigen'):
        A(f"color {hexcol(ROLE[r]['cmap'], ROLE[r]['base'])}, chain {roles[r]}")
    A("")
    # ground-truth residues as sticks, coloured by attribution
    interface_sels = []
    for r in ('heavy', 'light', 'antigen'):
        ch = roles[r]; cmap = ROLE[r]['cmap']
        nums = []
        for res in mapping['residues'][r]:
            n = res['resnum']; v = res[attr_key]
            if v is None:
                v = 0.0
            nums.append(n)
            sel = f"{obj} and chain {ch} and resi {n}"
            A(f"show sticks, ({sel}) and (sidechain or name CA)")
            A(f"show spheres, ({sel}) and name CA")
            A(f"color {hexcol(cmap, attr_to_t(v))}, {sel}")
            A(f'label ({sel}) and name CA, "%s%s" % (resn, resi)')
        if nums:
            numsel = "+".join(str(n) for n in nums)
            gsel = f"gt_{r}"
            A(f"select {gsel}, {obj} and chain {ch} and resi {numsel}")
            interface_sels.append(f"chain {ch} and resi {numsel}")
    A("")
    A(f"select interface, {obj} and ((" + ") or (".join(interface_sels) + "))")
    A("deselect")
    # lighting / view
    A("set light_count, 2")
    A("set specular, 0.15")
    A("set ambient, 0.45")
    A("set cartoon_side_chain_helper, 1")
    A("orient interface")        # focus on the binding interface so residue labels have room to read
    A("zoom interface, 6")
    A(f"set ray_trace_fog, 0")
    A(f"ray 2800, 2100")
    A(f'png {IMG_REL}/{pdb}_{model}.png, dpi=350, ray=0')
    A(f'save {IMG_REL}/{pdb}_{model}.pse')
    A(f'print "WROTE {pdb} {model}"')
    path = os.path.join(PMLDIR, f"{pdb}_{model}.pml")
    open(path, 'w', encoding='utf-8').write("\n".join(L) + "\n")
    return path

def main():
    paths = []
    for f in sorted(glob.glob(os.path.join(ROOT, 'mapping_*.json'))):
        mp = json.load(open(f))
        pdb = mp['pdb']
        for model in ('ours', 'twostream'):
            paths.append(build(pdb, model, mp))
    print("Wrote:")
    for p in paths:
        print("  ", os.path.relpath(p, HERE))

if __name__ == '__main__':
    main()
