"""
render_figures.py
=================
Compute the scene + a force-directed label layout for the ground-truth binding-residue
figures, and WRITE one reproducible .pml per (PDB, model) in both labelled and no-label
variants. Labels sit ON their residue; only labels that would collide are nudged apart and
linked back to the residue with a dotted leader line.

Geometry (camera projection) is computed in a headless pymol2 session; the actual PNGs are
produced by running the generated .pml files through the `pymol -cq` launcher (which does not
stamp the evaluation watermark).

Run:  python pymol_groundtruth/render_figures.py   # writes pml/*.pml
      then: for each pml -> pymol -cq pml/<f>.pml   # writes figures/*.png + *.pse
"""
import os, json, glob
import numpy as np
import matplotlib.cm as cm
import pymol2

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.join(HERE, 'pymol_groundtruth')
PMLDIR = os.path.join(ROOT, 'pml'); os.makedirs(PMLDIR, exist_ok=True)
IMG_REL = 'pymol_groundtruth/figures'
PDB_REL = 'pymol_groundtruth/pdb'
os.makedirs(os.path.join(ROOT, 'figures'), exist_ok=True)
PDBFILE = {'1VFB': 'pdb1vfb.ent', '4ETQ': 'pdb4etq.ent', '5GRJ': 'pdb5grj.ent', '5Y9J': 'pdb5y9j.ent'}

ROLE = {
    'heavy':   {'cmap': cm.get_cmap('Greens'),  'base': 0.26},
    'light':   {'cmap': cm.get_cmap('Blues'),   'base': 0.22},
    'antigen': {'cmap': cm.get_cmap('Purples'), 'base': 0.30},
}

def hexcol(cmap, t):
    c = cmap(float(np.clip(t, 0, 1)))
    return "0x{:02X}{:02X}{:02X}".format(int(c[0]*255), int(c[1]*255), int(c[2]*255))

def attr_to_t(v):
    return 0.42 + 0.58 * float(np.clip(v, 0, 1))

def static_scene(pdb, model, mp):
    """Command strings that build the cartoon + coloured GT sticks (no labels yet)."""
    attr_key = 'attr_three' if model == 'ours' else 'attr_two'
    roles = mp['chains']; obj = pdb
    keep = "+".join(roles[r] for r in ('heavy', 'light', 'antigen'))
    L = [
        f"# {pdb} -- {'Ours (All-CDR)' if model=='ours' else 'Two-stream'} -- GT binding residues coloured by IG attribution",
        "reinitialize",
        f"load {PDB_REL}/{PDBFILE[pdb]}, {obj}",
        "remove not polymer",
        f"remove {obj} and not chain {keep}",
        "hide everything",
        "bg_color white",
        "set ray_opaque_background, 1",
        "set ray_shadows, 0",
        "set ray_trace_mode, 1",
        "set ray_trace_color, grey50",
        "set antialias, 2",
        "set cartoon_fancy_helices, 1",
        "set cartoon_transparency, 0.25",
        "set cartoon_side_chain_helper, 1",
        "set stick_radius, 0.30",
        "set sphere_scale, 0.34",
        "set valence, 0",
        "set label_size, 15",
        "set label_font_id, 7",
        "set label_color, black",
        "set label_outline_color, white",
        "set label_bg_color, white",
        "set label_bg_transparency, 0.25",
        "set float_labels, 1",
        "set dash_color, grey30",
        "set dash_width, 1.6",
        "set dash_gap, 0.40",
        "set dash_length, 0.40",
        "set dash_radius, 0.045",
        "set dash_round_ends, 0",
        "show cartoon",
    ]
    for r in ('heavy', 'light', 'antigen'):
        L.append(f"color {hexcol(ROLE[r]['cmap'], ROLE[r]['base'])}, chain {roles[r]}")
    interface_sels = []
    for r in ('heavy', 'light', 'antigen'):
        ch = roles[r]; cmap = ROLE[r]['cmap']; nums = []
        for res in mp['residues'][r]:
            n = res['resnum']; v = res[attr_key] or 0.0
            nums.append(n)
            sel = f"{obj} and chain {ch} and resi {n}"
            L.append(f"show sticks, ({sel}) and (sidechain or name CA)")
            L.append(f"show spheres, ({sel}) and name CA")
            L.append(f"color {hexcol(cmap, attr_to_t(v))}, {sel}")
        if nums:
            interface_sels.append(f"chain {ch} and resi " + "+".join(str(n) for n in nums))
    L.append(f"select interface, {obj} and ((" + ") or (".join(interface_sels) + "))")
    L.append("deselect")
    L.append("set light_count, 2")
    L.append("set specular, 0.15")
    L.append("set ambient, 0.45")
    return L, obj, roles, attr_key

def declutter(labels, R, origin, view, fov, img_h=2100, img_w=2800):
    """Hybrid 2D label placement in camera space (verified convention:
    cam = R.T @ (model-origin);  model = R @ cam + origin).
    Isolated labels stay ON their residue; crowded labels are moved to left/right edge gutters
    and linked back with a dotted leader (conn=True). Keeps the structure prominent."""
    cam_dist = abs(float(view[11])) or 100.0
    vis_h = 2.0 * cam_dist * np.tan(np.radians(fov / 2.0))   # visible height (A) at the centre
    vis_w = vis_h * img_w / img_h
    apx = vis_h / img_h                                      # angstrom per pixel
    # rendered label footprint in pixels in a 2800x2100 ray (label_size 15, font_id 7) -- measured
    CHAR_W_PX, LINE_H_PX, PAD_PX = 47.5, 48.0, 16.0
    for d in labels:
        d['cam'] = R.T.dot(d['xyz'] - origin)
        d['anchor'] = d['cam'][:2].astype(float).copy()
        d['hw'] = (0.5 * len(d['text']) * CHAR_W_PX + PAD_PX) * apx
        d['hh'] = (0.5 * LINE_H_PX + PAD_PX) * apx

    def box(p, d):
        return (p[0] - d['hw'], p[0] + d['hw'], p[1] - d['hh'], p[1] + d['hh'])

    placed = []
    def hits(b):
        for q in placed:
            if not (b[1] < q[0] or b[0] > q[1] or b[3] < q[2] or b[2] > q[3]):
                return True
        return False

    # process from least-crowded to most-crowded so isolated labels claim their residue first
    for d in labels:
        d['nbr'] = sum(1 for e in labels if e is not d
                       and abs(d['anchor'][0] - e['anchor'][0]) < (d['hw'] + e['hw'])
                       and abs(d['anchor'][1] - e['anchor'][1]) < (d['hh'] + e['hh']))
    order = sorted(labels, key=lambda d: d['nbr'])
    overflow = []
    for d in order:
        b = box(d['anchor'], d)
        if not hits(b):
            d['pos'] = d['anchor'].copy(); d['conn'] = False; placed.append(b)
        else:
            overflow.append(d)

    # edge gutters for the overflow (kept fully inside the frame)
    fx = 0.40 * vis_w
    fy = 0.43 * vis_h
    for sign in (-1, +1):
        grp = [d for d in overflow if (d['anchor'][0] < 0) == (sign < 0)]
        grp.sort(key=lambda d: d['anchor'][1], reverse=True)
        k = len(grp)
        if k == 0:
            continue
        ys = [0.0] if k == 1 else list(np.linspace(fy, -fy, k))
        gx = sign * fx
        for d, yy in zip(grp, ys):
            d['pos'] = np.array([gx, yy]); d['conn'] = True; placed.append(box(d['pos'], d))

    for d in labels:
        cam_lab = np.array([d['pos'][0], d['pos'][1], d['cam'][2]])
        d['lpos'] = R.dot(cam_lab) + origin
    return labels

def collect_labels(cmd, obj, roles, mp):
    labels = []
    for r in ('heavy', 'light', 'antigen'):
        ch = roles[r]
        for res in mp['residues'][r]:
            n = res['resnum']
            casel = f"{obj} and chain {ch} and resi {n} and name CA"
            try:
                xyz = np.array(cmd.get_atom_coords(casel))
            except Exception:
                continue
            sp = []
            cmd.iterate(casel, "rn.append(resn)", space={'rn': sp})
            labels.append({'ca_sel': casel, 'text': f"{sp[0] if sp else ''}{n}",
                           'xyz': xyz, 'tag': f"{ch}{n}"})
    return labels

def set_view_block(view):
    return "set_view (\\\n" + ",\\\n".join(
        "   " + ", ".join(f"{view[i+j]:.7f}" for j in range(3)) for i in range(0, 18, 3)) + " )"

def write_pml(path, lines, dyn, set_view, png_rel, pse_rel):
    render = ['set ray_trace_fog, 0', 'ray 2800, 2100',
              f'png {png_rel}, dpi=350, ray=0', f'save {pse_rel}']
    pml = lines + [""] + (dyn + [""] if dyn else []) + [set_view, ""] + render
    open(path, 'w', encoding='utf-8').write("\n".join(pml) + "\n")

def build(cmd, pdb, model, mp):
    lines, obj, roles, attr_key = static_scene(pdb, model, mp)
    cmd.reinitialize()
    for ln in lines:
        if not ln.startswith('#'):
            cmd.do(ln)
    cmd.orient('interface')
    cmd.zoom('interface', 9)            # interface centred, structure prominent
    view = cmd.get_view()
    R = np.array(view[0:9]).reshape(3, 3)
    origin = np.array(view[12:15])
    fov = float(cmd.get('field_of_view'))
    set_view = set_view_block(view)
    labels = collect_labels(cmd, obj, roles, mp)
    labels = declutter(labels, R, origin, view, fov)
    # dynamic label commands
    dyn = []
    for d in labels:
        x, y, z = [round(float(v), 3) for v in d['lpos']]
        pa = f"lab_{d['tag']}"
        dyn.append(f'pseudoatom {pa}, pos=[{x}, {y}, {z}], label="{d["text"]}"')
        if d['conn']:
            dl = f"dl_{d['tag']}"
            dyn.append(f"distance {dl}, ({d['ca_sel']}), {pa}")
            dyn.append(f"hide labels, {dl}")
    dyn.append("set label_color, black, lab_*")
    # write both variants (same camera)
    write_pml(os.path.join(PMLDIR, f"{pdb}_{model}.pml"), lines, dyn, set_view,
              f"{IMG_REL}/{pdb}_{model}.png", f"{IMG_REL}/{pdb}_{model}.pse")
    write_pml(os.path.join(PMLDIR, f"{pdb}_{model}_nolabel.pml"), lines, [], set_view,
              f"{IMG_REL}/{pdb}_{model}_nolabel.png", f"{IMG_REL}/{pdb}_{model}_nolabel.pse")
    nconn = sum(1 for d in labels if d['conn'])
    print(f"WROTE pml {pdb} {model}  ({len(labels)} labels, {nconn} leader lines)")

def main():
    with pymol2.PyMOL() as p:
        cmd = p.cmd
        cmd.cd(HERE)
        for f in sorted(glob.glob(os.path.join(ROOT, 'mapping_*.json'))):
            mp = json.load(open(f))
            for model in ('ours', 'twostream'):
                build(cmd, mp['pdb'], model, mp)

if __name__ == '__main__':
    main()
