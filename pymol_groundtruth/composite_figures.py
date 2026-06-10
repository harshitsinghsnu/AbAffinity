"""
composite_figures.py
====================
Publication-style composite interface figures (in the style of the reference image):
  * central overview: antibody (salmon) + antigen (green), transparent surface + cartoon,
    dashed circles marking the interface sub-regions;
  * 3 zoom-in inset boxes (connected to the circles by black lines) showing the ground-truth
    interface residues as sticks, hydrogen bonds as yellow dashes, and RES-NUM labels;
  * residue sticks are shaded by THAT model's IG attribution within the partner colour
    (antibody = Reds ramp, antigen = Greens ramp), so ours vs two-stream still differ.

Phase 1 (pymol2): build the scene, cluster the GT residues into inset groups, detect H-bonds,
compute the camera + label layout, and write .pml scripts (overview + insets).
Phase 2: render the .pml via the `pymol -cq` launcher (no watermark).
Phase 3 (matplotlib): assemble overview + insets with connector lines + dashed circles.

Run:  python pymol_groundtruth/composite_figures.py
"""
import os, json, glob, subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Ellipse
from scipy.cluster.vq import kmeans2
import pymol2

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.join(HERE, 'pymol_groundtruth')
PMLDIR = os.path.join(ROOT, 'pml_composite'); os.makedirs(PMLDIR, exist_ok=True)
IMGDIR = os.path.join(ROOT, 'figures_composite'); os.makedirs(IMGDIR, exist_ok=True)
IMG_REL = 'pymol_groundtruth/figures_composite'
PDB_REL = 'pymol_groundtruth/pdb'
PDBFILE = {'1VFB': 'pdb1vfb.ent', '4ETQ': 'pdb4etq.ent', '5GRJ': 'pdb5grj.ent', '5Y9J': 'pdb5y9j.ent'}

# same chain-role colour scheme as pymol_groundtruth/figures (Heavy=Greens, Light=Blues, Antigen=Purples)
ROLE = {
    'heavy':   {'cmap': cm.get_cmap('Greens'),  'base': 0.26},
    'light':   {'cmap': cm.get_cmap('Blues'),   'base': 0.22},
    'antigen': {'cmap': cm.get_cmap('Purples'), 'base': 0.30},
}
OV_W, OV_H = 2400, 2000            # overview ray size
IN_W, IN_H = 1300, 1000            # inset ray size

def hexcol(cmap, t):
    c = cmap(float(np.clip(t, 0, 1)))
    return "0x{:02X}{:02X}{:02X}".format(int(c[0]*255), int(c[1]*255), int(c[2]*255))

def attr_t(v):
    return 0.42 + 0.58 * float(np.clip(v or 0.0, 0, 1))   # matches figures/ stick shading

# ---------------------------------------------------------------- scene
def scene_cmds(pdb, mp, model, surface):
    """surface=True -> overview (transparent surface + cartoon + all GT sticks).
       surface=False -> inset    (no surface, faint cartoon, sticks shown later per cluster)."""
    attr_key = 'attr_three' if model == 'ours' else 'attr_two'
    roles = mp['chains']; obj = pdb
    keep = "+".join(roles[r] for r in ('heavy', 'light', 'antigen'))
    ab_ch = [roles['heavy'], roles['light']]; ag_ch = roles['antigen']
    L = [
        "reinitialize",
        f"load {PDB_REL}/{PDBFILE[pdb]}, {obj}",
        "remove not polymer",
        f"remove {obj} and not chain {keep}",
        "hide everything",
        "bg_color white",
        "set ray_opaque_background, 1",
        "set ray_shadows, 0",
        "set antialias, 2",
        "set cartoon_side_chain_helper, 1",
        "set stick_radius, %s" % ("0.26" if surface else "0.32"),
        "set sphere_scale, 0.30",
        "set valence, 0",
        "set two_sided_lighting, 1",
        "set label_size, %d" % (17 if surface else 24),
        "set label_font_id, 7",
        "set label_color, black",
        "set label_outline_color, white",
        "set label_bg_color, white",
        "set label_bg_transparency, 0.2",
        "set float_labels, 1",
        "set dash_color, yellow",
        "set dash_gap, 0.30",
        "set dash_width, %s" % ("3.0" if surface else "4.0"),
        "set dash_radius, %s" % ("0.06" if surface else "0.09"),
        "set cartoon_transparency, %s" % ("0.0" if surface else "0.72"),
        "set ray_trace_mode, %s" % ("1" if surface else "0"),
        "set ray_trace_color, grey60",
        "show cartoon",
    ]
    # cartoon coloured by chain role (Heavy=green, Light=blue, Antigen=purple), matching figures/
    for r in ('heavy', 'light', 'antigen'):
        L.append(f"color {hexcol(ROLE[r]['cmap'], ROLE[r]['base'])}, chain {roles[r]}")
    # overview is cartoon-only (no surface), per request
    # attribution-shaded GT sticks (shown now for overview; for insets shown later per cluster)
    for r in ('heavy', 'light', 'antigen'):
        ch = roles[r]; cmap = ROLE[r]['cmap']
        for res in mp['residues'][r]:
            n = res['resnum']; sel = f"{obj} and chain {ch} and resi {n}"
            if surface:
                L.append(f"show sticks, ({sel}) and (sidechain or name CA)")
            L.append(f"color {hexcol(cmap, attr_t(res[attr_key]))}, {sel}")
    return L, obj, roles

# ---------------------------------------------------------------- geometry helpers
def gt_atoms(cmd, obj, roles, mp):
    """Return list of dicts per GT residue: chain, resi, role, ca(np3), text, attr."""
    out = []
    for r in ('heavy', 'light', 'antigen'):
        ch = roles[r]
        for res in mp['residues'][r]:
            n = res['resnum']; casel = f"{obj} and chain {ch} and resi {n} and name CA"
            try:
                ca = np.array(cmd.get_atom_coords(casel))
            except Exception:
                continue
            sp = []
            cmd.iterate(casel, "rn.append(resn)", space={'rn': sp})
            out.append({'ch': ch, 'resi': n, 'role': r, 'ca': ca,
                        'text': f"{sp[0] if sp else ''}-{n}", 'tag': f"{ch}{n}",
                        'sel': f"{obj} and chain {ch} and resi {n}", 'ca_sel': casel})
    return out

def hbonds(cmd, obj, gts):
    """Polar-atom (N/O) pairs < 3.5 A between different GT residues -> list of (sel1, sel2)."""
    # gather polar atoms of GT residues
    atoms = []
    gtsel = " or ".join(f"(chain {g['ch']} and resi {g['resi']})" for g in gts)
    space = {'A': []}
    cmd.iterate_state(1, f"({obj}) and ({gtsel}) and (elem N+O)",
                      "A.append((chain, resi, name, x, y, z))", space=space)
    A = space['A']
    pairs = []
    seen = set()
    for i in range(len(A)):
        ci, ri, ni, xi, yi, zi = A[i]
        for j in range(i + 1, len(A)):
            cj, rj, nj, xj, yj, zj = A[j]
            if (ci, ri) == (cj, rj):
                continue
            d = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
            if d < 3.5**2:
                key = tuple(sorted([(ci, ri), (cj, rj)]))
                if key in seen:
                    continue
                seen.add(key)
                pairs.append((f"{obj} and chain {ci} and resi {ri} and name {ni}",
                              f"{obj} and chain {cj} and resi {rj} and name {nj}",
                              (ci, ri), (cj, rj)))
    return pairs

def set_view_block(view):
    return "set_view (\\\n" + ",\\\n".join(
        "   " + ", ".join(f"{view[i+j]:.7f}" for j in range(3)) for i in range(0, 18, 3)) + " )"

def project_px(p3, R, origin, view, fov, w, h):
    cam = R.T.dot(np.asarray(p3) - origin)
    vis_h = 2.0 * abs(float(view[11])) * np.tan(np.radians(fov / 2.0))
    apx = vis_h / h
    px = w / 2.0 + cam[0] / apx
    py = h / 2.0 - cam[1] / apx
    return px, py

# ---------------------------------------------------------------- main build
def build(cmd, pdb, model, mp):
    ov_lines, obj, roles = scene_cmds(pdb, mp, model, surface=True)
    in_lines, _, _ = scene_cmds(pdb, mp, model, surface=False)
    cmd.reinitialize()
    for ln in ov_lines:
        cmd.do(ln)
    gts = gt_atoms(cmd, obj, roles, mp)
    hb = hbonds(cmd, obj, gts)
    hb_residues = set()
    for _, _, a, b in hb:
        hb_residues.add(a); hb_residues.add(b)
    # cluster GT residues into <=3 inset groups by CA position
    cas = np.array([g['ca'] for g in gts])
    k = min(3, len(gts))
    cen, lab = kmeans2(cas, k, seed=0, minit='++')
    for g, c in zip(gts, lab):
        g['clu'] = int(c)
    res2clu = {(g['ch'], str(g['resi'])): g['clu'] for g in gts}
    # H-bond distance commands; tag each pair with the cluster if both ends share one
    hb_cmds = []          # all (overview)
    hb_by_clu = {c: [] for c in range(k)}
    for idx, (s1, s2, a, b) in enumerate(hb):
        nm = f"hb{idx}"
        cmds = [f"distance {nm}, ({s1}), ({s2})", f"hide labels, {nm}", f"set dash_color, yellow, {nm}"]
        hb_cmds += cmds
        ca, cb = res2clu.get(a), res2clu.get(b)
        if ca is not None and ca == cb:
            hb_by_clu[ca] += cmds

    # ---- overview ----
    cmd.orient(obj)
    cmd.zoom(obj, 3)
    ov_view = cmd.get_view()
    R = np.array(ov_view[0:9]).reshape(3, 3); origin = np.array(ov_view[12:15])
    fov = float(cmd.get('field_of_view'))
    # cluster centroids -> overview pixels (for circles)
    circles = []
    for c in range(k):
        pts = cas[lab == c]
        cen3 = pts.mean(axis=0)
        circles.append(project_px(cen3, R, origin, ov_view, fov, OV_W, OV_H))
    ov_dyn = hb_cmds + ["set label_color, black"]
    write_pml(os.path.join(PMLDIR, f"{pdb}_{model}_overview.pml"),
              ov_lines, ov_dyn, set_view_block(ov_view),
              f"{IMG_REL}/{pdb}_{model}_overview.png", OV_W, OV_H)

    # ---- insets ----
    inset_files = []
    for c in range(k):
        grp = [g for g in gts if g['clu'] == c]
        gsel = " or ".join(f"(chain {g['ch']} and resi {g['resi']})" for g in grp)
        cmd.orient(f"{obj} and ({gsel})")
        cmd.zoom(f"{obj} and ({gsel})", 1.6)
        iv = cmd.get_view()
        Ri = np.array(iv[0:9]).reshape(3, 3); oi = np.array(iv[12:15]); fovi = float(cmd.get('field_of_view'))
        # labels for this cluster (declutter inside inset)
        lab_objs = [{'xyz': g['ca'], 'text': g['text'], 'ca_sel': g['ca_sel'], 'tag': g['tag']} for g in grp]
        lab_objs = declutter(lab_objs, Ri, oi, iv, fovi, IN_H, IN_W)
        # inset uses the no-surface scene; show only this cluster's sticks + its H-bonds
        dyn = list(hb_by_clu[c])
        dyn.append("hide labels")
        for g in grp:
            dyn.append(f"show sticks, ({g['sel']}) and (sidechain or name CA)")
            dyn.append(f"show spheres, ({g['sel']}) and name CA")
            dyn.append(f"set sphere_scale, 0.20, ({g['sel']})")
        for d in lab_objs:
            x, y, z = [round(float(v), 3) for v in d['lpos']]
            pa = f"lab_{d['tag']}"
            dyn.append(f'pseudoatom {pa}, pos=[{x}, {y}, {z}], label="{d["text"]}"')
            if d['conn']:
                dl = f"dl_{d['tag']}"
                dyn.append(f"distance {dl}, ({d['ca_sel']}), {pa}")
                dyn.append(f"hide labels, {dl}")
                dyn.append(f"set dash_color, grey50, {dl}")
        dyn.append("set label_color, black, lab_*")
        write_pml(os.path.join(PMLDIR, f"{pdb}_{model}_inset{c}.pml"),
                  in_lines, dyn, set_view_block(iv),
                  f"{IMG_REL}/{pdb}_{model}_inset{c}.png", IN_W, IN_H)
        inset_files.append(f"{pdb}_{model}_inset{c}.png")
    return {'pdb': pdb, 'model': model, 'k': k, 'circles': circles,
            'overview': f"{pdb}_{model}_overview.png", 'insets': inset_files}

def write_pml(path, lines, dyn, set_view, png_rel, w, h):
    render = ['set ray_trace_fog, 0', f'ray {w}, {h}', f'png {png_rel}, dpi=300, ray=0']
    pml = lines + [""] + dyn + ["", set_view, ""] + render
    open(path, 'w', encoding='utf-8').write("\n".join(pml) + "\n")

# reuse the validated declutter from render_figures
import importlib.util as _u
_spec = _u.spec_from_file_location('rf', os.path.join(ROOT, 'render_figures.py'))
_rf = _u.module_from_spec(_spec); _spec.loader.exec_module(_rf)
declutter = _rf.declutter

# ---------------------------------------------------------------- compositing
def composite(meta):
    pdb, model = meta['pdb'], meta['model']
    ov = plt.imread(os.path.join(IMGDIR, meta['overview']))
    insets = [plt.imread(os.path.join(IMGDIR, f)) for f in meta['insets']]
    k = meta['k']
    fig = plt.figure(figsize=(13, 8.5), dpi=300)
    # overview centred
    axo = fig.add_axes([0.26, 0.10, 0.48, 0.80]); axo.imshow(ov); axo.axis('off')
    H, W = ov.shape[0], ov.shape[1]
    # inset corner slots (fig coords): TL, TR, BR (image-style)
    slots = [[0.005, 0.55, 0.255, 0.40], [0.74, 0.55, 0.255, 0.40], [0.74, 0.05, 0.255, 0.40]]
    for i in range(k):
        sx, sy, sw, sh = slots[i]
        axi = fig.add_axes([sx, sy, sw, sh]); axi.imshow(insets[i])
        axi.set_xticks([]); axi.set_yticks([])
        for s in axi.spines.values():
            s.set_color('black'); s.set_linewidth(1.6)
        # dashed circle on overview at cluster centroid
        cx, cy = meta['circles'][i]
        axo.add_patch(Ellipse((cx, cy), W*0.085, H*0.085, fill=False, ls=(0, (4, 3)),
                              lw=2.2, ec='black', zorder=5))
        # connector line from circle to inset box (figure coords)
        x_fig = 0.26 + 0.48 * (cx / W)
        y_fig = 0.10 + 0.80 * (1 - cy / H)
        box_cx = sx + (sw if sx < 0.5 else 0.0)   # inner edge of inset
        box_cy = sy + sh / 2.0
        line = plt.Line2D([box_cx, x_fig], [box_cy, y_fig], transform=fig.transFigure,
                          color='black', lw=1.6, zorder=1)
        fig.add_artist(line)
    out = os.path.join(IMGDIR, f"{pdb}_{model}_composite")
    fig.savefig(out + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(out + '.pdf', bbox_inches='tight')
    plt.close(fig)
    print("composited", os.path.basename(out) + '.png')

# ---------------------------------------------------------------- driver
def main():
    metas = []
    with pymol2.PyMOL() as p:
        cmd = p.cmd; cmd.cd(HERE)
        for f in sorted(glob.glob(os.path.join(ROOT, 'mapping_*.json'))):
            mp = json.load(open(f))
            for model in ('ours', 'twostream'):
                metas.append(build(cmd, mp['pdb'], model, mp))
                print("wrote pml", mp['pdb'], model)
    json.dump(metas, open(os.path.join(ROOT, '_composite_meta.json'), 'w'))
    # render all pml via launcher
    for pml in sorted(glob.glob(os.path.join(PMLDIR, '*.pml'))):
        subprocess.run(['pymol', '-cq', os.path.relpath(pml, HERE)], cwd=HERE,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for meta in metas:
        composite(meta)

if __name__ == '__main__':
    main()
