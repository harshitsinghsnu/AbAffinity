reinitialize
load pymol_groundtruth/pdb/pdb5y9j.ent, 5Y9J
remove not polymer
remove 5Y9J and not chain H+L+A
hide everything
bg_color white
set ray_opaque_background, 1
set ray_shadows, 0
set antialias, 2
set cartoon_side_chain_helper, 1
set stick_radius, 0.32
set sphere_scale, 0.30
set valence, 0
set two_sided_lighting, 1
set label_size, 24
set label_font_id, 7
set label_color, black
set label_outline_color, white
set label_bg_color, white
set label_bg_transparency, 0.2
set float_labels, 1
set dash_color, yellow
set dash_gap, 0.30
set dash_width, 4.0
set dash_radius, 0.09
set cartoon_transparency, 0.72
set ray_trace_mode, 0
set ray_trace_color, grey60
show cartoon
color 0xC4E7BD, chain H
color 0xCBDEF0, chain L
color 0xCECEE5, chain A
color 0x157E3A, 5Y9J and chain H and resi 31
color 0x339C51, 5Y9J and chain H and resi 54
color 0x2E974E, 5Y9J and chain H and resi 55
color 0x278F48, 5Y9J and chain H and resi 57
color 0x76C477, 5Y9J and chain H and resi 101
color 0x6EC173, 5Y9J and chain H and resi 102
color 0x7AC67B, 5Y9J and chain H and resi 103
color 0x64BC6E, 5Y9J and chain H and resi 104
color 0x2777B8, 5Y9J and chain L and resi 95
color 0x1D6CB1, 5Y9J and chain L and resi 31
color 0x2575B7, 5Y9J and chain L and resi 28
color 0x2070B4, 5Y9J and chain L and resi 50
color 0x776BB0, 5Y9J and chain A and resi 162
color 0x501F8B, 5Y9J and chain A and resi 163
color 0x7E7AB8, 5Y9J and chain A and resi 206
color 0x684DA1, 5Y9J and chain A and resi 231
color 0x7D77B7, 5Y9J and chain A and resi 233
color 0x7F7CB9, 5Y9J and chain A and resi 240
color 0x7261AB, 5Y9J and chain A and resi 242
color 0x6A51A3, 5Y9J and chain A and resi 265
color 0x776CB1, 5Y9J and chain A and resi 266
color 0x6950A2, 5Y9J and chain A and resi 224
color 0x8F8BC1, 5Y9J and chain A and resi 222
color 0x6C56A6, 5Y9J and chain A and resi 225

distance hb10, (5Y9J and chain A and resi 240 and name O), (5Y9J and chain A and resi 242 and name N)
hide labels, hb10
set dash_color, yellow, hb10
hide labels
show sticks, (5Y9J and chain A and resi 206) and (sidechain or name CA)
show spheres, (5Y9J and chain A and resi 206) and name CA
set sphere_scale, 0.20, (5Y9J and chain A and resi 206)
show sticks, (5Y9J and chain A and resi 240) and (sidechain or name CA)
show spheres, (5Y9J and chain A and resi 240) and name CA
set sphere_scale, 0.20, (5Y9J and chain A and resi 240)
show sticks, (5Y9J and chain A and resi 242) and (sidechain or name CA)
show spheres, (5Y9J and chain A and resi 242) and name CA
set sphere_scale, 0.20, (5Y9J and chain A and resi 242)
pseudoatom lab_A206, pos=[15.014, -3.491, 63.517], label="TYR-206"
pseudoatom lab_A240, pos=[8.882, -13.03, 68.631], label="LEU-240"
pseudoatom lab_A242, pos=[12.437, -15.303, 63.926], label="ASN-242"
set label_color, black, lab_*

set_view (\
   0.4613189, 0.4955139, -0.7359694,\
   0.8272718, -0.5400040, 0.1549746,\
   -0.3206343, -0.6803394, -0.6590387,\
   0.0000000, 0.0000000, -60.9150391,\
   12.7574673, -9.3059654, 64.6723557,\
   48.0258789, 73.8041992, -20.0000000 )

set ray_trace_fog, 0
ray 1300, 1000
png pymol_groundtruth/figures_composite/5Y9J_ours_inset1.png, dpi=300, ray=0
