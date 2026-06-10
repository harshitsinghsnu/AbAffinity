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
set stick_radius, 0.26
set sphere_scale, 0.30
set valence, 0
set two_sided_lighting, 1
set label_size, 17
set label_font_id, 7
set label_color, black
set label_outline_color, white
set label_bg_color, white
set label_bg_transparency, 0.2
set float_labels, 1
set dash_color, yellow
set dash_gap, 0.30
set dash_width, 3.0
set dash_radius, 0.06
set cartoon_transparency, 0.0
set ray_trace_mode, 1
set ray_trace_color, grey60
show cartoon
color 0xC4E7BD, chain H
color 0xCBDEF0, chain L
color 0xCECEE5, chain A
show sticks, (5Y9J and chain H and resi 31) and (sidechain or name CA)
color 0x157E3A, 5Y9J and chain H and resi 31
show sticks, (5Y9J and chain H and resi 54) and (sidechain or name CA)
color 0x339C51, 5Y9J and chain H and resi 54
show sticks, (5Y9J and chain H and resi 55) and (sidechain or name CA)
color 0x2E974E, 5Y9J and chain H and resi 55
show sticks, (5Y9J and chain H and resi 57) and (sidechain or name CA)
color 0x278F48, 5Y9J and chain H and resi 57
show sticks, (5Y9J and chain H and resi 101) and (sidechain or name CA)
color 0x76C477, 5Y9J and chain H and resi 101
show sticks, (5Y9J and chain H and resi 102) and (sidechain or name CA)
color 0x6EC173, 5Y9J and chain H and resi 102
show sticks, (5Y9J and chain H and resi 103) and (sidechain or name CA)
color 0x7AC67B, 5Y9J and chain H and resi 103
show sticks, (5Y9J and chain H and resi 104) and (sidechain or name CA)
color 0x64BC6E, 5Y9J and chain H and resi 104
show sticks, (5Y9J and chain L and resi 95) and (sidechain or name CA)
color 0x2777B8, 5Y9J and chain L and resi 95
show sticks, (5Y9J and chain L and resi 31) and (sidechain or name CA)
color 0x1D6CB1, 5Y9J and chain L and resi 31
show sticks, (5Y9J and chain L and resi 28) and (sidechain or name CA)
color 0x2575B7, 5Y9J and chain L and resi 28
show sticks, (5Y9J and chain L and resi 50) and (sidechain or name CA)
color 0x2070B4, 5Y9J and chain L and resi 50
show sticks, (5Y9J and chain A and resi 162) and (sidechain or name CA)
color 0x776BB0, 5Y9J and chain A and resi 162
show sticks, (5Y9J and chain A and resi 163) and (sidechain or name CA)
color 0x501F8B, 5Y9J and chain A and resi 163
show sticks, (5Y9J and chain A and resi 206) and (sidechain or name CA)
color 0x7E7AB8, 5Y9J and chain A and resi 206
show sticks, (5Y9J and chain A and resi 231) and (sidechain or name CA)
color 0x684DA1, 5Y9J and chain A and resi 231
show sticks, (5Y9J and chain A and resi 233) and (sidechain or name CA)
color 0x7D77B7, 5Y9J and chain A and resi 233
show sticks, (5Y9J and chain A and resi 240) and (sidechain or name CA)
color 0x7F7CB9, 5Y9J and chain A and resi 240
show sticks, (5Y9J and chain A and resi 242) and (sidechain or name CA)
color 0x7261AB, 5Y9J and chain A and resi 242
show sticks, (5Y9J and chain A and resi 265) and (sidechain or name CA)
color 0x6A51A3, 5Y9J and chain A and resi 265
show sticks, (5Y9J and chain A and resi 266) and (sidechain or name CA)
color 0x776CB1, 5Y9J and chain A and resi 266
show sticks, (5Y9J and chain A and resi 224) and (sidechain or name CA)
color 0x6950A2, 5Y9J and chain A and resi 224
show sticks, (5Y9J and chain A and resi 222) and (sidechain or name CA)
color 0x8F8BC1, 5Y9J and chain A and resi 222
show sticks, (5Y9J and chain A and resi 225) and (sidechain or name CA)
color 0x6C56A6, 5Y9J and chain A and resi 225

distance hb0, (5Y9J and chain A and resi 162 and name N), (5Y9J and chain A and resi 163 and name N)
hide labels, hb0
set dash_color, yellow, hb0
distance hb1, (5Y9J and chain A and resi 162 and name N), (5Y9J and chain L and resi 95 and name OD1)
hide labels, hb1
set dash_color, yellow, hb1
distance hb2, (5Y9J and chain A and resi 162 and name O), (5Y9J and chain A and resi 265 and name N)
hide labels, hb2
set dash_color, yellow, hb2
distance hb3, (5Y9J and chain A and resi 162 and name O), (5Y9J and chain A and resi 266 and name N)
hide labels, hb3
set dash_color, yellow, hb3
distance hb4, (5Y9J and chain A and resi 206 and name OH), (5Y9J and chain H and resi 31 and name ND2)
hide labels, hb4
set dash_color, yellow, hb4
distance hb5, (5Y9J and chain A and resi 206 and name OH), (5Y9J and chain H and resi 101 and name OD2)
hide labels, hb5
set dash_color, yellow, hb5
distance hb6, (5Y9J and chain A and resi 222 and name OD2), (5Y9J and chain L and resi 50 and name NZ)
hide labels, hb6
set dash_color, yellow, hb6
distance hb7, (5Y9J and chain A and resi 224 and name O), (5Y9J and chain A and resi 225 and name N)
hide labels, hb7
set dash_color, yellow, hb7
distance hb8, (5Y9J and chain A and resi 225 and name N), (5Y9J and chain L and resi 28 and name O)
hide labels, hb8
set dash_color, yellow, hb8
distance hb9, (5Y9J and chain A and resi 231 and name NH1), (5Y9J and chain H and resi 103 and name O)
hide labels, hb9
set dash_color, yellow, hb9
distance hb10, (5Y9J and chain A and resi 240 and name O), (5Y9J and chain A and resi 242 and name N)
hide labels, hb10
set dash_color, yellow, hb10
distance hb11, (5Y9J and chain A and resi 265 and name O), (5Y9J and chain A and resi 266 and name N)
hide labels, hb11
set dash_color, yellow, hb11
distance hb12, (5Y9J and chain A and resi 265 and name NH1), (5Y9J and chain H and resi 101 and name OD1)
hide labels, hb12
set dash_color, yellow, hb12
distance hb13, (5Y9J and chain A and resi 266 and name OE1), (5Y9J and chain H and resi 57 and name N)
hide labels, hb13
set dash_color, yellow, hb13
distance hb14, (5Y9J and chain H and resi 54 and name N), (5Y9J and chain H and resi 55 and name N)
hide labels, hb14
set dash_color, yellow, hb14
distance hb15, (5Y9J and chain H and resi 101 and name N), (5Y9J and chain H and resi 102 and name N)
hide labels, hb15
set dash_color, yellow, hb15
distance hb16, (5Y9J and chain H and resi 101 and name O), (5Y9J and chain H and resi 103 and name N)
hide labels, hb16
set dash_color, yellow, hb16
distance hb17, (5Y9J and chain H and resi 101 and name O), (5Y9J and chain H and resi 104 and name N)
hide labels, hb17
set dash_color, yellow, hb17
distance hb18, (5Y9J and chain H and resi 102 and name N), (5Y9J and chain H and resi 103 and name N)
hide labels, hb18
set dash_color, yellow, hb18
distance hb19, (5Y9J and chain H and resi 102 and name O), (5Y9J and chain H and resi 104 and name N)
hide labels, hb19
set dash_color, yellow, hb19
distance hb20, (5Y9J and chain H and resi 103 and name N), (5Y9J and chain H and resi 104 and name N)
hide labels, hb20
set dash_color, yellow, hb20
set label_color, black

set_view (\
   -0.1056180, 0.6314777, 0.7681672,\
   -0.9737224, 0.0910630, -0.2087395,\
   -0.2017659, -0.7700282, 0.6052660,\
   -0.0000000, 0.0000000, -324.3650208,\
   34.5380936, 16.8643036, 62.8831329,\
   255.7318420, 392.9981995, -20.0000000 )

set ray_trace_fog, 0
ray 2400, 2000
png pymol_groundtruth/figures_composite/5Y9J_ours_overview.png, dpi=300, ray=0
