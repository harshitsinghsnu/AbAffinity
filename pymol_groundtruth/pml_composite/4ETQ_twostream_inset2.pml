reinitialize
load pymol_groundtruth/pdb/pdb4etq.ent, 4ETQ
remove not polymer
remove 4ETQ and not chain H+L+C
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
color 0xCECEE5, chain C
color 0x258D46, 4ETQ and chain H and resi 29
color 0x40AA5C, 4ETQ and chain H and resi 31
color 0x2E974E, 4ETQ and chain H and resi 32
color 0x50B264, 4ETQ and chain H and resi 34
color 0x56B567, 4ETQ and chain H and resi 51
color 0x90D18D, 4ETQ and chain H and resi 53
color 0x5EB96B, 4ETQ and chain H and resi 55
color 0x6BBF71, 4ETQ and chain H and resi 56
color 0x50B264, 4ETQ and chain H and resi 58
color 0x4CB062, 4ETQ and chain H and resi 60
color 0x329B51, 4ETQ and chain H and resi 75
color 0x61BA6C, 4ETQ and chain H and resi 102
color 0x6CC072, 4ETQ and chain H and resi 103
color 0x73C375, 4ETQ and chain H and resi 105
color 0x4B98C9, 4ETQ and chain L and resi 31
color 0x82BADB, 4ETQ and chain L and resi 91
color 0x88BEDC, 4ETQ and chain L and resi 92
color 0x65AAD3, 4ETQ and chain L and resi 94
color 0x7466AE, 4ETQ and chain C and resi 3
color 0x8C89BF, 4ETQ and chain C and resi 5
color 0x5F3D99, 4ETQ and chain C and resi 39
color 0x7F7CB9, 4ETQ and chain C and resi 40
color 0x674CA0, 4ETQ and chain C and resi 41
color 0x5E3A98, 4ETQ and chain C and resi 44
color 0x5B3595, 4ETQ and chain C and resi 108
color 0x674CA0, 4ETQ and chain C and resi 145
color 0x7160AB, 4ETQ and chain C and resi 174
color 0x6950A2, 4ETQ and chain C and resi 175
color 0x6C55A5, 4ETQ and chain C and resi 176
color 0x796FB3, 4ETQ and chain C and resi 177
color 0x776BB0, 4ETQ and chain C and resi 204
color 0x7567AE, 4ETQ and chain C and resi 205
color 0x6950A2, 4ETQ and chain C and resi 215
color 0x7364AD, 4ETQ and chain C and resi 217
color 0x684EA1, 4ETQ and chain C and resi 219
color 0x715FAA, 4ETQ and chain C and resi 220
color 0x53258E, 4ETQ and chain C and resi 221
color 0x7364AD, 4ETQ and chain C and resi 223
color 0x562C91, 4ETQ and chain C and resi 224

distance hb0, (4ETQ and chain C and resi 39 and name N), (4ETQ and chain C and resi 40 and name N)
hide labels, hb0
set dash_color, yellow, hb0
distance hb1, (4ETQ and chain C and resi 39 and name O), (4ETQ and chain H and resi 31 and name ND2)
hide labels, hb1
set dash_color, yellow, hb1
distance hb5, (4ETQ and chain C and resi 44 and name NH1), (4ETQ and chain H and resi 56 and name OE1)
hide labels, hb5
set dash_color, yellow, hb5
distance hb30, (4ETQ and chain H and resi 53 and name N), (4ETQ and chain H and resi 58 and name O)
hide labels, hb30
set dash_color, yellow, hb30
distance hb31, (4ETQ and chain H and resi 53 and name O), (4ETQ and chain H and resi 55 and name N)
hide labels, hb31
set dash_color, yellow, hb31
distance hb32, (4ETQ and chain H and resi 53 and name O), (4ETQ and chain H and resi 56 and name N)
hide labels, hb32
set dash_color, yellow, hb32
distance hb33, (4ETQ and chain H and resi 55 and name N), (4ETQ and chain H and resi 56 and name N)
hide labels, hb33
set dash_color, yellow, hb33
distance hb34, (4ETQ and chain H and resi 56 and name O), (4ETQ and chain H and resi 58 and name N)
hide labels, hb34
set dash_color, yellow, hb34
hide labels
show sticks, (4ETQ and chain H and resi 31) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 31) and name CA
set sphere_scale, 0.20, (4ETQ and chain H and resi 31)
show sticks, (4ETQ and chain H and resi 53) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 53) and name CA
set sphere_scale, 0.20, (4ETQ and chain H and resi 53)
show sticks, (4ETQ and chain H and resi 55) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 55) and name CA
set sphere_scale, 0.20, (4ETQ and chain H and resi 55)
show sticks, (4ETQ and chain H and resi 56) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 56) and name CA
set sphere_scale, 0.20, (4ETQ and chain H and resi 56)
show sticks, (4ETQ and chain H and resi 58) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 58) and name CA
set sphere_scale, 0.20, (4ETQ and chain H and resi 58)
show sticks, (4ETQ and chain H and resi 75) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 75) and name CA
set sphere_scale, 0.20, (4ETQ and chain H and resi 75)
show sticks, (4ETQ and chain C and resi 39) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 39) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 39)
show sticks, (4ETQ and chain C and resi 40) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 40) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 40)
show sticks, (4ETQ and chain C and resi 44) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 44) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 44)
show sticks, (4ETQ and chain C and resi 145) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 145) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 145)
pseudoatom lab_H31, pos=[-22.22, 7.35, 26.142], label="ASN-31"
distance dl_H31, (4ETQ and chain H and resi 31 and name CA), lab_H31
hide labels, dl_H31
set dash_color, grey50, dl_H31
pseudoatom lab_H53, pos=[-11.47, 9.249, 33.127], label="ASP-53"
pseudoatom lab_H55, pos=[-13.326, 8.957, 21.294], label="SER-55"
distance dl_H55, (4ETQ and chain H and resi 55 and name CA), lab_H55
hide labels, dl_H55
set dash_color, grey50, dl_H55
pseudoatom lab_H56, pos=[-14.615, 13.675, 36.059], label="GLU-56"
pseudoatom lab_H58, pos=[-2.236, 14.434, 42.257], label="GLU-58"
distance dl_H58, (4ETQ and chain H and resi 58 and name CA), lab_H58
hide labels, dl_H58
set dash_color, grey50, dl_H58
pseudoatom lab_H75, pos=[-19.696, 14.567, 28.345], label="ARG-75"
pseudoatom lab_C39, pos=[-20.025, 13.247, 51.676], label="THR-39"
distance dl_C39, (4ETQ and chain C and resi 39 and name CA), lab_C39
hide labels, dl_C39
set dash_color, grey50, dl_C39
pseudoatom lab_C40, pos=[-31.114, 7.905, 30.695], label="GLY-40"
distance dl_C40, (4ETQ and chain C and resi 40 and name CA), lab_C40
hide labels, dl_C40
set dash_color, grey50, dl_C40
pseudoatom lab_C44, pos=[-18.649, 12.022, 45.473], label="ARG-44"
pseudoatom lab_C145, pos=[-23.661, 13.981, 35.942], label="ASN-145"
set label_color, black, lab_*

set_view (\
   0.4582287, -0.8888344, 0.0003696,\
   0.1207343, 0.0618313, -0.9907573,\
   0.8805962, 0.4540381, 0.1356457,\
   0.0000000, 0.0000000, -65.9916382,\
   -16.6757660, 12.1625509, 36.3398933,\
   52.0283127, 79.9549637, -20.0000000 )

set ray_trace_fog, 0
ray 1300, 1000
png pymol_groundtruth/figures_composite/4ETQ_twostream_inset2.png, dpi=300, ray=0
