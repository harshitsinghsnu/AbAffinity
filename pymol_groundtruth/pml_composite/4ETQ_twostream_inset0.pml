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

distance hb4, (4ETQ and chain C and resi 41 and name NZ), (4ETQ and chain H and resi 32 and name O)
hide labels, hb4
set dash_color, yellow, hb4
distance hb15, (4ETQ and chain C and resi 204 and name N), (4ETQ and chain C and resi 205 and name N)
hide labels, hb15
set dash_color, yellow, hb15
distance hb16, (4ETQ and chain C and resi 217 and name OE1), (4ETQ and chain H and resi 102 and name OH)
hide labels, hb16
set dash_color, yellow, hb16
distance hb17, (4ETQ and chain C and resi 217 and name OE2), (4ETQ and chain H and resi 103 and name NH1)
hide labels, hb17
set dash_color, yellow, hb17
distance hb18, (4ETQ and chain C and resi 219 and name O), (4ETQ and chain C and resi 220 and name N)
hide labels, hb18
set dash_color, yellow, hb18
distance hb19, (4ETQ and chain C and resi 219 and name O), (4ETQ and chain H and resi 103 and name NE)
hide labels, hb19
set dash_color, yellow, hb19
distance hb21, (4ETQ and chain C and resi 220 and name NH1), (4ETQ and chain H and resi 102 and name O)
hide labels, hb21
set dash_color, yellow, hb21
distance hb27, (4ETQ and chain H and resi 29 and name O), (4ETQ and chain H and resi 32 and name N)
hide labels, hb27
set dash_color, yellow, hb27
distance hb35, (4ETQ and chain H and resi 102 and name N), (4ETQ and chain H and resi 103 and name N)
hide labels, hb35
set dash_color, yellow, hb35
hide labels
show sticks, (4ETQ and chain H and resi 29) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 29) and name CA
set sphere_scale, 0.20, (4ETQ and chain H and resi 29)
show sticks, (4ETQ and chain H and resi 32) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 32) and name CA
set sphere_scale, 0.20, (4ETQ and chain H and resi 32)
show sticks, (4ETQ and chain H and resi 34) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 34) and name CA
set sphere_scale, 0.20, (4ETQ and chain H and resi 34)
show sticks, (4ETQ and chain H and resi 102) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 102) and name CA
set sphere_scale, 0.20, (4ETQ and chain H and resi 102)
show sticks, (4ETQ and chain H and resi 103) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 103) and name CA
set sphere_scale, 0.20, (4ETQ and chain H and resi 103)
show sticks, (4ETQ and chain C and resi 41) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 41) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 41)
show sticks, (4ETQ and chain C and resi 204) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 204) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 204)
show sticks, (4ETQ and chain C and resi 205) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 205) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 205)
show sticks, (4ETQ and chain C and resi 215) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 215) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 215)
show sticks, (4ETQ and chain C and resi 217) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 217) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 217)
show sticks, (4ETQ and chain C and resi 219) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 219) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 219)
show sticks, (4ETQ and chain C and resi 220) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 220) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 220)
pseudoatom lab_H29, pos=[-22.122, 4.466, 28.784], label="SER-29"
pseudoatom lab_H32, pos=[-19.181, 12.338, 25.387], label="PHE-32"
distance dl_H32, (4ETQ and chain H and resi 32 and name CA), lab_H32
hide labels, dl_H32
set dash_color, grey50, dl_H32
pseudoatom lab_H34, pos=[-11.384, 5.307, 30.564], label="TRP-34"
pseudoatom lab_H102, pos=[-13.888, -0.766, 33.284], label="TYR-102"
pseudoatom lab_H103, pos=[-13.071, -4.114, 34.879], label="ARG-103"
pseudoatom lab_C41, pos=[-32.056, 6.531, 40.145], label="LYS-41"
distance dl_C41, (4ETQ and chain C and resi 41 and name CA), lab_C41
hide labels, dl_C41
set dash_color, grey50, dl_C41
pseudoatom lab_C204, pos=[-24.708, 6.036, 34.529], label="SER-204"
pseudoatom lab_C205, pos=[-25.434, 7.55, 31.091], label="LEU-205"
pseudoatom lab_C215, pos=[-25.004, 1.482, 37.758], label="ILE-215"
pseudoatom lab_C217, pos=[-15.424, -10.182, 46.713], label="GLU-217"
distance dl_C217, (4ETQ and chain C and resi 217 and name CA), lab_C217
hide labels, dl_C217
set dash_color, grey50, dl_C217
pseudoatom lab_C219, pos=[-15.372, -1.383, 42.024], label="TYR-219"
pseudoatom lab_C220, pos=[-1.655, -0.192, 36.262], label="ARG-220"
distance dl_C220, (4ETQ and chain C and resi 220 and name CA), lab_C220
hide labels, dl_C220
set dash_color, grey50, dl_C220
set label_color, black, lab_*

set_view (\
   0.7177585, 0.6805339, 0.1472969,\
   -0.5680652, 0.4499952, 0.6890618,\
   0.4026470, -0.5782542, 0.7095755,\
   0.0000000, 0.0000000, -65.6457748,\
   -17.3478794, 0.8668375, 35.8320084,\
   51.7556305, 79.5359192, -20.0000000 )

set ray_trace_fog, 0
ray 1300, 1000
png pymol_groundtruth/figures_composite/4ETQ_twostream_inset0.png, dpi=300, ray=0
