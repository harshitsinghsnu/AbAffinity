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
color 0x339C51, 4ETQ and chain H and resi 29
color 0x1D8641, 4ETQ and chain H and resi 31
color 0x258D46, 4ETQ and chain H and resi 32
color 0x68BE70, 4ETQ and chain H and resi 34
color 0x7BC77C, 4ETQ and chain H and resi 51
color 0x41AB5D, 4ETQ and chain H and resi 53
color 0x41AB5D, 4ETQ and chain H and resi 55
color 0x1D8641, 4ETQ and chain H and resi 56
color 0x117C38, 4ETQ and chain H and resi 58
color 0x2F984E, 4ETQ and chain H and resi 60
color 0x43AC5E, 4ETQ and chain H and resi 75
color 0x04702F, 4ETQ and chain H and resi 102
color 0x016E2C, 4ETQ and chain H and resi 103
color 0x167F3B, 4ETQ and chain H and resi 105
color 0x2878B8, 4ETQ and chain L and resi 31
color 0x2676B7, 4ETQ and chain L and resi 91
color 0x1D6CB1, 4ETQ and chain L and resi 92
color 0x0B559F, 4ETQ and chain L and resi 94
color 0x7B73B5, 4ETQ and chain C and resi 3
color 0x7C76B6, 4ETQ and chain C and resi 5
color 0x8480BB, 4ETQ and chain C and resi 39
color 0x664AA0, 4ETQ and chain C and resi 40
color 0x807DBA, 4ETQ and chain C and resi 41
color 0x705DA9, 4ETQ and chain C and resi 44
color 0x7A71B3, 4ETQ and chain C and resi 108
color 0x8F8BC1, 4ETQ and chain C and resi 145
color 0x8683BD, 4ETQ and chain C and resi 174
color 0x776BB0, 4ETQ and chain C and resi 175
color 0x776BB0, 4ETQ and chain C and resi 176
color 0x8E8AC0, 4ETQ and chain C and resi 177
color 0x8F8BC1, 4ETQ and chain C and resi 204
color 0x786EB2, 4ETQ and chain C and resi 205
color 0x5A3394, 4ETQ and chain C and resi 215
color 0x6F5BA8, 4ETQ and chain C and resi 217
color 0x9794C5, 4ETQ and chain C and resi 219
color 0x7E7AB8, 4ETQ and chain C and resi 220
color 0x66499F, 4ETQ and chain C and resi 221
color 0x9693C4, 4ETQ and chain C and resi 223
color 0x63459D, 4ETQ and chain C and resi 224

distance hb6, (4ETQ and chain C and resi 174 and name N), (4ETQ and chain C and resi 175 and name N)
hide labels, hb6
set dash_color, yellow, hb6
distance hb7, (4ETQ and chain C and resi 174 and name O), (4ETQ and chain C and resi 176 and name N)
hide labels, hb7
set dash_color, yellow, hb7
distance hb8, (4ETQ and chain C and resi 175 and name N), (4ETQ and chain C and resi 176 and name N)
hide labels, hb8
set dash_color, yellow, hb8
distance hb9, (4ETQ and chain C and resi 175 and name O), (4ETQ and chain C and resi 177 and name N)
hide labels, hb9
set dash_color, yellow, hb9
distance hb12, (4ETQ and chain C and resi 176 and name N), (4ETQ and chain C and resi 177 and name N)
hide labels, hb12
set dash_color, yellow, hb12
distance hb13, (4ETQ and chain C and resi 176 and name O), (4ETQ and chain H and resi 60 and name NH1)
hide labels, hb13
set dash_color, yellow, hb13
distance hb14, (4ETQ and chain C and resi 176 and name ND1), (4ETQ and chain L and resi 94 and name OH)
hide labels, hb14
set dash_color, yellow, hb14
distance hb24, (4ETQ and chain C and resi 223 and name O), (4ETQ and chain C and resi 224 and name N)
hide labels, hb24
set dash_color, yellow, hb24
distance hb25, (4ETQ and chain C and resi 223 and name OH), (4ETQ and chain H and resi 105 and name OD1)
hide labels, hb25
set dash_color, yellow, hb25
distance hb29, (4ETQ and chain H and resi 51 and name N), (4ETQ and chain H and resi 60 and name O)
hide labels, hb29
set dash_color, yellow, hb29
distance hb36, (4ETQ and chain L and resi 31 and name N), (4ETQ and chain L and resi 92 and name OG1)
hide labels, hb36
set dash_color, yellow, hb36
distance hb37, (4ETQ and chain L and resi 91 and name N), (4ETQ and chain L and resi 92 and name N)
hide labels, hb37
set dash_color, yellow, hb37
distance hb38, (4ETQ and chain L and resi 91 and name NE1), (4ETQ and chain L and resi 94 and name OH)
hide labels, hb38
set dash_color, yellow, hb38
hide labels
show sticks, (4ETQ and chain H and resi 51) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 51) and name CA
set sphere_scale, 0.20, (4ETQ and chain H and resi 51)
show sticks, (4ETQ and chain H and resi 60) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 60) and name CA
set sphere_scale, 0.20, (4ETQ and chain H and resi 60)
show sticks, (4ETQ and chain H and resi 105) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 105) and name CA
set sphere_scale, 0.20, (4ETQ and chain H and resi 105)
show sticks, (4ETQ and chain L and resi 31) and (sidechain or name CA)
show spheres, (4ETQ and chain L and resi 31) and name CA
set sphere_scale, 0.20, (4ETQ and chain L and resi 31)
show sticks, (4ETQ and chain L and resi 91) and (sidechain or name CA)
show spheres, (4ETQ and chain L and resi 91) and name CA
set sphere_scale, 0.20, (4ETQ and chain L and resi 91)
show sticks, (4ETQ and chain L and resi 92) and (sidechain or name CA)
show spheres, (4ETQ and chain L and resi 92) and name CA
set sphere_scale, 0.20, (4ETQ and chain L and resi 92)
show sticks, (4ETQ and chain L and resi 94) and (sidechain or name CA)
show spheres, (4ETQ and chain L and resi 94) and name CA
set sphere_scale, 0.20, (4ETQ and chain L and resi 94)
show sticks, (4ETQ and chain C and resi 3) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 3) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 3)
show sticks, (4ETQ and chain C and resi 5) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 5) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 5)
show sticks, (4ETQ and chain C and resi 108) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 108) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 108)
show sticks, (4ETQ and chain C and resi 174) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 174) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 174)
show sticks, (4ETQ and chain C and resi 175) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 175) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 175)
show sticks, (4ETQ and chain C and resi 176) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 176) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 176)
show sticks, (4ETQ and chain C and resi 177) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 177) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 177)
show sticks, (4ETQ and chain C and resi 221) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 221) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 221)
show sticks, (4ETQ and chain C and resi 223) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 223) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 223)
show sticks, (4ETQ and chain C and resi 224) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 224) and name CA
set sphere_scale, 0.20, (4ETQ and chain C and resi 224)
pseudoatom lab_H51, pos=[-1.046, 18.565, 26.435], label="MET-51"
distance dl_H51, (4ETQ and chain H and resi 51 and name CA), lab_H51
hide labels, dl_H51
set dash_color, grey50, dl_H51
pseudoatom lab_H60, pos=[-3.022, 10.98, 33.253], label="ARG-60"
pseudoatom lab_H105, pos=[-7.225, -4.125, 31.789], label="ASP-105"
pseudoatom lab_L31, pos=[2.295, -8.855, 36.641], label="SER-31"
pseudoatom lab_L91, pos=[4.827, -16.11, 32.679], label="TRP-91"
distance dl_L91, (4ETQ and chain L and resi 91 and name CA), lab_L91
hide labels, dl_L91
set dash_color, grey50, dl_L91
pseudoatom lab_L92, pos=[12.861, -14.957, 33.599], label="THR-92"
distance dl_L92, (4ETQ and chain L and resi 92 and name CA), lab_L92
hide labels, dl_L92
set dash_color, grey50, dl_L92
pseudoatom lab_L94, pos=[2.499, 3.913, 34.187], label="TYR-94"
pseudoatom lab_C3, pos=[10.369, 20.972, 37.287], label="GLN-3"
distance dl_C3, (4ETQ and chain C and resi 3 and name CA), lab_C3
hide labels, dl_C3
set dash_color, grey50, dl_C3
pseudoatom lab_C5, pos=[0.444, -15.907, 42.487], label="LEU-5"
distance dl_C5, (4ETQ and chain C and resi 5 and name CA), lab_C5
hide labels, dl_C5
set dash_color, grey50, dl_C5
pseudoatom lab_C108, pos=[-3.269, 19.57, 42.588], label="LYS-108"
distance dl_C108, (4ETQ and chain C and resi 108 and name CA), lab_C108
hide labels, dl_C108
set dash_color, grey50, dl_C108
pseudoatom lab_C174, pos=[-18.346, 17.633, 43.66], label="ILE-174"
distance dl_C174, (4ETQ and chain C and resi 174 and name CA), lab_C174
hide labels, dl_C174
set dash_color, grey50, dl_C174
pseudoatom lab_C175, pos=[-9.023, 5.628, 38.642], label="ASN-175"
pseudoatom lab_C176, pos=[-12.165, 18.097, 39.136], label="HIS-176"
distance dl_C176, (4ETQ and chain C and resi 176 and name CA), lab_C176
hide labels, dl_C176
set dash_color, grey50, dl_C176
pseudoatom lab_C177, pos=[-5.029, 9.292, 37.479], label="SER-177"
pseudoatom lab_C221, pos=[-10.248, -5.691, 40.628], label="ASN-221"
pseudoatom lab_C223, pos=[-13.74, -17.512, 46.183], label="TYR-223"
distance dl_C223, (4ETQ and chain C and resi 223 and name CA), lab_C223
hide labels, dl_C223
set dash_color, grey50, dl_C223
pseudoatom lab_C224, pos=[-6.367, -16.605, 45.16], label="LYS-224"
distance dl_C224, (4ETQ and chain C and resi 224 and name CA), lab_C224
hide labels, dl_C224
set dash_color, grey50, dl_C224
set label_color, black, lab_*

set_view (\
   0.1325913, 0.9381021, 0.3199751,\
   -0.9879616, 0.0991286, 0.1187665,\
   0.0796964, -0.3318705, 0.9399524,\
   0.0000000, 0.0000000, -96.8861465,\
   -3.2416897, 1.1523581, 37.1629906,\
   76.3857727, 117.3865204, -20.0000000 )

set ray_trace_fog, 0
ray 1300, 1000
png pymol_groundtruth/figures_composite/4ETQ_ours_inset1.png, dpi=300, ray=0
