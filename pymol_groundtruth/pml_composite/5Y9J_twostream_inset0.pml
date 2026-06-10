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
color 0x37A055, 5Y9J and chain H and resi 31
color 0x3CA659, 5Y9J and chain H and resi 54
color 0x30994F, 5Y9J and chain H and resi 55
color 0x4EB163, 5Y9J and chain H and resi 57
color 0x17803C, 5Y9J and chain H and resi 101
color 0x04702F, 5Y9J and chain H and resi 102
color 0x006629, 5Y9J and chain H and resi 103
color 0x006127, 5Y9J and chain H and resi 104
color 0x2C7CBB, 5Y9J and chain L and resi 95
color 0x4795C8, 5Y9J and chain L and resi 31
color 0x4896C8, 5Y9J and chain L and resi 28
color 0x8BC0DD, 5Y9J and chain L and resi 50
color 0x5B3595, 5Y9J and chain A and resi 162
color 0x5A3495, 5Y9J and chain A and resi 163
color 0x6E59A7, 5Y9J and chain A and resi 206
color 0x7A72B4, 5Y9J and chain A and resi 231
color 0x9C98C7, 5Y9J and chain A and resi 233
color 0x7568AF, 5Y9J and chain A and resi 240
color 0x7261AB, 5Y9J and chain A and resi 242
color 0x5D3897, 5Y9J and chain A and resi 265
color 0x776BB0, 5Y9J and chain A and resi 266
color 0x6A51A3, 5Y9J and chain A and resi 224
color 0x8F8BC1, 5Y9J and chain A and resi 222
color 0x776BB0, 5Y9J and chain A and resi 225

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
distance hb9, (5Y9J and chain A and resi 231 and name NH1), (5Y9J and chain H and resi 103 and name O)
hide labels, hb9
set dash_color, yellow, hb9
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
hide labels
show sticks, (5Y9J and chain H and resi 31) and (sidechain or name CA)
show spheres, (5Y9J and chain H and resi 31) and name CA
set sphere_scale, 0.20, (5Y9J and chain H and resi 31)
show sticks, (5Y9J and chain H and resi 54) and (sidechain or name CA)
show spheres, (5Y9J and chain H and resi 54) and name CA
set sphere_scale, 0.20, (5Y9J and chain H and resi 54)
show sticks, (5Y9J and chain H and resi 55) and (sidechain or name CA)
show spheres, (5Y9J and chain H and resi 55) and name CA
set sphere_scale, 0.20, (5Y9J and chain H and resi 55)
show sticks, (5Y9J and chain H and resi 57) and (sidechain or name CA)
show spheres, (5Y9J and chain H and resi 57) and name CA
set sphere_scale, 0.20, (5Y9J and chain H and resi 57)
show sticks, (5Y9J and chain H and resi 101) and (sidechain or name CA)
show spheres, (5Y9J and chain H and resi 101) and name CA
set sphere_scale, 0.20, (5Y9J and chain H and resi 101)
show sticks, (5Y9J and chain H and resi 102) and (sidechain or name CA)
show spheres, (5Y9J and chain H and resi 102) and name CA
set sphere_scale, 0.20, (5Y9J and chain H and resi 102)
show sticks, (5Y9J and chain H and resi 103) and (sidechain or name CA)
show spheres, (5Y9J and chain H and resi 103) and name CA
set sphere_scale, 0.20, (5Y9J and chain H and resi 103)
show sticks, (5Y9J and chain H and resi 104) and (sidechain or name CA)
show spheres, (5Y9J and chain H and resi 104) and name CA
set sphere_scale, 0.20, (5Y9J and chain H and resi 104)
show sticks, (5Y9J and chain L and resi 95) and (sidechain or name CA)
show spheres, (5Y9J and chain L and resi 95) and name CA
set sphere_scale, 0.20, (5Y9J and chain L and resi 95)
show sticks, (5Y9J and chain A and resi 162) and (sidechain or name CA)
show spheres, (5Y9J and chain A and resi 162) and name CA
set sphere_scale, 0.20, (5Y9J and chain A and resi 162)
show sticks, (5Y9J and chain A and resi 163) and (sidechain or name CA)
show spheres, (5Y9J and chain A and resi 163) and name CA
set sphere_scale, 0.20, (5Y9J and chain A and resi 163)
show sticks, (5Y9J and chain A and resi 231) and (sidechain or name CA)
show spheres, (5Y9J and chain A and resi 231) and name CA
set sphere_scale, 0.20, (5Y9J and chain A and resi 231)
show sticks, (5Y9J and chain A and resi 233) and (sidechain or name CA)
show spheres, (5Y9J and chain A and resi 233) and name CA
set sphere_scale, 0.20, (5Y9J and chain A and resi 233)
show sticks, (5Y9J and chain A and resi 265) and (sidechain or name CA)
show spheres, (5Y9J and chain A and resi 265) and name CA
set sphere_scale, 0.20, (5Y9J and chain A and resi 265)
show sticks, (5Y9J and chain A and resi 266) and (sidechain or name CA)
show spheres, (5Y9J and chain A and resi 266) and name CA
set sphere_scale, 0.20, (5Y9J and chain A and resi 266)
pseudoatom lab_H31, pos=[11.507, 3.962, 49.756], label="ASN-31"
distance dl_H31, (5Y9J and chain H and resi 31 and name CA), lab_H31
hide labels, dl_H31
set dash_color, grey50, dl_H31
pseudoatom lab_H54, pos=[18.256, 2.232, 68.565], label="MET-54"
pseudoatom lab_H55, pos=[20.409, -0.749, 69.633], label="PHE-55"
pseudoatom lab_H57, pos=[25.816, -0.297, 73.871], label="THR-57"
pseudoatom lab_H101, pos=[24.552, 1.425, 60.535], label="ASP-101"
pseudoatom lab_H102, pos=[27.397, -0.81, 59.303], label="LEU-102"
pseudoatom lab_H103, pos=[31.784, -10.184, 45.69], label="LEU-103"
distance dl_H103, (5Y9J and chain H and resi 103 and name CA), lab_H103
hide labels, dl_H103
set dash_color, grey50, dl_H103
pseudoatom lab_H104, pos=[25.493, -4.492, 46.773], label="LEU-104"
distance dl_H104, (5Y9J and chain H and resi 104 and name CA), lab_H104
hide labels, dl_H104
set dash_color, grey50, dl_H104
pseudoatom lab_L95, pos=[35.257, -0.775, 64.7], label="ASN-95"
pseudoatom lab_A162, pos=[30.253, -4.423, 66.745], label="SER-162"
pseudoatom lab_A163, pos=[38.911, -9.071, 75.088], label="TYR-163"
distance dl_A163, (5Y9J and chain A and resi 163 and name CA), lab_A163
hide labels, dl_A163
set dash_color, grey50, dl_A163
pseudoatom lab_A231, pos=[25.322, -9.423, 52.081], label="ARG-231"
pseudoatom lab_A233, pos=[14.996, -7.578, 50.304], label="ILE-233"
distance dl_A233, (5Y9J and chain A and resi 233 and name CA), lab_A233
hide labels, dl_A233
set dash_color, grey50, dl_A233
pseudoatom lab_A265, pos=[17.069, 1.811, 80.064], label="ARG-265"
distance dl_A265, (5Y9J and chain A and resi 265 and name CA), lab_A265
hide labels, dl_A265
set dash_color, grey50, dl_A265
pseudoatom lab_A266, pos=[27.682, -4.272, 77.755], label="GLU-266"
distance dl_A266, (5Y9J and chain A and resi 266 and name CA), lab_A266
hide labels, dl_A266
set dash_color, grey50, dl_A266
set label_color, black, lab_*

set_view (\
   -0.2814932, 0.8631999, -0.4191031,\
   -0.1302924, -0.4671030, -0.8745506,\
   -0.9506763, -0.1915741, 0.2439547,\
   -0.0000000, 0.0000000, -82.0670471,\
   25.0446472, -2.8968241, 62.5173645,\
   64.7022858, 99.4318085, -20.0000000 )

set ray_trace_fog, 0
ray 1300, 1000
png pymol_groundtruth/figures_composite/5Y9J_twostream_inset0.png, dpi=300, ray=0
