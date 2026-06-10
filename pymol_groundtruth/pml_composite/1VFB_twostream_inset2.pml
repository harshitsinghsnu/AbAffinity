reinitialize
load pymol_groundtruth/pdb/pdb1vfb.ent, 1VFB
remove not polymer
remove 1VFB and not chain B+A+C
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
color 0xC4E7BD, chain B
color 0xCBDEF0, chain A
color 0xCECEE5, chain C
color 0x80C97F, 1VFB and chain B and resi 30
color 0x82CA82, 1VFB and chain B and resi 31
color 0x90D18D, 1VFB and chain B and resi 32
color 0x319A50, 1VFB and chain B and resi 52
color 0x38A256, 1VFB and chain B and resi 53
color 0x41AB5D, 1VFB and chain B and resi 54
color 0x80C97F, 1VFB and chain B and resi 56
color 0x70C274, 1VFB and chain B and resi 58
color 0x51B365, 1VFB and chain B and resi 99
color 0x43AC5E, 1VFB and chain B and resi 100
color 0x2D964D, 1VFB and chain B and resi 101
color 0x4EB163, 1VFB and chain B and resi 102
color 0x2676B7, 1VFB and chain A and resi 30
color 0x8BC0DD, 1VFB and chain A and resi 32
color 0x3282BE, 1VFB and chain A and resi 49
color 0x2272B5, 1VFB and chain A and resi 50
color 0x7BB7D9, 1VFB and chain A and resi 52
color 0x519CCC, 1VFB and chain A and resi 53
color 0x509BCB, 1VFB and chain A and resi 91
color 0x60A6D1, 1VFB and chain A and resi 92
color 0x539DCC, 1VFB and chain A and resi 93
color 0x3585BF, 1VFB and chain A and resi 94
color 0x776BB0, 1VFB and chain C and resi 18
color 0x928EC2, 1VFB and chain C and resi 19
color 0x7C75B5, 1VFB and chain C and resi 22
color 0x664AA0, 1VFB and chain C and resi 24
color 0x583093, 1VFB and chain C and resi 26
color 0x7A71B3, 1VFB and chain C and resi 102
color 0x7364AD, 1VFB and chain C and resi 116
color 0x7D77B7, 1VFB and chain C and resi 121
color 0x63449C, 1VFB and chain C and resi 129

distance hb30, (1VFB and chain B and resi 100 and name OD1), (1VFB and chain C and resi 24 and name OG)
hide labels, hb30
set dash_color, yellow, hb30
distance hb31, (1VFB and chain B and resi 100 and name OD2), (1VFB and chain C and resi 22 and name O)
hide labels, hb31
set dash_color, yellow, hb31
distance hb35, (1VFB and chain C and resi 18 and name N), (1VFB and chain C and resi 19 and name N)
hide labels, hb35
set dash_color, yellow, hb35
distance hb36, (1VFB and chain C and resi 19 and name O), (1VFB and chain C and resi 22 and name N)
hide labels, hb36
set dash_color, yellow, hb36
distance hb37, (1VFB and chain C and resi 24 and name O), (1VFB and chain C and resi 26 and name N)
hide labels, hb37
set dash_color, yellow, hb37
hide labels
show sticks, (1VFB and chain B and resi 100) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 100) and name CA
set sphere_scale, 0.20, (1VFB and chain B and resi 100)
show sticks, (1VFB and chain C and resi 18) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 18) and name CA
set sphere_scale, 0.20, (1VFB and chain C and resi 18)
show sticks, (1VFB and chain C and resi 19) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 19) and name CA
set sphere_scale, 0.20, (1VFB and chain C and resi 19)
show sticks, (1VFB and chain C and resi 22) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 22) and name CA
set sphere_scale, 0.20, (1VFB and chain C and resi 22)
show sticks, (1VFB and chain C and resi 24) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 24) and name CA
set sphere_scale, 0.20, (1VFB and chain C and resi 24)
show sticks, (1VFB and chain C and resi 26) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 26) and name CA
set sphere_scale, 0.20, (1VFB and chain C and resi 26)
show sticks, (1VFB and chain C and resi 102) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 102) and name CA
set sphere_scale, 0.20, (1VFB and chain C and resi 102)
show sticks, (1VFB and chain C and resi 121) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 121) and name CA
set sphere_scale, 0.20, (1VFB and chain C and resi 121)
pseudoatom lab_B100, pos=[44.16, -7.7, -5.079], label="ASP-100"
pseudoatom lab_C18, pos=[48.156, -4.668, 5.971], label="ASP-18"
pseudoatom lab_C19, pos=[49.179, -3.643, 2.437], label="ASN-19"
pseudoatom lab_C22, pos=[50.823, -4.028, -1.856], label="GLY-22"
pseudoatom lab_C24, pos=[34.548, -5.889, 12.91], label="SER-24"
distance dl_C24, (1VFB and chain C and resi 24 and name CA), lab_C24
hide labels, dl_C24
set dash_color, grey50, dl_C24
pseudoatom lab_C26, pos=[43.135, -10.923, 3.83], label="GLY-26"
pseudoatom lab_C102, pos=[54.946, -7.513, -4.755], label="GLY-102"
pseudoatom lab_C121, pos=[33.822, -15.688, -6.017], label="GLN-121"
distance dl_C121, (1VFB and chain C and resi 121 and name CA), lab_C121
hide labels, dl_C121
set dash_color, grey50, dl_C121
set label_color, black, lab_*

set_view (\
   -0.9380272, 0.0444029, -0.3437054,\
   -0.2917644, 0.4340312, 0.8523441,\
   0.1870254, 0.8998029, -0.3941779,\
   0.0000000, 0.0000000, -70.2840576,\
   45.9740219, -6.2818604, 0.6911407,\
   55.4124870, 85.1556244, -20.0000000 )

set ray_trace_fog, 0
ray 1300, 1000
png pymol_groundtruth/figures_composite/1VFB_twostream_inset2.png, dpi=300, ray=0
