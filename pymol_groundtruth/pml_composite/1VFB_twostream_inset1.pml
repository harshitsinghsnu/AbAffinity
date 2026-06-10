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

distance hb15, (1VFB and chain B and resi 30 and name N), (1VFB and chain B and resi 31 and name N)
hide labels, hb15
set dash_color, yellow, hb15
distance hb16, (1VFB and chain B and resi 31 and name N), (1VFB and chain B and resi 32 and name N)
hide labels, hb16
set dash_color, yellow, hb16
distance hb17, (1VFB and chain B and resi 32 and name OH), (1VFB and chain C and resi 116 and name NZ)
hide labels, hb17
set dash_color, yellow, hb17
distance hb18, (1VFB and chain B and resi 52 and name N), (1VFB and chain B and resi 56 and name O)
hide labels, hb18
set dash_color, yellow, hb18
distance hb19, (1VFB and chain B and resi 52 and name O), (1VFB and chain B and resi 53 and name N)
hide labels, hb19
set dash_color, yellow, hb19
distance hb20, (1VFB and chain B and resi 52 and name O), (1VFB and chain B and resi 54 and name N)
hide labels, hb20
set dash_color, yellow, hb20
distance hb21, (1VFB and chain B and resi 52 and name NE1), (1VFB and chain B and resi 58 and name OD2)
hide labels, hb21
set dash_color, yellow, hb21
distance hb22, (1VFB and chain B and resi 53 and name N), (1VFB and chain B and resi 54 and name N)
hide labels, hb22
set dash_color, yellow, hb22
distance hb23, (1VFB and chain B and resi 54 and name OD1), (1VFB and chain B and resi 56 and name N)
hide labels, hb23
set dash_color, yellow, hb23
hide labels
show sticks, (1VFB and chain B and resi 30) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 30) and name CA
set sphere_scale, 0.20, (1VFB and chain B and resi 30)
show sticks, (1VFB and chain B and resi 31) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 31) and name CA
set sphere_scale, 0.20, (1VFB and chain B and resi 31)
show sticks, (1VFB and chain B and resi 32) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 32) and name CA
set sphere_scale, 0.20, (1VFB and chain B and resi 32)
show sticks, (1VFB and chain B and resi 52) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 52) and name CA
set sphere_scale, 0.20, (1VFB and chain B and resi 52)
show sticks, (1VFB and chain B and resi 53) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 53) and name CA
set sphere_scale, 0.20, (1VFB and chain B and resi 53)
show sticks, (1VFB and chain B and resi 54) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 54) and name CA
set sphere_scale, 0.20, (1VFB and chain B and resi 54)
show sticks, (1VFB and chain B and resi 56) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 56) and name CA
set sphere_scale, 0.20, (1VFB and chain B and resi 56)
show sticks, (1VFB and chain B and resi 58) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 58) and name CA
set sphere_scale, 0.20, (1VFB and chain B and resi 58)
show sticks, (1VFB and chain B and resi 99) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 99) and name CA
set sphere_scale, 0.20, (1VFB and chain B and resi 99)
show sticks, (1VFB and chain C and resi 116) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 116) and name CA
set sphere_scale, 0.20, (1VFB and chain C and resi 116)
pseudoatom lab_B30, pos=[42.976, -18.901, -12.451], label="THR-30"
pseudoatom lab_B31, pos=[45.126, -16.338, -10.518], label="GLY-31"
pseudoatom lab_B32, pos=[43.196, -13.097, -11.238], label="TYR-32"
pseudoatom lab_B52, pos=[24.461, -6.101, -5.919], label="TRP-52"
distance dl_B52, (1VFB and chain B and resi 52 and name CA), lab_B52
hide labels, dl_B52
set dash_color, grey50, dl_B52
pseudoatom lab_B53, pos=[27.03, -17.379, -8.558], label="GLY-53"
distance dl_B53, (1VFB and chain B and resi 53 and name CA), lab_B53
hide labels, dl_B53
set dash_color, grey50, dl_B53
pseudoatom lab_B54, pos=[29.83, -29.003, -9.647], label="ASP-54"
distance dl_B54, (1VFB and chain B and resi 54 and name CA), lab_B54
hide labels, dl_B54
set dash_color, grey50, dl_B54
pseudoatom lab_B56, pos=[32.588, -17.303, -8.598], label="ASN-56"
pseudoatom lab_B58, pos=[29.269, -11.142, -8.468], label="ASP-58"
pseudoatom lab_B99, pos=[44.751, -7.824, -8.85], label="ARG-99"
pseudoatom lab_C116, pos=[55.586, -11.711, -5.549], label="LYS-116"
distance dl_C116, (1VFB and chain C and resi 116 and name CA), lab_C116
hide labels, dl_C116
set dash_color, grey50, dl_C116
set label_color, black, lab_*

set_view (\
   -0.9662228, -0.2134501, -0.1444049,\
   -0.2412813, 0.9461299, 0.2159207,\
   0.0905375, 0.2434696, -0.9656736,\
   0.0000000, 0.0000000, -78.2184372,\
   41.0270195, -14.1246691, -8.9316015,\
   61.6680145, 94.7688599, -20.0000000 )

set ray_trace_fog, 0
ray 1300, 1000
png pymol_groundtruth/figures_composite/1VFB_twostream_inset1.png, dpi=300, ray=0
