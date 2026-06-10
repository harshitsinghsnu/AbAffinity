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
color 0xC4E7BD, chain B
color 0xCBDEF0, chain A
color 0xCECEE5, chain C
show sticks, (1VFB and chain B and resi 30) and (sidechain or name CA)
color 0x80C97F, 1VFB and chain B and resi 30
show sticks, (1VFB and chain B and resi 31) and (sidechain or name CA)
color 0x82CA82, 1VFB and chain B and resi 31
show sticks, (1VFB and chain B and resi 32) and (sidechain or name CA)
color 0x90D18D, 1VFB and chain B and resi 32
show sticks, (1VFB and chain B and resi 52) and (sidechain or name CA)
color 0x319A50, 1VFB and chain B and resi 52
show sticks, (1VFB and chain B and resi 53) and (sidechain or name CA)
color 0x38A256, 1VFB and chain B and resi 53
show sticks, (1VFB and chain B and resi 54) and (sidechain or name CA)
color 0x41AB5D, 1VFB and chain B and resi 54
show sticks, (1VFB and chain B and resi 56) and (sidechain or name CA)
color 0x80C97F, 1VFB and chain B and resi 56
show sticks, (1VFB and chain B and resi 58) and (sidechain or name CA)
color 0x70C274, 1VFB and chain B and resi 58
show sticks, (1VFB and chain B and resi 99) and (sidechain or name CA)
color 0x51B365, 1VFB and chain B and resi 99
show sticks, (1VFB and chain B and resi 100) and (sidechain or name CA)
color 0x43AC5E, 1VFB and chain B and resi 100
show sticks, (1VFB and chain B and resi 101) and (sidechain or name CA)
color 0x2D964D, 1VFB and chain B and resi 101
show sticks, (1VFB and chain B and resi 102) and (sidechain or name CA)
color 0x4EB163, 1VFB and chain B and resi 102
show sticks, (1VFB and chain A and resi 30) and (sidechain or name CA)
color 0x2676B7, 1VFB and chain A and resi 30
show sticks, (1VFB and chain A and resi 32) and (sidechain or name CA)
color 0x8BC0DD, 1VFB and chain A and resi 32
show sticks, (1VFB and chain A and resi 49) and (sidechain or name CA)
color 0x3282BE, 1VFB and chain A and resi 49
show sticks, (1VFB and chain A and resi 50) and (sidechain or name CA)
color 0x2272B5, 1VFB and chain A and resi 50
show sticks, (1VFB and chain A and resi 52) and (sidechain or name CA)
color 0x7BB7D9, 1VFB and chain A and resi 52
show sticks, (1VFB and chain A and resi 53) and (sidechain or name CA)
color 0x519CCC, 1VFB and chain A and resi 53
show sticks, (1VFB and chain A and resi 91) and (sidechain or name CA)
color 0x509BCB, 1VFB and chain A and resi 91
show sticks, (1VFB and chain A and resi 92) and (sidechain or name CA)
color 0x60A6D1, 1VFB and chain A and resi 92
show sticks, (1VFB and chain A and resi 93) and (sidechain or name CA)
color 0x539DCC, 1VFB and chain A and resi 93
show sticks, (1VFB and chain A and resi 94) and (sidechain or name CA)
color 0x3585BF, 1VFB and chain A and resi 94
show sticks, (1VFB and chain C and resi 18) and (sidechain or name CA)
color 0x776BB0, 1VFB and chain C and resi 18
show sticks, (1VFB and chain C and resi 19) and (sidechain or name CA)
color 0x928EC2, 1VFB and chain C and resi 19
show sticks, (1VFB and chain C and resi 22) and (sidechain or name CA)
color 0x7C75B5, 1VFB and chain C and resi 22
show sticks, (1VFB and chain C and resi 24) and (sidechain or name CA)
color 0x664AA0, 1VFB and chain C and resi 24
show sticks, (1VFB and chain C and resi 26) and (sidechain or name CA)
color 0x583093, 1VFB and chain C and resi 26
show sticks, (1VFB and chain C and resi 102) and (sidechain or name CA)
color 0x7A71B3, 1VFB and chain C and resi 102
show sticks, (1VFB and chain C and resi 116) and (sidechain or name CA)
color 0x7364AD, 1VFB and chain C and resi 116
show sticks, (1VFB and chain C and resi 121) and (sidechain or name CA)
color 0x7D77B7, 1VFB and chain C and resi 121
show sticks, (1VFB and chain C and resi 129) and (sidechain or name CA)
color 0x63449C, 1VFB and chain C and resi 129

distance hb0, (1VFB and chain A and resi 30 and name O), (1VFB and chain A and resi 32 and name N)
hide labels, hb0
set dash_color, yellow, hb0
distance hb1, (1VFB and chain A and resi 30 and name ND1), (1VFB and chain C and resi 129 and name OXT)
hide labels, hb1
set dash_color, yellow, hb1
distance hb2, (1VFB and chain A and resi 32 and name O), (1VFB and chain A and resi 91 and name N)
hide labels, hb2
set dash_color, yellow, hb2
distance hb3, (1VFB and chain A and resi 49 and name N), (1VFB and chain A and resi 53 and name O)
hide labels, hb3
set dash_color, yellow, hb3
distance hb4, (1VFB and chain A and resi 49 and name O), (1VFB and chain A and resi 50 and name N)
hide labels, hb4
set dash_color, yellow, hb4
distance hb5, (1VFB and chain A and resi 49 and name O), (1VFB and chain A and resi 52 and name N)
hide labels, hb5
set dash_color, yellow, hb5
distance hb6, (1VFB and chain A and resi 50 and name O), (1VFB and chain A and resi 52 and name N)
hide labels, hb6
set dash_color, yellow, hb6
distance hb7, (1VFB and chain A and resi 50 and name OH), (1VFB and chain C and resi 18 and name OD1)
hide labels, hb7
set dash_color, yellow, hb7
distance hb8, (1VFB and chain A and resi 52 and name N), (1VFB and chain A and resi 53 and name N)
hide labels, hb8
set dash_color, yellow, hb8
distance hb9, (1VFB and chain A and resi 53 and name OG1), (1VFB and chain C and resi 19 and name ND2)
hide labels, hb9
set dash_color, yellow, hb9
distance hb10, (1VFB and chain A and resi 91 and name N), (1VFB and chain A and resi 92 and name N)
hide labels, hb10
set dash_color, yellow, hb10
distance hb11, (1VFB and chain A and resi 91 and name O), (1VFB and chain C and resi 121 and name NE2)
hide labels, hb11
set dash_color, yellow, hb11
distance hb12, (1VFB and chain A and resi 92 and name N), (1VFB and chain A and resi 93 and name N)
hide labels, hb12
set dash_color, yellow, hb12
distance hb13, (1VFB and chain A and resi 93 and name N), (1VFB and chain C and resi 121 and name OE1)
hide labels, hb13
set dash_color, yellow, hb13
distance hb14, (1VFB and chain A and resi 93 and name O), (1VFB and chain A and resi 94 and name N)
hide labels, hb14
set dash_color, yellow, hb14
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
distance hb24, (1VFB and chain B and resi 99 and name N), (1VFB and chain B and resi 100 and name N)
hide labels, hb24
set dash_color, yellow, hb24
distance hb25, (1VFB and chain B and resi 99 and name N), (1VFB and chain B and resi 102 and name O)
hide labels, hb25
set dash_color, yellow, hb25
distance hb26, (1VFB and chain B and resi 99 and name O), (1VFB and chain B and resi 101 and name N)
hide labels, hb26
set dash_color, yellow, hb26
distance hb27, (1VFB and chain B and resi 99 and name NH2), (1VFB and chain C and resi 102 and name O)
hide labels, hb27
set dash_color, yellow, hb27
distance hb28, (1VFB and chain B and resi 100 and name N), (1VFB and chain B and resi 101 and name N)
hide labels, hb28
set dash_color, yellow, hb28
distance hb29, (1VFB and chain B and resi 100 and name O), (1VFB and chain B and resi 102 and name NH1)
hide labels, hb29
set dash_color, yellow, hb29
distance hb30, (1VFB and chain B and resi 100 and name OD1), (1VFB and chain C and resi 24 and name OG)
hide labels, hb30
set dash_color, yellow, hb30
distance hb31, (1VFB and chain B and resi 100 and name OD2), (1VFB and chain C and resi 22 and name O)
hide labels, hb31
set dash_color, yellow, hb31
distance hb32, (1VFB and chain B and resi 101 and name N), (1VFB and chain B and resi 102 and name N)
hide labels, hb32
set dash_color, yellow, hb32
distance hb33, (1VFB and chain B and resi 101 and name OH), (1VFB and chain C and resi 121 and name N)
hide labels, hb33
set dash_color, yellow, hb33
distance hb34, (1VFB and chain B and resi 102 and name NH1), (1VFB and chain C and resi 22 and name O)
hide labels, hb34
set dash_color, yellow, hb34
distance hb35, (1VFB and chain C and resi 18 and name N), (1VFB and chain C and resi 19 and name N)
hide labels, hb35
set dash_color, yellow, hb35
distance hb36, (1VFB and chain C and resi 19 and name O), (1VFB and chain C and resi 22 and name N)
hide labels, hb36
set dash_color, yellow, hb36
distance hb37, (1VFB and chain C and resi 24 and name O), (1VFB and chain C and resi 26 and name N)
hide labels, hb37
set dash_color, yellow, hb37
set label_color, black

set_view (\
   0.5467612, 0.3717805, -0.7502210,\
   -0.5598075, 0.8286186, 0.0026435,\
   0.6226297, 0.4185339, 0.6611819,\
   0.0000000, 0.0000000, -225.3664856,\
   42.8503799, -5.8561935, -5.8587456,\
   177.6806641, 273.0523071, -20.0000000 )

set ray_trace_fog, 0
ray 2400, 2000
png pymol_groundtruth/figures_composite/1VFB_twostream_overview.png, dpi=300, ray=0
