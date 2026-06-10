reinitialize
load pymol_groundtruth/pdb/pdb5grj.ent, 5GRJ
remove not polymer
remove 5GRJ and not chain H+L+A
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
color 0x2E974E, 5GRJ and chain H and resi 27
color 0x39A357, 5GRJ and chain H and resi 28
color 0x278F48, 5GRJ and chain H and resi 31
color 0x2B944B, 5GRJ and chain H and resi 33
color 0x56B567, 5GRJ and chain H and resi 52
color 0x41AB5D, 5GRJ and chain H and resi 53
color 0x38A256, 5GRJ and chain H and resi 54
color 0x50B264, 5GRJ and chain H and resi 57
color 0x53B365, 5GRJ and chain H and resi 59
color 0x59B769, 5GRJ and chain H and resi 101
color 0x7EC87E, 5GRJ and chain H and resi 102
color 0x70C274, 5GRJ and chain H and resi 103
color 0x76C477, 5GRJ and chain H and resi 104
color 0x7AC67B, 5GRJ and chain H and resi 105
color 0x78C67A, 5GRJ and chain H and resi 106
color 0x4694C7, 5GRJ and chain L and resi 32
color 0x509BCB, 5GRJ and chain L and resi 34
color 0x5EA5D1, 5GRJ and chain L and resi 93
color 0x85BCDB, 5GRJ and chain L and resi 97
color 0x7BB7D9, 5GRJ and chain L and resi 99
color 0x6BAED6, 5GRJ and chain L and resi 95
color 0x5A3394, 5GRJ and chain A and resi 54
color 0x8E8AC0, 5GRJ and chain A and resi 60
color 0x6950A2, 5GRJ and chain A and resi 61
color 0x684DA1, 5GRJ and chain A and resi 62
color 0x5A3495, 5GRJ and chain A and resi 63
color 0x562B90, 5GRJ and chain A and resi 66
color 0x6B54A4, 5GRJ and chain A and resi 56
color 0x6C56A6, 5GRJ and chain A and resi 58
color 0x5C3696, 5GRJ and chain A and resi 59
color 0x7F7CB9, 5GRJ and chain A and resi 69
color 0x5A3394, 5GRJ and chain A and resi 73
color 0x7363AC, 5GRJ and chain A and resi 75
color 0x66499F, 5GRJ and chain A and resi 113
color 0x7C75B5, 5GRJ and chain A and resi 121
color 0x6E59A7, 5GRJ and chain A and resi 123
color 0x5B3595, 5GRJ and chain A and resi 125
color 0x61409B, 5GRJ and chain A and resi 115
color 0x664AA0, 5GRJ and chain A and resi 117
color 0x7E79B8, 5GRJ and chain A and resi 78

distance hb4, (5GRJ and chain A and resi 58 and name O), (5GRJ and chain A and resi 59 and name N)
hide labels, hb4
set dash_color, yellow, hb4
distance hb6, (5GRJ and chain A and resi 59 and name N), (5GRJ and chain A and resi 60 and name N)
hide labels, hb6
set dash_color, yellow, hb6
distance hb7, (5GRJ and chain A and resi 59 and name N), (5GRJ and chain A and resi 62 and name O)
hide labels, hb7
set dash_color, yellow, hb7
distance hb8, (5GRJ and chain A and resi 59 and name O), (5GRJ and chain A and resi 61 and name N)
hide labels, hb8
set dash_color, yellow, hb8
distance hb9, (5GRJ and chain A and resi 60 and name N), (5GRJ and chain A and resi 61 and name N)
hide labels, hb9
set dash_color, yellow, hb9
distance hb10, (5GRJ and chain A and resi 60 and name O), (5GRJ and chain A and resi 62 and name N)
hide labels, hb10
set dash_color, yellow, hb10
distance hb11, (5GRJ and chain A and resi 60 and name OE1), (5GRJ and chain L and resi 34 and name OH)
hide labels, hb11
set dash_color, yellow, hb11
distance hb12, (5GRJ and chain A and resi 61 and name N), (5GRJ and chain A and resi 62 and name N)
hide labels, hb12
set dash_color, yellow, hb12
distance hb13, (5GRJ and chain A and resi 61 and name O), (5GRJ and chain H and resi 104 and name N)
hide labels, hb13
set dash_color, yellow, hb13
distance hb14, (5GRJ and chain A and resi 61 and name OD1), (5GRJ and chain H and resi 105 and name N)
hide labels, hb14
set dash_color, yellow, hb14
distance hb15, (5GRJ and chain A and resi 61 and name OD2), (5GRJ and chain H and resi 106 and name OG1)
hide labels, hb15
set dash_color, yellow, hb15
distance hb16, (5GRJ and chain A and resi 61 and name OD2), (5GRJ and chain L and resi 99 and name NH2)
hide labels, hb16
set dash_color, yellow, hb16
distance hb17, (5GRJ and chain A and resi 62 and name N), (5GRJ and chain A and resi 63 and name N)
hide labels, hb17
set dash_color, yellow, hb17
distance hb33, (5GRJ and chain H and resi 103 and name O), (5GRJ and chain H and resi 104 and name N)
hide labels, hb33
set dash_color, yellow, hb33
distance hb34, (5GRJ and chain H and resi 103 and name OG1), (5GRJ and chain H and resi 106 and name N)
hide labels, hb34
set dash_color, yellow, hb34
distance hb35, (5GRJ and chain H and resi 104 and name N), (5GRJ and chain H and resi 105 and name N)
hide labels, hb35
set dash_color, yellow, hb35
distance hb36, (5GRJ and chain H and resi 105 and name N), (5GRJ and chain H and resi 106 and name N)
hide labels, hb36
set dash_color, yellow, hb36
distance hb37, (5GRJ and chain H and resi 105 and name O), (5GRJ and chain L and resi 99 and name NE)
hide labels, hb37
set dash_color, yellow, hb37
distance hb38, (5GRJ and chain H and resi 106 and name OG1), (5GRJ and chain L and resi 99 and name NH2)
hide labels, hb38
set dash_color, yellow, hb38
distance hb39, (5GRJ and chain L and resi 32 and name O), (5GRJ and chain L and resi 34 and name N)
hide labels, hb39
set dash_color, yellow, hb39
distance hb40, (5GRJ and chain L and resi 32 and name OH), (5GRJ and chain L and resi 95 and name N)
hide labels, hb40
set dash_color, yellow, hb40
distance hb41, (5GRJ and chain L and resi 34 and name O), (5GRJ and chain L and resi 93 and name N)
hide labels, hb41
set dash_color, yellow, hb41
distance hb42, (5GRJ and chain L and resi 93 and name OH), (5GRJ and chain L and resi 97 and name OG)
hide labels, hb42
set dash_color, yellow, hb42
distance hb43, (5GRJ and chain L and resi 95 and name O), (5GRJ and chain L and resi 97 and name N)
hide labels, hb43
set dash_color, yellow, hb43
hide labels
show sticks, (5GRJ and chain H and resi 103) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 103) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 103)
show sticks, (5GRJ and chain H and resi 104) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 104) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 104)
show sticks, (5GRJ and chain H and resi 105) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 105) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 105)
show sticks, (5GRJ and chain H and resi 106) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 106) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 106)
show sticks, (5GRJ and chain L and resi 32) and (sidechain or name CA)
show spheres, (5GRJ and chain L and resi 32) and name CA
set sphere_scale, 0.20, (5GRJ and chain L and resi 32)
show sticks, (5GRJ and chain L and resi 34) and (sidechain or name CA)
show spheres, (5GRJ and chain L and resi 34) and name CA
set sphere_scale, 0.20, (5GRJ and chain L and resi 34)
show sticks, (5GRJ and chain L and resi 93) and (sidechain or name CA)
show spheres, (5GRJ and chain L and resi 93) and name CA
set sphere_scale, 0.20, (5GRJ and chain L and resi 93)
show sticks, (5GRJ and chain L and resi 97) and (sidechain or name CA)
show spheres, (5GRJ and chain L and resi 97) and name CA
set sphere_scale, 0.20, (5GRJ and chain L and resi 97)
show sticks, (5GRJ and chain L and resi 99) and (sidechain or name CA)
show spheres, (5GRJ and chain L and resi 99) and name CA
set sphere_scale, 0.20, (5GRJ and chain L and resi 99)
show sticks, (5GRJ and chain L and resi 95) and (sidechain or name CA)
show spheres, (5GRJ and chain L and resi 95) and name CA
set sphere_scale, 0.20, (5GRJ and chain L and resi 95)
show sticks, (5GRJ and chain A and resi 60) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 60) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 60)
show sticks, (5GRJ and chain A and resi 61) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 61) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 61)
show sticks, (5GRJ and chain A and resi 62) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 62) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 62)
show sticks, (5GRJ and chain A and resi 63) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 63) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 63)
show sticks, (5GRJ and chain A and resi 58) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 58) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 58)
show sticks, (5GRJ and chain A and resi 59) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 59) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 59)
pseudoatom lab_H103, pos=[20.146, -80.644, 26.349], label="THR-103"
pseudoatom lab_H104, pos=[29.455, -79.15, 32.999], label="VAL-104"
distance dl_H104, (5GRJ and chain H and resi 104 and name CA), lab_H104
hide labels, dl_H104
set dash_color, grey50, dl_H104
pseudoatom lab_H105, pos=[12.112, -65.357, 20.008], label="THR-105"
distance dl_H105, (5GRJ and chain H and resi 105 and name CA), lab_H105
hide labels, dl_H105
set dash_color, grey50, dl_H105
pseudoatom lab_H106, pos=[17.755, -75.588, 25.698], label="THR-106"
pseudoatom lab_L32, pos=[10.938, -82.49, 17.284], label="TYR-32"
pseudoatom lab_L34, pos=[5.825, -76.148, 16.992], label="TYR-34"
distance dl_L34, (5GRJ and chain L and resi 34 and name CA), lab_L34
hide labels, dl_L34
set dash_color, grey50, dl_L34
pseudoatom lab_L93, pos=[15.273, -77.012, 18.349], label="TYR-93"
pseudoatom lab_L97, pos=[20.87, -77.828, 15.363], label="SER-97"
pseudoatom lab_L99, pos=[16.999, -72.905, 18.658], label="ARG-99"
pseudoatom lab_L95, pos=[12.047, -80.071, 4.454], label="SER-95"
distance dl_L95, (5GRJ and chain L and resi 95 and name CA), lab_L95
hide labels, dl_L95
set dash_color, grey50, dl_L95
pseudoatom lab_A60, pos=[26.667, -95.36, 19.518], label="GLU-60"
distance dl_A60, (5GRJ and chain A and resi 60 and name CA), lab_A60
hide labels, dl_A60
set dash_color, grey50, dl_A60
pseudoatom lab_A61, pos=[27.038, -85.371, 29.638], label="ASP-61"
distance dl_A61, (5GRJ and chain A and resi 61 and name CA), lab_A61
hide labels, dl_A61
set dash_color, grey50, dl_A61
pseudoatom lab_A62, pos=[25.462, -91.129, 25.636], label="LYS-62"
distance dl_A62, (5GRJ and chain A and resi 62 and name CA), lab_A62
hide labels, dl_A62
set dash_color, grey50, dl_A62
pseudoatom lab_A63, pos=[23.852, -84.609, 24.74], label="ASN-63"
pseudoatom lab_A58, pos=[24.901, -85.692, 20.604], label="GLU-58"
pseudoatom lab_A59, pos=[21.576, -86.693, 19.166], label="MET-59"
set label_color, black, lab_*

set_view (\
   0.6776001, 0.0954967, -0.7292040,\
   -0.5356704, 0.7434639, -0.4003982,\
   0.5039001, 0.6619228, 0.5549260,\
   0.0000000, 0.0000000, -70.0373383,\
   18.5009518, -80.8461685, 20.4394722,\
   55.2179718, 84.8567047, -20.0000000 )

set ray_trace_fog, 0
ray 1300, 1000
png pymol_groundtruth/figures_composite/5GRJ_twostream_inset2.png, dpi=300, ray=0
