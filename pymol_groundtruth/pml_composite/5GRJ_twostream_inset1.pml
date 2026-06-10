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

distance hb0, (5GRJ and chain A and resi 54 and name N), (5GRJ and chain A and resi 117 and name O)
hide labels, hb0
set dash_color, yellow, hb0
distance hb1, (5GRJ and chain A and resi 56 and name N), (5GRJ and chain A and resi 115 and name O)
hide labels, hb1
set dash_color, yellow, hb1
distance hb2, (5GRJ and chain A and resi 56 and name OH), (5GRJ and chain H and resi 54 and name N)
hide labels, hb2
set dash_color, yellow, hb2
distance hb22, (5GRJ and chain A and resi 113 and name NH2), (5GRJ and chain H and resi 52 and name OH)
hide labels, hb22
set dash_color, yellow, hb22
distance hb25, (5GRJ and chain H and resi 52 and name N), (5GRJ and chain H and resi 57 and name O)
hide labels, hb25
set dash_color, yellow, hb25
distance hb26, (5GRJ and chain H and resi 52 and name O), (5GRJ and chain H and resi 53 and name N)
hide labels, hb26
set dash_color, yellow, hb26
distance hb27, (5GRJ and chain H and resi 52 and name O), (5GRJ and chain H and resi 54 and name N)
hide labels, hb27
set dash_color, yellow, hb27
distance hb28, (5GRJ and chain H and resi 53 and name N), (5GRJ and chain H and resi 54 and name N)
hide labels, hb28
set dash_color, yellow, hb28
hide labels
show sticks, (5GRJ and chain H and resi 52) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 52) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 52)
show sticks, (5GRJ and chain H and resi 53) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 53) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 53)
show sticks, (5GRJ and chain H and resi 54) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 54) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 54)
show sticks, (5GRJ and chain H and resi 57) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 57) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 57)
show sticks, (5GRJ and chain H and resi 59) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 59) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 59)
show sticks, (5GRJ and chain A and resi 54) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 54) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 54)
show sticks, (5GRJ and chain A and resi 66) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 66) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 66)
show sticks, (5GRJ and chain A and resi 56) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 56) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 56)
show sticks, (5GRJ and chain A and resi 69) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 69) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 69)
show sticks, (5GRJ and chain A and resi 113) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 113) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 113)
show sticks, (5GRJ and chain A and resi 121) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 121) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 121)
show sticks, (5GRJ and chain A and resi 123) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 123) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 123)
show sticks, (5GRJ and chain A and resi 125) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 125) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 125)
show sticks, (5GRJ and chain A and resi 115) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 115) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 115)
show sticks, (5GRJ and chain A and resi 117) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 117) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 117)
pseudoatom lab_H52, pos=[31.091, -75.729, 23.976], label="TYR-52"
pseudoatom lab_H53, pos=[35.732, -70.36, 33.855], label="PRO-53"
distance dl_H53, (5GRJ and chain H and resi 53 and name CA), lab_H53
hide labels, dl_H53
set dash_color, grey50, dl_H53
pseudoatom lab_H54, pos=[39.483, -76.901, 31.237], label="SER-54"
distance dl_H54, (5GRJ and chain H and resi 54 and name CA), lab_H54
hide labels, dl_H54
set dash_color, grey50, dl_H54
pseudoatom lab_H57, pos=[40.854, -72.96, 30.377], label="ILE-57"
distance dl_H57, (5GRJ and chain H and resi 57 and name CA), lab_H57
hide labels, dl_H57
set dash_color, grey50, dl_H57
pseudoatom lab_H59, pos=[28.94, -71.46, 16.64], label="PHE-59"
pseudoatom lab_A54, pos=[41.832, -83.732, 29.555], label="ILE-54"
distance dl_A54, (5GRJ and chain A and resi 54 and name CA), lab_A54
hide labels, dl_A54
set dash_color, grey50, dl_A54
pseudoatom lab_A66, pos=[36.11, -88.571, 33.322], label="GLN-66"
distance dl_A66, (5GRJ and chain A and resi 66 and name CA), lab_A66
hide labels, dl_A66
set dash_color, grey50, dl_A66
pseudoatom lab_A56, pos=[38.804, -91.672, 31.466], label="TYR-56"
distance dl_A56, (5GRJ and chain A and resi 56 and name CA), lab_A56
hide labels, dl_A56
set dash_color, grey50, dl_A56
pseudoatom lab_A69, pos=[40.729, -84.927, 29.324], label="HIS-69"
pseudoatom lab_A113, pos=[24.158, -91.348, 9.746], label="ARG-113"
distance dl_A113, (5GRJ and chain A and resi 113 and name CA), lab_A113
hide labels, dl_A113
set dash_color, grey50, dl_A113
pseudoatom lab_A121, pos=[38.171, -79.42, 15.546], label="ALA-121"
pseudoatom lab_A123, pos=[32.257, -82.593, 13.386], label="TYR-123"
pseudoatom lab_A125, pos=[26.693, -85.905, 11.519], label="ARG-125"
pseudoatom lab_A115, pos=[24.523, -69.326, 9.839], label="MET-115"
distance dl_A115, (5GRJ and chain A and resi 115 and name CA), lab_A115
hide labels, dl_A115
set dash_color, grey50, dl_A115
pseudoatom lab_A117, pos=[44.096, -79.607, 28.099], label="SER-117"
distance dl_A117, (5GRJ and chain A and resi 117 and name CA), lab_A117
hide labels, dl_A117
set dash_color, grey50, dl_A117
set label_color, black, lab_*

set_view (\
   0.5572788, -0.1334560, -0.8195303,\
   -0.0127383, 0.9855086, -0.1691466,\
   0.8302277, 0.1047013, 0.5475030,\
   0.0000000, 0.0000000, -71.4320068,\
   32.1758957, -80.3933945, 20.3103237,\
   56.3175392, 86.5464783, -20.0000000 )

set ray_trace_fog, 0
ray 1300, 1000
png pymol_groundtruth/figures_composite/5GRJ_twostream_inset1.png, dpi=300, ray=0
