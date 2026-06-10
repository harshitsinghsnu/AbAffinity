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

distance hb23, (5GRJ and chain H and resi 27 and name O), (5GRJ and chain H and resi 28 and name N)
hide labels, hb23
set dash_color, yellow, hb23
distance hb24, (5GRJ and chain H and resi 28 and name O), (5GRJ and chain H and resi 31 and name N)
hide labels, hb24
set dash_color, yellow, hb24
distance hb29, (5GRJ and chain H and resi 101 and name O), (5GRJ and chain H and resi 102 and name N)
hide labels, hb29
set dash_color, yellow, hb29
hide labels
show sticks, (5GRJ and chain H and resi 27) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 27) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 27)
show sticks, (5GRJ and chain H and resi 28) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 28) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 28)
show sticks, (5GRJ and chain H and resi 31) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 31) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 31)
show sticks, (5GRJ and chain H and resi 33) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 33) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 33)
show sticks, (5GRJ and chain H and resi 101) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 101) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 101)
show sticks, (5GRJ and chain H and resi 102) and (sidechain or name CA)
show spheres, (5GRJ and chain H and resi 102) and name CA
set sphere_scale, 0.20, (5GRJ and chain H and resi 102)
show sticks, (5GRJ and chain A and resi 73) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 73) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 73)
show sticks, (5GRJ and chain A and resi 75) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 75) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 75)
show sticks, (5GRJ and chain A and resi 78) and (sidechain or name CA)
show spheres, (5GRJ and chain A and resi 78) and name CA
set sphere_scale, 0.20, (5GRJ and chain A and resi 78)
pseudoatom lab_H27, pos=[26.997, -74.746, 38.999], label="PHE-27"
pseudoatom lab_H28, pos=[34.068, -81.724, 45.098], label="THR-28"
distance dl_H28, (5GRJ and chain H and resi 28 and name CA), lab_H28
hide labels, dl_H28
set dash_color, grey50, dl_H28
pseudoatom lab_H31, pos=[33.303, -67.651, 32.385], label="SER-31"
distance dl_H31, (5GRJ and chain H and resi 31 and name CA), lab_H31
hide labels, dl_H31
set dash_color, grey50, dl_H31
pseudoatom lab_H33, pos=[27.331, -75.146, 26.959], label="ILE-33"
pseudoatom lab_H101, pos=[16.439, -79.665, 21.891], label="LEU-101"
distance dl_H101, (5GRJ and chain H and resi 101 and name CA), lab_H101
hide labels, dl_H101
set dash_color, grey50, dl_H101
pseudoatom lab_H102, pos=[21.922, -83.02, 28.733], label="GLY-102"
pseudoatom lab_A73, pos=[31.638, -86.376, 32.168], label="ASP-73"
pseudoatom lab_A75, pos=[22.123, -94.805, 31.28], label="LYS-75"
distance dl_A75, (5GRJ and chain A and resi 75 and name CA), lab_A75
hide labels, dl_A75
set dash_color, grey50, dl_A75
pseudoatom lab_A78, pos=[21.197, -87.573, 29.085], label="HIS-78"
set label_color, black, lab_*

set_view (\
   0.5385166, -0.2121987, -0.8154581,\
   0.5807762, 0.7946427, 0.1767541,\
   0.6104907, -0.5687836, 0.5511681,\
   0.0000000, 0.0000000, -61.2534180,\
   26.2576294, -80.9123993, 32.8160706,\
   48.2926636, 74.2141724, -20.0000000 )

set ray_trace_fog, 0
ray 1300, 1000
png pymol_groundtruth/figures_composite/5GRJ_twostream_inset0.png, dpi=300, ray=0
