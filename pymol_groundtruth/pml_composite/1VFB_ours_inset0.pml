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
color 0x005321, 1VFB and chain B and resi 30
color 0x248C45, 1VFB and chain B and resi 31
color 0x077331, 1VFB and chain B and resi 32
color 0x3DA75A, 1VFB and chain B and resi 52
color 0x5BB769, 1VFB and chain B and resi 53
color 0x1C8540, 1VFB and chain B and resi 54
color 0x2C954C, 1VFB and chain B and resi 56
color 0x279048, 1VFB and chain B and resi 58
color 0x005B24, 1VFB and chain B and resi 99
color 0x005823, 1VFB and chain B and resi 100
color 0x005823, 1VFB and chain B and resi 101
color 0x17803C, 1VFB and chain B and resi 102
color 0x2171B5, 1VFB and chain A and resi 30
color 0x1E6DB2, 1VFB and chain A and resi 32
color 0x1967AD, 1VFB and chain A and resi 49
color 0x1967AD, 1VFB and chain A and resi 50
color 0x1765AB, 1VFB and chain A and resi 52
color 0x2272B5, 1VFB and chain A and resi 53
color 0x3D8DC3, 1VFB and chain A and resi 91
color 0x3A8AC1, 1VFB and chain A and resi 92
color 0x1866AC, 1VFB and chain A and resi 93
color 0x09539D, 1VFB and chain A and resi 94
color 0x7F7CB9, 1VFB and chain C and resi 18
color 0x7E79B8, 1VFB and chain C and resi 19
color 0x7C75B5, 1VFB and chain C and resi 22
color 0x776CB1, 1VFB and chain C and resi 24
color 0x6A52A3, 1VFB and chain C and resi 26
color 0x776BB0, 1VFB and chain C and resi 102
color 0x8F8BC1, 1VFB and chain C and resi 116
color 0x8380BB, 1VFB and chain C and resi 121
color 0x9B97C6, 1VFB and chain C and resi 129

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
distance hb8, (1VFB and chain A and resi 52 and name N), (1VFB and chain A and resi 53 and name N)
hide labels, hb8
set dash_color, yellow, hb8
distance hb10, (1VFB and chain A and resi 91 and name N), (1VFB and chain A and resi 92 and name N)
hide labels, hb10
set dash_color, yellow, hb10
distance hb12, (1VFB and chain A and resi 92 and name N), (1VFB and chain A and resi 93 and name N)
hide labels, hb12
set dash_color, yellow, hb12
distance hb14, (1VFB and chain A and resi 93 and name O), (1VFB and chain A and resi 94 and name N)
hide labels, hb14
set dash_color, yellow, hb14
distance hb32, (1VFB and chain B and resi 101 and name N), (1VFB and chain B and resi 102 and name N)
hide labels, hb32
set dash_color, yellow, hb32
hide labels
show sticks, (1VFB and chain B and resi 101) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 101) and name CA
set sphere_scale, 0.20, (1VFB and chain B and resi 101)
show sticks, (1VFB and chain B and resi 102) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 102) and name CA
set sphere_scale, 0.20, (1VFB and chain B and resi 102)
show sticks, (1VFB and chain A and resi 30) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 30) and name CA
set sphere_scale, 0.20, (1VFB and chain A and resi 30)
show sticks, (1VFB and chain A and resi 32) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 32) and name CA
set sphere_scale, 0.20, (1VFB and chain A and resi 32)
show sticks, (1VFB and chain A and resi 49) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 49) and name CA
set sphere_scale, 0.20, (1VFB and chain A and resi 49)
show sticks, (1VFB and chain A and resi 50) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 50) and name CA
set sphere_scale, 0.20, (1VFB and chain A and resi 50)
show sticks, (1VFB and chain A and resi 52) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 52) and name CA
set sphere_scale, 0.20, (1VFB and chain A and resi 52)
show sticks, (1VFB and chain A and resi 53) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 53) and name CA
set sphere_scale, 0.20, (1VFB and chain A and resi 53)
show sticks, (1VFB and chain A and resi 91) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 91) and name CA
set sphere_scale, 0.20, (1VFB and chain A and resi 91)
show sticks, (1VFB and chain A and resi 92) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 92) and name CA
set sphere_scale, 0.20, (1VFB and chain A and resi 92)
show sticks, (1VFB and chain A and resi 93) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 93) and name CA
set sphere_scale, 0.20, (1VFB and chain A and resi 93)
show sticks, (1VFB and chain A and resi 94) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 94) and name CA
set sphere_scale, 0.20, (1VFB and chain A and resi 94)
show sticks, (1VFB and chain C and resi 129) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 129) and name CA
set sphere_scale, 0.20, (1VFB and chain C and resi 129)
pseudoatom lab_B101, pos=[51.296, -3.265, -17.723], label="TYR-101"
distance dl_B101, (1VFB and chain B and resi 101 and name CA), lab_B101
hide labels, dl_B101
set dash_color, grey50, dl_B101
pseudoatom lab_B102, pos=[42.533, -2.714, -7.815], label="ARG-102"
pseudoatom lab_A30, pos=[37.711, 0.535, 5.256], label="HIS-30"
pseudoatom lab_A32, pos=[29.695, -0.724, 16.719], label="TYR-32"
distance dl_A32, (1VFB and chain A and resi 32 and name CA), lab_A32
hide labels, dl_A32
set dash_color, grey50, dl_A32
pseudoatom lab_A49, pos=[44.909, 1.972, -4.669], label="TYR-49"
pseudoatom lab_A50, pos=[44.11, 1.219, -0.969], label="TYR-50"
pseudoatom lab_A52, pos=[47.317, 6.079, -0.485], label="THR-52"
pseudoatom lab_A53, pos=[56.568, 9.41, 4.396], label="THR-53"
distance dl_A53, (1VFB and chain A and resi 53 and name CA), lab_A53
hide labels, dl_A53
set dash_color, grey50, dl_A53
pseudoatom lab_A91, pos=[36.4, -2.398, -3.25], label="PHE-91"
pseudoatom lab_A92, pos=[26.539, -6.092, 5.392], label="TRP-92"
distance dl_A92, (1VFB and chain A and resi 92 and name CA), lab_A92
hide labels, dl_A92
set dash_color, grey50, dl_A92
pseudoatom lab_A93, pos=[23.945, -12.508, -5.646], label="SER-93"
distance dl_A93, (1VFB and chain A and resi 93 and name CA), lab_A93
hide labels, dl_A93
set dash_color, grey50, dl_A93
pseudoatom lab_A94, pos=[30.324, -6.341, -3.696], label="THR-94"
pseudoatom lab_C129, pos=[39.649, -3.082, 12.4], label="LEU-129"
set label_color, black, lab_*

set_view (\
   -0.8488405, 0.2610745, -0.4596846,\
   -0.3453813, 0.3844351, 0.8561084,\
   0.4002270, 0.8854659, -0.2361538,\
   0.0000000, 0.0000000, -85.1851883,\
   40.6093559, -2.2060947, -0.4434462,\
   67.1606522, 103.2097244, -20.0000000 )

set ray_trace_fog, 0
ray 1300, 1000
png pymol_groundtruth/figures_composite/1VFB_ours_inset0.png, dpi=300, ray=0
