# 1VFB -- Ours (All-CDR) -- GT binding residues coloured by IG attribution
reinitialize
load pymol_groundtruth/pdb/pdb1vfb.ent, 1VFB
remove not polymer
remove 1VFB and not chain B+A+C
hide everything
bg_color white
set ray_opaque_background, 1
set ray_shadows, 0
set ray_trace_mode, 1
set ray_trace_color, grey50
set antialias, 2
set cartoon_fancy_helices, 1
set cartoon_transparency, 0.25
set cartoon_side_chain_helper, 1
set stick_radius, 0.30
set sphere_scale, 0.34
set valence, 0
set label_size, 15
set label_font_id, 7
set label_color, black
set label_outline_color, white
set label_bg_color, white
set label_bg_transparency, 0.25
set float_labels, 1
set dash_color, grey30
set dash_width, 1.6
set dash_gap, 0.40
set dash_length, 0.40
set dash_radius, 0.045
set dash_round_ends, 0
show cartoon
color 0xC4E7BD, chain B
color 0xCBDEF0, chain A
color 0xCECEE5, chain C
show sticks, (1VFB and chain B and resi 30) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 30) and name CA
color 0x005321, 1VFB and chain B and resi 30
show sticks, (1VFB and chain B and resi 31) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 31) and name CA
color 0x248C45, 1VFB and chain B and resi 31
show sticks, (1VFB and chain B and resi 32) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 32) and name CA
color 0x077331, 1VFB and chain B and resi 32
show sticks, (1VFB and chain B and resi 52) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 52) and name CA
color 0x3DA75A, 1VFB and chain B and resi 52
show sticks, (1VFB and chain B and resi 53) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 53) and name CA
color 0x5BB769, 1VFB and chain B and resi 53
show sticks, (1VFB and chain B and resi 54) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 54) and name CA
color 0x1C8540, 1VFB and chain B and resi 54
show sticks, (1VFB and chain B and resi 56) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 56) and name CA
color 0x2C954C, 1VFB and chain B and resi 56
show sticks, (1VFB and chain B and resi 58) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 58) and name CA
color 0x279048, 1VFB and chain B and resi 58
show sticks, (1VFB and chain B and resi 99) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 99) and name CA
color 0x005B24, 1VFB and chain B and resi 99
show sticks, (1VFB and chain B and resi 100) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 100) and name CA
color 0x005823, 1VFB and chain B and resi 100
show sticks, (1VFB and chain B and resi 101) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 101) and name CA
color 0x005823, 1VFB and chain B and resi 101
show sticks, (1VFB and chain B and resi 102) and (sidechain or name CA)
show spheres, (1VFB and chain B and resi 102) and name CA
color 0x17803C, 1VFB and chain B and resi 102
show sticks, (1VFB and chain A and resi 30) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 30) and name CA
color 0x2171B5, 1VFB and chain A and resi 30
show sticks, (1VFB and chain A and resi 32) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 32) and name CA
color 0x1E6DB2, 1VFB and chain A and resi 32
show sticks, (1VFB and chain A and resi 49) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 49) and name CA
color 0x1967AD, 1VFB and chain A and resi 49
show sticks, (1VFB and chain A and resi 50) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 50) and name CA
color 0x1967AD, 1VFB and chain A and resi 50
show sticks, (1VFB and chain A and resi 52) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 52) and name CA
color 0x1765AB, 1VFB and chain A and resi 52
show sticks, (1VFB and chain A and resi 53) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 53) and name CA
color 0x2272B5, 1VFB and chain A and resi 53
show sticks, (1VFB and chain A and resi 91) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 91) and name CA
color 0x3D8DC3, 1VFB and chain A and resi 91
show sticks, (1VFB and chain A and resi 92) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 92) and name CA
color 0x3A8AC1, 1VFB and chain A and resi 92
show sticks, (1VFB and chain A and resi 93) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 93) and name CA
color 0x1866AC, 1VFB and chain A and resi 93
show sticks, (1VFB and chain A and resi 94) and (sidechain or name CA)
show spheres, (1VFB and chain A and resi 94) and name CA
color 0x09539D, 1VFB and chain A and resi 94
show sticks, (1VFB and chain C and resi 18) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 18) and name CA
color 0x7F7CB9, 1VFB and chain C and resi 18
show sticks, (1VFB and chain C and resi 19) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 19) and name CA
color 0x7E79B8, 1VFB and chain C and resi 19
show sticks, (1VFB and chain C and resi 22) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 22) and name CA
color 0x7C75B5, 1VFB and chain C and resi 22
show sticks, (1VFB and chain C and resi 24) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 24) and name CA
color 0x776CB1, 1VFB and chain C and resi 24
show sticks, (1VFB and chain C and resi 26) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 26) and name CA
color 0x6A52A3, 1VFB and chain C and resi 26
show sticks, (1VFB and chain C and resi 102) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 102) and name CA
color 0x776BB0, 1VFB and chain C and resi 102
show sticks, (1VFB and chain C and resi 116) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 116) and name CA
color 0x8F8BC1, 1VFB and chain C and resi 116
show sticks, (1VFB and chain C and resi 121) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 121) and name CA
color 0x8380BB, 1VFB and chain C and resi 121
show sticks, (1VFB and chain C and resi 129) and (sidechain or name CA)
show spheres, (1VFB and chain C and resi 129) and name CA
color 0x9B97C6, 1VFB and chain C and resi 129
select interface, 1VFB and ((chain B and resi 30+31+32+52+53+54+56+58+99+100+101+102) or (chain A and resi 30+32+49+50+52+53+91+92+93+94) or (chain C and resi 18+19+22+24+26+102+116+121+129))
deselect
set light_count, 2
set specular, 0.15
set ambient, 0.45

pseudoatom lab_B30, pos=[42.976, -18.901, -12.451], label="THR30"
pseudoatom lab_B31, pos=[55.787, -28.261, -25.806], label="GLY31"
distance dl_B31, (1VFB and chain B and resi 31 and name CA), lab_B31
hide labels, dl_B31
pseudoatom lab_B32, pos=[46.903, -26.333, -24.743], label="TYR32"
distance dl_B32, (1VFB and chain B and resi 32 and name CA), lab_B32
hide labels, dl_B32
pseudoatom lab_B52, pos=[37.013, -14.326, -9.635], label="TRP52"
pseudoatom lab_B53, pos=[40.013, -16.626, -10.331], label="GLY53"
pseudoatom lab_B54, pos=[32.332, -29.232, -15.119], label="ASP54"
distance dl_B54, (1VFB and chain B and resi 54 and name CA), lab_B54
hide labels, dl_B54
pseudoatom lab_B56, pos=[32.588, -17.303, -8.598], label="ASN56"
pseudoatom lab_B58, pos=[29.269, -11.142, -8.468], label="ASP58"
pseudoatom lab_B99, pos=[44.751, -7.824, -8.85], label="ARG99"
pseudoatom lab_B100, pos=[44.16, -7.7, -5.079], label="ASP100"
pseudoatom lab_B101, pos=[41.403, -5.08, -5.045], label="TYR101"
pseudoatom lab_B102, pos=[51.309, 17.412, 6.571], label="ARG102"
distance dl_B102, (1VFB and chain B and resi 102 and name CA), lab_B102
hide labels, dl_B102
pseudoatom lab_A30, pos=[37.711, 0.535, 5.256], label="HIS30"
pseudoatom lab_A32, pos=[32.351, 15.361, 16.856], label="TYR32"
distance dl_A32, (1VFB and chain A and resi 32 and name CA), lab_A32
hide labels, dl_A32
pseudoatom lab_A49, pos=[56.251, 17.496, 4.476], label="TYR49"
distance dl_A49, (1VFB and chain A and resi 49 and name CA), lab_A49
hide labels, dl_A49
pseudoatom lab_A50, pos=[44.11, 1.219, -0.969], label="TYR50"
pseudoatom lab_A52, pos=[61.53, 16.868, 3.171], label="THR52"
distance dl_A52, (1VFB and chain A and resi 52 and name CA), lab_A52
hide labels, dl_A52
pseudoatom lab_A53, pos=[49.488, 3.702, -2.517], label="THR53"
pseudoatom lab_A91, pos=[36.4, -2.398, -3.25], label="PHE91"
pseudoatom lab_A92, pos=[34.588, -3.248, 0.032], label="TRP92"
pseudoatom lab_A93, pos=[15.206, -26.734, -11.486], label="SER93"
distance dl_A93, (1VFB and chain A and resi 93 and name CA), lab_A93
hide labels, dl_A93
pseudoatom lab_A94, pos=[30.324, -6.341, -3.696], label="THR94"
pseudoatom lab_C18, pos=[50.193, 9.247, 17.637], label="ASP18"
distance dl_C18, (1VFB and chain C and resi 18 and name CA), lab_C18
hide labels, dl_C18
pseudoatom lab_C19, pos=[49.179, -3.643, 2.437], label="ASN19"
pseudoatom lab_C22, pos=[68.461, 12.754, 5.737], label="GLY22"
distance dl_C22, (1VFB and chain C and resi 22 and name CA), lab_C22
hide labels, dl_C22
pseudoatom lab_C24, pos=[44.434, 10.886, 17.819], label="SER24"
distance dl_C24, (1VFB and chain C and resi 24 and name CA), lab_C24
hide labels, dl_C24
pseudoatom lab_C26, pos=[40.571, 8.525, 22.443], label="GLY26"
distance dl_C26, (1VFB and chain C and resi 26 and name CA), lab_C26
hide labels, dl_C26
pseudoatom lab_C102, pos=[54.946, -7.513, -4.755], label="GLY102"
pseudoatom lab_C116, pos=[41.079, -30.869, -16.504], label="LYS116"
distance dl_C116, (1VFB and chain C and resi 116 and name CA), lab_C116
hide labels, dl_C116
pseudoatom lab_C121, pos=[24.94, -30.456, -10.557], label="GLN121"
distance dl_C121, (1VFB and chain C and resi 121 and name CA), lab_C121
hide labels, dl_C121
pseudoatom lab_C129, pos=[30.663, 8.404, 26.58], label="LEU129"
distance dl_C129, (1VFB and chain C and resi 129 and name CA), lab_C129
hide labels, dl_C129
set label_color, black, lab_*

set_view (\
   -0.2379093, -0.9230943, -0.3021523,\
   -0.7699513, -0.0103998, 0.6380178,\
   -0.5920929, 0.3844330, -0.7082636,\
   -0.0000000, 0.0000000, -141.6154175,\
   41.7829208, -6.8829947, -2.9936199,\
   111.6506729, 171.5801697, -20.0000000 )

set ray_trace_fog, 0
ray 2800, 2100
png pymol_groundtruth/figures/1VFB_ours.png, dpi=350, ray=0
save pymol_groundtruth/figures/1VFB_ours.pse
