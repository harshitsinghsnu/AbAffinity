# 4ETQ -- Two-stream -- GT binding residues coloured by IG attribution
reinitialize
load pymol_groundtruth/pdb/pdb4etq.ent, 4ETQ
remove not polymer
remove 4ETQ and not chain H+L+C
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
color 0xC4E7BD, chain H
color 0xCBDEF0, chain L
color 0xCECEE5, chain C
show sticks, (4ETQ and chain H and resi 29) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 29) and name CA
color 0x258D46, 4ETQ and chain H and resi 29
show sticks, (4ETQ and chain H and resi 31) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 31) and name CA
color 0x40AA5C, 4ETQ and chain H and resi 31
show sticks, (4ETQ and chain H and resi 32) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 32) and name CA
color 0x2E974E, 4ETQ and chain H and resi 32
show sticks, (4ETQ and chain H and resi 34) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 34) and name CA
color 0x50B264, 4ETQ and chain H and resi 34
show sticks, (4ETQ and chain H and resi 51) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 51) and name CA
color 0x56B567, 4ETQ and chain H and resi 51
show sticks, (4ETQ and chain H and resi 53) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 53) and name CA
color 0x90D18D, 4ETQ and chain H and resi 53
show sticks, (4ETQ and chain H and resi 55) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 55) and name CA
color 0x5EB96B, 4ETQ and chain H and resi 55
show sticks, (4ETQ and chain H and resi 56) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 56) and name CA
color 0x6BBF71, 4ETQ and chain H and resi 56
show sticks, (4ETQ and chain H and resi 58) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 58) and name CA
color 0x50B264, 4ETQ and chain H and resi 58
show sticks, (4ETQ and chain H and resi 60) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 60) and name CA
color 0x4CB062, 4ETQ and chain H and resi 60
show sticks, (4ETQ and chain H and resi 75) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 75) and name CA
color 0x329B51, 4ETQ and chain H and resi 75
show sticks, (4ETQ and chain H and resi 102) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 102) and name CA
color 0x61BA6C, 4ETQ and chain H and resi 102
show sticks, (4ETQ and chain H and resi 103) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 103) and name CA
color 0x6CC072, 4ETQ and chain H and resi 103
show sticks, (4ETQ and chain H and resi 105) and (sidechain or name CA)
show spheres, (4ETQ and chain H and resi 105) and name CA
color 0x73C375, 4ETQ and chain H and resi 105
show sticks, (4ETQ and chain L and resi 31) and (sidechain or name CA)
show spheres, (4ETQ and chain L and resi 31) and name CA
color 0x4B98C9, 4ETQ and chain L and resi 31
show sticks, (4ETQ and chain L and resi 91) and (sidechain or name CA)
show spheres, (4ETQ and chain L and resi 91) and name CA
color 0x82BADB, 4ETQ and chain L and resi 91
show sticks, (4ETQ and chain L and resi 92) and (sidechain or name CA)
show spheres, (4ETQ and chain L and resi 92) and name CA
color 0x88BEDC, 4ETQ and chain L and resi 92
show sticks, (4ETQ and chain L and resi 94) and (sidechain or name CA)
show spheres, (4ETQ and chain L and resi 94) and name CA
color 0x65AAD3, 4ETQ and chain L and resi 94
show sticks, (4ETQ and chain C and resi 3) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 3) and name CA
color 0x7466AE, 4ETQ and chain C and resi 3
show sticks, (4ETQ and chain C and resi 5) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 5) and name CA
color 0x8C89BF, 4ETQ and chain C and resi 5
show sticks, (4ETQ and chain C and resi 39) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 39) and name CA
color 0x5F3D99, 4ETQ and chain C and resi 39
show sticks, (4ETQ and chain C and resi 40) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 40) and name CA
color 0x7F7CB9, 4ETQ and chain C and resi 40
show sticks, (4ETQ and chain C and resi 41) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 41) and name CA
color 0x674CA0, 4ETQ and chain C and resi 41
show sticks, (4ETQ and chain C and resi 44) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 44) and name CA
color 0x5E3A98, 4ETQ and chain C and resi 44
show sticks, (4ETQ and chain C and resi 108) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 108) and name CA
color 0x5B3595, 4ETQ and chain C and resi 108
show sticks, (4ETQ and chain C and resi 145) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 145) and name CA
color 0x674CA0, 4ETQ and chain C and resi 145
show sticks, (4ETQ and chain C and resi 174) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 174) and name CA
color 0x7160AB, 4ETQ and chain C and resi 174
show sticks, (4ETQ and chain C and resi 175) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 175) and name CA
color 0x6950A2, 4ETQ and chain C and resi 175
show sticks, (4ETQ and chain C and resi 176) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 176) and name CA
color 0x6C55A5, 4ETQ and chain C and resi 176
show sticks, (4ETQ and chain C and resi 177) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 177) and name CA
color 0x796FB3, 4ETQ and chain C and resi 177
show sticks, (4ETQ and chain C and resi 204) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 204) and name CA
color 0x776BB0, 4ETQ and chain C and resi 204
show sticks, (4ETQ and chain C and resi 205) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 205) and name CA
color 0x7567AE, 4ETQ and chain C and resi 205
show sticks, (4ETQ and chain C and resi 215) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 215) and name CA
color 0x6950A2, 4ETQ and chain C and resi 215
show sticks, (4ETQ and chain C and resi 217) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 217) and name CA
color 0x7364AD, 4ETQ and chain C and resi 217
show sticks, (4ETQ and chain C and resi 219) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 219) and name CA
color 0x684EA1, 4ETQ and chain C and resi 219
show sticks, (4ETQ and chain C and resi 220) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 220) and name CA
color 0x715FAA, 4ETQ and chain C and resi 220
show sticks, (4ETQ and chain C and resi 221) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 221) and name CA
color 0x53258E, 4ETQ and chain C and resi 221
show sticks, (4ETQ and chain C and resi 223) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 223) and name CA
color 0x7364AD, 4ETQ and chain C and resi 223
show sticks, (4ETQ and chain C and resi 224) and (sidechain or name CA)
show spheres, (4ETQ and chain C and resi 224) and name CA
color 0x562C91, 4ETQ and chain C and resi 224
select interface, 4ETQ and ((chain H and resi 29+31+32+34+51+53+55+56+58+60+75+102+103+105) or (chain L and resi 31+91+92+94) or (chain C and resi 3+5+39+40+41+44+108+145+174+175+176+177+204+205+215+217+219+220+221+223+224))
deselect
set light_count, 2
set specular, 0.15
set ambient, 0.45

pseudoatom lab_H29, pos=[-22.122, 4.466, 28.784], label="SER29"
pseudoatom lab_H31, pos=[-18.681, 8.097, 31.586], label="ASN31"
pseudoatom lab_H32, pos=[-17.908, 4.528, 32.707], label="PHE32"
pseudoatom lab_H34, pos=[-11.384, 5.307, 30.564], label="TRP34"
pseudoatom lab_H51, pos=[-5.758, 8.391, 29.324], label="MET51"
pseudoatom lab_H53, pos=[-11.47, 9.249, 33.127], label="ASP53"
pseudoatom lab_H55, pos=[-28.286, 29.316, 32.098], label="SER55"
distance dl_H55, (4ETQ and chain H and resi 55 and name CA), lab_H55
hide labels, dl_H55
pseudoatom lab_H56, pos=[-14.615, 13.675, 36.059], label="GLU56"
pseudoatom lab_H58, pos=[-23.697, 37.141, 32.117], label="GLU58"
distance dl_H58, (4ETQ and chain H and resi 58 and name CA), lab_H58
hide labels, dl_H58
pseudoatom lab_H60, pos=[-3.022, 10.98, 33.253], label="ARG60"
pseudoatom lab_H75, pos=[-25.586, 32.658, 26.197], label="ARG75"
distance dl_H75, (4ETQ and chain H and resi 75 and name CA), lab_H75
hide labels, dl_H75
pseudoatom lab_H102, pos=[-13.888, -0.766, 33.284], label="TYR102"
pseudoatom lab_H103, pos=[1.514, -30.224, 38.397], label="ARG103"
distance dl_H103, (4ETQ and chain H and resi 103 and name CA), lab_H103
hide labels, dl_H103
pseudoatom lab_H105, pos=[8.68, -18.884, 34.302], label="ASP105"
distance dl_H105, (4ETQ and chain H and resi 105 and name CA), lab_H105
hide labels, dl_H105
pseudoatom lab_L31, pos=[2.295, -8.855, 36.641], label="SER31"
pseudoatom lab_L91, pos=[-0.223, -2.826, 32.72], label="TRP91"
pseudoatom lab_L92, pos=[13.068, -10.775, 37.26], label="THR92"
distance dl_L92, (4ETQ and chain L and resi 92 and name CA), lab_L92
hide labels, dl_L92
pseudoatom lab_L94, pos=[2.499, 3.913, 34.187], label="TYR94"
pseudoatom lab_C3, pos=[21.894, 5.372, 42.447], label="GLN3"
distance dl_C3, (4ETQ and chain C and resi 3 and name CA), lab_C3
hide labels, dl_C3
pseudoatom lab_C5, pos=[17.294, -2.438, 42.582], label="LEU5"
distance dl_C5, (4ETQ and chain C and resi 5 and name CA), lab_C5
hide labels, dl_C5
pseudoatom lab_C39, pos=[-33.183, 21.923, 36.562], label="THR39"
distance dl_C39, (4ETQ and chain C and resi 39 and name CA), lab_C39
hide labels, dl_C39
pseudoatom lab_C40, pos=[-35.338, 17.815, 34.528], label="GLY40"
distance dl_C40, (4ETQ and chain C and resi 40 and name CA), lab_C40
hide labels, dl_C40
pseudoatom lab_C41, pos=[-37.729, 14.037, 35.915], label="LYS41"
distance dl_C41, (4ETQ and chain C and resi 41 and name CA), lab_C41
hide labels, dl_C41
pseudoatom lab_C44, pos=[-31.343, 26.475, 43.209], label="ARG44"
distance dl_C44, (4ETQ and chain C and resi 44 and name CA), lab_C44
hide labels, dl_C44
pseudoatom lab_C108, pos=[-1.976, 15.416, 42.673], label="LYS108"
pseudoatom lab_C145, pos=[-23.661, 13.981, 35.942], label="ASN145"
pseudoatom lab_C174, pos=[14.854, -6.146, 44.693], label="ILE174"
distance dl_C174, (4ETQ and chain C and resi 174 and name CA), lab_C174
hide labels, dl_C174
pseudoatom lab_C175, pos=[19.696, 1.323, 41.025], label="ASN175"
distance dl_C175, (4ETQ and chain C and resi 175 and name CA), lab_C175
hide labels, dl_C175
pseudoatom lab_C176, pos=[-5.21, 5.561, 38.353], label="HIS176"
pseudoatom lab_C177, pos=[24.39, 9.002, 39.523], label="SER177"
distance dl_C177, (4ETQ and chain C and resi 177 and name CA), lab_C177
hide labels, dl_C177
pseudoatom lab_C204, pos=[-42.143, 5.965, 33.341], label="SER204"
distance dl_C204, (4ETQ and chain C and resi 204 and name CA), lab_C204
hide labels, dl_C204
pseudoatom lab_C205, pos=[-39.614, 9.548, 29.927], label="LEU205"
distance dl_C205, (4ETQ and chain C and resi 205 and name CA), lab_C205
hide labels, dl_C205
pseudoatom lab_C215, pos=[-25.004, 1.482, 37.758], label="ILE215"
pseudoatom lab_C217, pos=[-47.068, -1.387, 38.218], label="GLU217"
distance dl_C217, (4ETQ and chain C and resi 217 and name CA), lab_C217
hide labels, dl_C217
pseudoatom lab_C219, pos=[-44.868, 2.659, 39.613], label="TYR219"
distance dl_C219, (4ETQ and chain C and resi 219 and name CA), lab_C219
hide labels, dl_C219
pseudoatom lab_C220, pos=[3.439, -25.792, 43.788], label="ARG220"
distance dl_C220, (4ETQ and chain C and resi 220 and name CA), lab_C220
hide labels, dl_C220
pseudoatom lab_C221, pos=[-10.248, -5.691, 40.628], label="ASN221"
pseudoatom lab_C223, pos=[10.382, -14.137, 42.962], label="TYR223"
distance dl_C223, (4ETQ and chain C and resi 223 and name CA), lab_C223
hide labels, dl_C223
pseudoatom lab_C224, pos=[5.728, -21.871, 43.882], label="LYS224"
distance dl_C224, (4ETQ and chain C and resi 224 and name CA), lab_C224
hide labels, dl_C224
set label_color, black, lab_*

set_view (\
   0.8574986, 0.5099646, -0.0680603,\
   -0.5031731, 0.8588616, 0.0957787,\
   0.1072981, -0.0478840, 0.9930732,\
   0.0000000, 0.0000000, -149.3691711,\
   -11.1778145, 3.5799732, 36.5161438,\
   117.7637939, 180.9745483, -20.0000000 )

set ray_trace_fog, 0
ray 2800, 2100
png pymol_groundtruth/figures/4ETQ_twostream.png, dpi=350, ray=0
save pymol_groundtruth/figures/4ETQ_twostream.pse
