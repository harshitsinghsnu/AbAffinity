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

set_view (\
   0.8574986, 0.5099646, -0.0680603,\
   -0.5031731, 0.8588616, 0.0957787,\
   0.1072981, -0.0478840, 0.9930732,\
   0.0000000, 0.0000000, -149.3691711,\
   -11.1778145, 3.5799732, 36.5161438,\
   117.7637939, 180.9745483, -20.0000000 )

set ray_trace_fog, 0
ray 2800, 2100
png pymol_groundtruth/figures/4ETQ_twostream_nolabel.png, dpi=350, ray=0
save pymol_groundtruth/figures/4ETQ_twostream_nolabel.pse
