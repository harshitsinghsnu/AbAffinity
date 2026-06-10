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

set_view (\
   -0.2379093, -0.9230943, -0.3021523,\
   -0.7699513, -0.0103998, 0.6380178,\
   -0.5920929, 0.3844330, -0.7082636,\
   -0.0000000, 0.0000000, -141.6154175,\
   41.7829208, -6.8829947, -2.9936199,\
   111.6506729, 171.5801697, -20.0000000 )

set ray_trace_fog, 0
ray 2800, 2100
png pymol_groundtruth/figures/1VFB_ours_nolabel.png, dpi=350, ray=0
save pymol_groundtruth/figures/1VFB_ours_nolabel.pse
