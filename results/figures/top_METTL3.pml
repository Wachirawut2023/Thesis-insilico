# Top hit: METTL3
# anchor: chain A Cys375
# binding_mode: monodentate
# covalent_score: 2.0
# functional_proximity_A: 3.45

load /home/runner/work/Thesis-insilico/Thesis-insilico/data/prepared/METTL3.clean.pdb, prot
bg_color white
hide everything
show cartoon, prot
color gray80, prot
set cartoon_transparency, 0.3

# Reactive anchor Cys (yellow)
select anchor_cys, chain A and resi 375 and resn CYS
show sticks, anchor_cys
color yellow, anchor_cys
label first (anchor_cys and name CA), "C%s" % resi

# Functional anchor residues (cyan)
select funcanc_377, resi 377 and (not anchor_cys)
show sticks, funcanc_377
color cyan, funcanc_377
select funcanc_395, resi 395 and (not anchor_cys)
show sticks, funcanc_395
color cyan, funcanc_395
select funcanc_396, resi 396 and (not anchor_cys)
show sticks, funcanc_396
color cyan, funcanc_396
select funcanc_397, resi 397 and (not anchor_cys)
show sticks, funcanc_397
color cyan, funcanc_397
select funcanc_398, resi 398 and (not anchor_cys)
show sticks, funcanc_398
color cyan, funcanc_398
select funcanc_511, resi 511 and (not anchor_cys)
show sticks, funcanc_511
color cyan, funcanc_511
select funcanc_538, resi 538 and (not anchor_cys)
show sticks, funcanc_538
color cyan, funcanc_538
select funcanc_549, resi 549 and (not anchor_cys)
show sticks, funcanc_549
color cyan, funcanc_549

# Camera + render
orient anchor_cys
zoom anchor_cys, 12
set ray_shadows, 0
set ray_opaque_background, 1
set label_size, 14
set label_color, black
ray 1200, 900
png /home/runner/work/Thesis-insilico/Thesis-insilico/results/figures/top_METTL3.png, dpi=150
save /home/runner/work/Thesis-insilico/Thesis-insilico/results/figures/top_METTL3.pse
