# Top hit: VIRMA
# anchor: chain A Cys1008
# binding_mode: bidentate
# covalent_score: 2.0
# functional_proximity_A: NA

load /home/runner/work/Thesis-insilico/Thesis-insilico/data/prepared/VIRMA.clean.pdb, prot
bg_color white
hide everything
show cartoon, prot
color gray80, prot
set cartoon_transparency, 0.3

# Reactive anchor Cys (yellow)
select anchor_cys, chain A and resi 1008 and resn CYS
show sticks, anchor_cys
color yellow, anchor_cys
label first (anchor_cys and name CA), "C%s" % resi

# Bidentate/tridentate partner Cys (orange)
select partner_A938, chain A and resi 938 and resn CYS
show sticks, partner_A938
color orange, partner_A938
label first (partner_A938 and name CA), "C%s" % resi

# Camera + render
orient anchor_cys
zoom anchor_cys, 12
set ray_shadows, 0
set ray_opaque_background, 1
set label_size, 14
set label_color, black
ray 1200, 900
png /home/runner/work/Thesis-insilico/Thesis-insilico/results/figures/top_VIRMA.png, dpi=150
save /home/runner/work/Thesis-insilico/Thesis-insilico/results/figures/top_VIRMA.pse
