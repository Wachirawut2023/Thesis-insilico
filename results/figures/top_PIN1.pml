# Top hit: PIN1
# anchor: chain A Cys113
# binding_mode: monodentate
# covalent_score: 2.0
# functional_proximity_A: 0.00

load /home/runner/work/Thesis-insilico/Thesis-insilico/data/prepared/PIN1.clean.pdb, prot
bg_color white
hide everything
show cartoon, prot
color gray80, prot
set cartoon_transparency, 0.3

# Reactive anchor Cys (yellow)
select anchor_cys, chain A and resi 113 and resn CYS
show sticks, anchor_cys
color yellow, anchor_cys
label first (anchor_cys and name CA), "C%s" % resi

# Functional anchor residues (cyan)
select funcanc_113, resi 113 and (not anchor_cys)
show sticks, funcanc_113
color cyan, funcanc_113

# Camera + render
orient anchor_cys
zoom anchor_cys, 12
set ray_shadows, 0
set ray_opaque_background, 1
set label_size, 14
set label_color, black
ray 1200, 900
png /home/runner/work/Thesis-insilico/Thesis-insilico/results/figures/top_PIN1.png, dpi=150
save /home/runner/work/Thesis-insilico/Thesis-insilico/results/figures/top_PIN1.pse
