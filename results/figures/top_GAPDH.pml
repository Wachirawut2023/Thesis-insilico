# Top hit: GAPDH
# anchor: chain O Cys152
# binding_mode: monodentate
# covalent_score: 2.0
# functional_proximity_A: 0.00

load /home/runner/work/Thesis-insilico/Thesis-insilico/data/prepared/GAPDH.clean.pdb, prot
bg_color white
hide everything
show cartoon, prot
color gray80, prot
set cartoon_transparency, 0.3

# Reactive anchor Cys (yellow)
select anchor_cys, chain O and resi 152 and resn CYS
show sticks, anchor_cys
color yellow, anchor_cys
label first (anchor_cys and name CA), "C%s" % resi

# Functional anchor residues (cyan)
select funcanc_152, resi 152 and (not anchor_cys)
show sticks, funcanc_152
color cyan, funcanc_152

# Camera + render
orient anchor_cys
zoom anchor_cys, 12
set ray_shadows, 0
set ray_opaque_background, 1
set label_size, 14
set label_color, black
ray 1200, 900
png /home/runner/work/Thesis-insilico/Thesis-insilico/results/figures/top_GAPDH.png, dpi=150
save /home/runner/work/Thesis-insilico/Thesis-insilico/results/figures/top_GAPDH.pse
