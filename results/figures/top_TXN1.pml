# Top hit: TXN1
# anchor: chain A Cys32
# binding_mode: bidentate
# covalent_score: 4.0
# functional_proximity_A: 0.00

load /home/runner/work/Thesis-insilico/Thesis-insilico/data/prepared/TXN1.clean.pdb, prot
bg_color white
hide everything
show cartoon, prot
color gray80, prot
set cartoon_transparency, 0.3

# Reactive anchor Cys (yellow)
select anchor_cys, chain A and resi 32 and resn CYS
show sticks, anchor_cys
color yellow, anchor_cys
label first (anchor_cys and name CA), "C%s" % resi

# Bidentate/tridentate partner Cys (orange)
select partner_A35, chain A and resi 35 and resn CYS
show sticks, partner_A35
color orange, partner_A35
label first (partner_A35 and name CA), "C%s" % resi

# Functional anchor residues (cyan)
select funcanc_32, resi 32 and (not anchor_cys)
show sticks, funcanc_32
color cyan, funcanc_32
select funcanc_35, resi 35 and (not anchor_cys)
show sticks, funcanc_35
color cyan, funcanc_35

# Camera + render
orient anchor_cys
zoom anchor_cys, 12
set ray_shadows, 0
set ray_opaque_background, 1
set label_size, 14
set label_color, black
ray 1200, 900
png /home/runner/work/Thesis-insilico/Thesis-insilico/results/figures/top_TXN1.png, dpi=150
save /home/runner/work/Thesis-insilico/Thesis-insilico/results/figures/top_TXN1.pse
