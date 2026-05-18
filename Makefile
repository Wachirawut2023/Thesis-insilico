.PHONY: help env ligand fetch prep cys pockets dock-nc dock-cov rank all smoke clean

PY := python
SCRIPTS := scripts

help:
	@echo "Targets:"
	@echo "  env       create conda env from env/environment.yml"
	@echo "  ligand    prepare As(OH)3 PDBQT"
	@echo "  fetch     download PDB + AlphaFold structures"
	@echo "  prep      clean, protonate, PDBQT-ify receptors"
	@echo "  cys       Cys reactivity landscape (Tier 1)"
	@echo "  pockets   fpocket + Cys-cluster docking boxes"
	@echo "  dock-nc   non-covalent Vina docking"
	@echo "  dock-cov  geometric covalent docking (primary readout)"
	@echo "  rank      composite ranking + heatmap + report"
	@echo "  all       fetch -> prep -> cys -> pockets -> dock-nc -> dock-cov -> rank"
	@echo "  smoke     run pipeline on a 3-target subset (METTL3, ALKBH5, YTHDF2)"
	@echo "  clean     remove derived data (keeps source PDBs)"

env:
	conda env create -f env/environment.yml || conda env update -f env/environment.yml

ligand:
	bash ligands/prepare_as3.sh

fetch:
	$(PY) $(SCRIPTS)/01_fetch_structures.py

prep:
	$(PY) $(SCRIPTS)/02_prepare_receptors.py

cys:
	$(PY) $(SCRIPTS)/03_cys_landscape.py

pockets:
	$(PY) $(SCRIPTS)/04_pocket_detection.py

dock-nc:
	$(PY) $(SCRIPTS)/05_dock_noncovalent.py

dock-cov:
	$(PY) $(SCRIPTS)/06_dock_covalent.py

rank:
	$(PY) $(SCRIPTS)/07_score_and_rank.py

all: ligand fetch prep cys pockets dock-nc dock-cov rank

smoke:
	SMOKE_GENES="METTL3,ALKBH5,YTHDF2" $(PY) $(SCRIPTS)/01_fetch_structures.py
	SMOKE_GENES="METTL3,ALKBH5,YTHDF2" $(PY) $(SCRIPTS)/02_prepare_receptors.py
	$(PY) $(SCRIPTS)/03_cys_landscape.py
	$(PY) $(SCRIPTS)/04_pocket_detection.py
	$(PY) $(SCRIPTS)/05_dock_noncovalent.py
	$(PY) $(SCRIPTS)/06_dock_covalent.py
	$(PY) $(SCRIPTS)/07_score_and_rank.py

clean:
	rm -rf data/prepared/* results/cys_table.tsv results/pocket_table.tsv \
	       results/docking_boxes.tsv results/docking_noncovalent.tsv \
	       results/docking_covalent.tsv results/composite_ranking.tsv \
	       results/protein_summary.tsv results/REPORT.md \
	       results/figures/*.png results/docking_poses
