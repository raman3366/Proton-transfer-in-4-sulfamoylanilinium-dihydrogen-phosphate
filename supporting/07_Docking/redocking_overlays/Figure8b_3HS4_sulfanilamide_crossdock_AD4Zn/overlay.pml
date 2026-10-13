# PyMOL: pymol overlay.pml (run in this folder)
load receptor.pdb
load crystal_ligand.sdf
hide everything
show cartoon, receptor
set cartoon_transparency, 0.6
show sticks, crystal_ligand
color grey70, crystal_ligand and elem C
load docked_validation_seed02.pdbqt, docked_validation_seed02
show sticks, docked_validation_seed02
color orange, docked_validation_seed02 and elem C
load docked_production_seed10.sdf, docked_production_seed10
show sticks, docked_production_seed10
color violet, docked_production_seed10 and elem C
show sticks, receptor and byres (receptor within 4.5 of crystal_ligand) and not name N+C+O
orient crystal_ligand
zoom crystal_ligand, 6
