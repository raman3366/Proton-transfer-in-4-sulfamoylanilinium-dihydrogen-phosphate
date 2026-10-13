# PyMOL: pymol overlay.pml (run in this folder)
load receptor.pdb
load crystal_ligand.sdf
hide everything
show cartoon, receptor
set cartoon_transparency, 0.6
show sticks, crystal_ligand
color grey70, crystal_ligand and elem C
load docked_production_seed03.sdf, docked_production_seed03
show sticks, docked_production_seed03
color orange, docked_production_seed03 and elem C
show sticks, receptor and byres (receptor within 4.5 of crystal_ligand) and not name N+C+O
orient crystal_ligand
zoom crystal_ligand, 6
