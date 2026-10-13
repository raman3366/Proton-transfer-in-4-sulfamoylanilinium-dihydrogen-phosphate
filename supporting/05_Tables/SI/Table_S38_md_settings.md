| Quantity | Value |
|:------------------------------|:-------------------------------------------------------------------------------------------|
| *System* |  |
| Structure | PDB 6RL9, chain A: 258 residues (His3 to Lys261, no residue 126), acetyl cap on His3; Glu106 as a carboxylate, His64 neutral (HIE) |
| Ligand | sulfonamidate anion of sulfanilamide (18 atoms), N2 bonded to the zinc; complex charge −1 |
| Start X | crystal pose (d1 3.90 Å, d2 7.74 Å) |
| Start N | docked near-native pose, 0.61 Å from the crystal ligand (d1 3.36 Å, d2 7.73 Å) |
| Start R | docked ring-rotated pose, 3.57 Å from the crystal ligand (d1 9.58 Å, d2 2.60 Å) |
| Atoms (X / N / R) | 40,902 / 40,911 / 40,857 |
| Water molecules (X / N / R) | 12,247 / 12,250 / 12,232 |
| Crystal waters kept (X / N / R) | 215 / 215 / 211 of 216 sites |
| Ions | 34 Na⁺ and 33 Cl⁻ (0.15 M counted in the solvent volume) |
| Cell | truncated octahedron, 10 Å buffer, edge 84.1 Å as built |
| *Force field* |  |
| Protein, water, ions | ff14SB, TIP3P, Joung–Cheatham |
| Ligand | GAFF2 2.2.30; equilibrium values of the bond S–N2 and of the angles C4–S–N2 and S–N2–H and the torsion about S–N2 (a series on H–N2–S–C4, H–N2–S–O set to zero) refitted to the DFT of this work |
| Zinc site | bonded model (MCPB.py, AmberTools 26.0): 4 Zn–N bonds, 14 angles, 39 torsions with zero force constant |
| Zinc set 2e (nine runs) | equilibrium values of 6RL9 (Zn–N2 1.931 Å; Zn–N(His) 2.010, 2.021, 2.018 Å); force constants 122.1 and 78.0, 73.3, 74.3 kcal mol⁻¹ Å⁻², 35 (N–Zn–N) and 50 kcal mol⁻¹ rad⁻² (Zn–N–C) |
| Zinc set 2s (one run) | Seminario force constants and the equilibrium values of the optimised model (Zn–N2 1.9665 Å, 102.4 kcal mol⁻¹ Å⁻²) |
| Terms at N2, as built | angles Zn–N2–S (50 kcal mol⁻¹ rad⁻², 117.55°, the angle of 6RL9), Zn–N2–H (50, 107.4°, the angle of the constructed hydrogen in start X) and S–N2–H (61.55, 118.1°), with the torsion H–N2–S–C4 of the ligand; the same in both zinc sets |
| Charges of the zinc site | RESP on the 101-atom capped model, B3LYP/6-31G* (pure d and f functions, density fitting); Zn +0.34 e |
| *Dynamics* |  |
| Program | OpenMM 8.6.1, CUDA, mixed precision |
| Non-bonded | particle-mesh Ewald, cut-off 1.0 nm, tolerance 0.0005, dispersion correction |
| Constraints, masses | bonds to hydrogen constrained, rigid water; hydrogen mass 3.024 u |
| Thermostat, barostat | Langevin (middle scheme), 300 K, 1 ps⁻¹; Monte Carlo, 1 bar |
| Minimisation | heavy atoms of the solute and crystal-water oxygens held, 1,000 kJ mol⁻¹ nm⁻² |
| Stage 1, heating | 100 ps, 2 fs, NVT, 50 to 300 K; heavy atoms of the solute held, 1,000 kJ mol⁻¹ nm⁻² |
| Stage 2 | 300 ps, 2 fs, NPT; heavy atoms of the solute held, 1,000 kJ mol⁻¹ nm⁻² |
| Stage 3 | 300 ps, 2 fs, NPT; heavy atoms of the solute held, 250 kJ mol⁻¹ nm⁻²; extended by 500 ps if no water lies in the place of the crystal ring (start R) |
| Stage 4, ligand released | 400 ps, 4 fs, NPT; protein backbone held, 250 kJ mol⁻¹ nm⁻² |
| Stage 5 | 400 ps, 4 fs, NPT; protein backbone held, 50 kJ mol⁻¹ nm⁻² |
| Stage 6 | 500 ps, 4 fs, NPT; nothing held |
| Production | 100 ns, NPT, 4 fs, no restraint and no bias |
| Frames written | solute every 10 ps (10,000 per run), all atoms every 100 ps (1,000 per run) |
| Runs | 9 with set 2e (3 per start), 2 sensitivity runs from start R, 2 validation runs of 5 ns; 1,110 ns of production in all |
| *Analysis* |  |
| Orientation | d1, N1 to Phe131 CZ; d2, N1 to Asn67 OD1. Crystal-like: d1 < 6.7 Å and d2 > 5.2 Å; rotated-like: d1 > 6.7 Å and d2 < 5.2 Å; intermediate otherwise |
| Verdict of a run | CONVERTED or STAYED if at least 90 % of the frames of the last 20 ns are in the state opposite to the start or in that of the start; MIXED otherwise |
| Time of a conversion | first frame of the first visit of at least 1 ns in the opposite state, counted from the release of the ligand |
| Ligand RMSD | 11 heavy atoms against the crystal pose after the fit of the 258 Cα atoms; smallest of the four equivalent labellings |
| Hydrogen bond to Thr199 | N2···OG1 at most 3.5 Å and N2–H···OG1 at least 135° |
| Place of a crystal water | occupied in an all-atom frame when a water oxygen lies within 1.4 Å of the crystal position after the fit of the Cα atoms |
| His64 deep in the site | Zn–NE2(His64) at most 7.0 Å with χ1 between 0° and 120°; the limit was chosen after looking at the data |
| Nearly planar nitrogen | signed volume of the unit vectors from N2 to Zn, H and S below 0.2 in size |

Table: Settings of the molecular dynamics of the carbonic anhydrase II–sulfanilamide complex. {#tbl:s-md-settings}

N1 is the aniline nitrogen and N2 the sulfonamide nitrogen of the ligand. Restraints are harmonic position restraints on the coordinates of the built system (*U* = *k* *d*²/2). The force constants of bonds and angles are those of the Amber energy function, *U* = *k*(*x* − *x*~0~)², which has no factor ½. The equilibrium Zn–N distances and N–Zn–N angles of set 2e are those of 6RL9, so the zinc geometry of the nine runs is put in and is no result. The force-field terms at N2 failed two quantum-chemical tests (@tbl:s-md-gates) and are used as built. The cell holds 0.12 M NaCl when the salt is counted in its whole volume.
