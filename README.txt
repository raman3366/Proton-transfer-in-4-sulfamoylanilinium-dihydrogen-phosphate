SAPH_Supplementary_structures.zip: machine-readable supplementary material

Proton transfer in 4-sulfamoylanilinium dihydrogen phosphate: a combined molecular and periodic DFT study with molecular docking

This file belongs to the supplementary material of the article and supplements its PDF document. Coordinates are in
angstrom.

molecular/
  Optimised structures of Table S2 and Section S2 (XYZ; the comment line gives the level of theory and
  the electronic energy). File names are the case names of Table S1:
  molecular/sfa_pair_pcm_d3.xyz  (27 atoms): Ion pair (M1), salt minimum, IEF-PCM (water)
  molecular/sfa_pair_pcm_d3_oform.xyz  (27 atoms): Ion pair (M1), neutral (O5-protonated) minimum, IEF-PCM (water)
  molecular/sfa_pair_gas_d3.xyz  (27 atoms): Ion pair (M1) optimised in the gas phase: neutral adduct
  molecular/sfa_pair_pcm_b3lyp.xyz  (27 atoms): Ion pair (M1), IEF-PCM (water), B3LYP without D3(BJ)
  molecular/sfa_pair_pcm_d3_tz.xyz  (27 atoms): Ion pair (M1), IEF-PCM (water), def2-TZVP
  molecular/cation_pcm_d3.xyz  (20 atoms): Cation (M2), IEF-PCM (water)
  molecular/cation_gas_d3.xyz  (20 atoms): Cation (M2), gas phase
  molecular/h2po4_pcm_d3.xyz  (7 atoms): Dihydrogen phosphate (M3), IEF-PCM (water)
  molecular/h2po4_gas_d3.xyz  (7 atoms): Dihydrogen phosphate (M3), gas phase
  molecular/h3po4_pcm_d3.xyz  (8 atoms): Phosphoric acid (M6), IEF-PCM (water)
  molecular/h3po4_gas_d3.xyz  (8 atoms): Phosphoric acid (M6), gas phase
  molecular/sulfanilamide_pcm_d3.xyz  (19 atoms): Sulfanilamide (M4), IEF-PCM (water)
  molecular/sulfanilamide_gas_d3.xyz  (19 atoms): Sulfanilamide (M4), gas phase
  molecular/sulfanilamide_pcm_d3_tz.xyz  (19 atoms): Sulfanilamide (M4), IEF-PCM (water), def2-TZVP
  molecular/sulfonamidate_pcm_d3.xyz  (18 atoms): Sulfonamidate (M5), IEF-PCM (water)
  molecular/urea_gas_d3.xyz  (8 atoms): Urea (M7), gas phase

periodic/
  The periodic Pc cell (54 atoms, PBE-D3(BJ), cell fixed at the X-ray values) as CIF (P1) and as extended XYZ (lattice
  in the comment line): relax_i, the salt relaxed from the X-ray structure; relax_ii and relax_iir1, relaxed from
  starts with H1A on O5 in both formula units and in one only (Table S13); relax_iirfix, the
  constrained configuration with O5-H1A held at 0.98 A; relax_prof110 to relax_prof170, the points of the
  constrained N1-H1A profile at 1.10-1.70 A (Table S15).

clusters/
  The five-pair clusters (135 atoms) cut from the relaxed periodic salt and from the cell with H1A held on O5
  (Table S22).

docking/
  receptors/<target>/: receptor.pdbqt, box.json and vina_config.txt of the four production targets (3TZF, 6RL9, 3HS4,
  7S2K), and for the hCA II targets the AutoDock4Zn receptor and grid parameter file (ad4zn/).
  ligands/: the prepared ligand files (PDBQT, SDF) with ligands.json.
  top_poses/<target>/: the top pose of every seed of every production job (SDF; score, seed and cluster as SDF
  properties), in the receptor frame.
  control/: the control receptor of Section S13 (Escherichia coli DNA gyrase B, 4DUH): receptors/<receptor>/ with
  receptor.pdbqt, box.json and the crystal inhibitor (ref_RLI.sdf) of the receptor as prepared (C1_4DUH) and of its
  variants (second loop build, built residues removed, ordered crystal water kept), and top_poses/<receptor>/ as
  above (Tables S33 and S34).
  tests_hCA_II/<receptor>/: the receptors of the two docking tests of Section S15 (Tables S36 and S37).
  Each is a prepared receptor of the article with one change and is otherwise identical to it line for line:
  receptor.pdbqt, variant.json (what was changed, with checksums), box.json and the AutoDock4Zn receptor and grid
  parameter file (ad4zn/); comparison.json (Glu106 test) or comparison_water.json (water test) holds every compared
  number, and rescore_fixed_poses*.json the stored poses scored in the changed receptor without a search.
    T3_6RL9_glu106GLU: Glu106 test: the 6RL9 receptor with Glu106 as the carboxylate (the proton HE2 removed,
      charges of the other Glu residues); stored Vina, Vinardo and AutoDock4Zn runs of sulfanilamide, its anion and
      acetazolamide repeated in it (comparison.json lists them)
    T4_3HS4_glu106GLU: Glu106 test: the same change and the same kind of runs in the 3HS4 receptor
    T3_6RL9_wat4: water test: the 6RL9 receptor as reported (Glu106 the neutral acid) with crystal waters 423, 516,
      553 and 568 kept, their hydrogens placed by the rule hbond (the orientation with the most hydrogen bonds)
    T3_6RL9_glu106GLU_wat4: water test: the 6RL9 receptor with Glu106 as the carboxylate and the same four waters
      (rule hbond)
    T3_6RL9_wat4rev: water test, second placement of the water hydrogens (rule rev: the best network that keeps no
      water–water hydrogen bond of hbond); the key run only
    T3_6RL9_wat4bare: water test, the four water oxygens without hydrogens (rule bare); the key run only
    T3_6RL9_wat4fixed: water test, every water in one fixed orientation (rule fixed); stored poses were scored in
      it, no docking search was run

molecular_dynamics/
  The simulations of the carbonic anhydrase II-sulfanilamide complex of Section S16 (Tables S38, S39 and S40).
  inputs/: protocol.json (every setting the driver read); system/ (the prepared protein, the zinc, the crystal
  water, the three starts X, N and R as PDB and SDF, residue_map.csv, and the protonation audit in audit/); ligand/
  (the ligand parameters, their repaired terms and the checks against the DFT scans of the article); zinc_site/
  (the files of the zinc site that the tleap inputs load: the five residue templates with their RESP charges and
  the zinc parameter files of sets 2e and 2s, with CA2_tleap.in and CA2_mcpbpy.pdb as MCPB.py wrote them); topology/<set>/
  <start>/ (tleap input and log, the added terms, the build record and the checksums of each built system, for the
  three parameter sets of the reported runs). The built topology and coordinate files are not included; they are
  rebuilt with scripts/build_topology.py, and build.json and SHA256SUMS give their checksums.
  scripts/: the system and topology builders, the driver run_md.py, the two analysis programs (analysis/), the
  programs that collect the results (production/, with the records of the zinc-site gate), the queue program and
  the task files of the runs (queue/), and the software environment (env/).
  results/series/: one gzip-compressed CSV file for each of the 11 runs of 100 ns, one line per 10 ps from the
  release of the ligand to the end (time, d1, d2, ligand RMSD, state, configuration at the sulfonamide nitrogen,
  hydrogen bond to Thr199); the columns are explained in its README.md.
  results/final_frames/: the last solute frame of each of these runs (gzip-compressed PDB), with a README.md.
  results/runs/<set>/<run>/: the summary of each of the two analysis programs and the short series of the first
  (equilibration, all-atom frames, C-alpha fluctuations); results/campaign/: the tables over all runs and the
  report of the comparison of the two programs; conversion_times.json, run_health.json (the driver's records and
  checksums), run_events.json, site_checks_solute.json and site_checks_full.json.
  quantum_chemistry/: the zinc-site models and their charges (zinc_site/); gate_G-R/ and gate_G-R2/, the two tests of
  the force-field terms at the sulfonamide nitrogen against B3LYP-D3(BJ)/6-31G* (geometries, energies and outcome;
  both tests failed, as Section S16 reports), with the two relaxations of the ligand in the zinc-site model in
  gate_G-R/relax/ (its README.md says which file is which).
  The trajectories (about 13 GB) are not included; they are part of the data described under "Data
  availability" in the article. Paths inside these files are those of the project (md/..., /root/...).

inputs/pyscf/
  case_definitions/: charge, spin, level of theory, solvent, CIF atom labels and fragments of every molecular case;
  start_geometries/: the start structures cut from the CIF; scripts/: the PySCF driver scripts that run every
  molecular calculation (python dft/pipeline.py <case>; the scan, dielectric ladder, functional check, TD-DFT,
  Delta-SCF, NLO, interaction-energy and torsion-scan scripts are called by case).

inputs/quantum_espresso/
  The pw.x input files of the 45 reported periodic runs: relaxations (relax_*, with the _pre coarse
  pre-relaxations), constrained profile points (relax_prof*), convergence tests (scf_*), and PBE and PBE0 single
  points (pbe_*, hybrid_*). The pseudopotential files are not included (SSSP efficiency 1.3.0 and SG15 ONCV v1.2).

tddft/
  Tables S25 and S26 as CSV (UTF-8): excitation energies, wavelengths, oscillator strengths and
  charge-transfer numbers of the lowest singlet states.

Software: PySCF 2.14.0 with geomeTRIC 1.1.1 and D3(BJ) dispersion (molecular DFT); Quantum ESPRESSO 7.5 (periodic DFT);
AutoDock Vina 1.2.7 (Vina and Vinardo) and AutoDock4Zn on AutoGrid 4.2.7 maps, with Meeko 0.8.0 (docking);
OpenMM 8.6.1 with AmberTools 26.0 (molecular dynamics; the full list is in
molecular_dynamics/scripts/env/).

Labels and comments inside the script and JSON files keep the file names of the project, which are not the printed
numbers of the article: "Table S27" there is Table S31 of the supplementary material and "Fig. 11" is
Fig. 10 of the article. Absolute paths in these files (pseudo_dir and outdir of the Quantum ESPRESSO inputs, the "file" entries of
ligands.json) are those of the machines that ran the calculations; set them to your own.
