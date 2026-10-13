# 05 Tables

Every table of the main paper (`main/`) and of the supplementary material (`SI/`), in two forms:

- `.md`: the source used to build the Word and PDF files (pandoc Markdown: caption after `Table:`, notes below);
- `.csv`: the same rows and columns as plain text, generated from the Markdown by the packaging script. The
  CSV files are UTF-8 with a byte-order mark, so Excel shows Å, ⁻¹ and Greek letters correctly. Formatting is dropped
  (*d* → d, pK~a~ → pKa, E~HOMO~ → E_HOMO), footnote marks become [a], [b], citations become the printed reference
  numbers, and a number printed with a minus sign is written with an ASCII minus so that it reads as a number.
  Captions and notes are only in the `.md` files.

Files are named by their printed number. The source name (as in the repository) is in the third column.

| Printed | File (.md, .csv) | Source file | Rows | Caption (first sentence) |
|---|---|---|---:|---|
| Table 1 | `main/Table_01_models` | `table01_models.md` | 11 | Models, their charges and sizes, and what each is used for. |
| Table 2 | `main/Table_02_proton` | `table02_proton.md` | 13 | Position of the proton in each environment. |
| Table 3 | `main/Table_03_hbonds` | `table04_hbonds.md` | 9 | Hydrogen bonds of the crystal. |
| Table 4 | `main/Table_04_geom` | `table03_geom.md` | 24 | Deviations of the optimised bond lengths (and of the heavy-atom bond angles, last block) from the X-ray structure, by bond class. |
| Table 5 | `main/Table_05_fmo` | `table05_fmo.md` | 7 | Frontier orbital energies and vertical ΔSCF reactivity descriptors (eV), B3LYP-D3(BJ)/6-311++G(d,p). |
| Table 6 | `main/Table_06_redock` | `table07_redock.md` | 11 | Redocking and cross-docking of the crystal ligands. |
| Table 7 | `main/Table_07_scores` | `table08_scores.md` | 23 | Production docking scores and ranks. |
| Table S1 | `SI/Table_S01_settings` | `tableS01_settings.md` | 18 | Settings of every calculation, with the case name and the folder of the data package (06_DFT_results/) that holds its results. |
| Table S2 | `SI/Table_S02_energies` | `tableS02_energies.md` | 16 | Energies of the optimised structures. |
| Table S3 | `SI/Table_S03_bonds` | `tableS03_bonds.md` | 26 | All bond lengths of the ion pair: X-ray and optimised (B3LYP-D3(BJ)/6-311++G(d,p) unless stated; periodic model PBE-D3(BJ), formula unit 1). |
| Table S4 | `SI/Table_S04_angles` | `tableS04_angles.md` | 22 | Deviations of the 22 heavy-atom bond angles from the X-ray structure (model − X-ray). |
| Table S5 | `SI/Table_S05_freq` | `tableS05_freq.md` | 75 | Harmonic frequencies and IR intensities of the ion pair in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p), finite-difference Hessian (step 0.005 Bohr), no imaginary frequency. |
| Table S6 | `SI/Table_S06_thermo` | `tableS06_thermo.md` | 14 | Harmonic frequencies and thermochemistry of the two minima of the ion pair in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p). |
| Table S7 | `SI/Table_S07_scan` | `tableS07_scan.md` | 8 | Relaxed proton scan along *d*(N1–H1A) in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p): chain A from the salt minimum, chain B from the neutral minimum, 0.10 Å steps. |
| Table S8 | `SI/Table_S08_eps` | `tableS08_eps.md` | 5 | Dielectric ladder: the ion pair optimised in IEF-PCM with *ε* = 1–78.4, B3LYP-D3(BJ)/6-311++G(d,p). |
| Table S9 | `SI/Table_S09_fcheck` | `tableS09_fcheck.md` | 5 | Functional check (kcal mol⁻¹). |
| Table S10 | `SI/Table_S10_zpe` | `tableS10_zpe.md` | 9 | Zero-point energy of the N1–H1A stretch against the proton-transfer barrier in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p). |
| Table S11 | `SI/Table_S11_sensitivity` | `tableS11_sensitivity.md` | 6 | The two contacts that close the R²₂(8) ring of the ion pair, and the agreement of every model with the X-ray geometry. |
| Table S12 | `SI/Table_S12_periodic_conv` | `tableS12_periodic_conv.md` | 7 | Convergence of the periodic PBE-D3(BJ) model at the X-ray geometry. |
| Table S13 | `SI/Table_S13_periodic_relax` | `tableS13_periodic_relax.md` | 4 | Periodic relaxations of the *Pc* cell and their energies relative to (i). |
| Table S14 | `SI/Table_S14_periodic_hbonds` | `tableS14_periodic_hbonds.md` | 9 | The nine hydrogen bonds of the crystal: deposited values and the relaxed periodic PBE-D3(BJ) model (i). |
| Table S15 | `SI/Table_S15_periodic_profile` | `tableS13b_periodic_profile.md` | 20 | Constrained proton-transfer profile in the crystal. |
| Table S16 | `SI/Table_S16_descriptors` | `tableS15_descriptors.md` | 10 | Global reactivity descriptors from vertical ΔSCF energies and from the orbital energies (Koopmans-type, *I* = −*E*~HOMO~, *A* = −*E*~LUMO~), B3LYP-D3(BJ)/6-311++G(d,p). |
| Table S17 | `SI/Table_S17_fukui` | `tableS_fukui.md` | 45 | Condensed Fukui functions of neutral sulfanilamide and the 4-sulfamoylanilinium cation in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p), each at its own IEF-PCM optimum. |
| Table S18 | `SI/Table_S18_charges` | `tableS16_charges.md` | 28 | Atomic charges of the ion pair in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p). |
| Table S19 | `SI/Table_S19_mep` | `tableS17_mep.md` | 14 | Extrema of the electrostatic potential on the 0.001 a.u. electron-density surface of the ion pair in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p). |
| Table S20 | `SI/Table_S20_qtaim` | `tableS18_qtaim.md` | 5 | QTAIM descriptors at the intermolecular bond critical points of the ion pair, B3LYP-D3(BJ)/6-311++G(d,p), with three hydrogen-bond energy estimates (Espinosa *V*/2; Emamian fits for neutral and charged complexes). |
| Table S21 | `SI/Table_S21_interaction` | `tableS19_interaction.md` | 3 | Interaction energies of the ion pair (kcal mol⁻¹), B3LYP-D3(BJ)/6-311++G(d,p), counterpoise-corrected; deformation energies relative to the relaxed fragments in the same environment (cation, then anion; sulfanilamide, then H₃PO₄). |
| Table S22 | `SI/Table_S22_cluster_qtaim` | `tableS18b_cluster_qtaim.md` | 9 | QTAIM descriptors of the nine hydrogen bonds of the crystal in the five-pair cluster. |
| Table S23 | `SI/Table_S23_hirshfeld` | `tableS20_hirshfeld.md` | 11 | Decomposition of the Hirshfeld surfaces of the cation and the anion into element-pair contacts. |
| Table S24 | `SI/Table_S24_hirshfeld_contacts` | `tableS21_hirshfeld_contacts.md` | 24 | Local minima of *d*~norm~ (red spots) on the Hirshfeld surfaces, with the interatomic distance of the contact (neutron-normalised H positions). |
| Table S25 | `SI/Table_S25_tddft_pair` | `tableS22_tddft_pair.md` | 20 | Lowest ten singlet excitations of the ion pair in IEF-PCM (water, non-equilibrium unless stated), TD-DFT/6-311++G(d,p) at the B3LYP-D3(BJ) geometry. |
| Table S26 | `SI/Table_S26_tddft_fragments` | `tableS23_tddft_fragments.md` | 30 | Lowest ten singlet excitations of the cation and sulfanilamide in IEF-PCM (water, non-equilibrium unless stated), TD-DFT/6-311++G(d,p) at the B3LYP-D3(BJ) geometry. |
| Table S27 | `SI/Table_S27_nlo` | `tableS24_nlo.md` | 4 | Finite-field NLO properties, B3LYP-D3(BJ)/6-311++G(d,p). |
| Table S28 | `SI/Table_S28_dock_validation` | `tableS25_dock_validation.md` | 38 | All validation and sensitivity docking runs: top-pose RMSD to the crystal ligand (mean, SD over seeds), success fractions, mean score and PoseBusters physics pass fraction. |
| Table S29 | `SI/Table_S29_dock_scores` | `tableS26_dock_scores.md` | 72 | Production docking scores of every species against every validated receptor. |
| Table S30 | `SI/Table_S30_torsion_scan` | `tableS26b_torsion_scan.md` | 29 | Relaxed torsional strain of neutral sulfanilamide and the sulfonamidate. |
| Table S31 | `SI/Table_S31_dock_torsion` | `tableS27_dock_torsion.md` | 30 | The sulfonamide torsion from the crystal to the docked poses, with its torsional strain. |
| Table S32 | `SI/Table_S32_dock_contacts` | `tableS28_dock_contacts.md` | 8 | Contacts of the representative top pose of the largest cluster for the sulfanilamide species (heavy-atom distances; PLIP/ProLIF/geometric). |
| Table S33 | `SI/Table_S33_control_receptors` | `tableS29_control_receptors.md` | 5 | Control receptor, the ATP site of *Escherichia coli* DNA gyrase B (4DUH, chain A). |
| Table S34 | `SI/Table_S34_control_scores` | `tableS30_control_scores.md` | 10 | Vina scores in the control receptor (4DUH) beside the scores of the same species in the four production structures. |
| Table S35 | `SI/Table_S35_admet` | `tableS_admet.md` | 27 | Drug-likeness descriptors and ADMET-AI predictions for neutral sulfanilamide, the species present at pH 7.4, with the 4-sulfamoylanilinium cation for reference, against experimental values for sulfanilamide. |
| Table S36 | `SI/Table_S36_dock_glu106_water` | `tableS31_dock_glu106_water.md` | 22 | Docking of sulfanilamide into carbonic anhydrase II (6RL9) with Glu106 as acid or carboxylate, without and with four crystal waters. |
| Table S37 | `SI/Table_S37_dock_glu106_other` | `tableS32_dock_glu106_other.md` | 10 | The other AutoDock4Zn runs of the carbonic anhydrase II receptors repeated with Glu106 as a carboxylate. |
| Table S38 | `SI/Table_S38_md_settings` | `tableS33_md_settings.md` | 43 | Settings of the molecular dynamics of the carbonic anhydrase II–sulfanilamide complex. |
| Table S39 | `SI/Table_S39_md_runs` | `tableS34_md_runs.md` | 11 | The eleven molecular dynamics runs of 100 ns. |
| Table S40 | `SI/Table_S40_md_gates` | `tableS35_md_gates.md` | 7 | Tests of the force field of the simulations, with those that failed. |

`SI/Section_S2_coordinates.md` holds the Cartesian coordinates of Section S2 of the supplementary material (not a
table; the same structures are the `opt.xyz` files in 06_DFT_results).

## Facts

- `facts.json` and `facts.md` list **every number printed in the paper and the supplementary material** (11111
  facts), each with its source file and key, unit, precision, caveats and any withdrawn or superseded value.
  `facts.md` is the readable version; the notes at its top (notes on the values, caveats, what was not computed)
  apply to the whole file.
- The literature values that enter computed numbers (pK~a~ values, half-life, IR scale factor, PLBD binding constants)
  are explained with their sources in 09_References/SOURCE_NOTES.md.
- `facts_rederivations/`: the independent re-derivations (`d01`–`d17`) that cross-check 3556 facts from the raw
  results (0 disagreements).
- Source paths in the facts are repository paths. `molecules/sfa/results/<case>/...` is `06_DFT_results/<case>/...`
  here (the scan, the dielectric ladder and the functional check are in their own folders; `periodic_qe/` is
  `periodic_QE/` and `hirshfeld_xray/` is `hirshfeld_surface/`), `dock/...` is in `07_Docking`, `md/...` is
  `11_Molecular_dynamics/md/...`, and
  `paper/audit/research/plbd_hca2_azm_sfm.json` is `10_Reproducibility/inputs/plbd_hca2_azm_sfm.json`.
  10_Reproducibility/source_map.csv maps every file of this package to its source.
- Rounding: round half to even on the full-precision value. A difference is computed from full precision and can
  differ by one unit in the last place from the difference of the two rounded end points (facts.md, "Rounding rule"
  at the top and section 30).
