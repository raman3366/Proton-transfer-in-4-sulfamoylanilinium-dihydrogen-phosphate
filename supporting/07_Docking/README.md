# 07 Docking

Neutral sulfanilamide, the species SAPH releases at pH 7.4, is the primary ligand; the sulfonamidate (zinc site) and
the 4-sulfamoylanilinium cation are sensitivity cases. The salt or ion pair was never docked. The three sulfanilamide
species were docked from their DFT geometries (06_DFT_results). Reference ligands: sulfamethoxazole (neutral and
anion), acetazolamide (neutral and anion), *p*-aminobenzoate and 7,8-dihydropteroate.

**Scores are empirical docking scores** in kcal mol⁻¹ by the programs' convention, reported as mean ± SD over seeds
(10 seeds per production job; 20 for `T3_6RL9/sfa_neutral__ad4zn__prod`; 20 for `T4_3HS4/sfa_neutral__ad4zn__prod`). They are not binding affinities and were never converted to *K*~i~
or *K*~d~. No free energy is computed. The complex of hCA II with sulfanilamide was also simulated by molecular
dynamics, from the crystal pose and from two docked poses of this folder (11_Molecular_dynamics; Section
S16 of the supplementary material); the simulations replace no docking number.

## Reports and summaries

| File | What |
|---|---|
| `VALIDATION.md` | Gate G2: redocking and cross-docking of the crystal ligands, per target, with the pass criteria. |
| `RESULTS.md` | Production docking: scores, ranks, poses, contacts, sensitivity and diagnostics. |
| `summary/validation_summary.json`, `validation_gate_G2.json`, `validation_contacts.json` | The validation runs analysed (RMSD, PoseBusters, site checks) and the gate decision. |
| `summary/production_summary.json`, `production_contacts.json` | The production runs analysed (scores, clusters, RMSD to the crystal ligands, contacts). |
| `summary/targets_inventory.json` | Inventory of the PDB entries used (T1–T5; 3TYZ as the source of the *p*-aminobenzoate control ligand; the control receptor C1, 4DUH, and 1KZN as the source of the clorobiocin pose): chain, resolution, hetero groups, gaps, metals, ligand and box. |
| `control_gyrase_B/CONTROL.md`, `control_summary.json` | The control receptor (below): what was run, the redocking, the scores and the receptor variants. |
| `GLU106.md`, `WATER4.md` | The two tests on hCA II (below): Glu106 as the carboxylate, and four ordered crystal waters kept in 6RL9. |

## Targets

| Folder | PDB | Protein | Crystal ligand | Validation |
|---|---|---|---|---|
| `T1_3TZF` | 3TZF | DHPS, Y. pestis (chain B) | sulfamethoxazole (08D) | passes the 2.0 Å criterion only (SMX anion, Vina, 1.92 ± 0.04 Å); poses qualitative |
| `T1_3TZF_seed12345` | 3TZF | variant of `T1_3TZF` | | second loop rebuild (receptor-preparation sensitivity) |
| `T2_1AJ0` | 1AJ0 | DHPS (E. coli) | sulfanilamide | **failed validation; not used** |
| `T2_1AJ0_noSO4` | 1AJ0 | variant of `T2_1AJ0` | | variant without the sulfate ion |
| `T3_6RL9` | 6RL9 | hCA II | sulfanilamide (SAN) | partially validated with AD4Zn (2.07 ± 1.45 Å, f≤1 Å 0.5, bimodal: fails the 1.0 Å gate) |
| `T3_6RL9_glu106GLU` | 6RL9 | variant of `T3_6RL9` | | Glu106 test: the 6RL9 receptor with Glu106 as the carboxylate (the proton HE2 removed, charges of the other Glu residues); stored Vina, Vinardo and AutoDock4Zn runs of sulfanilamide, its anion and acetazolamide repeated in it (`comparison.json` lists them) |
| `T3_6RL9_glu106GLU_wat4` | 6RL9 | variant of `T3_6RL9` | | water test: the 6RL9 receptor with Glu106 as the carboxylate and the same four waters (rule `hbond`) |
| `T3_6RL9_his64HID` | 6RL9 | variant of `T3_6RL9` | | His64 as the HID tautomer (sensitivity) |
| `T3_6RL9_wat4` | 6RL9 | variant of `T3_6RL9` | | water test: the 6RL9 receptor as reported (Glu106 the neutral acid) with crystal waters 423, 516, 553 and 568 kept, their hydrogens placed by the rule `hbond` (the orientation with the most hydrogen bonds) |
| `T3_6RL9_wat4bare` | 6RL9 | variant of `T3_6RL9` | | water test, the four water oxygens without hydrogens (rule `bare`); the key run only |
| `T3_6RL9_wat4fixed` | 6RL9 | variant of `T3_6RL9` | | water test, every water in one fixed orientation (rule `fixed`); stored poses were scored in it, no docking search was run |
| `T3_6RL9_wat4rev` | 6RL9 | variant of `T3_6RL9` | | water test, second placement of the water hydrogens (rule `rev`: the best network that keeps no water–water hydrogen bond of `hbond`); the key run only |
| `T4_3HS4` | 3HS4 | hCA II | acetazolamide (AZM) | validated cross-dock for sulfanilamide with AD4Zn (0.62 ± 0.03 Å, 10/10); the AZM pose itself is not reproduced |
| `T4_3HS4_glu106GLU` | 3HS4 | variant of `T4_3HS4` | | Glu106 test: the same change and the same kind of runs in the 3HS4 receptor |
| `T4_3HS4_his64HID` | 3HS4 | variant of `T4_3HS4` | | His64 as the HID tautomer (sensitivity) |
| `T5_7S2K` | 7S2K | Sul2 (sulfonamide-resistant DHPS) | 7,8-dihydropteroate (78H) | validated with Vina and Vinardo (78H redock 0.51 / 0.54 Å, 10/10 seeds within 1 Å) |

`receptors/<target>/`: `receptor.pdb` (prepared protein, H added) and `receptor.pdbqt` (docking input), `box.json`
and `vina_config.txt` (search box), `ref_*.sdf` (crystal ligand in the receptor frame; `ref_cross_*` the ligand of the
other hCA II entry), `cofactors.pdbqt`, `propka.pka`, `prep_report*.json` (what the preparation did), the
intermediate PDB files, and for the zinc targets `ad4zn/` (AutoDock4Zn receptor, grid parameter and log files; the
maps themselves are not included, 22 MB, and are regenerated with `dock/run_ad4zn.py --stage maps` on an x86 machine).
`pdb_sources/` has the deposited mmCIF files.

## Poses

| Folder | What |
|---|---|
| `top_poses/<target>/` | For every target × species × scoring function: `*__top_pose_per_seed.sdf` (the top pose of each seed, with its score, seed and cluster as SDF properties) and `*__best_seed*_all_modes.pdbqt` (all modes of the best seed). AD4Zn sulfanilamide runs with other input geometries carry `csn102`, `csn109` or `etkdg` in the name (input-geometry diagnostic). |
| `redocking_overlays/Figure8*/` | The 4 overlays of Figure 8: `receptor.pdb`, `crystal_ligand.sdf` and the docked pose of the median-RMSD seed, in the receptor frame (no superposition). `overlay.pml` opens the set in PyMOL (`pymol overlay.pml` in the folder). |
| `redocking_overlays/crystal_references_in_receptor_frames/` | Crystal ligands carried into the other receptors' frames (used for the cross-dock RMSD). |
| `poses_for_dft/` | The docked conformers of the three sulfanilamide species passed to the rigid DFT strain single points of burst 2 (`__Hrot`: polar-H rotors re-optimised, because Vina does not score H positions), with `manifest.json`. Those single points are retired; the torsional strain now comes from the relaxed torsion scans (06_DFT_results/torsion_scan/). |

Figure 8 overlays:

| Folder | Crystal ligand | Job | Seed | RMSD to the crystal ligand (Å) |
|---|---|---|---:|---:|
| `Figure8a_6RL9_sulfanilamide_AD4Zn` | sulfanilamide (SAN, 6RL9) | `sfa_neutral__ad4zn__prod` (production) | 3 | 0.610 |
| `Figure8b_3HS4_sulfanilamide_crossdock_AD4Zn` | sulfanilamide of 6RL9 carried into 3HS4 by C-alpha superposition | `sfa_neutral__ad4zn` (validation) | 2 | 0.620 |
| `Figure8b_3HS4_sulfanilamide_crossdock_AD4Zn` | sulfanilamide of 6RL9 carried into 3HS4 by C-alpha superposition | `sfa_neutral__ad4zn__prod` (production) | 10 | 3.213 |
| `Figure8c_3TZF_sulfamethoxazole_anion_Vina` | sulfamethoxazole (08D, 3TZF chain B) | `smx_anion__vina__prod` (production) | 3 | 1.933 |
| `Figure8d_7S2K_dihydropteroate_Vina` | 7,8-dihydropteroate (78H, 7S2K) | `dhp_anion__vina__prod` (production) | 3 | 0.516 |

## Raw runs

`runs/validation_and_sensitivity/<target>/<run>/` and `runs/production/<target>/<run>/` hold every docking job:
`results.json` and, per seed, `ligand.pdbqt` (input), `out.pdbqt` (all poses) and `vina.log`. Run names are
`<species>__<scoring>[__variant]`; production runs end in `__prod`. `logs/` has the job logs, `diagnostics/` the
input-geometry diagnostic, `failed_receptor_prep/` one receptor preparation that failed (with `WHY.txt`).

## Two tests on hCA II: Glu106 and four crystal waters

Section S15 of the supplementary material. Both tests change one thing in a prepared receptor of the paper
and repeat stored runs with byte-identical ligand inputs. No printed docking number of the main text comes from
them: they say how far two results of the paper depend on the receptor as prepared.

- **Glu106** (`GLU106.md`). The two hCA II receptors of the paper carry Glu106 as the neutral acid, because PROPKA was
  run on the protein without the zinc (p*K*~a~ 8.91; with the zinc in its input, computed for 6RL9 only, 5.33).
  `receptors/T3_6RL9_glu106GLU/` and `T4_3HS4_glu106GLU/` hold the same receptors with Glu106 as the carboxylate;
  `variant.json` lists every changed line. With the carboxylate no Vina or Vinardo verdict changes; with AutoDock4Zn
  the DFT geometry of sulfanilamide gives the crystal pose in 15 of 20 seeds in 6RL9, against 20 with the
  acid. Which ring orientation ranks first therefore depends on this protonation state as well as on the input
  geometry. The test does not establish the protonation state of Glu106.
- **Four crystal waters** (`WATER4.md`). The ring-rotated pose of sulfanilamide lies on the sites of waters 423, 516,
  553 and 568 of the 6RL9 crystal. `receptors/T3_6RL9_wat4/` and `T3_6RL9_glu106GLU_wat4/` hold the receptor with
  these four waters for the two states of Glu106; `T3_6RL9_wat4rev/`, `T3_6RL9_wat4bare/` and `T3_6RL9_wat4fixed/`
  place the water hydrogens by other rules. In the 8 AutoDock4Zn runs with the waters no mode of the rotated pose
  is written. A rigid receptor with waters fixed in place forbids their displacement by construction: the test is a
  diagnostic of the two orientations, not a validation of the protocol, and it gives no cost of displacing the waters.

Files: `comparison.json` (Glu106) and `comparison_water.json` (waters) in a receptor folder hold every compared
number; `rescore_fixed_poses*.json` the stored poses scored without a search; `baseline_check*.json` the stored
runs of the unchanged receptor repeated bit for bit. The runs are under `runs/validation_and_sensitivity/<receptor>/`
and `runs/production/<receptor>/` like all others (`T3_6RL9_wat4fixed` has none). The scripts are
`10_Reproducibility/scripts/dock/glu_protonation_variant.py`, `water_variant.py`, `run_glu106_sensitivity.sh` and
`run_water4_sensitivity.sh`.

## Control receptor

`control_gyrase_B/` holds the docking into the ATP site of *Escherichia coli* DNA gyrase B (PDB 4DUH, chain A), an
enzyme that is not a target of sulfanilamide, run with the protocol of the five receptors (Section S13 of the supplementary
material). It is kept apart from them: none of its runs enters gate G2 or the production summary. The control is
read for the level of the scores only; it does not show that sulfanilamide does not bind gyrase B.

| Path | What |
|---|---|
| `CONTROL.md` | The report: redocking of the crystal inhibitor, the scores of every species, the comparison with the targets, the receptor variants and the whole-fragment search. |
| `control_summary.json` | Every number of the report and of the paper's control tables. |
| `receptors/C1_4DUH/` | The control receptor as prepared (residues 109–116 of the ATP lid built), with `ref_RLI.sdf` (the crystal inhibitor), `ref_cross_X_1KZN_CBN.sdf` (the clorobiocin pose of 1KZN in the 4DUH frame), `box.json` and `box_whole.json` (the whole-fragment box). |
| `receptors/C1_4DUH_seed12345/`, `C1_4DUH_noloop/`, `C1_4DUH_wat/` | The variants: residues 109–116 built from a second seed, removed, and the receptor with the ordered water of the crystal put back (diagnostic). |
| `runs/<receptor>/<run>/` | Every docking job, as in `runs/` above; `__whole` is the whole-fragment search. Clorobiocin is not docked: its pose in 1KZN does not fit this receptor (`CONTROL.md`, section 3.3). |
| `top_poses/<receptor>/` | The top pose of every seed (SDF) and all modes of the best seed (PDBQT) for every job. |
| `logs/` | The run log and the job logs. |

## Opening the files

PDB, PDBQT and SDF open in PyMOL, UCSF ChimeraX, Discovery Studio Visualizer or Maestro. Load `receptor.pdb` with a
pose file; all coordinates are in the receptor frame.
