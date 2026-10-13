# Control receptor: the ATP site of *E. coli* DNA gyrase B (4DUH)

Generated 2026-10-05 00:53 IST by `python dock/control_gyrb.py analyse` and `report`. Every number in this file is read from `dock/control/control_summary.json`.

## 1. Why a control receptor, and why this one

The five receptors of the study (T1–T5, `dock/VALIDATION.md` and `dock/RESULTS.md`) all belong to sulfonamide pharmacology: DHPS, hCA II and Sul2, the DHPS variant that rejects sulfonamides. A docking score says something about a ligand only if it changes when the receptor is not a target of the ligand, and nothing in gate G2 or in the production run measures that. The control is the ATP site of DNA gyrase B of *Escherichia coli*, in the 24 kDa N-terminal fragment crystallised with a 4,5′-bithiazole inhibitor (PDB 4DUH, chain A, 1.50 Å; Brvar et al. 2012, doi:10.1021/jm300395d). This site binds ATP and the aminocoumarins. The antibacterial target of sulfanilamide is dihydropteroate synthase, and the fluoroquinolones act on gyrase at the GyrA–DNA cleavage complex, not here; neither sulfanilamide nor ciprofloxacin is a known ligand of this pocket. No measurement says that sulfanilamide does not bind it, and a fragment of 11 heavy atoms may well bind weakly, so the control is described as a receptor that is not a target, never as one that the ligand does not bind.

The control is kept apart from the five receptors by construction: its directories are named `C1_*`, which the `T*` globs of `analyse.py`, `production_analysis.py` and `contacts.py` never match, and its output is in `dock/control/`. Gate G2, the production summary and every statistic of the five receptors are unchanged.

## 2. Receptors, ligands and protocol

Protocol of the five receptors: AutoDock Vina 1.2.7, Vina and Vinardo scoring, a 20 Å cube on the heavy-atom centroid of the crystal ligand, exhaustiveness 32, 10 seeds, `--cpu 3`; receptor from one chain without water, PROPKA protonation at pH 7.4, restrained minimisation of the built atoms and hydrogens. The three sulfanilamide species are docked from their B3LYP-D3(BJ)/IEF-PCM geometries, every other ligand from one ETKDG/MMFF94s conformer per seed. New ligands (`dock/prep_ligands.py`, SMILES checked against the RCSB chemical components): the 4DUH inhibitor as its benzoate (`rli_anion`, component RLI) and ciprofloxacin as its zwitterion (`cfx_zwit`, CPF). Clorobiocin is not docked (section 3.3).

| Receptor directory | What it is | Rebuilt residues to the crystal inhibitor (Å) | Rebuilt heavy atoms in the box |
|---|---|---:|---:|
| `C1_4DUH` | as prepared (residues 109–116 built, PDBFixer seed 20260928) | 8.27 | 10 |
| `C1_4DUH_seed12345` | residues 109–116 built from a second seed (12345) | 9.79 | 7 |
| `C1_4DUH_noloop` | residues 109–116 removed (observed residues only) | – | – |
| `C1_4DUH_wat` | as prepared, with the ordered water of chain A put back (diagnostic) | as `C1_4DUH` | as `C1_4DUH` |

Chain A of 4DUH lacks residues 109–116 of the ATP lid. PDBFixer builds them, as it builds the internal gaps of the five receptors; a loop built without density is a guess, so it was built twice. The two builds differ by 2.24 Å RMS in the built residues (largest 5.46 Å) and by 0.002 Å in every other heavy atom. The deposited atoms moved by at most 0.218 Å in the restrained minimisation. PROPKA leaves every Asp and Glu charged; non-default protonation: A:HIS83=HIP. `C1_4DUH_noloop` is `C1_4DUH` with the 74 atoms of the built residues deleted. `C1_4DUH_wat` is `C1_4DUH` with 322 of the 323 waters of chain A added as oxygen atoms (37 inside the box; one water that overlaps a built atom is left out). The water is fixed and has no hydrogens: each oxygen is an acceptor-type atom that fills space. It is a diagnostic, not the protocol of the five receptors, which were docked without such water; its scores are given beside those of the receptor as prepared and are marked as such.

6 of the 29 focused jobs ran with `--cpu 1` on the spare core beside the queue (`C1_4DUH_noloop: cfx_zwit__vina__prod`, `C1_4DUH_noloop: rli_anion__vina`, `C1_4DUH_noloop: sfa_neutral__vina__prod`, `C1_4DUH_wat: sfa_anion__vina__prod`, `C1_4DUH_wat: sfa_cation__vina__prod`, `C1_4DUH_wat: smx_anion__vina__prod`); the others ran with `--cpu 3`. The result does not depend on `--cpu`: the output file of seed 1 of `C1_4DUH/runs/sfa_neutral__vina__prod` was reproduced byte for byte at `--cpu 1` (sha256 identical), as `dock/VALIDATION.md` found for `--cpu 2`.

## 3. Redocking of the crystal inhibitor

Same analysis as gate G2 (`analyse.analyse_job`): symmetry-corrected heavy-atom RMSD of the top pose in the receptor frame, PoseBusters, and the three polar contacts of the deposited pose as site checks. The near-native mode is the first mode of a seed within 1.0 Å of the crystal pose.

| Receptor | Scoring | RMSD (Å) | f ≤ 1 Å | f ≤ 2 Å | Score | Near-native mode: seeds | rank | its RMSD (Å) | its score − top | Waters under the top pose | 1.0 Å gate | 2.0 Å |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| `C1_4DUH` | vina | 6.683 ± 0.015 | 0.0 | 0.0 | -9.269 ± 0.051 | 10/10 | 2–3 | 0.45–0.53 | 0.478 to 0.765 | 3–3 | fail | fail |
| `C1_4DUH` | vinardo | 6.771 ± 0.015 | 0.0 | 0.0 | -7.159 ± 0.047 | 10/10 | 2–2 | 0.31–0.57 | 0.105 to 0.366 | 3–3 | fail | fail |
| `C1_4DUH_wat` | vina | 0.518 ± 0.016 | 1.0 | 1.0 | -9.035 ± 0.080 | 10/10 | 1–1 | 0.50–0.55 | 0.000 to 0.000 | 0–0 | PASS | pass |
| `C1_4DUH_seed12345` | vina | 6.670 ± 0.016 | 0.0 | 0.0 | -9.179 ± 0.064 | 10/10 | 2–4 | 0.44–0.54 | 0.417 to 0.617 | 3–3 | fail | fail |
| `C1_4DUH_noloop` | vina | 6.681 ± 0.015 | 0.0 | 0.0 | -8.991 ± 0.043 | 10/10 | 4–5 | 0.40–0.53 | 0.231 to 0.468 | 3–3 | fail | fail |

**The receptor fails the gate.** With Vina the top pose of every seed lies 6.65–6.71 Å from the crystal pose, one cluster (10/10 seeds), and every top pose passes the PoseBusters physics checks (1.0). The search is not at fault: 10 of 10 seeds hold a mode within 0.45–0.53 Å of the crystal pose, at rank 2–3, 0.478–0.765 kcal/mol behind the top pose (mean 0.616).

| Contact of the deposited pose | Crystal (Å) | Top poses, mean (Å) | Pass fraction |
|---|---:|---:|---:|
| amide N to Asp73 N/O | 2.82 | 2.90 | 1.0 |
| carboxylate O to Arg136 N/O | 3.06 | 15.60 | 0.0 |
| aryl-NH-aryl N to Gly101 O | 2.82 | 3.85 | 0.0 |

The top pose keeps the first contact and gives up the other two: the amide N–H stays on the Asp73 carboxylate, the bithiazole stays in place, and the benzoate leaves Arg136 and turns under the ATP lid. Residues within 4 Å of the representative top pose: ALA100 (2.80), ASP73 (2.85), GLY117 (3.19), GLY101 (3.26), THR165 (3.37), VAL71 (3.39), PHE104 (3.40), LYS103 (3.41), ASN46 (3.42), GLY77 (3.56). The top pose overlaps 3–3 crystal waters of chain A (oxygen within 2.4 Å of a heavy atom of the pose); the near-native mode overlaps 0–0, and the crystal ligand itself 0 (its nearest water is 2.81 Å away). Two things can favour the turned pose: neither Vina nor Vinardo has an electrostatic term, so the salt bridge between the benzoate and the Arg136 guanidinium group counts only as a hydrogen bond, and the protocol removes the water that fills the space under the lid in the crystal. Only the second was tested (section 3.1), with Vina, and with all the ordered water of chain A at once; the three waters under the lid were not tested alone.

### 3.1 The ordered water put back (diagnostic)

With the ordered water of chain A in the receptor the top pose of the inhibitor lies 0.52 ± 0.02 Å from the crystal pose (1.0 of the seeds within 1 Å, 1.0 within 2 Å), against 6.68 Å without it; score -9.035 against -9.269 kcal/mol.

| Ligand | Score without water | Score with water | Difference | In the inhibitor site, without / with water |
|---|---:|---:|---:|---|
| 4DUH inhibitor (crystal ligand, RLI) | -9.269 | -9.035 ± 0.080 | +0.234 | 1.0 / 1.0 |
| sulfanilamide | -5.960 | -5.932 ± 0.006 | +0.028 | 0.0 / 1.0 |
| ciprofloxacin | -8.294 | -7.768 ± 0.256 | +0.526 | 1.0 / 1.0 |
| 4-sulfamoylanilinium | -5.945 | -5.932 ± 0.011 | +0.013 | 0.0 / 1.0 |
| sulfonamidate | -5.955 | -5.936 ± 0.008 | +0.019 | 0.0 / 1.0 |
| sulfamethoxazole anion | -7.408 | -7.234 ± 0.028 | +0.174 | 1.0 / 1.0 |
| sulfamethoxazole | -7.531 | -7.237 ± 0.032 | +0.294 | 1.0 / 1.0 |
| p-aminobenzoate | -5.541 | -5.512 ± 0.011 | +0.029 | 1.0 / 1.0 |
| acetazolamide | -6.271 | -6.198 ± 0.120 | +0.073 | 1.0 / 1.0 |
| acetazolamide anion | -6.288 | -6.143 ± 0.111 | +0.145 | 1.0 / 1.0 |

### 3.2 The built residues: a second build, and none

| Ligand | Primary | Second build | Difference | Residues removed | Difference | Top pose within 4 Å of a built residue (primary) | Representative pose moved, second build / no loop (Å) |
|---|---:|---:|---:|---:|---:|---:|---|
| 4DUH inhibitor (crystal ligand, RLI) | -9.269 ± 0.051 | -9.179 | +0.090 | -8.991 | +0.278 | 0.1 | 0.08 / 0.04 |
| sulfanilamide | -5.960 ± 0.011 | -5.852 | +0.108 | -5.832 | +0.128 | 1.0 | 0.05 / 8.57 |
| ciprofloxacin | -8.294 ± 0.636 | -8.279 | +0.015 | -8.253 | +0.041 | 0.0 | 0.03 / 0.03 |

Jobs of the primary receptor whose top pose lies within 4 Å of a built residue: `rli_anion__vina`, `rli_anion__vinardo`, `sfa_anion__vina__prod`, `sfa_cation__vina__prod`, `sfa_neutral__vina__prod`, `smx_anion__vina__prod`, `smx_neutral__vina__prod`; the closest approach of any top pose to a built heavy atom is 3.47 Å.

### 3.3 Clorobiocin: no cross-docking test

Clorobiocin is crystallised with the same fragment (1KZN, 2.3 Å; Lafitte et al. 2002, doi:10.1021/bi0159837). 1KZN is not numbered like *E. coli* GyrB between residues 82 and 96, so its pose is carried into the 4DUH frame by a sequence-aligned superposition of 183 Cα atoms (RMSD 0.781 Å; 0.522 Å over 164 atoms after trimming). In 4DUH that pose does not fit: its closest heavy-atom contact with the prepared receptor is 1.65 Å (LYS103:NZ), with 11 pairs below 2.5 Å (LYS103 7, ARG76 1, ASP73 1, GLU50 1, GLY101 1) and 39 below 3.0 Å, in residues ASN46, GLU50, ASP73, ARG76, GLY77, ILE94, GLY101, LYS103, ARG136, THR165. 8 of the 11 pairs are with GLY101, LYS103, residues that are ordered in 4DUH and have no counterpart in the 1KZN model (residues of 4DUH without one: 83, 84, 85, 86, 101, 102, 103, 104, 105, 106, 107, 108); the other pairs are with residues ordered in both crystals and show that the two structures also differ locally. A cross-docking test of 4DUH with clorobiocin could therefore not pass, and clorobiocin is not docked: a score of a pose that the receptor cannot hold would say nothing.

## 4. Scores of every species in the control receptor

A pose is *in the cavity of the inhibitor* when its heavy-atom centroid lies within 4.0 Å of a heavy atom of the crystal inhibitor. This says where a pose is, not that it is right: the top pose of the inhibitor itself in the receptor as prepared counts in every seed although it lies 6.7 Å from the crystal pose. Scores are empirical docking scores in kcal/mol, mean ± SD over the seeds, never affinities.

### Vina, `C1_4DUH` (protocol of the targets)

| Rank | Ligand | Heavy atoms | Score | Best seed | LE | Clusters (seeds in the first) | Seeds in the cavity of the inhibitor | Centroid to the inhibitor (Å) | Seeds with N/O within 3.5 Å of Asp73 | Waters under the pose | To a built residue (Å) | Nearest residues |
|---:|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---|
| 1 | 4DUH inhibitor (crystal ligand, RLI) | 26 | -9.269 ± 0.051 | -9.337 | 0.357 | 1 (10) | 10/10 | 1.48 | 10/10 | 3.0 | 4.0 | ALA100, ASP73, GLY117, GLY101 |
| 2 | ciprofloxacin | 24 | -8.294 ± 0.636 | -9.225 | 0.346 | 3 (4) | 10/10 | 1.07 | 7/10 | 2.2 | 5.97 | GLY77, THR165, GLU50, ARG76 |
| 3 | sulfamethoxazole | 17 | -7.531 ± 0.013 | -7.548 | 0.443 | 1 (10) | 10/10 | 3.33 | 0/10 | 4.0 | 3.54 | ASN46, ALA100, GLY119, GLY117 |
| 4 | sulfamethoxazole anion | 17 | -7.408 ± 0.043 | -7.466 | 0.436 | 1 (10) | 10/10 | 3.21 | 0/10 | 5.0 | 3.73 | VAL120, ASN46, LYS103, ALA100 |
| 5 | acetazolamide anion | 13 | -6.288 ± 0.015 | -6.316 | 0.484 | 2 (2) | 10/10 | 1.28 | 8/10 | 1.8 | 5.37 | ASN46, HIS99, ALA100, PHE104 |
| 6 | acetazolamide | 13 | -6.271 ± 0.021 | -6.303 | 0.482 | 2 (7) | 10/10 | 1.43 | 7/10 | 2.2 | 5.45 | ASN46, VAL43, LYS103, ALA47 |
| 7 | sulfanilamide | 11 | -5.960 ± 0.011 | -5.978 | 0.542 | 1 (10) | 0/10 | 5.02 | 0/10 | 4.0 | 3.5 | ASN46, VAL120, ALA100, GLY117 |
| 8 | sulfonamidate | 11 | -5.955 ± 0.020 | -6.000 | 0.541 | 1 (10) | 0/10 | 5.04 | 0/10 | 4.0 | 3.48 | ASN46, VAL120, ALA100, GLY117 |
| 9 | 4-sulfamoylanilinium | 11 | -5.945 ± 0.009 | -5.959 | 0.540 | 1 (10) | 0/10 | 5.03 | 0/10 | 4.0 | 3.47 | ASN46, VAL120, ALA100, GLY117 |
| 10 | p-aminobenzoate | 10 | -5.541 ± 0.003 | -5.545 | 0.554 | 1 (10) | 10/10 | 0.62 | 10/10 | 2.0 | 6.58 | ASN46, THR165, ASP73, LYS103 |

Pearson r between the mean score and the number of heavy atoms over these 10 species (10–26 heavy atoms): -0.986; least-squares line -0.2126 kcal/mol per heavy atom.

### Vinardo, `C1_4DUH`

| Rank | Ligand | Heavy atoms | Score | Best seed | LE | Clusters (seeds in the first) | Seeds in the cavity of the inhibitor | Centroid to the inhibitor (Å) | Seeds with N/O within 3.5 Å of Asp73 | Waters under the pose | To a built residue (Å) | Nearest residues |
|---:|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---|
| 1 | 4DUH inhibitor (crystal ligand, RLI) | 26 | -7.159 ± 0.047 | -7.216 | 0.275 | 1 (10) | 10/10 | 1.55 | 10/10 | 3.0 | 3.83 | ALA100, ASP73, GLY117, GLY101 |
| 2 | ciprofloxacin | 24 | -5.805 ± 0.414 | -6.287 | 0.242 | 1 (10) | 10/10 | 1.3 | 9/10 | 0.9 | 7.15 | THR165, GLY101, GLY77, ARG76 |
| 3 | sulfanilamide | 11 | -4.974 ± 0.025 | -5.004 | 0.452 | 2 (6) | 4/10 | 3.54 | 4/10 | 3.2 | 4.5 | SER121, ALA100, HIS99, VAL120 |

Pearson r between the mean score and the number of heavy atoms over these 3 species (11–26 heavy atoms): -0.859; least-squares line -0.1163 kcal/mol per heavy atom.

### Vina, `C1_4DUH_wat` (ordered crystal water kept; the receptor that redocks the inhibitor)

| Rank | Ligand | Heavy atoms | Score | Best seed | LE | Clusters (seeds in the first) | Seeds in the cavity of the inhibitor | Centroid to the inhibitor (Å) | Seeds with N/O within 3.5 Å of Asp73 | Waters under the pose | To a built residue (Å) | Nearest residues |
|---:|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---|
| 1 | 4DUH inhibitor (crystal ligand, RLI) | 26 | -9.035 ± 0.080 | -9.129 | 0.348 | 1 (10) | 10/10 | 1.3 | 10/10 | 0.0 | 8.38 | GLY101, ARG76, HOH614, ARG136 |
| 2 | ciprofloxacin | 24 | -7.768 ± 0.256 | -8.167 | 0.324 | 1 (10) | 10/10 | 1.26 | 2/10 | 0.0 | 7.1 | HOH614, GLY101, ARG76, GLY77 |
| 3 | sulfamethoxazole | 17 | -7.237 ± 0.032 | -7.278 | 0.426 | 1 (10) | 10/10 | 1.36 | 9/10 | 0.0 | 7.77 | ASP73, ARG76, LYS103, ARG136 |
| 4 | sulfamethoxazole anion | 17 | -7.234 ± 0.028 | -7.274 | 0.426 | 1 (10) | 10/10 | 1.27 | 0/10 | 0.0 | 7.71 | ARG76, GLU50, ILE78, LYS103 |
| 5 | acetazolamide | 13 | -6.198 ± 0.120 | -6.289 | 0.477 | 1 (10) | 10/10 | 0.69 | 10/10 | 0.0 | 7.09 | HOH445, ASP73, VAL167, LYS103 |
| 6 | acetazolamide anion | 13 | -6.143 ± 0.111 | -6.295 | 0.473 | 1 (10) | 10/10 | 0.76 | 10/10 | 0.0 | 7.15 | HOH445, ASP73, VAL167, LYS103 |
| 7 | sulfonamidate | 11 | -5.936 ± 0.008 | -5.943 | 0.540 | 1 (10) | 10/10 | 0.34 | 10/10 | 0.0 | 7.13 | HOH432, ASP73, ASN46, LYS103 |
| 8 | sulfanilamide | 11 | -5.932 ± 0.006 | -5.944 | 0.539 | 1 (10) | 10/10 | 0.32 | 10/10 | 0.0 | 7.13 | ASN46, ASP73, LYS103, HOH445 |
| 9 | 4-sulfamoylanilinium | 11 | -5.932 ± 0.011 | -5.946 | 0.539 | 1 (10) | 10/10 | 0.33 | 10/10 | 0.0 | 7.11 | ASN46, ASP73, LYS103, HOH445 |
| 10 | p-aminobenzoate | 10 | -5.512 ± 0.011 | -5.527 | 0.551 | 1 (10) | 10/10 | 0.19 | 10/10 | 0.0 | 7.35 | ASP73, ASN46, HOH445, HOH432 |

Pearson r between the mean score and the number of heavy atoms over these 10 species (10–26 heavy atoms): -0.973; least-squares line -0.1886 kcal/mol per heavy atom.

## 5. The same species in the four production structures

Production means of `dock/results/production_summary.json` (read only): 3TZF (DHPS), 7S2K (Sul2, which rejects sulfonamides), 6RL9 and 3HS4 (hCA II; Vina fails the redocking gate there, having no metal term). Difference = control receptor minus production structure; positive = weaker score in the control receptor.

### Vina, `C1_4DUH` (protocol of the targets)

| Species | 4DUH | 3TZF | 7S2K | 6RL9 | 3HS4 | 4DUH − 3TZF | 4DUH − 7S2K |
|---|---:|---:|---:|---:|---:|---:|---:|
| sulfanilamide | -5.960 | -6.423 | -6.409 | -5.880 | -5.767 | +0.463 | +0.449 |
| 4-sulfamoylanilinium | -5.945 | -6.403 | -6.376 | -5.873 | -5.763 | +0.458 | +0.431 |
| sulfonamidate | -5.955 | -6.400 | -6.678 | -5.838 | -5.758 | +0.445 | +0.723 |
| sulfamethoxazole anion | -7.408 | -7.636 | -7.853 | -7.728 | -7.721 | +0.228 | +0.445 |
| sulfamethoxazole | -7.531 | -7.799 | -8.170 | -7.529 | -7.432 | +0.268 | +0.639 |
| p-aminobenzoate | -5.541 | -6.348 | -7.501 | -5.449 | -5.418 | +0.807 | +1.960 |
| acetazolamide | -6.271 | -5.582 | -6.971 | -6.653 | -6.409 | -0.689 | +0.700 |
| acetazolamide anion | -6.288 | -5.552 | -6.826 | -6.637 | -6.477 | -0.736 | +0.538 |

- against 3TZF (8 species): mean difference +0.155 kcal/mol (-0.736 to +0.807); Pearson r 0.743, Spearman ρ 0.381 (p = 0.3518).
- against 7S2K (8 species): mean difference +0.736 kcal/mol (+0.431 to +1.960); Pearson r 0.736, Spearman ρ 0.595 (p = 0.1195).
- against 6RL9 (8 species): mean difference +0.086 kcal/mol (-0.117 to +0.382); Pearson r 0.97, Spearman ρ 0.929 (p = 0.0009).
- against 3HS4 (8 species): mean difference -0.019 kcal/mol (-0.197 to +0.313); Pearson r 0.979, Spearman ρ 0.952 (p = 0.0003).

Size relation over the 8 species docked into every receptor (Pearson r between mean score and heavy atoms): 4DUH -0.992, 3TZF -0.678, 7S2K -0.757, 6RL9 -0.99, 3HS4 -0.994.

### Vinardo, `C1_4DUH`

| Species | 4DUH | 3TZF | 7S2K | 6RL9 | 3HS4 | 4DUH − 3TZF | 4DUH − 7S2K |
|---|---:|---:|---:|---:|---:|---:|---:|
| sulfanilamide | -4.974 | -5.310 | -4.991 | -5.077 | -4.987 | +0.336 | +0.017 |


### Vina, `C1_4DUH_wat` (ordered crystal water kept; the receptor that redocks the inhibitor)

| Species | 4DUH | 3TZF | 7S2K | 6RL9 | 3HS4 | 4DUH − 3TZF | 4DUH − 7S2K |
|---|---:|---:|---:|---:|---:|---:|---:|
| sulfanilamide | -5.932 | -6.423 | -6.409 | -5.880 | -5.767 | +0.491 | +0.477 |
| 4-sulfamoylanilinium | -5.932 | -6.403 | -6.376 | -5.873 | -5.763 | +0.471 | +0.444 |
| sulfonamidate | -5.936 | -6.400 | -6.678 | -5.838 | -5.758 | +0.464 | +0.742 |
| sulfamethoxazole anion | -7.234 | -7.636 | -7.853 | -7.728 | -7.721 | +0.402 | +0.619 |
| sulfamethoxazole | -7.237 | -7.799 | -8.170 | -7.529 | -7.432 | +0.562 | +0.933 |
| p-aminobenzoate | -5.512 | -6.348 | -7.501 | -5.449 | -5.418 | +0.836 | +1.989 |
| acetazolamide | -6.198 | -5.582 | -6.971 | -6.653 | -6.409 | -0.616 | +0.773 |
| acetazolamide anion | -6.143 | -5.552 | -6.826 | -6.637 | -6.477 | -0.591 | +0.683 |

- against 3TZF (8 species): mean difference +0.252 kcal/mol (-0.616 to +0.836); Pearson r 0.752, Spearman ρ 0.347 (p = 0.3993).
- against 7S2K (8 species): mean difference +0.832 kcal/mol (+0.444 to +1.989); Pearson r 0.714, Spearman ρ 0.635 (p = 0.0909).
- against 6RL9 (8 species): mean difference +0.183 kcal/mol (-0.098 to +0.494); Pearson r 0.968, Spearman ρ 0.898 (p = 0.0024).
- against 3HS4 (8 species): mean difference +0.078 kcal/mol (-0.178 to +0.487); Pearson r 0.977, Spearman ρ 0.874 (p = 0.0045).

Size relation over the 8 species docked into every receptor (Pearson r between mean score and heavy atoms): 4DUH -0.987, 3TZF -0.678, 7S2K -0.757, 6RL9 -0.99, 3HS4 -0.994.

Neutral sulfanilamide with Vina, every variant of the control receptor: `C1_4DUH` -5.960, `C1_4DUH_seed12345` -5.852, `C1_4DUH_noloop` -5.832, `C1_4DUH_wat` -5.932; production structures: 3TZF -6.423, 7S2K -6.409, 6RL9 -5.880, 3HS4 -5.767. Largest difference between any variant and any production structure: 0.591 kcal/mol.

## 6. Whole-fragment search

Box 53.5 × 56.0 × 57.0 Å around the whole fragment (21.35 times the volume of the 20 Å cube), exhaustiveness 683 (32 x volume ratio, rounded: the same number of Monte Carlo runs per unit volume as the focused search).

| Job | Seeds | Score | Top pose in the inhibitor site | Seeds with any mode in the site | Top pose to the reference (Å) |
|---|---:|---:|---:|---:|---|
| `C1_4DUH/runs/sfa_neutral__vina__whole` | 5 | -5.959 ± 0.010 | 0/5 | 5/5 | 0.08 ± 0.02 (≤ 2 Å in 5/5; representative pose of the focused search (sfa_neutral__vina__prod, seed 1)) |

Top poses outside the inhibitor site (`C1_4DUH/runs/sfa_neutral__vina__whole`), by seed, with their nearest residues: seed 1: ASN46, VAL120, ALA100, GLY117, LYS103, PHE104; seed 2: ASN46, VAL120, ALA100, GLY117, LYS103, PHE104; seed 3: ASN46, VAL120, ALA100, GLY117, LYS103, PHE104; seed 4: ASN46, VAL120, ALA100, GLY117, LYS103, PHE104; seed 5: ASN46, VAL120, ALA100, GLY117, LYS103, PHE104.

## 7. What the control supports

- The level of the scores, not the poses. The receptor as prepared does not pass the redocking gate, so no pose in it is read as a binding mode; the receptor that passes it holds the crystal water and is not the protocol of the five receptors.
- A pocket that is no target of sulfanilamide can be compared with the production structures only through the same scoring function, Vina (section 5); the control was not run with AutoDock4Zn, which is specific to the zinc site.
- It does not show that sulfanilamide does not bind gyrase B. It shows that the Vina score does not tell this pocket from DHPS, Sul2 or hCA II.
- The diagnostic with the crystal water, the two loop builds and the loop-free receptor say how much of each number depends on choices of the receptor model (section 3).

## 8. Reproduce

```
python dock/targets.py
python dock/prep_receptor.py C1_4DUH
python dock/prep_receptor.py C1_4DUH --seed 12345 --suffix _seed12345
python dock/prep_ligands.py --only rli_anion cfx_zwit
python dock/control_gyrb.py run        # builds C1_4DUH_wat and C1_4DUH_noloop, then every job (resumes)
python dock/control_gyrb.py analyse    # dock/control/control_summary.json, poses
python dock/control_gyrb.py report     # this file
```
