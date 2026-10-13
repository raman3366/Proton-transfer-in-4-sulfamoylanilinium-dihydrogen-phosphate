| Ligand, scoring, input | Glu106 | Four waters | RMSD (Å) | Top pose ≤ 1.0 Å | Modes of the rotated pose | Score (kcal mol⁻¹) | Zn–N (Å) | Strict gate | 2.0 Å criterion |
|:----------------------|:----------|:-----|----------:|--------:|-------:|-------:|-----:|:-----|:-------|
| Sulfanilamide, AD4Zn, DFT geometry, 20 seeds | acid | no | 0.61 ± 0.04 | 20 of 20 | 40 | −8.836 | 2.20 | – | – |
|  | carboxylate | no | 1.34 ± 1.30 | 15 of 20 | 44 | −8.882 | 2.21 | – | – |
|  | acid | yes | 0.51 ± 0.04 | 20 of 20 | 0 | −9.232 | 2.16 | – | – |
|  | carboxylate | yes | 0.51 ± 0.03 | 20 of 20 | 0 | −9.241 | 2.16 | – | – |
|  | acid | yes^a^ | 0.51 ± 0.03 | 20 of 20 | 0 | −9.242 | 2.16 | – | – |
|  | acid | yes^b^ | 0.52 ± 0.04 | 20 of 20 | 0 | −9.538 | 2.12 | – | – |
| Sulfanilamide, AD4Zn, ETKDG, 10 seeds | acid | no | 2.07 ± 1.45 | 5 of 10 | 20 | −8.862 | 2.22 | fail | fail |
|  | carboxylate | no | 2.70 ± 1.38 | 3 of 10 | 20 | −8.911 | 2.23 | fail | fail |
|  | acid | yes | 0.58 ± 0.03 | 10 of 10 | 0 | −9.286 | 2.14 | pass | pass |
|  | carboxylate | yes | 0.60 ± 0.03 | 10 of 10 | 0 | −9.291 | 2.15 | pass | pass |
| Sulfonamidate, AD4Zn, ETKDG, 10 seeds | acid | no | 1.12 ± 0.32 | 3 of 10 | 20 | −8.133 | 2.21 | fail | pass |
|  | carboxylate | no | 0.86 ± 0.33 | 7 of 10 | 20 | −7.764 | 2.22 | pass | pass |
|  | acid | yes | 0.61 ± 0.04 | 10 of 10 | 0 | −8.131 | 2.23 | pass | pass |
|  | carboxylate | yes | 0.58 ± 0.04 | 10 of 10 | 0 | −7.764 | 2.22 | pass | pass |
| Sulfanilamide, Vina, ETKDG, 10 seeds | acid | no | 1.10 ± 0.01 | 0 of 10 | 14 | −5.870 | 2.84 | fail | fail |
|  | carboxylate | no | 1.09 ± 0.00 | 0 of 10 | 19 | −5.871 | 2.83 | fail | fail |
|  | acid | yes | 1.10 ± 0.01 | 0 of 10 | 0 | −6.279 | 2.93 | fail | fail |
|  | carboxylate | yes | 1.10 ± 0.01 | 0 of 10 | 0 | −6.302 | 2.90 | fail | fail |
| Sulfanilamide, Vinardo, ETKDG, 10 seeds | acid | no | 1.13 ± 0.00 | 0 of 10 | 23 | −5.075 | 2.48 | fail | fail |
|  | carboxylate | no | 1.12 ± 0.00 | 0 of 10 | 19 | −5.078 | 2.48 | fail | fail |
|  | acid | yes | 1.12 ± 0.00 | 0 of 10 | 0 | −5.498 | 2.50 | fail | fail |
|  | carboxylate | yes | 1.12 ± 0.00 | 0 of 10 | 0 | −5.498 | 2.50 | fail | fail |

Table: Docking of sulfanilamide into carbonic anhydrase II (6RL9) with Glu106 as acid or carboxylate, without and with four crystal waters. {#tbl:s-dock-glu106-water}

The acid is the state of the receptors behind every docking number of the main text and of Section S13 for which the carboxylate is not named; the waters are 423, 516, 553 and 568 of 6RL9, kept as fixed receptor atoms with hydrogens placed for the largest number of hydrogen bonds. RMSD is the symmetry-corrected heavy-atom RMSD of the top pose from the crystal ligand (mean ± SD over the seeds), followed by the number of seeds whose top pose lies within 1.0 Å of it. A mode of the rotated pose lies beyond 2.0 Å from the crystal ligand and within 2.0 Å of the ring-rotated pose of the dry receptor (seed 14, mode 3 of the first row); it is counted over every mode of every seed. The score is the mean score of the top poses, and Zn–N the mean distance from the zinc to their sulfonamide nitrogen. The gates of Section S13 apply to the validation runs (ETKDG conformers) only. ^a^ The hydrogen-bond ring of the four waters reversed. ^b^ The waters as bare oxygen atoms, as in the control receptor of Section S13. Exhaustiveness 32, 20 Å box, ligand inputs identical within each block of rows; AD4Zn, AutoDock4Zn. The SD of the score over the seeds is at most 0.021 kcal mol⁻¹ with AD4Zn, 0.061 with Vina and 0.013 with Vinardo. The rows with the waters are a diagnostic and not results of the protocol: a rigid receptor cannot let a ligand displace water, so the rotated pose is excluded by construction, the four waters were always kept together, and a pass of a gate in such a receptor does not validate the protocol. Stored modes of the rotated pose scored without a search in the maps with the waters give +5,346 to +7,644 kcal mol⁻¹ over every placement of the water hydrogens tried; in the receptors of the table the best near-native mode of each seed changes by −0.089 (acid) and −0.084 kcal mol⁻¹ (carboxylate) on average.
