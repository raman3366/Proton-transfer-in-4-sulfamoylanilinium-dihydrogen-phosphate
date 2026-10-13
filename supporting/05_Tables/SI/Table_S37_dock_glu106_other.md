| Receptor | Ligand, input | Glu106 | RMSD (Å) | Top pose ≤ 1.0 Å | Score (kcal mol⁻¹) | Zn–N (Å) | Seeds with Zn–N ≤ 2.3 Å | Strict gate | 2.0 Å criterion |
|:-------|:----------------------|:----------|----------:|--------:|-------:|-----:|--------:|:-----|:-------|
| 3HS4 | Sulfanilamide, ETKDG, cross-dock, 10 seeds | acid | 0.62 ± 0.03 | 10 of 10 | −8.789 | 2.18 | 10 of 10 | pass | pass |
|  |  | carboxylate | 0.84 ± 0.72 | 9 of 10 | −8.842 | 2.13 | 10 of 10 | pass | pass |
| 3HS4 | Sulfanilamide, DFT geometry, cross-dock, 20 seeds | acid | 3.22 ± 0.02 | 0 of 20 | −8.817 | 2.32 | 6 of 20 | – | – |
|  |  | carboxylate | 3.09 ± 0.60 | 1 of 20 | −8.883 | 2.31 | 9 of 20 | – | – |
| 6RL9 | Acetazolamide, ETKDG, cross-dock, 10 seeds | acid | 2.68 ± 0.81 | 0 of 10 | −9.484 | 3.68 | 5 of 10 | fail | fail |
|  |  | carboxylate | 2.20 ± 0.71 | 0 of 10 | −9.520 | 2.43 | 9 of 10 | fail | fail |
| 3HS4 | Acetazolamide, ETKDG, redock, 10 seeds | acid | 2.70 ± 0.20 | 0 of 10 | −9.258 | 2.08 | 10 of 10 | fail | fail |
|  |  | carboxylate | 2.51 ± 0.39 | 0 of 10 | −9.306 | 2.10 | 10 of 10 | fail | fail |
| 3HS4 | Acetazolamide anion, ETKDG, redock, 10 seeds | acid | 3.65 ± 0.44 | 0 of 10 | −9.006 | 4.86 | 0 of 10 | fail | fail |
|  |  | carboxylate | 3.34 ± 0.04 | 0 of 10 | −8.832 | 4.81 | 0 of 10 | fail | fail |

Table: The other AutoDock4Zn runs of the carbonic anhydrase II receptors repeated with Glu106 as a carboxylate. {#tbl:s-dock-glu106-other}

Columns as in @tbl:s-dock-glu106-water. For sulfanilamide every top pose that is not within 1.0 Å of the reference lies beyond 2.0 Å, in the ring-rotated orientation. The reference of a cross-docking run is the crystal ligand of the other structure carried over by Cα superposition. The acid rows are the runs of the main text and of Section S13; ligand inputs are identical within each pair of rows. The SD of the score over the seeds is at most 0.075 kcal mol⁻¹. Stored modes of the four sulfanilamide runs scored without a search in both sets of maps change by −0.043 to −0.053 kcal mol⁻¹ (near-native modes) and −0.061 to −0.066 kcal mol⁻¹ (modes beyond 2.0 Å) with the carboxylate, so the ring-rotated orientation gains 0.013 to 0.018 kcal mol⁻¹, the size of the separation of the two orientations. The Vina and Vinardo scores of 80 stored top poses do not change (0.000 kcal mol⁻¹), and no verdict of the 21 repeated Vina and Vinardo runs changes. Not repeated with the carboxylate: the six input-geometry runs of sulfanilamide, the 20-seed and higher-exhaustiveness variants of the 6RL9 validation run, every run of sulfamethoxazole, *p*-aminobenzoate and 4-sulfamoylanilinium, the other Vina validation variants, and the His64 variants. The heavy atoms of the receptor were not moved, and the test does not establish the protonation state of Glu106.
