| Receptor | Ligand | Scoring function | Ligand geometry | RMSD (Å) | Success fraction ≤ 1 Å (≤ 2 Å) | Strict gate / 2.0 Å criterion |
|:----------------------|:------------------------------|:----------------|:---------------|------------:|------------------------------:|:-----------------------------|
| 7S2K (Sul2) | Dihydropteroate (redock) | Vina | ETKDG | 0.51 ± 0.03 | 1.0 (1.0) | pass / pass |
| 7S2K (Sul2) | Dihydropteroate (redock) | Vinardo | ETKDG | 0.54 ± 0.01 | 1.0 (1.0) | pass / pass |
| 3HS4 (hCA II) | Sulfanilamide (cross-dock) | AD4Zn | ETKDG | 0.62 ± 0.03 | 1.0 (1.0) | pass / pass |
| 6RL9 (hCA II) | Sulfanilamide (redock) | AD4Zn | ETKDG | 2.07 ± 1.45 | 0.5 (0.5) | fail / fail |
| 6RL9 (hCA II) | Sulfonamidate (redock) | AD4Zn | ETKDG | 1.12 ± 0.32 | 0.3 (1.0) | fail / pass |
| 6RL9 (hCA II) | Sulfanilamide (redock) | Vina | ETKDG | 1.10 ± 0.01 | 0.0 (1.0) | fail / fail |
| 3TZF (DHPS, chain B) | Sulfamethoxazole anion (redock) | Vina | ETKDG | 1.92 ± 0.04 | 0.0 (1.0) | fail / pass |
| 3HS4 (hCA II) | Acetazolamide (redock) | AD4Zn | ETKDG | 2.70 ± 0.20 | 0.0 (0.0) | fail / fail |
| 1AJ0 (DHPS, *E. coli*) | Sulfanilamide (redock) | Vina | ETKDG | 11.64 ± 0.05 | 0.0 (0.0) | fail / fail |
| 6RL9 (hCA II) | Sulfanilamide (redock) | AD4Zn | DFT (20 seeds) | 0.61 ± 0.04 | 1.0 (1.0) | – |
| 3HS4 (hCA II) | Sulfanilamide (cross-dock) | AD4Zn | DFT (20 seeds) | 3.22 ± 0.02 | 0.0 (0.0) | – |

Table: Redocking and cross-docking of the crystal ligands. {#tbl:redock}

RMSD is the symmetry-corrected heavy-atom RMSD of the top-ranked pose in the receptor frame (mean ± SD over seeds, 10 seeds unless stated); the success fraction is the fraction of seeds whose top pose lies within 1 Å (2 Å). Every run passed the PoseBusters physics checks (pass fraction 1.0). Exhaustiveness 32, 20 Å box; AD4Zn, AutoDock4Zn. Validation used ETKDG/MMFF94s conformers; the last two rows repeat the hCA II sulfanilamide runs with the DFT geometry used in production. The hCA II rows are for the receptors as prepared, with Glu106 as the neutral acid (Section 2.5); with Glu106 as a carboxylate the 6RL9 run from the DFT geometry gives 1.34 ± 1.30 Å with 15 of 20 seeds within 1 Å, and the sulfonamidate redock then passes the strict gate (Table S36). The 6RL9 AD4Zn result is bimodal in validation (near-native pose found in every seed but ranked first in half); in that run the two orientations differ by 0.008–0.039 kcal mol⁻¹ (over eight AD4Zn runs the near-native pose trails the top pose by 0.000–0.219; near-native: within 1.0 Å of the crystal pose), and which one ranks first depends on the input geometry and on the state of Glu106. The Vina run in 6RL9 misses the strict gate by its RMSD and fails the 2.0 Å criterion on the zinc contact (Zn–N 2.84 Å against the 2.3 Å of the site check), not on the distance from the crystal pose. The 3HS4 cross-dock reference is 6RL9 sulfanilamide carried over by Cα superposition (0.328 Å).
