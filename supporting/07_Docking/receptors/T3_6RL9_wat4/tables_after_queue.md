
#### The waters and their hydrogen bonds, by hydrogen rule

| Variant (rule) | Water | B (Å²) | Donates to: acceptor, H···A (Å), O–H···A (°) | Accepts from: donor, H···O (Å), D–H···O (°) | Hydrogens without a partner | Hydrogen pairs closer than 2.0 Å |
|:--|:--|--:|:--|:--|--:|:--|
| `T3_6RL9_wat4` (hydrogen-bond network) | 423 | 12.7 | HOH568 O, 1.79, 176.8; ASN67 OD1, 1.66, 176.7 | ASN62 ND2-H, 2.05, 155.2; HOH516 H1, 1.71, 180.0 | 0 | – |
|  | 516 | 17.6 | HOH423 O, 1.71, 180.0 | HOH553 H1, 1.70, 173.9 | 1 | – |
|  | 553 | 11.5 | HOH516 O, 1.70, 173.9; THR200 OG1, 2.06, 174.3 | HOH568 H1, 1.93, 180.0 | 0 | – |
|  | 568 | 12.0 | HOH553 O, 1.93, 180.0 | GLN92 NE2-H, 2.15, 144.7; HOH423 H1, 1.79, 176.8 | 1 | – |
| `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 423 | 12.7 | HOH568 O, 1.79, 176.8; ASN67 OD1, 1.66, 176.7 | ASN62 ND2-H, 2.05, 155.2; HOH516 H1, 1.71, 180.0 | 0 | – |
|  | 516 | 17.6 | HOH423 O, 1.71, 180.0 | HOH553 H1, 1.70, 173.9 | 1 | – |
|  | 553 | 11.5 | HOH516 O, 1.70, 173.9; THR200 OG1, 2.06, 174.3 | HOH568 H1, 1.93, 180.0 | 0 | – |
|  | 568 | 12.0 | HOH553 O, 1.93, 180.0 | GLN92 NE2-H, 2.15, 144.7; HOH423 H1, 1.79, 176.8 | 1 | – |
| `T3_6RL9_wat4rev` (network, ring reversed) | 423 | 12.7 | HOH516 O, 1.79, 149.8; ASN67 OD1, 1.75, 149.5 | ASN62 ND2-H, 2.05, 155.2; HOH568 H1, 1.79, 180.0 | 0 | – |
|  | 516 | 17.6 | HOH553 O, 1.69, 180.0 | HOH423 H1, 1.79, 149.8 | 1 | – |
|  | 553 | 11.5 | HOH568 O, 1.93, 174.7; THR200 OG1, 2.06, 174.9 | HOH516 H1, 1.69, 180.0 | 0 | – |
|  | 568 | 12.0 | HOH423 O, 1.79, 180.0 | GLN92 NE2-H, 2.15, 144.7; HOH553 H1, 1.93, 174.7 | 1 | H2···GLN92 HE22 1.90 |
| `T3_6RL9_wat4bare` (bare oxygen (precedent)) | 423 | 12.7 | no hydrogens | ASN62 ND2-H, 2.05, 155.2 | – | – |
|  | 516 | 17.6 | no hydrogens | – | – | – |
|  | 553 | 11.5 | no hydrogens | – | – | – |
|  | 568 | 12.0 | no hydrogens | GLN92 NE2-H, 2.15, 144.7 | – | – |
| `T3_6RL9_wat4fixed` (one fixed orientation) | 423 | 12.7 | – | ASN62 ND2-H, 2.05, 155.2 | 2 | H1···ASN62 HD22 1.84; H2···ASN67 HD21 1.96 |
|  | 516 | 17.6 | – | – | 2 | – |
|  | 553 | 11.5 | THR200 OG1, 2.21, 141.9 | – | 1 | – |
|  | 568 | 12.0 | – | GLN92 NE2-H, 2.15, 144.7 | 2 | H2···GLN92 HE22 1.66 |

#### Orientations: every mode of every seed

| Run | Glu106 | Receptor (hydrogen rule) | Seeds | Top pose: near-native / rotated pose / other beyond 2.0 Å / between | All modes: near-native / rotated pose / other beyond 2.0 Å / between | Seeds holding a mode of the rotated pose | Best score: near-native | Best score: rotated pose | Rotated pose behind the best near-native mode of its seed, min–max | Modes beyond 2.0 Å within 2 kcal mol⁻¹ of the seed best (of them the rotated pose) | Closest heavy atom of any mode to a water site (Å) |
|:--|:--|:--|--:|:--|:--|--:|--:|--:|:--|:--|--:|
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry: `T3_6RL9` | 20 | 20 / 0 / 0 / 0 | 41 / 40 / 22 / 35 | 20 | −8.843 | −8.823 | 0.014–0.070 | 57 (39) | 0.07 |
|  |  | water: `T3_6RL9_wat4` (hydrogen-bond network) | 20 | 20 / 0 / 0 / 0 | 39 / 0 / 1 / 11 | 0 | −9.238 | – | – | 0 (0) | 2.79 |
| sulfanilamide, AD4Zn, ETKDG (validation) | acid (reported) | dry: `T3_6RL9` | 10 | 5 / 5 / 0 / 0 | 20 / 20 / 5 / 15 | 10 | −8.874 | −8.888 | −0.039–0.020 | 23 (20) | 0.20 |
|  |  | water: `T3_6RL9_wat4` (hydrogen-bond network) | 10 | 10 / 0 / 0 / 0 | 20 / 0 / 0 / 6 | 0 | −9.295 | – | – | 0 (0) | 2.79 |
| sulfonamidate, AD4Zn, ETKDG (validation) | acid (reported) | dry: `T3_6RL9` | 10 | 3 / 0 / 0 / 7 | 9 / 20 / 8 / 21 | 10 | −8.138 | −8.068 | 0.061–0.097 | 28 (20) | 0.34 |
|  |  | water: `T3_6RL9_wat4` (hydrogen-bond network) | 10 | 10 / 0 / 0 / 0 | 18 / 0 / 0 / 4 | 0 | −8.142 | – | – | 0 (0) | 2.79 |
| sulfanilamide, Vina, ETKDG (validation) | acid (reported) | dry: `T3_6RL9` | 10 | 0 / 0 / 0 / 10 | 0 / 14 / 160 / 26 | 9 | – | −5.517 | – | 174 (14) | 0.23 |
|  |  | water: `T3_6RL9_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 0 / 0 / 185 / 15 | 0 | – | – | – | 155 (0) | 2.68 |
| sulfanilamide, Vinardo, ETKDG (validation) | acid (reported) | dry: `T3_6RL9` | 10 | 0 / 0 / 0 / 10 | 2 / 23 / 145 / 30 | 10 | −4.679 | −4.442 | 0.237–0.237 | 168 (23) | 0.18 |
|  |  | water: `T3_6RL9_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 5 / 0 / 179 / 15 | 0 | −5.160 | – | – | 139 (0) | 2.46 |
| sulfanilamide, Vina, DFT geometry (production) | acid (reported) | dry: `T3_6RL9` | 10 | 0 / 0 / 0 / 10 | 1 / 16 / 154 / 29 | 10 | −5.714 | −5.567 | 0.200–0.200 | 170 (16) | 0.22 |
|  |  | water: `T3_6RL9_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 0 / 0 / 182 / 16 | 0 | – | – | – | 137 (0) | 2.69 |
| sulfanilamide, Vinardo, DFT geometry (production) | acid (reported) | dry: `T3_6RL9` | 10 | 0 / 0 / 0 / 10 | 0 / 14 / 156 / 30 | 9 | – | −4.473 | – | 170 (14) | 0.29 |
|  |  | water: `T3_6RL9_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 1 / 0 / 184 / 15 | 0 | −5.078 | – | – | 128 (0) | 2.53 |
| sulfanilamide, AD4Zn, DFT geometry (production) | carboxylate | dry: `T3_6RL9_glu106GLU` | 20 | 15 / 5 / 0 / 0 | 42 / 44 / 22 / 29 | 20 | −8.886 | −8.897 | −0.017–0.050 | 61 (41) | 0.09 |
|  |  | water: `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 20 | 20 / 0 / 0 / 0 | 38 / 0 / 1 / 10 | 0 | −9.247 | – | – | 1 (0) | 2.78 |
| sulfanilamide, AD4Zn, ETKDG (validation) | carboxylate | dry: `T3_6RL9_glu106GLU` | 10 | 3 / 7 / 0 / 0 | 21 / 20 / 7 / 16 | 10 | −8.919 | −8.945 | −0.046–0.005 | 27 (20) | 0.19 |
|  |  | water: `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 10 | 10 / 0 / 0 / 0 | 17 / 0 / 0 / 7 | 0 | −9.300 | – | – | 0 (0) | 2.79 |
| sulfonamidate, AD4Zn, ETKDG (validation) | carboxylate | dry: `T3_6RL9_glu106GLU` | 10 | 7 / 0 / 0 / 3 | 16 / 20 / 30 / 14 | 10 | −7.770 | −7.708 | 0.030–0.095 | 30 (20) | 0.32 |
|  |  | water: `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 10 | 10 / 0 / 0 / 0 | 18 / 0 / 0 / 5 | 0 | −7.772 | – | – | 0 (0) | 2.79 |
| sulfanilamide, Vina, ETKDG (validation) | carboxylate | dry: `T3_6RL9_glu106GLU` | 10 | 0 / 0 / 0 / 10 | 0 / 19 / 156 / 25 | 10 | – | −5.517 | – | 175 (19) | 0.19 |
|  |  | water: `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 2 / 0 / 183 / 15 | 0 | −6.184 | – | – | 150 (0) | 2.68 |
| sulfanilamide, Vinardo, ETKDG (validation) | carboxylate | dry: `T3_6RL9_glu106GLU` | 10 | 0 / 0 / 0 / 10 | 1 / 19 / 148 / 31 | 10 | −4.582 | −4.457 | 0.151–0.151 | 167 (19) | 0.18 |
|  |  | water: `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 2 / 0 / 181 / 17 | 0 | −5.107 | – | – | 138 (0) | 2.52 |
| sulfanilamide, Vina, DFT geometry (production) | carboxylate | dry: `T3_6RL9_glu106GLU` | 10 | 0 / 0 / 0 / 10 | 0 / 14 / 157 / 29 | 10 | – | −5.547 | – | 171 (14) | 0.24 |
|  |  | water: `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 0 / 0 / 182 / 17 | 0 | – | – | – | 131 (0) | 2.70 |
| sulfanilamide, Vinardo, DFT geometry (production) | carboxylate | dry: `T3_6RL9_glu106GLU` | 10 | 0 / 0 / 0 / 10 | 1 / 20 / 153 / 26 | 10 | −4.668 | −4.443 | 0.245–0.245 | 173 (20) | 0.25 |
|  |  | water: `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 0 / 0 / 182 / 18 | 0 | – | – | – | 130 (0) | 2.49 |
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry: `T3_6RL9` | 20 | 20 / 0 / 0 / 0 | 41 / 40 / 22 / 35 | 20 | −8.843 | −8.823 | 0.014–0.070 | 57 (39) | 0.07 |
|  |  | water: `T3_6RL9_wat4rev` (network, ring reversed) | 20 | 20 / 0 / 0 / 0 | 39 / 0 / 1 / 7 | 0 | −9.249 | – | – | 0 (0) | 2.79 |
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry: `T3_6RL9` | 20 | 20 / 0 / 0 / 0 | 41 / 40 / 22 / 35 | 20 | −8.843 | −8.823 | 0.014–0.070 | 57 (39) | 0.07 |
|  |  | water: `T3_6RL9_wat4bare` (bare oxygen (precedent)) | 20 | 20 / 0 / 0 / 0 | 41 / 0 / 0 / 6 | 0 | −9.542 | – | – | 0 (0) | 2.51 |

#### Does the crystal pose still dock: top poses, contacts, gates

| Run | Glu106 | Receptor | RMSD top (Å) | best of seeds | f ≤ 1 Å (≤ 2 Å) | Score (kcal mol⁻¹) | Zn–N (Å) / pass | Thr199 OG1···N (Å) / pass | Best mode of each seed (Å), min–max | PoseBusters physics | Strict gate / 2.0 Å |
|:--|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry | 0.61 ± 0.04 | 0.56 | 1.00 (1.00) | −8.836 ± 0.004 | 2.20 / 1.00 | 2.55 / 1.00 | 0.47–0.66 | – | – |
|  |  | water, hydrogen-bond network | 0.51 ± 0.04 | 0.45 | 1.00 (1.00) | −9.232 ± 0.004 | 2.16 / 1.00 | 2.66 / 1.00 | 0.32–0.51 | – | – |
| sulfanilamide, AD4Zn, ETKDG (validation) | acid (reported) | dry | 2.07 ± 1.45 | 0.65 | 0.5 (0.5) | −8.862 ± 0.017 | 2.22 / 1.00 | 2.56 / 1.00 | 0.40–0.73 | 1.0 | fail / fail |
|  |  | water, hydrogen-bond network | 0.58 ± 0.03 | 0.53 | 1.0 (1.0) | −9.286 ± 0.005 | 2.14 / 1.00 | 2.68 / 1.00 | 0.48–0.58 | 1.0 | pass / pass |
| sulfonamidate, AD4Zn, ETKDG (validation) | acid (reported) | dry | 1.12 ± 0.32 | 0.61 | 0.3 (1.0) | −8.133 ± 0.006 | 2.21 / 1.00 | 2.60 / 1.00 | 0.60–1.29 | 1.0 | fail / pass |
|  |  | water, hydrogen-bond network | 0.61 ± 0.04 | 0.52 | 1.0 (1.0) | −8.131 ± 0.009 | 2.23 / 1.00 | 2.59 / 1.00 | 0.46–0.62 | 1.0 | pass / pass |
| sulfanilamide, Vina, ETKDG (validation) | acid (reported) | dry | 1.10 ± 0.01 | 1.09 | 0.0 (1.0) | −5.870 ± 0.012 | 2.84 / 0.00 | 4.76 / 0.00 | 1.09–1.11 | 1.0 | fail / fail |
|  |  | water, hydrogen-bond network | 1.10 ± 0.01 | 1.09 | 0.0 (1.0) | −6.279 ± 0.061 | 2.93 / 0.00 | 4.81 / 0.00 | 1.09–1.11 | 1.0 | fail / fail |
| sulfanilamide, Vinardo, ETKDG (validation) | acid (reported) | dry | 1.13 ± 0.00 | 1.12 | 0.0 (1.0) | −5.075 ± 0.010 | 2.48 / 0.00 | 4.45 / 0.00 | 0.35–1.13 | 1.0 | fail / fail |
|  |  | water, hydrogen-bond network | 1.12 ± 0.00 | 1.11 | 0.0 (1.0) | −5.498 ± 0.013 | 2.50 / 0.00 | 4.47 / 0.00 | 0.26–1.12 | 1.0 | fail / fail |
| sulfanilamide, Vina, DFT geometry (production) | acid (reported) | dry | 1.08 ± 0.00 | 1.08 | 0.0 (1.0) | −5.880 ± 0.012 | 2.75 / 0.00 | 4.68 / 0.00 | 0.33–1.08 | – | – |
|  |  | water, hydrogen-bond network | 1.08 ± 0.00 | 1.08 | 0.0 (1.0) | −6.348 ± 0.004 | 2.75 / 0.00 | 4.67 / 0.00 | 1.08–1.08 | – | – |
| sulfanilamide, Vinardo, DFT geometry (production) | acid (reported) | dry | 1.11 ± 0.00 | 1.11 | 0.0 (1.0) | −5.077 ± 0.004 | 2.40 / 0.00 | 4.35 / 0.00 | 1.10–1.12 | – | – |
|  |  | water, hydrogen-bond network | 1.10 ± 0.00 | 1.10 | 0.0 (1.0) | −5.537 ± 0.004 | 2.43 / 0.00 | 4.37 / 0.00 | 0.25–1.11 | – | – |
| sulfanilamide, AD4Zn, DFT geometry (production) | carboxylate | dry | 1.34 ± 1.30 | 0.41 | 0.75 (0.75) | −8.882 ± 0.007 | 2.21 / 1.00 | 2.55 / 1.00 | 0.41–0.68 | – | – |
|  |  | water, hydrogen-bond network | 0.51 ± 0.03 | 0.47 | 1.00 (1.00) | −9.241 ± 0.004 | 2.16 / 1.00 | 2.65 / 1.00 | 0.30–0.56 | – | – |
| sulfanilamide, AD4Zn, ETKDG (validation) | carboxylate | dry | 2.70 ± 1.38 | 0.67 | 0.3 (0.3) | −8.911 ± 0.021 | 2.23 / 1.00 | 2.57 / 1.00 | 0.66–0.73 | 1.0 | fail / fail |
|  |  | water, hydrogen-bond network | 0.60 ± 0.03 | 0.55 | 1.0 (1.0) | −9.291 ± 0.004 | 2.15 / 1.00 | 2.67 / 1.00 | 0.43–0.62 | 1.0 | pass / pass |
| sulfonamidate, AD4Zn, ETKDG (validation) | carboxylate | dry | 0.86 ± 0.33 | 0.58 | 0.7 (1.0) | −7.764 ± 0.003 | 2.22 / 1.00 | 2.60 / 1.00 | 0.45–0.71 | 1.0 | pass / pass |
|  |  | water, hydrogen-bond network | 0.58 ± 0.04 | 0.53 | 1.0 (1.0) | −7.764 ± 0.006 | 2.23 / 1.00 | 2.59 / 1.00 | 0.52–0.64 | 1.0 | pass / pass |
| sulfanilamide, Vina, ETKDG (validation) | carboxylate | dry | 1.09 ± 0.01 | 1.09 | 0.0 (1.0) | −5.871 ± 0.011 | 2.83 / 0.00 | 4.77 / 0.00 | 1.09–1.10 | 1.0 | fail / fail |
|  |  | water, hydrogen-bond network | 1.10 ± 0.01 | 1.09 | 0.0 (1.0) | −6.302 ± 0.054 | 2.90 / 0.00 | 4.80 / 0.00 | 0.33–1.11 | 1.0 | fail / fail |
| sulfanilamide, Vinardo, ETKDG (validation) | carboxylate | dry | 1.12 ± 0.00 | 1.12 | 0.0 (1.0) | −5.078 ± 0.010 | 2.48 / 0.00 | 4.45 / 0.00 | 0.46–1.13 | 1.0 | fail / fail |
|  |  | water, hydrogen-bond network | 1.12 ± 0.01 | 1.11 | 0.0 (1.0) | −5.498 ± 0.012 | 2.50 / 0.00 | 4.47 / 0.00 | 0.26–1.12 | 1.0 | fail / fail |
| sulfanilamide, Vina, DFT geometry (production) | carboxylate | dry | 1.08 ± 0.00 | 1.08 | 0.0 (1.0) | −5.878 ± 0.012 | 2.75 / 0.00 | 4.68 / 0.00 | 1.08–1.08 | – | – |
|  |  | water, hydrogen-bond network | 1.09 ± 0.00 | 1.08 | 0.0 (1.0) | −6.348 ± 0.004 | 2.74 / 0.00 | 4.66 / 0.00 | 1.08–1.09 | – | – |
| sulfanilamide, Vinardo, DFT geometry (production) | carboxylate | dry | 1.11 ± 0.00 | 1.11 | 0.0 (1.0) | −5.076 ± 0.003 | 2.40 / 0.00 | 4.35 / 0.00 | 0.37–1.11 | – | – |
|  |  | water, hydrogen-bond network | 1.11 ± 0.00 | 1.10 | 0.0 (1.0) | −5.541 ± 0.003 | 2.42 / 0.00 | 4.37 / 0.00 | 1.10–1.11 | – | – |
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry | 0.61 ± 0.04 | 0.56 | 1.00 (1.00) | −8.836 ± 0.004 | 2.20 / 1.00 | 2.55 / 1.00 | 0.47–0.66 | – | – |
|  |  | water, network, ring reversed | 0.51 ± 0.03 | 0.45 | 1.00 (1.00) | −9.242 ± 0.003 | 2.16 / 1.00 | 2.65 / 1.00 | 0.37–0.53 | – | – |
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry | 0.61 ± 0.04 | 0.56 | 1.00 (1.00) | −8.836 ± 0.004 | 2.20 / 1.00 | 2.55 / 1.00 | 0.47–0.66 | – | – |
|  |  | water, bare oxygen (precedent) | 0.52 ± 0.04 | 0.42 | 1.00 (1.00) | −9.538 ± 0.003 | 2.12 / 1.00 | 2.70 / 1.00 | 0.36–0.53 | – | – |

#### What the search finds besides the crystal pose: families of modes beyond 1.0 Å (2.0 Å clusters, all seeds pooled)

| Run | Glu106 | Receptor | Family | Modes (seeds; top poses) | Best score | Behind the best mode of the run | RMSD to the crystal ligand (Å) | RMSD to the rotated reference (Å) | Zn–N (Å) | Aniline N···Asn67 OD1 (Å) | Nearest water site (Å) |
|:--|:--|:--|:--|:--|--:|--:|--:|--:|--:|--:|--:|
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry | 1: rotated pose | 40 (20; 0) | −8.823 | 0.020 | 3.57 | 0.00 | 2.24 | 2.57 | 0.19 |
|  |  |  | 2: between | 39 (20; 0) | −8.404 | 0.439 | 1.17 | 2.75 | 2.19 | 5.79 | 1.49 |
|  |  |  | 3: other, beyond 2.0 A | 8 (7; 0) | −7.657 | 1.186 | 2.14 | 3.64 | 2.00 | 7.52 | 0.90 |
|  |  |  | 4: other, beyond 2.0 A | 9 (9; 0) | −7.635 | 1.208 | 3.10 | 4.99 | 1.95 | 10.42 | 1.18 |
|  |  |  | 5: other, beyond 2.0 A | 1 (1; 0) | −6.115 | 2.728 | 2.91 | 2.15 | 1.86 | 4.59 | 0.55 |
|  |  | water, hydrogen-bond network | 1: between | 11 (11; 0) | −8.857 | 0.381 | 1.57 | 3.61 | 2.28 | 7.12 | 2.79 |
|  |  |  | 2: between | 1 (1; 0) | −7.721 | 1.517 | 1.65 | 3.74 | 2.23 | 8.88 | 2.80 |
| sulfanilamide, AD4Zn, ETKDG (validation) | acid (reported) | dry | 1: rotated pose | 21 (10; 5) | −8.888 | 0.000 | 3.37 | 0.26 | 2.23 | 2.70 | 0.23 |
|  |  |  | 2: between | 14 (8; 0) | −8.400 | 0.488 | 1.15 | 2.77 | 2.13 | 5.86 | 1.61 |
|  |  |  | 3: other, beyond 2.0 A | 2 (2; 0) | −7.728 | 1.160 | 3.16 | 5.09 | 2.21 | 10.62 | 1.32 |
|  |  |  | 4: other, beyond 2.0 A | 3 (3; 0) | −7.443 | 1.445 | 2.17 | 3.72 | 2.00 | 7.74 | 0.91 |
|  |  | water, hydrogen-bond network | 1: between | 6 (5; 0) | −8.509 | 0.786 | 1.58 | 3.58 | 2.32 | 6.95 | 2.79 |
| sulfonamidate, AD4Zn, ETKDG (validation) | acid (reported) | dry | 1: between | 19 (9; 7) | −8.140 | 0.000 | 1.37 | 3.86 | 2.18 | 7.21 | 2.48 |
|  |  |  | 2: rotated pose | 20 (10; 0) | −8.068 | 0.072 | 3.32 | 0.32 | 2.24 | 2.81 | 0.41 |
|  |  |  | 3: between | 3 (3; 0) | −7.310 | 0.830 | 1.67 | 3.06 | 2.13 | 7.41 | 1.89 |
|  |  |  | 4: other, beyond 2.0 A | 4 (4; 0) | −6.937 | 1.203 | 3.12 | 5.00 | 1.99 | 10.45 | 1.31 |
|  |  |  | 5: other, beyond 2.0 A | 3 (3; 0) | −6.875 | 1.265 | 2.20 | 3.54 | 2.36 | 7.23 | 0.84 |
|  |  | water, hydrogen-bond network | 1: between | 3 (3; 0) | −7.656 | 0.486 | 1.59 | 3.55 | 2.28 | 7.00 | 2.79 |
|  |  |  | 2: between | 1 (1; 0) | −6.910 | 1.232 | 1.52 | 3.70 | 2.25 | 8.73 | 2.80 |
| sulfanilamide, AD4Zn, DFT geometry (production) | carboxylate | dry | 1: rotated pose | 45 (20; 5) | −8.897 | 0.000 | 3.51 | 0.08 | 2.24 | 2.58 | 0.16 |
|  |  |  | 2: between | 26 (18; 0) | −8.445 | 0.452 | 1.02 | 2.81 | 2.17 | 5.97 | 1.92 |
|  |  |  | 3: between | 1 (1; 0) | −8.021 | 0.876 | 1.85 | 3.28 | 2.11 | 8.00 | 1.98 |
|  |  |  | 4: other, beyond 2.0 A | 12 (12; 0) | −7.722 | 1.175 | 3.07 | 4.98 | 2.04 | 10.39 | 1.23 |
|  |  |  | 5: other, beyond 2.0 A | 11 (11; 0) | −7.644 | 1.253 | 2.17 | 3.61 | 2.12 | 7.44 | 0.84 |
|  |  | water, hydrogen-bond network | 1: between | 9 (9; 0) | −8.858 | 0.389 | 1.59 | 3.53 | 2.31 | 6.99 | 2.79 |
|  |  |  | 2: other, beyond 2.0 A | 1 (1; 0) | −7.770 | 1.477 | 2.16 | 4.88 | 2.21 | 9.11 | 3.01 |
|  |  |  | 3: between | 1 (1; 0) | −7.621 | 1.626 | 1.57 | 3.65 | 2.16 | 8.62 | 2.78 |
| sulfanilamide, AD4Zn, ETKDG (validation) | carboxylate | dry | 1: rotated pose | 21 (10; 7) | −8.945 | 0.000 | 3.41 | 0.22 | 2.24 | 2.71 | 0.20 |
|  |  |  | 2: between | 14 (10; 0) | −8.408 | 0.537 | 1.09 | 2.73 | 2.15 | 5.95 | 1.75 |
|  |  |  | 3: other, beyond 2.0 A | 5 (4; 0) | −7.669 | 1.276 | 3.16 | 5.08 | 2.20 | 10.61 | 1.33 |
|  |  |  | 4: other, beyond 2.0 A | 3 (3; 0) | −7.358 | 1.587 | 2.14 | 3.65 | 1.94 | 7.54 | 0.87 |
|  |  | water, hydrogen-bond network | 1: between | 7 (7; 0) | −8.483 | 0.817 | 1.59 | 3.51 | 2.31 | 6.84 | 2.79 |
| sulfonamidate, AD4Zn, ETKDG (validation) | carboxylate | dry | 1: between | 6 (4; 3) | −7.767 | 0.003 | 1.39 | 3.93 | 2.15 | 7.34 | 2.57 |
|  |  |  | 2: rotated pose | 21 (10; 0) | −7.708 | 0.062 | 3.28 | 0.38 | 2.24 | 2.95 | 0.39 |
|  |  |  | 3: between | 8 (8; 0) | −7.041 | 0.729 | 1.34 | 3.22 | 2.15 | 7.60 | 2.16 |
|  |  |  | 4: other, beyond 2.0 A | 5 (5; 0) | −6.584 | 1.186 | 3.11 | 4.99 | 2.03 | 10.40 | 1.33 |
|  |  |  | 5: other, beyond 2.0 A | 6 (5; 0) | −6.558 | 1.212 | 2.15 | 3.62 | 2.01 | 7.45 | 0.84 |
|  |  |  | 6: other, beyond 2.0 A | 18 (10; 0) | −4.869 | 2.901 | 5.51 | 4.48 | 6.12 | 8.38 | 0.57 |
|  |  | water, hydrogen-bond network | 1: between | 3 (3; 0) | −7.353 | 0.419 | 1.67 | 3.67 | 2.30 | 7.15 | 2.79 |
|  |  |  | 2: between | 2 (2; 0) | −7.034 | 0.738 | 1.10 | 3.65 | 2.01 | 8.40 | 3.17 |
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry | 1: rotated pose | 40 (20; 0) | −8.823 | 0.020 | 3.57 | 0.00 | 2.24 | 2.57 | 0.19 |
|  |  |  | 2: between | 39 (20; 0) | −8.404 | 0.439 | 1.17 | 2.75 | 2.19 | 5.79 | 1.49 |
|  |  |  | 3: other, beyond 2.0 A | 8 (7; 0) | −7.657 | 1.186 | 2.14 | 3.64 | 2.00 | 7.52 | 0.90 |
|  |  |  | 4: other, beyond 2.0 A | 9 (9; 0) | −7.635 | 1.208 | 3.10 | 4.99 | 1.95 | 10.42 | 1.18 |
|  |  |  | 5: other, beyond 2.0 A | 1 (1; 0) | −6.115 | 2.728 | 2.91 | 2.15 | 1.86 | 4.59 | 0.55 |
|  |  | water, network, ring reversed | 1: between | 8 (8; 0) | −8.822 | 0.427 | 1.59 | 3.58 | 2.29 | 7.10 | 2.79 |
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry | 1: rotated pose | 40 (20; 0) | −8.823 | 0.020 | 3.57 | 0.00 | 2.24 | 2.57 | 0.19 |
|  |  |  | 2: between | 39 (20; 0) | −8.404 | 0.439 | 1.17 | 2.75 | 2.19 | 5.79 | 1.49 |
|  |  |  | 3: other, beyond 2.0 A | 8 (7; 0) | −7.657 | 1.186 | 2.14 | 3.64 | 2.00 | 7.52 | 0.90 |
|  |  |  | 4: other, beyond 2.0 A | 9 (9; 0) | −7.635 | 1.208 | 3.10 | 4.99 | 1.95 | 10.42 | 1.18 |
|  |  |  | 5: other, beyond 2.0 A | 1 (1; 0) | −6.115 | 2.728 | 2.91 | 2.15 | 1.86 | 4.59 | 0.55 |
|  |  | water, bare oxygen (precedent) | 1: between | 4 (4; 0) | −9.144 | 0.398 | 1.60 | 3.59 | 2.25 | 7.08 | 2.79 |
|  |  |  | 2: between | 2 (2; 0) | −8.219 | 1.323 | 1.62 | 3.72 | 2.25 | 8.76 | 2.80 |

#### Stored modes of the dry runs scored in the maps with the waters (no search)

| Maps | Stored dry run | Near-native modes: n, score in the dry maps → in the water maps (change) | Modes of the rotated pose: n, score in the water maps, mean (min to max) |
|:--|:--|:--|:--|
| `T3_6RL9_wat4` | `T3_6RL9`, sulfanilamide, AD4Zn, DFT geometry (production) | 20, −8.836 → −8.925 (−0.100 to −0.070) | 20, +6840 (+6000 to +7432) |
| `T3_6RL9_glu106GLU_wat4` | `T3_6RL9_glu106GLU`, sulfanilamide, AD4Zn, DFT geometry (production) | 20, −8.878 → −8.962 (−0.098 to −0.067) | 20, +6818 (+6472 to +7360) |
| `T3_6RL9_wat4rev` | `T3_6RL9`, sulfanilamide, AD4Zn, DFT geometry (production) | 20, −8.836 → −8.926 (−0.102 to −0.073) | 20, +6867 (+5846 to +7209) |
| `T3_6RL9_wat4bare` | `T3_6RL9`, sulfanilamide, AD4Zn, DFT geometry (production) | 20, −8.836 → −9.175 (−0.358 to −0.313) | 20, +5877 (+5346 to +6858) |
| `T3_6RL9_wat4fixed` | `T3_6RL9`, sulfanilamide, AD4Zn, DFT geometry (production) | 20, −8.836 → −8.923 (−0.103 to −0.068) | 20, +7151 (+6410 to +7644) |
