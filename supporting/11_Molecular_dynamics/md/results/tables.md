# Results of the MD campaign: the long tables

Written with `md/RESULTS.md` by `md/production/assemble/make_results_md.py` on 11 Oct 2026, 11:50 AM IST; every number is read from the files named.

## A. The six N–Zn–N angles: mean of each run and its departure from 6RL9 (°)

Source: `gate_GV.GV3_zinc.criteria` of each `summary.json` (first program; the second gives the same). Reference: the angle in 6RL9 as the specification prints it (section 11).

| Run | N2–Zn–NE2(His94) | N2–Zn–NE2(His96) | N2–Zn–ND1(His119) | NE2(His94)–Zn–NE2(His96) | NE2(His94)–Zn–ND1(His119) | NE2(His96)–Zn–ND1(His119) | S–N2–Zn |
|---|---|---|---|---|---|---|---|
| 6RL9 | 108.9 | 114.3 | 116.8 | 103.8 | 113.6 | 98.4 | 117.6 |
| X1 | 105.9 (-3.0) | 115.7 (+1.4) | 118.4 (+1.6) | 102.5 (-1.3) | 114.7 (+1.1) | 98.2 (-0.2) | 113.7 |
| X2 | 105.7 (-3.2) | 115.7 (+1.4) | 118.6 (+1.8) | 102.8 (-1.0) | 114.5 (+0.9) | 98.3 (-0.1) | 114.0 |
| X3 | 105.8 (-3.1) | 115.9 (+1.6) | 118.5 (+1.7) | 102.5 (-1.3) | 114.6 (+1.0) | 98.2 (-0.2) | 114.1 |
| N1 | 105.3 (-3.6) | 116.0 (+1.7) | 118.5 (+1.7) | 102.8 (-1.0) | 114.6 (+1.0) | 98.4 (+0.0) | 114.0 |
| N2 | 105.3 (-3.6) | 116.1 (+1.8) | 118.6 (+1.8) | 102.8 (-1.0) | 114.6 (+1.0) | 98.2 (-0.2) | 113.7 |
| N3 | 105.4 (-3.5) | 115.3 (+1.0) | 118.8 (+2.0) | 102.9 (-0.9) | 114.4 (+0.8) | 98.8 (+0.4) | 114.5 |
| R1 | 105.7 (-3.2) | 116.1 (+1.8) | 118.2 (+1.4) | 102.4 (-1.4) | 114.8 (+1.2) | 98.4 (+0.0) | 113.9 |
| R2 | 105.5 (-3.4) | 116.0 (+1.7) | 118.6 (+1.8) | 102.4 (-1.4) | 114.5 (+0.9) | 98.5 (+0.1) | 114.0 |
| R3 | 105.5 (-3.4) | 116.0 (+1.7) | 118.5 (+1.7) | 102.6 (-1.2) | 114.6 (+1.0) | 98.3 (-0.1) | 114.1 |
| 2s/R1 | 113.6 (+4.7) | 107.6 (-6.7) | 119.0 (+2.2) | 101.6 (-2.2) | 114.4 (+0.8) | 96.8 (-1.6) | 109.7 |
| 2e_n2fit/R1 | 109.0 (+0.1) | 113.5 (-0.8) | 117.5 (+0.7) | 102.4 (-1.4) | 115.0 (+1.4) | 97.7 (-0.7) | 112.1 |
| 2e/X0 | 105.6 (-3.3) | 115.8 (+1.5) | 118.3 (+1.5) | 102.4 (-1.4) | 114.7 (+1.1) | 98.7 (+0.3) | 114.0 |
| 2s/X0 | 113.2 (+4.3) | 108.0 (-6.3) | 118.5 (+1.7) | 101.7 (-2.1) | 115.0 (+1.4) | 96.6 (-1.8) | 109.3 |

## B. Equilibration: steps, restraints and mean temperature of the state log (K)

Source: `results/run_health.json` (manifest of each run and its `log_<step>.csv`).

| Step | Length (ps) | Time step (fs) | Ensemble | Atoms held; force constant (kJ mol⁻¹ nm⁻²) | 2s/X0 | 2e/X0 | X1 | N1 | R1 | X2 | N2 | R2 | X3 | N3 | R3 | 2s/R1 | 2e_n2fit/R1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| eq01 (heat) | 100 | 2 | NVT | solute_heavy; 1000 | 185.7 | 186.0 | 185.3 | 186.0 | 186.0 | 185.8 | 185.7 | 185.4 | 185.5 | 185.7 | 186.1 | 185.9 | 186.0 |
| eq02 (density) | 300 | 2 | NPT | solute_heavy; 1000 | 300.6 | 300.2 | 299.8 | 300.7 | 300.0 | 300.4 | 300.9 | 300.5 | 300.7 | 300.3 | 300.7 | 300.3 | 300.5 |
| eq03 (loosen) | 300 | 2 | NPT | solute_heavy; 250 | 300.2 | 300.6 | 300.9 | 300.4 | 299.9 | 300.4 | 300.2 | 300.5 | 300.0 | 300.8 | 300.4 | 300.6 | 300.3 |
| eq03x1 (loosen, extension 1 of at most 4) | 500 | 2 | NPT | solute_heavy; 250 | not run | not run | not run | not run | 300.0 | not run | not run | not run | not run | not run | not run | not run | not run |
| eq04 (release_side_chains_zinc_ligand) | 400 | 4 | NPT | protein_backbone; 250 | 301.2 | 301.0 | 300.8 | 301.2 | 301.1 | 300.7 | 301.0 | 301.1 | 300.8 | 300.6 | 301.3 | 300.8 | 300.5 |
| eq05 (loosen_backbone) | 400 | 4 | NPT | protein_backbone; 50 | 300.5 | 300.8 | 300.8 | 300.7 | 300.5 | 300.8 | 300.7 | 300.7 | 300.8 | 301.0 | 301.1 | 301.0 | 301.1 |
| eq06 (free) | 500 | 4 | NPT | none; 0 | 300.7 | 301.2 | 300.9 | 300.7 | 300.8 | 301.0 | 300.7 | 301.3 | 300.7 | 300.9 | 300.5 | 301.0 | 300.9 |

## C. The production in blocks of 20 ns

Frames of each state (crystal-like / rotated-like / intermediate, of 2,000), frames with the sign at N2 negative, frames with the hydrogen bond N2–H···OG1(Thr199), mean ligand RMSD to the crystal pose (Å). Source: `results/run_events.json`, counted from the production series of both programs.

| Run | 0 to 20 ns | 20 to 40 ns | 40 to 60 ns | 60 to 80 ns | 80 to 100 ns |
|---|---|---|---|---|---|
| X1 | 1990 / 0 / 10; N2− 0; HB 1998; 0.57 | 1983 / 0 / 17; N2− 0; HB 1998; 0.52 | 1992 / 0 / 8; N2− 0; HB 1995; 0.53 | 1993 / 0 / 7; N2− 0; HB 1998; 0.56 | 1995 / 0 / 5; N2− 0; HB 1997; 0.60 |
| X2 | 1990 / 0 / 10; N2− 0; HB 1999; 0.54 | 1937 / 0 / 63; N2− 0; HB 1997; 0.56 | 1946 / 0 / 54; N2− 0; HB 1997; 0.53 | 1998 / 0 / 2; N2− 0; HB 1996; 0.52 | 1993 / 0 / 7; N2− 0; HB 1997; 0.55 |
| X3 | 1990 / 0 / 10; N2− 0; HB 1999; 0.53 | 1990 / 0 / 10; N2− 0; HB 1998; 0.53 | 1998 / 0 / 2; N2− 0; HB 1999; 0.51 | 1997 / 0 / 3; N2− 0; HB 1996; 0.55 | 1975 / 0 / 25; N2− 0; HB 1991; 0.61 |
| N1 | 1996 / 0 / 4; N2− 0; HB 1999; 0.53 | 1998 / 0 / 2; N2− 0; HB 1997; 0.52 | 1996 / 0 / 4; N2− 0; HB 1996; 0.52 | 1999 / 0 / 1; N2− 0; HB 1998; 0.51 | 1991 / 0 / 9; N2− 0; HB 1992; 0.57 |
| N2 | 1997 / 0 / 3; N2− 0; HB 1995; 0.52 | 1974 / 0 / 26; N2− 0; HB 2000; 0.58 | 1986 / 0 / 14; N2− 0; HB 1996; 0.58 | 1989 / 0 / 11; N2− 0; HB 1997; 0.53 | 1998 / 0 / 2; N2− 0; HB 1997; 0.50 |
| N3 | 1996 / 0 / 4; N2− 0; HB 1998; 0.55 | 1996 / 0 / 4; N2− 0; HB 1997; 0.55 | 1997 / 0 / 3; N2− 0; HB 2000; 0.54 | 1994 / 0 / 6; N2− 0; HB 1983; 0.82 | 1988 / 0 / 12; N2− 0; HB 1999; 0.58 |
| R1 | 1980 / 0 / 20; N2− 0; HB 1998; 0.56 | 1994 / 0 / 6; N2− 0; HB 2000; 0.60 | 1962 / 0 / 38; N2− 0; HB 1997; 0.55 | 1984 / 0 / 16; N2− 0; HB 1999; 0.50 | 1993 / 0 / 7; N2− 0; HB 1997; 0.49 |
| R2 | 1996 / 0 / 4; N2− 0; HB 1990; 0.66 | 1998 / 0 / 2; N2− 0; HB 1999; 0.52 | 1980 / 0 / 20; N2− 0; HB 1998; 0.53 | 1996 / 0 / 4; N2− 0; HB 1998; 0.51 | 1993 / 0 / 7; N2− 0; HB 1997; 0.51 |
| R3 | 1962 / 0 / 38; N2− 0; HB 1995; 0.58 | 1999 / 0 / 1; N2− 0; HB 2000; 0.52 | 1990 / 0 / 10; N2− 0; HB 1997; 0.52 | 2000 / 0 / 0; N2− 0; HB 1998; 0.52 | 2000 / 0 / 0; N2− 0; HB 1999; 0.53 |
| 2s/R1 | 1979 / 0 / 21; N2− 0; HB 1998; 0.70 | 1976 / 1 / 23; N2− 0; HB 1998; 0.74 | 1964 / 0 / 36; N2− 1; HB 1998; 0.92 | 1968 / 0 / 32; N2− 1; HB 1997; 0.87 | 1989 / 0 / 11; N2− 0; HB 1997; 0.78 |
| 2e_n2fit/R1 | 1820 / 23 / 157; N2− 120; HB 1961; 0.90 | 1911 / 11 / 78; N2− 44; HB 1987; 0.83 | 1796 / 50 / 154; N2− 134; HB 1959; 1.05 | 1708 / 126 / 166; N2− 169; HB 1967; 1.00 | 1487 / 183 / 330; N2− 235; HB 1955; 1.19 |

## D. The frames around the change of sign at N2 in R1, R2 and R3

From two frames before the first frame with the sign at N2 positive to two frames after the first crystal-like frame. Source: `results/run_events.json` (`frames_around_the_conversion`), from `timeseries_equil_solute.csv`.

**R1**

| Time from the release (ps) | State | Signed volume at N2 | N2–H···OG1 present | N2···OG1 (Å) | N2–H···OG1 (°) | Ligand RMSD (Å) | d1 (Å) | d2 (Å) |
|---|---|---|---|---|---|---|---|---|
| 240 | R | -0.564 | yes | 2.67 | 175 | 3.59 | 9.69 | 3.07 |
| 250 | R | -0.504 | yes | 2.80 | 158 | 3.81 | 10.24 | 2.99 |
| 260 | R | +0.532 | yes | 2.79 | 162 | 2.75 | 7.35 | 3.00 |
| 270 | C | +0.663 | yes | 2.73 | 164 | 1.36 | 4.78 | 5.85 |
| 280 | C | +0.792 | yes | 3.06 | 154 | 1.58 | 5.07 | 5.29 |
| 290 | C | +0.607 | yes | 2.87 | 150 | 1.31 | 4.38 | 5.83 |

**R2**

| Time from the release (ps) | State | Signed volume at N2 | N2–H···OG1 present | N2···OG1 (Å) | N2–H···OG1 (°) | Ligand RMSD (Å) | d1 (Å) | d2 (Å) |
|---|---|---|---|---|---|---|---|---|
| 340 | R | -0.623 | yes | 2.97 | 167 | 3.79 | 10.04 | 2.92 |
| 350 | R | -0.453 | yes | 2.89 | 161 | 3.94 | 10.36 | 2.80 |
| 360 | I | +0.627 | yes | 2.96 | 167 | 2.38 | 6.54 | 4.16 |
| 370 | C | +0.759 | yes | 2.95 | 163 | 1.22 | 3.89 | 5.96 |
| 380 | C | +0.710 | yes | 2.94 | 159 | 0.78 | 3.60 | 7.13 |
| 390 | C | +0.727 | yes | 2.86 | 166 | 0.77 | 4.12 | 7.21 |

**R3**

| Time from the release (ps) | State | Signed volume at N2 | N2–H···OG1 present | N2···OG1 (Å) | N2–H···OG1 (°) | Ligand RMSD (Å) | d1 (Å) | d2 (Å) |
|---|---|---|---|---|---|---|---|---|
| 120 | R | -0.486 | yes | 2.80 | 177 | 3.93 | 10.30 | 3.03 |
| 130 | R | -0.628 | yes | 2.99 | 167 | 3.86 | 9.63 | 3.08 |
| 140 | R | +0.106 | yes | 2.84 | 157 | 2.84 | 7.81 | 3.09 |
| 150 | I | +0.780 | yes | 2.84 | 166 | 1.80 | 5.40 | 4.75 |
| 160 | C | +0.840 | yes | 2.98 | 169 | 1.44 | 4.35 | 5.33 |
| 170 | C | +0.650 | yes | 2.87 | 162 | 1.14 | 4.11 | 5.75 |
| 180 | C | +0.684 | yes | 2.84 | 172 | 1.17 | 4.64 | 5.32 |

## E. Water in the place of the crystal ring through the equilibration of the runs from start R

All-atom frames every 100 ps; time from the start of the production (ns); count of water oxygens within 2.5 Å of the crystal positions of the six ring carbons (`n_water_ring_place`), distance of the nearest one (Å), and the state of the solute frame of the same time. Source: `results/run_events.json`.

**R1**: -2.4: 0 (3.35, R); -2.3: 2 (1.19, R); -2.2: 3 (1.86, R); -2.1: 1 (1.08, R); -2.0: 2 (2.35, R); -1.9: 2 (0.95, R); -1.8: 0 (3.29, R); -1.7: 3 (1.25, R); -1.6: 1 (1.58, R); -1.5: 4 (1.16, R); -1.4: 1 (1.78, R); -1.3: 2 (0.75, R); -1.2: 2 (0.54, R); -1.1: 3 (0.48, R); -1.0: 0 (4.12, C); -0.9: 1 (2.19, C); -0.8: 1 (1.84, C); -0.7: 2 (2.13, C); -0.6: 1 (2.18, C); -0.5: 1 (1.02, C); -0.4: 1 (2.04, C); -0.3: 0 (2.66, C); -0.2: 0 (3.17, C); -0.1: 0 (3.53, C); 0.0: 0 (3.60, C)

**R2**: -1.9: 1 (2.02, R); -1.8: 3 (0.88, R); -1.7: 3 (0.94, R); -1.6: 3 (0.35, R); -1.5: 2 (0.75, R); -1.4: 2 (0.68, R); -1.3: 2 (0.89, R); -1.2: 4 (0.62, R); -1.1: 2 (0.91, R); -1.0: 0 (2.90, R); -0.9: 0 (3.86, C); -0.8: 0 (2.97, C); -0.7: 0 (3.39, C); -0.6: 0 (3.16, C); -0.5: 0 (3.52, C); -0.4: 0 (3.59, C); -0.3: 0 (3.23, C); -0.2: 0 (3.41, C); -0.1: 0 (2.75, C); 0.0: 0 (3.35, C)

**R3**: -1.9: 1 (2.19, R); -1.8: 3 (1.15, R); -1.7: 2 (1.12, R); -1.6: 1 (0.82, R); -1.5: 2 (2.34, R); -1.4: 1 (1.42, R); -1.3: 1 (2.31, R); -1.2: 1 (2.04, R); -1.1: 0 (2.77, C); -1.0: 0 (3.89, C); -0.9: 0 (2.95, C); -0.8: 0 (3.06, C); -0.7: 0 (3.35, C); -0.6: 0 (3.26, C); -0.5: 0 (3.13, C); -0.4: 0 (3.40, C); -0.3: 0 (3.10, C); -0.2: 0 (3.50, C); -0.1: 0 (2.86, C); 0.0: 0 (3.25, C)

**2s/R1**: -1.9: 1 (1.66, R); -1.8: 1 (2.25, R); -1.7: 1 (0.86, R); -1.6: 1 (1.67, R); -1.5: 1 (1.15, R); -1.4: 2 (1.98, R); -1.3: 1 (2.43, R); -1.2: 0 (3.87, I); -1.1: 0 (2.75, C); -1.0: 1 (1.80, C); -0.9: 1 (2.15, C); -0.8: 0 (2.65, C); -0.7: 1 (1.99, I); -0.6: 1 (1.79, I); -0.5: 1 (2.29, I); -0.4: 2 (1.82, C); -0.3: 1 (2.11, C); -0.2: 0 (3.50, C); -0.1: 0 (3.61, C); 0.0: 0 (3.23, C)

**2e_n2fit/R1**: -1.9: 0 (2.58, R); -1.8: 3 (0.49, R); -1.7: 2 (1.18, R); -1.6: 1 (1.28, R); -1.5: 1 (2.16, R); -1.4: 3 (0.83, R); -1.3: 1 (2.49, R); -1.2: 0 (2.93, C); -1.1: 0 (2.91, C); -1.0: 0 (3.27, C); -0.9: 0 (3.94, C); -0.8: 1 (2.36, C); -0.7: 0 (3.09, C); -0.6: 0 (3.33, C); -0.5: 0 (2.93, C); -0.4: 0 (3.63, C); -0.3: 0 (3.12, C); -0.2: 0 (3.60, C); -0.1: 0 (3.02, C); 0.0: 0 (3.27, C)

## F. Hydration check of the driver and its reproduction by the two programs

Source: `equilibration.stages[].driver_hydration_check` of each `summary.json`; `count_agrees` and `nearest_agrees` are true in both programs for every record.

| Run | Step | Driver: count; nearest water oxygen (Å); applies; decision | Compared; count agrees; nearest agrees |
|---|---|---|---|
| 2s/X0 | eq03 | 0; 3.103; false; go on | true; true; true |
| 2e/X0 | eq03 | 0; 3.186; false; go on | true; true; true |
| X1 | eq03 | 0; 3.263; false; go on | true; true; true |
| N1 | eq03 | 0; 3.057; false; go on | true; true; true |
| R1 | eq03 | 0; 3.293; true; extend | true; true; true |
| R1 | eq03x1 | 2; 0.745; true; go on | true; true; true |
| X2 | eq03 | 0; 3.203; false; go on | true; true; true |
| N2 | eq03 | 0; 3.319; false; go on | true; true; true |
| R2 | eq03 | 2; 0.885; true; go on | true; true; true |
| X3 | eq03 | 0; 3.259; false; go on | true; true; true |
| N3 | eq03 | 0; 2.846; false; go on | true; true; true |
| R3 | eq03 | 1; 2.312; true; go on | true; true; true |
| 2s/R1 | eq03 | 1; 2.429; true; go on | true; true; true |
| 2e_n2fit/R1 | eq03 | 1; 2.493; true; go on | true; true; true |

## G. Frames within the tolerance of a threshold, per run

Own reading of the raw frames (`md/production/site_checks.py`; `results/site_checks_solute.json`, `results/site_checks_full.json`), equilibration and production together. Tolerances of `md/analysis/SPEC.md`, section 14: 0.001 Å, 0.01°, 0.001 for the signed volume. This reading, the first program and the second class every one of these frames alike.

| Run | Solute frames read | d1 at 6.7 Å | d2 at 5.2 Å | N2···OG1 at 3.5 Å | N2–H···OG1 at 135° | nearer sulfonyl O at 2.6 Å | signed volume at N2 below 0.001 in size | Solute frames, any | All-atom frames read | nearest oxygen at 2.6 Å from the zinc | water at the 2.5 Å of the place of the ring | water at 3.5 Å of N1 | water at 3.5 Å of the ring with N1 | All-atom frames, any |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| X1 | 10,200 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 1,020 | 0 | 0 | 8 | 6 | 14 |
| X2 | 10,200 | 2 | 0 | 0 | 0 | 2 | 0 | 4 | 1,020 | 0 | 0 | 7 | 10 | 17 |
| X3 | 10,200 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 1,020 | 0 | 1 | 2 | 9 | 12 |
| N1 | 10,200 | 0 | 0 | 1 | 0 | 0 | 0 | 1 | 1,020 | 0 | 0 | 6 | 6 | 11 |
| N2 | 10,200 | 0 | 0 | 0 | 0 | 2 | 0 | 2 | 1,020 | 0 | 0 | 3 | 9 | 12 |
| N3 | 10,200 | 0 | 0 | 0 | 1 | 2 | 0 | 3 | 1,020 | 0 | 1 | 7 | 10 | 18 |
| R1 | 10,250 | 1 | 0 | 0 | 0 | 1 | 0 | 2 | 1,025 | 0 | 0 | 3 | 6 | 9 |
| R2 | 10,200 | 0 | 0 | 0 | 0 | 1 | 0 | 1 | 1,020 | 0 | 0 | 5 | 11 | 15 |
| R3 | 10,200 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1,020 | 0 | 0 | 10 | 11 | 21 |
| 2s/R1 | 10,200 | 0 | 0 | 0 | 0 | 7 | 0 | 7 | 1,020 | 2 | 1 | 9 | 8 | 20 |
| 2e_n2fit/R1 | 10,200 | 6 | 2 | 0 | 0 | 8 | 5 | 21 | 1,020 | 0 | 2 | 6 | 8 | 16 |
| 2e/X0 | 700 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 70 | 0 | 0 | 1 | 0 | 1 |
| 2s/X0 | 700 | 0 | 1 | 0 | 0 | 0 | 0 | 1 | 70 | 0 | 0 | 1 | 0 | 1 |
| all 13 | 113,650 | 11 | 3 | 1 | 2 | 23 | 5 | 45 | 11,365 | 2 | 5 | 68 | 94 | 167 |

## H. The places of four crystal waters and the water at N2, whole production

Own reading of the raw all-atom frames (`results/site_checks_full.json`): a water oxygen within 1.4 Å of the place of the crystal water after the fit of the Cα atoms, per cent of the 1,000 all-atom frames of the production; a water donating a hydrogen bond to N2 (O within 3.5 Å, O–H···N2 at least 135°); the same two for the frames in which the solute frame of the same time is crystal-like, rotated-like or intermediate. Waters 423, 516, 553 and 568 are the four removed to build start R.

| Run | 423 | 516 | 553 | 568 | water donating to N2 | all-atom frames C / R / I | place 553 occupied in C / R / I frames (%) | water in the place of the crystal ring in C / R / I frames (%) |
|---|---|---|---|---|---|---|---|---|
| X1 | 32.8 | 59.0 | 92.0 | 86.3 | 92.6 | 990 / 0 / 10 | 92.1 / none / 80.0 | 2.2 / none / 0.0 |
| X2 | 27.0 | 75.6 | 95.9 | 82.1 | 94.2 | 979 / 0 / 21 | 96.1 / none / 85.7 | 2.6 / none / 14.3 |
| X3 | 25.3 | 28.3 | 94.8 | 86.9 | 95.6 | 991 / 0 / 9 | 94.9 / none / 88.9 | 1.3 / none / 0.0 |
| N1 | 20.3 | 79.5 | 95.8 | 84.7 | 94.4 | 999 / 0 / 1 | 95.8 / none / 100.0 | 2.0 / none / 0.0 |
| N2 | 29.0 | 79.0 | 94.4 | 84.0 | 92.8 | 998 / 0 / 2 | 94.4 / none / 100.0 | 3.3 / none / 0.0 |
| N3 | 24.2 | 55.2 | 86.8 | 71.5 | 89.7 | 997 / 0 / 3 | 86.8 / none / 100.0 | 2.7 / none / 0.0 |
| R1 | 17.1 | 10.6 | 96.7 | 86.5 | 95.9 | 989 / 0 / 11 | 96.9 / none / 81.8 | 1.5 / none / 18.2 |
| R2 | 14.8 | 17.5 | 93.2 | 87.4 | 92.4 | 1,000 / 0 / 0 | 93.2 / none / none | 1.8 / none / none |
| R3 | 12.4 | 21.4 | 96.9 | 88.2 | 96.4 | 997 / 0 / 3 | 97.0 / none / 66.7 | 1.2 / none / 0.0 |
| 2s/R1 | 41.4 | 49.9 | 78.3 | 65.4 | 80.4 | 986 / 0 / 14 | 78.9 / none / 35.7 | 11.0 / none / 57.1 |
| 2e_n2fit/R1 | 32.3 | 70.9 | 49.5 | 46.0 | 52.3 | 864 / 47 / 89 | 57.3 / 0.0 / 0.0 | 17.6 / 70.2 / 60.7 |
