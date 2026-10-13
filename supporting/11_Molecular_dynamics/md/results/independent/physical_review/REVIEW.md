# Physical review of the MD campaign (hCA II with the sulfanilamide anion on the zinc, 6RL9)

Written 8 Oct 2026, 4:53 AM to 5:45 AM IST, beside the two analysis programs of stage S6. Reviewer's own code on
the raw trajectories and logs. Nothing was written into a run folder, no file under `md/` was changed, the outputs
of `analysis/` and `analysis_recheck/` were not opened, the rented machine was not touched, nothing was committed.

## 0. What was read, and how

- The eleven runs of 100 ns: `bonded_2e/X1 X2 X3 N1 N2 N3 R1 R2 R3`, `bonded_2s/R1`, `bonded_2e_n2fit/R1`. For each:
  every solute frame of the equilibration and of the production (10 ps; 10,200 frames, 10,250 for R1 of set 2e), every
  all-atom frame (100 ps; 1,020 or 1,025), the 20 logs and the records of the steps. The two validation runs X0 were
  read for one comparison only (section 2.3).
- Atoms are chosen by residue number and name from `solute.pdb` and `full.pdb` of each run. `selections.json` was
  compared afterwards: no index differs in any run (`extract.log`, `sel_mismatch []`).
- Fit: CA atoms on `reference.pdb` (Kabsch, no reflection). Images: minimum image with a search of the 27 neighbouring
  cells. Hydrogen bond: donor to acceptor at most 3.5 Å and angle at the hydrogen at least 135° (the limits of the
  specification). States C, R, I from d1 and d2 with the midpoints 6.7 and 5.2 Å. V is the signed volume at N2
  (Zn, H, S), +0.724 in start X, -0.710 in start R.
- Positions of the crystal: start X as built (`crystal_values.json`), waters 423, 516, 553, 568 from
  `md/system/waters_all.pdb`, the rotated ring from `md/system/ligand_R.pdb`. A crystal-water position counts as
  occupied when a water oxygen lies within 1.4 Å of it after the CA fit (limit set before any frame was read; the
  count with 1.0 Å is beside it in `tables.txt`).
- One limit chosen after looking: the class "His64 deep in" (chi1 between 0 and 120° and Zn to NE2 at most 7.0 Å).
  It describes a two-peaked distribution seen in the data and is not a rule.
- CPU: one heavy process at a time, at niceness 10. Four short inspection commands (a few seconds each) were run
  beside the first extraction.

Files, all in this folder: `extract.py` (per-frame quantities, one `data_<model>_<run>.npz` per run), `crystal.py`,
`summ.py` (`tables.txt`, `tables.json`), `events.py` (`events.txt`), `n2state.py` (`n2state.txt`), `compare.py`
(`compare.txt`), `logs.py` (`logs_table.txt`), `rmsf.py` (`rmsf_table.txt`), `dssp_site.py` (`dssp_table.txt`),
`followup.py`, `followup2.py`, `followup3.py` (`followup*.txt`), `extra_numbers.txt`, `common.py`, `run_all.sh`.

## 1. Integrity: nothing wrong was found in any of the eleven runs

Protein (`tables.txt` section 1, `rmsf_table.txt`, `dssp_table.txt`):

| Run | CA RMSD to 6RL9, Å: 0-10 ns, 45-55 ns, 90-100 ns, largest | Radius of gyration, heavy atoms, Å (start 17.378) | Protein to its nearest image, Å, smallest |
|---|---|---|---|
| 2e/X1 | 0.80, 0.89, 0.91, 1.18 | 17.628 ± 0.046 | 21.6 |
| 2e/X2 | 0.85, 0.91, 0.90, 1.22 | 17.615 ± 0.045 | 21.6 |
| 2e/X3 | 0.84, 0.94, 1.02, 1.44 | 17.635 ± 0.066 | 22.0 |
| 2e/N1 | 0.82, 0.88, 0.92, 1.17 | 17.609 ± 0.043 | 22.6 |
| 2e/N2 | 0.88, 0.86, 0.91, 1.21 | 17.620 ± 0.044 | 21.4 |
| 2e/N3 | 0.81, 0.88, 0.96, 1.18 | 17.610 ± 0.047 | 22.3 |
| 2e/R1 | 0.84, 1.06, 1.01, 1.40 | 17.595 ± 0.042 | 22.9 |
| 2e/R2 | 0.83, 0.98, 0.99, 1.44 | 17.582 ± 0.052 | 22.8 |
| 2e/R3 | 0.80, 0.91, 1.00, 1.23 | 17.617 ± 0.042 | 22.1 |
| 2s/R1 | 0.84, 0.90, 0.93, 1.30 | 17.627 ± 0.050 | 23.3 |
| n2fit/R1 | 0.86, 1.01, 0.88, 1.48 | 17.634 ± 0.043 | 21.8 |

- No jump: the largest change of the CA RMSD between two frames 10 ps apart is 0.24 to 0.31 Å. Slope over the last
  50 ns between -0.29 and +0.24 Å per 100 ns. The radius of gyration is 1.2 to 1.5 % above the crystal and steady
  (17.40 to 17.86 Å over all runs).
- Secondary structure on one frame per ns: 93.6 to 95.0 % of the 258 residues keep the assignment of the start;
  strand residues 74.5 to 75.7 against 78, helix residues 32.7 to 35.7 against 40. Every residue at the site (62, 64,
  65, 67, 92, 94, 96, 106, 117, 119, 121, 131, 141, 143, 145, 198 to 202, 207, 209, 244, 246) keeps its assignment
  in 95 to 100 % of the frames. What is lost are three short helices (13 to 15, 155 to 157, 162 to 163) and strand
  ends (87, 170, 233).
- CA fluctuation: mean 0.61 to 0.67 Å; correlation with sqrt(3B/8π²) of 6RL9 0.70 to 0.83 (Pearson), 0.80 to 0.84
  (Spearman). The mean structure is 0.57 to 0.74 Å from the crystal.
- The two cis peptides stay cis: |ω| of Ser29-Pro30 averages 7.0 to 7.3° (largest 39°), of Pro201-Pro202 14.0 to
  15.2° (largest 63°).

Images and bonds (`tables.txt`, `followup.txt`):

- Solute frames: in no frame of any run did the zinc or the ligand have to be moved by a lattice vector (0 of
  112,250 frames); the largest raw Zn-N2 distance is 2.16 Å. All-atom frames: raw Zn to NE2(His94) at most 2.25 Å,
  raw N2 to Zn at most 2.14 Å. The largest distance between neighbouring CA atoms is 4.13 to 4.17 Å (limit of the
  specification 4.5 Å).
- Ligand, all eleven runs: S-N2 1.573 to 1.574 Å (1.462 to 1.676), S-C4 1.821 to 1.822 Å (1.672 to 2.011), S-O 1.460
  to 1.463 Å (1.365 to 1.574), C1-N1 1.355 to 1.356 Å, ring C-C 1.267 to 1.540 Å, N2-H 1.022 Å in every frame (the
  constraint), ring carbons at most 0.114 Å r.m.s. from their plane. The ligand is intact and on the zinc in every
  frame (Zn-N2 at most 2.16 Å).

Logs (`logs_table.txt`):

- 20 segments and 10,000 rows per run; the step is 4.000 fs in every row and every record; the last step is
  25,000,000 at 100,000 ps; every boundary continues by 2,500 steps and 10 ps.
- Temperature 300.96 to 301.00 K (sd 1.46 to 1.48), segment means 300.75 to 301.15 K, slopes within ±0.0012 K per ns.
  Density 1.02530 to 1.02555 g cm⁻³, segment means 1.02495 to 1.02596. Potential energy steady (slopes -0.28 to
  +0.84 kJ mol⁻¹ per ns on means of -542,259 to -543,038 kJ mol⁻¹).
- No jump at a segment boundary: the largest change of the potential energy across a boundary is 2,146 to 2,952
  kJ mol⁻¹, inside the row-to-row scatter within a segment (99th percentile 2,632 to 2,720, largest 3,637 to 4,404).
  The same holds for density and temperature.
- The 1.0 K above the thermostat at the 4 fs step is the known point of the plan (7 Oct, 12:10 AM, point 6); it is
  the same in all eleven runs.

## 2. The zinc site

### 2.1 Set 2e (nine runs), production

| Quantity | Nine runs of set 2e: range of the run means (sd within a run; extremes over all frames) | 6RL9 |
|---|---|---|
| Zn-N2 | 1.925 to 1.928 Å (0.048 to 0.049; 1.707 to 2.161) | 1.931 |
| Zn-NE2 His94 | 2.014 to 2.019 Å (0.058 to 0.059; 1.766 to 2.266) | 2.010 |
| Zn-NE2 His96 | 2.040 to 2.048 Å (0.060 to 0.062; 1.797 to 2.325) | 2.021 |
| Zn-ND1 His119 | 2.035 to 2.039 Å (0.060 to 0.061; 1.742 to 2.297) | 2.018 |
| N2-Zn-His94 | 105.3 to 105.9° | 108.9 |
| N2-Zn-His96 | 115.3 to 116.1° | 114.3 |
| N2-Zn-His119 | 118.2 to 118.8° | 116.8 |
| His94-Zn-His96 | 102.4 to 102.9° | 103.8 |
| His94-Zn-His119 | 114.4 to 114.8° | 113.6 |
| His96-Zn-His119 | 98.2 to 98.8° | 98.4 |
| Zinc out of the ring plane of His94, His96, His119 | -0.25 to -0.30, 0.11 to 0.15, 0.46 to 0.53 Å | -0.17, 0.06, 0.15 |
| Nearer sulfonyl oxygen to Zn | 2.968 to 2.992 Å | 2.985 |

- No excursion: the angles have sd 3.7 to 4.1° and stay between 82.5 and 135.5°; the distributions are single
  peaked. The zinc lies more than 1.0 Å out of the His119 plane in 234 to 316 frames of 10,000 per run.
- Oxygen at the zinc: a sulfonyl oxygen within 2.6 Å in 7 to 17 of 10,000 solute frames per run (0.07 to 0.17 %)
  and in 0 to 3 of 1,000 all-atom frames. No water oxygen and no protein oxygen within 2.6 Å in any all-atom frame
  of any run. A water oxygen within 3.0 Å in 0 to 2 frames of 1,000, except N3 with 9 (smallest 2.73 Å, five
  different waters between 1.7 and 90.3 ns; `followup2.txt` A1). The nearest protein oxygen is Thr199 OG1 (72 to
  84 % of the frames) or Glu106 OE1, never nearer than 3.31 Å.
- Glu106 stays off the zinc: OE1 3.26 Å at the nearest, OE2 3.32 Å. In X2, N2 and n2fit/R1 the carboxylate turns
  over for part of the run (the hydrogen bond from Thr199 goes to OE1 in 45.8, 83.3 and 72.5 % of the frames, to one
  of the two oxygens in 100 %).
- Thr199 OG1-H···Glu106: 99.9 to 100.0 % in every run. The hydroxyl hydrogen never turns to the ligand, so Thr199
  OG1 is the acceptor of N2-H throughout. Arg246 N-H···Glu106: 99.2 to 99.8 %.
- Second shell: His94 ND1-H···Gln92 OE1 98.9 to 99.6 %, His96 ND1-H···Asn244 O 97.9 to 98.6 %, His119 NE2-H···Glu117
  98.5 to 99.3 %.
- The zinc lies 0.28 to 0.33 Å from its crystal position after the CA fit (largest 0.99 Å).

### 2.2 His64 is not settled, and its inner position is not one of 6RL9

Built "out" (conformer A of 6RL9, occupancy 0.55; chi1 -55°, Zn to NE2 11.8 Å). In the runs it moves between out
and in on a scale of tens of ns (`tables.txt`, `followup.txt`, `followup2.txt` A6, `followup3.txt`):

| Run | chi1 gauche+ ("in"), % of the production | "Deep in" (Zn to NE2 at most 7.0 Å), % of the last 80 ns | Water 516 occupied, % of the all-atom frames of 20-100 ns (`compare.json`) |
|---|---|---|---|
| 2e/X1 | 37.2 | 0.7 | 62 |
| 2e/X2 | 8.2 | 6.1 | 74 |
| 2e/X3 | 89.4 | 52.9 | 30 |
| 2e/N1 | 4.2 | 0.0 | 82 |
| 2e/N2 | 0.0 | 0.0 | 78 |
| 2e/N3 | 35.2 | 35.3 | 50 |
| 2e/R1 | 84.2 | 96.8 | 1 |
| 2e/R2 | 84.6 | 87.0 | 7 |
| 2e/R3 | 80.2 | 89.3 | 8 |
| 2s/R1 | 83.5 | 13.3 | 48 |
| n2fit/R1 | 0.0 | 0.0 | 71 |

- In the deep position NE2 lies 6.0 to 6.3 Å from the zinc and on the place of crystal water 516 (1.40 to 1.64 Å
  from it on average in R1 to R3). Over the nine runs of set 2e, water 516 is occupied in 4 % of the deep frames
  (3,279 all-atom frames), in 52 % of the other gauche+ frames and in 76 % of the out frames. With His64 deep in,
  7.5 water oxygens within 8 Å of the zinc become 5.0.
- 6RL9 has no such conformer: its "in" conformers B (0.39) and C (0.05) have the nearer ring nitrogen 7.6 and 8.26 Å
  from the zinc, and water 516 (occupancy 1.0) is hydrogen bonded to ND1 of conformer B at 2.84 Å
  (`md/system/audit/waters.csv`, `histidines.json`). The simulation replaces that water by the ring.
- The aniline nitrogen never comes within 4.0 Å of the His64 ring in any run.

### 2.3 Sets 2s and n2fit (one run each)

- 2s/R1: Zn-N2 1.951 Å, Zn-His96 2.098 Å; N2-Zn-His94 113.6°, N2-Zn-His96 107.6° (crystal 108.9 and 114.3, so the
  two angles are exchanged in size); zinc 0.70 Å out of the His119 plane (more than 1.0 Å in 1,190 frames); nearer
  sulfonyl oxygen 2.855 Å, within 2.6 Å in 254 of 10,000 solute frames (2.54 %) and 26 of 1,000 all-atom frames
  (2.6 %). Its validation run X0 by the same code: 13 of 500 and 3 of 50, as the plan records. Over 100 ns the count
  is below the 5 % of the gate; the mean distance is 0.13 Å shorter than in 6RL9.
- n2fit/R1: Zn-N distances and N(His)-Zn-N(His) angles within 0.01 Å and 0.5° of those of set 2e; N2-Zn-His94 109.0°, N2-Zn-His96 113.5°;
  sulfonyl oxygen within 2.6 Å in 255 solute frames (2.55 %), 18 all-atom frames.

## 3. The ligand

Bound and intact in every run (section 1). Production (`tables.txt` section 3):

| | Nine runs of set 2e | 2s/R1 | n2fit/R1 | 6RL9 or start X |
|---|---|---|---|---|
| Ligand RMSD to the crystal pose, Å | 0.53 to 0.61 (sd 0.21 to 0.27; largest 1.75 to 2.28) | 0.80 ± 0.38 (2.68) | 0.99 ± 0.54 (4.05) | 0 |
| N2-H···OG1 Thr199, % | 99.8 to 99.9 | 99.9 | 98.3 | N···O 2.86 Å |
| Thr199 N-H···O1, % (O2: 0.0 in all) | 99.2 to 99.7 | 97.7 | 95.3 | N···O1 2.94 Å |
| C4-S-N2 | 109.3 to 109.5° | 110.2° | 110.6° | 108.8° |
| τ_fold | 129.6 to 135.0° | 130.7° | 149.1° | 132.1° |
| Zn-N2-S | 113.7 to 114.5° | 109.7° | 112.1° | 117.6° |
| Zn-N2-S-C4 | -121.0 to -123.2° (sd 6.3 to 7.3) | -116.5° (11.4) | -105.7° (16.2) | -117.0° |
| V at N2 | +0.689 to +0.703 (sd 0.074 to 0.079) | +0.715 | +0.489 (sd 0.256) | +0.724 |
| Sum of the three angles at N2 | 333.7 to 334.8° | 332.5° | 344.5° | 332.4° |

- The two sulfonyl oxygens never change places: O1 is the one at Thr199 N-H and O2 the one near the zinc in every
  run, also after a conversion.
- Tilt of the ring in state C (set 2e): the axis C4 to C1 is 6.2 to 7.9° (sd 3.5 to 4.4) from its crystal
  direction, the ring normal 9.5 to 15.1° from the crystal normal, the angle between Zn-N2 and the ring axis 53.8
  to 55.6° (crystal 57.2). In the rotated state before a conversion: 64 to 65° from the crystal axis, 116 to 118°
  between Zn-N2 and the axis (start R 123.8°).
- Aniline nitrogen in state C: no protein N or O within 3.5 Å in 96.9 to 98.2 % of the frames (Pro201 O 0.8 to
  1.8 %, N1 donating; Gln92 NE2 0.4 to 2.1 %). Water: 2.39 to 2.52 water oxygens within 3.5 Å; N1 donates to 1.35 to
  1.41 of them and accepts from 0.23 to 0.28. A chloride within 4.0 Å in 0 to 4 all-atom frames of 1,000; no ion
  nearer than 9.2 Å to the zinc.
- Aniline nitrogen in the rotated state (13 to 35 frames per run before the conversion): Asn67 OD1 in 97 to 100 %
  of the frames (N1 donating in 93 to 100 %), Asn62 ND2 in 83 to 93 % (donating to N1 in 37 to 64 %), no water in
  the all-atom frames of R1 and R3.
- Ring contacts in state C (smallest heavy-atom distance, means): Thr200 3.3 to 3.5, Leu198 3.6, Val121 3.9, Gln92
  3.8 to 4.0, His94 3.9 to 4.0, Phe131 4.4 to 4.7 Å. In the rotated state: Asn67 3.0, Asn62 3.1, His94 3.1, Thr200
  3.4, His96 3.7 to 3.8, Gln92 3.8, Ala65 3.9 Å.
- **Water 553 donates a hydrogen bond to N2 in the crystal-like state.** In 6RL9 it lies 2.93 Å from N2 (also 3.01
  Å from Thr200 OG1). In the runs of set 2e a water is at that position in 87 to 97 % of the all-atom frames, 2.86
  to 2.90 Å from N2, with O-H···N2 at least 135° in 96 to 98 % of those frames (`followup.txt`). N2 in the
  crystal-like state therefore has four partners: the zinc, the sulfur, its hydrogen (to Thr199) and this water.
  The rotated ligand lies on that position (water 553 is one of the four removed for start R) and has no water at
  N2.
- A tilted variant of the crystal-like pose occurs in every run: a water oxygen within 2.5 Å of a ring-carbon
  position of the crystal, the ring turned about S-C4 (τ_fold 154 to 168° against 129 to 135°), ligand RMSD 0.9 to
  1.15 Å against 0.52 to 0.59 Å. It holds 12 to 33 of 1,000 all-atom frames in set 2e, 116 in 2s/R1 and 239 in
  n2fit/R1 (`followup2.txt` A2). It is where the ligand RMSD above 1.0 Å comes from (2.9 to 10.0 % of the frames
  in set 2e; above 2.0 Å at most 0.03 %).

## 4. How the orientation changes

### 4.1 Counts by the rule (own count, from d1 and d2 only)

| Run | Production frames C / R / I, % | Last 20 ns C / R / I, % | Verdict by the rule | Changes of the sign of V from release | First long visit in the opposite state, from release |
|---|---|---|---|---|---|
| 2e/X1 | 99.53 / 0 / 0.47 | 99.8 / 0 / 0.2 | STAYED | 0 | none |
| 2e/X2 | 98.64 / 0 / 1.36 | 99.7 / 0 / 0.3 | STAYED | 0 | none |
| 2e/X3 | 99.50 / 0 / 0.50 | 98.8 / 0 / 1.2 | STAYED | 0 | none |
| 2e/N1 | 99.80 / 0 / 0.20 | 99.5 / 0 / 0.5 | STAYED | 0 | none |
| 2e/N2 | 99.44 / 0 / 0.56 | 99.9 / 0 / 0.1 | STAYED | 0 | none |
| 2e/N3 | 99.71 / 0 / 0.29 | 99.4 / 0 / 0.6 | STAYED | 0 | none |
| 2e/R1 | 99.13 / 0 / 0.87 | 99.7 / 0 / 0.3 | CONVERTED | 1 | 540 ps (stage 5) |
| 2e/R2 | 99.63 / 0 / 0.37 | 99.7 / 0 / 0.3 | CONVERTED | 1 | 370 ps (stage 4) |
| 2e/R3 | 99.51 / 0 / 0.49 | 100.0 / 0 / 0 | CONVERTED | 1 | 160 ps (stage 4) |
| 2s/R1 | 98.76 / 0.01 / 1.23 | 99.5 / 0 / 0.6 | CONVERTED | 5 | 890 ps (stage 6) |
| n2fit/R1 | 87.22 / 3.93 / 8.85 | 74.3 / 9.2 / 16.5 | MIXED | 979 | 60 ps (stage 4) |

The three rotated starts of the primary set converted during the equilibration, the six others stayed; by rule 3
that is the first outcome. The production of the nine runs holds no frame of state R and no frame within 1.64 Å
(ligand RMSD) of the rotated start (`extra_numbers.txt`).

### 4.2 Set 2e, runs R1, R2, R3: one inversion at N2 each, the hydrogen bond to Thr199 kept

Before the release (ligand heavy atoms held), in all three: V -0.55, sum of the angles at N2 344.0 to 344.5°,
Zn-N2-S 123.4 to 124.2°, N2-H···OG1 99 to 100 %, Thr199 N-H···O only 15 to 27 %, ligand 0.38 to 0.40 Å from the
rotated start, Zn-N2 1.94 to 1.96 Å (the 2.235 Å of the docked pose is gone after the minimisation). The nitrogen of
the rotated state is flatter than that of the crystal-like state (sum of angles 344° against 334°).

Sequence, times from the release (`events.txt`; 10 ps between frames):

| | R1 | R2 | R3 |
|---|---|---|---|
| Rotated state after release | 25 frames, 250 ps | 35 frames, 350 ps | 13 frames, 130 ps |
| Last frame before: V; sum of angles | -0.50; 347.6° | -0.45; 350.4° | -0.63; 339.4° |
| First frame with V > 0 | 260 ps: V +0.53, d1 7.35, d2 3.00 (still state R) | 360 ps: V +0.63, d1 6.54, d2 4.16 (state I) | 140 ps: V +0.11, sum 359.5°, d1 7.81, d2 3.09 (state R) |
| First frame in state C | 270 ps | 370 ps | 160 ps |
| Torsion His94 NE2-Zn-N2-S, mean before and mean afterwards | -8° to +78° | -5° to +78° | -5° to +78° |
| Torsion His94 NE2-Zn-N2-H, mean before and mean afterwards | -144° to -155° | -142° to -155° | -147° to -155° |
| N2-H···OG1 in the 20 frames around | every frame | every frame | every frame |
| Thr199 N-H···O1 | forms with the change | forms with the change | forms with the change |
| Water in the place of the crystal ring, last all-atom frame before | 3 (200 ps) | 0, nearest 2.90 Å (300 ps); 2 at 200 ps, 4 at 100 ps | 1 (100 ps) |
| Water at 423, 516, 553, 568 afterwards | all four within 1.3 Å at 1,100 ps | 553 and 568 at 500 ps | 423 at 200 ps, all four at 500 ps |
| Ligand RMSD to the crystal pose, window means | 1.07 to 1.42 Å until 1,000 ps, 0.52 Å from then | 0.42 Å from 400 ps | 0.99 Å at 200 to 300 ps, 0.55 Å from 300 ps |
| Rule time | 540 ps (stage 5) | 370 ps (stage 4) | 160 ps (stage 4) |

Reading:

1. **It is an inversion at N2.** The sign of V changes once per run, between two frames; the sulfur goes from one
   tetrahedral place at N2 to the other (its torsion about Zn-N2 changes by about 85°), the hydrogen stays turned
   to Thr199 (its torsion changes by about 10°, with a swing of 25 to 45° in the frames of the change itself), and
   the hydrogen bond N2-H···OG1 is present in every frame around the change. The ring follows within one or two
   frames (10 to 20 ps). It is not a rotation about Zn-N2 (the hydrogen would leave Thr199) and not one about S-C4
   (which does not move N1).
2. **The ring arrives in the tilted variant first** (section 3): τ(C3-C4-S-N2) 160 to 180°, d2 5.3 to 6.6 Å. In R1 one
   or two waters stay within 2.5 Å of the ring place in the all-atom frames from 400 to 900 ps and the ligand stays
   1.1 to 1.4 Å from the crystal pose until they have left; in R2 and R3 the place is empty and the ring seats within 100 to 200 ps.
3. **All three changes fall in equilibration stage 4**, in which the protein backbone is held to the crystal
   positions with 250 kJ mol⁻¹ nm⁻² (`logs_table.txt`). The primary set therefore holds no frame of the rotated
   pose with a free protein.
4. **The rule time of R1 (540 ps) is 280 ps later than the event.** One frame at 530 ps has d1 7.86 and d2 3.25 Å
   with V +0.30 and a nearly planar nitrogen (sum 356.0°): the ring swung back for less than 30 ps without an
   inversion, and that frame cuts the first crystal-like stretch (270 to 520 ps) below the 1 ns of a long visit.
   By the rule as worded the time is 540 ps; the inversion is at 260 ps.
5. R2: in the last all-atom frame before the change the place of the crystal ring holds no water. The hydration
   check of the plan was met at the end of the held part (2 waters), so the rule was kept; the concern behind it
   (an empty place beside the ligand) applies to that moment of R2. In R1 and R3 water was there and was displaced.
6. After the change the sign of V never changes again (101 ns each).

### 4.3 After the change: the ligand is that of the crystal-start runs; the water at 516 and 423 is not

`compare.txt`: production 20 to 100 ns in eight blocks of 10 ns; the 24 blocks of X1 to X3 give mean ± sd.

| Quantity | X runs (blocks) | R1 | R2 | R3 |
|---|---|---|---|---|
| Ligand RMSD to the crystal pose, Å | 0.55 ± 0.04 | 0.53 | 0.52 | 0.52 |
| d1; d2, Å | 4.71 ± 0.27; 7.76 ± 0.32 | 4.69; 7.64 | 4.65; 7.78 | 4.44; 7.74 |
| V; sum of angles at N2 | 0.70; 334.1 ± 0.3° | 0.70; 334.3 | 0.70; 334.1 | 0.70; 334.1 |
| Zn-N2-S-C4; τ_fold | -122.9 ± 0.7°; 131.3 ± 1.9° | -122.9; 128.9 | -123.3; 129.3 | -123.5; 129.5 |
| N2-H···OG1; Thr199 N-H···O, % | 99.83; 99.54 | 99.91; 99.62 | 99.90; 99.74 | 99.93; 99.74 |
| Ring axis against the crystal axis | 6.6 ± 0.5° | 6.1 | 6.2 | 6.1 |
| Water 553; 568 occupied, % | 94 ± 5; 86 ± 6 | 97; 89 | 98; 93 | 98; 91 |
| Water 423; 516 occupied, % | 29 ± 10; 55 ± 27 | 15; 1 | 12; 7 | 11; 8 |
| Water oxygens in the place of the rotated ring | 3.29 ± 0.32 | 2.58 | 2.58 | 2.61 |
| Water oxygens within 8 Å of the zinc | 7.14 ± 0.93 | 4.78 | 5.32 | 5.06 |
| His64 gauche+, % | 46 ± 47 (run means 10 to 93) | 98 | 100 | 100 |

- Every ligand quantity of R1, R2, R3 lies within 2 block standard deviations of the X runs (largest |z| 1.7). In
  ligand geometry, hydrogen bonds and contacts the converted runs cannot be told from the crystal-start runs.
- The water differs, and the difference is that of His64 (section 2.2): in the first 10 ns the runs from R have the
  hydration of the runs from X (7.1 to 7.8 water oxygens within 8 Å of the zinc, 3.2 to 3.6 in the place of the
  rotated ring; water 516 occupied in 48 to 74 % of the first 20 ns), then His64 moves deep in at 13.6, 16.4 and
  21.8 ns and stays, and water 516 is gone. Within the deep frames the runs from X and N have the same counts as
  the runs from R (5.3 against 4.9 water oxygens within 8 Å; `followup3.txt`).
- Whether three of three runs from R against at most one half in the six others is chance in a slow coordinate or
  a memory of the rotated start cannot be said: Asn62, Asn67, Gln92 and Thr200 have the same side-chain torsions in
  the first 10 ns of all runs (`followup2.txt` B), and no cause was found.

### 4.4 The runs that do not change (X1 to X3, N1 to N3): what holds them

- The nitrogen never flattens: V between +0.16 and +0.92 over 6 × 101.3 ns, below +0.2 in 0 to 1 frame of 10,000;
  the sum of the angles at N2 is above 355° in at most 0.07 % of the frames. No change of sign.
- N2-H···OG1 in 99.8 to 99.9 % of the frames, never absent in two frames running; Thr199 N-H···O1 99.2 to 99.7 %; the
  water at 553 on N2 in 87 to 96 % of the all-atom frames.
- The place where the rotated ring would lie holds 3.0 to 3.5 water oxygens (553 in 87 to 96 %, 568 in 72 to 87 % of
  the frames), and the ring is packed against Leu198, Thr200, Val121, Gln92 and His94.
- Largest excursions: d1 up to 9.29 Å (X2, X3) with d2 never below 4.34 Å; the two conditions of state R are never
  met together. Frames of state I: 21 to 136 of 10,130.

### 4.5 2s/R1: inversion after 30 ps, seating after 0.9 ns

- V changes sign between 20 and 30 ps after the release, again with N2-H···OG1 kept. The ring then stays between the
  two places for 0.9 ns (first 325 ps after release: C 34 %, R 19 %, I 47 %); the first crystal-like visit of 1 ns
  begins at 890 ps, in the free stage 6.
- Of its 18 frames of state R from the release, 16 have V > 0: they are a half-turned ring on a nitrogen of the
  crystal configuration (ring axis 38.5° from the crystal axis, ligand 1.35 Å from the rotated start), not the
  docked pose.
- Two single frames with V < 0 in the production (57.0 and 76.3 ns): in the first the hydrogen, not the sulfur, has
  changed side (its torsion about Zn-N2 -74° against -155°; N2-H···OG1 angle 96°), back within 10 ps.
- Afterwards crystal-like in 98.8 % of the production, but with a ligand RMSD of 0.80 Å, Thr199 N-H···O 97.7 %,
  water 553 in 78 % of the frames and the tilted variant in 11.6 % of the all-atom frames.

### 4.6 n2fit/R1: the nitrogen is nearly planar and the ring goes back and forth

- First change of sign of V 20 ps after the release, first crystal-like frame at 20 ps.
- Production: the sum of the angles at N2 is 344.5 ± 9.4°, above 355° in 17.0 % of the frames; V < 0 in 7.0 %; 490
  stretches with V < 0 from the release, median 1 frame, longest 12 frames (120 ps).
- 125 changes from C to R and 126 back in 101.3 ns; a stay in R lasts 30 ps in the median and 750 ps at the longest;
  no rotated-like visit reaches 1 ns. 111 of the 125 changes have a change of the sign of V within 50 ps.
- State R here is mostly not the docked pose: d1 7.80 ± 0.89 Å, d2 3.85 ± 0.63 Å, ring axis 38° from the crystal axis
  and 28° from that of the rotated start, sum of the angles 355°, V > 0 in 181 of 396 frames. The ligand is within
  1.0 Å of the rotated start in 80 production frames (0.8 %).
- The share of state C falls from 93.3 % in the first half to 81.1 % in the second (blocks of 10 ns: 89.9, 92.1, 98.3,
  92.8, 93.4, 86.2, 78.8, 92.0, 64.5, 84.2 %). Over the ten blocks the share outside C goes with water in the place
  of the crystal ring (correlation +0.82) and against water 553 (-0.73). In the frames of state R and I and in
  the frames with V < 0, water 553 is absent (0 %) and a water lies in the place of the crystal ring in 61 to 70 %.
- The candidate terms distort the crystal-like pose itself: Zn-N2-H 119.6° (as built 106.7 to 107.0°), τ_fold 149°
  (crystal 132.1°), ligand RMSD in state C 0.85 Å.

## 5. The force-field question

What `md/GR.md` and `md/N2.md` established for the terms at N2 as built, on the 64-atom zinc-site model in vacuum
(B3LYP-D3(BJ)/6-31G*; no Thr199, no Thr200, no water):

- (a) planar nitrogen on the constructed inversion path: 13.8 against 6.6 kcal mol⁻¹, the force field too high; at
  the one relaxed geometry near a planar nitrogen 8.0 against 5.6, 2.4 too high;
- (b) the hydrogen alone moved into the plane: 5.2 against 0.5, too high; the free minimum has Zn-N2-H 107.5°
  against 122.8°, so the force field holds the hydrogen too far out of the plane;
- (c) a turn of the ligand by 30° about Zn-N2: 1.8 against 4.2, too free;
- (d) the rotated-like placement against the crystal-like one on rigid geometries: 14.2 against 10.8, the rotated
  one 3.4 too high;
- (e) in both descriptions the rotated-like placement is not a minimum of the bare zinc site.

What the trajectories show, and which error each observation depends on:

1. **The conversion uses the path of error (a), against its direction.** Every conversion is an inversion at N2
   with the sulfur moving. The force field overstates the cost of that path, so error (a) works against a
   conversion; three of three runs converted all the same, within 140 to 360 ps. A lower cost would make it
   faster, not rarer. The observation "the rotated start converts" does not rest on error (a).
2. **Error (d) points the way the result points.** A rotated placement drawn 3.4 kcal mol⁻¹ too high gives the
   force field more reason to leave it and less to return than quantum chemistry has. How much of the 3.4, a number
   from rigid unrelaxed points without the protein, acts in the protein is not known. The one-way outcome of set 2e
   (rotated state for 130 to 350 ps, then never again in about 910 ns of crystal-like sampling) is therefore the
   outcome one expects from this error as well as from the molecule.
3. **Error (c) helps the event a little.** At the frame of the inversion the hydrogen swings by 25 to 45° about
   Zn-N2, a motion the force field makes about 2.4 kcal mol⁻¹ too cheap at 30°. In the settled states the torsions
   about Zn-N2 have sd 6 to 7°, so this error shows only at the event.
4. **Errors (a) and (b) are why set 2e shows no fluctuation at all.** With the terms as built N2 is within |V| < 0.2
   of planar in 0.00 to 0.01 % of the frames and the six crystal-like starts have no change of sign in 608 ns. That
   stiffness is a property of terms known to be too stiff, in the hydrogen by a factor of ten on the test point.
   "The crystal pose never visits the rotated side" is not a statement about the molecule.
5. **The bracket (n2fit/R1) shows the size of that dependence.** With terms that put the planar nitrogen 1.9 to
   2.7 kcal mol⁻¹ too low and remove the 3.4 in sample, the nitrogen is near planar in 17 % of the frames, the sign
   of V changes 979 times and the ring leaves the crystal-like state 125 times. Crystal-like still holds 87.2 % of
   the production, rotated-like 3.9 % (the docked pose itself 0.8 %), and the last 20 ns are MIXED by the rule
   (74.3 %). So: that the crystal orientation is the preferred one holds on both sides of the bracket; that the
   rotated one is absent holds only on the side of the terms as built. The share of 3.9 % is one run of a model that
   failed its gate and drifts over the run; it is not a population of the molecule.
6. **2s/R1 says that the zinc terms do not decide the direction**, and little more: same terms at N2, other bonds
   and angles at the zinc, inversion after 30 ps instead of 140 to 360 ps, crystal-like afterwards. One run, with a
   set that failed gate G-V narrowly and that holds the ligand less well (section 4.5). The earlier inversion goes
   with a flatter nitrogen in its held rotated state (sum of angles 348.3° against 344.0 to 344.5°).
7. **What the gates could not see and the trajectories show.** The 64-atom model has neither Thr199 nor water. In
   the protein the crystal-like nitrogen carries two more partners than in that model: its hydrogen bond to Thr199
   OG1, kept through every inversion, and a hydrogen bond from water 553 (2.93 Å in 6RL9, 87 to 97 % in the runs).
   The rotated placement has no water at N2 and holds Thr199 N-H···O in only 8 to 46 % of its frames, against 99.5 %.
   These terms (fixed charges of -0.914 on N2, TIP3P, ff14SB) favour the crystal-like side and were tested by no
   gate; the crystal structure shows both partners, which speaks for their direction and says nothing about their
   size.

Where the simulation cannot tell:

- **How fast** the two orientations interconvert, and the height of any barrier. The two parameter sets differ by
  orders of magnitude in the number of events (1 against 979 changes of sign in the same time). Three stays of 130,
  250 and 350 ps are three stays; they are not a rate (plan, 6 Oct, 9:15 PM, point 4).
- **By how much** the crystal orientation is preferred. Set 2e gives no number (one-way events); n2fit gives one
  that belongs to a rejected model. No free energy was computed, by rule 5.
- **Whether the rotated pose would last longer in a protein that can adapt.** All three conversions of the primary
  set happened with the backbone held to the crystal; the two sensitivity runs converted even earlier.
- **Whether quantum chemistry would hold the rotated placement in the protein.** Gate G-R says it is no minimum of
  the bare site; with Thr199, Asn62, Asn67 and without water 553 this was not computed.
- **Anything about the neutral molecule.** The docked rotated pose is of neutral sulfanilamide; the runs are of the
  anion with one hydrogen rebuilt on N2.

A wording that the eleven runs support: with the zinc-bound anion, in explicit water, the ring-rotated placement
did not persist in any of the five runs started from it (three with the primary parameters, two sensitivity runs);
it went over to the crystal placement by inversion at the sulfonamide nitrogen, within 0.4 ns, during the
equilibration; the crystal placement was kept in all six runs started from it; the terms at the nitrogen failed
two quantum-chemical tests, one of whose errors (the rotated placement 3.4 kcal mol⁻¹ too high) acts in the
direction of this result, and with refitted terms that err the other way on the planar nitrogen the crystal
orientation still dominates but the ring is no longer fixed.

## 6. What a referee who knows carbonic anhydrase would ask

1. **His64.** Why "out", which tautomer, and does it matter? Built as conformer A (0.55) with the hydrogen on NE2.
   It moves in during eight of the nine runs of set 2e (all but N2), to a position about 2 Å nearer the zinc than either inner
   conformer of 6RL9, and removes water 516 (section 2.2). It does not touch the ligand, but it is a slow
   coordinate that 100 ns do not converge, and the text should not call the water structure of the site reproduced.
   Other tautomers and the protonated form were not run.
2. **The water on the sulfonamide nitrogen** (water 553; sections 3 and 5). It is in the crystal structure and in
   the runs, and it is absent from the quantum-chemical model that judged the nitrogen.
3. **The four ordered waters under the rotated pose.** The docked rotated ring lies on waters 423, 516, 553 and 568
   (B 11.5 to 17.6 Å²). In the production of the runs from X and N, 553 is occupied in 87 to 96 % of the all-atom
   frames, 568 in 72 to 87 %, 516 in 28 to 79 % (it goes with His64) and 423 in 20 to 33 %. The docking that found the two orientations equal had no water there.
4. **The hydrophobic pocket.** With the ring away it holds 1.1 to 2.3 water oxygens on average (0 to 4 in a frame;
   held part of the five runs from R); with the ring in place none in 96.7 to 98.8 % of the all-atom frames of the
   production (set 2e).
5. **The Thr199 gate.** OG1-H to Glu106 in 99.9 to 100 %, OG1 accepting N2-H in 99.8 to 99.9 %, Thr199 N-H to a
   sulfonyl oxygen in 99.2 to 99.7 % (crystal-like state). Glu106 never binds the zinc.
6. **Is the zinc four-coordinate?** Yes in every frame of every run. The nearer sulfonyl oxygen is at 2.97 to 2.99 Å
   (6RL9 2.985) in set 2e, a scaled 1-4 interaction in this model, as the plan notes.
7. **Is the zinc geometry a result?** No for set 2e: its equilibrium values are those of 6RL9. Steady departures
   that no gate tests, confirmed over all nine runs: zinc 0.46 to 0.53 Å out of the His119 plane (0.15), N2-Zn-His94
   105.3 to 105.9° (108.9), Zn-N2-S 113.7 to 114.5° (117.6), S-C4 1.821 to 1.822 Å (1.770), S-N2 1.574 Å (1.596).
8. **Second shell.** Gln92, Asn244 and Glu117 hold their hydrogen bonds to the zinc histidines in 97.9 to 99.6 % of
   the frames.
9. **Independence of the runs.** Three velocity sets per start; 2s/R1 and n2fit/R1 use the seeds of R1. All starts
   share one protein conformation and one water placement.
10. **Sampling.** The ligand pose is the same in nine runs after 1 ns. His64, water 516 and the Glu106 carboxylate
    are not sampled to convergence.
11. **Chloride and sodium.** None within 9.2 Å of the zinc; a chloride at the aniline nitrogen in at most 4 frames
    of 1,000.
12. **Start R.** It is the docked pose of the neutral molecule turned into the anion, with four crystal waters
    deleted. The plan says so; the text has to.

## 7. Findings, by weight

No finding invalidates a run.

Limits a claim:

- L1. All three conversions of the primary set happen in equilibration stage 4, under backbone restraints, 140 to
  360 ps after the release; the production holds no rotated frame (sections 4.1, 4.2).
- L2. The event is an inversion at N2, so it runs through the terms that failed gates G-R and G-R2; one of their
  errors (rotated placement 3.4 kcal mol⁻¹ too high) acts in the direction of the result (section 5).
- L3. With the candidate terms the run is MIXED by the rule and the ring leaves the crystal-like state 125 times;
  "never returns" belongs to the terms as built (sections 4.6, 5).
- L4. Set 2e shows no fluctuation of the nitrogen at all; that is the known stiffness of its terms, not evidence
  about the molecule (section 5, point 4).

Must be reported:

- M1. Water 553 donates to N2 in the crystal-like state; it is not in the quantum-chemical model of the gates
  (section 3).
- M2. His64 moves to an inner position that 6RL9 does not have and displaces water 516; the runs from R spend 87 to
  97 % of their last 80 ns there, the runs from X and N 0 to 53 % (sections 2.2, 4.3).
- M3. After the conversion the ligand is indistinguishable from the crystal-start runs; the water at 516 and 423
  and the number of waters in the site are not (section 4.3).
- M4. The rule time of R1 is 540 ps; the inversion is at 260 ps (section 4.2, point 4).
- M5. R2: the place of the crystal ring was empty in the last all-atom frame before the conversion (section 4.2,
  point 5).
- M6. 2s/R1: sulfonyl oxygen within 2.6 Å of the zinc in 2.54 % of the solute frames, zinc 0.70 Å out of the
  His119 plane, ligand RMSD 0.80 Å (sections 2.3, 4.5).
- M7. State R of the specification does not imply the rotated configuration at N2 in the sensitivity runs (16 of
  18 frames in 2s/R1, 181 of 396 in n2fit/R1 have V > 0) (sections 4.5, 4.6).

Notes:

- N1. Integrity is clean in all eleven runs (section 1).
- N2. A tilted variant of the crystal-like pose with a water under the ring exists in every run (section 3).
- N3. Thr199 N-H to the sulfonyl oxygen is mostly missing in the rotated state (8 to 46 %).
- N4. No defect was found in a file that must not be edited; the two analysis programs and their outputs were not
  read.
