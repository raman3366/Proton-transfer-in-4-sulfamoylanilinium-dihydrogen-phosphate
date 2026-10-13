# Simulation topologies of the hCA II–sulfonamidate complex (stage S3, second part)

Work of 5 Oct 2026, 10:17 PM, to 6 Oct 2026, about 1:05 AM IST, on the Pi with two threads. Files:
`md/build_topology.py`, `md/topology/` (small records), this file; the built systems are in
`/root/cqc-sfa-md/systems/`. Nothing was committed and nothing else in the repository was changed. Every number
below is in a file under `md/topology/`; the file is named in brackets.

## 1. In short

1. **The restrained model is built for the three starts, and every check of the built files passes** (18 checks
   each; section 6). A second build gives the same six files, by sha256 (start X built twice with the final
   script; the five system files of all three starts were unchanged through every rebuild).
2. **The systems hold 40,902 (X), 40,911 (N) and 40,857 (R) atoms, 60 % more than the plan's 25,600.** A
   10 Å buffer of tleap's truncated octahedron leaves 32 Å between the protein and its image as built. The plan's
   number would need a buffer of 4.5 Å, which leaves 12 Å if the protein turns (section 3). The buffer is the
   plan's; GPU time and output size should grow by about the same 60 % (not measured). Someone has to accept
   this or name a smaller buffer (`--buffer-A`).
3. **Restraints: k = 500,000 kJ mol⁻¹ nm⁻² on each of the four Zn–N distances** (thermal spread 0.022 Å). With the
   200,000 first tried, the force field pulled Zn–N2 0.11 Å short (section 5).
4. **The driver minimises all three systems**, the energy falls and the four distances end within 0.07 Å of
   their targets. One check fails and is reported as failed: in start R a sulfonyl oxygen moved 0.62 Å, above the
   0.50 Å fixed beforehand (section 6).
5. **Builder, driver and analysis fit together on the real system**: start N went through the minimisation, the
   six equilibration stages shortened to 4 ps and two production segments of 1 ps, and `analyse_md.py` read the
   run without a warning (section 6).
6. **Warning about the restrained model.** In each of the three 4 ps tests a fifth atom was on the zinc at the
   end: the sulfonyl O2 (X, 2.50 Å), the carboxylate O of Glu106 (R, 1.80 Å), and Glu106 O and water 553
   together (N, 1.90 and 2.09 Å). The four distances hold; the site does not. This is 4 ps of warming up, not a
   gate, but it is what a nonbonded Zn²⁺ with distance restraints only is known to do, and it bears on rule 4 of
   the plan (section 5.4).
7. **The bonded model: the path is written and tested on files that MCPB.py itself wrote**, without the real
   quantum chemistry. The real output has not been built (it does not exist yet). Interface in section 7.
8. **Start R is built with no water in the place the crystal ring takes** (section 4).

## 2. How to run

```
export AMBERHOME=/root/cqc-sfa-tools/md
PY=/root/cqc-sfa-tools/md/bin/python
$PY md/build_topology.py --model restrained --start X        # also N, R; about 1 min each
$PY md/build_topology.py --model bonded --start X            # when md/qm/mcpb/ holds the output of MCPB.py step 4
$PY md/build_topology.py --model restrained --start X --min-check --short-md      # 20 to 25 min
$PY md/build_topology.py --model restrained --start X --k-scan                    # after the line above; 4 min
$PY md/build_topology.py --buffer-scan                                            # 2 min
$PY md/build_topology.py --selftest-bonded --start X [--mcpb-dir DIR]             # 1 min
```

Software: AmberTools 26.0 (tleap; GAFF2 2.2.30), ParmEd 4.3.1, NumPy 2.5.3, SciPy 1.18.1, Python 3.12.14, and for
the checks OpenMM 8.6.1 and MDAnalysis 2.10.0 [`build.json`, `software`]. No random numbers are used.

## 3. What the builder does

1. It writes one PDB file from the four files of `md/system/` (coordinates copied as text): protein (259 residues
   with the cap), zinc, ligand, crystal water, numbered 1 to N. The zinc is residue 260 and the ligand 261, as in
   `md/system/residue_map.csv` and in the input of MCPB.py.
2. **Pass 1 of tleap** (`leaprc.protein.ff14SB`, `leaprc.gaff2`, `leaprc.water.tip3p`; the ligand frcmod as
   assigned, then the repair frcmod; `sfa_anion.mol2`): hydrogens, then `solvateOct sys TIP3PBOX 10.0`.
3. **Ions.** Pairs = round(0.15 × N_water / 55.3), the concentration in the volume of the solvent (55.3 mol L⁻¹ is
   water at 300 K): 33 pairs for 12,314, 12,317 and 12,299 water molecules, and one Na⁺ for the charge −1 of the
   complex: **34 Na⁺ and 33 Cl⁻** in each system. The box volume as built would give 41 pairs, but that box is
   about 10 % too large: its density is 0.91 g cm⁻³ and it will shrink at 1 bar. With 33 pairs the concentration
   is 0.12 mol L⁻¹ in the box as built [`build.json`, `ions`]. This is my reading of "from the box volume"; the
   other number is one constant away.
4. **Where the ions go.** tleap's `addIonsRand` has no seed (its help text says the result differs from call to
   call) and can replace a crystal water; `addIons` puts every ion on the protein surface. The builder therefore
   places them itself, without chance: each ion replaces one water molecule added by `solvateOct`, never a
   crystal water. Candidates are the added molecules at least 6 Å from every atom of protein, zinc, ligand and
   crystal water (about 11,000). They are put in the order of sha256("cqc-sfa-md ion sites|20261005|k") and taken
   in that order if at least 6 Å from every site already taken; Na⁺ and Cl⁻ alternate. Smallest distances found:
   ion to solute 6.35 Å, ion to ion 6.16 Å [`build.json`, `ions.placement`; sites in `ions.pdb`].
5. **Pass 2 of tleap**: the same input, then `remove` for the 67 water molecules and the ions from `ions.pdb`. The
   builder checks that every other atom is exactly where pass 1 put it.
6. The date in the first line of `system.prmtop` is replaced by a fixed stamp, so that a rebuild gives the same
   file. The driver records the sha256 of the topology and will not continue a run on another one.
7. It writes `selections.json`, `reference.pdb`, `extra_forces.json` and `build.json` as in the shared conventions,
   reads the topology back with ParmEd and OpenMM, and runs the checks of section 6.

Some keys of `selections.json` go beyond the conventions; none of the fixed keys is changed:

- `crystal_water_O` (with `crystal_water_numbers_deposited`): the driver then holds the crystal-water oxygens in
  the minimisation, as its protocol foresees. This settles the open point of `md/RUN.md`, section 7, point 3.
- `asn67_ND2`: the analysis reads it if present and can then see a turn of the Asn67 amide.
- `ions_Na`, `ions_Cl`.

### The box against the plan [`buffer_scan_X.json`]

| Buffer of `solvateOct` (Å) | Atoms | Cell edge (Å) | Protein to its image, as built (Å) | The same if the protein turns (Å) |
|---:|---:|---:|---:|---:|
| 4.5 | 26,104 | 73.1 | 22.6 | 12.4 |
| 6.0 | 29,653 | 76.1 | 25.3 | 15.4 |
| 8.0 | 35,057 | 80.1 | 28.8 | 19.4 |
| **10.0 (built)** | **40,902** | **84.1** | **32.4** | **23.4** |
| 12.0 | 47,205 | 88.1 | 36.2 | 27.4 |

"If the protein turns" is the cell edge minus the largest extent of the solute (60.8 Å): the nearest a rigid
protein can come to its image. tleap fits the box to the bounding box of the solute, corners included, which is
why 10 Å gives so much. All values are for the box as built; at 1 bar the edge will be about 3 % shorter (my
estimate from the density of 0.91 g cm⁻³, not measured). The non-bonded cut-off is 10 Å. For the octahedron
tleap accepts one buffer only and always makes the cell isometric (tried here).

## 4. The three systems [`restrained/<start>/build.json`]

| | X | N | R |
|---|---:|---:|---:|
| Atoms | 40,902 | 40,911 | 40,857 |
| Over the plan's 25,600 | +59.8 % | +59.8 % | +59.6 % |
| Solute atoms (protein 4,075, zinc 1, ligand 18) | 4,094 | 4,094 | 4,094 |
| Water molecules (crystal + added) | 12,247 (215 + 12,032) | 12,250 (215 + 12,035) | 12,232 (211 + 12,021) |
| Na⁺ / Cl⁻ | 34 / 33 | 34 / 33 | 34 / 33 |
| Net charge (e) | 0.000000 | 0.000000 | 0.000000 |
| Cell edge (Å), angles 109.47° | 84.135 | 84.134 | 84.121 |
| Volume (nm³) | 458.5 | 458.4 | 458.2 |
| Density as built (g cm⁻³) | 0.912 | 0.913 | 0.912 |
| Protein to its image, as built (Å) | 32.4 | 32.4 | 32.4 |
| tleap errors | 0 | 0 | 0 |
| tleap warnings, final pass | 47 | 47 | 46 |
| Ligand RMSD to the crystal pose (Å); `starts.csv` | 0.0000; 0.0 | 0.6104; 0.61 | 3.5732; 3.573 |
| d1, aniline N to Phe131 CZ (Å); `starts.csv` | 3.8966; 3.897 | 3.3561; 3.356 | 9.5849; 9.585 |
| d2, aniline N to Asn67 OD1 (Å); `starts.csv` | 7.7443; 7.744 | 7.7255; 7.725 | 2.5968; 2.597 |
| `system.prmtop`, sha256 begins | f297a0ed8d6c | a2a489f9e995 | 6446066b509e |
| `system.inpcrd`, sha256 begins | abca2bda29fc | 52bf69d93574 | 7789f44966c8 |

- **tleap warnings.** 45, 45 and 44 pairs of atoms closer than 1.5 Å, every one with a hydrogen (shortest 0.90 Å,
  a lysine hydrogen against a water hydrogen), none with the zinc or the ligand and none between heavy atoms; and
  twice "the charge of the unit (−1) is not zero", before the ions are added. Nothing else. Pass 1 has one more
  charge warning and the note that `savepdb` writes CLYS as LYS. [`build.json`, `tleap`; `leap.log`]
- **Geometry.** Fitted on the 258 CA atoms of `reference.pdb`, every heavy atom of the input files is at its
  deposited or docked coordinates (largest deviation 0.000 Å). Fourteen numbers of `starts.csv` per start are
  recomputed from `system.inpcrd`; the largest difference is 0.0005 Å and 0.04°.
- **Water at the site.** `solvateOct` put no water into the active site: the nearest added water is 5.5 to 6.0 Å
  from the ligand and 10.5 to 11.3 Å from the zinc. In start R no water oxygen lies within 2.5 Å of any position of
  the crystal ligand: the place the crystal ring takes, which the rotated ligand leaves, is empty at the start, as
  `SYSTEM.md` (section 9, point 3) expected. In the 4 ps test the driver counted 1 or 2 water oxygens within 3.5 Å
  of the ligand throughout. Whether water enters during the 600 ps of held equilibration has to be looked at in
  the first real run. [`build.json`, `water_at_the_site`; `min_check.json`, `short_md.stages`]
- **Size on disk**: about 9 MB per system; the folder `build/` beside it (work files of tleap, about 11 MB) need
  not go to the GPU.

### Parameters read back from `system.prmtop` [`build.json`, `force_field`, `ligand_parameters`]

- **Repaired ligand terms are the ones in force**: bond N2–S 461.21, 1.5906 Å; angles C4–S–N2 101.19, 117.133° and
  H–N2–S 61.55, 118.059°; torsion H–N2–S–C4 0.2624 (n = 1, 180°) and 2.7030 (n = 2, 0°), nothing left of the
  terms assigned by analogy (6.8, n = 1; 0.5, n = 3); H–N2–S–O zero. All 18 bonds, 28 angles, 37 torsions and
  7 impropers of the ligand equal those of `md/ligand/sfa_anion_repaired.prmtop`.
- **Ligand charges**: AM1-BCC of the free anion, as in `sfa_anion.mol2`, sum −1.000000.
- **Zinc**: type Zn2+, charge +2, not bonded, 12-6 Lennard-Jones with R_min/2 = 1.271 Å and ε = 0.00330286
  kcal mol⁻¹. This is the line "CM set for TIP3P water from Li et al., JCTC, 2013, 9, 2733" of
  `frcmod.ions234lm_126_tip3p`, the Li/Merz compromise set, which `leaprc.water.tip3p` loads last. No 12-6-4 term.
- **Na⁺, Cl⁻**: Joung–Cheatham for TIP3P (1.369 Å, 0.0874393; 2.513 Å, 0.0355910), as `frcmod.ionsjc_tip3p` has them.
- **Histidines** as in `SYSTEM.md`: His94 and His96 HID, His119 HIE (residue names in `selections_table.txt`).

## 5. The restrained model

### 5.1 What it is

Zinc as a free ion (above), and four harmonic terms U = ½ k (r − r₀)² in `extra_forces.json`, which the driver
adds: zinc to N2 of the ligand, NE2 of His94, NE2 of His96 and ND1 of His119. Nothing acts on an angle, a
torsion or the ring. r₀ are the distances of the 1.0 Å structure, the same for the three starts: **1.9314, 2.0099,
2.0206 and 2.0183 Å**. In N and R the docked Zn–N2 distance (2.207 and 2.235 Å) is therefore pulled to 1.93 Å.

### 5.2 The force constant [`restrained/X/k_scan.json`; `build.json`, `zinc_restraints`]

The spring has to hold against the charge of the ion: zinc (+2) and N2 (−1.25) attract each other with nothing
in between. Offsets of the four distances from r₀ in start X at 0 K, the 502 atoms within 8 Å of the zinc free:

| k (kJ mol⁻¹ nm⁻²) | Zn–N2 (Å) | Zn–N(His), largest (Å) | Thermal spread √(kT/k) at 300 K (Å) | Period of a free N (14.0 amu) against Zn on the spring alone (fs) |
|---:|---:|---:|---:|---:|
| 100,000 | −0.141 | −0.080 | 0.050 | 67.5 |
| 200,000 | −0.111 | −0.055 | 0.035 | 47.7 |
| 300,000 | −0.086 | −0.040 | 0.029 | 39.0 |
| **500,000** | **−0.060** | **−0.027** | **0.022** | **30.2** |
| 700,000 | −0.046 | −0.021 | 0.019 | 25.5 |
| 1,000,000 | −0.033 | −0.015 | 0.016 | 21.3 |

**Chosen: 500,000 kJ mol⁻¹ nm⁻²** (1,195 kcal mol⁻¹ Å⁻² in the same form).

- Strong enough: offset plus two spreads is 0.10 Å for Zn–N2 and 0.07 Å or less for the histidines.
- Not too stiff for 4 fs: with the masses the driver sets, the period is 28.3 fs for Zn–N2 and 30.2 fs for
  Zn–N(His), 7.1 and 7.5 steps. The fastest bond already in the topology that is not constrained, CD1–CG of
  tryptophan, has 21.7 fs, 5.4 steps. The restraint is slower than bonds the 4 fs step carries anyway. At
  700,000 it would be nearly as fast as they are; that value is the alternative if 0.06 Å is too much.
- For scale: the empirical constants MCPB.py gives a *bonded* Zn–N are 73 to 122 kcal mol⁻¹ Å⁻² in Amber's form
  K (r − r₀)², against 598 here [`selftest_bonded_mock/selftest_X.json`]. A restraint on an unscreened ion has to
  be several times stiffer than a bond.

### 5.3 Does it hold? [`restrained/<start>/min_check.json`]

After the minimisation the four distances are 0.013 to 0.067 Å short of r₀ (section 6). In the 2.6 ps after zinc
and ligand are released (26 frames, 273 to 293 K):

| | X | N | R |
|---|---|---|---|
| Mean minus r₀: N2, His94, His96, His119 (Å) | −0.024, −0.011, −0.004, −0.020 | −0.033, −0.003, +0.005, +0.012 | −0.037, +0.015, +0.008, +0.005 |
| Standard deviation (Å); expected 0.022 | 0.031, 0.014, 0.022, 0.021 | 0.024, 0.020, 0.018, 0.025 | 0.019, 0.015, 0.019, 0.025 |
| Largest deviation of any distance in any frame (Å) | 0.081 | 0.098 | 0.079 |

The distances stay within 0.1 Å, and the run is stable at 4 fs for the 650 steps tried.

### 5.4 What it does not hold

Heavy atoms within 2.6 Å of the zinc at the end of the 4 ps (`short_md`):

| Start | On the zinc |
|---|---|
| X | N2 1.87, His94 2.00, His96 2.01, His119 2.01, **sulfonyl O2 2.50 Å** (crystal 2.99 Å; mean over the 26 frames 2.08 Å, all 26 within 2.6 Å) |
| N | N2 1.89, **Glu106 OE1 1.90**, His94 2.00, His119 2.02, His96 2.06, **water 553 2.09 Å**; mean N2–Zn–N(His96) angle 161° (crystal 114°) |
| R | **Glu106 OE1 1.80**, N2 1.88, His119 2.01, His94 2.02, His96 2.04 Å; mean N–Zn–N angles 86 to 124° (crystal 98 to 117°) |

Already in the minimisation, with all heavy atoms held, the N–Zn–N angles move by up to 7.6° (X) and 10.5° (R).
In the crystal Glu106 OE1 is 4.0 Å from the zinc and takes the hydrogen bond of Thr199 (`SYSTEM.md`, section 2). A
free ion of charge +2 draws in the nearest negative oxygen within picoseconds, and four distance restraints
cannot stop it.

This is a short test during warming up, at a density still below 1 g cm⁻³, and no gate was evaluated. But it
says that fallback (b) of the plan, as worded ("the four Zn–N distances only"), will very probably not keep the
site that the 1.0 Å structure shows, with Glu106 on the zinc. A repeat of the rotated start with this model
(rule 4) would then test another active site, not the same one with other zinc terms. The plan's owner should
decide before production. What I can say from here: `MCPB.py -s 2e` (fallback (a), bonded terms from empirical
constants) ran on this machine without quantum chemistry; `md/ZINC.md` describes it.

## 6. Checks

**Of the built files** [`build.json`, `checks`]: 18 checks per start, all passed for X, N and R.
tleap without error and without an unexpected warning; protein chain unbroken; net charge zero; ion numbers as
computed; pass 2 equal to pass 1 except for the ions; input heavy atoms at their coordinates; selections resolve
to the right atoms; geometry agrees with `starts.csv`; zinc type, charge and Lennard-Jones values; zinc not
bonded; Na⁺ and Cl⁻ values; repaired ligand terms in force; ligand terms equal to the free-anion topology; ligand
charges those of the mol2 file; ions keep their distances; no crystal water replaced.

`selections_table.txt` prints residue name, deposited number and atom name for every index (for example
`his119_ND1 1817 HIE 119 ND1`), with the distance from the same atom of the `md/system` files after the fit:
0.0000 Å for the 22 heavy atoms and the ligand hydrogen; the two hydrogens of Thr199, which tleap builds, sit on
OG1 and on N.

**The topology loads** in ParmEd and in OpenMM with the settings of `md/protocol.json`: 40,902 particles and
38,761 constraints for X; the energy of the built coordinates is finite [`openmm_check.json`].

**Minimisation with `md/run_md.py`** [`min_check.json`]. Test protocol: `md/protocol.json` with at most 300
L-BFGS iterations per pass (the production protocol allows 5,000) and output every 0.1 ps; CPU platform, two
threads. Limits were fixed before the first run (`MIN_CHECK` in the script).

| | X | N | R |
|---|---:|---:|---:|
| Potential energy before (kJ mol⁻¹) | −447,644 | −447,990 | −445,147 |
| after | −620,755 | −623,676 | −623,568 |
| Iterations per pass | 300, 300, 300, 183 | 300, 300, 300, 300 | 300, 300, 300, 300, 105 |
| Wall time (s) | 572 | 489 | 554 |
| Zn–N2 before → after (Å); r₀ 1.931 | 1.931 → 1.864 | 2.207 → 1.917 | 2.235 → 1.881 |
| Zn–N(His94, His96, His119) after (Å); r₀ 2.010, 2.021, 2.018 | 1.990, 2.008, 1.997 | 1.966, 1.969, 1.979 | 1.981, 1.980, 1.993 |
| Largest distance minus r₀ (Å); limit 0.10 | −0.067 | −0.052 | −0.051 |
| Zinc moved (Å) | 0.116 | 0.122 | 0.208 |
| Ligand heavy atom that moved most (Å) | N2 0.153 | N1 0.435 | O1 0.625 |
| Limit for zinc and ligand heavy atoms (Å) | 0.30 | 0.50 | 0.50 |
| Ligand hydrogen that moved most (Å) | 0.345 | 1.068 | 1.137 |
| Ligand RMSD to the crystal pose, before → after (Å) | 0 → 0.083 | 0.610 → 0.336 | 3.573 → 3.403 |
| d1 before → after (Å) | 3.897 → 3.971 | 3.356 → 3.633 | 9.585 → 9.567 |
| d2 before → after (Å) | 7.744 → 7.460 | 7.725 → 7.570 | 2.597 → 2.846 |
| S–N2 before → after (Å) | 1.596 → 1.574 | 1.678 → 1.527 | 1.677 → 1.564 |
| Held atoms (2,289; 2,285 in R): r.m.s. and largest move (Å) | 0.240, 1.431 | 0.247, 1.368 | 0.245, 1.383 |
| Energy lower | passed | passed | passed |
| Four distances within 0.10 Å of r₀ | passed | passed | passed |
| Zinc and ligand heavy atoms within the limit | passed | passed | **failed** |

- **The failed check (R).** O1 moved 0.62 Å and O2 0.55 Å; the limit of 0.50 Å was my estimate before the run. The
  cause is the docked geometry relaxing: Zn–N2 goes from 2.24 to 1.88 Å, S–N2 from 1.68 to 1.56 Å, and the SO₂
  group turns with it (Zn–N2–S–C4 from −32° to −51°). The ligand stays in the rotated orientation (d1 9.57 Å,
  d2 2.85 Å). I left the limit as it was.
- The held atoms that moved most are side-chain ends and crystal waters (Asp41 OD2, Glu26 OE1, Lys39 NZ, waters
  472 and 565), up to 1.4 Å against a restraint of 1,000 kJ mol⁻¹ nm⁻².
- **The production minimisation was not run to its end.** With 5,000 iterations per pass it needs more than an
  hour per system on the shared Pi; I stopped two attempts. `--min-check --min-full` runs it.

**Six equilibration stages shortened to 4 ps** (`--equil-scale 0.002`; 0.2, 0.6, 0.6 ps at 2 fs, then 0.8, 0.8,
1.0 ps at 4 fs): all ran for X, N and R. Temperature at the end 292.6, 292.7, 292.7 K; density from 0.912 to
0.976, 0.990, 0.980 g cm⁻³; 0.4 to 0.9 ns per day on the Pi. A smoke test, not an equilibration.

**End to end, start N**: two production segments of 1 ps followed, and `md/analysis/analyse_md.py` read the run
folder: exit code 0, 20 solute frames and 2 full frames, 0 warnings (its guards compare atom names at every
mapped index and the residue names of the CA atoms with `reference.pdb`) [`restrained/N/min_check.json`,
`end_to_end`]. The numbers of 2 ps mean nothing and are not quoted.

## 7. Interface for the bonded model

**Note of 6 Oct 2026, 11:40 AM IST** (after the first build from the real MCPB.py output; account in
`md/test/REHEARSAL_2e.md`, change in `md/test/rehearsal_2e/build_topology_6oct.diff`). The text below describes the
builder as it was on 6 Oct, 1:00 AM. What differs now:

- The bonded model is built for a named parameter set: `--model bonded --set 2e` (or `2s`) `--start X`; the folders
  are `bonded_2e` and `bonded_2s`, not `bonded`. The MCPB.py folder is `--mcpb-dir`, by default the first of
  `/root/cqc-sfa-md/qm/gr/mcpb_<set>` and `md/qm/mcpb` that holds `CA2_tleap.in`; the builder stops if the frcmod
  there was not made by the step the set names.
- The terms at N2 (Zn–N2–S, Zn–N2–H, S–N2–H, H–N2–S–C4) come from one file loaded after every file of MCPB.py:
  `--n2-frcmod FILE`, or for set 2e without it a file the builder writes with the terms as built before gate G-R
  (Zn–N2–H at 107.41°, not the 120° of step 2e). It is kept as `md/topology/bonded_<set>/<start>/n2_terms.frcmod`.
- **The six "one sided connection" warnings are no longer expected and no longer occur.** The renamed histidines
  now get head, tail and connect atoms, and their six peptide `bond` lines are not used: with the `bond` lines
  tleap had built the amide hydrogen of Phe95, Trp97 and Leu120 before those bonds existed, 0.4 to 1.2 Å out of
  place. The two mock tests below were made with that fault.
- More is read back and checked (zinc terms against the parameter files, charges against the templates, the
  pairs across the zinc, every solute atom against the restrained build): `build.json`, keys `zinc_site_bonded`,
  `terms_at_N2`, `against_the_restrained_system_of_this_start`.

`--model bonded` reads the output of `MCPB.py -s 4` from `md/qm/mcpb/` (or `--mcpb-dir`), which
`md/qm/after_qm.sh` will leave there:

| File | What the builder takes |
|---|---|
| `CA2_tleap.in` | the `addAtomTypes` block (M1, Y1 to Y4); the lines `X = loadmol2 X.mol2`; the `loadamberparams` lines in their order; the ten `bond` lines (four to the zinc, six peptide bonds); the name of the PDB file |
| `CA2_mcpbpy.pdb` | residue names only: residues 93, 95, 118, 260 and 261 become HD1, HD2, HE1, ZN1 and SN1 in the builder's own PDB file. Its coordinates are **not** used, so the same parameters serve X, N and R. Residue k there must be residue k here (CA atoms within 0.5 Å), or the build stops |
| `HD1.mol2`, `HD2.mol2`, `HE1.mol2`, `ZN1.mol2`, `SN1.mol2` | templates of the five residues of the site with their charges and atom types |
| `CA2_mcpbpy.frcmod`, `sfa_anion_gaff2.frcmod`, `sfa_anion_repair.frcmod` | parameters, loaded from that folder; the builder notes whether the two ligand files equal those of `md/ligand/` |

Not taken from MCPB.py: its `source` lines (the builder sources `leaprc.protein.ff14SB`, `leaprc.gaff2` and
`leaprc.water.tip3p`; MCPB.py writes `oldff/leaprc.ff14SB` and no GAFF2 line), `solvatebox`, `addions`, and
**`frcmod.ions1lm_126_tip3p`**, which would replace the Joung–Cheatham Na⁺ and Cl⁻ by the Li/Merz ones. The two
models should differ in the zinc site only.

Everything after that is the same code as for the restrained model: solvation, ions, selections, reference,
checks. `extra_forces.json` then holds an empty list. Read back and checked: the zinc has exactly the four
bonds; the chain is unbroken around the renamed residues; the repaired ligand terms are in force under the new
type of N2 (reported, not required, because gate G-R may refit them); bonds, angles and torsions with zinc are
listed in `build.json` (`zinc_site_bonded`), with the charge of the site. The six "one sided connection"
warnings of tleap are expected: the renamed residues have no connect atoms, and the peptide bonds come from the
`bond` lines.

**Tested** [`selftest_bonded_mock/`, `selftest_bonded_synthetic/`]:

1. On a mock made here with MCPB.py itself: its steps 2e and 4 on a copy of the step-1 files of `md/qm/mcpb/`,
   and mol2 files with the unfitted charges of ff14SB and of the free anion. Built for X and R: all checks pass;
   zinc bonded to NE2 of His94 and His96, ND1 of His119 and N2; 6 angles centred on the zinc; 39 torsions through
   it, all zero; the six repaired terms in force under type Y4.
2. On the folder that `md/qm/test_after_qm.py` makes, where MCPB.py steps 2, 3 and 4 ran on synthetic
   quantum-chemical files: the mol2 files are those MCPB.py step 3 writes, the frcmod in the folder is that of
   step 2e. Built for X and R: all checks pass, net charge 0.000000.

Both results are marked `MOCK_NOT_A_MODEL` in `/root/cqc-sfa-md/systems/_selftest_bonded_*/`; they must never be
simulated. Without `CA2_tleap.in` the builder stops with exit code 2 and says which file is missing.

## 8. Not tested, and open

1. **The bonded model with the real parameters**: not built, not minimised, gate G-Z not evaluated. The first
   build should be followed by `--min-check --short-md`.
2. **Decision, box**: 40,900 atoms against 25,600 (section 3). `md/RUN.md` estimated output sizes for 25,600.
3. **Decision, sensitivity model**: the restrained model lets Glu106 or a sulfonyl oxygen onto the zinc within
   4 ps (section 5.4).
4. **Decision, ions**: 33 pairs by the solvent volume; 41 by the box volume as built.
5. **Empty pocket in start R** (section 4): to be watched in the first real equilibration.
6. **The minimisation with the production limit**, any dynamics beyond 4 ps, the density at 1 bar, and
   everything on a GPU: not run here.
7. **Two-thread check runs do not repeat to the last digit** (`md/RUN.md`, test 0); the numbers of sections 5.3
   and 6 would differ slightly in a second run. The built systems do repeat.
8. **The second self-test reads another agent's throwaway folder** (`…/scratchpad/mdprep/zinc/after_qm_test/mcpb`).
   If it is gone, `md/qm/test_after_qm.py` makes it again; the checksums of the files used are in
   `selftest_bonded_synthetic/bonded/<start>/build.json`.
9. `md/ZINC.md` appeared while I worked; I read its parts on file names and `md/qm/after_qm.sh`, and changed
   nothing under `md/qm/`.

## 9. Files

`md/build_topology.py`: the builder and its checks.

`md/topology/`:

| Path | Content |
|---|---|
| `restrained/<start>/build.json` | the record of the build (copy of the one beside the system): inputs and outputs with sha256, software, composition, box, ions with their sites, parameters read back, warnings, geometry, checks |
| `restrained/<start>/leap_stage1.in`, `leap.in`, `leap.log`, `ions.pdb` | the two tleap inputs, the log of the final pass, the ion positions |
| `restrained/<start>/selections_table.txt` | what every index of `selections.json` points at |
| `restrained/<start>/SHA256SUMS` | the six files of the built system |
| `restrained/<start>/openmm_check.json` | the topology in OpenMM |
| `restrained/<start>/min_check.json` | minimisation, 4 ps test and, for N, the end-to-end test |
| `restrained/X/k_scan.json` | offsets against the force constant |
| `buffer_scan_X.json` | box size against the buffer |
| `selftest_bonded_mock/`, `selftest_bonded_synthetic/` | the two tests of the bonded path |

`/root/cqc-sfa-md/systems/` (not in git): `restrained/X`, `N`, `R` with `system.prmtop`, `system.inpcrd`,
`selections.json`, `extra_forces.json`, `reference.pdb`, `build.json` and `build/`; `_checks/` (the run folders
of the tests, 170 MB); `_selftest_bonded_mock/` and `_selftest_bonded_synthetic/` (mocks).
