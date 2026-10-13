# Production: parameter sets, built systems, validation runs and the queue

Written on 6 Oct 2026 (IST) by the integration task of the evening; state of 6 Oct 2026, 10:45 PM IST.
The contract is the plan's entry of 6 Oct 2026, 9:15 PM IST (`md/PLAN.md`): production with the terms at N2 as
built, by the user's decision. Every number below is read from a file named beside it; this file is assembled
from those files by a script, not typed. Times are IST (UTC+5:30). What was added on 7 Oct 2026 is marked where
it stands and listed at the end of the file.

## 1. Outcome in short

| | |
|---|---|
| Primary zinc set | **2e** (empirical constants, MCPB.py step 2e), systems `bonded_2e` |
| Why not 2s | Set 2s (Seminario) passed gate G-Z, but its validation run X0 **failed gate G-V**: part d of GV3, an oxygen within 2.6 Å of the zinc in 3 of 50 all-atom frames (6 %; at most 5 % allowed). Every other number of the gate is met (section 5) |
| Set 2e | passed gate G-Z and gate G-V with every number (sections 4 and 5; what parts a to c of GV3 can show for this set: section 5.3, paragraph added on 7 Oct 2026) |
| The two analysis programs | agree on every value of both runs: 42,850 values compared for X0 of 2s and 42,845 for X0 of 2e, 0 disagreements, 0 threshold frames, exit status 0 both times |
| Queued | nine runs X1 N1 R1, X2 N2 R2, X3 N3 R3 with `bonded_2e` (100 ns each), then two sensitivity runs of start R: `R_zn` (zinc set 2s) and `R_n2` (set 2e with the candidate terms at N2 that gate G-R2 did not accept) |
| Terms at N2 | as built, the same file in both zinc sets (sha256 `113e3a2ebdbe99a4389e2ed497b20caaf0e64aebcbc8d547cf39823e9ee2cd00`); gates G-R and G-R2 on these terms failed (`md/GR.md`, `md/N2.md`) |

The plan's rule (entry of 9:15 PM, point 1): "the Seminario set (`bonded_2s`) is the primary one, provided it
passes gate G-Z and its validation run X0 passes gate G-V; if it fails either, that is reported and the empirical
set (`bonded_2e`) … is primary, again only after its own X0." That rule was applied as written. The failing part
of set 2s is narrow (three frames where two are allowed; the same event is seen in
13 of 500 of the ten times denser solute frames,
2.6 %); it is reported with all its numbers in section 5.3 and the rule
was not read again after the result.

One thing was done ahead of the rule: **X0 of both sets was put into the queue at once** (X0 of 2s first, X0 of
2e behind it), before the result of 2s was known. `md/RUN.md` gives every simulated zinc model its own
validation run, set 2e is simulated in either case (as the primary set or as the sensitivity set), and the two
runs take 12 minutes each. Nothing of production was queued before both gates were read.

## 2. Parameter sets

### 2.1 Where they come from

| Set | MCPB.py step | Work folder (Pi) | `CA2_mcpbpy.frcmod` sha256 | Made |
|---|---|---|---|---|
| 2e | 2e, empirical zinc constants | `/root/cqc-sfa-md/qm/gr/mcpb_2e` | `38221d1d52d8b0af047145619d29bdeee11fdaf3010c452503facd01ec8f29c9` | 6 Oct 2026, 7:10 AM IST (date of the folder), unchanged since |
| 2s | 2s, Seminario's method on the Hessian of the optimised small model | `/root/cqc-sfa-md/qm/gr/mcpb_2s` | `a10b145de977a167e5d56071aa721c0e454dc154edb4f9d7c6c443fb54dd8518` | 6 Oct 2026, 9:18 PM IST, this task |

Set 2s was made with `bash md/qm/after_qm.sh 2s /root/cqc-sfa-md/qm/gr/mcpb_2s <report folder>` in a copy of the
step-1 files of `md/qm/mcpb/` with `CA2_small_opt.fchk` (sha256 `cd7b5ae7ff949e4353dead159f5f6475815b8bcfcb5a03b2697f56cbbf97fc8d`, the
file the machine wrote at 8:47 PM IST). `md/qm/after_qm.sh` (sha256 `ee4df1172773207e…`) and the
guard `md/qm/mcpb_step2s.py` (sha256 `332290ab477652e1…`) were used as they are; **neither
needed a change**. The guard's line: "36 values passed to MCPB.py's averaging; 18 were of complex type, largest
|imaginary|/|value| 0.0e+00 (real part taken)". The four hard checks of the script passed
(`/root/cqc-sfa-md/qm/gr/after_qm_report_2s/after_qm_check_2s.json`):

- H1 ligand terms with the ligating nitrogen arrive under the new type: passed
- H2 tleap: no error, 4,094 atoms, net charge -1: passed
- H3 four Zn-N bonds, 14 angles with zinc, torsions through zinc all zero: passed
- H4 charges of the site: five residues sum to +1, ligand charges refitted: passed

The two sets differ in the four zinc bonds and the fourteen angles with the zinc in `CA2_mcpbpy.frcmod`, and in
nothing else: the residue templates (`HD1.mol2`, `HD2.mol2`, `HE1.mol2`, `ZN1.mol2`, `SN1.mol2`, with the
equal-charge RESP set), `CA2_tleap.in` and `CA2_mcpbpy.pdb` of the two work folders are identical byte for byte
(sha256 in `md/topology/bonded_2e/X/build.json` and `md/topology/bonded_2s/X/build.json`, key `mcpb.files_sha256`).
Step 2ms (modified Seminario) is tabulated beside them and not used (plan, 6 Oct 2026, 3:02 AM, point 3).

### 2.2 The zinc terms of both sets (and 2ms beside them)

Amber form: bonds U = k (r − r0)², k in kcal mol⁻¹ Å⁻²; angles U = k (θ − θ0)², k in kcal mol⁻¹ rad⁻².
The columns 2e and 2s are read back from the built topologies of start X and equal the parameter files
(`build.json`, key `zinc_site_bonded`, check `zinc_bonds_and_angles_as_the_parameter_files`); the column 2ms is
from the frcmod that step 2ms writes (`md/production/seminario_blocks.json`). "In 6RL9" is the value in the
crystal start X.

**Bonds to the zinc**

| Bond | Types | In 6RL9 (Å) | 2e: k | 2e: r0 (Å) | 2s: k | 2s: r0 (Å) | 2ms: k | 2ms: r0 (Å) |
|---|---|---|---|---|---|---|---|---|
| Zn-N2 | M1-Y4 | 1.9314 | 122.1 | 1.9310 | 102.4 | 1.9665 | 102.4 | 1.9665 |
| Zn-His94 NE2 | Y1-M1 | 2.0099 | 78.0 | 2.0100 | 73.6 | 2.0258 | 73.6 | 2.0258 |
| Zn-His96 NE2 | Y2-M1 | 2.0206 | 73.3 | 2.0210 | 51.2 | 2.0840 | 51.2 | 2.0840 |
| Zn-His119 ND1 | Y3-M1 | 2.0183 | 74.3 | 2.0180 | 65.4 | 2.0423 | 65.4 | 2.0423 |

**Angles with the zinc** (2ms has the equilibrium values of 2s)

| Angle | Types | In 6RL9 (°) | 2e: k | 2e: θ0 (°) | 2s: k | 2s: θ0 (°) | 2s θ0 − 6RL9 (°) | 2ms: k |
|---|---|---|---|---|---|---|---|---|
| His94 NE2-Zn-His96 NE2 | Y1-M1-Y2 | 103.82 | 35.00 | 103.82 | 40.22 | 101.65 | -2.17 | 34.35 |
| His94 NE2-Zn-His119 ND1 | Y1-M1-Y3 | 113.59 | 35.00 | 113.59 | 36.78 | 113.97 | +0.38 | 26.95 |
| His94 NE2-Zn-N2 | Y1-M1-Y4 | 108.92 | 35.00 | 108.92 | 47.10 | 118.19 | +9.27 | 34.67 |
| His96 NE2-Zn-His119 ND1 | Y2-M1-Y3 | 98.41 | 35.00 | 98.41 | 29.79 | 102.43 | +4.02 | 25.21 |
| His96 NE2-Zn-N2 | Y2-M1-Y4 | 114.28 | 35.00 | 114.28 | 51.42 | 100.64 | -13.64 | 42.66 |
| His119 ND1-Zn-N2 | Y3-M1-Y4 | 116.83 | 35.00 | 116.83 | 30.27 | 116.09 | -0.74 | 22.24 |
| His94 CD2-His94 NE2-Zn | CV-Y1-M1 | 124.21 | 50.00 | 124.21 | 54.00 | 128.47 | +4.26 | 27.02 |
| His94 CE1-His94 NE2-Zn | CR-Y1-M1 | 127.89 | 50.00 | 127.89 | 51.01 | 124.95 | -2.94 | 25.53 |
| His96 CD2-His96 NE2-Zn | CV-Y2-M1 | 122.91 | 50.00 | 122.91 | 65.13 | 132.81 | +9.90 | 32.62 |
| His96 CE1-His96 NE2-Zn | CR-Y2-M1 | 129.31 | 50.00 | 129.31 | 64.74 | 120.86 | -8.45 | 32.43 |
| His119 CE1-His119 ND1-Zn | CR-Y3-M1 | 121.05 | 50.00 | 121.05 | 54.65 | 121.77 | +0.72 | 27.77 |
| His119 CG-His119 ND1-Zn | CC-Y3-M1 | 131.22 | 50.00 | 131.22 | 51.88 | 130.31 | -0.91 | 26.52 |

**The two angle lines at N2 that MCPB.py writes are not in force** in any system: the file of the terms at N2 is
loaded last and replaces them (section 2.5). For the record:

| Types | 2e: k, θ0 | 2s: k, θ0 | 2ms: k, θ0 |
|---|---|---|---|
| M1-Y4-hn | 50.00, 120.00 | 43.30, 122.76 | 31.59, 122.76 |
| M1-Y4-sy | 50.00, 117.55 | 74.18, 107.74 | 52.91, 107.74 |

Torsions through the zinc: 39 in each set, every one with zero force constant;
no improper torsion holds the zinc (`torsions_through_zinc` in both `build.json`).

### 2.3 What step 2s read from the Hessian

The Hessian is a partial one: the six atoms that were displaced are
HID93@NE2 (12), HID95@NE2 (27), HIE118@ND1 (38), ZN260@ZN (46), SAN261@N2 (54), SAN261@S (57) of the
64-atom small model (MCPB.py's numbering: residues 93, 95 and 118 are His94, His96 and His119 of 6RL9, 260 is the
zinc and 261 the ligand), so 732 of the 4096 3 × 3 blocks are computed
and the others are zero. `md/production/seminario_blocks.py` runs MCPB.py's own step 2s (and 2ms) in a
temporary copy with the reading of the Hessian instrumented and records every block that is read
(`md/production/seminario_blocks.json`):

- step 2s reads 4 bonds and 14 angles; for a bond a–b the block
  [a, b], for an angle a–b–c the blocks [a, b] and [c, b] (the averaging over both directions is off);
- **every block read for every bond and every angle is a computed one: true**
  (each has a displaced atom on at least one side; none is a zero block); the same holds for step 2ms
  (true);
- the constants that the instrumented run gives equal the frcmod of the work folder
  (true), and a second calculation with code of
  its own (same blocks) gives the same numbers to the printed digits.

So nothing blocks set 2s on the side of the Hessian.

### 2.4 Are the constants of 2s physical?

In size, yes: bonds 51.2 to 102.4 kcal mol⁻¹ Å⁻², angles at the zinc
29.79 to 51.42 and Zn–N–C
51.01 to 65.13 kcal mol⁻¹ rad⁻²; no negative
or near-zero constant, no complex eigenvalue with an imaginary part, and MCPB.py's rules of thumb hold (bond
constants below 200 and lengths below 2.8 Å; angle constants below 100 and angles above 100°;
`after_qm_check_2s.json`, key `mcpb_rules_of_thumb`). Three things have to be said plainly:

1. **The equilibrium values of 2s are those of the small model optimised in vacuum, not of the crystal.** Two of
   them are far from 6RL9: N2–Zn–NE2 His96 has θ0 = 100.64° against 114.28° in the crystal (−13.64°), and
   N2–Zn–NE2 His94 118.19° against 108.92° (+9.27°); the two Zn–N–C angles of His96 are 132.81° and 120.86°
   against 122.91° and 129.31°. The gate asks for means within 10° of 6RL9, so set 2s was at risk at that
   angle before any run; in X0 the mean came out at 107.96° (-6.34°, inside the window).
2. **The bond Zn–NE2 His96 is soft and long in 2s**: k = 51.2 with r0 = 2.0840 Å (2e: 73.3, 2.0210 Å; crystal
   2.0206 Å). In X0 its mean is 2.097 Å, inside the window that ends at 2.12 Å.
3. **Seminario's method is not unique for Zn–N2**: with the block [Zn, N2], which MCPB.py reads, k = 102.44;
   with the block [N2, Zn] 99.53 (mean 100.99). For the three histidine bonds the two blocks agree within 0.3.

What followed in the run is in section 5.3: with set 2s the sulfonyl oxygen O2 is nearer the zinc than with
set 2e.

### 2.5 The terms at N2

Both zinc sets are built with **the same file of the terms at N2, as built**
(`md/topology/bonded_2e/<start>/n2_terms.frcmod`, sha256 `113e3a2ebdbe99a4389e2ed497b20caaf0e64aebcbc8d547cf39823e9ee2cd00`), loaded after every
file of MCPB.py. Read back from the topology of each of the six systems (`build.json`, key
`terms_at_N2.in_the_topology.entries_of_the_file`; check `terms_at_N2_are_those_of_the_n2_file`):

| Term | Types | k | Equilibrium or phase |
|---|---|---|---|
| angle Zn–N2–S | M1-Y4-sy | 50.00 kcal mol⁻¹ rad⁻² | 117.55° |
| angle Zn–N2–H | M1-Y4-hn | 50.00 kcal mol⁻¹ rad⁻² | 107.407° (its value in start X, not the 120° of step 2e) |
| angle S–N2–H | hn-Y4-sy | 61.550 kcal mol⁻¹ rad⁻² | 118.059° |
| torsion H–N2–S–C4 | ca-sy-Y4-hn | 0.2624 kcal mol⁻¹, n = 1 | 180° |
| | | 2.703 kcal mol⁻¹, n = 2 | 0° |

These terms did not pass gate G-R (`md/GR.md`) and the refit did not pass gate G-R2 (`md/N2.md`); they are used
by the user's decision of 6 Oct 2026 (plan, entry of 9:15 PM).

**The candidate of the sensitivity run** (`md/qm/n2/n2_terms_new_charges_as_built.frcmod`, sha256
`2244a06c6e760ca079ffdd525eac9314dba9c3989e7b116d23e5545101cb7f29`; its first line says "NOT ACCEPTED"): read back from the topology of
`bonded_2e_n2fit/R`, 12 of 12 entries of the file are in force as written:

| Entry of the file | In force (k and equilibrium; for torsions k, n, phase) |
|---|---|
| ANGLE M1-Y4-sy | 20.423 kcal mol⁻¹ rad⁻², 113.229° |
| ANGLE M1-Y4-hn | 23.041 kcal mol⁻¹ rad⁻², 118.486° |
| ANGLE hn-Y4-sy | 51.074 kcal mol⁻¹ rad⁻², 117.127° |
| DIHE ca-sy-Y4-hn | 2.124 kcal mol⁻¹, n = 1, 180°; 2.8826 kcal mol⁻¹, n = 2, 0°; 1.2233 kcal mol⁻¹, n = 3, 0° |
| DIHE M1-Y4-sy-ca | 3.1094 kcal mol⁻¹, n = 1, 180°; 0.1183 kcal mol⁻¹, n = 2, 0°; 0.5797 kcal mol⁻¹, n = 3, 0° |
| DIHE Y1-M1-Y4-sy | 0.3739 kcal mol⁻¹, n = 3, 292.439° |
| DIHE Y2-M1-Y4-sy | 0.3739 kcal mol⁻¹, n = 3, 292.439° |
| DIHE Y3-M1-Y4-sy | 0.3739 kcal mol⁻¹, n = 3, 292.439° |
| DIHE Y1-M1-Y4-hn | 0.4154 kcal mol⁻¹, n = 3, 8.243° |
| DIHE Y2-M1-Y4-hn | 0.4154 kcal mol⁻¹, n = 3, 8.243° |
| DIHE Y3-M1-Y4-hn | 0.4154 kcal mol⁻¹, n = 3, 8.243° |
| IMPROPER M1-hn-Y4-sy | 0.5301 kcal mol⁻¹, n = 2, 180° (N2 the central atom) |

### 2.6 The charges of the site

The equal-charge RESP set (`md/qm/resp_equiv.py`; the six pairs of chemically equal ligand atoms made equal), the
same in every system of this file: the residue templates of the two work folders are identical, and the charges
of the five site residues in each topology equal the templates (`build.json`, key `zinc_site_bonded.charges`;
checks `site_charges_are_those_of_the_templates`, `site_charge_plus_one`, `ligand_equal_atoms_have_equal_charges`).
Record of the fit: `/root/cqc-sfa-md/qm/gr/after_qm_report_2s/resp_equiv_2s.json`, sha256 `5392f6c6bef6f28f7f667d502ae80676c86bce23db31e14800c6f91739788c70`.

| Residue | Sum of charges (e) | Largest difference from the unconstrained fit (e) |
|---|---|---|
| HD1 (number 94 of 6RL9) | +0.334448 | 0.023369 |
| HD2 (number 96 of 6RL9) | +0.413171 | 0.023943 |
| HE1 (number 119 of 6RL9) | +0.447207 | 0.074304 |
| ZN1 (number 301 of 6RL9) | +0.335600 | 0.043853 |
| SN1 (number 304 of 6RL9) | -0.530426 | 0.049095 |
| all five | +1.0000 | |

| Atom | Charge (e) |
|---|---|
| Zn | +0.335600 |
| NE2 of His94 | -0.001770 |
| NE2 of His96 | +0.047689 |
| ND1 of His119 | +0.122659 |
| C1 of the ligand | +0.464944 |
| C2 of the ligand | -0.274695 |
| C3 of the ligand | -0.153047 |
| C4 of the ligand | +0.080205 |
| C5 of the ligand | -0.153047 |
| C6 of the ligand | -0.274695 |
| N1 of the ligand | -0.954927 |
| N2 of the ligand | -0.914224 |
| O1 of the ligand | -0.488550 |
| O2 of the ligand | -0.488550 |
| S of the ligand | +0.790275 |
| H2 of the ligand | +0.146108 |
| H3 of the ligand | +0.155431 |
| H5 of the ligand | +0.155431 |
| H6 of the ligand | +0.146108 |
| HN11 of the ligand | +0.417016 |
| HN12 of the ligand | +0.417016 |
| HN21 of the ligand | +0.398775 |

## 3. The built systems

Builder: `md/build_topology.py`, sha256 `a9b8ef8d699dfce68abfddcd21186ebdb939900d67b4c549b92692bf39d38d70` (the same file built all eight systems:
`script_sha256` in each `build.json`). Built on the Pi on 6 Oct 2026 between
9:27 PM and 9:31 PM IST with

```bash
cd /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate
export AMBERHOME=/root/cqc-sfa-tools/md OMP_NUM_THREADS=2 OPENMM_CPU_THREADS=2
PY="nice -n 10 /root/cqc-sfa-tools/md/bin/python"
for s in X N R; do
  $PY md/build_topology.py --model bonded --set 2e --start $s --note "…"
  $PY md/build_topology.py --model bonded --set 2s --start $s --n2-frcmod md/topology/bonded_2e/$s/n2_terms.frcmod \
      --compare-with bonded_2e --note "…"
done
for b in 2e 2s; do
  $PY md/build_topology.py --model bonded --set $b --variant n2fit --start R \
      --n2-frcmod md/qm/n2/n2_terms_new_charges_as_built.frcmod --compare-with bonded_$b --note "…"
done
```

(the notes are in `STATUS_NOTE.txt` of each system). Systems: `/root/cqc-sfa-md/systems/<system>/<start>/`
(`system.prmtop`, `system.inpcrd`, `selections.json`, `reference.pdb`, `extra_forces.json` with an empty list,
`build.json`, `STATUS_NOTE.txt`); records: `md/topology/<system>/<start>/` (`build.json`, `leap.in`,
`n2_terms.frcmod`, `SHA256SUMS`, `openmm_check.json`, `compare_with_<other>.json`).

| System | Start | Atoms | Water | Na+ / Cl− | Checks true / n.a. / false | system.prmtop sha256 | system.inpcrd sha256 | builder sha256 |
|---|---|---|---|---|---|---|---|---|
| bonded_2e | X | 40,902 | 12,247 | 34 / 33 | 26 / 7 / 0 | `a6c57a7c20aaea633317173b691837de56301feeed6f5cf38d687e91c9085aad` | `27acb79e664525a7f2c7d7f245d7e8bb4092fb0b2a2f67b1d1ee6c6cd0220d19` | `a9b8ef8d699dfce6…` |
| bonded_2e | N | 40,911 | 12,250 | 34 / 33 | 26 / 7 / 0 | `98e46802e34ae08c96b5ca4cb187a061bed3c89202ffeeafaa512be43a7e4f16` | `9ce16b4c588d2235467e12ebbf6591fec5faa5e1b599437ff5bc8b948f532725` | `a9b8ef8d699dfce6…` |
| bonded_2e | R | 40,857 | 12,232 | 34 / 33 | 26 / 7 / 0 | `7ff894cfda73a2bb8bb68cc3b1198328b7f78a413fb15768de0b82f02395d732` | `fe5bc95b9f90444fc2dcd9aa1ba619e991dea0ebe500716b2ced28ee1ec4e360` | `a9b8ef8d699dfce6…` |
| bonded_2s | X | 40,902 | 12,247 | 34 / 33 | 27 / 6 / 0 | `5c1e5543950f05d4db0a707f6e739657ff506c7b1166047e50af3fdc9674b319` | `27acb79e664525a7f2c7d7f245d7e8bb4092fb0b2a2f67b1d1ee6c6cd0220d19` | `a9b8ef8d699dfce6…` |
| bonded_2s | N | 40,911 | 12,250 | 34 / 33 | 27 / 6 / 0 | `06f1e9820d156b66d51c572a15e43e3bbaa6654bba0350d182dab9568edabc56` | `9ce16b4c588d2235467e12ebbf6591fec5faa5e1b599437ff5bc8b948f532725` | `a9b8ef8d699dfce6…` |
| bonded_2s | R | 40,857 | 12,232 | 34 / 33 | 27 / 6 / 0 | `7cd9f3d1da06d2703f5a567801ba0cc3c6ac4802370e3d68adbe24e28bcc7296` | `fe5bc95b9f90444fc2dcd9aa1ba619e991dea0ebe500716b2ced28ee1ec4e360` | `a9b8ef8d699dfce6…` |
| bonded_2e_n2fit | R | 40,857 | 12,232 | 34 / 33 | 26 / 7 / 0 | `d4999d63b48ef7c43faebb5b108b45bd87d6b9cc7cc6d6f69c22fddc296d8673` | `fe5bc95b9f90444fc2dcd9aa1ba619e991dea0ebe500716b2ced28ee1ec4e360` | `a9b8ef8d699dfce6…` |
| bonded_2s_n2fit | R | 40,857 | 12,232 | 34 / 33 | 26 / 7 / 0 | `b4116282903bd85fa1abbf8730f97e5c710a6d366b5818bf20fb8a88b3989bd5` | `fe5bc95b9f90444fc2dcd9aa1ba619e991dea0ebe500716b2ced28ee1ec4e360` | `a9b8ef8d699dfce6…` |

No check failed in any system (`checks_failed` is empty in all eight `build.json`). The checks that do not apply
to a bonded system are counted as n.a. The three systems of `bonded_2e` have the checksums of the five system
files that the rehearsal ran with (`md/test/REHEARSAL_2e.md`): only `build.json` and `STATUS_NOTE.txt` are new
(the note and the new checks).

**Crystal water** (added on 7 Oct 2026 after the independent check of the set-up; plan, entry of 7 Oct 2026,
12:10 AM, point 2). The systems hold 215 of the 216 water sites of chain A of 6RL9: site 510 (occupancy 0.32) is
dropped in all starts. Start R holds 211, because four ordered waters (423, 516, 553 and 568; all fully occupied,
B 11.5 to 17.6 Å² against a median of 22.5 Å² for all water) lie on the rotated ligand,
0.28 Å (N1), 1.43 Å (C6), 0.55 Å (C5) and 0.19 Å (C2) from its heavy atoms, and were removed when the start was built, by
the rule that a water whose oxygen is within 2.2 Å of a heavy atom of the ligand is removed (`md/SYSTEM.md`,
section 4; `md/system/decisions.json`, keys `water` and `starts.R.waters_removed`; `composition.crystal_water_molecules`
in the eight `build.json`: 215 in X and N, 211 in R). So the runs from X and N start with these four waters in the place
that the rotated ring would take, and the runs from R start with the ligand there and without them. The docking
test `dock/WATER4.md` is built on the same four waters.

**Each start has the same coordinates, box, ions and selections in every set**: `system.inpcrd`,
`selections.json` and `reference.pdb` have the same sha256 as in `bonded_2e` (table, and
`compare_with_*.json`, keys `system_inpcrd_same_sha256`, `selections_json_same_sha256`,
`reference_pdb_same_sha256`).

**What differs between the sets, by read-back of the two topologies** (`compare_with_<other>.json` in each
record folder; every bond, angle, torsion and improper term compared atom by atom, with the atom names, types,
charges, masses, Lennard-Jones parameters, the excluded-atom lists and the 1–4 pairs with their scaling):

| System | Compared with | Terms that differ | Any other difference | 1–4 pairs | Sections of the prmtop that differ as text |
|---|---|---|---|---|---|
| bonded_2s/X | bonded_2e | 6 × angle N-Zn-N, 6 × angle Zn-N-C, 4 × bond Zn-N | none (atoms 0, classes not allowed 0) | 10,878 = 10,878, identical with their scaling: true | ANGLE_EQUIL_VALUE, ANGLE_FORCE_CONSTANT, BOND_EQUIL_VALUE, BOND_FORCE_CONSTANT |
| bonded_2s/N | bonded_2e | 6 × angle N-Zn-N, 6 × angle Zn-N-C, 4 × bond Zn-N | none (atoms 0, classes not allowed 0) | 10,878 = 10,878, identical with their scaling: true | ANGLE_EQUIL_VALUE, ANGLE_FORCE_CONSTANT, BOND_EQUIL_VALUE, BOND_FORCE_CONSTANT |
| bonded_2s/R | bonded_2e | 6 × angle N-Zn-N, 6 × angle Zn-N-C, 4 × bond Zn-N | none (atoms 0, classes not allowed 0) | 10,878 = 10,878, identical with their scaling: true | ANGLE_EQUIL_VALUE, ANGLE_FORCE_CONSTANT, BOND_EQUIL_VALUE, BOND_FORCE_CONSTANT |
| bonded_2e_n2fit/R | bonded_2e | 3 × angle at N2, 1 × improper at N2, 2 × torsion about N2-S, 6 × torsion about Zn-N2 | none (atoms 0, classes not allowed 0) | 10,878 = 10,878, identical with their scaling: true | ANGLE_EQUIL_VALUE, ANGLE_FORCE_CONSTANT, DIHEDRALS_INC_HYDROGEN, DIHEDRALS_WITHOUT_HYDROGEN, DIHEDRAL_FORCE_CONSTANT, DIHEDRAL_PERIODICITY, DIHEDRAL_PHASE, POINTERS, SCEE_SCALE_FACTOR, SCNB_SCALE_FACTOR |
| bonded_2s_n2fit/R | bonded_2s | 3 × angle at N2, 1 × improper at N2, 2 × torsion about N2-S, 6 × torsion about Zn-N2 | none (atoms 0, classes not allowed 0) | 10,878 = 10,878, identical with their scaling: true | ANGLE_EQUIL_VALUE, ANGLE_FORCE_CONSTANT, DIHEDRALS_INC_HYDROGEN, DIHEDRALS_WITHOUT_HYDROGEN, DIHEDRAL_FORCE_CONSTANT, DIHEDRAL_PERIODICITY, DIHEDRAL_PHASE, POINTERS, SCEE_SCALE_FACTOR, SCNB_SCALE_FACTOR |

So `bonded_2s` differs from `bonded_2e` in the 4 zinc bonds and the 12 zinc angles of section 2.2 only, and each
`n2fit` system differs from its base set only in the terms at N2: the 3 angles, the 8 torsion quadruples (2 about
N2–S: H–N2–S–C4 and Zn–N2–S–C4; 6 about Zn–N2: each histidine nitrogen with S and with H) and the 1 improper at
N2; the 1–4 pairs are identical. (In the `n2fit` systems the pointers and the 1–4 scale-factor tables differ as
text because the new torsion types add entries to them; the pairs themselves and their scaling are compared one
by one and are identical.)

The builder's checks were made **set-specific, not weaker**: `torsions_through_zinc_all_zero` holds for the sets
as built and is replaced, only in a variant, by `torsions_through_zinc_zero_or_as_the_n2_file` (a torsion through
the zinc is zero unless the file of the terms at N2 carries it, and then it is the file's);
`n2_file_is_of_the_class_of_the_set` refuses a file that holds anything but the terms of its class;
`terms_at_N2_are_those_of_the_n2_file` reads every entry of the file back from the topology;
`differs_from_the_compared_set_only_where_the_sets_differ` is the comparison above.

## 4. Gate G-Z

The plan (stage S2; entries of 5 Oct 2026, 9:41 PM, and 6 Oct 2026, 3:02 AM): "after minimisation and 100 ps",
(a) Zn–N2 within 0.10 Å of 1.93 Å, (b) each Zn–N(His) within 0.10 Å of 2.01–2.02 Å, (c) each mean N–Zn–N angle
within 10° of its value in 6RL9, (d) an oxygen (sulfonyl, water or protein) within 2.6 Å of the zinc in at most
5 % of the frames; Zn···Glu106 OE1 reported by name. These are parts a to d of GV3 of `md/analysis/SPEC.md`,
section 11, and their numbers are taken as printed there.

**Which frames.** The gate is read from the validation run X0 of each set, on **the first 100 ps after the
minimisation**: step `eq01` (heating from 50 K to 300 K, NVT), frames 0 to 9 of the equilibration series of the
analysis (10 solute frames at 10, 20, … 100 ps after the minimisation; −1.99 to −1.90 ns on the time axis of the
series) and frame 0 of the all-atom series (at 100 ps). Program: `md/production/gate_GZ.py` (sha256
`680c974ea76ac5ee…`), which reads `timeseries_equil_solute.csv` and `timeseries_equil_full.csv`
of the analysis; records `md/production/gate_GZ_bonded_2s.json` and `md/production/gate_GZ_bonded_2e.json`. The two
analysis programs give the same numbers on these frames (144 values compared in each run, largest
difference 0; the two series files are identical byte for byte).

| Part | Quantity | Must lie in | 2s: minimised | 2s: mean of the 10 frames | 2s: frame at 100 ps | 2s | 2e: minimised | 2e: mean of the 10 frames | 2e: frame at 100 ps | 2e |
|---|---|---|---|---|---|---|---|---|---|---|
| a | Zn–N2 (Å) | 1.83 to 2.03 | 1.949 | 1.911 | 1.898 | met | 1.927 | 1.904 | 1.951 | met |
| b | Zn–NE2 His94 (Å) | 1.91 to 2.12 | 1.993 | 1.998 | 2.053 | met | 2.008 | 2.011 | 2.039 | met |
| b | Zn–NE2 His96 (Å) | 1.91 to 2.12 | 2.105 | 2.109 | 2.110 | met | 2.034 | 2.018 | 1.970 | met |
| b | Zn–ND1 His119 (Å) | 1.91 to 2.12 | 2.056 | 2.048 | 1.999 | met | 2.042 | 2.024 | 2.005 | met |
| c | N2–Zn–NE2 His94 (°) | 98.9 to 118.9 | 111.845 | 113.080 | 114.053 | met | 104.331 | 104.119 | 106.049 | met |
| c | N2–Zn–NE2 His96 (°) | 104.3 to 124.3 | 109.855 | 110.451 | 110.191 | met | 115.722 | 117.204 | 114.036 | met |
| c | N2–Zn–ND1 His119 (°) | 106.8 to 126.8 | 117.737 | 116.390 | 110.970 | met | 118.769 | 118.422 | 115.926 | met |
| c | NE2 His94–Zn–NE2 His96 (°) | 93.8 to 113.8 | 102.797 | 103.436 | 108.361 | met | 102.677 | 101.979 | 103.616 | met |
| c | NE2 His94–Zn–ND1 His119 (°) | 103.6 to 123.6 | 115.094 | 114.933 | 116.148 | met | 114.822 | 115.073 | 117.658 | met |
| c | NE2 His96–Zn–ND1 His119 (°) | 88.4 to 108.4 | 97.181 | 96.046 | 95.433 | met | 99.753 | 99.236 | 99.002 | met |

| | 2s | 2e |
|---|---|---|
| d: sulfonyl O within 2.6 Å of the zinc, solute frames | 0 of 10 (nearest 2.642 Å) | 0 of 10 (nearest 2.810 Å) |
| d: any O within 2.6 Å, all-atom frame | 0 of 1 (nearest: sulfonyl 2.895, water 3.517, protein 3.906 Å) | 0 of 1 (nearest: sulfonyl 2.973, water 3.648, protein 3.937 Å) |
| d | met | met |
| Zn···Glu106 OE1 (reported, no criterion): mean, smallest | 4.049, 3.874 Å; within 2.6 Å in 0 of 10 | 3.928, 3.780 Å; within 2.6 Å in 0 of 10 |
| **Gate G-Z** | **PASSED** | **PASSED** |

The columns "minimised" (the structure after the minimisation, from `state_min.xml`) and "frame at 100 ps" are
for information; the verdict is taken from the mean of the 10 frames and from part d.

**What this gate can and cannot show.** In `eq01` the heavy atoms of the solute, the zinc and the ligand
included, are held to the crystal positions (U = ½ k d², k = 1000 kJ mol⁻¹ nm⁻² = 2.39 kcal mol⁻¹ Å⁻²). That is
far weaker than a zinc bond term, but not negligible beside an angle term at the zinc, so the angles of these
frames lean towards the crystal. The same four parts for every stage of the equilibration, for information and
without a verdict (stages 4 to 6 hold no atom of the zinc site):

| Stage | Held atoms, k | 2s: Zn–N2 | 2s: Zn–N His94/96/119 | 2s: six angles | 2s: O within 2.6 Å, solute / all-atom | 2s a–d | 2e: Zn–N2 | 2e: Zn–N His94/96/119 | 2e: six angles | 2e: O within 2.6 Å | 2e a–d |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 (heat, 100 ps) | solute_heavy, 1000 | 1.911 | 1.998 / 2.109 / 2.048 | 113.1 / 110.5 / 116.4 / 103.4 / 114.9 / 96.0 | 0 of 10 / 0 of 1 | met | 1.904 | 2.011 / 2.018 / 2.024 | 104.1 / 117.2 / 118.4 / 102.0 / 115.1 / 99.2 | 0 of 10 / 0 of 1 | met |
| 2 (density, 300 ps) | solute_heavy, 1000 | 1.926 | 1.991 / 2.103 / 2.023 | 112.6 / 110.1 / 117.2 / 102.3 / 115.8 / 95.7 | 2 of 30 / 0 of 3 | **NOT MET** (d: fifth ligand) | 1.912 | 1.998 / 2.052 / 2.002 | 105.0 / 115.2 / 118.7 / 102.8 / 115.2 / 98.8 | 0 of 30 / 0 of 3 | met |
| 3 (loosen, 300 ps) | solute_heavy, 250 | 1.962 | 2.006 / 2.078 / 2.041 | 111.0 / 110.2 / 115.8 / 103.0 / 117.2 / 97.3 | 1 of 30 / 0 of 3 | met | 1.935 | 2.023 / 2.020 / 2.040 | 104.0 / 116.6 / 116.9 / 104.1 / 115.3 / 99.3 | 0 of 30 / 0 of 3 | met |
| 4 (release_side_chains_zinc_ligand, 400 ps) | protein_backbone, 250 | 1.949 | 2.007 / 2.106 / 2.078 | 112.7 / 107.6 / 119.0 / 100.0 / 115.6 / 98.0 | 0 of 40 / 0 of 4 | met | 1.933 | 2.020 / 2.039 / 2.039 | 104.1 / 116.4 / 118.5 / 101.9 / 115.6 / 99.1 | 0 of 40 / 0 of 4 | met |
| 5 (loosen_backbone, 400 ps) | protein_backbone, 50 | 1.951 | 2.008 / 2.086 / 2.063 | 112.1 / 108.5 / 119.7 / 102.1 / 114.8 / 96.0 | 1 of 40 / 0 of 4 | met | 1.917 | 2.025 / 2.041 / 2.037 | 105.3 / 117.1 / 118.5 / 102.7 / 114.3 / 97.7 | 0 of 40 / 0 of 4 | met |
| 6 (free, 500 ps) | none, 0 | 1.969 | 2.029 / 2.106 / 2.057 | 112.2 / 108.4 / 118.5 / 102.6 / 114.9 / 96.9 | 1 of 50 / 0 of 5 | met | 1.927 | 2.016 / 2.039 / 2.040 | 104.8 / 115.9 / 118.4 / 103.0 / 114.9 / 98.7 | 0 of 50 / 0 of 5 | met |

The six angles are in the order of the table above. For set 2s one stage (the second, 30 frames) has the
sulfonyl oxygen within 2.6 Å in 2 of 30 solute frames, more than 5 %; it is not one of the frames of the gate,
and it is the same behaviour that fails part d of gate G-V in the 5 ns (section 5.3).

## 5. Validation run X0 and gate G-V

### 5.1 The runs

| | X0 of 2s | X0 of 2e |
|---|---|---|
| Folder (machine and, pulled, Pi) | `/root/cqc-sfa-md/runs/bonded_2s/X0` | `/root/cqc-sfa-md/runs/bonded_2e/X0` |
| System, `system.prmtop` sha256 | `5c1e5543950f05d4db0a707f6e739657ff506c7b1166047e50af3fdc9674b319` | `a6c57a7c20aaea633317173b691837de56301feeed6f5cf38d687e91c9085aad` |
| Driver `md/run_md.py` sha256 | `6d797f080ff75679…` | `6d797f080ff75679…` |
| Run seed (from the run id X0) | 1793718340 | 1793718340 |
| Started, ended | 6 Oct 2026, 9:37:14 PM IST; 6 Oct 2026, 9:49:27 PM IST | 6 Oct 2026, 9:49:30 PM IST; 6 Oct 2026, 10:01:50 PM IST |
| Minimisation | 21.64 s, L-BFGS iterations per pass [906, 3312, 6, 7] | 26.2 s, [879, 4362, 6, 6] |
| Equilibration, 6 stages, 2 ns | 238.6 s; no extension step | 239.0 s; no extension step |
| Hydration check at the end of stage 3 | skipped, the ligand itself is at the site (0.282 Å from it on average) | skipped (0.228 Å) |
| Production, 5 ns in one segment | 464.4 s, **930.15 ns per day**; 500 solute and 50 all-atom frames | 467.0 s, **925.03 ns per day**; 500 and 50 |
| Checksums of the pulled files against the driver's records | ok (watchdog log, "driver checksums … 0 bad") | ok |

Both runs have the same seeds (the seed follows from the run id), so they start with the same velocities and
differ in the zinc terms alone. Full protocol of `md/protocol.json` (sha256 `decda42f5340d06f…`), on the L40S, CUDA, mixed precision.

### 5.2 The analysis and the comparison of the two programs

```bash
nice -n 10 /root/cqc-sfa-tools/md/bin/python md/analysis/recheck_md.py --compare /root/cqc-sfa-md/runs/bonded_2s/X0 --run-primary
nice -n 10 /root/cqc-sfa-tools/md/bin/python md/analysis/recheck_md.py --compare /root/cqc-sfa-md/runs/bonded_2e/X0 --run-primary
```

run on the Pi after the watchdog had pulled each run; `analyse_md.py` sha256 `f75e6155de700700…`,
`recheck_md.py` sha256 `5b4c0595c067ff57…`, specification 1.3. Both programs analysed each run in this call.

| | X0 of 2s | X0 of 2e |
|---|---|---|
| Values compared | 42,850 | 42,845 |
| Disagreements beyond tolerance | 0 runs; 0 values | 0 runs; 0 values |
| Threshold frames; verdicts that differ through them | 0; 0 | 0; 0 |
| Warnings (first / second program) | 0 / 0 | 0 / 0 |
| Exit status | 0 | 0 |
| Gate verdict, first / second program | FAILED / FAILED | PASSED / PASSED |
| Report | `runs/bonded_2s/X0/analysis_recheck/recheck_report.md` | `runs/bonded_2e/X0/analysis_recheck/recheck_report.md` |

**The two programs agree on every value of both runs.**

### 5.3 Gate G-V, every number

Gate G-V as `md/analysis/SPEC.md`, section 11, version 1.3, words it (the plan's entries of 5 Oct 2026, 9:41 PM);
means over the 500 solute frames of the 5 ns as written to `summary.json` (key `gate_GV`); the verdict is the
program's own (`gate_GV.verdict`; role `validation`, not provisional). In brackets: the deviation from the
reference.

| Criterion | Quantity | Must be | 2s (X0) | 2s | 2e (X0) | 2e |
|---|---|---|---|---|---|---|
| GV1 | mean heavy-atom RMSD of the ligand from the crystal pose (Å) | ≤ 1.0 | 0.712205 | met | 0.550160 | met |
| GV2 | frames with the hydrogen bond N2–H···OG1 Thr199 | ≥ 80 % | 500 of 500 (100.0 %) | met | 499 of 500 (99.8 %) | met |
| GV3 a | mean Zn–N2 (Å) | 1.83 to 2.03 | 1.952221 | met | 1.924440 | met |
| GV3 b | mean Zn–NE2 His94 (Å) | 1.91 to 2.12 | 2.017209 | met | 2.016962 | met |
| GV3 b | mean Zn–NE2 His96 (Å) | 1.91 to 2.12 | 2.096659 | met | 2.047101 | met |
| GV3 b | mean Zn–ND1 His119 (Å) | 1.91 to 2.12 | 2.054899 | met | 2.034849 | met |
| GV3 c | mean N2–Zn–NE2 His94 (°); 6RL9 108.9 | 98.9 to 118.9 | 113.239917 (+4.339917) | met | 105.617480 (-3.282520) | met |
| GV3 c | mean N2–Zn–NE2 His96 (°); 6RL9 114.3 | 104.3 to 124.3 | 107.957375 (-6.342625) | met | 115.783846 (+1.483846) | met |
| GV3 c | mean N2–Zn–ND1 His119 (°); 6RL9 116.8 | 106.8 to 126.8 | 118.493403 (+1.693403) | met | 118.314403 (+1.514403) | met |
| GV3 c | mean NE2 His94–Zn–NE2 His96 (°); 6RL9 103.8 | 93.8 to 113.8 | 101.700238 (-2.099762) | met | 102.413849 (-1.386151) | met |
| GV3 c | mean NE2 His94–Zn–ND1 His119 (°); 6RL9 113.6 | 103.6 to 123.6 | 115.046929 (+1.446929) | met | 114.710536 (+1.110536) | met |
| GV3 c | mean NE2 His96–Zn–ND1 His119 (°); 6RL9 98.4 | 88.4 to 108.4 | 96.633566 (-1.766434) | met | 98.718929 (+0.318929) | met |
| GV3 d | solute frames with a sulfonyl O within 2.6 Å of the zinc | ≤ 5 % | 13 of 500 (2.6 %) | met | 1 of 500 (0.2 %) | met |
| GV3 d | all-atom frames with any O (sulfonyl, water, protein) within 2.6 Å of the zinc | ≤ 5 %, at least one frame | 3 of 50 (6.0 %) | **NOT MET** | 0 of 50 (0.0 %) | met |
| GV3 | all parts | | | **NOT MET** | | met |
| GV4 a | mean C4–S–N2 (°), s.d. | 105.8 to 111.8 | 109.968094, 2.532340 | met | 109.562265, 2.506503 | met |
| GV4 b | mean τ_fold (°), s.d. | 112.1 to 152.1 | 126.777612, 14.613415 | met | 131.375949, 13.634488 | met |
| **Gate G-V** | `gate_GV.verdict` | | | **FAILED** | | **PASSED** |

**What parts a to c of GV3 show for set 2e** (added on 7 Oct 2026 after the independent check of the set-up;
plan, entry of 7 Oct 2026, 12:10 AM, point 3). For set 2e the equilibrium values of the four Zn–N bonds and the
six N–Zn–N angles are those of 6RL9 itself (MCPB.py step 2e; section 2.2; `md/ZINC.md`, section 8). Parts a to c of GV3 therefore show only
that the rest of the force field does not move the site away from these values (largest shifts of the means of
X0 from the equilibrium values: 0.026 Å, Zn–NE2 His96, and 3.30°, N2–Zn–NE2 His94; the brackets of the table are taken
against the references of the gate, which are rounded to 0.1°); they are not an independent reproduction of the zinc geometry.
The equilibrium value of the angle Zn–N2–S in the file of the terms at N2 is the value of 6RL9 as well (117.55°, section 2.5);
its mean in X0, which is no number of the gate, lies 3.56° below it.
The parts of the gate that are not put in are GV1, GV2, GV3 d and GV4.

Reported beside the gate, no criterion: Zn···Glu106 OE1 mean 4.412 Å, smallest 3.557 Å, within
2.6 Å in 0 of 500 frames (2s); mean 4.279 Å, smallest 3.405 Å, within 2.6 Å in 0 of 500 (2e).
Zn···Glu106 OE2: mean 6.225 and 6.201 Å, smallest 5.168 and 5.093 Å.
Orientation of both runs: STAYED (2s: 498 of 500 frames in the crystal state, 0 rotated) and STAYED (2e: 499 of 500, 0 rotated).

**The failure of set 2s, in full.** Part d of GV3 has two counts, each of which may be at most 5 %: the solute
frames (every 10 ps; the sulfonyl oxygens) and the all-atom frames (every 100 ps; any oxygen).

| | 2s | 2e |
|---|---|---|
| Solute frames with a sulfonyl O within 2.6 Å of the zinc | 13 of 500 (2.6 %): met | 1 of 500 (0.2 %): met |
| All-atom frames with any O within 2.6 Å | **3 of 50 (6 %): not met** (2 of 50 would be the most allowed) | 0 of 50: met |
| The all-atom frames concerned: time, distance of the sulfonyl O | 1.6 ns, 2.582 Å; 4.6 ns, 2.595 Å; 4.9 ns, 2.568 Å | none |
| Which oxygen | the sulfonyl oxygen O2 every time (O1 never: 0 frames); no water and no protein oxygen within 2.6 Å in any all-atom frame (nearest water O 3.180 Å, nearest protein O 3.703 Å) | nearest water O 3.154 Å, nearest protein O 3.405 Å |
| Nearest sulfonyl O in the 500 solute frames: mean, s.d., smallest | 2.835, 0.124, 2.464 Å | 2.973, 0.123, 2.547 Å |
| Solute frames with it within 2.7 Å; within 2.8 Å | 74; 200 | 6; 45 |
| Mean angle S–N2–Zn | 109.34° | 113.99° (start X: 117.55°) |

Read plainly: the three all-atom frames are three of the thirteen solute frames, so the two counts measure the
same event, a sulfonyl oxygen that comes within 2.6 Å of the zinc for single frames; no water and no protein
oxygen is involved, and the zinc does not gain a fifth ligand that stays. With 50 frames the limit of 5 % is
passed by the third frame. By the rule as it was fixed before the run this is a failed gate, and it is counted
as one. It is also a real difference between the sets and not only a matter of counting: with set 2s the
sulfonyl oxygen is on average 0.14 Å nearer the zinc than with set 2e
(2.835 against 2.973 Å) and within 2.7 Å in 74 frames against 6,
and it goes with smaller mean angles N2–Zn–NE2 His96 (107.96° against 115.78°) and S–N2–Zn
(109.34° against 113.99°), the direction of the equilibrium angles of section 2.4.

## 6. The queue as added

Machine `cqc-md-gpu-01`, task folder `/root/cqc-sfa-md/queue/tasks.d/`; one GPU task at a time; each task waits
for the one named in `depends_on`. The files are in `md/gpu/` and have the same sha256 in the task folder.

| File | sha256 | Added (IST) | Tasks |
|---|---|---|---|
| `md/gpu/tasks_x0.json` | `0740b84bf9c2c8b7a1b0f4c59f7edfaa0686487c463fe8a5276bb0df77e25d42` | 6 Oct 2026, 9:36:57 PM | `X0_2s`, then `X0_2e` (both done) |
| `md/gpu/tasks_production.json` | `14d243f3727544ae79e1a91355573d3e4970541414d007fca14e93495ddb377a` | 6 Oct 2026, 10:14:10 PM | `X1`, `N1`, `R1`, `X2`, `N2`, `R2`, `X3`, `N3`, `R3` |
| `md/gpu/tasks_sensitivity.json` | `98ccd1715cb46d514643222e425fffd327b24d6e5161ed501a9054a8209c362e` | 6 Oct 2026, 10:14:10 PM | `R_zn`, `R_n2` |

| Order | Task | Waits for | System | Run folder on the machine | Run id | Length |
|---|---|---|---|---|---|---|
| 1 | `X1` | X0_2e, preflight_speed, preflight_interrupt | `systems/bonded_2e/X` | `/root/cqc-sfa-md/runs/bonded_2e/X1` | X1 | 100 ns in segments of 5 ns |
| 2 | `N1` | X1 | `systems/bonded_2e/N` | `/root/cqc-sfa-md/runs/bonded_2e/N1` | N1 | 100 ns in segments of 5 ns |
| 3 | `R1` | N1 | `systems/bonded_2e/R` | `/root/cqc-sfa-md/runs/bonded_2e/R1` | R1 | 100 ns in segments of 5 ns |
| 4 | `X2` | R1 | `systems/bonded_2e/X` | `/root/cqc-sfa-md/runs/bonded_2e/X2` | X2 | 100 ns in segments of 5 ns |
| 5 | `N2` | X2 | `systems/bonded_2e/N` | `/root/cqc-sfa-md/runs/bonded_2e/N2` | N2 | 100 ns in segments of 5 ns |
| 6 | `R2` | N2 | `systems/bonded_2e/R` | `/root/cqc-sfa-md/runs/bonded_2e/R2` | R2 | 100 ns in segments of 5 ns |
| 7 | `X3` | R2 | `systems/bonded_2e/X` | `/root/cqc-sfa-md/runs/bonded_2e/X3` | X3 | 100 ns in segments of 5 ns |
| 8 | `N3` | X3 | `systems/bonded_2e/N` | `/root/cqc-sfa-md/runs/bonded_2e/N3` | N3 | 100 ns in segments of 5 ns |
| 9 | `R3` | N3 | `systems/bonded_2e/R` | `/root/cqc-sfa-md/runs/bonded_2e/R3` | R3 | 100 ns in segments of 5 ns |
| 10 | `R_zn` | R3 | `systems/bonded_2s/R` | `/root/cqc-sfa-md/runs/bonded_2s/R1` | R1 | 100 ns in segments of 5 ns |
| 11 | `R_n2` | R_zn | `systems/bonded_2e_n2fit/R` | `/root/cqc-sfa-md/runs/bonded_2e_n2fit/R1` | R1 | 100 ns in segments of 5 ns |

Every command is that of `md/RUN.md` and `md/gpu/GPU_RUN.md`, section 5:
`run_md.py --system … --out … --run-id … --platform CUDA --stage all --prod-ns 100 --segment-ns 5 --resume`; the
test of a task is the run's `manifest.json` with `status=complete`, `production.complete=true`,
`production.completed_ns=100` and `test_only=false`. Each run has a time limit of 6 hours for one attempt
(`timeout_s`; a run takes about 2.7 hours): a run that hangs is ended and continued by the next attempt, three
attempts in all.

**The two sensitivity runs** (plan, entry of 9:15 PM, point 3; each once, start R, 100 ns, after the nine; in no
count of rule 3):

- `R_zn`: the other zinc set, here **2s**, which did not pass gate G-V (section 5.3); terms at N2 as built.
  Folder `runs/bonded_2s/R1`.
- `R_n2`: zinc set 2e as in the nine runs, with the terms at N2 of
  `md/qm/n2/n2_terms_new_charges_as_built.frcmod`, the candidate that gate G-R2 did **not accept**; it brackets
  the terms as built (plan: at the relaxed near-planar nitrogen the terms as built are 2.4 kcal mol⁻¹ too high
  and this candidate 1.9 to 2.7 too low). Folder `runs/bonded_2e_n2fit/R1`.

Both have run id R1 in a folder of their own model, so the analysis reads them as start R, replica 1, of the
models `bonded_2s` and `bonded_2e_n2fit` (the model is taken from the manifest), never as runs of `bonded_2e`;
and both have the seeds of run R1 of the primary set, so the three runs R1 start with the same velocities.
Each task writes a file `SENSITIVITY_RUN.txt` into its run folder before the run starts; its text is in
`md/gpu/tasks_sensitivity.json`. The system `bonded_2s_n2fit/R` is built and on the machine but in no task (it
would have been the system of `R_n2` had 2s been the primary set).

**Not queued and not decided here:** what follows if a sensitivity run and the nine runs disagree; that is the
plan's rule 4 and the user's to read.

## 7. The start of production, the speed and the end

Run X1 (start X, replica 1, `bonded_2e`) started on the machine at 10:14:22 PM IST on 6 Oct 2026 and was watched
from the Pi (the driver's log, copied to the Pi at 6 Oct 2026, 10:44 PM IST):

| Step | Driver's line |
|---|---|
| 10:14:51 PM IST | min: potential energy after -659872.1 kJ/mol (change -214652.6), 5077 calls, 25 s; L-BFGS iterations per pass 1066, 3998, 6, 7 |
| 10:15:05 PM IST | eq01 (heat): done in 13 s, 652.0 ns/day; T 294.4 K, density 0.9123 g/cm3 |
| 10:15:57 PM IST | eq02 (density): done in 51 s, 504.1 ns/day; T 298.8 K, density 1.0288 g/cm3 |
| 10:16:50 PM IST | eq03 (loosen): done in 52 s, 495.5 ns/day; T 303.1 K, density 1.0304 g/cm3 |
| 10:16:50 PM IST | eq03: hydration check: 0 water oxygen(s) within 2.5 A of the crystal ring site (nearest 3.26 A); ligand 0.26 A from the site on average; check skipped: the ligand itself is at the site; next step eq04 |
| 10:17:29 PM IST | eq04 (release_side_chains_zinc_ligand): done in 37 s, 922.3 ns/day; T 299.9 K, density 1.0252 g/cm3 |
| 10:18:07 PM IST | eq05 (loosen_backbone): done in 37 s, 922.7 ns/day; T 302.4 K, density 1.0275 g/cm3 |
| 10:18:55 PM IST | eq06 (free): done in 47 s, 921.8 ns/day; T 300.4 K, density 1.0252 g/cm3 |

| Segment | Ended (IST) | Wall time | ns per day |
|---|---|---|---|
| seg0001 (to 5 ns) | 6 Oct 2026, 10:26:43 PM IST | 467 s | 924.8 |
| seg0002 (to 10 ns) | 6 Oct 2026, 10:34:34 PM IST | 470 s | 918.9 |
| seg0003 (to 15 ns) | 6 Oct 2026, 10:42:23 PM IST | 468 s | 922.6 |

The hydration check does not apply to a run from start X (the ligand itself is at the site); it applies to the
runs from start R, where the held part of the equilibration is extended by up to 4 × 500 ps if no water is at
the place of the crystal ring (`md/RUN.md`, section 1A).

**Pulled to the Pi.** The watchdog's run that ended at 10:44 PM IST (pulls start every 15 minutes) fetched 17 files (111.9 MB); its line says "driver checksums: 158 ok, 0 bad". On the Pi, `/root/cqc-sfa-md/runs/bonded_2e/X1/` then held the minimisation, the six equilibration stages and segments seg0001, seg0002 (8 files of these segments; the sha256 of each was computed again here and equals the driver's record in the pulled manifest: true).

**The status program.** `python3 ops/pi/cqc_gpu_status.py` at 6 Oct 2026, 10:44 PM IST shows the progress of every run
and an end time (Queue (as pulled; the runner's state of 6 Oct, 10:40 PM): 22 done, 1 running, 10 waiting, 0 failed; the three lines of the rehearsal folders are left out
here, see section 9, point 6):

```
**MD runs** (speed now: 924.3 ns per day, from the last 5 finished production segment(s))
| Run | State | ns done / planned | ns per day (measured) |
|---|---|---|---|
| bonded_2e/X1 | **running** | 15 / 100 | 922.1 |
| bonded_2e/N1 | waiting | 0 / 100 |  |
| bonded_2e/R1 | waiting | 0 / 100 |  |
| bonded_2e/X2 | waiting | 0 / 100 |  |
| bonded_2e/N2 | waiting | 0 / 100 |  |
| bonded_2e/R2 | waiting | 0 / 100 |  |
| bonded_2e/X3 | waiting | 0 / 100 |  |
| bonded_2e/N3 | waiting | 0 / 100 |  |
| bonded_2e/R3 | waiting | 0 / 100 |  |
| bonded_2s/R1 | waiting | 0 / 100 |  |
| bonded_2e_n2fit/R1 | waiting | 0 / 100 |  |
| bonded_2s/X0 | done | 5 / 5 | 930.1 |
| bonded_2e/X0 | done | 5 / 5 | 925 |
| all 16 | 2 done | 25 / 1125 | |
Test and speed-test folders (not counted): bench, preflight_a, preflight_b.
**Production ends: about 8 Oct 2026, 3:35 AM IST**, if the machine runs through and the GPU does nothing else (1085 ns of production and 27 ns-equivalents of equilibration, less 1.48 ns of the step now running = 28.8 h at 924.3 ns per day).
Cost at that end: about US$36.14 of the cap US$60.00.
Not in this estimate: extensions of the held equilibration of the R runs (0 to 2 ns each), pre-emptions and their waiting times.
```

**Speed and end.** Measured: 3 finished segment(s) of X1 at 924.8, 918.9, 922.6 ns per day (1405 s for 15 ns, that is 922.4 ns per day); the two validation runs gave 930.15 and 925.03. One run of 100 ns then takes
9,640 s = 2.68 h (9,367 s of production at 922.4 ns per day, 237 s of equilibration,
36 s of minimisation and start). From the start of X1 at 10:14:22 PM IST on 6 Oct 2026:

| | Estimated end (IST) |
|---|---|
| the nine runs (R3 done) | **7 Oct 2026, 10:20 PM IST** |
| `R_zn` done | 8 Oct 2026, 1:01 AM IST |
| `R_n2` done, the queue empty | **8 Oct 2026, 3:41 AM IST** |

The status program's own estimate at that time: "Production ends: about 8 Oct 2026, 3:35 AM IST". Not in these times: extension steps of
the equilibration of the R runs (at most 4 × 500 ps, about 87 s each), a pre-emption of the spot machine
and the wait for its restart, and the time the Pi needs to pull the last files (the files of one segment
are 55.9 MB on the Pi, so a run of 20 segments is about 1.1 GB).

Cost: the status program gives about US$36.14 at its end time, against the cap of US$60.00.

## 8. How to follow and how to resume

On the Pi:

```bash
python3 ops/pi/cqc_gpu_status.py            # machine, cost, queue, every run with ns done, speed and the end time
tail -n 3 /root/cqc-sfa-md/status/watchdog.log
ls /root/cqc-sfa-md/runs/bonded_2e/         # the pulled runs; a finished segment arrives within about 15 minutes
```

On the machine (read only; `$SSH` as in `md/gpu/GPU_RUN.md`, section 2):

```bash
$SSH "sudo python3 /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md/gpu/queue.py status"
$SSH "sudo tail -n 3 /root/cqc-sfa-md/runs/bonded_2e/X1/run_md.log"     # last finished step and its ns per day
$SSH "sudo tail -n 5 /root/cqc-sfa-md/logs/X1.log"                      # the task's own output
```

**Nothing has to be done to resume.** A machine that was stopped or pre-empted is started again by the Pi's
watchdog, the runner starts by itself and starts the interrupted task again, and the run goes on after its last
complete segment (`--resume`; at most 5 ns are computed twice). This path was tested on the real machine at
9:21 PM IST on 6 Oct 2026 with the queue empty (plan, entry of 9:22 PM). A run whose command fails is tried
three times; after that the task is `failed` and the chain behind it waits: read `/root/cqc-sfa-md/logs/<task>.log`
and the run's `run_md.log`, then `sudo python3 …/md/gpu/queue.py retry <task>` on the machine. To take the
sensitivity runs out before they start: `queue.py cancel R_zn R_n2`.

**After the runs** (on the Pi, when every folder is pulled; `md/analysis/RECHECK.md`):

```bash
nice -n 10 /root/cqc-sfa-tools/md/bin/python md/analysis/recheck_md.py --compare-campaign /root/cqc-sfa-md/runs --run-primary
```

The campaign is evaluated per zinc model (`md/analysis/SPEC.md`, section 12), and the model of a run is read
from its manifest: `bonded_2s` and `bonded_2e_n2fit` are therefore models of their own there, with one run each
(R1) and an outcome that is marked as not final. They are the sensitivity runs and belong in no count of rule 3
of `bonded_2e`. The two validation runs X0 are left out of the campaign by the programs.

The file `/root/cqc-sfa-md/status/KEEP_RUNNING` (last touched at 10:39 PM IST on 6 Oct 2026, not by this task, which
touched it last at 10:13 PM) holds back the watchdog's stops for "all done" and "idle" for 6 hours after the last
touch (`ops/pi/GPU_WATCH.md`). It was not set by this task and is not removed by it; it runs out long before the
queue is empty, so the watchdog stops the machine by itself when every task is done and pulled, unless someone
touches the file again.

## 9. Defects found and fixed, and things to know

1. **Builder: a refused build could destroy a system that existed** (the accident of `md/N2.md`, section 5.6):
   the folder was emptied before the inputs were judged. Now every refusal comes before anything is touched,
   the system is written to `<start>/build/new/` and moved into place with atomic renames only when every step
   has ended; tested with a refused build in a scratch folder, which left the system as it was.
2. **Builder: the checks did not know a second zinc set or a variant of the terms at N2.** Added `--variant`,
   `--compare-with`, the reading of `IMPROPER` entries, the set-specific checks of section 3 and the
   term-by-term comparison of two systems. `--mcpb-dir` is made absolute; `--set 2s` refuses to run without
   `--n2-frcmod`, because MCPB.py's own lines at N2 differ between 2e and 2s (section 2.2).
3. **`md/qm/after_qm.sh`: a relative report folder is taken relative to the work folder**, because the script
   changes into the work folder first. The reports of set 2s first landed in
   `/root/cqc-sfa-md/qm/gr/mcpb_2s/md/qm/results/`; they were moved to
   `/root/cqc-sfa-md/qm/gr/after_qm_report_2s/` beside those of 2e. The script was not changed (give it an
   absolute path).
4. **Set 2s fails gate G-V** by one all-atom frame in part d of GV3 (section 5.3). Not a defect of the pipeline;
   listed here because it decides the primary set.
5. **Gate G-Z as worded can hardly fail on angles**: in the first 100 ps the zinc site is held to the crystal
   positions (section 4). The numbers of the free stages and of the 5 ns are given beside it for that reason.
6. **The status table on the Pi shows the queue as of the last pull** (up to 15 minutes old), and it lists the
   three rehearsal runs of `runs_test/` as "done, pull pending" with 0 of 5 ns, because the task file of the
   rehearsal is still in the queue and `runs_test/` is not pulled. These three lines are not production and the
   totals of the table include their 15 ns as planned. Taking `tasks_rehearsal_2e.json` out of the queue (and
   only then deleting the `runs_test` folders on the machine) removes them; a done task is tested again at every
   start of the runner, so the folders must not be deleted while the file is in the queue.

## 10. Files

| What | Where |
|---|---|
| This file | `md/PRODUCTION.md` |
| Builder | `md/build_topology.py` |
| Records of the eight systems | `md/topology/bonded_2e/`, `bonded_2s/` (X, N, R), `bonded_2e_n2fit/R/`, `bonded_2s_n2fit/R/` |
| Systems | `/root/cqc-sfa-md/systems/` on the Pi and on the machine (copied with a checksum list: `/root/cqc-sfa-md/status/push_production.sha256` on the machine, 57 files, all verified) |
| Step 2s | `/root/cqc-sfa-md/qm/gr/mcpb_2s/`, reports `/root/cqc-sfa-md/qm/gr/after_qm_report_2s/` |
| Parameter files of the two zinc sets in the repository | `md/qm/gr/charges/CA2_mcpbpy_2e.frcmod` and, added on 7 Oct 2026, `md/qm/gr/charges/CA2_mcpbpy_2s.frcmod`: copies of `CA2_mcpbpy.frcmod` of the two work folders, with the sha256 of section 2.1 (compared each time this file is assembled) |
| Blocks of the Hessian, table 2e / 2s / 2ms | `md/production/seminario_blocks.py`, `md/production/seminario_blocks.json` |
| Gate G-Z | `md/production/gate_GZ.py`, `md/production/gate_GZ_bonded_2s.json`, `md/production/gate_GZ_bonded_2e.json` |
| Validation runs and their analyses | `/root/cqc-sfa-md/runs/bonded_2s/X0/`, `/root/cqc-sfa-md/runs/bonded_2e/X0/` (`analysis/`, `analysis_recheck/`) |
| Task files | `md/gpu/tasks_x0.json`, `md/gpu/tasks_production.json`, `md/gpu/tasks_sensitivity.json` |
| How this file was assembled | `md/production/assemble/`: `frags.py` and `gates_md.py` (tables from the records, written to `frags.md` and `gates.md`), `x1_watch.py` (section 7, from `x1/run_md.log.machine`, the copy of the machine's log of X1, and `x1/status.txt`, the saved status table), `make_production_md.py` (this file; it stops if one of its own statements does not hold in the records), `make_tasks.py` (the two task files), `push/` (the log and the checksum lists of the copies to the machine) |
| Documents brought up to date | `md/RUN.md` (section 2; the path in section 4; section 7, points 4 and 11), `md/gpu/GPU_RUN.md` (sections 3 and 5); on 7 Oct 2026 `md/RUN.md` again, to the hydration check as it is run (six ring carbons, 2.5 Å), to the runs made on the GPU with their drivers and protocols (the outcome of the pre-flight test among them) and to the analysis specification as fixed (its header and state paragraph; section 1A; section 5, rows "Hydration check" and "GPU precision"; the first lines of section 6.2; section 7, points 1, 2, 6 and 10) |

This file as first written (state of 6 Oct 2026, 10:45 PM IST) is in commit `d517798` (6 Oct 2026, 10:47 PM IST),
which also holds the files of section 10 that lie under `md/` in their state of 6 Oct 2026, all but
`md/qm/gr/charges/CA2_mcpbpy_2s.frcmod`. Added on 7 Oct 2026, after the independent check of the set-up
(`paper/audit/research/md_production_check_2026-10-07.json`; plan, entry of 7 Oct 2026, 12:10 AM IST): the file
`md/qm/gr/charges/CA2_mcpbpy_2s.frcmod` and its row in section 10, the paragraph on crystal water in section 3,
the paragraph after the table of section 5.3, the pointer to it in the row "Set 2e" of section 1, the second
half of the row "Documents brought up to date" (the changes of `md/RUN.md` are named there), the sentence on
these additions at the top, and this paragraph, which replaces the line "Nothing of this is committed".
`md/production/assemble/make_production_md.py` was changed to write and to check these additions.
