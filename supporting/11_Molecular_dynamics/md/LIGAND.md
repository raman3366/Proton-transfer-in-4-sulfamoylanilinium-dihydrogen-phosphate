# Ligand parameters for the simulation (stage S1): what was found

Work of 5 Oct 2026, 8:05 PM to 10:15 PM IST, on the Pi (2 threads). Stage S1 of `PLAN.md`: the GAFF2 terms of
sulfanilamide at sulfur, checked against the DFT of the paper and repaired where gate G-L fails. Two species: the
sulfonamidate anion (charge −1, one H on the sulfonamide N; the species that will sit on the zinc) and the neutral
molecule (used here only as a guard against over-fitting). Nothing in `paper/`, `dock/` or `molecules/` was changed.
Every number below is in a file of `md/ligand/`; the file is named in brackets.

## 1. In short

1. **GAFF2 as assigned fails G-L for both species**, in vacuum and in implicit water (section 5).
2. **The scouting value "C–S–N 96°" was a tool fault.** In this AmberTools 26 build `parmchk2 -s gaff2` opens
   `gaff.dat` (GAFF 1.81) although the help text lists that spelling for GAFF2; only `parmchk2 -s 2` opens
   `gaff2.dat`. The frcmod of the first spelling, loaded after `leaprc.gaff2`, replaces 10 (anion) and 17 (neutral)
   native GAFF2 terms. With `-s 2` no native term is replaced (86 and 89 terms checked) [`gaff2_as_assigned.json`].
   True GAFF2 gives C–S–N 107.2° (anion) and 103.0° (neutral) in vacuum, not 96°. It still fails.
3. **The largest fault is the rotation about S–N, not the angle.** GAFF2 has no torsion term for H–N–S–X of these
   nitrogen types; parmchk2 substitutes one by analogy (anion: penalty score 160.5). The paper holds 36 DFT single
   points on rigidly rotated geometries. On the 19 of them where only the sulfonamide is turned, GAFF2 is wrong by
   3.9–4.1 kcal mol⁻¹ (r.m.s., implicit water) and 7.5–8.0 (vacuum); the repaired set by 0.5–0.8 (implicit water).
   GAFF2 puts the N–H hydrogen on or beside a sulfonyl oxygen, and the conformer of the DFT minimum is not a
   minimum at all [`rigid_single_points.csv`, `repair.json`].
4. **Repair.** Anion: three equilibrium values (S–N bond, C–S–N, S–N–H) and the S–N rotor (two Fourier terms on
   H–N–S–C4). Neutral: five equilibrium values, the rotor, and the C4–S torsion on the nitrogen path. Every repaired
   term contains the sulfonamide nitrogen; the shared sulfonyl terms stay GAFF2 for both species (section 4).
5. **After the repair G-L passes in implicit water for both species.** In vacuum the anion passes narrowly and the
   neutral molecule fails on the scan (2.19 against 0.5 kcal mol⁻¹). The DFT reference is in PCM water, so implicit
   water is the like-for-like comparison; see section 4.1 for why and section 7 for what this leaves open.
6. **G-L does not reach the coordinate that decides the simulation.** The two ring orientations of the docking
   differ in the rotation about N–S as seen from the zinc: Zn–N2–S–C4 is −117° in the crystal, −108° to −104° in the
   crystal-like docked poses and −35° to −32° in the rotated ones; the ring torsion C3–C4–S–N2 that G-L tests is
   the same in both (141–154°) [`bound_rotamer.csv`]. In the planned zinc model this rotation is governed by the
   H–N–S–C4 term (fitted here on the free anion), by the angles at N and by 1–4 non-bonded terms to the zinc. None of
   that is validated by a test of the free ligand. It needs a check on the zinc model before production (section 6).

## 2. Reference data (from the repository; no new quantum chemistry)

| Quantity | DFT anion, PCM | DFT neutral, PCM | DFT neutral, gas | Crystal 6RL9, SAN 304 |
|---|---:|---:|---:|---:|
| C4–S–N2 (°) | 109.83 | 109.45 | 107.39 | 108.80 |
| O1–S–O2 (°) | 113.75 | 119.46 | 122.95 | 116.83 |
| O–S–N2, mean (°) | 111.46 | 105.34 | 105.15 | 109.64 |
| C4–S–O, mean (°) | 104.94 | 108.46 | 107.66 | 105.77 |
| S–N2 (Å) | 1.580 | 1.677 | 1.686 | 1.596 |
| C4–S (Å) | 1.806 | 1.773 | 1.785 | 1.771 |
| S–O, mean (Å) | 1.498 | 1.469 | 1.458 | 1.436 |
| C3–C4–S–N2 (°) | 101.7 | 87.2 | 89.6 | 132.1 |
| C4–S–N2–H (°) | −68.7 | 61.1, −61.5 | 63.7, −63.7 | not seen by X-rays |

[`reference_geometry.csv`]. DFT: B3LYP-D3(BJ)/6-311++G(d,p), `molecules/sfa/results/sulfonamidate_pcm_d3/`,
`sulfanilamide_pcm_d3/`, `sulfanilamide_gas_d3/`. Crystal: heavy atoms of SAN 304 of `dock/pdb/6RL9.cif.gz`.

- **Relaxed scans** (`<case>/torsion_scan/`): C3–C4–S1–N2 held at 90…180° (18 points anion, 21 neutral on the C3
  branch), every other coordinate relaxed, each point started from the minimum with the sulfur side turned rigidly;
  IEF-PCM water. Anion: flat, 0.098 kcal mol⁻¹ between lowest (170°) and highest (110°) point. Neutral: minimum at
  90°, 1.88 kcal mol⁻¹ at 180° [`dft_scan_anion.csv`, `dft_scan_neutral.csv`]. Only the rotation about C4–S was
  scanned; no relaxed scan about S–N exists.
- **Rigid single points** (`<case>/dock_strain.json`, made for Table S27): the DFT minimum with the three rotors
  C4–S, S–N and C1–N1 turned rigidly, 16 geometries for the anion and 20 for the neutral molecule, same level, PCM.
  With the aniline as in the minimum they give: anion, N–H *anti* to C4 (176–179°), 4.64–6.00 kcal mol⁻¹ above the
  minimum; neutral, both N–H eclipsing the oxygens (±117–121°), 3.33–5.15 [`rigid_single_points.csv`]. These are the
  only DFT energies on the S–N rotor.
- **Hessian**: none stored. `sulfanilamide_gas_d3/freq.json` holds frequencies and intensities only, nothing for the
  anion. Angle force constants therefore stay GAFF2.

## 3. GAFF2 as assigned

antechamber (GAFF2 types, AM1-BCC charges, on the DFT minimum) [`sfa_anion.mol2`, `sfa_neutral.mol2`]:

- Anion: N2 gets type **n2** because antechamber draws S=N and one S–O⁻; the two oxygens get different charges
  (−0.736, −0.780). Neutral: N2 is **n8**. Aniline N: nv. Sulfur: sy.
- The AM1-BCC charges summed to −0.997 (anion) and +0.001 (neutral), because sqm prints the Mulliken charges to
  three decimals. The residual was spread evenly over all atoms; the delivered mol2 files sum to −1 and 0 exactly.
- Heavy-atom names are those of SAN in 6RL9 (sulfur is `S`); hydrogens as in the chemical component (H2, H3, H5, H6,
  HN11, HN12, HN21, HN22; the anion keeps HN21).

Terms at sulfur, anion [`terms_gaff2_anion.csv`; neutral in `terms_gaff2_neutral.csv`]:

| Term | Types | k | Equilibrium | Source |
|---|---|---:|---:|---|
| C4–S | ca-sy | 209.82 | 1.7913 Å | gaff2.dat, native |
| S–N2 | n2-sy | 461.21 | 1.5440 Å | gaff2.dat, native |
| S–O | o-sy | 614.68 | 1.4656 Å | gaff2.dat, native |
| C4–S–N2 | ca-sy-n2 | 101.19 | 112.23° | gaff2.dat, native (4 structures) |
| C4–S–O | ca-sy-o | 120.43 | 108.34° | gaff2.dat, native |
| N2–S–O | n2-sy-o | 132.73 | 122.18° | gaff2.dat, native (7 structures) |
| O–S–O | o-sy-o | 157.75 | 121.38° | gaff2.dat, native |
| S–N2–H | hn-n2-sy | 61.55 | 111.27° | gaff2.dat, native (1 structure) |
| C–C4–S–N2, C–C4–S–O | X-ca-sy-X | 1.30, n = 2, phase 180° | | gaff2.dat, generic |
| H–N2–S–C4, H–N2–S–O | hn-n2-sy-X | 6.8 (n = 1) + 0.5 (n = 3), phase 180° | | **parmchk2 by analogy with X-ne-sy-X, penalty 160.5** |

The six equilibrium angles at sulfur add up to 694.7°; a tetrahedral centre allows about 657°. The set is strained
by construction, and the C–S–N angle at the minimum is 5° below its own equilibrium value. Neutral: S–N 300.51 /
1.692 Å, C–S–N 102.01 / 103.18°, N–S–O 134.01 / 106.96°, S–N–H 56.36 / 109.66° (all native); H–N–S–X by analogy
with X-n3-sy-X (3.13, n = 2, phase 0°; penalty 0.0).

Minimum from 145 starts (the DFT minimum, and 144 starts with C4–S, S–N and the aniline pyramid turned and a 0.02 Å
jitter, seed 20261005), OpenMM Reference platform [`geometry_before.csv`, `minima_before.csv`]:

| | C–S–N (°) | S–N (Å) | C4–S–N2–H (°) | S–N–H (°) | Scan, largest deviation (kcal mol⁻¹) |
|---|---:|---:|---:|---:|---:|
| Anion, DFT | 109.83 | 1.580 | −68.7 | 107.5 | |
| Anion, GAFF2, vacuum | 107.25 | 1.509 | 116.6 | 95.7 | 0.25 |
| Anion, GAFF2, implicit water | 103.76 | 1.533 | 0.0 | 99.3 | 0.40 |
| Neutral, DFT | 109.45 | 1.677 | 61.1, −61.5 | 111.0 | |
| Neutral, GAFF2, vacuum | 102.96 | 1.593 | −127.1, 127.1 | 95.5 | 2.23 |
| Neutral, GAFF2, implicit water | 100.05 | 1.619 | 118.1, 11.7 | 97.4 | 2.10 |

- In every case the N–H hydrogens fall towards the sulfonyl oxygens (S–N–H 96–99° against 107–111°) and the S–N
  bond is shortened by 0.05–0.08 Å, in the neutral molecule to 0.10 Å below its own equilibrium value. The cause is
  the 1–4 attraction between H(N) and O with nothing in the torsion term to hold against it.
- The scans of GAFF2 are not like for like: the conformer of the DFT minimum does not exist in GAFF2, so the MM scan
  runs in another one, and for the neutral molecule it jumps between two (the two starts of a point differ by up to
  7.3 kcal mol⁻¹). This is the jump the synthesis saw in the scouting scan [`scan_before_*.csv`].
- Control, the scouting route (`-s gaff2` frcmod over gaff2.dat): C–S–N 89.3° (anion) and 102.3° (neutral) at the
  global minimum [`gaff2_as_assigned.json`, `control_s_gaff2`].

## 4. The repair

### 4.1 Environment of the fit

The DFT is in PCM water. The MM was therefore fitted and gated in implicit water (generalised Born OBC2, dielectric
78.3553 as in the PCM, no surface term). Fitting a vacuum model to energies in a dielectric would put the screening
by the solvent into the torsion terms, and the explicit water of the simulation would then count it a second time.
The size of that effect was measured: the cos φ coefficient of the anion's S–N rotor comes out as −0.26 kcal mol⁻¹ in
implicit water and −2.90 in vacuum, a difference of about 5 kcal mol⁻¹ between *syn* and *anti*
[`rotor_models.json`]. Vacuum results are reported beside. The plan did not fix the environment; this choice is mine and is listed in section 7.

### 4.2 What was changed

Anion [`sfa_anion_repair.frcmod`, `repair.json`]:

| Term | Types | GAFF2 | Repaired | DFT value met at the minimum |
|---|---|---:|---:|---:|
| S–N2 bond, r₀ | n2-sy | 1.5440 Å | **1.5906 Å** | 1.580 Å |
| C4–S–N2, θ₀ | ca-sy-n2 | 112.23° | **117.133°** | 109.83° |
| S–N2–H, θ₀ | hn-n2-sy | 111.27° | **118.059°** | 107.47° |
| H–N2–S–C4 | hn-n2-sy-ca | by analogy (6.8, n = 1; 0.5, n = 3) | **0.2624 (n = 1, phase 180°) + 2.7030 (n = 2, phase 0°)** | |
| H–N2–S–O | hn-n2-sy-o | by analogy (the same) | **0** | |

Force constants are unchanged. N2–S–O was fitted too, moved by 0.4° only, and went back to GAFF2 (122.18°); the
angle then misses DFT by 0.2°. The C4–S torsion passes with the GAFF2 term and was not touched.

Neutral [`sfa_neutral_repair.frcmod`]: n8-sy 1.692 → 1.7375 Å; ca-sy-n8 103.18 → 113.293°; n8-sy-o 106.96 →
109.767°; hn-n8-sy 109.66 → 118.443°; hn-n8-hn 106.38 → 111.247°; hn-n8-sy-ca 2.3128 (n = 1, 180°) + 2.2133 (n = 2,
0°); hn-n8-sy-o 0; ca-ca-sy-n8 1.4733 (n = 2, 180°) in place of the generic 1.30.

How:

- **S–N rotor**: linear least squares on the rigid DFT single points (MM single points on the same geometries; the
  bond and angle terms cancel because all geometries share one frame). Two coefficients on H–N–S–C4. A third
  coefficient absorbs the error of the aniline rotor in the 7 (anion) and 10 (neutral) geometries where the aniline
  is turned as well; it is reported in section 7 and **not** put into the force field. Fit quality, anion: r.m.s.
  0.75, largest 1.56, leave-one-out 1.02 kcal mol⁻¹ over 16 points; neutral: 0.61, 1.33, 0.87 over 20. A third
  term (n = 3) raised the leave-one-out error of the anion from 1.02 to 1.91 and lowered that of the neutral
  molecule from 0.68 to 0.51; the two-term form was kept for both [`rotor_models.json`].
- **Equilibrium values**: each moved until the minimum of the DFT conformer has the DFT value (to 0.01° and
  0.0001 Å). The repaired θ₀ of C–S–N is 117° to give 109.8°: it is an effective number inside a strained set of six
  angles, not a geometry.
- **C4–S torsion**: refitted only where the scan missed (neutral, 0.70 → 0.05 kcal mol⁻¹), against the relaxed scan
  (weight 1/0.1) and the rigid points (weight 1/0.5) together, the scan relaxed again after each fit until the
  coefficients stopped changing (4 cycles).
- The frcmod was written, the topology rebuilt from it with tleap, and all "after" numbers come from the rebuilt
  topology. Fitted model and rebuilt topology agree to 0.0001 kcal mol⁻¹. OpenMM and sander agree to 0.005
  kcal mol⁻¹ on relative energies of all four topologies, vacuum and implicit water [`engine_check.json`].

### 4.3 Over-fitting guard

No repaired term is shared between the two species (all carry n2 or n8), so the anion's repair cannot be tested on
the neutral molecule term by term. What can be said: the same procedure, with the same functional form, works on
both; the shared sulfonyl terms (o-sy-o, ca-sy-o, ca-sy, o-sy) were left alone and both species still come out
within 2° on O–S–O and C–S–O. One set of **nitrogen** terms cannot serve both species: the DFT values differ too much
(S–N 1.580 against 1.677 Å, O–S–N 111.5 against 105.3°).

## 5. Gate G-L

G-L: at the MM minimum C–S–N within 2° and S–N within 0.02 Å of DFT; scan within 0.5 kcal mol⁻¹ (largest deviation,
each profile relative to its lowest point). Minimum = lowest of 145 starts. Verdict taken in implicit water.

| Species | Force field | Environment | C–S–N, MM − DFT (°) | S–N, MM − DFT (Å) | Scan (kcal mol⁻¹) | G-L |
|---|---|---|---:|---:|---:|---|
| Anion | GAFF2 as assigned | implicit water | −6.07 | −0.047 | 0.40 | **fails** |
| Anion | GAFF2 as assigned | vacuum | −2.58 | −0.071 | 0.25 | fails |
| Anion | repaired | implicit water | +1.37 | +0.000 | 0.09 | **passes** |
| Anion | repaired | vacuum | −0.08 | −0.019 | 0.48 | passes, narrowly |
| Neutral | GAFF2 as assigned | implicit water | −9.40 | −0.059 | 2.10 | **fails** |
| Neutral | GAFF2 as assigned | vacuum | −6.49 | −0.084 | 2.23 | fails |
| Neutral | repaired | implicit water | −0.01 | −0.000 | 0.05 | **passes** |
| Neutral | repaired | vacuum | +1.03 | −0.020 | 2.19 | fails (scan) |

[`gaff2_as_assigned.json`, `repair.json`, `geometry_before.csv`, `geometry_after.csv`, `scan_before_*.csv`,
`scan_after_*.csv`; figure `torsion_profiles.png`.]

- Anion, repaired, implicit water: four minima within 0.19 kcal mol⁻¹ (the ring torsion is flat, as in DFT); the
  lowest has C3–C4–S–N2 160° and C–S–N 111.20°; the one reached from the DFT geometry has 99° and 109.82°. Other
  angles at sulfur: O–S–O +0.8°, O–S–N +0.7°, C–S–O −1.9° from DFT.
- The N–H torsion of the repaired anion is ±91° at the minimum, DFT −69°. The repaired rotor is 0.5–1.5 kcal mol⁻¹
  above its minimum at 75–60° [`sn_rotor_profile.csv`]. The two-term series cannot place it better.
- Neutral, repaired, implicit water: N–H torsions ±61.3° (DFT 61.1, −61.5), angle sum at N2 331.99° (331.99),
  barrier 1.885 (1.876).

## 6. For the zinc-site step

**Terms that carry the ligating nitrogen** (MCPB.py renames its type; all must arrive under the new type):

| Kind | Types (n2 → new type) | Value to arrive |
|---|---|---|
| bond | n2-sy | 461.21, 1.5906 Å (repaired) |
| angle | ca-sy-n2 | 101.19, 117.133° (repaired) |
| angle | hn-n2-sy | 61.55, 118.059° (repaired) |
| angle | n2-sy-o | 132.73, 122.18° (GAFF2 native) |
| torsion | hn-n2-sy-ca | 0.2624, phase 180°, n = 1 and 2.7030, phase 0°, n = 2 (repaired) |
| torsion | hn-n2-sy-o | 0 (repaired) |
| torsion | ca-ca-sy-n2 | generic X-ca-sy-X, 1.30, phase 180°, n = 2 (GAFF2; MCPB.py leaves it to the wildcard) |
| bond | hn-n2 | 454.59, 1.022 Å (GAFF2 native) |

How to carry them, read from the installed code (`pymsmt/lib.py`, `pymsmt/mcpb/gene_pre_frcmod_file.py`), not run:

1. Use `parmchk2 -s 2`, never `-s gaff2`.
2. MCPB.py input: `gaff 2` and `frcmod_files sfa_anion_gaff2.frcmod sfa_anion_repair.frcmod`, in this order. MCPB.py
   copies bond, angle and torsion parameters to the renamed type from gaff2.dat and these files; a later file
   replaces an earlier entry only under the same spelling of the types, which is why the repair frcmod uses the order
   of gaff2.dat. Each section of the file ends with one blank line, as MCPB.py requires.
3. `amber_modeling.py`, line 240, prints `source leaprc.gaff2` to the screen and not into the tleap input when
   `gaff` is 2. Add the line to the generated tleap input by hand.
4. After MCPB.py, compare the terms of the new type in its frcmod with the table above.

**Two checks on the zinc model that S1 cannot replace:**

- *Geometry.* The repaired θ₀ and r₀ are effective values that balance non-bonded pulls of the AM1-BCC charges of the
  free anion. With the RESP charges of the zinc site the balance shifts. After MCPB.py, minimise the small model in
  MM and compare C–S–N, S–N, O–S–N, S–N–H, the angle sum at N and Zn–N–S–C4 with its QM optimum; if C–S–N or S–N
  miss by more than the G-L tolerances, move the equilibrium values there with the iteration of `ligand_fit.py`
  (`fit_equilibria`).
- *Rotation about N–S.* Crystal: Zn–N2–S–C4 −117.0°, zinc eclipsing O2 (Zn–N2–S–O2 −3.0°, Zn···O2 2.99 Å), and the
  direction from N2 to Thr199 OG1 at +116.9° from C4, over the other oxygen. Rotated docked poses: Zn–N2–S–C4 −35° to −32°, zinc between C4 and an
  oxygen (80–82° from the nearer oxygen), Zn–N2–S 124–126° against 118° [`bound_rotamer.csv`, `bound_rotamer.json`].
  The step between the orientations is therefore a rotation about N–S (or an inversion at N), and its energy in
  the model comes from the H–N–S–C4 term, the angles at N and the 1–4 terms between zinc and O, C4. GAFF2 and the
  repaired rotor differ along H–N–S–C4 by up to 5.3 kcal mol⁻¹ on the free anion (at 0°: GAFF2 has its minimum
  there, the repaired set its maximum) [`sn_rotor_profile.csv`]. Either choice, carried over unchecked, could decide
  the outcome by construction. Proposed: about 8–12 QM single points of the small zinc model with the SO₂–aryl part
  turned rigidly about N–S, at the level of the Hessian, compared with MM single points on the same geometries
  (the method of section 4.2), and the H–N–S–C4 coefficients refitted on the zinc model if they miss. I suggest this
  as an additional gate of S2.

## 7. What is uncertain

1. **Environment.** Fitted in implicit water against PCM. In vacuum the repaired neutral molecule has a barrier of
   4.06 against 1.88 kcal mol⁻¹; no gas-phase DFT scan exists to say how much of that is a real solvent effect. The
   one gas-phase DFT geometry (neutral) is not reproduced by the vacuum MM: C–S–N 110.48 against 107.39°, S–N 1.657
   against 1.686 Å. The bound anion sits on zinc, half buried: neither environment is its own.
2. **S–N rotor data.** Rigid single points, not relaxed ones; with the aniline unrotated the anion has only three
   distinct N–H positions (−69°, 46°, 176–179°). Leave-one-out error about 1 kcal mol⁻¹. The DFT minimum itself may
   not be the global one: microwave spectroscopy finds the *eclipsed* NH₂ rotamer for gas-phase sulfanilamide
   (Vigorito et al., *Molecules* 2022, 27, 2820, doi:10.3390/molecules27092820; abstract read only), while the
   paper's gas-phase minimum is the staggered one (N–H at ±63.7°). The repaired MM puts the eclipsed rotamer 3.63
   kcal mol⁻¹ above the staggered one in implicit water and 0.59 below it in vacuum [`minima_after.csv`], which agrees
   with both observations, but no DFT energy of a relaxed eclipsed rotamer exists to confirm it. This touches the
   paper's M4 gas geometry and is for the main session to judge.
3. **Transfer to the bound form.** See section 6. The nitrogen of the free anion carries two lone pairs; on zinc it
   carries one. The fitted rotor is evidence about the first, not the second.
4. **Charges.** AM1-BCC on the free anion. With the two oxygen charges set to their mean the repaired anion keeps
   C–S–N 109.82°, S–N 1.5798 Å and a scan deviation of 0.09 kcal mol⁻¹ [`repair.json`, `charge_sensitivity`]. No
   test with RESP charges was possible here.
5. **Force constants** are GAFF2's; no Hessian is stored.
6. **Left as GAFF2, with a miss:** S–O 1.468 Å against DFT 1.498 (the 1.0 Å structure has 1.436, so GAFF2 lies
   between them); C–S–O −1.9° (anion); sulfur lies in the ring plane (0.1–0.4° against 2.1–3.3° in DFT and 1.9° in
   the crystal; parmchk2 gives the improper at C4 10.5 kcal mol⁻¹ by analogy, penalty 41.9).
7. **Aniline group, not repaired, outside the gate.** The aniline N is planar in MM (angle sum 359.1–360.0°; DFT
   343.8° anion, 349.1° neutral), and the rigid points ask for about 0.9 kcal mol⁻¹ more on each of the four
   ca-ca-nv-hn terms (GAFF2 1.05, by analogy), i.e. the NH₂ rotation barrier of GAFF2 is about half of the DFT one.
   This can matter for the aniline hydrogen bond to Asn67 in the rotated pose. A repair would be one term; it needs
   a decision, because it is outside what S1 was asked to change.
8. **The DFT reference itself** moves with the basis set: C–S–N of the neutral molecule is 109.45° here and 108.85°
   with def2-TZVP (`paper/main/data/facts.json`, `frag.sulfanilamide_pcm_d3_tz.ang_CSN`).

## 8. Files and how to run again

Scripts (`md/`), each writes the same files on a second run (checked: `SHA256SUMS`; the first line of a prmtop
carries the build date):

| Script | Does | Time |
|---|---|---|
| `ligand_params.py` | reference tables, antechamber, parmchk2, tleap, term table, minima and scans of GAFF2 | 8 min |
| `ligand_fit.py` | the repair, rebuilt topologies, "after" tables, rigid points, charge check | 20 min |
| `ligand_plot.py` | `torsion_profiles.png` | seconds |
| `ligand_engine_check.py` | OpenMM against sander | seconds |
| `ligand_sn_rotor.py` | relaxed MM profile about S–N, both force fields | 15 min |
| `ligand_rotor_models.py` | the rotor fit in vacuum and implicit water, with two and three terms | seconds |
| `ligand_bound_rotamer.py` | crystal and docked poses seen from the zinc (docking environment) | seconds |
| `ligand_common.py`, `ligand_eval.py` | helpers: geometry, readers, MM engine, minima, scan, gate | |

Run with `/root/cqc-sfa-tools/md/bin/python` from the repository root (`ligand_bound_rotamer.py` with
`/root/cqc-sfa-tools/dock/bin/python`). Software: AmberTools 26.0, gaff2.dat 2.2.30 (sha256 in the JSON files),
ParmEd 4.3.1, OpenMM 8.6.1 (Reference platform), NumPy 2.5.3, SciPy 1.18.1. Scratch of antechamber, sqm and tleap:
`/root/cqc-sfa-md/ligand/`.

Deliverables in `md/ligand/`: `sfa_anion.mol2`, `sfa_anion_gaff2.frcmod`, `sfa_anion_repair.frcmod`,
`sfa_anion_gaff2.prmtop`, `sfa_anion_repaired.prmtop` (and `.inpcrd`), `sfa_anion_repaired_min.xyz`; the same for
`sfa_neutral`; `terms_gaff2_*.csv`, `terms_repaired_*.csv`; `reference_geometry.csv`, `dft_scan_*.csv`;
`geometry_before.csv`, `geometry_after.csv`, `minima_before.csv`, `minima_after.csv`; `scan_before_*.csv` (vacuum;
`_obc2` implicit water), `scan_after_*.csv` (implicit water; `_vacuum`); `rigid_single_points.csv`;
`sn_rotor_profile.csv`, `rotor_models.json`; `bound_rotamer.csv`, `bound_rotamer.json`; `gaff2_as_assigned.json`,
`repair.json`, `engine_check.json`; `torsion_profiles.png`; `SHA256SUMS`.
