# MD system of the hCA II–sulfanilamide complex: protonation, water, the three starts

Stage S3 of `md/PLAN.md`, first part (structure files; no solvation, no equilibration, no parameters). Written on
5 Oct 2026, about 8:30 PM IST. Nothing in `paper/`, `dock/`, `lit/` or `molecules/` was changed.

**Gates.** tleap loads the protein with the kept crystal water of each start: **passed**. The two orientation
distances of section 5 of the plan: **confirmed** for the crystal start and for the rotated start (section 6 below).

## How it is made

| Step | Command | Writes |
|---|---|---|
| 1 | `/root/cqc-sfa-tools/dock/bin/python md/build_system.py` | `md/system/`: `protein_md.pdb`, `zinc.pdb`, `waters_all.pdb`, `waters_X/N/R.pdb`, `ligand_X/N/R.sdf` and `.pdb`, `starts.csv`, `residue_map.csv`, `decisions.json` |
| 2 | `/root/cqc-sfa-tools/dock/bin/python md/system_audit.py` | `md/system/audit/`: the evidence (histidines, Glu106, PROPKA with and without zinc, conformers, water, all docked modes) |
| 3 | `AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/check_tleap.py` | `md/system/tleap/`: leap inputs and logs, `protein_H.pdb`, `tleap_check.json`; test topologies in `/root/cqc-sfa-md/system/tleap/` |
| 4 | `/root/cqc-sfa-tools/dock/bin/python md/verify_system.py` | `md/system/verify.json`: independent re-derivation and a second run of steps 1 and 2 |

`md/system_common.py` holds the shared code (reading of the entry, the conformer rule, geometry). No step uses
random numbers; a second run of steps 1 and 2 gives byte-identical files (`verify.json`). Software: gemmi 0.7.5,
RDKit 2026.03.1, NumPy 1.26.4, PROPKA 3.5.1 (docking environment); AmberTools 26.0 (tleap), ParmEd 4.3.1 (MD
environment). Source checksums are in `decisions.json`.

## 1. What the entry holds (`audit/deposited.json`)

6RL9, 1.0 Å, one chain. 258 residues are modelled, His3 to Lys261 (there is no residue 126; Thr125 and Lys127 are
joined by a normal peptide bond). Seven residues of the construct before His3 are not modelled (GSPEF of the tag,
Met1, Ser2). The entry carries 2,369 riding hydrogens on the protein, none on the ligands or the water. Hetero
groups: Zn 301; sulfanilamide 304 on the zinc (occupancy 1.00, one conformer, B 5.7–12.8 Å²); a second
sulfanilamide 305 (0.68) on the surface, 14 Å away; a third, 306 (0.85), at the rim of the active site, 5.7 Å from
the zinc-bound one; three glucose molecules; mercury on Cys206 (4-mercuribenzoate 310, component BE7, 0.76; Hg 302, 0.24)
and Hg 303 (0.05) at His64; 216 water sites, 10 of them with a conformer label. 56 residues have alternate
conformers. Distances and neighbours of every hetero group are in `audit/deposited.json`.

Zinc site as deposited: Zn–N(sulfonamide) 1.931 Å, Zn–NE2(His94) 2.010, Zn–NE2(His96) 2.021, Zn–ND1(His119)
2.018 Å; N–Zn–N angles 98.4–116.8°; Zn···O2 2.99 Å, Zn···O1 4.05 Å.

## 2. Protonation at pH 7.4

Evidence used: the geometry of 6RL9 (`audit/histidines.json`, `audit/glu106.json`), the riding hydrogens of the
entry, PROPKA 3.5.1 run on the docking's heavy-atom protein without and with the zinc (`audit/propka.json`), and
four papers read in full text: NMR and neutron proton occupancies of the 12 histidines (Fisher, Raum and Weininger
2025, DOI 10.1002/cbic.202400930), the NMR study of His64 (Raum, Fisher and Weininger 2023, DOI
10.1007/s00018-023-04936-z), the neutron structure at pH 7.8 (Fisher et al. 2011, DOI 10.1021/bi201487b, author
manuscript) and the neutron structure with acetazolamide (Fisher et al. 2012, DOI 10.1021/ja3068098). All four DOIs
resolve on Crossref. The NMR values quoted, with the protonated fraction at pH 7.4 computed from them, are in
`decisions.json` (`protein.histidines_nmr`).

### Histidines

| Residue | Here | Docking | Evidence |
|---|---|---|---|
| His94 | HID | HID | binds zinc through NE2 (2.01 Å); ND1–H donates to Gln92 OE1 (2.79 Å). NMR and neutron: N(δ)–H only |
| His96 | HID | HID | binds zinc through NE2 (2.02 Å); ND1–H donates to the carbonyl O of Asn244 (2.80 Å). NMR and neutron: N(δ)–H only |
| His119 | HIE | HIE | binds zinc through ND1 (2.02 Å); NE2–H donates to Glu117 OE2 (2.72 Å). NMR and neutron: N(ε)–H only |
| His64 | HIE, conformer A (out) | HIE, conformer A | **assumption**, see below |
| His107 | **HID** | HIE | buried. ND1–H donates to Glu117 OE2 (2.76 Å), NE2 accepts from Tyr194 OH (2.63 Å). NMR: N(δ)–H only, pKa < 3.5. The entry carries HD1 |
| His122 | **HID** | HIE | buried. ND1–H donates to the carbonyl O of Ala142 (3.02 Å), NE2 accepts from Tyr51 OH (2.74 Å). NMR: N(δ)–H only, pKa < 3.5. The entry carries HD1 |
| His4 | **HID** | HIE | ND1 is 2.70 Å from the carbonyl O of Trp5 (a hydrogen there: 1.84 Å, 142°); NE2 has no partner; the entry carries HD1. NMR in solution: 60–72 % N(ε)–H, pKa 6.2. Surface, B 29 Å², 15 Å from the zinc |
| His3 | HIE | HIE | not decided by the structure (NE2···water 2.86 Å). NMR: 71–85 % N(ε)–H, pKa 6.1 |
| His10 | HIE | HIE | conformer A, no partner. NMR: 79–93 % N(ε)–H, pKa 5.6 |
| His15 | HIE | HIE | conformer A: ND1 accepts from Lys18 NZ (2.94 Å), NE2–H donates to the carbonyl O of Lys9 (2.74 Å). NMR 76–90 % N(ε)–H |
| His17 | HIE | HIE | NE2–H donates to the carbonyl O of Glu187 (2.86 Å). NMR finds both tautomers |
| His36 | HIE | HIE | conformer A: ND1···water 2.59 Å. NMR: 57–69 % N(ε)–H, pKa 6.94, so about 26 % is protonated at pH 7.4; the neutral majority is taken. 25 Å from the zinc |

**His64.** The entry has three conformers: A, pointing away from the zinc (out, occupancy 0.55, χ1 −55°, 11.0 Å
from the zinc), and B and C, pointing towards it (in, 0.39 and 0.05, 7.6 and 8.3 Å). Conformer A is taken, as in
the docking. In solution (Raum 2023) the two positions are populated about 1:1 at every pH and exchange at
10⁷–10⁹ s⁻¹; the pKa is 6.25 for the inward and 7.60 for the outward position (6.27 and 7.43 with acetazolamide
bound), and the two neutral tautomers are about 1:1. From these values 7 % of the inward and 52–61 % of the outward
form are protonated at pH 7.4, about 30 % overall. One fixed state has to be chosen for a simulation: the neutral
form is the majority, the tautomer is not determined (ND1···water 2.73 Å), and HIE is kept as in the docking
receptor and in the entry. This is the least certain protonation choice of the system.

### Glu106 (changed) and Glu117

Glu106 is **deprotonated** here; the docking receptors have it protonated.

- PROPKA gave 8.91 because the zinc was not in its input. On the same coordinates with the zinc as a +2 ion it gives
  **5.33**; the run without the zinc reproduces 8.91 (`audit/propka.json`).
- OE1 accepts a hydrogen bond from Thr199 OG1 (2.62 Å; the riding hydrogen of the entry points at OE1, 1.81 Å,
  162°). OE2 accepts from the backbone N–H of Arg246 (2.79 Å) and from water 415 (2.65 Å). The zinc is 4.00 Å from
  OE1.
- Thr199 OG1 has one hydrogen. It accepts the N–H of the sulfonamide (N···OG1 2.86 Å), so its own hydrogen must go
  to Glu106. This is the scheme seen by neutron diffraction with acetazolamide: the sulfonamide N–D donates to
  Thr199, which donates to Glu106 (Fisher 2012).
- The two C–O bonds are equal (1.246 and 1.249 Å). This is weak evidence: the geometry restraints of the
  refinement pull every carboxylate towards 1.25 Å (mean difference over 21 fully occupied Asp and Glu 0.011 Å,
  largest 0.029 Å).

What the structure alone does not exclude: a hydrogen on OE2 pointing at water 415. I found no direct measurement
of the Glu106 charge in a text I could read. Stage S3b tests what the change does to the docking.

Glu117 is deprotonated (PROPKA 0.73 without and 3.60 with the zinc); it accepts four hydrogen bonds (His119 NE2–H,
His107 ND1–H, the N–H of Glu106 and of His107).

### Everything else

- Asp and Glu: all deprotonated; highest PROPKA pKa apart from Glu106 is 4.83 (Glu187).
- Lys and Arg: all protonated; lowest PROPKA pKa 9.49 (Lys154) and 10.38 (Arg246).
- Tyr: all neutral; lowest PROPKA pKa 10.43 (Tyr40). Tyr7 (PROPKA 13.5–14.7) is protonated at pH 7.8 in the
  neutron structure and deprotonated only at pH 10 (Fisher 2011).
- Cys206, the only cysteine: free thiol (PROPKA 10.7). In the crystal it carries mercury from an additive; the
  mercury compounds are removed and conformer A (0.76) is taken. Same as in the docking.
- PROPKA flags, within one unit of pH 7.4: Glu106 (above), His36 (6.72; treated above), His107 (6.43 with the
  zinc; NMR says below 3.5) and the amino group of His3 (7.78; capped, see below).
- Asn and Gln: the 21 amide orientations of the entry were checked against their neighbours; none is suspicious
  (`audit/amides.json`). Asn67 OD1, used in the orientation coordinate, accepts from water 423 (2.64 Å).

Formal charge of the protein: 19 Asp + 13 Glu + the C-terminus against 24 Lys + 7 Arg, **−2**. With the zinc (+2)
and the anion (−1) the complex is −1.

## 3. Conformers, termini, atoms without density

**Rule** (the rule of the docking, with coupling added): residues with alternates form a group when they are
neighbours in the sequence or when atoms of two different labels come closer than 2.4 Å; a group takes one label,
the one with the most modelled atoms, then the highest occupancy, then A. This gives 37 groups for the 56 residues
(`audit/altlocs.csv`, `decisions.json`). 32 groups are decided by occupancy, 4 by the number of atoms, 1 by the tie
rule.

Differences from the docking receptor: Ser73, Gln74, Asp75, Ala77, Val78, Leu79, Lys80, Gly81 and Arg89 take B
(docking: A), and Phe176 takes A (docking: B). In the loop 73–81 the two conformers differ by up to 2.05 Å in the
backbone, and the docking had taken Lys76 from B between neighbours from A; Asp75 (A) and Arg89 (B) come within
2.29 Å, so Arg89 follows the loop. All ten residues are at least 10.4 Å from the ligand. Tyr7 takes B (0.61), His64
A (0.55), Cys206 A (0.76), as in the docking.

**N-terminus: acetyl cap on His3** (docking: NH₃⁺ on His3). The chain continues before His3, so a charged amino
group there would be a charge the protein does not have. The three cap atoms are placed by geometry, with φ of His3
(−100°) chosen in 10° steps to keep them farthest from all other atoms (nearest atom 3.31 Å). **C-terminus:**
carboxylate on Lys261, the true end of the chain (OXT is in the entry).

**22 side-chain atoms without density** are not written; tleap builds them from its templates: Lys9 (CG, CD, CE,
NZ), Lys45 (CE, NZ), Asp52 (OD1, OD2), Lys80 (CE, NZ), Asp165 (CG, OD1, OD2), Lys168, Lys170, Lys213 (NZ), Lys228
(CE, NZ), Lys261 (CG, CD, CE, NZ). All are on the surface. After tleap the closest contact of a built atom is
Asp165 OD1 at 2.36 Å from water 419; no built atom is closer than 2.3 Å to anything.

**Removed:** the second and third sulfanilamide, the three glucose molecules and the mercury compounds. The third
sulfanilamide (306, occupancy 0.85) lies at the rim of the active site, its sulfonyl O 3.07 Å from Asn67 ND2, and
touches a neighbouring molecule of the crystal. The simulation holds one ligand, so its place is left to the
solvent.

## 4. Water

All crystal water of chain A is kept, with these rules (`decisions.json`, `audit/waters.csv`):

- 206 sites have no label and full occupancy: kept.
- Of the 10 labelled sites the copy is kept that does not clash (2.4 Å) with a conformer taken; of several free
  copies the one with the highest occupancy. Kept: 410 A, 412 A, 414 A, 417 A, 425 B, 481 B, 505 A, 513 A, 570 B.
  Dropped: 510 (its only copy, B, 0.32, is 1.78 Å from Glu214 OE2 of conformer A).
- Result: **215 of 216** sites. Two pairs of kept waters are closer than 2.4 Å (427–569, 2.37 Å; 502–580, 2.39 Å);
  both are in the entry at full occupancy and are left as they are.
- For each start a water whose O is within 2.2 Å of a ligand heavy atom is removed.

| Start | Waters removed | Waters kept | Nearest kept water to the ligand |
|---|---|---|---|
| X | none | 215 | 553, 2.93 Å from N2 |
| N | none | 215 | 553, 3.15 Å from C5 |
| R | **423, 516, 553, 568** | 211 | 582, 2.67 Å from C2 |

The rotated pose lies on four water sites of the crystal: its ring atoms are 0.19 Å (568), 0.55 Å (553) and 1.43 Å
(516) from them, and its aniline N is 0.28 Å from water 423, the water that bridges Asn67 OD1 (2.64 Å) and Asn62
ND2 (3.00 Å). All four are fully occupied with B 11.5–17.6 Å² (median of all water 22.5 Å²). The docking receptor
held no water, so the docking could not see this.

**Deep water.** There is none: no water lies within 4 Å of Thr199 N, and the sulfonyl O1 of the ligand takes that
place (Thr199 N···O1 2.94 Å). No water is in the first shell of the zinc; the nearest are 553 (3.69 Å) and 516
(4.89 Å).

**Aniline end.** Crystal pose: one water within 4 Å of the aniline N (582, 3.80 Å, B 32 Å²). Near-native pose:
none. Rotated pose, after the removal: 568 would be 2.56 Å and 516 2.86 Å away (both removed), 582 is 3.29 Å away.

tleap gives every water its two hydrogens in one standard orientation. It reports 45 pairs of hydrogens closer than
1.5 Å (shortest 0.90 Å), nearly all with a water hydrogen, and no close contact that involves a heavy atom. The
minimisation has to relax them.

## 5. The three starts

Species: the sulfonamidate anion ArSO₂NH⁻, charge −1, 18 atoms, residue SAN 304. Heavy-atom names are those of
SAN in 6RL9; the paper's C3–C4–S1–N2 is C3–C4–S–N2 here. Hydrogens: H2, H3, H5, H6, HN11, HN12, HN21. The atom
order of the files is that of the paper's DFT files and of the ligand parameters of stage S1
(`md/ligand/sfa_anion.mol2`).

**Selection rule** (`decisions.json`, `ligand.pose_rule`). Both docked starts come from the production job of the
paper (AutoDock4Zn, neutral sulfanilamide from its DFT geometry, 20 seeds, 138 modes;
`dock/T3_6RL9/runs/sfa_neutral__ad4zn__prod`). A mode must pass the paper's site check (Zn–N ≤ 2.3 Å, Thr199
OG1···N ≤ 3.3 Å).

- **N**: the best-scored mode within 1.0 Å of the crystal pose.
- **R**: the best-scored mode more than 2.0 Å away with the aniline N within 3.5 Å of Asn67 OD1.

This is the rule proposed for the task. I wrote it down after all modes had been tabulated, so it was not fixed
blind. The clause on Asn67 is there because 24 of the 62 modes beyond 2.0 Å are other orientations; it does not
change the pick, the best-scored mode beyond 2.0 Å is the same one.

In this job the near-native pose is ranked first in all 20 seeds and the rotated pose third in all 20, 0.014–0.070
kcal mol⁻¹ behind the best near-native mode of its seed (`audit/poses_summary.json`). N and R therefore share one
input geometry and differ by 0.020 kcal mol⁻¹ in score.

| | X (crystal) | N (near-native) | R (rotated) |
|---|---|---|---|
| Source | SAN 304 of 6RL9 | seed 3, mode 1 | seed 14, mode 3 |
| File | `dock/pdb/6RL9.cif.gz` | `…/seed03/out.pdbqt` | `…/seed14/out.pdbqt` |
| AutoDock4Zn score (kcal mol⁻¹) | – | −8.843 | −8.823 |
| RMSD to the crystal pose (Å) | 0 | 0.61 | 3.57 |
| Zn–N(sulfonamide) (Å) | 1.93 | 2.21 | 2.24 |
| Zn···O(sulfonyl), nearer / farther (Å) | 2.99 / 4.05 | 2.78 / 4.24 | 3.91 / 4.53 |
| N2···OG1(Thr199) (Å) | 2.86 | 2.55 | 2.59 |
| H···OG1 (Å), N–H···OG1 (°) | 1.87, 164 | 1.62, 149 | 1.72, 140 |
| Thr199 N···O(sulfonyl) (Å) | 2.94 | 2.98 | 3.02 |
| C3–C4–S–N2 / C5–C4–S–N2 (°) | 132.1 / −50.1 | 148.6 / −28.9 | 142.7 / −34.9 |
| τ_fold (°) | 132.1 | 148.6 | 142.7 |
| C4–S–N2 (°) | 108.8 | 109.5 | 109.4 |
| S–N2 (Å) | 1.596 | 1.677 | 1.677 |
| aniline N to Phe131 CZ (Å) | 3.90 | 3.36 | 9.58 |
| aniline N to Asn67 OD1 (Å) | 7.74 | 7.73 | 2.60 |
| aniline N to the His64 ring, conformer A / B / C (Å) | 9.4 / 6.8 / 6.9 | 10.2 / 7.4 / 7.5 | 5.4 / 3.3 / 3.5 |

All numbers are in `starts.csv` and `decisions.json`; distances are to the deposited coordinates. In the frame of
the docking receptor (heavy atoms moved by its minimisation) the two orientation distances differ by about
0.02 Å at most. `verify_system.py` recomputes the table from the written files with its own code (distances within
0.002 Å, angles within 0.06°, the effect of three-decimal coordinates).

**Hydrogens.** The entry has no ligand hydrogens. In all three starts the single sulfonamide hydrogen is placed by
one construction: N2 is taken as tetrahedral with S and Zn as two substituents; of the two remaining directions the
one nearer to Thr199 OG1 gets the hydrogen (N–H 1.021 Å, the value of the paper's DFT anion), the other is the lone
pair. In the crystal pose this gives a good hydrogen bond (1.87 Å, 164°).

In the two docked poses the neutral molecule has two sulfonamide hydrogens, and **neither donates to Thr199**:

| Pose | Hydrogen | H···Zn (Å) | Zn–N–H (°) | H···OG1 (Å) | N–H···OG1 (°) | Fate |
|---|---|---|---|---|---|---|
| N | 14 | 1.88 | 57.9 | 2.21 | 97.5 | points at the zinc: deleted |
| N | 15 | 3.08 | 143.0 | 2.07 | 106.7 | replaced by the construction (moved 1.02 Å) |
| R | 15 | 2.16 | 72.6 | 2.14 | 104.2 | points at the zinc: deleted |
| R | 14 | 2.93 | 124.4 | 2.20 | 100.5 | replaced by the construction (moved 1.27 Å) |

The plan says the hydrogen that donates to Thr199 is kept. No docked hydrogen meets a hydrogen-bond criterion
(H···O ≤ 2.5 Å and N–H···O ≥ 120°), so the remaining one was replaced by the same construction as in X. The other
hydrogens: ring hydrogens by RDKit; aniline hydrogens by RDKit for X (planar NH₂) and as docked for N and R
(pyramidal NH₂ of the DFT geometry).

What N and R carry from the docking and the force field must relax: the bond lengths of the neutral molecule (S–N
1.677 Å against 1.596 Å in the crystal and 1.580 Å in the DFT anion), a Zn–N distance of 2.2 Å and a short
N2···OG1 contact (2.55–2.59 Å). With the hydrogens tleap gives the protein, the shortest contacts of a ligand
hydrogen are 1.66 Å (N, H2 to a methyl hydrogen of Val121) and 1.69 Å (R, an aniline hydrogen to HD21 of Asn67);
no ligand hydrogen is closer than 2.08 Å to a kept water oxygen (`verify.json`, for information, not a gate).

## 6. The numbers of section 5 of the plan

| Distance | Plan, crystal | X here | Plan, rotated | R here |
|---|---|---|---|---|
| aniline N to Phe131 CZ | 3.9 Å | 3.90 Å | 9.3–9.8 Å | 9.58 Å |
| aniline N to Asn67 OD1 | 7.7 Å | 7.74 Å | 2.6–2.7 Å | 2.60 Å |

**They match.** Two remarks for the rule, which I did not change:

- The plan gives ranges for the rotated pose, so the midpoints are fixed only to about ±0.13 Å and ±0.03 Å. From
  the X and R starts they are **6.74 Å** and **5.17 Å**. The near-native start (3.36 and 7.73 Å) is on the crystal
  side of both.
- The ranges of the plan are not the range of all rotated modes: over the 38 rotated modes of the production job
  the distances are 8.72–10.32 Å and 2.56–3.20 Å. All of them are on the rotated side of both midpoints.

## 7. Files and the tleap test

`protein_md.pdb` has 2,040 heavy atoms (258 residues and the cap), deposited coordinates unchanged, deposited residue
numbers, Amber residue names (5 HID, 7 HIE, no HIP, GLU for Glu106, CYS for Cys206). `residue_map.csv` gives the
index tleap assigns: deposited number − 1 up to Thr125, − 2 from Lys127 on (ACE is 1, His94 is 93, His96 95,
His119 118, Thr199 197, Lys261 259); water follows in file order.

Test (`tleap/tleap_check.json`), with leaprc.protein.ff14SB and leaprc.water.tip3p, zinc and ligand left out:

| Load | Atoms | Waters | Net charge | Heavy atoms added | Errors |
|---|---|---|---|---|---|
| protein alone | 4,075 | 0 | −2.000 | 22, as listed | 0 |
| protein + water of X | 4,720 | 215 | −2.000 | 22, as listed | 0 |
| protein + water of N | 4,720 | 215 | −2.000 | 22, as listed | 0 |
| protein + water of R | 4,708 | 211 | −2.000 | 22, as listed | 0 |

No atom name was unknown to the templates, the input heavy atoms kept their coordinates (largest shift 0.000 Å),
and tleap's warnings are the non-zero charge and the hydrogen contacts named above. After tleap: Thr199 has HG1,
and tleap happens to place it towards Glu106 OE1 (1.71 Å, 158°); Glu106 and Glu117 have no carboxyl hydrogen; the
twelve histidines have the hydrogens decided; Cys206 has HG; His3 has one amide hydrogen; Lys261 has OXT.
`tleap/protein_H.pdb` is the protein with these hydrogens and the deposited numbers, for the stages that need a
protonated file.

A throwaway test at 8:23 PM IST: tleap loads the three ligand files against `md/ligand/sfa_anion.mol2` as it then
was, all 18 atoms matched by name, charge −1.

## 8. Differences from the docking receptor

1. Glu106 deprotonated (docking: protonated).
2. His107 and His122 HID, His4 HID (docking: HIE for all three). The docking default left the bare ND1 of His107
   against the Glu117 carboxylate, 8 Å from the zinc.
3. Acetyl cap on His3 (docking: NH₃⁺). With point 1 the protein is −2 here and 0 in the docking receptor.
4. Ten residues take another conformer (section 3).
5. 215 crystal waters (docking: none).
6. Atoms without density built by tleap (docking: PDBFixer, then minimised).
7. Deposited heavy-atom coordinates unchanged (docking: minimised with restraints, shifts up to 0.175 Å).
8. The ligand is the anion (AutoDock4Zn docking: the neutral molecule).

Same as in the docking: the zinc histidines, His64 (HIE, conformer A), Tyr7 conformer B, Cys206, all other
histidines HIE, the removal of the other hetero groups.

## 9. What is uncertain

1. **His64.** One fixed state stands for a residue that is about 30 % protonated at pH 7.4, has two tautomers 1:1
   and two positions 1:1. In the rotated start the aniline N is 3.3 Å from the ring of the inward conformers, so
   the state of His64 can bear on the question. A repeat of the rotated start with HIP64 would test it; the plan
   has no such run.
2. **Glu106.** Deprotonated on indirect evidence (zinc in PROPKA, hydrogen-bond scheme, neutron structure of
   another sulfonamide). A hydrogen on OE2 towards water 415 is not excluded by geometry.
3. **The rotated start has an empty pocket.** Its ring replaces four waters, and the space of the crystal ring at
   the hydrophobic wall holds neither ligand nor crystal water. The solvation step may not fill a closed pocket;
   the equilibration should check this before the restraints are released.
4. **Asn67.** In the crystal its ND2 faces the third sulfanilamide, which is removed. If the amide turns over in
   solution, the distance to OD1 changes by about 0.7 Å in the rotated pose (ND2 is 3.31 Å from the aniline N) and
   1.2 Å in the crystal pose (ND2 at 6.55 Å). Neither crosses a midpoint.
5. **Minor protonation choices far from the site:** His36 (26 % protonated), His4 (tautomer mixed in solution),
   His3 and His10 (tautomer not decided by the structure), Cys206 (side chain taken from the mercurated crystal).
6. **Start geometry of N and R:** neutral bond lengths, Zn–N 2.2 Å, a constructed hydrogen. The first minimisation
   changes all of this; the starts are docked orientations, not relaxed structures.
7. **Hydrogens of water and hydroxyl groups** are placed by tleap without regard to hydrogen bonds.

Not done here, by the task: parameters of the ligand and the zinc, solvation, ions, minimisation, MD.
