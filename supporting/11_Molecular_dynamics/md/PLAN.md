# Molecular dynamics of the hCA II–sulfanilamide complex: plan

Decided by the user on 5 Oct 2026 ("Let's go with B"), after the scoping of the same day
(`paper/audit/research/md_scoping_2026-10-05.md`). This file fixes the question, the design, the gates and the rules
for reporting **before** any production run. A later change to it is made in a dated line at the end, never silently.

## 1. Scope

One system: human carbonic anhydrase II with sulfanilamide on the zinc (PDB 6RL9, 1.0 Å). No molecular dynamics of
DHPS, Sul2 or gyrase B: a fragment of 11 heavy atoms in pockets built by mobile loops, beside Mg²⁺ and a pyrophosphate,
gives runs that cannot be read cleanly in either direction (scoping, protocol report).

## 2. The question

The docking of the paper could not separate two orientations of the ring at the zinc: the crystal one and one rotated
about the zinc anchor (AutoDock4Zn ranks them within 0.000–0.219 kcal mol⁻¹; Section 3.6). The question for the
simulation is:

> In explicit water, with a zinc site that leaves the ring free, does the rotated orientation persist, or does it
> relax to the crystal orientation? And does the crystal pose itself stay where the 1.0 Å structure has it?

What the simulation cannot show, and will not be said to show: that sulfanilamide binds (the crystal structure shows
it), how strongly it binds (no MM-GBSA or MM-PBSA is computed), or anything about the Zn–N bond, which every
defensible zinc model keeps by construction.

## 3. Design

| Item | Choice |
|---|---|
| Protein | 6RL9 chain A, all crystal water kept, protonation re-audited (Glu106, His64, the three zinc histidines) |
| Ligand | the sulfonamidate anion, ArSO₂NH⁻, on the zinc (the species that binds); net site charge +1 |
| Starts | X: crystal pose (SAN 304). N: docked near-native pose. R: docked ring-rotated pose. The docked poses are of the neutral molecule; the hydrogen that points at the zinc is removed, the one that donates to Thr199 is kept |
| Runs | 3 starts × 3 independent runs (new velocities) × 100 ns = 900 ns, plus equilibration and a 5 ns validation run |
| Force field | protein ff14SB, water TIP3P, 0.15 M NaCl, ligand GAFF2 with RESP charges and repaired terms at sulfur |
| Zinc site | bonded model after the MCPB.py procedure: bonds and angles to zinc from a B3LYP/6-31G* Hessian (Seminario), RESP charges from a larger capped model; torsions through the metal are zero, so the ring is free |
| Settings | rhombic dodecahedron, 10 Å shell (about 25,600 atoms), PME, 300 K, 1 bar, 4 fs step with hydrogen mass repartitioning, staged release of restraints over about 2 ns |
| Engine | OpenMM 8.6 (already cited in the paper); parameters with AmberTools 26; analysis with MDAnalysis and ProLIF |

The run length is 100 ns and not the 50 ns of the first sketch: no docking paper of the journal with a known length
has less, and the difference is about US$4 of GPU time.

## 4. Stages and gates

Each stage ends in a gate with numbers fixed here. A failed gate stops the work and is reported; it is not tuned away.

| Stage | Work | Where | Gate |
|---|---|---|---|
| S0 | Software environment (AmberTools 26, ParmEd, OpenMM, MDAnalysis) in its own prefix; exact package list recorded | Pi | imports and versions recorded |
| S1 | Ligand: GAFF2 terms checked against the paper's own DFT (geometry of the anion and of the neutral molecule, relaxed torsion scans); angles at sulfur and the S–N and C–S torsions refitted where they miss | Pi | G-L: C–S–N within 2° and S–N within 0.02 Å of DFT at the MM minimum; torsion profile within 0.5 kcal mol⁻¹ of the scan |
| S2 | Zinc site: small model (zinc, three imidazoles, ligand) optimised with its Hessian, large capped model for the electrostatic potential, both at B3LYP/6-31G* with PySCF; Seminario and RESP through MCPB.py | Pi | G-Z: after minimisation and 100 ps, Zn–N(sulfonamide) within 0.10 Å of 1.93 Å, Zn–N(His) within 0.10 Å of 2.01–2.02 Å, N–Zn–N angles 100–120°, no fifth ligand within 2.6 Å |
| S3 | System builder: protonation, crystal water, the three starts, solvation, equilibration protocol; one script per step, seeds recorded | Pi | topology loads; energy minimises; density and temperature settle |
| S3b | Docking sensitivity: hCA II redocked with Glu106 deprotonated (the docking receptors have it protonated) | Pi | reported as found |
| S4 | Validation run, 5 ns from the crystal pose | Pi or GPU | G-V: ligand heavy-atom RMSD to the crystal pose at most 1.0 Å on average; Thr199 OG1···N hydrogen bond present in at least 80 % of frames; zinc geometry as in G-Z; τ_fold and C–S–N near the crystal values (132.1°, 108.8°) |
| S5 | Production: 9 runs of 100 ns; a 10-minute speed test first; checkpoints; checksums; results pulled to the Pi; the machine deleted | rented GPU | all runs complete; checksums verified; cost within the cap |
| S6 | Analysis, each number re-derived by a second, independent script | Pi | the two scripts agree |
| S7 | Text: a section of the supplementary material with one figure and one or two tables; a few sentences in the main text; references through `lit.py`; facts | Pi | shown to the user before it enters the manuscript |
| S8 | Verification round as before submission (independent readers), rebuild, package, sync | Pi | all gates of the project pass |

Fallbacks for S2, in this order, each stated openly in the text if used: (a) empirical zinc force constants
(MCPB.py step 2e; Yu, Li and Merz 2018) with the same RESP charges; (b) a nonbonded zinc with harmonic restraints on
the four Zn–N distances only. No restraint or collective variable is ever put on the orientation of the ring.

## 5. Rules fixed before production

1. **Validation first.** If the crystal-start runs leave the crystal orientation (ligand RMSD above 2 Å for the last
   half of the run) in two or three of three runs, the force field contradicts the 1.0 Å structure. The simulation is
   then not used as support for anything and is reported as a failed validation.
2. **Orientation coordinate.** The orientation is read from two distances fixed here: aniline N to Phe131 CZ (3.9 Å
   in the crystal pose, 9.3–9.8 Å in the rotated pose) and aniline N to Asn67 OD1 (7.7 Å against 2.6–2.7 Å). A run
   "converts" when it stays on the other side of the midpoint of both for the last 20 ns.
3. **Reading the outcome.**
   - Rotated starts convert in two or three of three runs, crystal and near-native starts stay: the rotated
     orientation is not stable outside the scoring function.
   - No run converts: reported as "both orientations persisted for 100 ns"; the simulation does not decide, and the
     text says so.
   - Mixed results are reported run by run; no average is taken over runs that disagree.
4. **Sensitivity.** If the outcome is a conversion, the rotated start is repeated once with the fallback zinc model
   (b), to show that the result does not come from the bonded terms.
5. **No end-point free energies.** MM-GBSA and MM-PBSA are not computed.
6. **Everything is reported.** All nine runs enter the supplementary material as time series, whichever way they
   fall.

## 6. Compute and cost

- Measured on the Pi (OpenMM, CPU): 4–5 ns per day for this system. The Pi does set-up, short tests and analysis.
- Production needs one rented GPU for about a day: 930 ns at about 25,600 atoms is 16–26 GPU-hours on an RTX 4090 or
  L40S class card. These speeds are published benchmarks, not measured here; a 10-minute test on the rented card
  comes first.
- Estimated cost: US$8–14 (Rs 750–1,300) on a marketplace RTX 4090; US$19–34 on the IBM L40S spot profile; US$76–134
  on IBM on demand. Tax and card mark-up are not included.
- **Nothing is rented before the user names the route and the ceiling.** Proposed ceiling: US$40 (Rs 3,900).
- Trajectories (a few GB) go to the Pi archive and R2, not into git; inputs, parameters, final frames and analysis
  tables go into git and the ZIP file of the supplementary material.

## 7. Changes to the paper (held until the results are shown)

Supplementary material: one new section (methods, validation, results), one figure, one or two tables, about ten
references. Main text: one sentence in the methods, one or two in Section 3.6, the limitation in Section 4.3 and one
sentence of the Conclusions; about 70–100 words. Cover letter, highlights if the result warrants it, README of the
submission folder, and the facts file. No text is changed before the user has seen the results.

## 8. Risks

- The zinc-site route has not been run on this machine; it is the least certain step.
- GAFF2 puts the C–S–N angle near 102–103°, the value that flips the orientation in the paper's own docking; without
  the repair of S1 the simulation would repeat the docking artefact.
- A null result (both orientations persist) is plausible.
- A molecular-dynamics study of this complex is published (Nguyen et al. 2025, cited in the Introduction); the gain
  here is the test of the docking degeneracy with a bonded zinc model and independent runs, not a first simulation.
- The receptor for the simulation (Glu106 deprotonated, crystal water kept) differs from the docking receptor; S3b
  covers the first difference.

## 9. Needed from the authors

1. The route for the GPU and the ceiling (section 6).
2. Whether the main text may grow by about 100 words (the 8,000-word limit is the project's own; the journal has
   none), or cuts must pay for it.
3. The full text of the 2017 article, if at hand (independent of this plan).

## 10. Time

Set-up and validation on the Pi: 3–5 days. Production: about 1 day once a GPU is there. Analysis, text, rebuild and
verification: 2–3 days. Planned total: 10–14 days; less if the zinc-site route works at the first attempt.

## Changes to this plan

- 5 Oct 2026: first version.
- 5 Oct 2026, 8:40 PM IST, before any run: three points made operational, none changed in substance.
  (1) Section 5, rule 2: the midpoints are 6.7 Å (aniline N to Phe131 CZ) and 5.2 Å (aniline N to Asn67 OD1), halfway
  between the crystal start and the rotated start as built (3.90 and 9.59 Å; 7.74 and 2.60 Å; `md/system/starts.csv`).
  A run has converted if at least 90 % of the frames of its last 20 ns are in the state opposite to its start, has
  stayed if at least 90 % are in the state of its start, and is mixed otherwise. (2) Section 5, rule 1: "ligand RMSD
  above 2 Å for the last half of the run" is read as the mean RMSD over the last half of the run. (3) Section 3: the
  box is a truncated octahedron with a 10 Å shell, because tleap, which writes the topology, offers no rhombic
  dodecahedron; the atom count is reported with the built systems. The full definitions are in
  `md/analysis/SPEC.md`, which is fixed before production.
- 5 Oct 2026, 9:41 PM IST, before any run: four points of the gates corrected or made operational after the analysis
  specification was written (`md/analysis/SPEC.md`, version 1.0, which raised them).
  (1) Gates G-Z and G-V, N–Zn–N angles: the window "100–120°" was wrong, because the 1.0 Å structure itself has one
  angle of 98.4° (His96 NE2–Zn–His119 ND1; the other five are 103.8–116.8°). The criterion is: each mean N–Zn–N
  angle within 10° of its value in 6RL9.
  (2) "No fifth ligand within 2.6 Å" is read as: an oxygen (ligand sulfonyl, water or protein) within 2.6 Å of the
  zinc in at most 5 % of the frames. The fraction is reported in any case.
  (3) Gate G-V, "τ_fold and C–S–N near the crystal values": the mean C4–S–N2 angle within 3° of 108.8°, and the mean
  τ_fold within 20° of 132.1° (the torsion profile of the anion is flat within 0.10 kcal mol⁻¹, so the pocket sets
  this angle and the tolerance is wide).
  (4) The 5 ns validation run is run X0 of its zinc model and is analysed on its own; it is not one of the nine
  production runs. Section 3 named ProLIF for the analysis; the analysis is the project's own code after the
  specification, checked by a second, independent program, and ProLIF is not used.
- 5 Oct 2026, 10:15 PM IST, before any run: after stage S1 (`md/LIGAND.md`).
  (1) Section 8 said that GAFF2 puts the C–S–N angle near 102–103°. Measured with GAFF2 2.2.30 and AM1-BCC charges:
  anion 107.2° in vacuum and 103.8° in implicit water (DFT 109.8°); neutral molecule 103.0° and 100.1° (DFT 109.5°).
  The earlier "96°" of the scoping came from a tool fault (`parmchk2 -s gaff2` reads GAFF 1.81 in this build; only
  `-s 2` reads GAFF2).
  (2) Gate G-L is judged in implicit water (generalised Born), like for like with the DFT reference in IEF-PCM water;
  the vacuum values are reported beside it. With the repaired terms it passes for the anion (C–S–N 111.2° against
  109.8°, S–N 1.580 Å against 1.580 Å, scan within 0.09 kcal mol⁻¹).
  (3) New gate G-R in stage S2, because G-L does not test the coordinate that separates the two orientations. They
  differ in the rotation about N–S as seen from the zinc (Zn–N2–S–C4 about −117° in the crystal pose and −35° to
  −32° in the rotated docked poses), which the force field sets through the H–N–S–C4 term, the angles at the
  nitrogen and the 1–4 terms to the zinc. Gate G-R: on 8 to 12 geometries of the small zinc model in which the
  SO₂–aryl part is turned rigidly about N–S, the force-field energies follow B3LYP/6-31G* single points within
  1.0 kcal mol⁻¹ (root mean square, relative to the lowest point), and the difference between the two orientations
  within 1.0 kcal mol⁻¹. If it fails, the H–N–S–C4 terms are refitted on the zinc model and the gate is run again;
  the refit and its residuals are reported.
- 6 Oct 2026, 1:03 AM IST, before any run: after the systems were built (`md/TOPOLOGY.md`).
  (1) Size. With the 10 Å buffer of tleap's truncated octahedron the systems hold 40,902 (X), 40,911 (N) and 40,857
  (R) atoms, not the 25,600 of section 3, which came from another program's box. The built size is kept: it leaves
  32 Å between the protein and its image (23 Å if the protein turns). The GPU time of section 6 grows by about
  60 %: roughly 26–42 GPU-hours, US$13–22 on a marketplace RTX 4090, US$30–55 on the IBM spot profile (estimates
  from the atom count, not measured). The proposed ceiling of US$40 on a marketplace card stands.
  (2) Ions: 33 ion pairs and one Na⁺, 0.15 M counted in the volume of the solvent.
  (3) Start R and water. The rotated start has no water where the ring of the crystal pose lies (the docking put
  the ring elsewhere and the solvation adds no water inside the site). An empty pocket beside the ligand could
  draw the ring back and make a conversion look easier than it is. Check added to stage S3 and to every R run: at
  the end of the held part of the equilibration (stages 1 to 3, 700 ps, solute heavy atoms held, water free) the
  number of water oxygens within 3.5 Å of the positions that the ring atoms C1 to C6 and N1 have in the crystal
  pose is counted; if it is zero, the held part is extended in steps of 500 ps, up to 2 ns more, before anything
  is released. The count is reported for every run.
- 6 Oct 2026, 3:02 AM IST, before any run and before any charge or force constant of the zinc site has been seen
  (the potential of the large model is still being computed): after the independent review of the set-up
  (`paper/audit/research/md_setup_review_2026-10-06.json`) and stage S2 as far as it has gone (`md/ZINC.md`).
  (1) **What separates the two starts.** The crystal start and the rotated start are not two rotamers about N–S.
  With the Zn–N2 bond and the hydrogen bond N2–H···OG1(Thr199) both kept, the sulfur sits on one or the other of
  the two remaining tetrahedral sites of the sulfonamide nitrogen (signed volume of Zn, H and S about N2: +0.72 in
  X, +0.79 in N, −0.71 in R). Going from R to X needs an inversion at N2 or the loss of that hydrogen bond. In a
  bonded model the inversion barrier comes from the three angle terms at N2 (Zn–N2–S, Zn–N2–H, S–N2–H), so these
  terms and the H–N2–S–C4 series decide whether a run converts.
  (2) **Gate G-R, extended.** B3LYP/6-31G* single points of the small zinc model, with the settings of S2, on:
  (a) 5 to 7 geometries along the out-of-plane coordinate of N2, from the configuration of X through a planar
  nitrogen to the other configuration; (b) the ligand as it lies in starts X and R; (c) the turn about N–S in each
  configuration, at least four points each. The force field with its final charges must follow these energies
  within 1.0 kcal mol⁻¹ (root mean square, relative to the lowest point), and give the difference between (b) X
  and (b) R within 1.0 kcal mol⁻¹. If not, the equilibrium values and constants of the three angles at N2 and the
  H–N2–S–C4 series are fitted on these points, the residuals are reported, and the fitted terms stand for all runs.
  The 120° that MCPB.py step 2e would set for Zn–N2–H is not used.
  (3) **The two bonded parameter sets.** The terms at N2 come from (2) in both. For the Zn–N bonds and the N–Zn–N
  and Zn–N–C angles there are two sets: Seminario's method on the optimised small model (MCPB.py step 2s, the
  method named in section 3) and the empirical constants (step 2e). The Seminario set is the primary one if its
  calculation is finished when production starts (about 1 to 3.5 days on the Pi, or a few hours on a rented
  machine); otherwise the empirical set is primary. The other set is the sensitivity set of rule 4; if only one
  set exists, the repeat of rule 4 uses the terms at N2 moved by the residual of the fit in (2). Which set is
  primary thus depends on computing time only, and is fixed here. Values of the modified Seminario method (step
  2ms) are tabulated beside and not used.
  (4) **Fallback (b) is dropped.** The model with a free Zn²⁺ and four distance restraints does not keep the site:
  in 4 ps tests Glu106 OE1 moved onto the zinc (1.80 and 1.90 Å; 4.00 Å in the crystal) or a sulfonyl oxygen came
  to 2.50 Å (`md/TOPOLOGY.md`, section 5.4). Rule 4 uses the second bonded set of (3). The distance from the zinc
  to Glu106 OE1 is reported by name in gates G-Z and G-V.
  (5) **Charges.** The potential is that of pose X. In the RESP fit the chemically equal atoms of the ligand are
  made equal (O1 and O2; C2 and C6; C3 and C5; H2 and H6; H3 and H5; HN11 and HN12), so that the charges do not
  favour the pose they were fitted in. The unconstrained set is reported beside.
  (6) **Gate G-L was a calibration** of three equilibrium values on the free anion with AM1-BCC charges; the scan
  criterion could not fail, because the DFT profile is flat. The binding check is made on the zinc model with the
  final charges: at the force-field minimum of the small model, C–S–N within 2°, S–N within 0.02 Å and Zn–N2–S
  within 3° of the quantum-chemical optimum of that model, or of 6RL9 (108.8°, 1.596 Å, 117.6°) if the
  optimisation has not been run.
  (7) **Corrections of what earlier sections say.** Section 3: in neither docked pose does a sulfonamide hydrogen
  donate to Thr199 (N–H···OG1 98–107°); both hydrogens were removed and one was rebuilt on N2 by the same
  tetrahedral construction as in start X. Stage S2 as run: pure d and f functions, density fitting, a medium
  integration grid for the small model, force constants from central differences of analytic gradients for the
  zinc, its four nitrogens and the sulfur only, hydrogens of the large model as built. Sections 6 and 10: the Pi
  runs the built systems of 40,900 atoms at about 2 to 3 ns per day, and stage S2 alone needs 1 to 3.5 days on
  it, so the time of section 10 holds only if a rented machine also does S2.
  (8) **Empty pocket of start R** (line of 1:03 AM, point 3): if no water has entered after the extra 2 ns, the run
  goes on as built, and the text states that an unsolvated pocket may favour a conversion. No water is placed by
  hand.
- 6 Oct 2026, 5:00 AM IST: **the GPU route is decided by the user** ("Go with spot on ibm01. Start, run and report every
  2 hours on status and eta. Monitor remaining credits. ibm01 is not used for any other projects now.").
  One IBM Cloud spot machine with one NVIDIA L40S (profile gx3-24x120x1l40s, region us-south, US$0.808 per hour by
  the IBM catalog as of 1 Oct 2026; the standard class at US$3.232 is not used). Ceiling for the whole job:
  US$60, enforced by a watchdog on the Pi that stops the machine. Credit of the account when decided: US$183.85,
  in one offer that expires on 21 Oct 2026, 1:02 PM IST. The machine also does the quantum chemistry of stage S2
  that the Pi would need days for (optimisation and force constants of the small zinc model, and the points of gate
  G-R that the Pi has not finished), so both bonded parameter sets of the entry of 3:02 AM can exist before
  production. Order on the machine: speed test, interruption test, quantum chemistry on its CPUs beside GPU tests,
  charges and gate G-R, the bonded systems with gate G-Z, validation run X0 with gate G-V, then the nine runs and
  the repeat of rule 4. A pre-empted machine is restarted by the watchdog and the runs resume from their last
  complete segment. No gate is relaxed because the machine is rented.
- 6 Oct 2026, 8:01 AM IST, before any energy of the test set of gate G-R has been read and before any run: the gate as
  built (`md/GR.md`) and two points of the specification and driver updates.
  (1) **Test set of G-R as built**: 26 geometries of the small zinc model with the zinc and the imidazoles fixed.
  Seven along the inversion at N2 from start X to an R-like placement (N2 has to move 0.25 Å on the sphere about
  the zinc, because with N2 fixed the ring of the rotated placement runs into His96); starts X, R and N as built;
  an R-like and an N-like placement of the crystal ligand, which have the crystal bond lengths (Rx, Nx); ten turns
  about N2–S, all of which the rigid protein would not allow and which are flagged for that; six points in which
  only the hydrogen on N2 is moved, added so that a fit can tell the Zn–N2–H and S–N2–H terms apart.
  (2) **Verdict of the second criterion**: X against Rx, within 1.0 kcal mol⁻¹. R as built differs from X by 0.08 Å
  in S–N and 0.30 Å in Zn–N (bond lengths of the docked neutral molecule), so X against R as built is reported
  beside and does not decide. Reference of the first criterion: the lowest quantum-chemical point; root mean square
  over the other points of the plan's set (the six hydrogen points, Rx, Nx and N enter a fit, not this number).
  (3) **Before a fit** the Zn–N2–H angle has the value of start X, 107.41°. The capping groups of the force-field
  copy of the small model are charged by the scheme "methyl"; the gate numbers are also given for two other schemes,
  and the spread between them (up to 0.5 kcal mol⁻¹ with placeholder charges) is stated as the uncertainty of the
  comparison. Equal charges on equal ligand atoms are imposed in the first stage of the RESP fit, so that the other
  charges adapt; the unconstrained set is reported beside.
  (4) **Fit, if the gate fails**: one of the two linear models of `md/GR.md`, section 6.2 (5 or 8 free parameters),
  chosen by the rule fixed there before any energy was seen, with a leave-one-out check; R and N as built are
  predicted, not fitted.
  (5) **Binding check** (entry of 3:02 AM, point 6): made on the vacuum minimum of the small model with all atoms
  free, and given also with zinc and imidazoles fixed.
  (6) **Analysis specification 1.2** (`md/analysis/SPEC.md`, fixed 7:18 AM IST): reported quantities only, no verdict,
  gate, rule or tolerance changed: the configuration at N2 per frame (starts +0.724, +0.786, −0.710), the
  N2–H···OG1 fraction before and after a first passage, the distances from the zinc to Glu106 OE1 and OE2, the
  water in the place of the crystal ring, and the equilibration frames as a separate series.
  (7) **Hydration check, corrected.** As worded at 1:03 AM the check could not act: in start R as built one water of
  the solvent lies 2.84 Å from the place of the aniline nitrogen N1, which points into the solvent, while no water
  is within 3.5 Å of the place of any ring carbon (nearest 4.15 Å). The count is therefore taken on the positions
  of the six ring carbons C1 to C6 only; the rest of the rule stands (applies when the ligand's own ring is more
  than 2.0 Å away from those positions; if the count is zero at the end of the held part, extension in steps of
  500 ps, at most four; then the run goes on and the text says what was found).
  (8) Output of the campaign: about 12 GB of trajectories (1.2 GB per run of 100 ns), not the "few GB" of section 6.
- 6 Oct 2026, 8:56 AM IST: **the machine exists** (IBM spot, one L40S, us-south-1, created 8:36 AM IST; set-up complete
  8:55 AM IST; quantum chemistry on it reproduces the Pi's test energy to 2 × 10⁻¹¹ Hartree). The same entry settles,
  before any further energy of gate G-R is read (read so far, by the reviewer only: X and R as built), what the
  independent review of the gate found (`paper/audit/research/md_gr_review_2026-10-06.json`):
  (1) **Dispersion.** The quantum-chemical points are B3LYP/6-31G* without a dispersion term, while the force field
  has one in its Lennard-Jones terms; between points the missing term differs by up to 10 kcal mol⁻¹. The D3(BJ)
  correction is added to every quantum-chemical energy of the gate (no new calculation); the numbers without it are
  reported beside.
  (2) **Two parts of the gate, read separately.** Part A, the terms at N2: the inversion path, the points in which
  only the hydrogen on N2 moves, and added points (turns of the ligand about Zn–N2, the bend of Zn–N2–S in its
  plane, the hydrogen points in the other configuration). Criterion as before, 1.0 kcal mol⁻¹ root mean square. If
  it fails, a fit is accepted only if its residual and its leave-one-out error are at most 1.0 kcal mol⁻¹, the
  gate number falls, the torsion coefficients stay below 5 kcal mol⁻¹ and the equilibrium angles within 10° of the
  structure; otherwise the terms stay as built and the gate is reported as failed, split by energy part.
  Part B, contacts: X against Rx and the turns about N2–S are almost wholly ring-to-imidazole contact and Zn···O
  electrostatics at unrelaxed geometries, which no term at N2 can repair. They are reported with the split by
  energy part and are not fitted. Added for part B: for each configuration the ligand relaxed with the zinc and the
  imidazoles fixed, in the force field and in quantum chemistry; the difference between the two relaxed
  configurations in the two descriptions is reported, and a mismatch above 3 kcal mol⁻¹ is stated in the text as a
  limit of the model. Part B is not a pass-or-fail gate.
  (3) **R and N as built** (docked geometries of the neutral molecule) are reported only; neither criterion uses them.
  (4) **Binding check** (entry of 3:02 AM, point 6): with the real charges the vacuum force-field minimum misses the
  crystal values on all three quantities; crystal against vacuum is not like for like. The check is judged against
  the quantum-chemical optimum of the small model, which the machine is computing; the result decides whether the
  ligand terms at sulfur are re-tuned on the zinc model with the tolerances of G-L.
  (5) **Starting orbitals** of the test points: the atomic guess (orbitals taken over from another geometry cost
  twice the cycles).
  (6) **Hydration check, threshold.** A water oxygen within 3.5 Å of a ring-carbon position is found in the crystal
  structure with the ring still in place, so that threshold did not show that the place of the ring holds water.
  The count is taken within 2.5 Å of the positions of the six ring carbons (a water inside the place of the ring);
  the analysis specification follows (version 1.3).
- 6 Oct 2026, 11:40 AM IST, before the validation run and any production run: after a REHEARSAL on the GPU
  (`md/test/REHEARSAL_2e.md`; the empirical zinc set with the terms at N2 as built, which gate G-R has not yet
  passed; three runs of 5 ns; a test of the pipeline, not a result).
  (1) **What the rehearsal showed.** The pipeline works end to end (builder, driver at 832–845 ns per day with all
  output written, copy back, both analysis programs in agreement on every value). The bonded zinc site holds: the
  four Zn–N distances and the six angles stay at the crystal values, no oxygen comes within 3.2 Å of the zinc,
  Glu106 OE1 stays at 4.2 Å, and the hydrogen bond to Thr199 is present throughout. The rotated start converted
  to the crystal-like pose 70 ps after the ligand was released, during the equilibration, by an inversion at N2
  with the hydrogen bond kept; the near-native start visited the rotated state for three frames. With terms at N2
  that are not validated this says nothing about the molecule; it shows that these terms decide the outcome.
  **No production run starts before gate G-R is settled.**
  (2) **A conversion can happen before time zero.** The time of a conversion is counted from the release of the
  ligand (the start of equilibration stage 4). A rotated run whose production holds no rotated frame is reported as
  "converted during the equilibration, t ps after release", with the series of the equilibration shown; it counts
  as converted in rule 3. The same holds, mirrored, for a crystal or near-native run.
  (3) **Folders.** The built systems of the two bonded sets are `bonded_2e` and `bonded_2s`; the terms at N2 are one
  file, the same for both sets. The validation run of a set is `runs/<set>/X0`.
  (4) The displacement limit of 0.50 Å in the minimisation, set for the dropped restrained model, is exceeded by
  the bonded model for the docked starts (0.501 Å in N, 0.528 Å in R), whose bonds to the zinc shorten from 2.2 to
  1.94 Å; it is reported and not used as a gate.
- 6 Oct 2026, 12:19 PM IST: **outcome of gate G-R** (`md/GR.md`; all 38 single points, with D3(BJ); read by the rules of the
  entries of 8:01 and 8:56 AM). **Part A has failed**: the force field misses the quantum-chemical energies by 3.08
  kcal mol⁻¹ root mean square against the criterion of 1.0 (2.85 without the dispersion term). No fit is accepted
  under the rule fixed beforehand (the 5-parameter model leaves 1.25 residual and 2.35 leave-one-out; the
  8-parameter model fits, but with a negative force constant and an angle 34° from the structure). The terms at N2
  therefore stay as built. The direction of the failure: along the constructed inversion path the planar nitrogen
  lies 6.6 kcal mol⁻¹ above X in quantum chemistry and 13.8 in the force field, so the force field resists the
  inversion about twice as strongly; with only the hydrogen moved into the plane the cost is 0.5 against 5.2; a
  turn of the ligand by 30° about Zn–N2 costs 4.2 against 1.8, so that rotation is too free in the force field.
  Part B (reported, not a gate): X against the R-like placement +3.4; the ten turns about N2–S 2.0 root mean
  square. **Relaxed references:** in the force field the R-like placement has no minimum of its own in the small
  model; in quantum chemistry the relaxation from the R-like placement, with zinc and imidazoles fixed, flattened
  at N2 and crossed to the other configuration at step 16, 18.8 kcal mol⁻¹ below its start. Both descriptions thus
  give one minimum for the zinc-site model. By section 4 a failed gate stops the stage and is reported: the
  decision on what follows is with the user (asked at this time). Not affected and going on: the quantum chemistry
  of the Seminario set, the relaxation from X, the binding check.
  Operational details of the gate that were fixed while it was built (`md/GR.md`, section 8): the two-body D3(BJ)
  term; part A read in two families, each against its own lowest point; twelve added points; a relaxation that
  crosses is stopped and reported as crossed; gate tasks on 7 threads.
- 6 Oct 2026, 1:09 PM IST: **correction of the entry of 12:19 PM**, after the independent check of the corrected gate
  (`paper/audit/research/md_gr_check_2026-10-06.json`), which reproduced every gate number by its own route (part A
  3.08 with D3(BJ), 2.85 without; 3.50 against a single reference; 2.20 without the inversion path; no fit
  acceptable).
  (1) **Withdrawn:** the sentence "Both descriptions thus give one minimum for the zinc-site model". The
  quantum-chemical relaxation from the R-like placement was stopped by its own rule when the configuration at N2
  changed sign, at a planar nitrogen (signed volume +0.097), unconverged, 3.5 kcal mol⁻¹ above and 3.3 Å away from
  where the relaxation from X then stood. What is shown: the R-like placement with a pyramidal nitrogen is not a
  minimum of the zinc-site model in quantum chemistry. Where the relaxation ends is not known; it is continued to
  convergence, and a run counts as crossed only when the signed volume has passed zero and the run has converged.
  (2) **What part A measures.** Nine of its 23 comparisons move more than the atoms at N2: on the inversion path
  most of the force-field energy is electrostatics between the sulfonyl group and the zinc site, and the turns
  about Zn–N2 are contacts. The failure concerns the force field's whole description of the nitrogen, fixed charges
  and 1–4 terms included, not the three angle terms alone. The statement that the force field resists the inversion
  about twice as strongly as quantum chemistry rests on the constructed path (6.6 against 13.8 kcal mol⁻¹) and on
  the hydrogen moved into the plane (0.5 against 5.2).
  (3) **Binding check** (entry of 8:56 AM, point 4): with the optimisation of the small model nearly converged, the
  vacuum force-field minimum misses the quantum-chemical optimum on two of three quantities (provisional). The
  reference itself, B3LYP/6-31G* in vacuum, misses the 1.0 Å structure on S–N by twice the tolerance, so the ligand
  terms are not re-tuned to it. The like-for-like test of the ligand geometry is the validation run in the protein
  against 6RL9 (gate G-V: C4–S–N2 within 3°, τ_fold within 20°, ligand RMSD), which stands as defined.
  (4) **Fit rule, two gaps closed for any later reading of the gate** (for example with the Seminario terms): an
  accepted fit must also bring part A to at most 1.0, otherwise it is reported as "terms improved, gate still
  failed" and not used; for the two angles at the hydrogen the 10° window is taken about the quantum-chemical
  optimum of the small model, not about the constructed hydrogen. A second reading writes its files under its own
  names and does not overwrite the record of this one.
- 6 Oct 2026, 1:34 PM IST: **the user's decision after the failed gate G-R: "3"**, that is, of the three courses put to him
  (proceed and report the failed gate; stop the simulation; build a better model of the nitrogen first), the third.
  What "better" means is fixed here, before any new energy is computed.
  (1) **Aim.** A force-field description of the sulfonamide nitrogen of the zinc-bound anion that follows
  B3LYP-D3(BJ)/6-31G* on the zinc-site model within 1.0 kcal mol⁻¹, including the energy of the R-like configuration
  against the crystal-like one, which the terms as built miss by 3.4 kcal mol⁻¹ (too high) while they overstate the
  cost of a planar nitrogen.
  (2) **What may be changed.** Only terms of the standard AMBER forms that carry N2 or act across it, with the
  zinc and histidine terms, the Lennard-Jones terms and the functional forms untouched: the three harmonic angles
  at N2 (force constants 10 to 300 kcal mol⁻¹ rad⁻², equilibrium angles free within 15° of the quantum-chemical
  optimum of the small model); Fourier torsions about N2–S (H–N2–S–C4/O and Zn–N2–S–C4/O) and about Zn–N2
  (N(His)–Zn–N2–S and N(His)–Zn–N2–H), which MCPB.py sets to zero, with periodicities 1 to 3 and each coefficient at
  most 5 kcal mol⁻¹; one improper term at N2. If these do not suffice, one further step is allowed and declared
  here: charges fitted to the electrostatic potentials of two placements (X and the R-like one) of the large model
  instead of one, with the ligand's equal-atom constraints kept, followed by the same fit. Nothing else.
  (3) **Data.** The 38 points of gate G-R and new points of the same construction (zinc and imidazoles fixed, ligand
  rigid inside) that sample the directions the terms act on (the place of the sulfur and of the hydrogen about N2,
  the turn about N2–S, the turn about Zn–N2), in both configurations and between them, within the ranges that the
  rehearsal trajectories visited widened by half, and without points whose force-field energy as built lies more
  than 20 kcal mol⁻¹ above the lowest. The new points are split into a fitting set and a test set (one in four) by a
  seed fixed before any of their energies is computed; the test set is used once, at the end.
  (4) **Acceptance, all of it** (gate G-R2): on the test set the force field follows quantum chemistry within 1.0
  kcal mol⁻¹ root mean square against ONE reference for all points (the lowest quantum-chemical point of the whole
  set), so that the energy of one configuration against the other is tested too; the largest single deviation at
  most 2.5; the 23 comparisons of part A of gate G-R, read as on 6 Oct (two families), at most 1.0, with the
  single-reference number reported; every parameter within the bounds of (2); leave-one-out error on the fitting
  set at most 1.0; and the force-field minimum of the small model within 5° of its quantum-chemical optimum in the
  three angles at N2. A model that misses any of these is not used.
  (5) **After acceptance**: the systems are rebuilt with the accepted terms, and the plan goes on unchanged
  (gate G-Z, validation run X0 with gate G-V, the nine runs). The terms as built become the second set of rule 4
  for the rotated start. The failed gate G-R and this rebuild are reported in full.
  (6) **Limit in time and money.** The work is given until 11:00 PM IST on 6 Oct 2026. If no model is accepted by
  then, it stops, and the user chooses between the two remaining courses. The ceiling of US$60 stands: the nine
  runs, the validation run and the repeat of rule 4 need about US$30 of it, so the model work may use the machine
  for about a day at most; the report warns when the projected total passes US$50.
- 6 Oct 2026, 3:51 PM IST: gate G-R, the relaxed references (`md/GR.md`, sections 6.6 and 8.2). Fixed before any step of
  the continued run existed: a relaxation is never stopped at a change of sign at N2; it counts as crossed only when
  it has converged with the sign opposite to its start; the frozen atoms are held exactly; quantum-chemical
  geometries are accepted by the force-field side within 0.01 Å on the frozen atoms; a second reading of the gate
  writes to its own folder. Results so far: the relaxation from X converged at 3:11 PM IST with the nitrogen
  pyramidal in the configuration of X, 2.29 Å from the crystal placement and 0.26 Å from the force-field minimum
  reached from X, so both descriptions agree on where the crystal-like ligand settles when only the zinc site holds
  it. The relaxation from the R-like placement, continued past the planar nitrogen, is moving towards that same
  place (signed volume +0.36 at step 24, 4.6 Å and 4.8 kcal mol⁻¹ from the end point of X) and is not converged;
  nothing is concluded from it yet. No number of the failed gate changed.
- 6 Oct 2026, 8:41 PM IST: **outcome of gate G-R2: NOT ACCEPTED** (`md/N2.md`, sections 6 and 9; `md/qm/n2/n2_gate_result.json`).
  188 new single points were computed (141 for fitting, 47 for testing; all converged). No model of the allowed
  class meets the conditions on the fitting data alone: of 30 candidates (five models, three ridge values, with the
  charges as built and with charges from two placements, the declared second step) the best leave-one-out errors
  are 1.29 and 1.18 kcal mol⁻¹ against the limit of 1.0, and the force-field minimum of the small model misses the
  quantum-chemical Zn–N2–S angle by 8° (limit 5°). Because no result on the test set could change this, the test
  set was not read, by the builder or by the independent examiner. The force field stays as built. What the
  attempt showed: the allowed terms bring the 23 comparisons of the old part A from 3.08 to 0.74 in-sample (charges
  as built), so the description of the nitrogen can be improved but not to the standard set; the misfit sits in
  points with hydrogen contacts shorter than any the rehearsal visited. By point 6 of the entry of 1:34 PM the user
  chooses between the two remaining courses.
  **Relaxed reference, now complete** (`md/qm/gr/relax/qm_Rx.json`): the quantum-chemical relaxation from the R-like
  placement, continued past the planar nitrogen with the frozen atoms held exactly, converged at step 66 in the
  minimum reached from X: 0.12 Å from that end point (heavy atoms of the ligand), 0.01 kcal mol⁻¹ from its energy,
  signed volume at N2 +0.62 against +0.63. In the zinc-site model in vacuum the R-like placement therefore has no
  minimum of its own in quantum chemistry; the force field as built gives the same. This concerns the bare zinc
  site, without Thr199, the pocket wall or water.
- 6 Oct 2026, 9:15 PM IST: **the user's decision after gate G-R2 was not accepted: "Go with option 1"**, that is, production
  with the terms at N2 as built, both failed gates (G-R and G-R2) reported in full with the direction of the
  error. Fixed here before the validation run and any production run:
  (1) **Zinc set.** The force constants of the small model were finished at 8:47 PM IST, before production starts,
  so by the rule of the entry of 3:02 AM (point 3) the Seminario set (`bonded_2s`) is the primary one, provided it
  passes gate G-Z and its validation run X0 passes gate G-V; if it fails either, that is reported and the
  empirical set (`bonded_2e`), whose rehearsal met every number of G-V, is primary, again only after its own X0.
  (2) **Runs.** Validation run X0 (5 ns) of the primary set; then X1 N1 R1, X2 N2 R2, X3 N3 R3 of 100 ns each, in
  that order; rules 1 to 6 of section 5 and the entries on how a conversion before time zero is counted stand.
  (3) **Sensitivity runs of the rotated start (rule 4), each once, 100 ns, after the nine runs:** with the other
  zinc set; and with the best candidate of the refit of the terms at N2 with the charges as built
  (`md/qm/n2/n2_terms_new_charges_as_built.frcmod`), which was not accepted by gate G-R2 and is used only to bracket:
  at the relaxed near-planar nitrogen the terms as built are 2.4 kcal mol⁻¹ too high and this candidate 1.9 to 2.7
  too low. They are reported as sensitivity runs and enter no count of rule 3.
  (4) **What the text may and may not claim** is unchanged: no statement on rates or on the height of a barrier;
  the crystal-start runs are the validation against the 1.0 Å structure; all runs are reported whichever way they
  fall.
  (5) **Before the long unattended part** the recovery after a stop is tested once on the real machine (a shutdown
  from inside while the queue is empty; the watchdog must see it and start the machine again, and the task runner
  must come back by itself).
- 6 Oct 2026, 9:22 PM IST: the recovery after a stop was tested on the real machine with the queue empty (entry of
  9:15 PM, point 5). A shutdown from inside at 9:16:00 PM was reported by IBM as "stopped" within 26 s; the watchdog
  saw it at its next run (9:20 PM), wrote the stop to the ledger and started the machine; at 9:21 PM it was running
  at the same address, and the task runner and the switch had come back by themselves. About 5 minutes in all.
  This is the path a pre-emption of the spot machine takes, and also what the switch on the machine relies on.
- 6 Oct 2026, 10:50 PM IST: **zinc set and queue.** Gate G-Z passed for both zinc sets (first 100 ps of X0). Validation
  run X0: set 2s (Seminario) FAILED gate G-V in one part (GV3 d: a sulfonyl oxygen within 2.6 Å of the zinc in 3 of
  50 all-atom frames, 6 % against the 5 % allowed; 13 of 500 solute frames; every other number met); set 2e
  (empirical) PASSED with every number; the two analysis programs agree on every value of both runs. By point 1 of
  the entry of 9:15 PM the empirical set (`bonded_2e`) is therefore primary; the rule was applied as written and not
  read again after the result. Queued at 10:14 PM IST: X1 N1 R1, X2 N2 R2, X3 N3 R3 with `bonded_2e`; then the two
  sensitivity runs of start R, with `bonded_2s` (`runs/bonded_2s/R1`) and with `bonded_2e_n2fit`
  (`runs/bonded_2e_n2fit/R1`), run id R1 in each, so with the seeds of R1 of the primary set. Details and every
  number: `md/PRODUCTION.md`. Three things said openly:
  (a) X0 of both sets was put into the queue together, before the result of set 2s was known (each simulated zinc
  model has its own validation run by `md/RUN.md`, and set 2e is simulated either way); no production task was
  queued before both gates were read.
  (b) Gate G-Z as worded is read on frames in which the zinc site is still held to the crystal positions, so it
  tests little about the angles; the free stages are tabulated beside it in `md/PRODUCTION.md`, and the test that
  counts is gate G-V.
  (c) The zinc sensitivity run uses a set that failed gate G-V. It stays in the queue as the entry of 9:15 PM fixed
  it (point 3, "with the other zinc set"); it is reported as a sensitivity run of a set that did not pass
  validation, with the failed number beside it, and enters no count of rule 3.
- 7 Oct 2026, 12:10 AM IST: **after an independent check of the production set-up**, made while X1 was running and
  before any run from start N or R had begun (four readers: the gates, the queue, the systems, the first output;
  one sceptic per reader; record `paper/audit/research/md_production_check_2026-10-07.json`). No finding stops or
  changes a run; the queue is as it was queued. Fixed or corrected here:
  (1) **How the time of a conversion is taken (completes the entry of 6 Oct, 11:40 AM, point 2; written before R1
  starts).** The two analysis programs count the first passage from the start of the production only
  (specification 1.3, section 4), so the time "from the release" is derived as follows, for every run alike. The
  solute frames of equilibration stages 4, 5 and 6 (400 + 400 + 500 ps; a frame every 10 ps; extension steps of
  stage 3 come before the release and do not count) and the solute frames of the production are joined in order
  of time. On the joined series the assigned states and the visits are formed by the rules of section 4 of the
  specification, unchanged (the state before the first frame is the state of the start; a long visit has at least
  100 frames, 1 ns). The time of the conversion is the time of the first frame of the first long visit in the
  opposite state, counted from the start of stage 4. If that frame lies in the equilibration the run is reported
  as "converted during the equilibration, t ps after release"; if it lies in the production both times are given
  (in the production, and from the release = that time + 1.3 ns). A run without such a visit has no time. The
  number is taken from the `state` of each frame as written by both programs (`timeseries_equil_solute.csv`,
  `timeseries_solute.csv` of `analysis/` and of `analysis_recheck/`), and the two must give the same frame. The
  verdict of a run and rule 3 are untouched (they read the last 20 ns), and the time describes one run; it is not
  a rate (entry of 6 Oct, 9:15 PM, point 4).
  (2) **Water in the systems (corrects "all crystal water kept" of section 3 and completes the entry of 6 Oct,
  1:03 AM, point 3).** The systems hold 215 of the 216 water sites of chain A of 6RL9 (site 510, occupancy 0.32,
  dropped in all starts). Start R holds 211: waters 423, 516, 553 and 568 lie 0.28, 1.43, 0.55 and 0.19 Å from
  heavy atoms of the rotated ligand and were removed when the start was built, by the 2.2 Å rule of `md/SYSTEM.md`
  section 4 (`md/system/decisions.json`, `starts/R/waters_removed`). So in the runs from X and N the place of the
  rotated ring holds four ordered waters at the start, and in the runs from R the ligand lies there and the
  place of the crystal ring is empty until solvent fills it (the hydration check). These are the four waters of the
  docking test `dock/WATER4.md`. Nothing is rebuilt.
  (3) **What gate G-V tests for set 2e.** In set 2e the equilibrium values of the four Zn–N bonds and of the six
  N–Zn–N angles are those of 6RL9 itself (`md/ZINC.md`, section 8.1), so parts a to c of GV3 show only that the
  rest of the force field does not move the site away from them; they are not an independent reproduction of the
  zinc geometry. The parts of the gate that are not put in are GV1, GV2, GV3 d and GV4. "Passed with every number"
  in the entry of 6 Oct, 10:50 PM, is to be read with this, and the text of the paper must say it.
  (4) **How narrow the failure of set 2s is.** The verdict stands as the rule was worded before the run, and so
  does the choice of the primary set. For the report: set 2s missed the limit by one frame of 50 (3 against the 2
  allowed), while the ten times denser series gives 13 of 500 (2.6 %), and of the ten ways to take every tenth
  frame of that series two give 3 or more. What is solid is the mean distance of the nearer sulfonyl oxygen from
  the zinc: 2.835 Å with 2s, 2.973 Å with 2e, 2.986 Å in 6RL9 (third reading, in the record above). Set 2s is not
  to be described as five-coordinate or as an invalid model. In the bonded model the sulfonyl oxygens and the zinc
  are three bonds apart, so their interaction is the scaled 1–4 one; this is said where the count is reported.
  (5) **Gate G-Z.** Which frames it is read on was not laid down before the run. The reading used (the mean of the
  ten solute frames of stage 1, with part d) was written into `md/production/gate_GZ.py` when X0 of set 2s had
  started and before its files were analysed; those frames are held to the crystal positions. Both sets pass
  under five other readings tried by the checker, so nothing hangs on the choice, and it is reported as a check on
  held frames.
  (6) **Temperature.** The kinetic temperature of the log averages 301.0 K in every production segment at the
  4 fs step (300.2 to 300.6 K in the stages at 2 fs), with the thermostat at 300 K; it is steady, the same in both
  zinc sets and in the rehearsal, and density and geometry are steady. The pattern is that of the kinetic-energy
  estimate at a long step; the cause was not established by a test. The runs go on as set, and both numbers go
  into the methods.
  (7) **Correction of the entry of 6 Oct, 11:40 AM, point 1:** "no oxygen comes within 3.2 Å of the zinc" holds for
  water and protein oxygens; the nearer sulfonyl oxygen is within 3.2 Å in nearly every frame, as in the crystal.
  (8) **Steady departures from the crystal that no gate tests** (set 2e; X0 and the first 20 ns of X1): the zinc
  lies 0.5 Å out of the plane of the His119 ring (0.17 Å in 6RL9); S–C4 is 1.824 Å (1.771 Å); S–N2–Zn 114.0°
  (117.55°); O1–S–C4 101.6° (106.6°). They are named in the text, which does not call the zinc site or the
  sulfonyl group as a whole reproduced.
  (9) **Reading the campaign.** No analysis folder is written into a production run folder before the queue is
  empty and the last pull is complete (a look in between goes to a scratch folder with `--out`); the comparison of
  the two programs over the campaign is run after that pull. The campaign table prints an outcome line for the
  two sets with one run each (`bonded_2s`, `bonded_2e_n2fit`); these lines are not outcomes of rule 3 and are not
  quoted (entry of 6 Oct, 9:15 PM, point 3).
  (10) **While the queue runs:** nothing is pushed to the machine and no system is rebuilt (the driver compares
  the build record of a system when a run is continued); the file `KEEP_RUNNING` is not touched (with tasks in the
  queue it is not needed, and at the end it would hold back the stop).
- 7 Oct 2026, 1:00 AM IST: **the guard's check at the end of the work: a fix is prepared and reviewed, and NOT
  installed.** The production check (queue part, finding F1, confirmed by its sceptic) found that the final
  verification compares, among all pulled files, three files that the machine itself rewrites all the time
  (`status/heartbeat.json`, `status/queue.json`, `status/deadman.json`), so that a try can fail by timing alone
  with every result identical on both sides. A change of `final_check` in `ops/pi/cqc_gpu_watchdog.py` that
  leaves these three out of the comparison and names them in the record was written on a copy, tested against
  the faked machine (the existing scenarios unchanged, 216 of 216; a new scenario 37 of 37, 16 of 37 on the file
  as installed) and accepted by an independent reviewer without a change to the code
  (`paper/audit/research/md_production_check_2026-10-07/guard_fix_proposed_not_installed/`). The installation was
  declined by the session's permission check, because it changes what a guard verifies; it is therefore left to
  the user. The installed guard is the file of 6 Oct, 11:41 AM IST (sha256 `99bf04b9ab19…`), unchanged.
  **What follows without the fix:** no result is at risk. At the end a final try may fail although nothing
  differs; the tries repeat every 5 minutes; after 3 hours the guard stops the machine all the same, with
  `ALERT_DONE_UNVERIFIED` (about US$2.50 of idle time), and `cqc_gpu.py delete` then refuses. In that case the
  record `final_verification.json` is read first: if `pi_against_machine.differ` names only these three files and
  the driver's checksums of every run file are right, the results are the same on both sides, and the deletion of
  the machine is put to the user.
- 7 Oct 2026, 7:25 AM IST: **the fix of the guard's final check is installed** (entry of 1:00 AM). The user ran the
  install script on the Pi at 7:16:48 AM IST; the installed file is the reviewed one (sha256 `15e51515cc6e…`). The
  new end-of-work scenario passes on the installed file in the faked world (37 of 37), and the first cron run of
  the new file, at 7:20 AM IST, logged as the runs before it. The changed lines act only in the final try at the
  end of the work; that will be their first run on the real machine. Record: `ops/pi/GPU_WATCH.md`, section 14.
- 8 Oct 2026, 4:45 AM IST: **stage S5 is complete.** The nine runs of 100 ns and the two sensitivity runs ended in
  the order queued, each at its first attempt; no task failed, none was repeated, and IBM did not pre-empt the
  machine once (the one restart in the ledger is the test of 6 Oct, 9:16 PM). The last run ended at 3:45 AM IST. The
  guard's final pull and verification passed at the first try at 4:12 AM IST: 3,128 files have on the Pi the sha256
  they have on the machine (the machine's three running-state files named as not compared), and 1,284 files of
  333 finished steps match the driver's own records (`/root/cqc-sfa-md/status/final_verification.json`, verified
  true). The guard then stopped the machine; it is not deleted. Cost at the stop: US$37.86 of the US$60 cap
  (43.60 h of running). Gate of S5: all runs complete, checksums verified, cost within the cap: met.
  Stage S6 (analysis by the two programs, with a blind third reading and a physical review beside them) started at
  4:42 AM IST. The stopped machine (US$0.03 per hour) is kept until the two programs agree on every run, so that
  a run could be repeated if a file proved unusable; its deletion is then put to the user.
- 8 Oct 2026, 7:00 AM IST: **stage S6 is complete; its gate is met.** The two analysis programs agree on every value
  of all runs and on the campaign rules (the documented command, exit status 0; 6,721,008 values compared, no
  disagreement). A blind third reading from the raw frames agrees with the record in every number, and an auditor's
  own reader agrees with all three. Outcome by the rules as worded: rule 1 is not failed (0 of 3 crystal-start runs
  left the crystal orientation); rule 3 gives its first case (R1, R2, R3 converted; X1 to X3 and N1 to N3 stayed).
  All three conversions lie in the equilibration after the release, with the protein backbone still held. The
  sensitivity runs, in no count: converted with zinc set 2s; MIXED with the candidate terms at N2. The record with
  every number, the limits and the corrections after two audits and a physical review is `md/RESULTS.md`
  (`md/results/`; agents' returns in `paper/audit/research/md_campaign_check_2026-10-08/`). No run has to be
  repeated. Next: the results are put to the user; no manuscript text before that (stage S7).
- 8 Oct 2026, 8:30 AM IST: **the rented machine is deleted, by the user's word** ("Delete the vm and related
  resources. Ensure nothing is billing in ibm01."). `cqc_gpu.py delete --yes` at 8:20 AM IST removed the instance
  `cqc-md-gpu-01` with its boot volume, the floating IP, the subnet, the VPC and the key (5 of 5; the verification
  record said verified, nothing left on the machine). Read back from IBM afterwards: no resource named `cqc-md-*`
  is left in the account (lists of the region and the account-wide search). The two guard lines were taken out of
  the Pi's cron table (nothing else in it changed), and the session's two-hourly report job was cancelled.
  Cost: US$37.97 by the project's ledger (43.60 h of running, the boot volume and address for 47.7 h, 14.24 GB
  pulled), against the cap of US$60; IBM's usage record for October names 43.46 instance-hours of the GPU profile,
  US$35.12. Stage S5 is closed. What still bills in the account belongs to the other project there (its stopped
  machine `pulsed-p00` keeps a 250 GB boot volume and a floating IP); this project did not touch it, and it is
  put to the user.
- 11 Oct 2026, 11:50 AM IST: **stage S7 is complete: the results are in the paper.** After the user's word of
  8 Oct 2026 ("All approved. Prepare the full manuscript ready to submit in the journal with all recent results.")
  the supplementary material received Section S15 (the two docking tests: Glu106 as the carboxylate, four ordered
  waters kept) and Section S16 (these simulations: methods, the tests of the force field with the two that failed,
  every run, the limits), Tables S36 to S40 and Figs. S9 and S10, and the main text was changed where the results
  bear on it. What the text may and may not say is fixed in
  `paper/audit/research/md_paper_2026-10-08/DESIGN.md` (D4). The documents were read by two rounds of four
  independent readers, each round followed by a repair (`paper/audit/FINAL.md`, section 12). Rule 6 of section 5 is
  met by Table S39, Figs. S9 and S10 and the time series of every run in `md/results/series/`.
  Corrections to this record found on the way; the entries above stay as written:
  (a) entry of 7 Oct, 12:10 AM, point 8: the zinc lies 0.15 Å out of the plane of the His119 ring in 6RL9 (0.154 Å
  from the prepared structure); the 0.17 Å printed there is the value for His94. The paper prints 0.15 Å.
  (b) same entry: "300.2 to 300.6 K in the stages at 2 fs" are the two stage means of the validation run X0 of
  set 2e; over the 27 stage means at 2 fs of all 13 runs the range is 299.8 to 300.9 K, which the paper prints.
  (c) Section S16 says in plain words what the entries of 6 Oct (1:34 PM and 9:15 PM) record: by section 4 a failed
  gate was to stop the work, and going on with the terms as built after the second failure was a departure from
  that rule, decided by the user and written down here before the validation runs and the production.
  (d) Section S16 also says that the rule "a change during the equilibration counts as a conversion" was laid down
  after the rehearsal of 6 Oct, in which the run from start R had changed in that way (entry of 11:40 AM), and
  that the definition of its time was completed while X1 was running (entry of 7 Oct, 12:10 AM).
  The data of the campaign on the Pi (`/root/cqc-sfa-md/`, 14.5 GiB in 6,773 files, the trajectories among them)
  were copied to the project's storage (`md-archive/` under the R2 prefix of the project) on 11 Oct at 10:42 AM
  IST and compared by size, with no difference; the paper offers the trajectories on request. Stage S8 (submission
  folder, package, GitHub, R2) follows.
