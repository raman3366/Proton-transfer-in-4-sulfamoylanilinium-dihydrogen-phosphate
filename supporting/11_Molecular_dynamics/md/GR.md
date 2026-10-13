# Gate G-R: the force-field terms at the sulfonamide nitrogen, tested against quantum chemistry

Stage S2 of `md/PLAN.md`. First version: 6 Oct 2026, 3:05 to 5:05 AM IST. Second version: 9:00 to 11:30 AM IST, after
the independent review of the gate (`paper/audit/research/md_gr_review_2026-10-06.json`) and the plan entry of 8:56 AM
IST; it is in commit e27c380. **This version: 6 Oct 2026, 1:15 to 3:55 PM IST**, after the independent check of the
corrected gate (`paper/audit/research/md_gr_check_2026-10-06.json`) and the plan entry of 1:09 PM IST, which it
implements: what a crossing of a relaxation is, the amended rule for accepting a fit, a second reading of the gate
under its own names, and the record brought up to date. Changed: sections 1, 4.4, 6.3, 6.4, 6.6, 7, 8, 9 and 10, and
one sentence each in sections 3.4 and 6.2. The test set, the 38 single points and every number of parts A and B are
those of 11:30 AM (compared file by file: no number differs). `md/PLAN.md`, `paper/` and `dock/` were not changed by
this work, and this version is not committed. Every number is in a file named in brackets; files without a folder are
in `md/qm/gr/`.

> **Note added on 11 Oct 2026 (nothing above or below was rewritten).** This version describes the state of 6 Oct
> 2026, 3:50 PM IST, when the quantum-chemical relaxation from Rx had made 23 steps. The run finished later that
> day: it **converged at step 66 (6 Oct 2026, 6:49:52 PM IST)** in the minimum reached from X. The finished record is
> in git since 11 Oct 2026 as `relax/qm_Rx_converged_step66.json` (`relax/README.md` says which file is which, with
> checksums); `relax/qm_Rx.json`, `relax_qm.json` and `gr_gate.json` are kept as they stood at 3:48 PM IST. The
> places of this file that call the run unconverged carry a dated note: sections 1, 6.6, 9 and 10. The plan has the
> result in its entry "Relaxed reference, now complete" (`md/PLAN.md`).

## 1. The gate in plain words

**What is compared.** The 64-atom small model of the zinc site (zinc, three ethyl-imidazoles, the sulfonamidate
anion; charge +1), in vacuum, with the zinc and the imidazoles fixed, at 38 geometries of the ligand. For each
geometry two energies: quantum chemistry, B3LYP/6-31G* plus the D3(BJ) dispersion energy, and the force field of the
production runs on the same 64 atoms. If the two follow each other, the force field's description of the sulfonamide
nitrogen N2 can be trusted with the question of the simulation: whether the ring-rotated start, which differs from
the crystal pose by the configuration at N2, converts.

**The gate is read in two parts** (plan, 8:56 AM, point 2).

- **Part A. This is the gate.** 25 geometries made from X or from Rx by moving the atoms at N2: the inner points of
  the inversion path, the hydrogen moved alone (in both configurations), the bend of Zn–N2–S, and rigid turns of the
  ligand about Zn–N2. Criterion: the force field follows quantum chemistry within 1.0 kcal mol⁻¹ (root mean square
  over 23 comparisons, in two families; the number against a single reference is quoted beside it and does not
  decide). If it fails, the three angle terms at N2 and the H–N2–S–C4 torsion may be fitted, and the fit is accepted
  only under the rule of section 6.3; otherwise the terms stay as built and the gate is reported as failed.
- **What part A measures: the force field's whole description of the nitrogen, fixed charges and 1–4 terms included,
  not the three angle terms alone** (check of 6 Oct 2026, finding 7; plan, 1:09 PM, point 2). In nine of the 23
  comparisons more moves than the atoms at N2: the five inner points of the path, the three turns about Zn–N2 that
  are not the reference, and X itself, which is compared against the turned point e_p120, the reference of its
  family. Along the path the sulfonyl group and the ring move with the nitrogen, and most of the force-field energy
  is non-bonded energy between the ligand and the fixed atoms: at the planar point a4 the force field is
  13.76 kcal mol⁻¹ above X, of which 7.61 are non-bonded energy between ligand and fixed atoms
  (electrostatic 8.00) and 4.39 the terms a fit could change; at a6, 11.56, 8.35 and 1.14. In the
  four turns about Zn–N2 the fitted terms do not change at all (0.00); the whole difference is contact and
  electrostatics with the fixed atoms. And where only the hydrogen moves, the fitted terms are still not the whole
  of it: with the hydrogen in the plane (d1) the force field costs 5.20 against 0.46 in quantum chemistry, and
  5.17 of the 5.20 are non-bonded energy inside the ligand (the fixed charge and the 1–4 terms of that
  hydrogen), 0.33 the fitted terms [`gr_gate.json`, `part_A_composition`]. **The failure does not rest on the nine
  comparisons**: the same number is 2.20 without the path (18 comparisons), 2.71 without the turns about Zn–N2
  (19), 2.70 on the hydrogen points and bends alone (14), and 3.17 with X and Rx themselves as the
  references (23) [`gr_gate.json`, `verdict_part_A.subsets`; each subset read like part A, against the lowest
  quantum-chemical point of each family within the subset]. So a failed part A says that this force field (angle
  terms, torsion, fixed charges, 1–4 terms) does not describe the nitrogen as quantum chemistry does. It does not say
  that the three angle terms are the fault, and the statement that the force field resists the inversion about twice
  as strongly as quantum chemistry rests on the constructed path (6.6 against 13.8) and on the hydrogen moved
  into the plane (0.5 against 5.2).
- **Part B, contacts. Reported, not a gate, never fitted.** The crystal placement against the ring-rotated one (X
  against Rx), the turns about N2–S, and the two placements after relaxation of the ligand. These energies are
  almost wholly contact between the ring and the imidazoles and electrostatics between the sulfonyl oxygens and the
  zinc, at geometries that are not relaxed; no term at N2 can repair them. They are given with the force-field
  energy split by part. For the relaxed placements a mismatch above 3 kcal mol⁻¹ between the two descriptions is to
  be stated in the text as a limit of the model. A quantum-chemical relaxation counts as crossed to the other
  configuration only when it has converged there (section 6.6).
- **Reported only:** starts R and N as built. They are docked geometries of the neutral molecule (S–N 0.08 Å and
  Zn–N 0.3 Å too long); their mismatch is bond and angle strain and tests nothing at N2.
- **Beside the gate, the binding check:** at the force-field minimum of the small model (all atoms free) C4–S–N2,
  S–N2 and Zn–N2–S are compared with the quantum-chemical optimum of the same model (tolerances 2°, 0.02 Å, 3°). Its
  result is reported; by the plan entry of 1:09 PM, point 3, the ligand terms are not re-tuned to that optimum
  (section 4.4).

**Where it stands.**

As of 3:50 PM IST on 6 Oct 2026. **38 of 38 single points are finished** (3 on the Pi, 35 on the rented machine). All
numbers in kcal mol⁻¹, force field minus quantum chemistry [`gr_gate.json`, `gr_table.csv`; written 3:48 PM IST].

| | with D3(BJ) | without D3 |
|---|---|---|
| **Part A**, r.m.s. over 23 of 23 comparisons, two families (criterion 1.0) | **3.08** | 2.85 |
| Part A, against one reference for all points (does not decide) | 3.50 | 4.01 |
| Part A, largest deviation | +6.62 at a4 | +6.55 at a4 |
| Part A, without the points flagged "protein" (14 comparisons) | 2.65 | 2.82 |
| Part A, without the inversion path (18 comparisons) | 2.20 | |
| Part A, without the turns about Zn–N2 (19) | 2.71 | |
| Part A, hydrogen points and bends only (14) | 2.70 | |
| Part A, X and Rx themselves as the references (23) | 3.17 | |
| Part B, X against Rx: quantum chemistry / force field / difference | +10.84 / +14.22 / **+3.38** | +19.93 / +14.22 / −5.71 |
| Part B, the 10 turns about N2–S, r.m.s. of the difference | 1.97 | 6.53 |
| Part B, Nx against X | +1.20 | +0.65 |
| Reported only, R as built against X | +19.48 | +13.30 |
| Reported only, N as built against X | +13.33 | +13.77 |

- **Part A has failed with the terms as built**: 3.08 in two families (3.50 against a single reference), criterion
  1.0. References: e_p120 in family X, g3 in family Rx. Under the other two cap schemes the number is 3.08 (3.56) and
  2.93 (3.25).
- The same number for the two other force-field models that are evaluated beside: 3.19 (3.53 against a single
  reference) with Zn–N2–H as step 2e wrote it (120°), 3.24 (3.81) with the charges fitted without the ligand
  equivalences [`gr_asis_gate.json`, `gr_free_gate.json`].
- Deviations of part A, each against the reference of its family: X itself: X −0.51, Rx −2.72; path (a): a4 +6.62, a3
  +4.34, a5 +6.57, a2 +1.25, a6 +4.76; hydrogen on X (d): d1 +4.23, d3 +2.83, d4 −1.03, d5 −1.62, d6 +1.36, d2 +2.14;
  turns about Zn–N2 (e): e_p015 −1.93, e_p030 −2.94, e_m120 −1.26; bend of Zn–N2–S (f): f_p08 +0.04, f_m08 +0.22;
  hydrogen on Rx (g): g1 +0.01, g2 −2.84, g4 −2.15, g5 −4.03, g6 −0.90.
- **Part B, X against Rx.** The force field puts Rx 14.22 above X, of which 11.63 are electrostatic energy between
  ligand and fixed atoms, 1.00 Lennard-Jones between them, 0.70 non-bonded inside the ligand, 1.35 the three angles at
  N2, 0.95 other angles and −1.41 the torsions about N2–S. Quantum chemistry puts it 10.84 above X with D3 and 19.93
  without: the dispersion term alone is worth −9.09, and the sign of the mismatch depends on it.
- **The fit: no model is accepted, so no parameter file is written and the terms at N2 stay as built; gate G-R is
  reported as failed.** The outcome is the same under the rule as it stood at 11 AM and under the rule as amended at
  1:09 PM (section 6.3; the fitted parameters are identical, only the list of conditions has grown). Part A before the
  fit 3.08 (3.50 against a single reference). M1: residual r.m.s. 1.25, leave-one-out 2.35, part A 1.33 (1.60) after;
  equilibrium angles 126.9, 113.4, 116.4°, c₁ −2.00, c₂ 2.25; conditions not met: residual r.m.s. at most 1.0,
  leave-one-out r.m.s. at most 1.0, part A with the fitted terms at most 1.0. M2: residual r.m.s. 0.63, leave-one-out
  0.94, part A 1.13 (1.39) after; equilibrium angles 105.4, 141.7, 117.6°, force constants −6.6, 7.6, 46.7, c₁ −1.69,
  c₂ 1.97; conditions not met: force constants in range, standard errors small, equilibrium angles within 10° of their
  reference, part A with the fitted terms at most 1.0 [`fit/fit.json`, with the split of every deviation by energy
  part; outcome "no fit accepted"].
- **What follows is not decided here.** By the plan a failed gate stops the stage and is reported, and it is not tuned
  away (section 4 of the plan; entry of 8:56 AM, point 2; the decision is with the user, plan entry of 12:19 PM). The
  terms at N2 are those of `mm_2e_main.json`, `fit_terms_in_topology`; nothing was changed after the numbers were
  seen.
- **Part B, relaxed references** (zinc and imidazoles fixed; section 6.6). Force field: the relaxation from Rx leaves
  the ring-rotated configuration (signed volume at N2 −0.68 → +0.78) and ends in the minimum reached from X; there is
  no second relaxed configuration [`relax_mm.json`]. Quantum chemistry, as of 3:50 PM IST [`relax_qm.json`,
  `relax/qm_X.json`, `relax/qm_Rx.json`]:
  - from X: **converged** after 55 steps (3:11 PM IST), 13.44 kcal mol⁻¹ below its first step; signed volume at N2
    +0.724 → +0.629, angles at N2 104.1, 125.4 and 109.8° (Zn–N2–S, Zn–N2–H, S–N2–H; sum 339.3°), Zn–N2–S–C4 −148.8°,
    heavy atoms of the ligand 2.29 Å from the start; pyramidal, configuration of the start; crossed by the rule: no.
  - from Rx: **not converged**, 23 steps so far (last step 3:45 PM IST), 19.30 kcal mol⁻¹ below its first step, r.m.s.
    gradient on the ligand 7.4e-04 Hartree bohr⁻¹; signed volume at N2 −0.678 → +0.273, angles at N2 117.9, 128.6 and
    110.1° (Zn–N2–S, Zn–N2–H, S–N2–H; sum 356.5°), Zn–N2–S–C4 −80.1°, heavy atoms of the ligand 1.44 Å from the start.
  - **Note added on 11 Oct 2026: the run from Rx has converged** [`relax/qm_Rx_converged_step66.json`, `final`]. 66
    steps (6 Oct 2026, 6:49:52 PM IST), 24.30 kcal mol⁻¹ below its first step, r.m.s. gradient on the ligand 1.3e-04
    Hartree bohr⁻¹; signed volume at N2 −0.678 → +0.619 (largest on the way +0.774), angles at N2 104.6, 125.5 and
    109.9° (Zn–N2–S, Zn–N2–H, S–N2–H; sum 340.0°), Zn–N2–S–C4 −147.8°, heavy atoms of the ligand 5.26 Å from the
    start; pyramidal, the configuration opposite to the start; crossed by the rule: yes. Against the end point of the
    run from X: 0.12 Å (heavy atoms of the ligand, r.m.s., no fitting) and E(from Rx) − E(from X) = −0.01 kcal mol⁻¹
    with D3(BJ), −0.05 without (from `e_total` and `e_scf` of the two records). In the bare zinc-site model the
    R-like placement therefore has no minimum of its own in quantum chemistry either. The three bullets that follow
    (the distance and the two energy differences "between the two geometries as they are", and the queue status)
    describe step 23 and are kept as written; they are not results.
  - Between the two geometries as they are (NOT both converged, so this is the distance between two runs on their way,
    not a result): 4.70 Å (heavy atoms of the ligand, r.m.s., no fitting), largest distance of one atom 8.22 Å; E(from
    Rx) − E(from X) = +4.98 kcal mol⁻¹ with D3(BJ), +6.38 without [`relax_qm.json`, `end_points`].
  - The force field at the same two geometries (single points of the cut-down model, as section 4): E(from Rx) −
    E(from X) = +5.48 kcal mol⁻¹, of which electrostatic energy between ligand and fixed atoms +8.98, Lennard-Jones
    between them −0.34, the three angles at N2 −0.94, other angles −2.03, torsions about N2–S +1.10, other torsions
    +0.86, bonds −0.71, non-bonded inside the ligand −1.44 [`relax_mm.json`,
    `Rx_minus_X.force_field_at_the_quantum_chemical_geometries`]; **not a result while a run is unconverged**.
  - On the machine (queue status of 3:47 PM IST): `gr_relax_X` done at 3:11 PM; `gr_relax_Rx_cont`, the continuation
    from Rx, running on 7 CPUs since 3:11 PM (attempt 1), beside the single points of gate G-R2 (`n2_points`) on 16.
    Nothing else of gate G-R is running; `gr_binding_final` waits for `qm_write`.
  - **Where the run from X ended.** 2.29 Å (heavy atoms) from the crystal placement, with the ring turned by 32° about
    N2–S (Zn–N2–S–C4 −148.8° against −117.0°) and by 11° about Zn–N2, and the nitrogen still pyramidal in the
    configuration of X (signed volume 0.63, with the hydrogen opened to Zn–N2–H 125.4°). That is 0.26 Å from the
    force-field minimum reached from X: on where the crystal-like ligand settles when only the zinc site holds it, the
    two descriptions agree. They differ in the geometry at N2 and at sulfur there (Zn–N2–H 125.4° against 106.9°, S–N2
    1.634 against 1.571 Å, C4–S–N2 106.4° against 111.3°), which is why the force field puts the quantum-chemical end
    point 19.2 kcal mol⁻¹ above its own minimum; they are the differences the binding check shows with all atoms free
    (section 4.4).
- **Binding check: failed (provisional)** [`binding_check.json`, written 3:48 PM IST; section 4.4]. Force field
  against the quantum-chemical optimum of the small model: C4–S–N2 111.34 against 107.10° (+4.24, tolerance 2), S–N2
  1.572 against 1.635 Å (−0.063, tolerance 0.02), Zn–N2–S 109.91 against 107.74° (+2.17, within 3). The optimum itself
  lies +0.039 Å in S–N2 and −9.86° in Zn–N2–S from the 1.0 Å structure, both outside the same tolerances, so by the
  plan the ligand terms are not re-tuned to it; the like-for-like test is the validation run in the protein (gate
  G-V).

**What the small model in vacuum cannot tell about the protein.**

1. It has no Thr199, no Thr200, no Leu198, no water. The hydrogen bond N2–H···OG1(Thr199) and the wall of the
   pocket are in neither energy. The gate compares two descriptions of the same 64 atoms; it does not give the
   barrier of the inversion in the protein, nor which configuration the protein prefers.
2. Left to itself, the ligand of the small model does not stay where the protein holds it, in either description.
   In the force field, with zinc and imidazoles fixed, the minimum reached from X lies 2.25 Å (heavy atoms, r.m.s.)
   from the crystal placement, and **the ring-rotated placement has no minimum of its own**: from Rx a descent that
   only goes downhill, in steps of at most 0.02 Å, passes through a planar nitrogen (15.3 kcal mol⁻¹ below the start)
   and ends in the minimum reached from X [`relax_mm.json`]. In quantum chemistry the relaxation from X has moved the
   ligand 2.29 Å from the crystal placement and has converged, and the pyramidal R-like placement is not a minimum
   either: the relaxation from Rx ran downhill at every one of its first 16 steps to a planar nitrogen (section 6.6).
   So a "difference between the two relaxed configurations" may not exist in this model at all; what keeps the
   rotated placement in the protein, if anything does, is outside the model.
3. Of the 38 geometries, 19 would run into protein atoms that the model does not contain (flag "protein" in sections
   3.3 and 3.4). They test the two descriptions against each other; they are not places the complex visits.
4. Two errors of the reference are not quantified: the basis-set superposition error at 6-31G*, which acts against
   the dispersion term, and the use of charges fitted on the large model in the small one (section 9).

**Which coordinates are still untested.**

- **The direction and length of the Zn–N2 bond against the three histidines** (the N(His)–Zn–N2 angles and Zn–N2
  itself). They change only along the inversion path, together with everything else, and in R and N as built. Their
  terms are those of the zinc site (step 2e or Seminario) and are judged by gate G-Z, not here.
- **The rotation about Zn–N2 and the bend of Zn–N2–S in the R-like configuration.** Both were added from X only
  (points e and f). Until 9:20 AM the rotation about Zn–N2 was not tested at all.
- **Nothing is relaxed in part A.** The path is one constructed path (section 3.1); the lowest path of the model,
  and so the barrier of the model, is not known.

**Decisions, and when they were taken.**

| When (IST) | Where | What |
|---|---|---|
| 5 Oct, 10:15 PM | plan | Gate G-R is created: the force field follows B3LYP/6-31G* within 1.0 kcal mol⁻¹ r.m.s. on rigid turns about N–S |
| 6 Oct, 3:02 AM | plan, points 2 and 6 | The starts differ by the configuration at N2, not by a turn. Test set: inversion path, X and R, turns in both configurations. If it fails: fit of the three angles at N2 and of H–N2–S–C4. The 120° of step 2e for Zn–N2–H is not used. Binding check on the zinc model |
| 6 Oct, 8:01 AM | plan, points 1 to 5 | The 26 points as built; second criterion on X against Rx; reference = lowest quantum-chemical point; Zn–N2–H at 107.41° before a fit; cap scheme "methyl"; fit models M1 and M2 |
| 6 Oct, 8:56 AM | plan, points 1 to 5 | **D3(BJ) is added** to every quantum-chemical energy, numbers without it beside. **Two parts**, A decides, B is reported and never fitted. **Rule for accepting a fit.** R and N as built reported only. Binding check judged against the quantum-chemical optimum. Atomic guess for the orbitals. Read by then: the energies of X and R |
| 6 Oct, 9:00 to 9:30 AM | here, section 8 | How the entry of 8:56 AM is made operational: the two-body D3(BJ) term; the families of part A; the 12 added points; the fit with one constant per family; the relaxations stop when the configuration crosses. Read by then: X, R, Rx and one r.m.s. that contained a4. Entered in the plan at 12:19 PM, after the numbers; the check found no sign that anything was tuned to a result (its finding 10) |
| 6 Oct, 12:19 PM | plan | **Outcome: part A has failed, no fit is accepted, the terms at N2 stay as built.** What follows is with the user |
| 6 Oct, 1:09 PM | plan, points 1 to 4 | After the independent check: "one minimum in both descriptions" is **withdrawn**; a relaxation counts as crossed only when the signed volume has passed zero **and** the run has converged, and the run from Rx is continued; what part A measures; the binding check is reported and the ligand terms are **not** re-tuned to the optimum; **fit rule amended** for any later reading (part A with the fitted terms at most 1.0; window of the two hydrogen angles about the quantum-chemical optimum); a second reading writes under its own names |
| 6 Oct, 1:15 to 3:50 PM | here, section 8 | How the entry of 1:09 PM is made operational: what "has passed zero and has converged" is in the program; the frozen atoms held exactly in the continued run; the tolerance of 0.01 Å for the fixed atoms of a quantum-chemical geometry; where the optimum for the window is read from; the folder of a second reading. Read by then: everything of the entry of 12:19 PM; no step of the continued relaxation |

## 2. How it is made

| Step | Command (from the repository root; `MD` = `/root/cqc-sfa-tools/md/bin/python`, `QM` = `/root/cqc-sfa-tools/env/bin/python`) | Writes |
|---|---|---|
| 1 | `MD md/qm/gr_points.py` | `points.json`, `points.csv`, `scan_ns_contacts.csv`, `scan_znn_contacts.csv`, `geom/<id>.xyz` and `.pdb` |
| 2 | `QM md/qm/gr_queue.py d3 CONFIG` (`d3-extract` once) | `d3.json`: D3(BJ) of every geometry; `d3_reference.json` |
| 3 | `MD md/qm/gr_mm.py run --params DIR --work WORK --out OUT.json --label TEXT` | force-field energies of all points with their split, three cap schemes, the two minima |
| 4 | `MD md/qm/gr_mm.py relax …`, `MD md/qm/gr_mm.py binding …` | `relax_mm.json`, `relax/mm_*.xyz`; `binding_check.json` |
| 5 | `MD md/qm/resp_equiv.py WORK` (called by `md/qm/after_qm.sh`) | constrained and unconstrained charges |
| 6 | `MD md/qm/gr_fit.py gate` / `fit` / `selftest` | `gr_table.csv`, `gr_gate.json`; `fit/`; `selftest/fit_selftest.json` |
| 7 | `QM md/qm/gr_queue.py run CONFIG` (the queue: steps 2 to 6 after every single point) | `qm/<id>.json`, everything above; bulky files in `/root/cqc-sfa-md/qm/gr/` |
| 8 | `QM md/qm/gr_queue.py relax CONFIG X` or `Rx` | `relax/qm_<start>.json`, `.xyz`; `relax_qm.json` |
| 9 | `QM md/qm/gr_queue.py evaluate CONFIG SRC` (on the Pi) | takes over results of the rented machine from SRC, then steps 2 to 6 |

`CONFIG` is `queue_prod.json` on the Pi and `queue_vm.json` on the rented machine (7 threads, nothing waited for;
otherwise the same). No step uses random numbers except the self-test of the fit (seed 20261006). A second run of
step 1 gives an identical `points.json`. Software: NumPy, OpenMM 8.6.1 (Reference platform), ParmEd,
AmberTools 26 (MD environment); PySCF 2.14.0, geomeTRIC 1.1.1 (quantum-chemistry environment); on the Pi also dftd3
1.6.0 and pyscf-dispersion 1.5.0, used as checks only. The versions are written into every result file. Nothing was
installed, on either machine.

## 3. The test set [`points.json`, `points.csv`]

Sections 3.1 to 3.3 describe the 26 points of 5 AM; section 3.4 the 12 points added at 9:20 AM. No geometry of
the first 26 was changed by the addition (`gr_points.py` checks every file against `SHA256SUMS`).

**Frame and fixed part.** The small model of MCPB.py (`md/qm/mcpb/CA2_small.pdb`, frame of 6RL9): three
ethyl-imidazoles and the zinc, 46 atoms, identical in all points; the 18 ligand atoms change. Thr199 OG1 is not in
the model; its position (17.496, −1.747, 17.667 Å) is read from the structure.

**Names.** z, u, h: unit vectors from N2 to Zn, S and H. φ: out-of-plane angle of N2–S, the angle between u and the
plane through z and h (+49.35° in X, 0 for a planar nitrogen, −44.81° in the R-like placement). τ: torsion
Zn–N2–S–C4. Body B: the 17 ligand atoms other than HN21 with the internal geometry of the crystal pose. Construction
T of the hydrogen: the one used for all three starts (`md/SYSTEM.md`, section 5; N–H 1.021 Å, of the two free
tetrahedral directions the one nearer to OG1).

### 3.1 Definitions

- **(b) X, R, N**: the ligand as it lies in the built starts (`md/system/ligand_X/R/N.pdb`). X is the geometry
  MCPB.py wrote (identical, checked). R and N are docked poses of the neutral molecule: S–N 1.677 Å, Zn–N2 2.235
  and 2.207 Å; against X, R also has S–O longer by 0.029 and 0.036 Å, O–S–N smaller by 4.6 and 4.0°, C3–C4–S–N2
  142.7° against 132.1°, and the amino group turned (C2–C1–N1–HN11 13.9° in R and −99.2° in N against 180.0° in X;
  N–H 1.010 against 1.030 Å) [review of 6 Oct 2026, findings 2 and 11]. Both are **reported only** (section 6.2).
- **(b) Rx, Nx ("R-like", "N-like")**: body B fitted on the 11 heavy atoms of the start (least squares, rotation and
  translation; r.m.s. 0.125 and 0.175 Å), then moved along the line Zn→N2 until Zn–N2 is 1.931 Å; hydrogen by T.
  Only the rigid-body placement of the crystal ligand differs from X. The same recipe applied to X returns X (largest
  shift 0.0004 Å, the rounding of the file), so the "X-like" point is X itself.
- **(a) a1 … a7, the inversion**, from X (= a1) to Rx (= a7). Zn–N2, N2–S, N2–H and B are kept. One parameter s:
  z(s) and h(s) are interpolated on great circles between their directions in X and in Rx (7.5° and 29.8° apart; N2
  moves 0.25 Å on the sphere about the zinc; the hydrogen keeps pointing to OG1: H···OG1 1.80–1.87 Å, N–H···O
  157–176°); φ, the in-plane direction ψ of N2–S, and τ run linearly from their values in X (49.35°, −135.24°,
  −117.03°) to those in Rx (−44.81°, −141.70°, −21.34°). The seven points are φ = φ_X·(1, ⅔, ⅓), 0, φ_Rx·(⅓, ⅔, 1),
  with s = (φ_X − φ)/(φ_X − φ_Rx). a4 has the planar nitrogen (angle sum 360.0°).
- **(c) c1_…, c2_…, the turn about N2–S**: B without N2 turned rigidly about the axis N2→S by δ, from X (c1) and from
  Rx (c2); N2, H and the zinc stay. A turn is "possible" when no ligand atom is closer to an imidazole atom than
  2.8 Å (heavy–heavy), 2.1 Å (heavy–hydrogen) or 1.9 Å (hydrogen–hydrogen). Taken: the five possible multiples of
  30° nearest to the reference, filled up with multiples of 15°. Result: c1 δ = −120, −90, −60, −30, +30°; c2 δ =
  −180, −150, −120, −30, −15°. The whole turn in 5° steps is in `scan_ns_contacts.csv`. Possible ranges: c1 from
  −150° to +45°; c2 from −180° to −95°, from −50° to +5°, and a window at +130° to +135°. The turn of R (τ = −32°)
  cannot be reached in the configuration of X: at δ = +85° the ring stands in His94 (C5···CE1 1.69 Å)
  [`points.json`, `turn_ranges_5deg_scan`].
- **(d) d1 … d6, only the hydrogen moved**, heavy atoms as in X. d1: H in the plane Zn–N2–S (planar nitrogen by
  moving H). d2: H on the other side of that plane. d3, d4: H turned about N2→S until Zn–N2–H is 12° larger,
  smaller (S–N2–H kept). d5, d6: H turned about N2→Zn until S–N2–H is 12° larger, smaller (Zn–N2–H kept). Added
  because without d3 to d6 a fit cannot tell the Zn–N2–H term from the S–N2–H term (section 6.5). In d1 and d2 the
  hydrogen does not point to OG1.

**Order of the queue**: X, R, Rx, a4, a3, a5, a2, a6, d1, d3, d4, d5, d6, d2, Nx, then the (c) points with the
smallest turns first and the two configurations in turn, N, then the added points of section 3.4. (Until 9 AM IST
on 6 Oct 2026, d2 and Nx came after the turns; review finding 11.) **Sets**: section 6.2.

### 3.2 Geometry at N2 (Å, degrees)

| # | Point | Zn–N2 | Zn–N2–S | Zn–N2–H | S–N2–H | Sum | φ | Zn–N2–S–C4 | H–N2–S–C4 | H···OG1, N–H···O |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | X | 1.931 | 117.6 | 107.4 | 107.4 | 332.4 | +49.4 | −117.0 | 121.8 | 1.87, 164 |
| 2 | R | 2.235 | 119.5 | 106.9 | 106.9 | 333.3 | −47.9 | −32.1 | 89.3 | 1.72, 140 |
| 3 | Rx | 1.931 | 123.8 | 105.8 | 105.8 | 335.4 | −44.8 | −21.3 | 100.6 | 1.83, 157 |
| 4 | a4 | 1.931 | 138.6 | 108.2 | 113.2 | 360.0 | 0.0 | −66.9 | 113.1 | 1.80, 174 |
| 5 | a3 | 1.931 | 135.0 | 108.3 | 113.1 | 356.5 | +16.5 | −83.6 | 120.8 | 1.81, 176 |
| 6 | a5 | 1.931 | 137.4 | 107.7 | 111.8 | 357.0 | −14.9 | −51.7 | 105.3 | 1.80, 169 |
| 7 | a2 | 1.931 | 127.4 | 108.1 | 111.2 | 346.7 | +32.9 | −100.3 | 123.9 | 1.84, 170 |
| 8 | a6 | 1.931 | 132.1 | 106.9 | 109.3 | 348.3 | −29.9 | −36.5 | 100.6 | 1.81, 163 |
| 9 | d1 | 1.931 | 117.6 | 121.2 | 121.2 | 360.0 | 0.0 | −117.0 | 63.0 | 2.29, 114 |
| 10 | d3 | 1.931 | 117.6 | 119.4 | 107.4 | 344.4 | +40.5 | −117.0 | 104.9 | 1.89, 157 |
| 11 | d4 | 1.931 | 117.6 | 95.4 | 107.4 | 320.4 | +54.8 | −117.0 | 137.0 | 1.94, 148 |
| 12 | d5 | 1.931 | 117.6 | 107.4 | 119.4 | 344.4 | +36.3 | −117.0 | 110.0 | 1.85, 170 |
| 13 | d6 | 1.931 | 117.6 | 107.4 | 95.4 | 320.4 | +58.5 | −117.0 | 130.1 | 1.98, 143 |
| 14 | c1_m030 | 1.931 | 117.6 | 107.4 | 107.4 | 332.4 | +49.4 | −147.0 | 91.8 | 1.87, 164 |
| 15 | c2_m015 | 1.931 | 123.8 | 105.8 | 105.8 | 335.4 | −44.8 | −36.3 | 85.6 | 1.83, 157 |
| 16 | c1_p030 | 1.931 | 117.6 | 107.4 | 107.4 | 332.4 | +49.4 | −87.0 | 151.8 | 1.87, 164 |
| 17 | c2_m030 | 1.931 | 123.8 | 105.8 | 105.8 | 335.4 | −44.8 | −51.3 | 70.6 | 1.83, 157 |
| 18 | c1_m060 | 1.931 | 117.6 | 107.4 | 107.4 | 332.4 | +49.4 | −177.0 | 61.8 | 1.87, 164 |
| 19 | c2_m120 | 1.931 | 123.8 | 105.8 | 105.8 | 335.4 | −44.8 | −141.3 | −19.4 | 1.83, 157 |
| 20 | c1_m090 | 1.931 | 117.6 | 107.4 | 107.4 | 332.4 | +49.4 | 153.0 | 31.8 | 1.87, 164 |
| 21 | c2_m150 | 1.931 | 123.8 | 105.8 | 105.8 | 335.4 | −44.8 | −171.3 | −49.4 | 1.83, 157 |
| 22 | c1_m120 | 1.931 | 117.6 | 107.4 | 107.4 | 332.4 | +49.4 | 123.0 | 1.8 | 1.87, 164 |
| 23 | c2_m180 | 1.931 | 123.8 | 105.8 | 105.8 | 335.4 | −44.8 | 158.7 | −79.4 | 1.83, 157 |
| 24 | d2 | 1.931 | 117.6 | 107.4 | 107.4 | 332.4 | −49.4 | −117.0 | 4.1 | 3.22, 61 |
| 25 | N | 2.207 | 105.9 | 110.4 | 110.3 | 326.6 | +56.9 | −107.5 | 133.1 | 1.62, 149 |
| 26 | Nx | 1.931 | 117.1 | 107.5 | 107.5 | 332.2 | +49.7 | −105.8 | 133.1 | 1.79, 156 |

N(His)–Zn–N2 for His94, His96, His119: X 108.9, 114.3, 116.8°; R 116.3, 112.6, 110.5°; Rx 114.2, 116.5, 109.3°;
N 117.3, 111.6, 110.3°; Nx 114.7, 115.0, 110.1°. The planar point a4 is not a symmetric one: the hydrogen stays on
OG1, so Zn–N2–S opens to 138.6° while Zn–N2–H stays at 108°. d1 is the other planar nitrogen, with both hydrogen
angles at 121.2°.

### 3.3 Closest contacts (Å)

In the model: ligand atom to imidazole atom. Protein: ligand heavy atom to a protein heavy atom that is not in the
model (`md/qm/mcpb/CA2_H.pdb`), split into N/O···N/O (hydrogen bonds are 2.5–3.0 Å) and pairs with a carbon or the
sulfur. Flag "protein": below 2.4 Å (N/O···N/O) or 3.0 Å (others), which the rigid protein of 6RL9 would not allow.
No point is below 2.0 Å between heavy atoms in the model; none was dropped.

| Point | heavy–heavy | heavy–H | H–H | Zn···O1, O2 | protein N/O···N/O | protein, other heavy | flag |
|---|---|---|---|---|---|---|---|
| X | 3.21 N2···His94 NE2 | 2.62 O2···His94 HE1 | 2.72 H3···His94 HE1 | 4.05, 2.99 | 2.86 N2···Thr199 OG1 | 3.17 C5···Thr200 OG1 | |
| R | 3.09 C5···His94 NE2 | 2.79 H5···His96 CE1 | 2.12 H5···His96 HE1 | 4.53, 3.91 | 2.58 N2···Thr199 OG1 | 3.19 C2···Gln92 NE2 | |
| Rx | 2.88 C3···His94 CE1 | 2.72 C3···His94 HE1 | 2.05 H5···His96 HE1 | 4.20, 3.81 | 2.44 N1···Asn67 OD1 | 3.12 C2···Gln92 NE2 | |
| a4 | 3.26 N2···His94 NE2 | 3.06 O2···His94 HE1 | 2.95 H3···His94 HE1 | 4.41, 3.67 | 2.51 O1···Thr199 N | 2.52 C5···Thr200 OG1 | protein |
| a3 | 3.24 N2···His94 NE2 | 3.02 O2···His94 HE1 | 2.88 H3···His94 HE1 | 4.37, 3.51 | 2.52 O1···Thr199 N | 2.62 C5···Thr200 OG1 | protein |
| a5 | 3.27 N2···His119 ND1 | 3.04 O2···His94 HE1 | 2.95 H3···His94 HE1 | 4.40, 3.75 | 2.61 O1···Thr199 N | 2.62 C5···Thr200 OG1 | protein |
| a2 | 3.23 N2···His94 NE2 | 2.86 O2···His94 HE1 | 2.78 H3···His94 HE1 | 4.26, 3.27 | 2.65 O1···Thr199 N | 2.82 C5···Thr200 OG1 | protein |
| a6 | 3.25 N2···His119 ND1 | 2.99 C3···His94 HE1 | 2.47 H5···His96 HE1 | 4.33, 3.79 | 2.76 N1···Asn62 ND2 | 3.02 C5···Thr200 OG1 | |
| d1 to d6 | 3.21 N2···His94 NE2 | 2.62 O2···His94 HE1 | 2.72 H3···His94 HE1 | 4.05, 2.99 | 2.86 N2···Thr199 OG1 | 3.17 C5···Thr200 OG1 | |
| c1_m030 | 2.95 O2···His94 CE1 | 2.17 O2···His94 HE1 | 3.39 HN21···His119 HE1 | 3.77, 3.11 | 2.86 N2···Thr199 OG1 | 1.81 N1···Leu198 CD1 | protein |
| c2_m015 | 3.19 C3···His94 CE1 | 2.72 O2···His94 HE1 | 1.94 H5···His96 HE1 | 4.27, 3.66 | 2.79 N2···Thr199 OG1 | 2.95 C5···Thr200 CG2 | protein |
| c1_p030 | 3.19 O2···His119 ND1 | 2.27 C3···His94 HE1 | 2.56 H3···His94 HE1 | 4.23, 3.07 | 1.86 N1···Gln92 NE2 | 2.29 C1···Gln92 NE2 | protein |
| c2_m030 | 3.12 O2···His94 CE1 | 2.45 O2···His94 HE1 | 2.17 H5···His96 HE1 | 4.32, 3.51 | 2.79 N2···Thr199 OG1 | 2.26 C6···Thr200 CG2 | protein |
| c1_m060 | 3.02 O2···His94 CE1 | 2.23 O2···His94 HE1 | 3.39 HN21···His119 HE1 | 3.44, 3.38 | 2.86 N2···Thr199 OG1 | 0.68 C1···Leu198 CB | protein |
| c2_m120 | 3.13 O2···His94 NE2 | 3.17 O2···His94 HE1 | 3.15 HN21···His96 HE1 | 3.92, 3.23 | 0.98 N1···Thr200 O | 1.19 N1···Thr200 C | protein |
| c1_m090 | 3.21 N2···His94 NE2 | 2.74 O2···His94 HE1 | 3.39 HN21···His119 HE1 | 3.15, 3.71 | 1.11 N1···Leu198 N | 0.38 C1···Leu198 CA | protein |
| c2_m150 | 3.22 N2···His119 ND1 | 2.62 O1···His94 HE1 | 3.15 HN21···His96 HE1 | 3.62, 3.46 | 2.21 N1···Leu198 O | 0.91 C6···Leu198 CA | protein |
| c1_m120 | 3.21 N2···His94 NE2 | 2.83 O1···His94 HE1 | 2.91 H5···His119 HE1 | 3.01, 4.01 | 0.67 N1···Ser197 O | 1.05 C1···Ser197 O | protein |
| c2_m180 | 2.81 O1···His94 CE1 | 2.25 O1···His94 HE1 | 3.15 HN21···His96 HE1 | 3.35, 3.75 | 1.60 N1···Leu198 N | 0.90 C2···Leu198 CB | protein |
| N | 2.93 O2···His119 ND1 | 2.74 O2···His119 HB3 | 3.23 H3···His94 HE1 | 4.24, 2.78 | 2.55 N2···Thr199 OG1 | 3.19 O1···Trp209 CZ2 | |
| Nx | 3.07 O2···His119 ND1 | 2.79 O2···His119 HB3 | 2.77 H3···His94 HE1 | 4.13, 2.98 | 2.76 N2···Thr199 OG1 | 3.33 O1···Leu198 CA | |

Two things to read from the last columns.

- **Every turn of 30° or more runs the ring through the pocket wall** (Leu198, Thr200, Gln92, Ser197). In the rigid
  protein only turns from −15° to +15° (c1) and from −10° to 0° (c2) stay clear of the flag [`points.json`,
  `turn_ranges_5deg_scan`]. The (c) points test the torsion term of the force field against quantum chemistry over a
  wide range; they are not geometries the complex visits. They are flagged, not dropped, because the model does not
  contain these protein atoms and the comparison is valid there. They belong to part B of the gate (section 6.2).
- **The middle of the inversion path presses the ring on Thr200 OG1** (2.52 Å at the planar point, 3.17 Å in X). In
  the protein the inversion needs room there; the small model does not see this.

### 3.4 Points added after the review (9:20 AM IST, 6 Oct 2026; before any of their energies existed)

Asked for by the plan entry of 8:56 AM, point 2, and by finding 6 of the review: the first set never turned the
ligand about Zn–N2 on its own, varied Zn–N2–S only along the path together with everything else, and moved the
hydrogen only with the heavy atoms of X.

- **(e) e_p015, e_p030, e_m120, e_p120: the whole ligand of X, hydrogen included, turned rigidly about the axis
  Zn→N2** by +15, +30, −120 and +120° (sense: that of the torsion NE2(His94)–Zn–N2–S, 74.2° in X). Bond lengths,
  the three angles at N2, the N(His)–Zn–N2 angles and every torsion inside the ligand stay; what changes are the
  torsions through the zinc, **which the force field sets to zero**, and the contacts. This is the test of the zero
  torsions. Rule for the choice: the possible multiples of 15° (limits of section 3.1) with |δ| ≤ 30°, and +120° and
  −120°, the two other staggered places against the three histidines. Possible multiples of 15°: +15, +30, +105,
  +120, +135, −120, −135 [`points.json`, `added_points`; the whole turn in 5° steps: `scan_znn_contacts.csv`]. In
  e_m120 and e_p120 the hydrogen bond to OG1 is given up, and the protein would not allow either (below).
- **(f) f_p08, f_m08: Zn–N2–S 8° larger and smaller than in X**, in its plane: body B without N2 turned about the
  normal of the plane Zn–N2–S through N2; the hydrogen put where Zn–N2–H and S–N2–H keep their values of X, on the
  same side of that plane. Of the three angles at N2 only Zn–N2–S changes; Zn–N2–S–C4 is kept.
- **(g) g1 … g6: the six hydrogen moves of (d), made on Rx** instead of X: the hydrogen terms in the other
  configuration.

| # | Point | Zn–N2–S | Zn–N2–H | S–N2–H | Sum | Signed volume at N2 | N94–Zn–N2–S | H–N2–S–C4 | H···OG1, N–H···O |
|---|---|---|---|---|---|---|---|---|---|
| 27 | e_p015 | 117.6 | 107.4 | 107.4 | 332.4 | +0.72 | 89.2 | 121.8 | 1.85, 172 |
| 28 | e_p030 | 117.6 | 107.4 | 107.4 | 332.4 | +0.72 | 104.2 | 121.8 | 1.92, 151 |
| 29 | e_m120 | 117.6 | 107.4 | 107.4 | 332.4 | +0.72 | −45.8 | 121.8 | 3.51, 44 |
| 30 | e_p120 | 117.6 | 107.4 | 107.4 | 332.4 | +0.72 | −165.8 | 121.8 | 3.25, 59 |
| 31 | f_p08 | 125.6 | 107.4 | 107.4 | 340.4 | +0.62 | 74.2 | 115.4 | 1.85, 172 |
| 32 | f_m08 | 109.6 | 107.4 | 107.4 | 324.4 | +0.81 | 74.2 | 126.6 | 1.89, 157 |
| 33 | g1 | 123.8 | 118.1 | 118.1 | 360.0 | 0.00 | −11.5 | 158.7 | 2.12, 121 |
| 34 | g2 | 123.8 | 105.8 | 105.8 | 335.4 | +0.68 | −11.5 | −143.3 | 3.06, 65 |
| 35 | g3 | 123.8 | 117.8 | 105.8 | 347.4 | −0.51 | −11.5 | 119.2 | 1.80, 163 |
| 36 | g4 | 123.8 | 93.8 | 105.8 | 323.4 | −0.77 | −11.5 | 84.4 | 1.96, 137 |
| 37 | g5 | 123.8 | 105.8 | 117.8 | 347.4 | −0.51 | −11.5 | 114.9 | 1.77, 175 |
| 38 | g6 | 123.8 | 105.8 | 93.8 | 323.4 | −0.77 | −11.5 | 90.5 | 1.98, 135 |

The signed volume at N2 is the triple product of the unit vectors from N2 to Zn, H and S (`md/analysis/SPEC.md`,
section 9A): +0.724 in X, +0.785 in N, −0.710 in R, −0.678 in Rx, 0 for a planar nitrogen.

Closest contacts, as in section 3.3 (Å):

| Point | heavy–heavy | heavy–H | H–H | Zn···O1, O2 | protein N/O···N/O | protein, other heavy | flag |
|---|---|---|---|---|---|---|---|
| e_p015 | 3.09 O2···His119 ND1 | 2.97 O2···His119 HB3 | 3.21 HN21···His96 HE1 | 4.05, 2.99 | 2.60 O1···Thr199 N | 2.92 N1···Phe131 CZ | protein |
| e_p030 | 2.87 O2···His119 ND1 | 3.21 O2···His119 HB3 | 3.03 HN21···His96 HE1 | 4.05, 2.99 | 2.27 O1···Thr199 N | 2.54 C1···Leu198 CD2 | protein |
| e_m120 | 2.97 O2···His94 NE2 | 2.90 O2···His96 HE1 | 2.86 H3···His96 HE1 | 4.05, 2.99 | 1.94 N1···Thr199 O | 0.61 C6···Thr200 N | protein |
| e_p120 | 3.18 O2···His119 CE1 | 2.36 O2···His119 HE1 | 3.31 HN21···His96 HE1 | 4.05, 2.99 | 1.19 O2···Thr199 OG1 | 1.29 C1···Ser197 O | protein |
| f_p08 | 3.21 N2···His94 NE2 | 2.80 O2···His94 HE1 | 3.09 H3···His94 HE1 | 4.14, 3.20 | 2.65 O1···Thr199 N | 2.81 C5···Thr200 OG1 | protein |
| f_m08 | 3.08 O2···His119 ND1 | 2.48 O2···His94 HE1 | 2.34 H3···His94 HE1 | 3.94, 2.76 | 2.86 N2···Thr199 OG1 | 3.30 C3···Val121 CG2 | |
| g1 to g6 | 2.88 C3···His94 CE1 | 2.72 C3···His94 HE1 | 2.05 H5···His96 HE1 | 4.20, 3.81 | 2.44 N1···Asn67 OD1 | 3.12 C2···Gln92 NE2 | |

The two large turns about Zn–N2 put the ring through Thr199 and Thr200 (0.6 and 1.2 Å): like the turns about
N2–S they are geometries of the model, not of the complex. In the rigid protein the rotation about Zn–N2 has room
for less than 15°: e_p015 is itself flagged (N1 2.92 Å from Phe131 CZ, limit 3.0 Å).

## 4. The force-field model of the small model [`md/qm/gr_mm.py`]

### 4.1 How it is built

MCPB.py writes no topology for its small model. tleap builds one from the same templates and parameter files as the
dry bonded complex of `md/qm/after_qm.sh`: the mol2 files of the ligand and the zinc unchanged; each histidine cut to
its side chain, CA and HA, with N and C replaced by the capping hydrogens of the quantum-chemical model (atom type
H1 on the CX carbon; ff14SB has all terms); the four bonds to the zinc; `leaprc.protein.ff14SB`, `leaprc.gaff2`, the
two ligand frcmod files, `CA2_mcpbpy.frcmod`, then any further frcmod. The atom order of the topology is that of the
quantum-chemical input (checked). Energies: OpenMM, Reference platform, no cut-off, vacuum.

**Checked against the dry complex** (on the folder of the zinc agent's known-content test, step 2e): 58 atoms are
common, 6 are caps; atom types and Lennard-Jones parameters identical for all 58; 61 bonds, 102 angles and 170
torsions among common atoms compared, none missing, one different: Zn–N2–H, 107.41° here against the 120.00° of
step 2e, by intent (below). With the cap scheme "group" the charges of all 58 common atoms are identical too. The
same check passes on a folder made with step 2s [`selftest/mm_placeholder.json`,
`selftest/after_qm_with_equivalence.json`].

**Zn–N2–H.** The plan says the 120° of step 2e is not used. The model evaluated first ("main") takes the value of
start X, 107.41°, with the force constant of the frcmod file; the model as MCPB.py wrote it ("asis") is evaluated
beside. This matters for the numbers before a fit and for the terms that stand if no fit is accepted.

### 4.2 What cancels and what does not

Only energies relative to point X are used and the 46 fixed atoms never move. **Cancels exactly:** every bonded and
non-bonded term among fixed atoms: the imidazoles, the caps (so the bonded parameters of the capping hydrogens do
not matter), the Zn–N(His) bonds, the N(His)–Zn–N(His) angles. **Does not cancel:** the ligand's own terms where its
geometry changes (in R and N as built: all of them), the three angles at N2 and the H–N2–S–X torsions, Zn–N2 (only R
and N), the three N(His)–Zn–N2 angles (path, Rx, Nx, R, N), and all non-bonded terms between ligand and fixed
atoms: imidazoles, zinc (including the 1–4 pairs through the zinc), and the caps.

### 4.3 The cap test: the comparison does depend on the caps, by up to 0.5 kcal mol⁻¹

The RESP charges come from the large model, which has the backbone. In the small model a methyl group stands in its
place, and the charge of the removed atoms N, H, C, O (−0.11 to −0.12 e per histidine) has to go somewhere. Three
schemes, all keeping the charge of each residue and of the model (+1):

| Scheme | Cap charges, His94 (e) | |
|---|---|---|
| methyl | carbon −0.271, the three hydrogens +0.087 each | a methyl group as the quantum-chemical model has one; **used** |
| group | H at the place of N −0.146, H at the place of C +0.029, CA +0.019, HA +0.087 | each capping hydrogen carries the charge of the group it replaces; nearest to the real system |
| spread | capping hydrogens −0.008, the rest spread over the 15 atoms | an extreme |

Largest difference of a relative energy over the 26 points, placeholder charges [`selftest/mm_placeholder.json`,
`cap_test`]: **0.49 kcal mol⁻¹ methyl against group** (r.m.s. 0.26), 0.50 methyl against spread, 0.48 group against
spread. The difference is systematic: −0.20 to +0.14 for the points in the configuration of X, −0.37 to −0.49 for
Rx and the c2 points. **The energy of X against the R-like point changes by 0.37 kcal mol⁻¹ with the scheme**, a
third of the tolerance of the gate. No cap atom is nearer than 4.8 Å to a ligand atom in any point [`points.json`,
`closest_cap_atom_over_all_points`]; the effect is the dipole of the cap acting on an anion with large partial
charges.

I could not show that the comparison is independent of the caps; it is not. What is done instead: `gr_fit.py`
writes the number of part A for all three schemes, and the verdict uses "methyl", because the quantum-chemical
model has methyl groups there. **With the real charges** [`mm_2e_main.json`, `cap_test`]: largest difference of a
relative energy 0.44 kcal mol⁻¹ (methyl against group, at R as built; r.m.s. 0.22 over the 38 points). A cleaner
way exists and is not done (section 9, point 2).

### 4.4 Binding check (plan, 3:02 AM, point 6; 8:56 AM, point 4; 1:09 PM, point 3) [`binding_check.json`]

`gr_mm.py run` minimises the small model with zinc and imidazoles fixed and with all atoms free. `gr_mm.py binding`
compares C4–S–N2, S–N2 and Zn–N2–S at the minimum with all atoms free with the **quantum-chemical optimum of the same
model** (tolerances 2°, 0.02 Å, 3°). It reads `md/qm/mcpb/CA2_small_opt.fchk` (on the Pi:
`/root/cqc-sfa-md/from_gpu/repo_md/qm/mcpb/`), written by task `qm_write` of the rented machine after the force
constants; before that, the geometry of the finished optimisation (`/root/cqc-sfa-md/qm/prod/opt_final.npy`, which the
watchdog brings to the Pi), and the verdict is then marked provisional. "Provisional" says only that the fchk file
does not exist yet: the geometry is the end of the optimisation (task `qm_opt`, done 12:46 PM IST). Until one of the
two exists the verdict is "pending". The queue runs the check after every point and writes the verdict into its log
and status; task `gr_binding_final` of the machine writes the final one after `qm_write`. `qm_write` waits for
`qm_hess`, which the work on gate G-R2 cancelled at 2:35 PM IST to give its single points the CPUs and which has to be
given new attempts (`queue.py retry qm_hess`); until then the verdict stays provisional. The geometry will be the
same.

**Result, as of 3:48 PM IST: failed (provisional), on two of three quantities.**

| | Force field, all atoms free | Quantum-chemical optimum | Difference | Tolerance | |
|---|---|---|---|---|---|
| C4–S–N2 | 111.34° | 107.10° | +4.24° | 2° | outside |
| S–N2 | 1.572 Å | 1.635 Å | −0.063 Å | 0.02 Å | outside |
| Zn–N2–S | 109.91° | 107.74° | +2.17° | 3° | within |

With all atoms free the ligand of the force-field minimum lies 1.12 Å (heavy atoms, r.m.s.) from start X; with zinc
and imidazoles fixed the three values are 111.30°, 1.571 Å and 107.50° [`mm_2e_main.json`].

**How far the reference itself is from the 1.0 Å structure** [`binding_check.json`, `optimum_against_crystal`]. The
optimum is B3LYP/6-31G* in vacuum with all atoms free and **without a dispersion term** (the single points of the gate
carry D3(BJ); this geometry was optimised without it), and it has no Thr199. Against 6RL9 (108.8°, 1.596 Å, 117.6°) it
lies −1.70° in C4–S–N2 (within 2°), **+0.039 Å in S–N2 (about twice the tolerance of 0.02 Å)** and −9.86° in Zn–N2–S
(more than three times the tolerance of 3°). In the optimum the nitrogen is pyramidal (signed volume +0.63) with
Zn–N2–H 122.8° and S–N2–H 109.1°, and Zn–N2 is 1.966 Å. The force field against the same crystal values: +2.54°,
−0.024 Å and −7.69°; in S–N2 the crystal value lies between the force field (1.572 Å) and the optimum (1.635 Å). A
re-tuning to the optimum would therefore move S–N2 past the value of the structure and away from the 1.580 Å that
stage S1 was calibrated on (plan, 5 Oct, 10:15 PM).

**What follows** (plan, 1:09 PM, point 3, which replaces the consequence named at 8:56 AM): the result is reported,
and **the ligand terms at sulfur (S–N `Y4-sy`, C–S–N `ca-sy-Y4`) are not re-tuned to this optimum**. The like-for-like
test of the ligand geometry is the validation run in the protein against 6RL9 (gate G-V: C4–S–N2 within 3°, τ_fold
within 20°, ligand RMSD), which stands as defined. A fit of part A could not have made the check pass in any case,
because two of the three quantities depend on terms that fit does not touch.

### 4.5 The force-field energy split by part [`mm_2e_main.json`, `points.<id>.split`]

For every point: bonds; the three angle terms at N2 and all other angles; every torsion term about N2–S (the
H–N2–S–C4 series; the terms through the zinc and H–N2–S–O are zero in this force field) and all other torsions;
the non-bonded energy among the 18 ligand atoms; the Lennard-Jones energy between ligand and fixed atoms (with its
r⁻¹² and r⁻⁶ parts) and the electrostatic energy between them, the scaled 1–4 pairs through the zinc included. The
non-bonded parts are summed in `gr_mm.py` from the particle and exception parameters of OpenMM's own system; their
sum equals OpenMM's non-bonded energy and the parts add up to the total within 10⁻⁵ kcal mol⁻¹ for every point
(asserted). The table `gr_table.csv` carries the split relative to X for every point; `gr_gate.json` carries it for
every comparison of part B.

## 5. Charges with the equal atoms made equal [`md/qm/resp_equiv.py`]

**Mechanism.** MCPB.py's two resp inputs are kept, with one change: in the first stage (all atoms free, restraint
weight 0.0005) the second atom of each pair is declared equivalent to the first: O2 = O1, C6 = C2, C5 = C3, H6 = H2,
H5 = H3, HN12 = HN11. The second stage (weight 0.001) refits only CH₂ and CH₃ groups of the protein residues and
keeps the ligand atoms, so the pairs stay equal. All constraints MCPB.py wrote stay in force (total charge, neutral
caps, backbone atoms at their ff14SB charges), so the fit itself keeps the total charge. This is the first stage and
not the second, as the review suggested: with the equivalence in the first stage all other charges adapt to it.

**Both sets.** The constrained charges go into the five mol2 files (used by everything after); the files MCPB.py
wrote are kept unchanged in `charges_unconstrained/`; `resp_equiv.json` holds both, their difference, the pair
average of the unconstrained set for comparison, and the checks. A second run gives the same files.

**Change in `md/qm/after_qm.sh`** (12 lines): the script is called right after MCPB.py step 3 and stops the pipeline
if one of its checks fails. Tested through the whole of `after_qm.sh` for step 2e and for step 2s on the synthetic
files: all four hard checks of `after_qm_check.py` pass both times [`selftest/after_qm_with_equivalence.json`]. I
did not run `md/qm/test_after_qm.py` again, so `md/qm/test_after_qm.json` is as its author left it; a new run would
report the zinc charge of the constrained set (1.838) where `md/ZINC.md` quotes 1.878.

**Self-test on synthetic potentials** [`selftest/resp_equiv_selftest.json`; passes]: the potential of known charges
on the 101 atoms of the large model at its 20,357 fit points, through MCPB.py step 3 and this script.

| Test | Result |
|---|---|
| Known charges with equal pairs; constrained fit without restraint | comes back to 7 × 10⁻⁶ e |
| Known charges with O1 ≠ O2 (−0.736, −0.780); unconstrained fit without restraint | comes back to 7 × 10⁻⁶ e |
| The same; constrained fit, first stage, against my own constrained least squares on the `.esp` file | 5 × 10⁻⁷ e |
| The same; pairs in the constrained set | equal to the last digit; site charge 1.000 in both sets |
| The same with the restraint (the path that will run): what the constraint changes | ligand by up to 0.054 e (O2 −0.757 → −0.719, S 1.441 → 1.388), other atoms by up to 0.071 e, zinc 1.878 → 1.838; relative r.m.s. of the fit 0.0023 → 0.0037 |
| Second run in the same folder | identical files |

Note the last but one line: forcing O1 = O2 on a potential that wants them 0.044 e apart moves the zinc and its
neighbours too, by more than a pair average would. With the real potential the size of this shift is the measure
of how much the charges were favouring pose X.

## 6. Reading the energies [`md/qm/gr_fit.py`, `md/qm/gr_queue.py`]

### 6.1 The quantum-chemical energies and the dispersion term

**Single points.** B3LYP (libxc 402, the B3LYP of Gaussian), 6-31G* with pure d and f functions, density fitting
with def2-universal-JKFIT, PySCF grid level 3, convergence 10⁻⁹ Hartree, at most 60 cycles and then the second-order
solver: the settings of stage S2 for the small model (`md/qm/config_prod.json`). Each result is in `qm/<id>.json`
with its cycles, times, starting orbitals, machine and versions. X, R and Rx were computed on the Pi before 9:05 AM
(X from the orbitals of the timing job, 3 cycles; R and Rx from orbitals taken over from another point, 32 and 23
cycles); every other point on the rented machine from the atomic guess (14 to 16 cycles and 2.4 to 3.1 minutes with 7
threads; d2 was stopped for the restart of the queue at 9:39 AM and continued from its own orbitals).

**The two machines agree.** X computed again on the rented machine from the atomic guess: 15 cycles, 166 s, and an
energy equal to the Pi's within 9 × 10⁻¹³ Hartree [`check_vm/agreement.json`].

**Dispersion** (plan, 8:56 AM, point 1). B3LYP has no dispersion term; the force field has one in its Lennard-Jones
terms. The D3 correction with Becke–Johnson damping for B3LYP is added to every quantum-chemical energy of the gate:
the two-body term (s6 = 1, s8 = 1.9889, a1 = 0.3981, a2 = 4.4211 bohr), without the three-body term. It needs no new
calculation, only the geometry. Every number of the gate is also written without it (`without_d3`); the verdicts use
the numbers with it.

- **Which program gives the numbers.** A NumPy routine of the project (`md/qm/gr_queue.py`, class `D3`), because the
  rented machine has no dispersion library and the relaxations of section 6.6 need the term there. Its reference
  data (reference coordination numbers and C6 coefficients, covalent radii, ⟨r⁴⟩/⟨r²⟩ for H, C, N, O, S and Zn) were
  read out of the source of the s-dftd3 library, tag v1.2.1, by `gr_queue.py d3-extract` [`d3_reference.json`, with
  the checksums of the three source files].
- **Checked on the Pi** [`d3.json`, `checks`]: against the dftd3 library 1.6.0 for all 38 geometries, largest
  difference 3.9 × 10⁻¹³ Hartree; the library's own parameter set for "b3lyp" gives the same numbers as the four
  parameters above (difference 0); the gradient (central differences, step 0.001 bohr) against the library's
  analytic gradient at X and Rx, largest difference 1.2 × 10⁻⁸ Hartree bohr⁻¹; against PySCF's own dispersion
  interface (pyscf-dispersion 1.5.0) at X, 3.4 × 10⁻¹³ Hartree.
- **Size** [`d3.json`, `rel_X_kcal_mol`]: relative to X, −9.09 kcal mol⁻¹ for Rx, −8.84 to −9.25 for the (g) points,
  −6.18 for R as built, −5.42 for c2_m015, −3.23 for a6, between −3.2 and +1.7 for the other turns, and within
  ±1.5 for everything made from X by moving atoms at N2. The correction is almost wholly between the ligand and
  the fixed atoms (`interaction_rel_X_kcal_mol`).
- **Left out:** the three-body term would change X against Rx and the (g) points by +0.30 kcal mol⁻¹ and no other
  point by more than 0.3 [`d3.json`, `library_d3bj_with_three_body`]. The basis-set superposition error at 6-31G*,
  which acts in the opposite direction to dispersion, is not quantified (section 9).

### 6.2 The two parts

With E_QM = B3LYP energy + D3(BJ) energy and E_MM the force-field energy (cap scheme "methyl"):

**Part A** [`gr_gate.json`: `with_d3.part_A`, `verdict_part_A`]. 25 points in two families, by the geometry a point
was made from:

- family X, 18 points: X, a2 to a6, d1 to d6, the four turns about Zn–N2 (e), the two bends (f);
- family Rx, 7 points: Rx and g1 to g6.

In each family the reference r is its point of lowest E_QM; d_i = [E_MM(i) − E_MM(r)] − [E_QM(i) − E_QM(r)]; the
number of part A is the root mean square of d_i over the 23 points that are not a reference. **Criterion: at most
1.0 kcal mol⁻¹.** The verdict exists when all 25 points are finished (`complete`).

Why two references. The plan words part A as one set "relative to the lowest point". The (g) points differ from X
by the whole difference between X and Rx, which is contact energy and belongs to part B; against a single reference
that difference would enter part A six times. Against their own base Rx they test the hydrogen terms, as (d) does
against X. This was fixed at 9:20 AM, before any energy of the (e), (f) and (g) points existed. The number with one
reference for all 25 points is written beside (`rms_single_reference`) and does not decide. Also beside: the number
without the points flagged "protein", and for the other two cap schemes (`part_A_by_cap_scheme`).

**Part B** [`with_d3.part_B`]. For each comparison: E_QM, E_MM, their difference, and the split of the force-field
energy (section 4.5). X against Rx; Nx against X; each turn about N2–S against the point it was turned from (c1
against X, c2 against Rx; also given against X); the relaxed references (section 6.6). The r.m.s. over the ten turns
is given as a summary. No criterion, no verdict, no fit.

**Reported only** [`with_d3.reported_only`]: R and N as built against X, with the split.

The numbers of the first version of the gate (sets `gate_plan` and `gate_rigid` of the plan entries of 3:02 and 8:01
AM, without D3) are kept in `gr_gate.json` under `first_version_not_a_verdict`, so that the earlier log lines can be
reproduced. They are not verdicts. (Until 1:30 PM IST the block was written twice, as `first_version_without_d3` and
as `gate`; `gate` was read by the controller of before 9 AM, which no longer runs, and its pass flags could be taken
for the verdict. The verdict is `verdict_part_A`.)

### 6.3 The fit, and the rule for accepting it

Run only when part A is complete and has failed. Everything but the terms at N2 is fixed: y_i = E_QM(i) − rest(i),
where rest is the force-field energy without the three angle terms at N2 and without the H–N2–S–C4 terms (computed two
ways; they agree to 10⁻¹³ kcal mol⁻¹). Model: y ≈ c_family + Σ k_a (θ_a − θ0_a)² + c₁ cos φ + c₂ cos 2φ, a = Zn–N2–S,
Zn–N2–H, S–N2–H, φ = H–N2–S–C4, with one free constant per family. Both variants are linear least-squares problems;
the points are those of part A.

| Model | Free parameters | |
|---|---|---|
| M1 | 5: three equilibrium angles, two torsion coefficients | force constants kept |
| M2 | 8: the same and the three force constants | needs at least 16 points |

**A model is acceptable when all of this holds** (plan, 8:56 AM, point 2, with the bounds of the first version;
conditions 4 and 5 as amended by the plan entry of 1:09 PM, point 4, for every reading of the gate from then on):

1. every force constant in 10–300 kcal mol⁻¹ rad⁻²; every standard error of an equilibrium angle at most 5° and (M2)
   of a force constant at most 30 %;
2. residual r.m.s. at most 1.0 kcal mol⁻¹, and leave-one-out r.m.s. at most 1.0 kcal mol⁻¹;
3. |c₁| and |c₂| below 5 kcal mol⁻¹;
4. every equilibrium angle within 10° of its reference. **Zn–N2–S: the crystal value, 117.55° (start X). Zn–N2–H
   and S–N2–H: their values in the quantum-chemical optimum of the small model with all atoms free, 122.76° and
   109.08°** (the geometry task `qm_opt` delivers; read from `binding_check.json`, key `quantum_chemical_optimum`,
   which
   names the file it was read from), and no longer the constructed hydrogen of start X (107.41° for both). **Until
   that optimum exists the condition cannot be met: it says so and no fit is accepted** [`fit/fit.json`,
   `window_centres`, `window_centres_record`];
5. the number of part A with the fitted terms lower than before **and at most 1.0**.

M2 is taken when it is acceptable and its leave-one-out r.m.s. is at least 20 % below that of M1; otherwise M1 when it
is acceptable. **Otherwise no parameter file is written, the terms stay as built, and gate G-R is reported as
failed**, with the split by energy part of every deviation of part A [`fit/fit.json`: `failed`, `part_A_split`]. A
model that meets every condition except "at most 1.0" (part A falls but stays above the criterion) is reported as
**"terms improved, gate still failed"** [`fit/fit.json`, `outcome`]; it is not used and no parameter file is written.

For an accepted model the queue rebuilds the topology with `fit/CA2_n2_fit.frcmod` (the terms `M1-Y4-sy`, `M1-Y4-hn`,
`hn-Y4-sy`, `ca-sy-Y4-hn`, loaded after `CA2_mcpbpy.frcmod`) and evaluates the gate again from that topology
[`gr_fit_gate.json`]. That number is the one that counts: if it does not fall, or falls and stays above 1.0, the file
is set aside as rejected and the log says which ("terms improved, gate still failed").

**Why the two additions** (check of 6 Oct 2026, findings 4 and 5). (a) Residuals carry a free constant per family,
part A is measured against one reference point per family, so a fit can have residual and leave-one-out below 1.0 and
leave part A above it. X and the four turns about Zn–N2 have the same three angles and the same H–N2–S–C4, so no
fitted term can change their deviations against each other (X −0.51, e_p015 −1.93, e_p030 −2.94, e_m120 −1.26 against
e_p120), and the reference of family X is one of them. (b) The window of condition 4 sat on a constructed hydrogen:
quantum chemistry has Zn–N2–H at 122.8°, 15.3° from the construction, so a fit that put the angle where quantum
chemistry has it would have been refused, while the term as built for S–N2–H (118.06°, GAFF2) lay 10.65° outside the
old window and lies 8.98° from the new centre.

**Today's reading is not changed by the amendment** [`fit/fit.json`, written 3:48 PM IST;
`selftest/fit_rule_cases.json`]. The fitted parameters of M1 and M2 are those of 11 AM to the last digit. M1 failed on
residual and leave-one-out and now also on "part A at most 1.0" (1.33; 1.60 against a single reference); M2 failed on
its force constants, its standard errors and its angles (Zn–N2–S 105.4° against 117.6°, Zn–N2–H 141.7° against 122.8°)
and now also on "part A at most 1.0" (1.13; 1.39). No model is chosen under any of the three cap schemes.

### 6.4 Self-test of gate and fit with the amended rule [`selftest/fit_selftest.json`, written 1:40 PM IST; passes]

Synthetic energies on the 38 real geometries, with the real force-field model (`mm_2e_main.json`); noise where stated:
0.3 kcal mol⁻¹, 200 draws, seed 20261006. The equilibrium angles of the known parameter sets are given as offsets from
the centres of the window of condition 4, so the test follows the rule; with the centres of now (117.55, 122.76,
109.08°) the set of T1 is 121.00, 127.35, 111.67° (it was 121, 112, 110° under the centres of 11 AM, with the same
offsets). **The numbers of T1 to T3 with noise therefore differ from those of the 11:30 AM version because the test
parameters moved with the window**; T0 to T0c are unchanged.

| Test | Energies given to the code | Right answer | Result |
|---|---|---|---|
| T0 | the force field itself | gate passes, nothing fitted | part A 0.000; with noise the gate passes in 200 of 200 draws, no fit |
| T0b, the reviewer's test | force field minus D3(BJ), read without the D3 term | "change nothing" | part A 0.945, passes, nothing fitted. (Against a single reference the same data give 5.185.) Read with the D3 term: 0.000 |
| T0c | force field + 3.4 on every point made from Rx + 2.0 on every turn about N2–S | part A passes | part A 0.000 (single reference: 1.84) |
| T1 | known other equilibrium angles (centres + 3.45, 4.59, 2.59°) and torsion coefficients (0.8, −1.5), force constants of the topology | fit, accepted, parameters back | part A 7.91 before (6.36 single reference); M1 and M2 return the set to 10⁻¹²; M1 taken; part A 0.000 after. With noise: a set is released in 200 of 200 draws (M1 195, M2 5); standard deviation of the angles 0.29–0.58°, of the coefficients 0.14–0.15 |
| T2 | the same with other force constants (38, 28, 45) | M2, accepted | M1 is refused (leave-one-out 1.28); M2 returns the set and is taken. With noise: M2 in 190 draws, none in 10; angles ± 0.5–2.5°, force constants ± 3.3–4.1 |
| T3 | known parameters outside the rule (Zn–N2–S 13.45° from its centre, c₂ = −6.5) | no parameter file | both models reproduce the data exactly and both are refused (conditions 3 and 4); with noise: none in 200 of 200 draws |
| T4, new | the data of T1, read as if the quantum-chemical optimum did not exist | no parameter file, for that reason alone | both models refused on condition 4 only |
| T5, new | parameters inside the rule; one of the five points with identical fitted coordinates (e_p120) lowered by 2.4 | "terms improved, gate still failed", no parameter file | part A 8.23 before; M1 residual 0.46, leave-one-out 0.51, part A 2.00 after; refused on "at most 1.0" alone, outcome "terms improved, gate still failed" |
| T6, new | exact data with Zn–N2–H 11° from the optimum | no parameter file | both models refused on condition 4 only |
| T7, new | exact data with both hydrogen angles at the optimum (15.35° and 1.67° from the constructed hydrogen) | accepted, parameters back | M1 taken, part A 0.000 after; the rule of 11 AM would have refused it |

The model form itself was checked first: with the parameters of the topology it reproduces the force field's own terms
at N2 for all 38 points up to a constant (spread 7 × 10⁻⁶ kcal mol⁻¹).

**The nine constructed cases of the independent check, run again on the program of now**
[`selftest/fit_rule_cases.json`, written 1:41 PM IST; the checker's cases and seed 20261007, its script
`mdgr2/grcheck/c3_fit.py` in the session's scratch]. (I) Exactly as the checker built them and read with the centres
it had: 5 refused and 4 accepted, the same choice as the checker's in all nine, and no residual, leave-one-out value,
part A number or angle differs (largest difference 0). (II) The same nine with every equilibrium angle moved by the
shift of its window centre and read with the centres of now: the same five refused and the same four accepted (in A4
the other of the two acceptable models is taken). (III) The nine unmoved, read with the centres of now: none is
accepted, and the four that were accepted are refused because their Zn–N2–H of 110° is 12.76° from the optimum. The
checker's gap case (its G1, accepted by the rule of 11 AM with part A at 2.00 after the fit) is now refused as "terms
improved, gate still failed" (part A 2.00). The two conditions on the number of part A alone: worse, equal, and "falls
to 1.19" are refused; "falls to 1.00" and "falls to 0.95" are accepted.

### 6.5 Why the points (d) [self-test of the first version, 5 AM, placeholder charges]

A move that changes both hydrogen angles together, as d1 or the path, fixes only a combination of the Zn–N2–H and
S–N2–H terms. Without d3 to d6, which move one angle each, the force constants of model M2 came back with standard
deviations of 1,200 and 460 kcal mol⁻¹ rad⁻² in a synthetic test; with them every parameter is determined.

### 6.6 The relaxed references of part B [`relax_mm.json`, `relax_qm.json`, `relax/`]

For each of the two configurations at N2 the ligand is relaxed with the zinc and the imidazoles fixed, from X and from
Rx, in both descriptions. At every step the signed volume at N2 (+0.72 in X, −0.68 in Rx) is recorded, so that a
relaxation cannot cross to the other configuration unnoticed.

**Force field** (`gr_mm.py relax`; cap scheme "methyl", final charges, terms as built; computed again at 3:48 PM IST
with the same numbers as at 10:25 AM):

| Start | Energy, start → minimum (kcal mol⁻¹) | Heavy atoms moved (r.m.s.) | Signed volume at N2, start → end |
|---|---|---|---|
| X | 72.77 → 62.68 | 2.25 Å | +0.724 → +0.783 |
| Rx | 86.99 → 62.68 | 5.38 Å | −0.678 → **+0.783** |

Both relaxations end in the same minimum (Zn–N2–S 107.5°, C4–S–N2 111.3°, S–N2 1.571 Å, Zn–N2 1.905 Å, Zn–N2–S–C4
−154.1°). The minimiser (L-BFGS) could have stepped over a barrier, so the relaxation from Rx was repeated as a
descent that only accepts a lower energy and moves no atom by more than 0.02 Å in a step: the signed volume changes
sign 15.3 kcal mol⁻¹ below the start and the descent arrives at the same minimum (62.79 after 39,220 steps, still
falling). **In the force field of the small model the ring-rotated placement is not a minimum; it inverts at N2
without a barrier.** There is therefore no force-field difference between two relaxed configurations, and
`relax_mm.json` says so (`configurations_kept`: false). On the force-field side a crossing is read with a threshold:
the two configurations are kept only if the minimum from X has a signed volume above +0.3 and the minimum from Rx one
below −0.3.

**Quantum chemistry** (`gr_queue.py relax`; tasks `gr_relax_Rx`, `gr_relax_X` and `gr_relax_Rx_cont` of the rented
machine): B3LYP/6-31G* with the settings of the single points (pure functions, density fitting, grid 3) and the D3(BJ)
energy and gradient of class `D3` in every step, geomeTRIC 1.1.1 with its own convergence criteria, atoms 1 to 46
frozen.

**What "crossed" means** (plan, 1:09 PM, point 1; check, finding 1). Until 1:09 PM the program stopped a relaxation at
the first step at which the signed volume at N2 had changed sign, called it "crossed", and section 6.6 of the version
of 11:30 AM read that as "the two descriptions agree that the small model has one minimum". That reading is
**withdrawn**: a run cut at a planar nitrogen, unconverged, shows where it was, not where it ends. The rule now, in
the program since 1:33 PM IST:

- a relaxation is **never stopped at a change of sign**; it runs until geomeTRIC's criteria are met;
- it **counts as crossed only when it has converged with a signed volume of the sign opposite to that of its start**
  (`crossed` in `relax/qm_<start>.json` and `relax_qm.json`); while it runs, `crossed` is false or empty, and the step
  of the first change of sign is recorded beside (`sign_first_changed_at_step`);
- a record with the old state "crossed" is not an end: the run goes on from its last geometry, and the summary names
  that state "stopped at the change of sign by the old rule, NOT converged";
- the end is described by numbers, with one word for the nitrogen from the threshold the force-field side already
  used: planar or nearly planar below 0.3 in the absolute signed volume, pyramidal above (`nitrogen_at_the_end`);
- `relax_qm.json` gives E(from Rx) − E(from X) when both have converged, whether or not one has crossed, and says
  `configurations_kept` only then; before that it gives the two runs as far as they are, marked as not a result.

Every step records the energy with and without D3(BJ), the signed volume, the three angles at N2 and their sum,
Zn–N2–S–C4, C4–S–N2, the distance from the start, the largest shift of a frozen atom and, since 1:33 PM, the distance
of the ligand (heavy atoms, r.m.s., no fitting, because the frame is common) and the energy from where the relaxation
from the other start has got to.

**"Zinc and imidazoles fixed" in the quantum-chemical runs** (check, finding 6). geomeTRIC holds atoms named under
`$freeze` only approximately unless its keyword `enforce` is set. Measured on the trajectories of the machine
[`relax/fixed_atom_shift.json`, 3:48 PM IST]: in the 16 steps from Rx before the continuation the frozen atoms drifted
by up to 0.0047 Å (0.0032 Å at step 16); in the 7 continued step(s), with `enforce`, by at most 7e-08 Å. In the 55
steps from X, a finished run without `enforce`, by up to 0.0117 Å (at step 39; 0.0001 Å at its last step). The
continued run sets `enforce` 0.1 (`queue_vm_relax2.json`: the frozen atoms are held exactly once they are within 0.1
bohr of their places) and puts the frozen atoms back exactly at every start. Tested on the Pi before use
[`selftest/relax_continue_and_enforce_test.json`; not a result of the gate]. On a small case with real forces
(formamide and water, HF/STO-3G, a distorted start, three atoms of the formamide frozen, so that frozen atoms are
bonded to free ones as in the zinc model): without `enforce` a frozen atom moves by up to 0.0052 Å during the run,
with it by 1e-07 Å; both converge (63 and 69 steps) to the same energy within 9e-08 Hartree. (With the water frozen
instead, a molecule of its own, geomeTRIC holds the atoms exactly with and without the keyword.) Then with a surrogate
energy in place of the quantum chemistry and the real 64 atoms, constraints and program: without `enforce` the frozen
atoms move by up to 0.024 Å during the run (0.00015 Å at its end), with it by 1e-07 Å; both runs converge (99 and 98
steps) to the same end (largest distance of an atom 0.0039 Å). The same test showed that a run which the old program
had stopped at the sign change (step 25 of the surrogate) is taken up from that geometry, converges (113 steps) and is
recorded as crossed only at its end. `gr_mm.py relax` accepts a quantum-chemical geometry whose frozen atoms lie
within 0.01 Å of X, writes the shift into `relax_mm.json`, and gives the force-field single point at the geometry as
it is and with the frozen atoms put back (their difference is the third column of the table at the end of this
section).

**State of both runs, as of 3:50 PM IST** [`relax_qm.json`, written 3:48 PM IST; `relax/qm_X.json`,
`relax/qm_Rx.json`]:

- from X: **converged** after 55 steps (3:11 PM IST), 13.44 kcal mol⁻¹ below its first step; signed volume at N2
  +0.724 → +0.629, angles at N2 104.1, 125.4 and 109.8° (Zn–N2–S, Zn–N2–H, S–N2–H; sum 339.3°), Zn–N2–S–C4 −148.8°,
  heavy atoms of the ligand 2.29 Å from the start; pyramidal, configuration of the start; crossed by the rule: no.
- from Rx: **not converged**, 23 steps so far (last step 3:45 PM IST), 19.30 kcal mol⁻¹ below its first step, r.m.s.
  gradient on the ligand 7.4e-04 Hartree bohr⁻¹; signed volume at N2 −0.678 → +0.273, angles at N2 117.9, 128.6 and
  110.1° (Zn–N2–S, Zn–N2–H, S–N2–H; sum 356.5°), Zn–N2–S–C4 −80.1°, heavy atoms of the ligand 1.44 Å from the start.
- **Note added on 11 Oct 2026: the run from Rx has converged** at step 66 (6 Oct 2026, 6:49:52 PM IST), 0.12 Å
  (heavy atoms of the ligand) and −0.01 kcal mol⁻¹ from the end point of the run from X, with the nitrogen pyramidal
  in the configuration of X (signed volume +0.619); the numbers are in the note of section 1 and in
  `relax/qm_Rx_converged_step66.json`. The bullets below on the two geometries "as they are" and the table of the
  steps (which ends at step 23) describe the state of 3:50 PM IST and are kept as written.
- Between the two geometries as they are (NOT both converged, so this is the distance between two runs on their way,
  not a result): 4.70 Å (heavy atoms of the ligand, r.m.s., no fitting), largest distance of one atom 8.22 Å; E(from
  Rx) − E(from X) = +4.98 kcal mol⁻¹ with D3(BJ), +6.38 without [`relax_qm.json`, `end_points`].
- The force field at the same two geometries (single points of the cut-down model, as section 4): E(from Rx) − E(from
  X) = +5.48 kcal mol⁻¹, of which electrostatic energy between ligand and fixed atoms +8.98, Lennard-Jones between
  them −0.34, the three angles at N2 −0.94, other angles −2.03, torsions about N2–S +1.10, other torsions +0.86, bonds
  −0.71, non-bonded inside the ligand −1.44 [`relax_mm.json`,
  `Rx_minus_X.force_field_at_the_quantum_chemical_geometries`]; **not a result while a run is unconverged**.
- On the machine (queue status of 3:47 PM IST): `gr_relax_X` done at 3:11 PM; `gr_relax_Rx_cont`, the continuation
  from Rx, running on 7 CPUs since 3:11 PM (attempt 1), beside the single points of gate G-R2 (`n2_points`) on 16.
  Nothing else of gate G-R is running; `gr_binding_final` waits for `qm_write`.
- **Where the run from X ended.** 2.29 Å (heavy atoms) from the crystal placement, with the ring turned by 32° about
  N2–S (Zn–N2–S–C4 −148.8° against −117.0°) and by 11° about Zn–N2, and the nitrogen still pyramidal in the
  configuration of X (signed volume 0.63, with the hydrogen opened to Zn–N2–H 125.4°). That is 0.26 Å from the
  force-field minimum reached from X: on where the crystal-like ligand settles when only the zinc site holds it, the
  two descriptions agree. They differ in the geometry at N2 and at sulfur there (Zn–N2–H 125.4° against 106.9°, S–N2
  1.634 against 1.571 Å, C4–S–N2 106.4° against 111.3°), which is why the force field puts the quantum-chemical end
  point 19.2 kcal mol⁻¹ above its own minimum; they are the differences the binding check shows with all atoms free
  (section 4.4).

What the first 16 steps from Rx do show, and it is all they show: the energy fell at every step (checked from the
record), 18.82 kcal mol⁻¹ in all, from a pyramidal nitrogen (signed volume −0.678) to a planar one (+0.097, angle sum
359.6°, Zn–N2–S–C4 −76.4°), with a gradient still 2.4 times geomeTRIC's criterion. **The R-like placement with a
pyramidal nitrogen is not a minimum of the zinc-site model in quantum chemistry either.** Whether the run ends in the
minimum reached from X, or in a minimum of its own with a flat nitrogen and the ring still turned, is what the
continuation decides.

**Where the two runs are going** [`relax/qm_turns.json`, written 3:48 PM IST, from the trajectories on the machine].
Two torsions describe it: the turn of the ligand about Zn–N2 (NE2(His94)–Zn–N2–S: 74.2° in X, −11.5° in Rx) and the
turn about N2–S (Zn–N2–S–C4: −117.0° and −21.3°). The run from X, at step 55: 85.2° about Zn–N2 and −148.8° about
N2–S, with N1 3.82 Å from its place in the crystal pose and the ligand 0.26 Å (heavy atoms) from the force-field
minimum, which has 92.3° and −154.1°; the quantum-chemical optimum of the model with all atoms free has 79.0° and
−144.2°. So with the frame fixed and nothing around it the ligand of the crystal pose turns about N2–S; in the protein
a turn of 30° there runs the ring into the pocket wall (section 3.3). The run from Rx, at step 23: 29.0° about Zn–N2
and −80.1° about N2–S, signed volume +0.273, 4.70 Å (heavy atoms) from the last geometry of the run from X.

**Estimate of the end of the continued run** (an estimate, not a result). Made at 2:14 PM IST from the 16 steps on
file, before the continuation had started. In its last three steps the run from Rx turned by 0.4, 2.0 and 4.6° about
Zn–N2 and by 0.2, 1.8 and 4.5° about N2–S, both towards the values of the run from X, while the energy fell by 0.05,
0.08 and 0.16 kcal mol⁻¹ a step; its gradient was near geomeTRIC's criterion of 3.0e-04 at steps 13 to 15 (4.0e-04 to
4.8e-04) and rose again at step 16 (7.3e-04). It has about 52° left about Zn–N2 and 52 to 68° about N2–S, and stands
5.2 kcal mol⁻¹ above the run from X. The run from X shows how flat this surface is: its gradient was 1.7e-04 to
2.2e-04 at steps 15 to 17, and it then moved on by another 0.9 Å in 19 steps. Both ends are therefore open. If the run
from Rx slides on at the pace of its last steps, it ends in the basin reached from X after some 40 to 80 further
steps, three to five hours at four minutes a step. If it converges near where it stands, it does so within 10 to 20
steps, as a flat nitrogen with the ring still turned, about 5 kcal mol⁻¹ above the end of the run from X. I take the
first as the more likely, because the gradient was growing and both turns were moving faster when the run was cut;
this is a guess from 16 steps, to be replaced by the run itself. **Revised at 3:50 PM IST from the first 7 steps of
the continuation** (steps 17 to 23, 3:17 to 3:45 PM; table below): the run is not settling at the planar nitrogen. The
signed volume has grown from +0.099 to +0.273, towards the configuration of X; the distance from the end point of the
run from X has fallen from 4.903 to 4.697 Å, by 0.002, 0.006, 0.021, 0.038, 0.055 and 0.084 Å in successive steps, so
the pace is still growing; the energy above that end point has fallen from 5.45 to 4.98 kcal mol⁻¹. With 4.7 Å still
to go, and if the pace grows to the 0.1 to 0.2 Å a step that the run from X reached on the same kind of slope and the
approach is followed by a tail of 10 to 15 steps as there, the end comes after some 35 to 65 further steps: about
three to five hours at 4.6 minutes a step, between about 6:30 and 9:00 PM IST, and probably in the basin reached from
X. That it ends there is an expectation from seven steps, not a result; a shallow minimum on the way is not excluded.

**The continued relaxation from Rx, step by step** (7 steps since the stop at step 16; energies in kcal mol⁻¹, angles
in degrees; "X run" is the end point of the relaxation from X when it has converged, else its last step):

| Step | Time (IST) | E − E(step 16) | Signed volume | Zn–N2–S | Zn–N2–H | S–N2–H | Sum | Zn–N2–S–C4 | From the X run (Å) | E − E(X run) | r.m.s. gradient |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 17 | 3:17 PM | −0.01 | +0.099 | 117.4 | 131.9 | 110.3 | 359.5 | −76.6 | 4.903 | +5.45 | 7.3e-04 |
| 18 | 3:22 PM | +0.00 | +0.116 | 116.1 | 132.2 | 111.1 | 359.4 | −77.8 | 4.901 | +5.47 | 7.5e-04 |
| 19 | 3:26 PM | −0.08 | +0.116 | 117.0 | 131.8 | 110.6 | 359.4 | −77.5 | 4.895 | +5.38 | 3.9e-04 |
| 20 | 3:31 PM | −0.14 | +0.131 | 117.5 | 131.2 | 110.5 | 359.2 | −77.7 | 4.874 | +5.33 | 3.6e-04 |
| 21 | 3:35 PM | −0.20 | +0.161 | 118.1 | 130.3 | 110.4 | 358.8 | −78.2 | 4.836 | +5.27 | 6.3e-04 |
| 22 | 3:40 PM | −0.31 | +0.206 | 118.3 | 129.5 | 110.3 | 358.0 | −78.9 | 4.781 | +5.16 | 6.6e-04 |
| 23 | 3:45 PM | −0.49 | +0.273 | 117.9 | 128.6 | 110.1 | 356.5 | −80.1 | 4.697 | +4.98 | 7.4e-04 |

**The force field at the quantum-chemical geometries** [`relax_mm.json`, `at_other_geometries`; kcal mol⁻¹ and Å]:

| Geometry | Largest shift of a frozen atom | Energy with the frozen atoms put back, minus as it is | Above the force-field minimum | Ligand from start X | Ligand from the force-field minimum |
|---|---|---|---|---|---|
| end point from X | 0.0001 | −0.003 | 19.17 | 2.29 | 0.26 |
| from Rx, last step (not converged) | 0.0000 | +0.000 | 24.64 | 3.06 | 4.74 |

> **Note added on 11 Oct 2026.** The second row is the geometry of step 23. The run converged at step 66 (note in
> section 1); `relax_mm.json` in git is that of 3:48 PM IST of 6 Oct 2026 and holds no force-field single point at
> the converged geometry. The force field against quantum chemistry along the whole finished run is in
> `md/qm/n2/final_test/behaviour.json` (`along_the_quantum_chemical_relaxations.runs.Rx`) and in `md/N2.md`.

**In one sentence.** So far the relaxations show that in the bare zinc-site model the crystal-like ligand leaves the
crystal placement in both descriptions and settles in the same place in both, with a pyramidal nitrogen, and that the
R-like placement with a pyramidal nitrogen is a minimum in neither; they do not yet show where the R-like placement
ends in quantum chemistry, and they cannot show what the ligand does in the protein, because the small model in vacuum
has no Thr199, no pocket wall and no water.

> **Note added on 11 Oct 2026.** Where the R-like placement ends in quantum chemistry is now known: in the minimum
> reached from X (0.12 Å, −0.01 kcal mol⁻¹; `relax/qm_Rx_converged_step66.json`). The second half of the sentence
> stands as written: the small model in vacuum has no Thr199, no pocket wall and no water, and this relaxation does
> not show what the ligand does in the protein.

**How the two descriptions are read together** (as fixed at 8:56 AM, with the correction of 1:09 PM). If quantum
chemistry keeps the ring-rotated configuration as a minimum and the force field does not, that is the mismatch the
plan wants stated as a limit of the model, and it bears directly on the simulation: a conversion of the rotated start
in the protein could then come from the force field. If the quantum-chemical run converges in the minimum reached from
X, the two descriptions agree that the small model has one minimum, and the model says nothing about the rotated
placement beyond that. If it converges somewhere else (a flat nitrogen with the ring still turned), the two
descriptions differ in where the R-like placement goes, and the energy difference between the two end points is given
in both, the force-field one as single points at the quantum-chemical geometries. In every case the model has no
Thr199, no pocket wall and no water, so none of the three says what the ligand does in the protein.

## 7. The queue and the two machines [`md/qm/gr_queue.py`, `queue_prod.json`, `queue_vm.json`, `queue_vm_relax2.json`, `md/gpu/tasks_gr.json`, `md/gpu/tasks_gr3.json`]

**What ran where.**

- **Pi.** The queue of the first version ran from 4:55 AM (PID 250515), computed X, R and Rx, and was stopped cleanly
  by the file `/root/cqc-sfa-md/qm/gr/STOP` at 9:04:31 AM, one second after it had begun a4; a4 was abandoned after
  one cycle and its files were moved out of the queue's folder. The STOP file is still there, so the Pi queue cannot
  start by accident. The Pi now only evaluates (`gr_queue.py evaluate`).
- **Rented machine** (`cqc-md-gpu-01`; tasks of `md/gpu/tasks_gr.json` in its runner `md/gpu/queue.py`):
  `gr_check_x` (X again, done 9:10 AM), `gr_points` (the same queue program with `queue_vm.json`, running since
  9:10 AM, restarted at 9:39 AM with the code of this version and the 12 added points), then `gr_relax_Rx` and
  `gr_relax_X`. Each uses 7 CPUs and 7 threads, beside the 16 of the optimisation of the small model. Not 8: the
  runner counts one CPU for a GPU task, and 16 + 1 + 8 is more than the 24 of the machine, so a task of 8 would have
  waited for every GPU task and, once running, kept every GPU task waiting. `gr_relax_Rx` ended at 11:57 AM at the
  change of sign (old rule); `gr_relax_X` has run since then. `gr_binding` (the provisional binding check, 12:46 PM)
  and `gr_binding_final` (after `qm_write`) use one CPU.
- **Rented machine, third file** (`md/gpu/tasks_gr3.json`, added 1:34 PM IST): `gr_relax_Rx_cont`, the relaxation
  from Rx continued from step 16 to convergence with `queue_vm_relax2.json` (section 6.6). 7 CPUs. **It waits for
  `gr_relax_X`**: the Hessian of the small model took over the 16 CPUs of the optimisation at 12:46 PM (`qm_hess`;
  37 gradients of seven to eight minutes each; it was cancelled at 2:35 PM IST by the work on gate G-R2 (`md/N2.md`;
  plan entry of 1:34 PM), whose single points, `n2_points`, hold those 16 CPUs since 2:37 PM), so there are 7 CPUs for
  one relaxation only (16 + 7 + 7 = 30 of 24). The
  dependency is written into the task and not left to the CPU count, because the runner starts whatever fits in the
  order of the list, this file comes before `tasks_initial.json`, and after a restart of the runner two relaxations
  would otherwise start first and keep `qm_hess` waiting. **It also waits for `gr_relax_Rx`** (done; dependency
  added 1:52 PM IST), for the case that the runner is started again while the continuation runs (section 9, point
  7). Its test holds only when the run has converged (state `done`).

**After every point** the queue computes what is missing of `d3.json`, runs the binding check and, when needed, the
force-field relaxation, evaluates the gate for the three force-field models (`gr_*`, `gr_asis_*`, `gr_free_*`) and
writes one line to `queue.log` and the block `gate_so_far` to `status.json`: part A (r.m.s. with and without D3,
number of comparisons, largest deviation, complete, passed), part B (X against Rx with and without D3, the turns,
the relaxed references in both descriptions), and the verdict of the binding check. When all points are finished
and part A has failed, it runs the fit (section 6.3). Since 1:33 PM IST the line quotes the single-reference number
beside the two-family number of part A, and a step beside the single points (D3, binding check, force-field
relaxation, summary, gate, fit) that ends with a return code other than its good ones is written into `queue.log`
and into `side_failures` of the status, where before only a traceback in `gate.log` showed it (check, finding 6).
The force-field relaxation is made again only when the model has other inputs; the single points at the
quantum-chemical geometries whenever one of those geometries has changed (`gr_mm.py relax --reuse`).

**Follow it.**

```
SSH="ssh -i /root/.ssh/aitf_pi -o IdentitiesOnly=yes -o BatchMode=yes ubuntu@<address>"
$SSH "sudo python3 /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md/gpu/queue.py status"
$SSH "sudo tail -3 /root/cqc-sfa-md/qm/gr/queue.log"                       # one line per point, both parts
$SSH "sudo grep '^relax' /root/cqc-sfa-md/logs/gr_relax_X.log | tail -3"        # one line per step of a relaxation
$SSH "sudo grep '^relax' /root/cqc-sfa-md/logs/gr_relax_Rx_cont.log | tail -3"  # the continued run from Rx
# on the Pi, after the watchdog has pulled (every 15 minutes), or after fetching md/qm/gr of the machine by hand
# (rsync -a --rsync-path="sudo rsync" of its qm/, relax/ and check_vm/ into a scratch folder); one thread, nice:
OMP_NUM_THREADS=1 OPENMM_CPU_THREADS=1 nice -n 10 /root/cqc-sfa-tools/env/bin/python md/qm/gr_queue.py evaluate md/qm/gr/queue_prod.json /root/cqc-sfa-md/from_gpu/repo_md/qm/gr
column -s, -t md/qm/gr/gr_table.csv | less -S
```

On the Pi of 6 Oct 2026, 1:39 PM IST, `evaluate` took 41 s, the force-field relaxation included.

`evaluate` never replaces a result the Pi has; a different energy for the same point on the two machines is
reported.

**Restart or end a task of the gate on the machine.** `sudo touch /root/cqc-sfa-md/qm/gr/STOP`: the running SCF or
relaxation stops at its next cycle or step and keeps its orbitals or geometry; the runner starts the task again
after 20 to 30 s, and the task removes the file. This is how new code or new points are taken up. To end a task:
`queue.py cancel gr_points` (or `gr_relax_Rx`, `gr_relax_X`, `gr_relax_Rx_cont`). The three relaxation tasks share
that STOP file: touched while two of them run, it stops both, and the first to start again removes it.

**Copying to the machine.** The files of the gate were copied with rsync to the same absolute paths and checked
there with `sha256sum -c` (three lists, 5, 87 and 5 files, in the morning). For this version five files went the
same way, with `rsync --update --checksum` and a checksum check on the machine: `md/qm/gr_queue.py`, `gr_mm.py`,
`gr_fit.py`, `md/qm/gr/queue_vm_relax2.json` and `md/gpu/tasks_gr3.json` (1:33 PM IST; `gr_queue.py` again at 1:40 PM
and `gr_fit.py` at 2:03 PM after changes outside the relaxation code; the task file again at 1:52 PM with its second
dependency).
Before the copy the three programs on the machine were those of commit 29d9e5d (same checksums). The relaxation from X
had started before the copy and runs the program of the morning to its end;
a task that starts later (the continued run, `gr_binding_final`) runs the program of now. The warning of the 11:30 AM
version about `md/gpu/push_to_vm.sh` is out of date (check, finding 11): since 11:42 AM the script copies with
`rsync --update` and lists what it left alone, and `md/protocol.json` is the same on both machines.

**Time.** A single point takes 2.4 to 3.1 minutes on the machine with 7 threads (12.6 to 62.6 minutes for the three
points of the Pi with 4). A step of a relaxation is one SCF from the orbitals of the step before and one gradient; the
time per step is written into `relax/qm_<start>.json` (`seconds_per_step_mean`). The first eight steps of the
relaxation from Rx on the machine took 264 to 321 s each (mean 287 s) with 7 threads [`relax/qm_Rx.json`].

## 8. Decisions taken here that need a dated line in `md/PLAN.md`

### 8.1 Of 9:00 to 9:30 AM: in the plan since 12:19 PM

All were taken between 9:00 and 9:30 AM IST on 6 Oct 2026. Energies read by then: X, R, Rx, and one r.m.s. of the
first version that contained a4. None of the points (e), (f), (g) existed. The plan entry of 12:19 PM names them
(two-body D3(BJ); two families; twelve added points; the stop rule of the relaxations; 7 threads); point 5 and the
stop rule of point 6 were changed by the entry of 1:09 PM (sections 6.3 and 6.6).

1. **D3(BJ) is the two-body term**, parameters of section 6.1, from the project's routine checked against two
   libraries; the three-body term (+0.30 kcal mol⁻¹ on X against Rx) is left out.
2. **Part A is read in two families**, each against its own lowest point (section 6.2); the number against a single
   reference is reported beside. Points of part A: X, a2 to a6, d1 to d6, e (4), f (2), Rx with g1 to g6.
3. **Part B**: X against Rx, Nx against X, each turn about N2–S against the point it was turned from, the relaxed
   references.
4. **The 12 added points** and the rule that chose the turns about Zn–N2 (section 3.4).
5. **The fit**: on the points of part A, one constant per family; M2 from 16 points; the conditions of section 6.3,
   in which the earlier range of 80–160° for the equilibrium angles is replaced by 10° around start X. To be
   settled by the plan's owner: whether "the structure" for the two hydrogen angles is the tetrahedral construction
   (107.41°), given that the force field as built has S–N2–H at 118.06°. **Settled at 1:09 PM**: the
   quantum-chemical optimum of the small model, and part A with the fitted terms at most 1.0 (section 6.3).
6. **Relaxed references**: D3(BJ) in the optimisation; a relaxation that crosses to the other configuration is
   stopped there and reported as crossed; geomeTRIC's own convergence criteria. **Changed at 1:09 PM**: a run is not
   stopped at a change of sign and counts as crossed only when it has converged there (section 6.6).
7. **Binding check**: a provisional verdict from the geometry of the finished optimisation, the final one from
   `CA2_small_opt.fchk`.
8. **The machine**: 7 threads for the tasks of the gate instead of 8 (section 7).

### 8.2 Of 1:15 to 3:50 PM: how the entry of 1:09 PM is made operational; not yet in the plan

Taken before any step of the continued relaxation existed and before any second reading of the gate. Read by then:
everything of the entry of 12:19 PM, the 16 steps from Rx and the steps from X up to that time.

1. **"Has passed zero and has converged"** is, in the program: the run has met geomeTRIC's criteria and the signed
   volume at N2 of its final geometry has the sign opposite to that of its start. A run that changes sign and comes
   back is not crossed. While a run is unconverged it is neither crossed nor "configuration kept". One word for the
   nitrogen at the end: planar or nearly planar below 0.3 in the absolute signed volume (the threshold the
   force-field side has used since the morning), pyramidal above.
2. **The continued run** (`queue_vm_relax2.json`, task `gr_relax_Rx_cont`): same level; geomeTRIC's `enforce` 0.1;
   the frozen atoms put back exactly at every start (0.0032 Å at most for the geometry of step 16); at most 150 steps
   an attempt and four attempts; it waits for `gr_relax_X` (section 7).
3. **Fixed atoms of a quantum-chemical geometry** in `gr_mm.py relax`: accepted within 0.01 Å of X; the force-field
   single point that is compared with the quantum-chemical energy is the one at the geometry as it is; the one with
   the frozen atoms put back is given beside.
4. **Fit rule**: the optimum for the window of the two hydrogen angles is read from `binding_check.json`
   (`quantum_chemical_optimum`): the geometry of the finished optimisation until `CA2_small_opt.fchk` exists, then
   that file (the same geometry). After an accepted fit the number of part A that must fall and be at most 1.0 is
   the one from the rebuilt topology. The outcome names are "accepted", "terms improved, gate still failed" and
   "no fit accepted".
5. **A second reading** is named by its `step2` (or by the key `reading` of its configuration) and writes into
   `md/qm/gr/reading_<name>/` and scratch folders of its own; the reading of today is "2e" and keeps its files where
   they are. The programs refuse to write one reading over another (section 9, point 8).
6. **`gr_gate.json`**: the key `gate` is gone, the first version's numbers are under `first_version_not_a_verdict`;
   new blocks `verdict_part_A.subsets`, `part_A_composition`, `reading`, `step2`, `written`.
7. **Binding check**: `binding_check.json` carries the consequence of the plan entry of 1:09 PM, point 3, and the
   optimum against the crystal with the same tolerances (`optimum_against_crystal`).

## 9. What is uncertain or not done

1. **Basis-set superposition** at 6-31G* is not quantified. It favours close contacts, as dispersion does, so the
   reference with D3 and without a counterpoise correction may overbind the ring-rotated placement; the true value
   of X against Rx probably lies between the numbers with and without D3 (an expectation: no counterpoise correction
   was computed). The review asks for X, Rx and a4 in a larger basis with a counterpoise estimate before X against Rx
   is trusted to 1 kcal mol⁻¹; not done.
2. **Caps** (section 4.3): up to 0.44 kcal mol⁻¹ on a relative energy. A cleaner treatment would fit the charges of
   the twelve cap atoms to the quantum-chemical potential of the small model itself; not done.
3. **The charges are fitted to the large model and used in the small one.** How well they reproduce the potential of
   the small model was not measured; an error there enters every comparison and, in a fit, the terms at N2.
4. **The intrinsic barrier of the three angle terms as built is small.** Their equilibrium values sum to 343.0°
   (332.4° in X); the lowest planar arrangement of the three terms alone costs 1.6 kcal mol⁻¹ (review, finding 10).
   The plan's sentence that the inversion barrier comes from the three angle terms does not hold for these
   parameters: at the planar point a4 the force field is 13.76 kcal mol⁻¹ above X, of which 8.00 are electrostatic
   energy between ligand and fixed atoms and 5.10 the three angle terms [`gr_table.csv`].
5. **The path (a) is one path**; four of its five inner points press on Thr200 in the rigid protein.
6. **The relaxation of the rotated placement in the force field crosses** (section 6.6). The plan's comparison of
   "the two relaxed configurations" has then no force-field side. A relaxation with the hydrogen bond to Thr199
   modelled (for example the hydrogen held, or OG1 added to the model) would give one; it is not in the plan and was
   not done.
7. **The quantum-chemical relaxations run on the rented machine only.** Their mechanics (frozen atoms held exactly,
   restart from the last geometry, a run stopped by the old rule taken up again, no stop at a change of sign,
   convergence, the record against the other run) were tested on the Pi with a surrogate energy in place of the SCF
   [`selftest/relax_continue_and_enforce_test.json`]. As of this version the relaxation from X has converged and the
   continued relaxation from Rx is running (step 23 at 3:45 PM IST); its end is not known. **To finish:** when
   `queue.py status` shows `gr_relax_Rx_cont` done, fetch the machine's `md/qm/gr` and run `evaluate` (section 7);
   `relax_qm.json` then holds where the run ended, the configuration at N2 and E(from Rx) − E(from X), `relax_mm.json`
   the force-field single points at the two end points with the split, and `gr_gate.json` both under
   `with_d3.part_B.relaxed`; sections 1 and 6.6 are then written again from those files.
   **Note added on 11 Oct 2026.** `gr_relax_Rx_cont` finished on 6 Oct 2026, 6:49:52 PM IST (converged, step 66). The
   watchdog brought the machine's `md/qm/gr` to `/root/cqc-sfa-md/from_gpu/repo_md/qm/gr/` on the Pi, with
   `relax_qm.json` of 6:49:52 PM IST (`both_converged` true, `Rx_minus_X_with_d3` −0.011) and `gr_gate.json` of
   8:47:30 PM IST. Sections 1 and 6.6 were not written again: the result was entered in the plan (`md/PLAN.md`,
   "Relaxed reference, now complete"), and the finished record of the run was put into git on 11 Oct 2026 as
   `relax/qm_Rx_converged_step66.json`. The dated notes in sections 1 and 6.6 stand in for the rewriting.
   If the machine is pre-empted while `gr_relax_Rx_cont` runs, the runner tests the finished task `gr_relax_Rx` again
   at its next start; that test asks for state `done` or `crossed`, the file then says `running`, and the runner
   runs `gr_relax_Rx` once more. It continues the same run from the same files, with `queue_vm.json` and therefore
   without `enforce`; the frozen atoms are still put back at the start, and the record names the settings of every
   part (`restarts`). `gr_relax_Rx_cont` waits for it (its second dependency; without it the two tasks would run the
   same relaxation at once, on 14 CPUs) and then finds the run done.
8. **Step 2s** (Seminario terms for the zinc site) is supported by every script; the gate was evaluated with step
   2e. With Seminario terms the same single points are reused: `evaluate` with `"step2": "2s"` in a copy of the
   configuration. **Such a second reading writes under its own names and cannot overwrite the record of today's**
   (plan, 1:09 PM, point 4; check, finding 8): everything that depends on the force-field model (`charges/`,
   `mm_2s_*.json`, `gr_*gate.json`, `gr_*table.csv`, `fit/`, `relax_mm.json`, `relax/mm_*.xyz`, `binding_check.json`)
   goes to `md/qm/gr/reading_2s/`, the scratch files to `mm_2s/` and `evaluate_status_2s.json`; the files of the
   reading "2e" stay where they are, with the mark `reading.json` beside them. Shared, because they do not depend on
   the force field: points, geometries, `qm/`, `d3.json`, the quantum-chemical relaxations. Tested in a scratch copy
   of `md/qm/gr`: a second reading wrote 16 files, all in its own folder, and changed none of the first reading's
   (only the shared summary `relax_qm.json` was written again); `evaluate` with the name "2e" and step 2s, and
   `gr_fit.py gate` or `fit` with another reading pointed at the files of "2e", stop with "must not overwrite" and
   write nothing. The fit of a second reading follows the amended rule of section 6.3. Not tested: step 2s itself,
   which needs `CA2_small_opt.fchk`.
9. **Interface to the system builder.** If a fit is accepted, `fit/CA2_n2_fit.frcmod` must be loaded after
   `CA2_mcpbpy.frcmod` in the builder; if none is, the terms as built stand, with Zn–N2–H at 107.41° (the frcmod
   line `n2h_structure.frcmod` of `gr_mm.py`), not at the 120° of step 2e. Neither line is written in the builder
   yet.

## 10. Files

`md/qm/`: `gr_points.py`, `gr_mm.py` (`run`, `relax`, `binding`), `gr_fit.py` (`gate`, `fit`, `selftest`),
`gr_queue.py` (`run`, `point`, `status`, `complete`, `evaluate`, `d3`, `d3-extract`, `relax`, `relax-summary`),
`resp_equiv.py`. `md/gpu/tasks_gr.json`: the six tasks of the gate on the rented machine; `md/gpu/tasks_gr3.json`:
the continued relaxation from Rx.

`md/qm/gr/`: `points.json`, `points.csv`, `scan_ns_contacts.csv`, `scan_znn_contacts.csv`, `geom/` (38 `.xyz` with
eight decimals, used by both programs, and 38 `.pdb`); `queue_prod.json` (Pi), `queue_vm.json` (machine),
`queue_vm_relax2.json` (machine, the continued relaxation), `queue_vm_check.json`; `d3_reference.json`, `d3.json`;
`charges/`; `mm_2e_main.json`, `mm_2e_asis.json`, `mm_2e_free.json`; `qm/<id>.json`; `check_vm/` (X on the machine
against the Pi); `gr_table.csv`, `gr_gate.json` (and `gr_asis_*`, `gr_free_*`); `binding_check.json`;
`relax_mm.json`, `relax_qm.json`, `relax/` (`mm_<start>.xyz`; `qm_<start>.json` with every step; `qm_<start>.xyz`
when converged, `qm_<start>_last.xyz` before; `qm_Rx_before_crossing.xyz`, the geometry of step 15, written by the
old rule and kept; `fixed_atom_shift.json` and `qm_turns.json`, measured on the Pi on the trajectories fetched from
the machine); after a fit `fit/`, and for an accepted one `mm_2e_fit.json`, `gr_fit_table.csv`,
`gr_fit_gate.json`; `reading.json` (the mark of the reading "2e"); a second reading in `reading_<name>/`;
`selftest/` (`fit_selftest.json`, `fit_rule_cases.json` and `relax_continue_and_enforce_test.json` of this version;
the other files are those of 5 AM); `SHA256SUMS` (scripts, point files, configurations, reference data; the lines
of the three programs changed and of `queue_vm_relax2.json` are of this version).

> **Note added on 11 Oct 2026.** `relax/` also holds `qm_Rx_converged_step66.json`, the finished record of the
> relaxation from Rx (66 steps), and `README.md`, which says which file is which and gives the checksums. The
> converged geometry `qm_Rx.xyz` of that run is on the Pi only
> (`/root/cqc-sfa-md/from_gpu/repo_md/qm/gr/relax/`); `relax/qm_Rx_last.xyz` in git is the geometry of step 23.

`/root/cqc-sfa-md/qm/gr/` on each machine (not for git): `queue.pid`, `status.json`, `queue.log`, `logs/` (PySCF
output per point and per relaxation), `chk/` (orbitals per point), `relax/` (steps and geometries of the
relaxations), `mcpb_2e/`, `mm/`, `charge_step.log`, `mm.log`, `gate.log`, `fit.log`, `points.log`; on the Pi also
`evaluate_status.json`; for a second reading `mm_<name>/` and `evaluate_status_<name>.json`. The watchdog brings the
machine's `md/qm/gr/` to
`/root/cqc-sfa-md/from_gpu/repo_md/qm/gr/` every 15 minutes.
