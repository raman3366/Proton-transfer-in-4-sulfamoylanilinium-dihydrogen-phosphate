# Rehearsal of the bonded zinc model, parameter set 2e, on the rented GPU

Work of 6 Oct 2026, 9:55 AM to 11:40 AM IST. **A REHEARSAL: a test of the pipeline and a first look at the bonded
site. It is NOT a result for the paper**, and nothing of it may be quoted as one, for three reasons that hold for
every number below:

1. **The terms at the ligand nitrogen N2 are not validated.** They are the terms as built before gate G-R
   (`md/PLAN.md`, 6 Oct 2026, 3:02 AM and 8:56 AM). Part A of that gate stood at a root mean square of 3.32
   kcal mol⁻¹ against a criterion of 1.0, with 28 of 38 points, when it was last evaluated before these runs
   (`/root/cqc-sfa-md/qm/gr/queue.log`, 10:25 AM IST; the gate is another agent's work and was not finished).
   These terms decide whether the ring changes its orientation.
2. **One run of 5 ns from each start**, not three of 100 ns.
3. **Set 2e only** (empirical zinc constants, MCPB.py step 2e); the Seminario set 2s did not exist.

I committed nothing, `md/PLAN.md` was not edited, no file of another agent was edited. **Some of these files are in
a commit all the same**: commit `2fd9cdc` of 10:23 AM IST ("GPU tooling: all 14 safety findings …", made by another
session while this work ran) took in `md/build_topology.py` as it then stood, which is its final state,
`md/gpu/tasks_rehearsal_2e.json` and `md/topology/bonded_2e/` without the three `min_check.json`; its message does
not name them. Every number is in a file named in brackets; files without a folder are in `md/test/rehearsal_2e/`
(its `README.md` lists them).

## 1. In short

1. **The builder makes the bonded model from the real MCPB.py output**: `md/build_topology.py --model bonded
   --set 2e --start X|N|R` writes `/root/cqc-sfa-md/systems/bonded_2e/<start>/`. All 24 checks that apply pass
   for X, N and R; a second build gives the same five system files (sha256). The built systems are, atom for
   atom, the restrained systems with other zinc terms: same water, same ions, same box (section 2).
2. **One defect of the bonded path was found and fixed before any run**: tleap built the amide hydrogen of the
   residue after each of the three zinc histidines out of the peptide plane, by 0.4 to 1.2 Å (section 8, D1).
3. **Three complete runs on the L40S** with the protocol of `md/protocol.json` unchanged, 14.2 to 14.5 minutes
   each: TX0, TR0, TN0 in `/root/cqc-sfa-md/runs_test/bonded_2e/` on the machine. Minimisation converged, no
   failure at 2 fs or 4 fs, temperature and density settle, 832 to 845 ns per day in production (section 5).
4. **The zinc site holds in all three**: the four Zn–N distances, the six N–Zn–N angles, no fifth ligand, Glu106
   and water stay off the zinc, the hydrogen bond to Thr199 is there in 1,499 of 1,500 frames. Every criterion
   of gates G-Z and G-V as now worded is met by the numbers of each run (no verdict is given: none of them is
   run X0) (section 6).
5. **TR0, the rotated start, changed to the crystal orientation during the equilibration**, 60 to 70 ps after
   the ligand was released (stage 4), by an inversion at N2 with the hydrogen bond to Thr199 kept; that is
   1.23 ns before production time zero, so its production shows the crystal orientation in all 500 frames.
   **TN0 visited the rotated orientation for three frames (30 ps)** at 1.38 to 1.40 ns and came back. In this
   model, with these N2 terms, the two orientations interconvert within tens of picoseconds (section 7).
   Whether that is the molecule or the unvalidated terms is exactly what gate G-R has to say.
6. **Both analysis programs ran on the three runs and agree** in all 128,018 values compared; one warning per
   run, the same in both, from a known difference between the driver's water count and the analysis of
   version 1.2, which the first program of version 1.3 no longer gives (section 8, D4).
7. **The quantum chemistry on the machine ran 10 % slower** per optimisation step while a GPU task ran (215 s
   before, 238 s during, 210 s after).
8. **Judgement**: builder, driver, queue, copy and both analysis programs work together on the real system; the
   pipeline is ready for run X0 as soon as the final parameters exist. What must change first is in section 10.

## 2. What was built

### 2.1 Changes of the builder [`build_topology_6oct.diff`: 489 lines added, 32 removed]

| Change | Why |
|---|---|
| `--set 2e` or `2s`, required with `--model bonded`. The folder is `bonded_<set>` everywhere (systems, `md/topology/`, and by the name of the folder the runs and the analysis); `build.json` carries `model` (the folder name), `model_kind` and `parameter_set` | plan, 6 Oct 2026, 3:02 AM, point 3: two bonded sets |
| The builder stops if the frcmod of the MCPB.py folder was not made by the step the set names (MCPB.py writes "Created by empirical method" or "… Seminario method" behind every zinc term) | a 2e folder must not be built as 2s |
| Default MCPB.py folder of a set: the first of `/root/cqc-sfa-md/qm/gr/mcpb_<set>` and `md/qm/mcpb` that holds `<group>_tleap.in`; `--mcpb-dir` as before | the real output of step 2e, 3 (equal-charge RESP) and 4 is in `/root/cqc-sfa-md/qm/gr/mcpb_2e/`, where the queue of gate G-R ran `md/qm/after_qm.sh 2e` at 7:07 AM IST; `md/qm/mcpb/` holds step 1 and the potential only |
| **The terms at N2 are one file, `n2_terms.frcmod`, loaded after every file of MCPB.py**: `--n2-frcmod FILE` (what gate G-R settles), or, for set 2e without that option, written by the builder with the terms as built before the gate | the plan fixes the N2 terms for both sets alike; `md/GR.md`, section 9, point 9, had left this line unwritten |
| `--set 2s` without `--n2-frcmod` is refused | the N2 terms of 2s are not Seminario's; give the gate's file or `md/topology/bonded_2e/<start>/n2_terms.frcmod` |
| The renamed histidines get head, tail and connect atoms; their peptide `bond` commands are left out | defect D1, section 8 |
| Read-back for the bonded model (section 2.3), new checks, the comparison with the restrained system | asked for this rehearsal |
| `--note TEXT` → `build.json` (`status_note`) and `STATUS_NOTE.txt` beside the system | every built folder says that it is a rehearsal |
| `--min-check` works for the bonded model (r₀ = the equilibrium lengths of the four zinc bonds); the site report gained the signed volume at N2 and the three angles at N2 | section 3 |

Tested besides the three builds, all in scratch and marked as tests: `--set 2s` on the 2e folder and without a
2s folder is refused (exit code 2); `--model bonded` without `--set` is refused; a build with a **synthetic**
file given by `--n2-frcmod` has exactly the angles and the torsion series of that file, and the check of the
repaired ligand terms follows the file where it replaces one [`look/n2_frcmod_test.json`]; on a copy of the 2e
folder in which only the word "empirical" was replaced by "Seminario", `--set 2s` is refused without
`--n2-frcmod` and builds with `md/topology/bonded_2e/X/n2_terms.frcmod`, all checks passed, and the topology
it writes is byte for byte that of `bonded_2e/X`, as it must be [`look/set_2s_test.json`]. **Not tested: a
real 2s folder** (none exists). One wart, as before this work: a build that is refused has already made its
empty folders (`md/topology/bonded_2s/X` appeared in the first test and was removed).

### 2.2 The three systems [`md/topology/bonded_2e/<start>/build.json`, `openmm_check.json`, `SHA256SUMS`]

Built 10:14 to 10:17 AM IST with `md/build_topology.py` sha256 `b858d1d1fe1b…` (the file as it stands).

| | X | N | R |
|---|---:|---:|---:|
| Atoms | 40,902 | 40,911 | 40,857 |
| Water molecules (crystal + added) | 12,247 (215 + 12,032) | 12,250 (215 + 12,035) | 12,232 (211 + 12,021) |
| Na⁺ / Cl⁻ (33 pairs and one Na⁺) | 34 / 33 | 34 / 33 | 34 / 33 |
| Net charge (e); solute with crystal water | 0.000000; −1.000000 | 0.000000; −1.000000 | 0.000000; −1.000000 |
| Cell edge (Å), density as built (g cm⁻³) | 84.135, 0.9124 | 84.134, 0.9126 | 84.121, 0.9118 |
| tleap: errors, close contacts with a hydrogen, unexpected warnings | 0, 45, 0 | 0, 45, 0 | 0, 44, 0 |
| Checks: passed, not applicable (restrained model only), failed | 24, 5, 0 | 24, 5, 0 | 24, 5, 0 |
| Against the restrained system of the start: atoms, water, ions | the same | the same | the same |
| … cell edge (Å): this, restrained | 84.135219, 84.1352196 | 84.1337263, 84.1337268 | 84.1208941, 84.1208946 |
| … largest coordinate difference of any atom (Å) | 0.0006 | 0.0006 | 0.0006 |
| OpenMM: particles, constraints | 40,902; 38,761 | 40,911; 38,770 | 40,857; 38,716 |
| OpenMM: energy as built (kJ mol⁻¹), finite | −445,219.6 | −445,770.1 | −443,130.2 |
| `system.prmtop`, sha256 begins | a6c57a7c20aa | 98e46802e34a | 7ff894cfda73 |
| `system.inpcrd`, sha256 begins | 27acb79e6645 | 9ce16b4c588d | fe5bc95b9f90 |

All crystal water is kept; solvation (10 Å, truncated octahedron) and ion rules are the code of the restrained
build, and the ions replace the same water molecules as there. The 0.0006 Å are hydrogens of the three
histidines, built from templates with four decimals. `selections.json`, `reference.pdb`, `extra_forces.json`
(empty list) and `build.json` as in the conventions; the CA residue names of `reference.pdb` are HD1, HD2, HE1
where the topology has them (the driver's guard finds 0 of 258 names different).

### 2.3 Parameters read back from `system.prmtop` with ParmEd [`build.json`: `zinc_site_bonded`, `terms_at_N2`, `ligand_parameters`]

The same for X, N and R. Amber form: K (r − r₀)², K (θ − θ₀)², kcal mol⁻¹.

| Term | Value | Source it was compared with |
|---|---|---|
| Zn–N2 (ligand) | 122.1, 1.931 Å | `CA2_mcpbpy.frcmod`, step 2e |
| Zn–NE2 (His94), Zn–NE2 (His96), Zn–ND1 (His119) | 78.0, 2.010 Å; 73.3, 2.021 Å; 74.3, 2.018 Å | the same |
| Six N–Zn–N, each 35.0 | His94–His96 103.82°, His94–His119 113.59°, His94–N2 108.92°, His96–His119 98.41°, His96–N2 114.28°, His119–N2 116.83° | the same (the angles of 6RL9) |
| Six Zn–N–C, each 50.0 | His94: CD2 124.21°, CE1 127.89°; His96: CD2 122.91°, CE1 129.31°; His119: CE1 121.05°, CG 131.22° | the same |
| Zn–N2–S | 50.0, 117.55° | `n2_terms.frcmod` (line copied from step 2e) |
| **Zn–N2–H** | **50.0, 107.407°** (start X; step 2e wrote 120.00°, which is not used) | `n2_terms.frcmod` |
| S–N2–H | 61.55, 118.059° (repaired, `md/LIGAND.md`) | `n2_terms.frcmod` (line copied) |
| H–N2–S–C4 | 0.2624 (n = 1, 180°) and 2.703 (n = 2, 0°); H–N2–S–O zero | `n2_terms.frcmod`; repair frcmod |
| Sum of the three equilibrium angles at N2 | 343.02° (332.4° in the geometry of start X) | |
| Torsions through the zinc | 39, every one zero; no improper with the zinc | |
| Ligand | 18 bonds, 28 angles, 37 torsions, 7 impropers, all equal to `md/ligand/sfa_anion_repaired.prmtop`; the six repaired terms in force under the type Y4 | |
| Zinc, non-bonded | type M1, R_min/2 1.395 Å, ε 0.014917 kcal mol⁻¹ (MCPB.py's frcmod), no 12-6-4 term | |
| Na⁺, Cl⁻ | Joung–Cheatham for TIP3P, as in the restrained systems (`frcmod.ions1lm_126_tip3p` of MCPB.py is not loaded) | |

The gate's own small model ("main" of `md/qm/gr_mm.py`) has Zn–N2–H at 107.406°; the two differ in the third
decimal because each rounds the angle of start X (107.4065°) for itself.

**Charges** (RESP of the large model, equal ligand atoms made equal; the five templates are byte for byte those of
`md/qm/gr/charges/`):

| Residue | Charge (e) | Atoms of interest (e) | Largest difference from the unconstrained fit (e) |
|---|---:|---|---:|
| Zinc (ZN1) | +0.335600 | | 0.044 |
| His94 (HD1) | +0.334448 | NE2 −0.001770 | 0.023 |
| His96 (HD2) | +0.413171 | NE2 +0.047689 | 0.024 |
| His119 (HE1) | +0.447207 | ND1 +0.122659 | 0.074 |
| Ligand (SN1) | −0.530426 | N2 −0.914224, S +0.790275, O1 = O2 −0.488550, N1 −0.954927, C1 +0.464944, C2 = C6 −0.274695, C3 = C5 −0.153047, C4 +0.080205, H2 = H6 +0.146108, H3 = H5 +0.155431, HN11 = HN12 +0.417016, HN21 +0.398775 | 0.049 |
| The five together | +1.000000 | | |

Every atom of the three histidines is listed in `build.json` (`zinc_site_bonded.charges.residues`).

**Pairs across the zinc**, in the topology (excluded-atom list and the 1–4 flags of the torsions) and again in
OpenMM (exceptions of the NonbondedForce): 18 pairs one or two bonds apart are excluded (the zinc with its four
nitrogens and with the eight atoms next to them, among them S and HN21; the six N···N pairs); 39 pairs three
bonds apart are 1–4 pairs with electrostatics divided by 1.2 and Lennard-Jones by 2.0, one term each; 80 pairs
four bonds apart interact in full. **Zn···O1 and Zn···O2 of the sulfonyl group are 1–4 pairs** (Zn–N2–S–O), so
their attraction is scaled; N(His)···O is in full.

## 3. Local checks on the Pi [`md/topology/bonded_2e/<start>/min_check.json`, `logs/mincheck_*.log`]

The topology loads in OpenMM (table above). Capped minimisation with `md/run_md.py` (the protocol with at most
300 L-BFGS iterations per pass, marked as a test; CPU, two threads, nice 10). The Pi was heavily loaded: the job
of X got 0.44 of a core when I looked.

| | X | N | R |
|---|---:|---:|---:|
| Energy before → after (kJ mol⁻¹) | −445,220 → −623,629 | −445,770 → −622,190 | −443,130 → −619,955 |
| Wall time (s) | 1,445 | 501 | 1,735 |
| Zn–N2 before → after (Å); bond r₀ 1.931 | 1.931 → 1.926 | 2.207 → 1.940 | 2.235 → 1.944 |
| Largest Zn–N minus r₀ after (Å); limit 0.10 | +0.017 | +0.012 | +0.040 |
| S–N2 before → after (Å) | 1.596 → 1.570 | 1.678 → 1.581 | 1.677 → 1.583 |
| Ligand RMSD to the crystal pose (Å) | 0 → 0.150 | 0.610 → 0.310 | 3.573 → 3.537 |
| d1, N1 to Phe131 CZ (Å) | 3.897 → 3.815 | 3.356 → 3.596 | 9.585 → 9.759 |
| d2, N1 to Asn67 OD1 (Å) | 7.744 → 7.716 | 7.725 → 7.659 | 2.597 → 2.823 |
| Signed volume at N2 | +0.724 → +0.740 | +0.786 → +0.758 | −0.710 → −0.551 |
| Zn–N2–S–C4 (°) | −117.0 → −126.2 | −107.5 → −125.2 | −32.1 → −37.7 |
| Zinc or ligand heavy atom that moved most (Å); limit | N2 0.212; 0.30 | N1 0.501; 0.50 | O1 0.528; 0.50 |
| Energy lower; four distances within 0.10 Å | passed; passed | passed; passed | passed; passed |
| Zinc and ligand heavy atoms within the limit | passed | **failed** (by 0.001 Å) | **failed** |

The limits are those fixed for the restrained model before its first minimisation (`MIN_CHECK` in the script);
I left them. In N and R the bond pulls Zn–N2 in by 0.27 and 0.29 Å; **both docked starts keep their orientation**
and their configuration at N2. The same after the full minimisation on the GPU [`look/min_site_<run>.json`]:
R 3.573 → 3.568 Å, d1 9.82, d2 2.84, signed volume −0.556, largest move O1 0.525 Å; N 0.610 → 0.380 Å, largest
move N2 0.411 Å; X 0 → 0.198 Å, largest move C6 0.257 Å.

**Order.** The check of X was still running (it needed 24 minutes) when I added the GPU tasks at 10:31 AM IST,
so the first GPU run started before a local minimisation had finished; the build checks and the OpenMM load had
passed. X finished at 10:43 AM, R at 11:12 AM, N at 11:28 AM (its first attempt was ended by the time limit of
my background command and started again). The GPU does the full minimisation in 26 to 35 s.

## 4. Copy to the machine and the tasks [`push/push_rehearsal.log`, `push/push_rehearsal_2e.sha256`]

**I did not run `md/gpu/push_to_vm.sh --yes`.** Its dry run against the machine showed that it would also have
sent other people's files that were being edited: `md/qm/gr_queue.py` (the program of the gate task then
running on the machine), `md/gpu/tasks_gr.json`, `tasks_initial.json`, `GPU_RUN.md`, `push_to_vm.sh`, the
deadman files, `md/analysis/recheck_md.py`, and nine result files of the gate whose copies on the Pi are newer
than the machine's originals. Copied instead, with the same mechanism (rsync `-a --update` through `sudo -n
rsync`, a dry run first that found nothing newer on the machine, then `sha256sum -c` on the machine): the 21
files of `systems/bonded_2e/` (without `build/`; 28.3 MB), `md/gpu/tasks_rehearsal_2e.json`, and **`md/run_md.py`
and `md/protocol.json`**. The machine still had the driver of before the 9:55 AM fix (sha256 `b0f3e2b2…`) and
the protocol with the 3.5 Å count of the hydration check; it now has `6d797f08…` and `decda42f…`, the Pi's.
24 files, verified at 10:21:35 AM IST.

Tasks [`md/gpu/tasks_rehearsal_2e.json`; checked with `queue.py check` on the Pi and on the machine; added
10:31:11 AM IST]: `reh2e_TX0`, `reh2e_TR0`, `reh2e_TN0`, one after another, the ready-made entry of
`md/gpu/GPU_RUN.md` with `RUNS=/root/cqc-sfa-md/runs_test`, `MODEL=bonded_2e`, 5 ns, a time limit of one hour
and a marker file `REHEARSAL_NOT_A_RESULT.txt` in every run folder. All three done at the first attempt: 852,
869 and 869 s (10:31:11 to 10:45:23, to 10:59:52, to 11:14:21 AM IST).

**The quantum chemistry** [`look/qm_speed.json`]: seconds per optimisation step of `qm_opt` (16 threads; a
7-thread gate task beside it throughout): 215.3 s before (13 steps, 196.5 to 231.5), **237.6 s while a GPU task
ran** (10 steps, 221.6 to 252.6), 209.9 s in the one step after. Ten per cent slower; not badly.

**The watchdog on the Pi does not pull `runs_test/`** (its roots are `runs`, `qm`, `systems`, `analysis`). I
pulled the three runs with rsync into the scratch folder (`scripts/pull_run.sh`; 124.5 MB per run, 3.5 to 5
minutes each at the speed of this line) and checked every file of every step against the sha256 the driver
recorded (`scripts/check_sums.py`, the check of `md/RUN.md`, section 4): 29 of 29 files in each run.

## 5. The runs: what the driver did [`look/<run>.json`, `runs/<run>/run_md.log`, `manifest.json`]

CUDA, mixed precision, deterministic forces set, NVIDIA L40S; driver `6d797f08…`, protocol `decda42f…`; the
manifests say `test_only: false` because protocol and equilibration are those of production.

| | TX0 (start X) | TR0 (start R) | TN0 (start N) |
|---|---:|---:|---:|
| Minimisation: energy before → after (kJ mol⁻¹) | −445,219.5 → −659,032.7 | −443,130.2 → −658,616.0 | −445,770.1 → −656,464.7 |
| L-BFGS iterations per pass (limit 5,000); none stopped by it | 1,054, 3,613, 6, 6 | 1,061, 3,607, 6, 7 | 1,017, 2,606, 6, 6 |
| Wall time of the minimisation (s) | 33 | 35 | 26 |
| Stage 1, 2 fs, NVT 50 → 300 K: temperature at the end (K) | 299.7 | 297.4 | 299.8 |
| Stage 2, 2 fs, NPT: mean temperature (K); density at the end (g cm⁻³) | 300.1; 1.0248 | 300.4; 1.0280 | 300.5; 1.0299 |
| Stage 3, 2 fs: mean temperature; density, second half | 300.3; 1.0271 | 300.7; 1.0262 | 299.9; 1.0265 |
| Stage 4, **first at 4 fs**, zinc and ligand released: mean temperature (K); density, second half | 300.7; 1.0263 | 301.2; 1.0260 | 301.2; 1.0273 |
| Stage 5, 4 fs: mean temperature; density, second half | 300.7; 1.0248 | 301.1; 1.0249 | 301.1; 1.0261 |
| Stage 6, 4 fs, free: mean temperature; density, second half | 301.1; 1.0259 | 300.7; 1.0258 | 300.6; 1.0256 |
| Production, 5 ns: temperature, mean and standard deviation (K); expected spread | 300.98, 1.45; 1.46 | 300.95, 1.44; 1.47 | 301.01, 1.46; 1.46 |
| Production: density, mean and standard deviation (g cm⁻³) | 1.02561, 0.00239 | 1.02578, 0.00239 | 1.02555, 0.00233 |
| Production: volume, mean (nm³); slope (nm³ per ns) | 407.84; +0.073 | 407.34; +0.012 | 407.96; +0.049 |
| Speed, ns per day: stages at 2 fs; at 4 fs; production | 455 to 635; 690 to 847; 845.1 | 455 to 643; 859 to 867; 832.3 | 382 to 581; 805 to 840; 838.0 |
| Wall time of all steps (s) | 812 | 811 | 818 |
| Frames: solute, all atoms | 500, 50 | 500, 50 | 500, 50 |

- **No failure at 2 fs or at 4 fs**: every energy of every state log is finite, the highest temperature of any
  log row is 305.9 K (TR0, stage 5).
- **Density**: 1.025 to 1.027 g cm⁻³ over the second half of every stage from stage 2 on; **the box does not drift** in production
  (the slopes are 0.003 to 0.02 % of the volume per ns, inside its fluctuation of 0.95 nm³).
- **The kinetic temperature reads 301.0 K at 4 fs** in all three productions and 300.6 to 301.2 K in the 4 fs
  stages, against 299.9 to 300.7 K in the 2 fs stages. The same offset in three runs is not noise. I did not
  look into it; it is for the owner of the driver to say whether one kelvin at 4 fs is accepted.
- **Speed**: 5 ns in 511 to 519 s with every frame written and the quantum chemistry beside it; the speed test
  without output gave 866 ns per day. 100 ns would take 2.84 to 2.88 hours.

**The hydration check** (count of water oxygens within 2.5 Å of the places of the six ring carbons of the
crystal pose, at the end of stage 3):

| | TX0 | TR0 | TN0 |
|---|---|---|---|
| Ligand, mean distance from those places (Å); limit 2.0 | 0.63 | 3.94 | 0.32 |
| Check | skipped (the ligand is there) | **applies** | skipped |
| Count; nearest water oxygen (Å) | 0; 3.45 | **3**; 1.22 (the three at 1.22, 1.32 and 2.45 Å) | 0; 3.19 |
| Decision; extension steps | go on; 0 | go on (water at the site); 0 | go on; 0 |

It acted as designed in the two branches these runs reached. In TR0 the driver's own counts at the end of each
step were 0 after the minimisation (nearest 4.48 Å), 0 after stage 1 (3.74 Å), 3 after stage 2 (1.96 Å) and 3
after stage 3: the water came into the place of the ring during the 300 ps of stage 2. **The branch that
extends the held part (count zero) did not occur** and has been run only in the driver's own tests on stand-in
systems (`md/RUN.md`, section 6.2). The CA fit of the check deviates by 0.31 Å in TR0.

## 6. The runs: the analysis [`look/<run>.json`; `runs/<run>/analysis/summary.json`; `comparison/recheck_report.md`]

`md/analysis/analyse_md.py` (sha256 `d7fea4c8…`, specification 1.2) and `md/analysis/recheck_md.py`
(`f229865d…`), then the one command of `md/analysis/RECHECK.md`, `recheck_md.py --compare-campaign <scratch>/runs_test
--run-primary`: **the two programs agree on all three runs** (42,662, 42,681 and 42,675 values, 0
disagreements, 0 threshold frames, 0 keys in one program only; the campaign rules agree and both leave the
three runs out of the campaign table, as their names ask). The start comes from `build_json.start` of the
manifest, the model from `build_json.model` (`bonded_2e`); the role is `production`, so **no verdict of gate
G-V is given**; "met" below is `passed` of each part, which the programs compute for every run.

| Quantity (production, 5 ns; 6RL9 or the limit) | TX0 | TR0 | TN0 |
|---|---:|---:|---:|
| Zn–N2, mean ± s.d. (Å); 1.83 to 2.03 | 1.923 ± 0.046 | 1.926 ± 0.049 | 1.930 ± 0.048 |
| Zn–NE2 His94 (Å); 1.91 to 2.12 | 2.019 ± 0.058 | 2.015 ± 0.059 | 2.020 ± 0.059 |
| Zn–NE2 His96 (Å) | 2.043 ± 0.056 | 2.043 ± 0.060 | 2.041 ± 0.065 |
| Zn–ND1 His119 (Å) | 2.038 ± 0.062 | 2.039 ± 0.059 | 2.041 ± 0.062 |
| N2–Zn–His94 (°); 108.9 ± 10 | 105.26 | 105.79 | 105.57 |
| N2–Zn–His96 (°); 114.3 | 116.37 | 115.47 | 115.84 |
| N2–Zn–His119 (°); 116.8 | 118.55 | 118.92 | 118.86 |
| His94–Zn–His96 (°); 103.8 | 102.50 | 102.99 | 102.92 |
| His94–Zn–His119 (°); 113.6 | 114.54 | 114.25 | 114.09 |
| His96–Zn–His119 (°); 98.4 | 98.34 | 98.09 | 98.24 |
| Largest deviation of the six from 6RL9 (°) | 3.64 | 3.11 | 3.33 |
| Sulfonyl oxygen within 2.6 Å of the zinc: solute frames of 500; at most 25 | 1 | 0 | 0 |
| Any oxygen within 2.6 Å: all-atom frames of 50 | 0 | 0 | 0 |
| Nearer sulfonyl oxygen to the zinc, mean; smallest (Å); 2.99 in 6RL9 | 2.980; 2.591 | 2.992; 2.604 | 2.990; 2.619 |
| **Zn to Glu106 OE1**, mean; smallest (Å); 4.00 in 6RL9 | 4.236; 3.552 | 4.174; 3.345 | 4.183; 3.502 |
| Zn to Glu106 OE2, mean; smallest (Å); 5.47 | 6.180; 5.305 | 5.977; 5.168 | 5.915; 4.925 |
| **Zn to the nearest water oxygen**, mean; smallest (Å); 3.69 in 6RL9 (HOH 553) | 3.571; 3.278 | 3.755; 3.201 | 3.793; 3.206 |
| Zn to the nearest protein oxygen, mean; smallest (Å) | 3.852; 3.514 | 3.863; 3.474 | 3.889; 3.359 |
| **N2–H···OG1 (Thr199)**: frames of 500; at least 400 | 500 | 500 | 499 |
| N2···OG1, mean (Å); N2–H···OG1 angle, mean (°) | 2.933; 166.7 | 2.943; 166.3 | 2.944; 166.3 |
| Thr199 N–H···O (sulfonyl); Thr199 OG1–H···OE (Glu106): frames of 500 | 498; 500 | 493; 500 | 497; 500 |
| **Ligand RMSD to the crystal pose**, mean ± s.d.; largest (Å); at most 1.0 | 0.523 ± 0.185; 1.264 | 0.572 ± 0.266; 1.563 | 0.543 ± 0.299; 3.311 |
| **C4–S–N2**, mean ± s.d. (°); 108.8 ± 3 | 109.24 ± 2.59 | 109.42 ± 2.56 | 109.39 ± 2.68 |
| **τ_fold**, mean ± s.d. (°); 132.1 ± 20 | 129.70 ± 9.40 | 132.78 ± 15.66 | 134.60 ± 14.57 |
| S–N2, mean (Å); 1.596. Zn–N2–S, mean (°); 117.6 | 1.573; 114.12 | 1.573; 114.58 | 1.576; 114.17 |
| Criteria GV1, GV2, GV3, GV4 met | all four | all four | all four |
| d1, mean; smallest; largest (Å); midpoint 6.7 | 4.60; 3.36; 6.54 | 4.40; 3.06; 6.40 | 4.67; 3.13; 9.74 |
| d2, mean; smallest; largest (Å); midpoint 5.2 | 7.23; 5.80; 9.37 | 7.86; 5.38; 10.19 | 7.60; 2.74; 9.93 |
| Frames crystal-like, rotated-like, intermediate | 500, 0, 0 | 500, 0, 0 | 496, 3, 1 |
| Verdict of the orientation (window: the whole 5 ns) | STAYED | **CONVERTED**, first passage 0.01 ns | STAYED |
| Sign at N2: frames +, − (start: +, −, +) | 500, 0 | 500, 0 | 498, 2 |
| Signed volume at N2, mean; smallest | 0.697; 0.426 | 0.690; 0.455 | 0.686; −0.431 |
| First production frame: state, sign, ligand RMSD (Å) | C, +, 0.519 | **C, +, 0.974** | C, +, 0.253 |
| Ring turned by 180° (labelling other than identity): frames | 0 | 0 | 0 |
| CA RMSD, mean; largest (Å). Largest CA RMSF (Å), residue | 0.841; 1.083. 1.56, 3 | 0.830; 1.026. 1.64, 154 | 0.786; 0.953. 1.70, 3 |
| Asn67 ND2 nearer to N1 than OD1: frames of 500 | 451 | 206 | 474 |
| Warnings (the same one in both programs; section 8, D4) | 1 | 1 | 1 |

**Orientation through the equilibration** (solute frames every 10 ps; state, sign at N2):

| Step | TX0 | TR0 | TN0 |
|---|---|---|---|
| Stages 1 to 3 (heavy atoms held), 70 frames | C, + in all | **R, − in all**; RMSD 3.42 to 3.99 Å | C, + in all |
| Stage 4 (zinc, ligand, side chains free), 40 frames | C, + in all | **6 frames R, −; then 34 frames C, +** | C, + in all |
| Stage 5, 40 frames | C, + in all | 39 C and 1 intermediate; + in all | C, + in all |
| Stage 6 (free), 50 frames | C, + in all | C, + in all | C, + in all |

## 7. What I saw in the series, said plainly

1. **TR0 leaves the rotated orientation 60 to 70 ps after the ligand is released.** Frames of stage 4
   [`runs/TR0/analysis/timeseries_equil_solute.csv`]: at −1.25 ns (50 ps after the release) still rotated (d1
   10.47 Å, d2 2.83 Å, signed volume −0.30, RMSD 3.93 Å); at −1.24 ns the nitrogen is flat (signed volume
   −0.017, d1 7.84 Å, RMSD 2.76 Å); at −1.23 ns crystal-like (d1 4.87 Å, d2 6.56 Å, signed volume +0.69, RMSD
   1.13 Å); 10 ps later RMSD 0.55 Å and then 0.13 to 0.9 Å. Zn–N2–S–C4 goes from −40° to −131°, Zn···O2 from
   3.8 to 2.9 Å. **It is an inversion at N2 with the hydrogen bond N2–H···OG1 present in every one of the 40
   frames of the stage**, the route the plan's entry of 3:02 AM names. During the held stages the nitrogen was
   already flatter than as built (signed volume −0.71 as built, −0.56 after the minimisation, −0.54 on average
   in stage 3, up to −0.36).
2. **The pocket was not empty when that happened.** Three water oxygens were in the place of the ring at the
   end of stage 3; in stage 4 the ring displaced them (nearest water oxygen to a crystal ring position 0.50 Å
   in the last all-atom frame of stage 3, 2.78 Å in the first of stage 4).
3. **TN0 visits the rotated orientation for three frames**: at 1.38, 1.39 and 1.40 ns d1 is 8.16, 8.01, 9.74 Å,
   d2 4.00, 3.11, 2.74 Å, the signed volume −0.43, +0.23, −0.34, the RMSD 2.50, 2.73, 3.31 Å; at 1.41 ns it is
   crystal-like again (RMSD 1.31 Å). Before it, from 1.27 ns, τ_fold stood at 157° to 176°. One more frame is
   intermediate (3.15 ns, d1 6.97 Å). **So in this model the two orientations exchange in both directions
   within tens of picoseconds**; the rotated one was held for 30 ps (TN0) and for 60 to 70 ps after release
   (TR0), the crystal-like one for nanoseconds. With the N2 terms not validated this says nothing yet about
   the molecule; it says that the answer to the plan's question will be set by these terms, and that it may
   be given inside the equilibration.
4. **In the crystal-like state the aniline nitrogen sits farther from Phe131 CZ than in the crystal**: d1 is
   4.4 to 4.7 Å on average against 3.90 Å, with single frames up to 6.5 Å (TX0: 9 frames above 6.0 Å), close
   to the midpoint of 6.7 Å. The ligand RMSD stays at 0.5 to 0.6 Å. In TX0 the mean of d1 per
   nanosecond goes 4.53, 4.47, 4.61, 4.66, 4.73 Å.
5. **Asn67**: in TR0 ND2 is the atom nearer to N1 in 206 of 500 frames, against 451 and 474 in TX0 and TN0. By
   the reading of `md/analysis/SPEC.md`, section 4, the amide of Asn67 is turned in most frames of TR0, the run
   that started with N1 on OD1. Reported; the states do not use it.
6. **Nothing at the zinc that should not be there.** No water oxygen or protein oxygen came nearer than 3.20 Å
   in any of the 150 all-atom frames of the productions or the 60 of the equilibrations; Glu106 OE1 never
   below 3.35 Å (10 ps frames); a sulfonyl oxygen within 2.6 Å in one frame of 1,500 (2.591 Å, TX0). The
   nearest water oxygen stands 3.2 to 4.0 Å from the zinc, where the crystal has HOH 553 (3.69 Å); in TX0 it is
   one and the same added water molecule in all 50 frames, in the other two runs several molecules take turns
   [`look/nearest_oxygens.json`]. The nearest protein oxygen is OG1 of Thr199 or OE1 of Glu106.
7. **The sulfonamide hydrogen stays on Thr199**: absent in one frame of 1,500 (TN0, 3.60 ns).
8. **No change of orientation while the solute was held** (stages 1 to 3) in any run, and no 180° turn of the
   ring anywhere (labelling "identity" in all 1,500 production and 600 equilibration frames).
9. **No drifting box** (section 5). **The temperature reads one kelvin high at 4 fs** (section 5).
10. **The quantum-chemistry optimisation beside the runs had not converged at step 40** and its gradient had
    risen again (r.m.s. 7.7 × 10⁻⁵ at step 24, 1.15 × 10⁻³ at step 39 [`logs/qm_steps_and_queue_log.txt`]).
    Not part of this work; it bears on whether set 2s exists tonight.

## 8. Defects

| # | Where | What | Done |
|---|---|---|---|
| D1 | `md/build_topology.py`, bonded path | tleap builds hydrogens when it loads the PDB file; the renamed histidines (mol2 units of MCPB.py) had no connect atoms, and their peptide bonds were made afterwards by `bond` commands. The amide hydrogen of **Phe95, Trp97 and Leu120** was therefore built without its C–N bond: C–N–H 99°, 118°, 83° and the improper C–CA–N–H 125°, 157°, 101° instead of 180°; 0.83, 0.41 and 1.17 Å from where the restrained build has them. The box and the ion sites then differed from the restrained system as well (cell edge 84.1340 against 84.1352 Å). The mock tests of `md/TOPOLOGY.md` had the same fault; no check saw it | **fixed** (I may edit this file): the units get head, tail and connect atoms, the peptide `bond` commands are left out, the six "one sided connection" warnings are gone and none is expected any more. New check: every solute atom within 0.01 Å of the restrained build (largest 0.0007 Å). After the GPU minimisation the three hydrogens are in plane (C–N–H 118° to 119°, improper 168° to 177°) [`look/min_site_<run>.json`] |
| D2 | `md/build_topology.py` | no way to name the parameter set, no folder for the second set, the 120° of step 2e for Zn–N2–H would have been built, no place for the terms of gate G-R | **done** (section 2.1) |
| D3 | the machine | it had the driver and the protocol of before 8:56 and 9:55 AM | **copied** the Pi's two files (section 4). No edit of either |
| D4 | `md/analysis/analyse_md.py`, `recheck_md.py` (version 1.2) against the driver | the driver counts water within 2.5 Å of the six ring carbons (protocol of 8:56 AM), the analysis within 3.5 Å of seven positions; its guard compares the two and warns in every run: "the driver counted 0 … the analysis 3" (TX0), "3 … 5" (TR0), "0 … 6" (TN0) | **not touched** (another agent's files; specification 1.3 is that work). `analyse_md.py` and `SPEC.md` changed to 1.3 at 11:23 and 11:26 AM IST, after my analysis. I ran that newer first program (sha256 `f75e6155…`) on the three runs into scratch: **0 warnings**, its count in the place of the ring equals the driver's in all three (0, 3, 0), and every key of gate G-V, of the orientation, of the configuration at N2 and of the fractions has the value of version 1.2 [`look/analysis_v13_against_v12.json`]. `recheck_md.py` was still at 1.2, so the two programs have not been compared at 1.3 |
| D5 | `md/gpu/push_to_vm.sh` | a full push sends every file of `md/` that is newer on the Pi, also files another agent is editing while their task runs on the machine | **not run**; targeted copy instead (section 4). Proposed: an option that takes a list of paths, or a rule that only the owner of a file pushes it |
| D6 | `ops/pi/cqc_gpu_watchdog.py` | `runs_test/` is not pulled, and its estimate lists the three runs as real runs "done, pull pending" with 15 ns planned and 0 done [`/root/cqc-sfa-md/status/gpu_status.json`, 11:20 AM IST]. No alert came of it and the estimate of time left is unaffected | **not touched**. Proposed: treat a run folder outside `runs/` as a test, or add `runs_test` to the roots |
| D7 | `md/build_topology.py --min-check` | the limit of 0.50 Å for the docked starts was set for the restrained model; the bonded model moves N1 by 0.501 Å (N) and O1 by 0.528 Å (R) | reported as failed (section 3); the limit is the plan's owner's |
| D8 | `md/RUN.md`, `md/gpu/GPU_RUN.md`, `md/TOPOLOGY.md` section 7 | they name the model folder `bonded` and `--model bonded --start X`; it is `bonded_2e` and `--model bonded --set 2e --start X` now | a dated note in `md/TOPOLOGY.md`, section 7; the other two are for their owners |

The fix of D1, as the diff has it:

```
+        for unit in chain_units:
+            L += [f"set {unit} head {unit}.1.N", f"set {unit} tail {unit}.1.C",
+                  f"set {unit}.1 connect0 {unit}.1.N", f"set {unit}.1 connect1 {unit}.1.C"]
 …
         for u1, r1, a1, u2, r2, a2 in mcpb["bonds"]:
+            if chain_units and a1 == "C" and a2 == "N" and int(r2) == int(r1) + 1:
+                continue                             # a peptide bond: made by loadpdb through the head and tail atoms
             L.append(f"bond sys.{r1}.{a1} sys.{r2}.{a2}")
```

Not a defect, known beforehand: the DCD header is flushed when a file is opened (driver of 9:55 AM); the runs
used that driver.

## 9. What this shows and what it does not

**Shows**: the builder turns the real MCPB.py output into three systems that pass every check; the driver runs
them on the GPU through the full protocol without a failure; the bonded site with set 2e keeps the crystal zinc
geometry for 5 ns in three runs, where the restrained model had lost it within 4 ps (`md/TOPOLOGY.md`, section
5.4); the copy back and the two analysis programs work on real runs and agree.

**Does not show**:

- anything about which orientation sulfanilamide prefers. The conversion of TR0 and the visit of TN0 come from
  a model whose terms at N2 failed the gate as far as it had been evaluated; the sum of the equilibrium angles
  at N2 is 343°, flatter than any of the three starts (327° to 333°), which may by itself make the inversion easy;
- gate G-V: none of the runs is X0, and X0 has to be run with the final terms;
- the behaviour of set 2s, of a run longer than 5 ns, or of a second run from the same start;
- the extension branch of the hydration check on the real system;
- a resumed run on this model (none of the three was interrupted).

## 10. Is the pipeline ready for run X0? What must change first

Yes, in my judgement, for the parts tested here. Before X0:

1. **The terms at N2 from gate G-R**, then all systems built again with `--n2-frcmod` (for 2e and 2s alike);
   run X0 belongs to those systems, not to the ones of this rehearsal. If the gate keeps the terms as built,
   the file is `md/topology/bonded_2e/X/n2_terms.frcmod`.
2. **The second analysis program at specification 1.3** and the comparison run again; the first program of
   1.3 already reads the three runs without a warning (D4).
3. **For set 2s**: `md/qm/after_qm.sh 2s` in a folder the builder finds (`/root/cqc-sfa-md/qm/gr/mcpb_2s` or
   `md/qm/mcpb`), then `--set 2s --n2-frcmod …`; a first build of real 2s output is still to be made and read.
4. **X0 goes to `runs/bonded_2e/X0`** (or `bonded_2s`), where the watchdog pulls it; `MODEL=bonded_2e` in the
   ready-made entry of `md/gpu/GPU_RUN.md`.
5. **Who pushes what** to the machine (D5).

**For a dated line in `md/PLAN.md`** (I did not edit it):

- The folders of the two bonded sets are `bonded_2e` and `bonded_2s`; there is no folder `bonded`.
- The terms at N2 are one file, loaded last, the same for both sets; as built on 6 Oct 2026: Zn–N2–S 50.0,
  117.55°; Zn–N2–H 50.0, 107.407°; S–N2–H 61.55, 118.059°; H–N2–S–C4 0.2624 (n = 1, 180°) and 2.703 (n = 2, 0°).
- **A run can change its orientation before time zero.** In the rehearsal the rotated start did so 770 ps into
  the equilibration, 70 ps after the release of the ligand. As rule 2 is worded such a run is CONVERTED with a
  first passage at the first frame, and its production holds no frame of its start. The plan should say
  whether the time of a conversion is counted from the release of the ligand (start of stage 4), and what is
  reported for a rotated run without a rotated production frame. If the final terms still let the two
  orientations exchange in tens of picoseconds, the question becomes one of populations, and rule 5 of the
  plan forbids reading those as free energies.
- No production before gate G-R is settled: the rehearsal shows that its terms decide the outcome.
- The limit of 0.50 Å of the minimisation check was exceeded by the bonded model in N (0.501 Å) and R (0.528 Å).

## 11. Files

`md/build_topology.py` (changed); `md/gpu/tasks_rehearsal_2e.json` (new); `md/topology/bonded_2e/X|N|R/`
(`build.json`, `leap.in`, `leap_stage1.in`, `leap.log`, `ions.pdb`, `n2_terms.frcmod`, `selections_table.txt`,
`openmm_check.json`, `min_check.json`, `SHA256SUMS`, `STATUS_NOTE.txt`); `md/TOPOLOGY.md` (a dated note in
section 7); this file; `md/test/rehearsal_2e/` (7 MB, 146 files; its `README.md`).

Not in the repository: `/root/cqc-sfa-md/systems/bonded_2e/` (the systems, with `build/`) and
`/root/cqc-sfa-md/systems/_checks/bonded_2e/` (the local minimisations) on the Pi; on the machine
`/root/cqc-sfa-md/systems/bonded_2e/`, `/root/cqc-sfa-md/runs_test/bonded_2e/TX0|TR0|TN0` (377 MB, every folder
with `REHEARSAL_NOT_A_RESULT.txt`), `/root/cqc-sfa-md/status/push_rehearsal_2e.sha256`,
`/root/cqc-sfa-md/queue/tasks.d/tasks_rehearsal_2e.json`. The trajectories are on the Pi only in the session's
scratch folder (`…/scratchpad/rehearsal/runs_test/`); to keep them:
`rsync -az --exclude '*.chk' -e "ssh -i /root/.ssh/aitf_pi -o IdentitiesOnly=yes" --rsync-path="sudo -n rsync"
ubuntu@<ip>:/root/cqc-sfa-md/runs_test/bonded_2e/ <a folder outside runs/>`. They go when the machine is deleted.

**Before the run folders on the machine are deleted**: the runner tests every finished task again at each of
its starts (`md/gpu/GPU_RUN.md`, section 4), and the test of the three tasks is the manifest of the run. Without
the folders it would run the three rehearsals again after the next restart of the machine (15 minutes of GPU
each). Take `tasks_rehearsal_2e.json` out of the queue first.

**Data pulled by hand**: 285.7 MB (124.6, 80.5 and 80.6 MB received), about US$0.03 at the catalog price the
watchdog uses; its ledger counts only its own pulls.

On the machine I copied files, added three tasks, read status and logs. I stopped, started or cancelled
nothing, touched no systemd unit, installed nothing and ran no ibmcloud command. `md/run_md.py` was compiled
there once to see that the copy was whole (`python -m py_compile`), which leaves a cache file in
`md/__pycache__/`.
