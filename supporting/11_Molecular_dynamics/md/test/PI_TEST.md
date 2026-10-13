# Pipeline test on the Raspberry Pi: builder, driver and analysis on the real system

**This is a test of the pipeline with the restrained-zinc model (free Zn²⁺, four Zn–N distance restraints).
It is not a result for the paper.** The plan dropped this model on 6 Oct 2026, 3:02 AM IST (`md/PLAN.md`, point 4),
while the test was running.

Written 6 Oct 2026, 8:45 AM IST. Work of 1:01 to 8:45 AM IST. Nothing was committed by this test and nothing outside
`md/test/` was changed in the repository. Every number below is in a file under `md/test/results/`, named in brackets.

## 1. In short

1. **The test was cut short from outside.** Run X1 (start X) finished the minimisation with the production limit and
   the whole shortened equilibration (100 ps). At 8:04 AM IST the main session killed it about 57 ps (by the speed; 50 ps are on disk) into the first
   production segment and stopped the chain (`/root/cqc-sfa-md/runs_test/STOPPED_README.txt`: too slow beside the
   quantum-chemical jobs, the model is dropped, the same test takes minutes on the GPU). I did not continue it.
   **Not run: the second half of X1 (150 ps), run R1, the speed test with 1 and 2 threads.**
2. **Builder, driver and analysis work together on the real system of 40,902 atoms** as far as the run went: 7 of
   9 steps complete, all 9 file checks of the complete steps pass, one could not be made (no all-atom frame was
   written before the kill) [`X1_check.json`].
3. **The kill was a real interruption, though not mine and in segment 1, not 2.** It left whole frames only (5 solute
   frames, 5 log rows, an empty all-atom file) and the lock of the dead process. On two copies of the folder
   `--resume` took the lock, marked the killed invocation, verified the last saved state, removed the three cut-off
   files and stopped before repeating the segment; this passes with the driver the run was made with and with the
   driver as it is now (6 of 6 checks each) [`X1_resume_check.json`]. **The repeated segment itself, and so the frame
   count after a resume, was not run** (about 4 hours here).
4. **Speed with 2 threads beside the 4-thread quantum-chemical job: 0.61 to 0.71 ns per day at 4 fs** (1.8 to 2.1
   steps per second). The two threads received 0.76 to 1.0 of a core; all four cores were busy the whole time (load
   average 9 to 13) [`pi_test_speed.json`].
5. **The restrained model does not keep the zinc site** (the plan's reason for dropping it, now seen over 150 ps
   instead of 4 ps). Crystal water 553 is on the zinc after the first 5 ps; the carboxylate of Glu106 joins when the
   side chains are released; in all 5 production frames Glu106 OE1 is 1.86 to 2.01 Å from the zinc and N2–Zn–NE2(His96)
   is 166° on average (114° in 6RL9). The four restrained distances hold; no histidine leaves [`X1_look.json`].
6. **One defect of `md/run_md.py`, still in the current version:** a new DCD file is empty until its first frame,
   and `analyse_md.py` on a run in progress warns "not readable as a DCD file". Seen on the real run; the one-line
   fix is tried on a copy (section 8). Not fixed by me: it does not stop anything.

## 2. What was run

`md/protocol.json` unchanged (sha256 `de0adcef…`), `md/run_md.py` sha256 `45e4cf28…` (both as committed in
`a98bd3f`), system `/root/cqc-sfa-md/systems/restrained/X` (prmtop `f297a0ed…`, inpcrd `abca2bda…`), run seed
154924513 (from the run id), CPU platform, 2 threads, nice 10 [`X1_check.json`, `inputs_sha256`].

```
bash md/test/pi_test_run.sh X     # run_md.py --equil-scale 0.05 --prod-ns 0.2 --segment-ns 0.1 --resume
```

| Planned | Done |
|---|---|
| X1: minimisation, equilibration ×0.05 (100 ps), production 200 ps in 2 segments | minimisation and equilibration complete; segment 1 cut off after 5 frames (50 ps on disk) |
| X1: kill by PID inside segment 2, then `--resume`, frame counts | killed by the main session inside segment 1; recovery checked on copies; segment not repeated |
| R1: equilibration ×0.03 (60 ps), production 100 ps | **not run** |
| `analyse_md.py` on both | run on X1 as it is (version of spec 1.2, sha256 `d7fea4c8…`), exit code 0, 1 warning |
| Speed test with 1 and 2 threads (`pi_test_bench.sh`) | **not run** |

R1 would have had one segment of 100 ps, not two of 50 ps: the driver refuses a segment that is not a whole multiple
of the 100 ps interval of the all-atom file (exit code 2, by design).

Timeline, 6 Oct 2026, IST [`X1_check.json`, `steps`; `/root/cqc-sfa-md/runs_test/restrained/X1.console`]: start
1:04:28 AM; minimisation done 1:32:09; equilibration stages done 1:44:37, 2:37:43, 3:36:00, 4:22:42, 5:10:03,
6:08:39; segment 1 from 6:08:41; killed 8:04:05 AM (exit code 137).

## 3. Speed and load [`pi_test_speed.json`; processor time from `cpu_samples.csv`, one sample a minute]

| Step | Steps | Wall time (s) | Steps per s | ns per day | Cores received | Cores used by everything else | Load average |
|---|---:|---:|---:|---:|---:|---:|---:|
| Minimisation | 4,636 iterations | 1,653 | 0.36 s per iteration | | 0.97 | | |
| eq01 heat, NVT, 2 fs | 2,500 | 744 | 3.36 | 0.58 | 1.00 | 2.98 | 9.0 |
| eq02 density, NPT, 2 fs | 7,500 | 3,183 | 2.36 | 0.41 | 0.89 | 3.11 | 9.8 |
| eq03 loosen, NPT, 2 fs | 7,500 | 3,493 | 2.15 | 0.37 | 0.86 | 3.14 | 9.6 |
| eq04 release, NPT, 4 fs | 5,000 | 2,798 | 1.79 | 0.62 | 0.78 | 3.22 | 11.5 |
| eq05, NPT, 4 fs | 5,000 | 2,838 | 1.76 | 0.61 | 0.78 | 3.22 | 10.4 |
| eq06 free, NPT, 4 fs | 6,250 | 3,512 | 1.78 | 0.62 | 0.76 | 3.24 | 12.0 |
| seg0001, NPT, 4 fs (to its last log row) | 12,500 of 25,000 | 6,085 | 2.05 | 0.71 | 0.86 | 3.0 to 3.3 | 9.1 to 13.1 |

- 2 threads, 40,902 atoms. "Cores received" is the processor time of the process divided by the wall time.
- The others were the quantum-chemical driver with 4 threads (2.0 to 2.3 cores; it ended at 7:07 AM and another
  began at 7:23 AM), docking jobs (vina, 0.3 to 0.8 core) and the agents themselves (`ps`, noted during the run).
- With 5 to 7 busy threads on 4 cores each of the two MD threads got about 0.4 of a core, less than a thread of the
  other jobs got. I did not measure why; two threads that wait for each other are a likely cause.
- Per second of processor time the run made 2.3 to 2.7 steps (pooled 2.45). This is the number that carries over
  to another load.
- The times of segment 1 come from the cumulative speed column of its log (two decimals, about 1 %).
- Minimisation with the production limit of 5,000 iterations per pass: passes of 1,005, 3,619, 6 and 6, none stopped
  by the limit; −447,644 → −661,131 kJ/mol. This had not been run to its end before (`md/TOPOLOGY.md`, section 6).

## 4. Files, frame counts, the interruption

**Complete steps** [`X1_check.json`]. Expected numbers are computed in `pi_test_common.py` from the protocol in the
manifest, independently of the driver.

| Step | Length (ps) | Log rows, solute frames, all-atom frames: expected = found |
|---|---:|---|
| eq01 | 5 | 0, 0, 0 |
| eq02, eq03 | 15 each | 1, 1, 0 |
| eq04, eq05 | 20 each | 2, 2, 0 |
| eq06 | 25 | 2, 2, 0 |
| seg0001 (cut off) | 100 planned | 5, 5, 0 on disk; 10, 10, 1 planned |

| Check | Result |
|---|---|
| Frames and log rows by file size, DCD header and record | passed, 7 of 9 steps complete |
| MDAnalysis reads the same numbers of frames | passed |
| Log times and DCD headers continue from step to step | passed, 6 steps |
| sha256 and size of every file equal the record | passed, 25 files |
| Each step started from the recorded state of the step before | passed, 6 links |
| Last solute frame of a step = its saved state | passed for eq04 and eq05, 3.7 × 10⁻⁶ Å (the other stages do not end on a frame) |
| The cut-off step holds whole frames only | passed: 246,316 bytes = header + 5 frames, 0 bytes beyond |
| Left behind by the kill: the lock and the three files of the cut-off step, nothing else | passed |
| Solute atoms of an all-atom frame = the solute frame of the same time | **not checked: the run holds no all-atom frame** |

With `--equil-scale` the driver counts frames inside each stage, so the scaled equilibration writes 8 solute frames
and no all-atom frame. The all-atom writer was therefore **not exercised by this test**; `md/TOPOLOGY.md` (start N,
output every 1 ps) is the only test of it on the real system.

**The interruption** [`X1_interrupt_snapshot.json`, `X1_resume_check.json`, the two `X1_resume_check_*.console`].
SIGKILL by the main session at 8:04:05 AM IST. `bash md/test/pi_test_after_stop.sh` copies the folder as the kill
left it and runs `run_md.py --stage min --resume` on each copy, which does everything a resume does up to the point of
repeating the segment:

| Check | Driver of the run (`45e4cf28…`) | Driver now (`b0f3e2b2…`) |
|---|---|---|
| Exit code 0, status `incomplete`, nothing simulated | passed | passed |
| Killed invocation marked "no end recorded (killed or crashed)" | passed | passed |
| Lock of the dead process taken and released | passed | passed |
| `full_seg0001.dcd`, `log_seg0001.csv`, `solute_seg0001.dcd` removed | passed | passed |
| Files of the complete steps unchanged (sha256) | passed, 8 files | passed, 8 files |
| "steps complete up to eq06; continuing with seg0001" | passed | passed |

The driver as it is now warned that the script had changed and kept the protocol of the run (no hydration check), as
it is meant to. **Not shown: that the repeated segment has 10 and 1 frames.** The original folder is untouched and
still holds the lock and the five frames.

## 5. Start R

**Not run.** What exists is the minimisation of `md/TOPOLOGY.md` with its test protocol (at most 300 iterations per
pass, 4 of 5 passes stopped by that limit), read again with `pi_test_look.py` [`R_min300_of_TOPOLOGY_look.json`]:

| | Docked start | Minimised (300 per pass) |
|---|---:|---:|
| Zn–N2 (Å); target 1.931 | 2.235 | 1.881 |
| S–N2 (Å) | 1.677 | 1.564 |
| d1, N1 to CZ Phe131 (Å) | 9.585 | 9.567 |
| d2, N1 to OD1 Asn67 (Å) | 2.597 | 2.846 |
| State by the plan's midpoints | R | R |
| Ligand RMSD to the crystal pose (Å) | 3.573 | 3.403 |
| Zn–N2–S–C4 (°) | −32.1 | −50.7 |
| N2–H···OG1 Thr199: distance (Å), angle (°) | 2.585, 140.2 | 2.742, 152.4 |

Ligand heavy-atom RMSD of the minimised start against the docked pose: **0.383 Å** after the fit on CA (0.384 Å
without a fit). The atom that moved most is O1, 0.625 Å. The ring is not thrown around: the zinc distance closes by
0.35 Å and the sulfonyl group turns by 19° about N–S. The minimisation with the production limit and any dynamics
of R with this driver beyond the 4 ps of `md/TOPOLOGY.md` remain untested.

## 6. The run against the gates, as far as 150 ps say anything

Equilibration: 12 snapshots (8 solute frames, 6 saved states with water, two of them at the time of a frame).
Production: 5 solute frames, 10 to 50 ps, no all-atom frame. `pi_test_look.py` is written independently of
`analyse_md.py`; on the 5 production frames the two agree within 0.00005 Å and 0.0005°, with identical states and
hydrogen-bond flags [`X1_look.json`, `comparison_with_analyse_md`].

| Quantity | Built | Minimised | End of equilibration | Production, 5 frames (mean ± s.d.) | Gate window |
|---|---:|---:|---:|---:|---|
| Zn–N2 (Å) | 1.931 | 1.874 | 1.911 | 1.904 ± 0.007 | 1.83–2.03 |
| Zn–NE2 His94 (Å) | 2.010 | 1.987 | 2.026 | 2.031 ± 0.021 | 1.91–2.12 |
| Zn–NE2 His96 (Å) | 2.021 | 1.992 | 2.076 | 2.031 ± 0.013 | 1.91–2.12 |
| Zn–ND1 His119 (Å) | 2.018 | 2.005 | 2.075 | 2.037 ± 0.024 | 1.91–2.12 |
| N2–Zn–NE2 His94 (°) | 108.9 | 106.9 | 94.2 | 95.7 ± 2.7 | 98.9–118.9 |
| N2–Zn–NE2 His96 (°) | 114.3 | 123.7 | 176.3 | 166.3 ± 9.3 | 104.3–124.3 |
| N2–Zn–ND1 His119 (°) | 116.8 | 110.5 | 94.8 | 94.4 ± 6.8 | 106.8–126.8 |
| His94–Zn–His96 (°) | 103.8 | 97.5 | 85.7 | 87.2 ± 0.9 | 93.8–113.8 |
| His94–Zn–His119 (°) | 113.6 | 112.9 | 100.5 | 95.9 ± 2.7 | 103.6–123.6 |
| His96–Zn–His119 (°) | 98.4 | 104.9 | 88.8 | 98.2 ± 4.1 | 88.4–108.4 |
| Zn–OE1 Glu106 (Å) | 4.00 | 3.48 | 1.92 | 1.95 ± 0.05, within 2.6 Å in 5 of 5 | fifth ligand: at most 5 % |
| Zn–nearest sulfonyl O (Å) | 2.99 | 3.13 | 3.94 | 3.63 to 4.04, within 2.6 Å in 0 of 5 | at most 5 % |
| Zn–nearest water O (Å) | 3.69 | 3.77 | 2.08 (water 553) | no all-atom frame | at most 5 % |
| N2–H···OG1 Thr199 | present | present | present | 5 of 5; 2.91 ± 0.07 Å, 159 ± 5° | at least 80 % |
| Ligand RMSD to crystal (Å) | 0.000 | 0.121 | 1.35 | 1.04 ± 0.17 | mean at most 1.0 |
| d1 (Å); crystal side below 6.7 | 3.90 | 4.02 | 6.19 | 5.88 ± 0.43 | |
| d2 (Å); crystal side above 5.2 | 7.74 | 7.27 | 5.68 | 6.29 ± 0.72 | |
| C4–S–N2 (°) | 108.8 | 111.6 | 107.8 | 111.5 ± 2.4 | 105.8–111.8 |
| τ_fold (°) | 132.1 | 129.6 | 134.6 | 145.5 ± 4.8 | 112.1–152.1 |

- **Zinc distances and their spread:** inside the windows, 0.007 to 0.024 Å spread (0.022 Å expected for the spring).
  Zn–N2 is 0.03 Å short of its target. No histidine leaves.
- **N–Zn–N angles:** one of six within 10° of 6RL9. The site is no longer tetrahedral.
- **Fifth ligand: yes, two.** Water 553 of the crystal (3.69 Å from the zinc as built) is at 2.05 Å at the end of
  the first 5 ps, while every heavy atom of protein, zinc and ligand is still held; it is within 2.6 Å in 5 of the 6
  saved states of the equilibration (2.05 to 2.31 Å; 2.66 Å in the other). Glu106 OE1 moves from 4.00 to 3.48 Å
  already in the minimisation, arrives during the stage that releases zinc, ligand and side chains (1.88 Å at its
  end) and stays. At the end of the equilibration the zinc has
  six neighbours: N2 1.91, Glu106 OE1 1.92, His94 2.03, His96 2.08, His119 2.08, water 553 2.08 Å.
- **Thr199:** the ligand hydrogen does not flip away (11 of 12 equilibration snapshots, 5 of 5 production frames);
  Thr199 OG1–H···OE(Glu106) is present in all 17, although OE1 sits on the zinc.
- **Ligand:** it has moved (RMSD 0.75 to 1.21 Å in production) but all 17 snapshots are on the crystal side of both
  midpoints; d2 came within 0.05 Å of its midpoint once (5.25 Å at 10 ps). With a site distorted like this the
  drift says nothing about the force field.
- **Density and temperature.** There are no "last 100 ps" of production; over the 50 ps written: 300.4 ± 1.1 K (5
  rows; 1.5 K expected for one row) and 1.0245 ± 0.0024 g/cm³. The density is flat from 30 ps of the equilibration
  on: 1.024, 1.029, 1.026, 1.026, 1.028, 1.026, 1.023 g/cm³, then 1.028, 1.024, 1.023, 1.022, 1.026 in production
  [`X1_log.csv`]. It is 1.025 and not 1.00 because 12 % of the mass is protein and salt (251,898 amu in all, 12,247
  water molecules; manifest of the run); no drift of the box.
- **No blow-up:** every energy finite, 4 fs stable for 28,750 steps with the restraints at k = 500,000 kJ mol⁻¹ nm⁻².
- `analyse_md.py` on the unfinished run: "0.05 ns seen of 0.2 planned, NOT FINISHED", verdict STAYED (provisional),
  Zn–Glu106 OE1 1.949 Å and within 2.6 Å in 5 of 5 frames, 1 warning (section 8) [`X1_analysis_summary.json`].

Figure: `results/X1_look.png` (all 19 snapshots); table of every plotted value: `results/X1_look.csv`.

## 7. How long the plan's validation run would take on this Pi (estimates)

The plan's equilibration is 675,000 steps (700 ps at 2 fs, 1,300 ps at 4 fs), the 5 ns validation run 1,250,000
steps, the minimisation 1,611 s of processor time [`pi_test_speed.json`, `estimate_validation_run_on_the_pi`].
Time = steps / (2.45 steps per processor-second × cores received).

| Case | ns per day at 4 fs | Equilibration, 2 ns | Validation, 5 ns | All, with minimisation |
|---|---:|---:|---:|---:|
| **2 threads beside the 4-thread quantum-chemical job, as measured here** (0.83 core received) | 0.71 | 92 h | 170 h | **10.9 days** |
| The same, if the two threads received 1.33 cores (an even share per thread) | 1.13 | 57 h | 106 h | 6.8 days |
| The same, if they received 2 cores | 1.70 | 38 h | 71 h | 4.6 days |
| **4 threads, nothing else running**, 15 % lost to threading (my assumption) | 2.9 | 23 h | 42 h | **2.7 days** |
| 4 threads, nothing else running, no loss (upper bound) | 3.4 | 19 h | 35 h | 2.3 days |

Only the first line is a measurement. Nothing was run with 4 threads or on an idle machine: the limit of this test
was 2 threads, and the speed test with 1 and 2 threads was not run. The scoping measured 3.34 to 4.07 ns per day
with 4 threads at 36,805 atoms on the shared machine (`md_scoping_2026-10-05.md`, section 6), which is 3.0 to 3.7
at 40,902 atoms by the atom count; the two 4-thread lines agree with that.

## 8. Defects and remarks

**Defect 1, `md/run_md.py`, `DcdWriter.__init__` (present in `45e4cf28…` and in the current `b0f3e2b2…`).**
The header of a new DCD file stays in the buffer until the first frame, so the file is empty for the first 10 ps
(solute) or 100 ps (all atoms) of a step, 2 % of every 5 ns segment. `analyse_md.py` on a run in progress then
writes a warning into `summary.json` and counts it in `n_warnings`. Seen on X1: "full_seg0001.dcd: not readable as a
DCD file (Reading DCD header failed: premature EOF found in DCD file); taken as holding no frame"
[`X1_analyse_md.console`]. Tried on a copy of the N test of `md/TOPOLOGY.md` with a third, unfinished segment: 2
warnings with empty files, 0 warnings with header-only files. Not fixed (it stops nothing, and X1 was running).
Suggested:

```diff
         head += struct.pack("<4i", 164, 4, n_atoms, 4)
         self.f.write(head)
+        self.f.flush()        # a file without a frame is then a valid DCD file with zero frames
```

**No defect found** in `md/build_topology.py` or `md/analysis/analyse_md.py` by this test. `--resume` accepted the
built system, so its six files are unchanged since 1:00 AM.

**Remarks, not defects.**

1. With `--equil-scale` a stage shorter than an output interval writes nothing (eq01: no log row, its record has
   `mean_over_state_log` null) and the equilibration writes no all-atom frame. The 5-minute pre-flight test of
   `md/RUN.md` uses `--equil-scale 0.05` and so has the same gap; its production part writes all-atom frames.
2. A scaled heating stage does not reach its target: 259.5 K at the end of the 5 ps ramp to 300 K; 300.1 K after
   the next stage.
3. The barostat seems to start with large volume steps in every new step of the chain. In the dense-output N
   test of `md/TOPOLOGY.md` (one barostat attempt per log row in the 4 fs steps) the volume changed in 0 to 3 of
   the 7 to 9 intervals of a step, by 0.04 to 0.96 % each [`pi_test_speed.json`,
   `barostat_volume_changes_in_the_dense_output_test_of_TOPOLOGY`]. I did not test this further. For 5 ns
   segments it should not matter; in X1 the density was at its plateau after 30 ps.
4. Water 553 reaches the zinc with the solute held, because crystal water is free after the minimisation. With
   the bonded model the same should be looked for in the first states of gate G-Z.
5. My own mistakes: I edited `pi_test_run.sh` while a copy of it was running; the running copy printed a syntax
   error after the kill (`X1.launcher.out`), with no effect on the run. The interim version of this file said the
   ring was "in the intermediate range" at the end of the equilibration; by the plan's midpoints it was still on
   the crystal side (d1 6.19 < 6.7, d2 5.68 > 5.2).

## 9. Not done, not verified

- Run R1; segment 2 of X1; the repeated segment after a resume and its frame counts; the speed test with 1 and 2
  threads; the campaign path of `analyse_md.py` on two runs.
- No all-atom frame was written, so nothing about water in production, and the all-atom writer is untested here.
- 5 production frames are 5 numbers. The gate table says what the site looks like, not what a 5 ns run would give.
- Nothing ran on a GPU. Whether a repeated segment is bit for bit the lost one was not tested (2 threads on the
  CPU platform do not repeat themselves, `md/RUN.md`, test 0).
- The estimates of section 7 for 4 threads rest on a rate measured with 2 threads on a loaded machine.

## 10. Files

Scripts, `md/test/`:

| File | What it does |
|---|---|
| `pi_test_run.sh` | starts or continues one test run (X or R) |
| `pi_test_cpu_sampler.sh` | once a minute, the processor time the run and the machine used |
| `pi_test_after_stop.sh` | **what produced the results here**: snapshot of the killed folder, recovery on copies, checks, analysis, look, speed |
| `pi_test_check.py` | `snapshot`, `check [--incomplete]`, `recovery` |
| `pi_test_look.py` | every snapshot from the built start to the last frame; figure; comparison with `analyse_md.py` |
| `pi_test_report.py` | speed, cores received, estimates |
| `pi_test_common.py` | own readers (DCD, inpcrd, saved state) and geometry; nothing imported from the code under test |
| `pi_test_chain.sh`, `pi_test_interrupt.sh`, `pi_test_bench.sh`, `pi_test_analyse.sh` | the plan of the test as it was (kill in segment 2, resume, R1, speed test, analysis); started, stopped at 8:04 AM, **not completed** |

Results, `md/test/results/`: `X1_check.json`, `X1_interrupt_snapshot.json`, `X1_resume_check.json` with two console
files, `X1_look.csv`, `X1_log.csv`, `X1_look.json`, `X1_look.png`, `X1_analysis_*` (copies of the output of
`analyse_md.py`), `R_min300_of_TOPOLOGY_look.*` and `N_min300_of_TOPOLOGY_*` (the test runs of `md/TOPOLOGY.md`, read
again), `pi_test_speed.json`.

Not in git, `/root/cqc-sfa-md/runs_test/` (68 MB): `restrained/X1/` (the run as the kill left it, with `analysis/`),
`restrained/X1.console`, `cpu_samples.csv`, `resume_check/` (the two copies), `driver_of_the_run/` (driver and
protocol of `a98bd3f`), `STOPPED_README.txt` (the main session's note).

Software: OpenMM 8.6.1, Python 3.12.14, NumPy 2.5.3, MDAnalysis 2.10.0 (recorded in every JSON file). No package
was missing.
