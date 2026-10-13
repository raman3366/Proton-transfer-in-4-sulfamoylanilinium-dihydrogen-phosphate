# Second analysis of the MD campaign: the recheck

Stage S6 of `md/PLAN.md` asks that every number is derived twice, by two programs written by different
people, and that the two agree. The first program is `md/analysis/analyse_md.py`. The second is
`md/analysis/recheck_md.py`, written from `md/analysis/SPEC.md` alone, now for **version 1.3** of that text.
This file gives the command that makes the comparison (the gate of S6), what the two programs gave on
6 Oct 2026, where they differ, what each needs, and which sentences of the specification allow two readings.
The file of version 1.2 is in the history of the repository (commit `2fd9cdc`).

Status, 6 Oct 2026, 12:41 PM IST: **on 37 runs the two programs agree in 1,547,845 of 1,547,851 values.**
The six that differ are times of one test run with frames every 0.1 ps, written with three decimals
(section 3.2); it cannot happen with the 10 ps of the campaign. The water count of version 1.3, its settings
and its comparison with the driver's own record agree in every run and in all 43 folders built to probe them.
No run of the campaign exists; nothing here is a result.

## 1. The gate: one command

For the whole campaign:

```
nice -n 10 /root/cqc-sfa-tools/md/bin/python md/analysis/recheck_md.py \
     --compare-campaign /root/cqc-sfa-md/runs --run-primary
```

For one run folder (the validation run, for example):

```
nice -n 10 /root/cqc-sfa-tools/md/bin/python md/analysis/recheck_md.py \
     --compare /root/cqc-sfa-md/runs/bonded/X0 --run-primary
```

| Exit status | Meaning |
|---|---|
| 0 | the two programs agree: every value within the tolerance of section 14 of the specification |
| 1 | at least one disagreement beyond tolerance (a value, a count, a flag, a line, a key, the number of warnings, a campaign rule) |
| 2 | no disagreement found, but something could not be compared: no output of the first program, output written after another version of the specification, or a run the second program cannot analyse |
| 3 | no disagreement, but a verdict, a part of gate G-V, rule 1 or a role differs between the two programs through frames at a threshold or a mean at a limit. Section 14 allows that and asks for both readings; they are in the report and a person has to read them |

Only 0 passes the gate.

What the command does, for every folder below the one named that holds a `solute_segNNNN.dcd`:

1. `recheck_md.py` analyses the run and writes `<run>/analysis_recheck/` (the five series, `summary.json`,
   `recheck_aux.json`). The folder `<run>/analysis/` of the first program is never written by it.
2. With `--run-primary`, `analyse_md.py` is called as a separate process where `<run>/analysis/summary.json`
   is missing, older than a trajectory file or the manifest, or written after another version of the
   specification. Without the option the first program's output must be there.
3. The two outputs are compared as section 14 says: the five CSV files line by line, the blocks of
   `summary.json` key by key (`water_ring_place` and the driver's records among them).
4. For a campaign the rules are evaluated by the second program and compared with
   `<DIR>/analysis_campaign/campaign.json` (`rules_per_model`, the folders left out, the order of the runs).
5. Report: `recheck_report.md` and `.json` in `<DIR>/analysis_campaign_recheck/` (one run:
   `<run>/analysis_recheck/`). The last line printed is the exit status with its meaning.

`--out DIR` puts all output of the second program below `DIR`; give it for a folder that belongs to someone
else. `--no-hash-trajectories` leaves out the checksums of the trajectory files.

Every folder with a trajectory is compared, also a test that lies beside the campaign: today
`/root/cqc-sfa-md/runs` holds `preflight_a` and `preflight_b` (compared, and left out of the campaign table by
both programs) and `bench` (no trajectory, passed over). With `--run-primary` the first program writes its
`analysis/` folder into each of them.

"Agree" means: every number within its tolerance (distances 0.001 Å, angles 0.01°, signed volume 0.001,
times 0.000001 ns), every count, flag, state, labelling, sign, verdict, role and truth value identical, the
site atoms, the cut-off and the source of the water count identical, the same lines in every series, the same
keys in the compared blocks, the same number of warnings. A frame whose value lies within the tolerance of a
threshold may be classed differently by two correct programs. Such frames are listed with both readings and
not counted. What follows from them is listed in the same way only where it can follow: a difference under a
key is put down to threshold frames only if a column that the key depends on has such frames in that run, a
count may differ by at most their number, and the statistics of a water count by at most what these frames
can make. Texts that the specification does not fix are listed where they differ and not counted.

Time. The second program needs Python and NumPy only. On the shared Pi it takes 0.16 ms per solute frame and
10.7 ms per all-atom frame of the real system (40,902 atoms), about 12 s of arithmetic for a run of 100 ns,
before reading and hashing its 1.2 GB (`timing.json`; frames of the real system, not a whole run). On the
rented machine the first program took 2.5 to 2.9 s and the second with the comparison 1.0 to 1.1 s for a
rehearsal run of 5 ns with its equilibration (`vm_rehearsal/vm_rehearsal.json`).

The command itself was run with `--run-primary` on two folders of links (`gate_checks.sh`): a campaign-shaped
copy of six synthetic runs (13,734 solute frames, the first program called seven times, figures drawn; 47 s,
status 0, 790,869 values) and a copy of `/root/cqc-sfa-md/runs` as it is today (the two pre-flight tests; 12 s,
status 0).

## 2. How the second program was made

- **From the text alone.** The changes of version 1.3 (sections 9D, 9E, 13, 14 and appendix A: the count that
  the driver acts on, its settings from the manifest, the wider count under its new name, the comparison with
  the driver's record) were written on top of the second program of version 1.2 without opening
  `analyse_md.py`, `md/run_md.py` or the test scripts of the first program (the review of gate G-R, which was
  read first, quotes one line of the first program of version 1.2, its list of seven site atoms). Before any
  comparison the new numbers were checked against the values that section 16.2 of the specification quotes
  for the runs `place/R1`, `default/R1`, `ext/R1` and `conv/R1`; all were returned.
- **Frozen before the comparison.** Program and output for 34 runs were frozen on 6 Oct 2026, 11:39:47 AM IST
  (sha256 of the program `5ec4c56ec3ded9d6…`; `frozen_before_compare/`). The first program was run on the same
  folders from 11:40:05 AM, the first comparison was made at 11:41:51 AM, and the code of the first program
  was read only after it and after the first round of edge cases (from 11:50 AM).
- **What changed afterwards** is in section 4. For the 34 runs the output is the same before and after: the
  168 series and `summary.json` files frozen are identical to those of the final program (series as text;
  `summary.json` under every key except the checksum of the program, its list of constants, the checksums of
  the input files and the folder name).
- **Another route** wherever the specification leaves one: own DCD and PDB readers instead of MDAnalysis
  (frames counted from the file size); the fit by Horn's quaternion method; angles as
  2 atan2(|û − v̂|, |û + v̂|); torsions from the two plane normals; minimum images through fractional
  coordinates, with a search of the 27 neighbouring cells where a distance is not below half the smallest
  width of the cell; the two water counts as minima over one table of distances per crystal position.
- **Checks of its own.** Geometry against plain formulas within 7 × 10⁻¹³ (`unit/unit_geometry.json`). For
  version 1.3 (`unit/unit_v13.json`): the settings of the count in 32 manifests typed by hand, 32 as expected;
  the minimum image in a truncated octahedron against a search over 1331 cells for 1538 points up to three
  cells away, largest deviation 4 × 10⁻¹⁴ Å; the crystal positions carried into a rotated and shifted frame,
  7 × 10⁻¹⁴ Å. (The first form of the minimum-image check searched too few cells and reported a false
  deviation of 105 Å; the check was wrong, not the routine.)
- **The generator's planted values against the second program.** `make_test_trajectory.py --compare` was run
  on the output of the second program in place of the first: **7339 of 7339 planted values are returned**
  (`planted_check.sh`). Two adaptations of form were needed and none of a value: an alias for one entry of
  `run.segment_overlap_check.boundaries[]` and the list `skipped` of `campaign.json` as folder paths; the
  specification names neither form.
- **The comparison itself was tested**: 75 cases, 75 as expected (`edge/compare_selftest.json`). Copies of the
  first program's output with one planted change each (a distance by 0.002 Å, a flag, a water count, the
  source or the order of the site atoms, a truth value of a driver record, a missing key or line, and so on)
  are reported as disagreements; changes within tolerance are not. A planted water oxygen on the cut-off is
  listed as a threshold frame; a verdict changed beside it is still a disagreement, because a water count
  cannot explain a verdict. A verdict changed beside a state at the midpoint, a gate part with its mean at
  the limit and rule 1 with its mean at 2.0 Å give status 3. The five exit statuses were run through the
  command's own code.

## 3. What was found on 6 Oct 2026

Programs: `recheck_md.py` sha256 `5b4c0595c067ff57…`, `analyse_md.py` sha256 `f75e6155de700700…` (the one
named in section 16 of the specification), `SPEC.md` sha256 `6b5b5f3215971690…`, driver `md/run_md.py`
`6d797f080ff75679…`, `md/protocol.json` `decda42f5340d06f…`. Every number below is in `recheck_today.json` of
the scratch folder (section 8); `run_all.sh` there writes it again (last run 12:33 to 12:39 PM IST).

### 3.1 The 37 runs

| Set | Runs | Agree | Solute / all-atom frames | Equilibration frames | Values compared | Disagreements | Threshold items | Keys in one program only | Largest difference: series (Å, °), `summary.json` (Å, °, volume) | Warning counts equal |
|---|---|---|---|---|---|---|---|---|---|---|
| synthetic runs of `make_test_trajectory.py` | 23 | 23 | 22,084 / 2,203 | 1,058 / 105 | 1,372,322 | 0 | 0 | 0 | 0.0001, 0.001; 0.000013, 0.000325, 0.000002 | 23 of 23 |
| test runs of the real driver, on the Pi | 11 | 10 | 333 / 136 | 178 / 102 | 46,962 | 6 | 0 | 0 | 0.0001, 0.000; 0.000005, 0.000000, 0.000000 | 11 of 11 |
| rehearsal runs on the rented machine | 3 | 3 | 1,500 / 150 | 600 / 60 | 128,567 | 0 | 0 | 0 | 0.0001, 0.000; 0.000003, 0.000000, 0.000000 | 3 of 3 |
| **all** | **37** | **36** | 23,917 / 2,489 | 1,836 / 267 | **1,547,851** | **6** | 0 | 0 | | 37 of 37 |

The differences in the series are one unit of the last decimal written. The campaign rules of the synthetic
campaign agree (23 values; outcome `validation_failed`, run `X0` left out by both).

The runs, one by one:

| Run | What it is | Frames, solute / all-atom: production; equilibration | Values | Disagreements | Driver records: compared, count agrees, nearest agrees | Warnings (first / second) |
|---|---|---|---|---|---|---|
| test/X1, X2, X3, R1 | synthetic, version 1.0 | 3000 / 300 each | 172,265; 172,237; 172,265; 172,244 | 0 | – | 0 / 0 |
| test/N1, test/X0, test/X9 | synthetic, version 1.0 and 1.1 | 1234 / 123; 500 / 50; 250 / 20 | 71,801; 30,057; 15,767 | 0 | – | 0 / 0; 0 / 0; 4 / 4 |
| angle, csn, edge, sulfonyl, tau, water /X0 | synthetic validation runs | 500 / 50 each | 31,083 each | 0 | – | 0 / 0 |
| long/X0, short/X0, glu/X0, gate/R1 | synthetic | 600 / 60; 300 / 30; 500 / 50; 400 / 40 | 36,963; 19,323; 31,089; 25,210 | 0 | – | 2 / 2; 0 / 0; 0 / 0; 0 / 0 |
| conv/R1, flip/X1, scaled/N1 | synthetic, with equilibration | conv 2200 / 220 and 200 / 20; flip 300 / 30 and 200 / 20; scaled 200 / 20 and 8 / 0 | 142,966; 31,115; 14,032 | 0 | – | 0 / 0 |
| ext/R1 | synthetic, two extension steps, settings of before 8:01 AM in its manifest | 200 / 20; 300 / 30 | 31,058 | 0 | 3 of 3, 3, 3 | 1 / 1 |
| place/R1 | synthetic, three extension steps, settings of the plan | 200 / 20; 350 / 35 | 33,983 | 0 | 4 of 4, 4, 4 | 0 / 0 |
| default/R1 | synthetic, manifest of the plain keys | 200 / 20 | 13,449 | 0 | – | 0 / 0 |
| `runs_test/restrained/X1` | Pi pipeline test: hCA II, restrained-zinc model, production cut off | 5 / 0; 8 / 0 | 2,549 | 0 | – | 1 / 1 |
| `runs/preflight_a`, `preflight_b` | pre-flight tests of the rented machine, pulled by the watchdog: hCA II, restrained-zinc model, shortened equilibration, hydration check with six atoms and 3.5 Å | 40 / 4; 8 / 0 each | 4,622; 4,627 | 0 | 0 of 1 each | 1 / 1 |
| driver as it stands, `dry` and `wet` | two runs of `run_md.py` `6d797f08…` on its stand-in systems, made by the maintainer of the specification; every step ends on an all-atom frame; 4 extension steps in `dry` | 4 / 4; 40 / 40 and 20 / 20 | 4,428; 2,849 | 0 | 5 of 5, 5, 5; 1 of 1, 1, 1 | 0 / 0 |
| driver of 5:28 AM: `dry`, `arrive`, `wet`, `skip`, `late` | test runs of the driver's maintainer, seven site atoms and 3.5 Å (11.2 Å in `arrive`); no step ends on an all-atom frame | 4 / 2; 8 to 16 / 3 to 7 | 2,450; 2,292; 1,818 × 3 | 0 | 0 of 5, 4, 1, 1, 1 | 7 / 7; 1 / 1; 7 / 7 × 3 |
| `test_driver/runs/long` | driver test: 22 segments, frames every 0.1 and 0.2 ps | 220 / 110; 40 / 20 | 17,691 | **6** | – | 0 / 0 |
| `runs_test/bonded_2e/TX0`, `TR0`, `TN0` | rehearsal on the rented machine: hCA II, bonded set 2e, GPU, 5 ns each, full equilibration, driver and protocol as they stand | 500 / 50; 200 / 20 each | 42,845; 42,864; 42,858 | 0 | 1 of 1, 1, 1 each | 0 / 0 |

None of these is a run of the campaign: synthetic frames, a system that is not hCA II, the zinc model that
the plan dropped, and a rehearsal marked on the machine as no result (parameter set not validated by gate
G-R). They test the two programs and the reading of the driver's files.

### 3.2 The one disagreement beyond tolerance

Six values, all `time_ns` of the solute series of `test_driver/runs/long`, each by 0.001 ns. That run has a
frame every 0.1 ps, and a frame at 4.5 ps cannot be written in ns with three decimals: the first program
writes `0.005` (it multiplies the frame number by the interval in ns), the second `0.004` (it divides the
time in ps by 1000). Both are rounding noise at a half-way case; neither is wrong by the text. Every other
value of that run agrees (largest difference 0.0001 Å in a series, one unit of the last decimal, and
0.000001 Å in `summary.json`). With the 10 ps and 100 ps of the campaign every
time is a whole number of the last decimal and the two programs write the same text. The sentence that
removes it is point 1 of section 6.

### 3.3 The count that the driver acts on, set beside the driver's own record

The 37 runs hold 30 records of the driver's hydration check. **16 are compared and all 16 agree in the count
and in the nearest distance** in both programs: 7 in the synthetic runs, 6 in the two runs of the driver as it
stands, 3 in the rehearsal. The other 14 are not compared because the step does not end on an all-atom frame
(a shortened equilibration: 12 in the test runs of 5:28 AM, 2 in the pre-flight tests); both programs write
`compared` false there.

The rehearsal is the first case with the real system, the real driver and the protocol of the plan:

| Run | Driver, end of `eq03`: count, nearest water oxygen (Å), `applies`, decision | All-atom frame at 700 ps, `n_water_ring_place` and nearest (Å): first program | second program | `compared`, `count_agrees`, `nearest_agrees` (both) |
|---|---|---|---|---|
| TX0 | 0, 3.451, false, go on | 0, 3.451016 | 0, 3.451017 | true, true, true |
| TR0 | 3, 1.216, true, go on | 3, 1.215560 | 3, 1.215559 | true, true, true |
| TN0 | 0, 3.188, false, go on | 0, 3.188018 | 0, 3.188020 | true, true, true |

These are numbers of a rehearsal and no result. In the same three runs both programs find the settings in the
manifest (six ring carbons, 2.5 Å, the defaults), the step table in the manifest (2000 ps, no extension step),
every step file with the frames expected, the DCD headers contiguous, the checksums of the system files right
and Zn–N2 of the two trajectories equal; neither gives a warning.

### 3.4 Run folders built to probe the text: 56 of 71 agree

71 folders and 3 small campaigns (`edge/make_edge_cases.py`; results `edge/edge_results.json`; 3,442,681
values compared in the 70 folders that both programs analyse). In all 70 no value of a series differs beyond
tolerance; five water counts of one folder are threshold frames (below).

**Version 1.3: 43 of 43 agree.** Run `place/R1` with other entries in its manifest:

| What was changed | Cases | Both programs |
|---|---|---|
| no `settings`, no `hydration_check`, `settings` null, `hydration_check` a text | 4 | the defaults, no warning, the four records still agree |
| `settings` a list or an empty object; a site atom the ligand lacks (C9), a hydrogen, small letters, a name twice, no name, names as one text; cut-off zero, negative, a text, true, missing | 13 | the defaults, one warning |
| the six ring carbons in another order; cut-off 2.5000000001 | 2 | from the manifest, `settings_are_the_defaults` true, no warning |
| cut-off 2.50000001; a whole number 3; six atoms and 3.5 Å; seven atoms and 3.5 Å; O1 and S with 4.0 Å | 5 | from the manifest, not the defaults, one warning, and one warning for each record that no longer agrees (0, 3, 3, 4, 4) |
| cut-off 2.4 Å, which lies on a planted water oxygen | 1 | the two programs class that water differently in 5 all-atom frames; listed as threshold frames with what follows (6 items); no disagreement |
| the driver's count as 1.0, as a text, as true | 3 | no record for that step |
| the count wrong; the driver's distance set to 2.4015, 2.402 or 2.403 Å beside the frame's 2.399999 Å; both wrong; no distance; the distance a text; no `applies`; a record on `eq05` | 9 | `count_agrees` and `nearest_agrees` as expected (2.4015 Å agrees, 2.402 and 2.403 Å do not; none without a number), one warning for a step that disagrees |
| a step that does not end on an all-atom frame; its all-atom file cut to the header or missing; no `full.pdb` | 4 | `compared` false |
| a tag listed twice in the manifest; `ext/R1` with one count changed from 1 to 2 | 2 | the last entry counts; `count_agrees` true, true, false |

**The 15 folders and 3 campaigns that differ.** In 9 of the 15 folders the only difference is the number of
warnings or one key, and in none does a value of a series differ; none of the 18 cases can occur in a complete
run of the campaign made by the driver.

| Case | First program | Second program | Whose reading |
|---|---|---|---|
| manifest `replica` 0.5, `"0"` or `false` with start X (three folders) | role `validation`, gate G-V with a verdict | role `production`, no verdict | first program departs from the text: section 1 says a manifest value that is not an integer never gives the role `validation` |
| manifest `replica` 0.0 | `validation` | `production` | open (section 6, point 14) |
| a `manifest.json` that is JSON but a list | stops with a Python error, no output | analyses the run as without a manifest | first program: a fault; the text covers only a file that "cannot be read as JSON" (point 9) |
| the all-atom file of a segment that the manifest lists is missing | 0 warnings | 1 warning (manifest and files) | first program: the guard "manifest" is not applied to an all-atom file that is absent |
| files of a step that is in no step table: `eq07` in both kinds; the six files of three extension steps under a manifest that cannot be read; a step record without its end | 2, 9 and 2 warnings (one per kind of file) | 1, 10 and 1 (one per step) | open (point 11) |
| a segment file that is empty (0 bytes), two cases | 1 warning more (not a DCD file) | no such warning | open (point 11) |
| an all-atom segment file with a header only and no other all-atom frame | 0 warnings | 1 warning (no all-atom frame was read) | open (point 11) |
| trajectories without a unit cell | 2 warnings | 0; `frames_without_box` is recorded | open (point 11) |
| all equilibration files cut to their headers | `equilibration.available` false | true | open (point 15) |
| planned length 5.005 ns, 500 frames | complete | not complete (500 < 500.5) | the text compares N with the quotient; the first program rounds it first (point 16) |
| campaign with runs X2 and X10 | order X10, X2 | order X2, X10 | first program: the replica is sorted as text |
| campaign of two R runs, both CONVERTED, no X or N run | `mixed` | `rotated_converts` | open (point 17) |
| campaign with four runs of start X, all complete | `final` false | `final` true | open (point 17) |

The 13 older folders that agree in everything, warnings included: the Pi test run without its manifest (18
warnings each); equilibration with extension steps, and a shortened one, read without the manifest (8 and 17);
a stage file cut to its header and one emptied (3); the driver's count changed (2); no all-atom trajectory,
and all-atom files without `full.pdb` (1); `--planned-ns 2` on a validation run; `--start R` on an X run;
folder `X0` with replica 2 in the manifest; frame intervals of 20 / 200 ps and 4 / 40 ps on the command line
(3 and 39); seven saved states of the Pi test run packed as all-atom frames of the real system.

### 3.5 What was done on the rented machine

The watchdog pulls `runs/` and not `runs_test/`, so the rehearsal runs were compared in place
(`vm_rehearsal.sh`, 12:39 to 12:40 PM IST with the final program; two earlier passes, at 12:02 and 12:12 PM
with the comparison as it stood then, gave the same counts). Each time the two programs were copied to
`/tmp/cqc_recheck` on the machine, each run was analysed by both with one thread at low priority (4 s per
run), all output went to that folder, the three reports and six `summary.json` files (57 kB packed, 0.4 MB
unpacked) were fetched, and the folder was removed. No task was added to the queue, no file of the
repository, of the run folders or of a service there was written, and the series of the rehearsal were not
fetched. The machine held the driver and the protocol that the specification names, by their checksums.

## 4. What was corrected in the second program after the first comparison

Three points of the analysis, none of which changes a number of the 37 runs, and four of the comparison:

1. `driver_hydration_check.nearest_water_O_A` is copied from the manifest as it stands. The second program
   first wrote none for a value that is no number; section 9E says "copied". (Case: a distance written as a
   text.)
2. `nearest_agrees` is none when the frame has no water oxygen at all, since there is no distance to compare.
   The text is silent; it cannot occur in the campaign (`edge/nowater/`: both programs now write the same).
3. A tag listed twice in `equilibration.stages` is one step and its last entry counts. The second program
   would have read the file twice. The text is silent; the first program reads it the same way.
4. Comparison: when the two programs class a water oxygen at the cut-off differently, the mean, standard
   deviation and block means of that count differ too. The comparison first reported these as disagreements;
   it now lists them with the frames, as long as they stay within what the listed frames can make, and
   reports a disagreement beyond that.
5. Comparison: a difference under a key is put down to threshold frames only if a column that the key depends
   on has such frames (before: if any column of the run had one). A verdict, a gate part, rule 1 or a role
   that differs in this way gives exit status 3 instead of 0. Rule 1 with its mean within 0.001 Å of 2.0 Å is
   treated like a gate part at its limit.
6. Comparison: output of the first program written after another version of the specification is not compared
   (exit status 2); with `--run-primary` it is made again.
7. Comparison of a campaign: the folders left out are compared as paths relative to the campaign folder.
   Before, a run folder directly below the campaign folder (`runs/preflight_a` today) was named differently
   for the two programs and reported as a disagreement of the campaign. Found by running the command on a
   linked copy of `/root/cqc-sfa-md/runs` as it is today: exit status 1 before, 0 after. With the layout
   `<campaign>/<model>/<run>` alone it did not show.

The correction of version 1.2 stands: the window of the verdict, the dwell and the blocks follow from 20 ns,
1 ns and 10 ns when the frame interval is not the default.

## 5. What the first program needs

`md/analysis/analyse_md.py`, sha256 `f75e6155de700700…`. **For version 1.3 itself nothing**: the two counts,
the settings, the step table with extension steps and the comparison with the driver's record are the same in
both programs in every case tried. Four points remain, three of them from version 1.2; none touches a run of
the campaign as planned:

1. Line 658, `rep_int = int(replica)`: a manifest replica of 0.5, `"0"` or `false` becomes 0 and a run of
   start X becomes the validation run, with a verdict of gate G-V. Take the role `validation` only for a JSON
   integer 0 (not a truth value, not a text, not a number with a fraction).
2. Line 1952, the campaign is sorted by `str(replica)`: replica 10 stands before replica 2. Sort by the
   number. No effect on replicas 1 to 3.
3. Lines 975 to 989, `manifest_segment_check` goes through the files that exist and adds only a missing
   solute file: a segment listed with all-atom frames whose `full_segNNNN.dcd` is missing is not reported.
4. New. Line 624 onward: a `manifest.json` whose top level is not an object (a list, a number) stops the
   program with `AttributeError: 'list' object has no attribute 'get'`. Treat it as section 1 treats a file
   that cannot be read as JSON: one warning, then a manifest without any key.

Line 1326 rounds the planned number of frames; the text compares N with the quotient itself. It shows only
when the planned length is not a whole number of frames, and the text may as well be changed (point 16).

## 6. Sentences of the specification that allow two readings

None changes a number of a complete run of the campaign as planned (folders `X0`, `X1` to `R3`, the driver's
manifest, 10 ps and 100 ps, files whole). "Both" means that the two programs read it the same way today and
the text does not say it. Points 1 to 10 came up with version 1.3; points 11 to 18 are open since version 1.2
and were run again today.

| # | Where | What is open | First program | Second program | Sentence that removes it |
|---|---|---|---|---|---|
| 1 | section 13 | `time_ns` of the production series has three decimals; a frame that is not at a whole number of ps falls on a half-way case (6 disagreements today, section 3.2) | 0.005 for 4.5 ps | 0.004 | "`time_ns` is written with six decimals in all four series." |
| 2 | section 14 | counts "may differ by at most their number" at threshold frames; the statistics of such a count are not named | writes what it counts | the same; the comparison now bounds them | "The mean, standard deviation and block means of a count, and its mean in the table row, may then differ by what the listed frames can make: at most D / n for a mean over n frames, D being the sum of the differences in the frames listed." |
| 3 | section 9E | `nearest_water_O_A` is "copied"; what if the manifest holds no number there | as it stands | as it stands (since section 4, point 1) | "`nearest_water_O_A` is the value of the manifest as it stands, a number with six decimals, any other value unchanged, none if the key is missing." |
| 4 | section 9E | `nearest_agrees` uses `d_water_ring_place_min_A` "of that frame": as computed, or as written with six decimals | as written | as computed | "… evaluated on the distance as computed, before it is rounded." The two can differ only within 5 × 10⁻⁷ Å of 0.002 Å, inside the 0.0001 Å that section 14 already lists. |
| 5 | sections 9D, 9E | a system without any atom in `water_O` | counts 0, distances `nan` in the series and none in `summary.json`, `nearest_agrees` none | the same (since section 4, point 2); run today on `place/R1` with such a system: the two agree in all 33,983 values | "If `water_O` is empty the two counts are 0, the two distances are `nan` in the series and none in `summary.json`, and `nearest_agrees` is none." |
| 6 | section 9E | "a whole number (a JSON integer)" for `n_water_O` | both: 1.0, `"1"` and `true` are not | | "A JSON integer is a number written without fraction and exponent; 1.0, a text and a truth value are not." The same sentence settles point 14. |
| 7 | section 9D | `hydration_check.settings` that is not null and not an object | both: warning (i), the defaults | | "A value of `hydration_check.settings` other than null that is not an object fails (a) and (b)." |
| 8 | section 9E | a tag listed twice in `equilibration.stages`; an entry without a start or an end | both: the last entry counts; an entry without both numbers is no step, and its files are files of no step | | "If a tag occurs more than once, the last entry counts. An entry without both numbers is no step." |
| 9 | section 1 | a `manifest.json` that is JSON but not an object | stops with an error | no key, no warning | "A `manifest.json` that cannot be read as JSON, or whose top level is not an object, gives one warning and is then treated as a manifest without any key." |
| 10 | section 9E, guard "header" | the driver rounds the first step of the header to a whole number of time steps, so a right time axis can get the warning when a step does not start on one | both warn alike (6 warnings in each of four test runs) | | "… each within 0.1 % of Δt or half a time step of the file, whichever is larger." Not needed for the campaign. |
| 11 | sections 1, 9E, 16.4 | `table_row.n_warnings` is compared, but there is no closed list of warnings with how many each gives. Open: (a) files of a step in no step table: one warning per kind of file or per step; (b) a file that is no DCD file; (c) no all-atom frame read although files are there; (d) frames without a unit cell; (e) an all-atom file of a listed segment missing; (f) all-atom frames beyond the last solute frame; (g) one warning per repeated boundary or per run | (a) per kind; (b) warns; (c) no; (d) warns, per trajectory; (e) no; (f) warns; (g) per boundary | (a) per step; (b) no; (c) warns; (d) no; (e) warns; (f) no; (g) per run | A table "every warning, its condition, how many": (a) one per kind of file, naming the files; (b) one per file; (c) one whenever no all-atom frame was read; (d) one per trajectory; (e) part of the one manifest warning; (f) one; (g) one per run. The three warnings of version 1.3 are already written this way. |
| 12 | sections 4, 9A, 10 | W, D and the block length are given as frames and as times; with another frame interval they differ | both: times (20 ns, 1 ns, 10 ns) | | "W, D and the block length are 20 ns, 1 ns and 10 ns: round(20000 / Δt), round(1000 / Δt) and round(10000 / Δt of the trajectory) frames, each at least 1." |
| 13 | section 1, guard "contiguous headers" | a header that starts at one interval "in every file" tells nothing; what if only a later file does | per boundary | every file must follow the one before it | "Each boundary is judged on its own: contiguous, or the second file starts at one interval (tells nothing), or neither." |
| 14 | section 1 | which manifest value of `replica` is "an integer" | everything `int()` accepts | a JSON integer | the sentence of point 6 |
| 15 | section 9E | `equilibration.available`: files are there, or frames were read | frames were read | a file of a step of the table is there | "`available` is true if at least one equilibration frame was read." |
| 16 | section 1 | "complete if N ≥ planned length / Δt" when the quotient is not whole | rounds the quotient | compares with the quotient | "… if N ≥ the planned length / Δt rounded to the nearest whole number of frames." |
| 17 | section 12 | `rotated_converts` without a run of start X or N; `final` with more than three runs of a start; "ordered by … replica" | needs one such run; exactly three of each; as text | "every run" of none is true; at least three complete of each; as a number | "… and the model has at least one run of start X or N, all of them STAYED." "`final` is true if the model has exactly the nine runs X1 to X3, N1 to N3, R1 to R3 and all are complete." "… and replica, as a number." |
| 18 | sections 10, 13 | texts that are compared key by key but not defined: `stats.<column>.kind` for a distance and an angle, `n2_configuration.order`, the form of the entries of `skipped` | `A`, `deg`; a formula; folder paths | `dist`, `ang`; "Zn, HN21, S"; folder and reason | "`kind` is one of A, deg, tors, vol, count. An entry of `skipped` is the path of the folder." (or: these texts are not compared; the second program lists them and does not count them) |

Settled by version 1.3: the unit of the step length in a DCD header and the form of the cell angles
(appendix A.2). Run today: points 1 to 3, 5 to 10, 11 (a to e), 12 and 14 to 18. Read in the two programs
and not run: points 4, 11 (f, g) and 13.

**The seven points the maintainer of the specification left for this program.** (1) Warnings per file, step or
run: point 11. (2) Empty `water_O`: point 5; both programs now write the same. (3) Site atoms and cut-off are
taken from the top-level settings only: read the same way; a step record with other settings shows as a
disagreeing count, as the text says. (4) Guard "header" on a shortened test run: both programs give the same
six warnings; point 10. (5) Three decimals of `time_ns`: point 1; it is the one disagreement of today.
(6) `nearest_agrees` on the rounded distance: point 4. (7) The count the driver notes at the end of every
step, production segments included (`monitor.crystal_ring_site`), is read by neither program, so no number
can differ. Reading it would add a like-for-like check every 5 ns of the production; it would be a new key
and a new version of the specification, to be decided before a run exists.

## 7. What the comparison does not show

- Both programs follow the same specification. A definition that is wrong in the text is wrong in both.
- No run of the campaign exists. The longest trajectories of the real system compared are the three rehearsal
  runs of 5 ns; a run of 100 ns (10,000 and 1,000 frames) was not compared, and no run with an extension step
  made by the driver on the real system exists (in the rehearsal the driver went on after `eq03` in all
  three; extension steps were compared on the driver's stand-in system and on synthetic runs).
- In the 37 runs the two programs classed no frame differently at a threshold (no threshold item). The
  threshold rules of section 14 were exercised by one built folder (a cut-off on a planted water) and by
  planted changes of the comparison's own test.
- Exit status 3 was reached only through planted outputs, not through a run.
- The figures of the first program are not checked by anything here.

## 8. Files

In the repository: `md/analysis/recheck_md.py` (sha256 `5b4c0595…`) and this file. Nothing was committed by
this work. Two commits made by other sessions while it ran hold earlier states of the program: `e27c380`
(11:40:44 AM IST) has it exactly as frozen before the first comparison (sha256 `5ec4c56e…`), `8f7a3a6`
(12:16:29 PM IST) has it one correction before the final one (`bb485115…`, without point 7 of section 4).

Throwaway, under
`/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/mdgr2/recheck/`:

| File | What |
|---|---|
| `run_all.sh`, `run_all.log` | everything of today on the Pi again, one thread, 6 minutes |
| `recheck_today.json`, `recheck_today.md` | every number of sections 2 to 4 (`collect.py`) |
| `frozen_before_compare/` | the second program and its output as they were at 11:39:47 AM IST, with checksums |
| `run_second.sh`, `run_second_real.sh`, `run_primary.sh`, `run_compare.sh` and their logs | the two programs on the 34 runs of the Pi; `run_compare_1.log` is the first comparison |
| `synth/test/` | the synthetic runs (seed of the generator 20261005), with `analysis/` of the first and `analysis_recheck/` of the second program in each |
| `real/` | both programs' output for the 11 test runs of the driver (the run folders themselves were only read) |
| `vm_rehearsal.sh`, `vm_rehearsal_collect.py`, `vm_rehearsal/` | the check on the rented machine: log, the fetched reports and `summary.json` files, `vm_rehearsal.json` |
| `planted_check.sh`, `planted_check/` | the generator's planted values against the second program |
| `edge/` | the built folders (`make_edge_cases.py`, `run_edge_cases.py`, `edge_results.json`), the system without water (`nowater.sh`) and the test of the comparison (`test_compare.py`, `compare_selftest.json`) |
| `unit/` | `unit_geometry.py`, `unit_v13.py` and their results |
| `gate_checks.sh`, `gate_demo/`, `gate_today/` | the one command with `--run-primary` on a campaign-shaped copy of six synthetic runs and on a copy of `/root/cqc-sfa-md/runs` as it is today |
| `timing.py`, `timing.json` | time per frame |
| `v12/` | the second program and this file as they were for version 1.2 |

Software: Python 3.12.14, NumPy 2.5.3 (`/root/cqc-sfa-tools/md`); the first program also uses MDAnalysis
2.10.0. No random numbers in either program; the unit checks use the seed 20261006. Each `summary.json` of
the second program records the sha256 of the program and of every file it read.

## 9. Proposed dated line for `md/PLAN.md`

"6 Oct 2026, 12:41 PM IST: stage S6, before any validation or production run. The second program
(`md/analysis/recheck_md.py`, written from the text of specification 1.3 alone and frozen before the first
program's output was seen) and the first agree in 1,547,845 of 1,547,851 values on 23 synthetic runs, 11 test
runs of the driver and the three rehearsal runs of the rented machine; the six that differ are times of a
test run with frames every 0.1 ps written with three decimals. The driver's hydration record is reproduced by
both programs in all 16 cases that can be compared, the three rehearsal runs among them. The gate of S6 is
`recheck_md.py --compare-campaign /root/cqc-sfa-md/runs --run-primary` ending with exit status 0; status 3
(a verdict that differs through a frame at a threshold) is read by a person. Open sentences of the
specification and four points of the first program: `md/analysis/RECHECK.md`, sections 5 and 6. No verdict,
gate, rule or tolerance changed."
