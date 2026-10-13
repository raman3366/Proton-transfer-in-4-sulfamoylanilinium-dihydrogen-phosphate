# Analysis of the MD campaign: specification

**Version 1.3**, fixed on 6 Oct 2026, 11:26 AM IST, **before any validation or production run exists**
(no folder of a run of the campaign is on this machine: `/root/cqc-sfa-md/runs/` held only 3 folders of tests of the rented machine, `bench`, `preflight_a` and `preflight_b`, at the start and at the end of the tests of this version, 6 Oct 2026, 9:56 AM and 11:22 AM IST). Version 1.0 was fixed on 5 Oct 2026, about 9:30 PM IST, version 1.1 the same evening, 10:55 PM IST,
version 1.2 on 6 Oct 2026, 7:18 AM IST.
This file makes sections 4 and 5 of `md/PLAN.md` operational: every quantity, threshold, formula, file
name and JSON key of the analysis (stage S6) is given here, so that the results cannot shape the definitions and a
second person can implement the analysis from this text alone.

- Primary implementation: `md/analysis/analyse_md.py` (MDAnalysis, NumPy, Matplotlib).
- Test with planted values: `md/analysis/make_test_trajectory.py` (section 16).
- The second, independent script of stage S6 is written from this text, not from the code of the first. What it
  has to know of the files that the driver writes is in appendix A, so that it need not open `md/run_md.py` either.
- A change is a dated line in section 17 and a new version number; the outputs carry the version.

## Changes against version 1.2

Made on 6 Oct 2026, 11:26 AM IST, still before any validation or production run: no trajectory of the campaign
exists, so no result could shape a definition. The folder `/root/cqc-sfa-md/runs/` did not exist on this machine when version 1.2 was fixed, nor when the work on this version began. It exists since the watchdog copies the results of the rented machine to this one (first copy with files at 9:01 AM IST, by its log). At the start and at the end of the tests of this version (9:56 AM and 11:22 AM IST) it held `bench`, `preflight_a` and `preflight_b`: the speed test (`bench`: a minimisation and a timing, no trajectory) and the pre-flight tests of that machine (run id `preflight`: the crystal start of the restrained-zinc model, which the plan has dropped, with the equilibration shortened to 5 % and marked as a test by the driver; their files are dated 9:06 AM IST and later, after the plan's line of 8:56 AM IST that this version carries out). None is a validation or a production run: no folder is named after a start and a replica, and no manifest has such a run id (`runs_state.py`; its output is in `run_all.log`). Of these folders only the fields run id, test flags, `equil_scale`, system folder, status and host of each `manifest.json` were read, and the dates of their files; no trajectory, no log and no number of them was opened. Beside `runs/`, the folder `/root/cqc-sfa-md/runs_test/` holds pipeline tests, which are no data of the campaign. While this version was tested, a rehearsal was run on the rented machine (`md/gpu/tasks_rehearsal_2e.json`: three runs of 5 ns, run ids TX0, TR0 and TN0, with a bonded parameter set that gate G-R has not validated, written to `runs_test/bonded_2e/` and marked there as no result): by the queue of that machine it ran from 10:31 to 11:14 AM IST (`vm_check.txt`, a read-only look at 11:18 AM IST). The two scripts of this version have had the checksums of section 16 since 9:40 AM IST (`scripts_since.txt`), before the rehearsal began, so the definitions were fixed by then; the wording of this text and the staged run of the broken copies (section 16.3) were finished afterwards. No file of the rehearsal was opened for this version. At 11:18 AM IST the machine held the driver and the protocol named below, by their checksums.

The changes carry out what the independent review of gate G-R and of version 1.2 found
(`paper/audit/research/md_gr_review_2026-10-06.json`, the finding "Specification 1.2 against the driver (water at
the ring positions)") and what the plan then fixed (`md/PLAN.md`, "Changes to this plan", 6 Oct 2026, 8:01 AM
point 7 and 8:56 AM point 6). **Only the water count of section 9D and what follows from it change: no verdict, no
gate, no rule and no tolerance of version 1.2 is changed.**

| # | Where | Version 1.2 | Version 1.3 | Why |
|---|---|---|---|---|
| 1 | section 9D; sections 13, 14 | one count: the water oxygens within 3.5 Å of the crystal positions of C1 to C6 and N1 | **the count that the driver acts on**: the water oxygens within 2.5 Å of the crystal positions of the six ring carbons C1 to C6 (`n_water_ring_place`) and the distance of the nearest one (`d_water_ring_place_min_A`), in every full frame of the equilibration and of the production. Site atoms and cut-off are read from the manifest of the run (`hydration_check.settings`) where it gives them, so that the analysis counts what the driver of that run counted; otherwise they are the six ring carbons and 2.5 Å. A run made with other settings is reported with its own and gets a warning | plan, 8:01 AM, point 7: the aniline nitrogen points into the solvent, so a water near its place says nothing of the place of the ring; plan, 8:56 AM, point 6: a water oxygen within 3.5 Å of a ring-carbon position is found in the crystal structure with the ring still in place; review, finding on version 1.2 |
| 2 | section 9D; section 13 | `n_water_ring_site_3p5A`, `d_water_ring_site_min_A`; in the table row `n_water_ring_site_3p5A_mean`, `d_water_ring_site_min_min_A` | the count of version 1.2 stays as a second, reported quantity with the same definition and the same values, under the names `n_water_ring_wide_3p5A`, `d_water_ring_wide_min_A`, `n_water_ring_wide_3p5A_mean`, `d_water_ring_wide_min_min_A`. The four names with `ring_site` no longer exist | "site" is the driver's word for its own positions, which are now six; one word for two sets of positions would mislead. No file of the campaign was ever written with the old names |
| 3 | section 9E | the driver's count at the end of a step compared with the count of version 1.2 | compared with `n_water_ring_place` and `d_water_ring_place_min_A` of the full frame at the end of that step; the rule for "at the end of the step" written with the step table (S + m Δt = E) | with the comparison of version 1.2 every run of the driver as it stands would have carried a false disagreement (review) |
| 4 | sections 1, 2, 13, 14, 15 | | the texts that name the manifest keys read, the columns, keys and figures, the agreement of two programs and the limits follow changes 1 to 3. New: two columns at the end of `timeseries_full.csv` and of `timeseries_equil_full.csv`, two entries at the end of the table row, the key `water_ring_place` of `summary.json`, and four numbers instead of two in `last_full_frame` of the equilibration steps | |
| 5 | new appendix A | the reader had to open `md/run_md.py` for the layout of its files | "Reading the driver's files": the file names, the DCD layout as far as a reader needs it, the manifest keys used, the time axis, and the hydration check in brief | so that the second program is written from this text alone |
| 6 | sections 1, 2, 9D, 9E | | wordings made exact where a second reader could have read them in two ways (list below); none changes what the primary implementation computes | |
| 7 | section 16 | 6232 planted values, 74 broken copies | the tests of version 1.2 kept, with the new columns planted in every run; two new synthetic runs; run ext/R1 with the settings of its hydration check in its manifest; the broken copy "the six ring carbons only" reversed and nine new ones; two short runs of the driver as it stands, with extension steps, read and compared count by count | |

Not changed: every other quantity, statistic, verdict, gate, rule and tolerance of version 1.2, and the names and
the order of its columns and keys but for the four renamed in change 2. On the 21 synthetic runs of
version 1.2 the two versions write the same values in every column and under every key that version 1.2 had
(section 16.4); one run of them, ext/R1, has one warning more, because its manifest now names the settings of the
driver as it stood before 8:01 AM.

**The driver that this version follows**: `md/run_md.py` with sha256 `6d797f080ff7567923469dabd86af71ca34183f5c482c14d55bf8038eb9d1822`
(6 Oct 2026, 9:53 AM IST. It is the file that version 1.2 named, sha256 `b0f3e2b23b3eedf5…`, with one line added while this version was written: the header of a new DCD file is written to the disk at once, so that the file of a step that has just begun is a DCD file without a frame and not an empty file. No name, no key and no layout changed (`driver_diff.txt` of the tests)) together with `md/protocol.json` with sha256
`decda42f5340d06f5f695b4c56b59a7bb796e5872b3a8af1c1a5b426f9c09d9a` (6 Oct 2026, 8:56 AM IST: site atoms C1 to C6, cut-off 0.25 nm). For this version its code
was read again where it writes what the analysis reads: the record of a step (`rec_<step>.json`, copied into
`equilibration.stages[]` and `production.segments[]` of the manifest), the hydration check with its records and
its summary in the manifest, the names of the files of the extension steps, the DCD writer and the PDB writer.
The keys that version 1.2 read are the keys this driver writes; what was added is the pair
`hydration_check.settings.site_atoms` and `hydration_check.settings.water_cutoff_A` (section 9D). Appendix A is
the result of that reading, and section 16.4 has the check on two runs made by this driver.

**Wordings made exact** (change 6):

- Section 1: the keys of the driver's manifest are named for the driver as it stands, with the new pair; no other
  key of a manifest is read; what happens with a manifest that cannot be read.
- Section 2: the minimum-image distance of the full trajectory is given by a formula, with the lattice vectors of
  appendix A.
- Section 9D: the site positions carried into the frame are given by a formula; the conditions under which the
  settings of a manifest are used, and that both come from the same place, are spelled out.
- Section 9E: for which steps `driver_hydration_check` is written; what `compared` means; which frame is compared.
- Section 9E, guard "header": what the header of the driver holds, and in which test runs it can differ from the
  step table by half a time step (the rule and its tolerance are unchanged).
- Section 14: frames that two correct programs may class differently, for the new count and for the comparison
  with the driver.

One wording of the task was made operational: "read the site atoms and the cut-off from the manifest's
hydration_check settings when present" is read as: both from the manifest if both are usable there, else both
from the defaults; never one of each.

## Changes of version 1.2 against version 1.1 (kept for the record)

Made on 6 Oct 2026, 7:18 AM IST, still before any run: no trajectory of the campaign exists, so no result could
shape a definition. The only run folder on this machine is a pipeline test of another agent
(`/root/cqc-sfa-md/runs_test/restrained/X1`, shortened equilibration, restrained-zinc model), which held no
production frame when the tests began. Its manifest was looked at for the layout of the stage records of a shortened
equilibration (a file format); no coordinate and no number of it was read. While this version was written the
driver was changed by its maintainer (`md/run_md.py` of 6 Oct 2026, 5:28 AM IST, sha256 `b0f3e2b23b3eedf5…`: the
hydration check of the plan, with extension steps `eq03x1` to `eq03x4` of the held part); section 9E reads the files
and the manifest of that driver. The changes carry out the reported quantities that the
independent review of the set-up (`paper/audit/research/md_setup_review_2026-10-06.json`, finding 1) and the plan
(`md/PLAN.md`, "Changes to this plan", 6 Oct 2026, 1:03 AM point 3 and 3:02 AM points 1, 4 and 8) ask of the analysis.
**Every new quantity is reported only: no verdict, no gate, no rule and no tolerance of version 1.1 is changed.**

| # | Where | Version 1.1 | Version 1.2 | Why |
|---|---|---|---|---|
| 1 | new section 9A; sections 13, 14 | | the configuration at the sulfonamide nitrogen N2 of every frame: the signed volume `n2_volume` of the unit vectors from N2 to Zn, to HN21 and to S (in this order) and its sign; frames of each sign over the run, over the last 20 ns and per orientation state; the sign changes that last at least 1 ns, with their times | review, finding 1: the crystal start and the rotated start differ in this configuration (+0.72 for X, +0.79 for N, −0.71 for R), so a conversion needs an inversion at N2 or the loss of the hydrogen bond to Thr199 |
| 2 | new section 9B; section 13 | the fraction of frames with N2–H···OG1(Thr199) over the run | also this fraction before the first passage and from it on, for every run that has a first passage | the same finding; plan, 3:02 AM, point 1 |
| 3 | new section 9C; sections 11, 13 | Glu106 seen only as "a protein oxygen" of the full trajectory, every 100 ps | the distances from the zinc to OE1 and to OE2 of Glu106 by name, every 10 ps, with mean, minimum and the fraction of frames within 2.6 Å; printed beside gate G-V; no criterion | plan, 3:02 AM, point 4: Glu106 OE1 moved onto a free zinc in the tests of the dropped fallback |
| 4 | new section 9D; section 13 | | the number of water oxygens within 3.5 Å of the positions that C1 to C6 and N1 have in the crystal pose, and the distance of the nearest one, in every full frame | plan, 1:03 AM point 3 and 3:02 AM point 8: the empty pocket of start R |
| 5 | new section 9E; sections 1, 13 | equilibration files not read | `solute_eqNN.dcd` and `full_eqNN.dcd`, and those of the extension steps of the held part (`eqNNxM`), if present, are read into a time series of their own with negative times; the driver's own count of water at the ring positions is set beside the count of change 4; whether the first production frame differs from the start as built (state, sign at N2) | plan, 3:02 AM: with the ligand free from stage 4 on, a change of pose or an inversion can happen before time zero; plan, 1:03 AM, point 3: the held part may be extended |
| 6 | sections 1, 2, 10 to 15 | | the texts that name the files read, the statistics, the output files, columns and keys, the agreement of two programs and the limits of the analysis follow changes 1 to 5; new columns stand after the columns of version 1.1, which keep their names and their order | |
| 7 | section 16 | 2917 planted values, 40 broken copies | the tests of version 1.1 kept; five new synthetic runs, four of them with equilibration files, one with extension steps; planted values and broken copies for every new quantity; the built systems as one-frame runs; the driver's own test run with its equilibration files | |

Not changed: every quantity, statistic, verdict, gate, rule and tolerance of version 1.1 (sections 2 to 12), and the
names and the order of its columns and keys. On the sixteen synthetic runs of version 1.1 the two versions write the
same values in every column and under every key that version 1.1 had (section 16.4). Two keys were added inside
blocks of version 1.1: `orientation.first_passage_frame` (the frame that belongs to `first_passage_ns`) and
`gate_GV.reported_zn_glu106` (numbers beside the gate; they are no part of `gate_GV.passed`).

One wording of the task was made operational: "the fraction of frames below 2.6 Å" for Glu106 is counted with the
limit included (distance ≤ 2.6 Å), as the fifth ligand of section 6 is, so that the two counts cannot differ by a
frame that sits on the limit.

## Changes of version 1.1 against version 1.0 (kept for the record)

Made on 5 Oct 2026, 10:55 PM IST, still before any run: no trajectory of the campaign exists, so no result
could shape a definition. The changes carry out the four decisions that the owner of the plan entered in
`md/PLAN.md` on 5 Oct 2026, 9:41 PM IST ("Changes to this plan"), after version 1.0 had raised them, and one
change of the conventions.

| # | Where | Version 1.0 | Version 1.1 | Why |
|---|---|---|---|---|
| 1 | section 11, GV3 c | each mean N–Zn–N angle in 100–120° | each mean N–Zn–N angle within 10° of its value in 6RL9 (six references, printed with the atom names) | the window was wrong: the 1.0 Å structure itself has 98.4° for His96 NE2–Zn–His119 ND1 (plan, point 1) |
| 2 | sections 6 and 11, GV3 d | no frame with an oxygen within 2.6 Å of the zinc | an oxygen within 2.6 Å in at most 5 % of the frames, counted as two fractions, one per trajectory, each at most 5 % | plan, point 2; the sulfonyl oxygens are written every 10 ps, water and protein oxygens every 100 ps only |
| 3 | section 11, GV4 | τ_fold and C4–S–N2 reported, no pass or fail | mean C4–S–N2 within 3° of 108.8° and mean τ_fold within 20° of 132.1°; both enter the verdict | plan, point 3 |
| 4 | sections 1, 5, 11, 12 | gate G-V computed with `--planned-ns 5`; "applies" to every run of start X | the validation run is run `X0` of its zinc model: recognised as start X with replica 0 (folder name `X0`, or the manifest), planned length 5 ns without an option, the only run with a verdict of gate G-V, never in the campaign table, not counted by rule 1; the production runs get the same numbers without a verdict | plan, point 4 |
| 5 | sections 1, 4, 13 | a turn of the Asn67 amide could not be seen | optional key `asn67_ND2` of `selections.json`: the distance N1 to ND2 and a flag per frame where ND2 is the nearer atom; the states still use OD1 only | conventions of the system builder |
| 6 | section 13 | | new keys and columns: `run.role`, `gate_GV` (new layout), `asn67_ND2`, two optional columns of `timeseries_solute.csv`, three columns of the campaign table | follow from 1 to 5 |
| 7 | section 16 | 856 planted values, 14 broken copies | the tests of version 1.0 kept; a validation run that passes the gate, runs in which each new criterion fails alone, broken copies for the new rules | the crystal geometry now meets the criteria |

Not changed: the definition of every quantity of a frame (sections 2 to 9), the statistics (section 10), the
orientation states, the verdict of a run and the rules of the campaign; on the test runs of version 1.0 the two
versions write the same time series and statistics (section 16.4). Section 1 now spells out the order in which
start and replica are found, because the role of a run depends on it; that order is the one of version 1.0.
The two remarks of version 1.0 on the plan (section 15) are settled by changes 1 and 2 and are removed.

## 0. Notation

- Coordinates and distances in Å (as stored in DCD and PDB files), angles and torsions in degrees, times in ns,
  frame intervals in ps. Force-field units do not occur.
- Atom names are those of the built system: the ligand SAN 304 has the heavy atoms C1–C6, N1 (aniline N), N2
  (sulfonamide N, on the zinc), O1, O2, S, and the sulfonamide hydrogen HN21 (`lig_HN2`); C1 carries N1, C4 carries
  S, C2 and C6 are next to C1, C3 and C5 next to C4. The paper's S1 is S here.
- Residue numbers are those of the deposited structure 6RL9 (His94, His96, His119, Thr199, Glu106, Phe131, Asn67),
  mapped to atom indices by `selections.json`.
- |a − b| is the distance of two atoms. ∠(a, b, c) is the angle at b, computed as atan2(|u × v|, u · v) with
  u = a − b, v = c − b. τ(a, b, c, d) is the torsion in (−180°, 180°] with the IUPAC sign: with b₁ = c − b
  normalised, v = (a − b) − ((a − b) · b₁) b₁ and w = (d − c) − ((d − c) · b₁) b₁, τ = atan2((b₁ × v) · w, v · w).
- "At most", "at least" and "within" include the limit (≤, ≥); "below" and "above" do not (<, >). "x within t
  of r" is |x − r| ≤ t.
- A **mean** without another word is the arithmetic mean over all frames read of the trajectory concerned
  (section 10: `stats.<column>.mean`).
- N is the number of frames of a trajectory, n = 0 … N − 1 the frame index through the whole run. N, n and every
  statistic belong to the production (the files `*_segNNNN.dcd`); the equilibration frames of section 9E are counted
  apart and enter nothing else.
- u × v = (u_y v_z − u_z v_y, u_z v_x − u_x v_z, u_x v_y − u_y v_x); u · v = u_x v_x + u_y v_y + u_z v_z.

## 1. Input and what is assumed of it

**Run folder** `/root/cqc-sfa-md/runs/<model>/<start><replica>/` (model `bonded` or `restrained`; start X, N or R;
replica 1, 2, 3 for the nine production runs of a model; `X0` for the validation run of the model, see "Role of a
run" below):

| File | Use |
|---|---|
| `solute_segNNNN.dcd` | the "solute" atoms (protein, zinc, ligand), one frame every 10 ps; segments in ascending NNNN |
| `full_segNNNN.dcd` | all atoms, one frame every 100 ps |
| `solute.pdb`, `full.pdb` | atom names and residue names, in the atom order of the trajectory of the same kind |
| `manifest.json` | optional; see below |
| `solute_eqNN.dcd`, `full_eqNN.dcd`, `solute_eqNNxM.dcd`, `full_eqNNxM.dcd` | optional: the frames of equilibration stage NN and of extension step M of that stage, same atoms and intervals as the production files of the same kind; read into a series of their own (section 9E) |

Other files of the folder (logs, saved states, checkpoints, the records `rec_<step>.json` of the steps) are not
read. Appendix A describes the files of the driver: their names, the layout of the DCD files, the manifest keys
used and the time axis.

**System folder** (`selections.json`, `reference.pdb`): taken from the command line (`--system`); else
`system_dir` of the manifest (a relative path counts from the run folder) if that folder is on this machine; else
`<run folder>/../../../systems/<model>/<start>/`. The last case is the normal one after the runs have been
pulled from the rented machine.

- `selections.json` holds zero-based atom indices of the topology. Keys read: `zn`, `lig_heavy` (eleven names),
  `lig_HN2`, `lig_all`, `his94_NE2`, `his96_NE2`, `his119_ND1`, `thr199_OG1`, `thr199_HG1`, `thr199_N`,
  `thr199_H`, `glu106_OE1`, `glu106_OE2`, `phe131_CZ`, `asn67_OD1`, `protein_CA`, `protein_heavy`, `water_O`,
  `solute`, `residue_numbers_deposited_for_CA`. One key is **optional**: `asn67_ND2` (section 4); a key that is
  missing or null is absent, and the analysis runs without it.
- **Index rule.** In the full trajectory an atom is found at its index. In the solute trajectory it is found at
  the position of its index in the list `solute`.
- `reference.pdb` (crystal coordinates, 6RL9 frame): the atoms named `CA`, in file order, are the reference of
  the atoms `protein_CA`, in list order; the ligand reference atom of a name is the one atom of the file with that
  name.

**Who the run is.** Start, replica, model, planned length and the two frame intervals are taken, in this order,
from the command line, from `manifest.json`, from the folder name (`<start><replica>`, parent folder = model), from
the defaults (100 ns, 10 ps, 100 ps). Two forms of manifest are read: the plain keys `start`, `replica`, `model`,
`planned_ns`, `solute_interval_ps`, `full_interval_ps`, `system_dir`; and the keys written by the driver
`md/run_md.py` (appendix A gives each with its meaning): `run_id`, `system_dir`, `build_json.start`,
`build_json.model`, `production.target_ns`, `protocol.production.length_ns`,
`protocol.output.solute_frame_interval_ps`, `protocol.output.full_frame_interval_ps`,
`production.segments[].{tag, time_start_ps, frames.solute, frames.full}`,
`inputs.{selections.json, reference.pdb}.sha256`, for the equilibration (section 9E)
`equilibration.stages[].{tag, time_start_ps, time_end_ps, hydration_check.{n_water_O, nearest_water_O.distance_A,
applies, decision}}`, and for the water count of section 9D `hydration_check.settings.{site_atoms,
water_cutoff_A}`. No other key of a manifest is read. A `manifest.json` that cannot be read as JSON gives a warning
and is then treated as a manifest without any key.
`summary.json` records where each value came from (`run.sources`). The first source that gives a value is used:

| Item | Sources, in this order |
|---|---|
| start | `--start`; manifest `start`; manifest `build_json.start` if it is X, N or R; the letter of manifest `run_id`; the letter of the folder name |
| replica | manifest `replica`; the number of manifest `run_id`; the number of the folder name; else none |
| model | `--model`; manifest `model`; manifest `build_json.model`; the name of the parent folder |
| planned length | production run: `--planned-ns`; manifest `planned_ns`; manifest `production.target_ns`; manifest `protocol.production.length_ns`; 100 ns. Validation run: `--planned-ns`; else **5 ns** (the manifest is not asked) |
| interval of the solute frames | `--solute-interval-ps`; manifest `solute_interval_ps`; manifest `protocol.output.solute_frame_interval_ps`; 10 ps |
| interval of the full frames | `--full-interval-ps`; manifest `full_interval_ps`; manifest `protocol.output.full_frame_interval_ps`; 100 ps |

A `run_id` or a folder name gives a letter and a number only if it is one of the letters X, N, R followed by digits
and nothing else (`X1`, `R3`, `X0`); the number is read as an integer. For the role the replica is compared as an
integer; a replica that is none, or a manifest value that is not an integer, never gives the role `validation`.

**Role of a run** (`run.role`). With start and replica found as above:

- `validation`: start X and replica 0. This is the 5 ns validation run of stage S4 of the plan, "run X0" of the
  zinc model of its parent folder (`/root/cqc-sfa-md/runs/<model>/X0`). There is one per zinc model.
- `production`: every other case (replica 1, 2, …, no replica, or replica 0 with start N or R, which the plan
  does not have).

What follows from the role `validation`, without any option on the command line:

1. the planned length is 5 ns (500 solute frames, 50 full frames), whatever the manifest says;
2. gate G-V is evaluated with a verdict (section 11); for every other run the same numbers are written without one;
3. rule 1 of the plan is not applied to the run (section 5), and the run is in no campaign table (section 12).

Three warnings belong to the role. (i) Validation run: `--planned-ns` is not given and the manifest names another
length (the value of `planned_ns`, or without it of `production.target_ns`, differs from 5 ns by more than
10⁻⁹ ns). (ii) Validation run: the run holds more solute frames than planned; all frames read are used, as for
every run, and nothing is cut off. (iii) Any run whose folder name has the form above: the letter and the number
of the folder name alone would give the other role (validation for X with 0, production for every other pair)
than the sources gave; the role found from the sources stands.

**Time.** The segments of a kind, in ascending NNNN, are one continuous trajectory. Frame n has the time
t = (n + 1) Δt, with Δt = 10 ps (solute) or 100 ps (full): the first frame is written one interval after the
start of the production, and time zero is the end of the equilibration. Full frame j therefore belongs to solute
frame 10 (j + 1) − 1. The length seen is N Δt. A run is **complete** if N ≥ planned length / Δt; the analysis
works on an incomplete run and says how many ns it saw. Only complete frames are read; a frame cut off at the end
of a file and a file without a frame are left out.

**What the driver must keep** (the analysis cannot repair it): no frame is written twice and none is missing
between segments; the atom order of each PDB file is that of its trajectory; every frame carries its unit cell;
the protein is written whole (not broken across the cell).

**Guards.** They are computed for every run and written to `summary.json` (`run.*`, `warnings`).

| Guard | Rule | If it fails |
|---|---|---|
| atom counts | atoms per frame = atoms of the PDB file; `solute.pdb` has as many atoms as the list `solute` | stop |
| system files | where the manifest records the checksums of `selections.json` and `reference.pdb` (`inputs.<file>.sha256`), the files used have them | warning |
| atom names | the names at the mapped positions are ZN, NE2, NE2, ND1, OG1, HG1, N, H, OE1, OE2, CZ, OD1, the eleven ligand names, a name starting with H for `lig_HN2`, CA for every `protein_CA`, a name starting with O for every `water_O`, ND2 for `asn67_ND2` if that key is present; the residue names of the CA atoms equal those of `reference.pdb` | stop (`--skip-name-check` turns it into a warning) |
| interval | the interval in each DCD header equals Δt within 0.1 % | warning |
| contiguous headers | where the header counts the steps of the whole run, each segment starts one interval after the last frame of the one before it; a header that starts at one interval in every file (OpenMM's own reporter) tells nothing and is recorded as "none" | warning |
| manifest | where the manifest lists the segments, its frame count of a segment equals the frames read and its start time equals the frames before it × Δt | warning |
| repeated frame | the first frame of a segment equals no frame of the segment before it (CA coordinates, largest difference below 0.001 Å = repeated) | warning |
| two trajectories | Zn–N2 of full frame j and of solute frame 10 (j + 1) − 1 differ by at most 0.01 Å | warning |
| whole protein | CA atoms of neighbouring residues (numbers i and i + 1, and 125–127, which are joined) are at most 4.5 Å apart in every frame | warning |

A run with a warning of these kinds has the wrong atoms, a wrong time axis or wrong coordinates; its numbers are
not used until the cause is found.

## 2. Frames, periodic images and the fit

Everything of sections 3 to 9C that is not marked "full trajectory" comes from the solute trajectory, every frame.

1. **Images** (solute trajectory, frames with a unit cell). The zinc is moved by the lattice vector that makes
   |Zn − NE2(His94)| smallest. Then all atoms of the ligand (`lig_all`) are moved together by the lattice vector
   that makes |N2 − Zn| smallest. The protein is not touched. Without a unit cell nothing is moved, and this is
   recorded.
2. **Fit.** The CA atoms of the frame (`protein_CA`, 258 atoms) are superposed on those of `reference.pdb` by least
   squares with equal weights (no masses), rotation without reflection (Kabsch). The same rotation and translation
   are applied to all atoms of the frame. Nothing else is fitted: the ligand is never fitted on itself.

Distances, angles and torsions do not depend on the fit; the ligand RMSD, the CA RMSD and the RMSF do.

**Full trajectory.** No fit. Every distance is the minimum-image distance in the unit cell of the frame: with the
lattice vectors a, b, c of the frame (appendix A.2), the distance of two points p and q is the smallest value of
|p − q + i a + j b + k c| over all whole numbers i, j, k. The cell of the campaign is a truncated octahedron, in
which rounding the three fractional coordinates one after the other does not always give the smallest value; an
implementation searches the neighbouring cells as well, or shows that it need not. The lattice vector of step 1
above is found in the same way.
The two water counts of section 9D need the fit of step 2 on the CA atoms of the full frame; it is described there
and touches no other quantity.

## 3. Ligand RMSD to the crystal pose

After the fit, with x_a the position of ligand heavy atom a in the frame and r_a its position in `reference.pdb`,
for a labelling σ (a permutation of the eleven names):

RMSD_σ = sqrt( (1/11) Σ_a |x_a − r_σ(a)|² )

The heavy-atom graph of the molecule has four equivalent labellings:

| Number | Name | σ |
|---|---|---|
| 0 | `identity` | every atom with its own name |
| 1 | `ring_flip` | C2 ↔ C6 and C3 ↔ C5 |
| 2 | `oxygen_swap` | O1 ↔ O2 |
| 3 | `ring_flip+oxygen_swap` | both |

`lig_rmsd_A` is the smallest of the four values, `lig_rmsd_labelling` the number of the labelling that gave it
(the lowest number if two are equal), `lig_rmsd_identity_A` the value of labelling 0.

**Check against `md/system/starts.csv`.** With this definition the three starts give 0.0000, 0.6104 and 3.5732 Å;
`starts.csv` has 0.0, 0.61 and 3.573 Å. The definition is the one the system builder used (symmetry-corrected,
no superposition of the ligand), and the fit changes nothing there because the protein of the starts is the
crystal structure. Table: `analyse_md.py --check-starts` (section 16.1).

## 4. Orientation of the ring (plan, section 5, rule 2)

- d1 = |N1 − CZ(Phe131)| (`d1_N1_F131CZ_A`); d2 = |N1 − OD1(Asn67)| (`d2_N1_N67OD1_A`).
- Midpoints, fixed: **m1 = 6.7 Å, m2 = 5.2 Å**, halfway between the crystal start (d1 3.897, d2 7.744 Å) and the
  rotated start (9.585, 2.597 Å) of `starts.csv` (exact midpoints 6.741 and 5.171 Å). The near-native start
  (3.356, 7.725 Å) is on the crystal side of both.

**State of a frame** (`state`):

| State | Condition |
|---|---|
| C, crystal-like | d1 < m1 and d2 > m2 |
| R, rotated-like | d1 > m1 and d2 < m2 |
| I, intermediate | everything else |

**State of the start:** C for X and N, R for R. The **opposite state** is the other of C and R.

**Verdict of a run.** Let W = min(N, 2000) frames, the last 20 ns, and c_S the number of frames of state S among
them (counted from `state`, intermediate frames count for neither):

- **CONVERTED** if 10 c_opposite ≥ 9 W;
- else **STAYED** if 10 c_start ≥ 9 W;
- else **MIXED**.

The comparison is made with integers. The verdict of a run that is not complete is **provisional**; if the run is
shorter than 20 ns the window is the whole run (`window_complete` false).

**Visits, first passage, transitions.**

1. Assigned state a_n (`state_assigned`): the state of frame n if it is C or R; otherwise a_(n−1); before the first
   frame the state of the start. (An intermediate frame belongs to the state last visited.)
2. A **visit** is a maximal stretch of frames with the same assigned state. It is **long** if it has at least
   D = 100 frames (1 ns).
3. **First passage** (`first_passage_ns`): the time of the first frame of the first long visit in the opposite
   state; none if there is no such visit.
4. **Transitions** (`n_transitions`): the number of changes of state in the sequence (state of the start, states of
   the long visits in order of time). Short visits are ignored.

For each long visit `fraction_strict` is the fraction of its frames whose `state` equals its assigned state.
Also given: the numbers and fractions of C, R and I frames over the whole run, the last half and the window.

**Asn67 ND2 (optional).** Only if `selections.json` has the key `asn67_ND2`, two more quantities are written for
every solute frame:

- `d_N1_N67ND2_A` = |N1 − ND2(Asn67)|;
- `asn67_ND2_nearer_N1` = 1 if `d_N1_N67ND2_A` < `d2_N1_N67OD1_A`, else 0 (equal distances give 0).

They are reported and nothing else: the state of a frame, the verdict of a run and every rule use d2 to OD1 as
defined above, also where ND2 is the nearer atom. Without the key the two quantities, their statistics and their
columns are not written (`asn67_ND2.available` false), and everything else is unchanged.

How the flag reads depends on the pose. Distances of N1 from OD1 and from ND2 in the starts as built
(`md/system/starts.csv`): rotated start 2.597 and 3.310 Å; crystal start 7.744 and 6.551 Å;
near-native start 7.725 and 6.557 Å.
In the rotated pose OD1 is the nearer atom, so in a rotated-like frame a flag of 1 means that the amide has turned
(OD1 and ND2 exchanged). In the crystal and the near-native pose ND2 is the nearer atom with the amide as the
crystal has it, so in a crystal-like frame the flag is 1 without any turn, and a flag of 0 means the turn. The flag
is therefore counted per state of the frame (`asn67_ND2.by_state`, by the column `state`); no single "turned" flag
is defined.

## 5. Validation quantities (gate G-V; plan, section 5, rule 1)

- **Hydrogen bond N2–H···OG1(Thr199)** (`hb_N2H_T199OG1`): present in a frame if
  |N2 − OG1| ≤ 3.5 Å and ∠(N2, HN21, OG1) ≥ 135° (angle at the hydrogen). Reported with the distance
  `d_N2_T199OG1_A`, the angle `ang_N2_HN2_T199OG1_deg` and the fraction of frames.
- **Rule 1.** "Ligand RMSD above 2 Å for the last half of the run" is: the mean of `lig_rmsd_A` over the last
  half of the run is above 2.0 Å (`rule1_left_crystal`). It is computed for every run and applies to the
  production runs of start X (`rule1_applies_to_this_run`: start X and role `production`); the validation run X0
  is not one of the crystal-start runs that rule 1 counts.

## 6. Zinc site

Solute trajectory:

- the four distances |Zn − N2|, |Zn − NE2(His94)|, |Zn − NE2(His96)|, |Zn − ND1(His119)|
  (`d_Zn_N2_A`, `d_Zn_H94NE2_A`, `d_Zn_H96NE2_A`, `d_Zn_H119ND1_A`);
- the six angles ∠(A, Zn, B) over the pairs of these four nitrogen atoms (`ang_<A>_Zn_<B>_deg` with A, B in the
  order N2, H94NE2, H96NE2, H119ND1);
- ∠(S, N2, Zn) (`ang_S_N2_Zn_deg`);
- |Zn − O1|, |Zn − O2| and the smaller of the two (`d_Zn_Osulfonyl_min_A`).

Full trajectory (every 100 ps, minimum image):

- `d_Zn_Owater_min_A`: the smallest distance from the zinc to an atom of `water_O`, with the index of that atom;
- `d_Zn_Oprotein_min_A`: the smallest distance to a protein oxygen, i.e. an atom of `protein_heavy` whose name
  in `full.pdb` starts with O, with its index;
- `d_Zn_Osulfonyl_min_full_A`: the smaller of |Zn − O1| and |Zn − O2| in that frame;
- `d_Zn_O_any_min_A`: the smallest of the three.

If two atoms are equally near, the first in the order of the list is named.

**Fifth ligand:** an oxygen within 2.6 Å of the zinc (distance ≤ 2.6 Å). Flags of a frame, 1 or 0:

| Flag | Trajectory | Is 1 if |
|---|---|---|
| `fifth_sulfonylO` | solute, every 10 ps | `d_Zn_Osulfonyl_min_A` ≤ 2.6 |
| `fifth_water` | full, every 100 ps | `d_Zn_Owater_min_A` ≤ 2.6 |
| `fifth_proteinO` | full | `d_Zn_Oprotein_min_A` ≤ 2.6 |
| `fifth_any` | full | `d_Zn_O_any_min_A` ≤ 2.6, that is: a sulfonyl, a water or a protein oxygen is within 2.6 Å in that full frame |

**The two fractions of the fifth ligand** (plan: "an oxygen (ligand sulfonyl, water or protein) within 2.6 Å of the
zinc in at most 5 % of the frames"). The sulfonyl oxygens are in both trajectories, water and protein oxygens only
in the full one. The 5 % is therefore counted twice, once in each trajectory, over all frames read of the run:

- **solute fraction** c_s / N_s: N_s is the number of solute frames read, c_s the number of them with
  `fifth_sulfonylO` = 1 (`fractions.fifth_sulfonylO`: `count`, `n`, `fraction`);
- **full fraction** c_f / N_f: N_f is the number of full frames read, c_f the number of them with
  `fifth_any` = 1 (`fractions.fifth_any`).

The criterion "no fifth ligand" (section 11, GV3 d) is met if **both** fractions are at most 5 %, in integers:
20 c_s ≤ N_s **and** N_f ≥ 1 **and** 20 c_f ≤ N_f. The two trajectories are never pooled: a solute frame is not
counted in the full fraction, a full frame is not counted in the solute fraction (the solute frame of the same time
is), and no weight is given to either. A sulfonyl oxygen near the zinc in a full frame counts in the full fraction
through `fifth_any`. Without a full trajectory (N_f = 0) the full fraction does not exist and the criterion is not
met. Both fractions are reported for every run, whether or not a gate depends on them.

## 7. Other contacts

- **Thr199 N–H···O(sulfonyl)** (`hb_T199NH_Osulfonyl`): present if for O1 or for O2 |N(Thr199) − O| ≤ 3.5 Å and
  ∠(N, H, O) ≥ 135°. Both distances, both angles and the smaller distance are written.
- **Thr199 OG1–HG1···OE(Glu106)** (`hb_T199OG1H_E106OE`): present if for OE1 or for OE2 |OG1 − OE| ≤ 3.5 Å and
  ∠(OG1, HG1, OE) ≥ 135°. Both distances, both angles and the smaller distance are written.
- **Water at the aniline N** (full trajectory): `n_water_N1_3p5A`, the number of atoms of `water_O` within 3.5 Å
  of N1 (minimum image).

## 8. Ligand internal coordinates

- `tors_C3_C4_S_N2_deg` = τ(C3, C4, S, N2) and `tors_C5_C4_S_N2_deg` = τ(C5, C4, S, N2);
- `tau_fold_deg` = max(|τ(C3, C4, S, N2)|, |τ(C5, C4, S, N2)|), the τ_fold of the paper (Section 3.8: "the larger
  of |C3–C4–S1–N2| and |C5–C4–S1–N2|"), which does not depend on which ring carbon is called C3; the crystal ligand
  has 132.1°;
- `ang_C4_S_N2_deg` = ∠(C4, S, N2) (crystal 108.8°);
- `d_S_N2_A` = |S − N2| (crystal 1.596 Å).

The two signed torsions exchange their values when the ring turns by 180°; τ_fold does not. Their means are
circular (section 10).

## 9. Protein

- `ca_rmsd_A`: after the fit, sqrt( (1/M) Σ_i |x_i − r_i|² ) over the M CA atoms, against `reference.pdb`.
- **RMSF** per residue: with x_i(n) the fitted position of CA atom i in frame n and ⟨x_i⟩ its mean over all frames
  of the run, RMSF_i = sqrt( (1/N) Σ_n |x_i(n) − ⟨x_i⟩|² ). The frames are those fitted to `reference.pdb`; there
  is no second fit to the mean structure. Residues carry the deposited numbers (3–125, 127–261).

## 9A. Configuration at the sulfonamide nitrogen N2 (reported)

The crystal start and the rotated start are not two rotamers about N–S: with the bond Zn–N2 and the hydrogen bond
N2–H···OG1(Thr199) kept, the sulfur sits on one or the other of the two remaining tetrahedral sites of N2 (plan,
change of 6 Oct 2026, 3:02 AM, point 1). The quantity below tells the two configurations apart in every frame.

**Signed volume.** For a solute frame, after step 1 of section 2 (the fit plays no part), with HN21 the atom
`lig_HN2`:

- e_Zn = (Zn − N2) / |Zn − N2|, e_H = (HN21 − N2) / |HN21 − N2|, e_S = (S − N2) / |S − N2|;
- `n2_volume` = V = e_Zn · (e_H × e_S).

The order is Zn, HN21, S. V is the determinant of the 3 × 3 matrix whose first row is e_Zn, whose second row is
e_H and whose third row is e_S. It has no unit and lies in −1 … 1. It does not change when the frame is rotated or
shifted; it changes its sign when two of the three atoms are exchanged or the frame is mirrored; V = 0 means that
the three directions lie in one plane (a planar nitrogen).

- `n2_sign` = +1 if V > 0, −1 if V < 0, 0 if V = 0. The sign is taken from V as computed in double precision from
  the stored coordinates, not from the rounded number in the CSV file.
- **Sign of the start** (constant): **+1 for the starts X and N, −1 for the start R.** The starts as built have
  V = +0.7240 (X), +0.7855 (N) and −0.7104 (R); the review gives +0.72, +0.79 and −0.71 (section 16.1).

**Counts.** From `n2_sign`, three numbers each time, the frames with +1, with −1 and with 0:

- over all N frames of the run, with the fractions count / N;
- over the window of the verdict, the last W = min(N, 2000) frames (section 4), with the fractions count / W;
- per orientation state: over the frames whose `state` (the state of the frame, section 4; not `state_assigned`)
  is C, over those with R and over those with I, each with the number of such frames n and the fractions count / n
  (none if n = 0).

**Sign changes that last at least 1 ns.**

1. Assigned sign s_n (`n2_sign_assigned`): `n2_sign` of frame n if it is +1 or −1; otherwise s_(n−1); before the
   first frame the sign of the start. (A frame with V = 0 belongs to the sign before it.)
2. A **stretch** is a maximal run of consecutive frames with the same assigned sign. It is **long** if it has at
   least D = 100 frames (1 ns, as the dwell of section 4).
3. Take the sequence (sign of the start, signs of the long stretches in order of time). Every long stretch whose
   sign differs from the element before it in this sequence is a **lasting sign change**. Its frame is the first
   frame f of that stretch, its time (f + 1) Δt, `from_sign` the element before, `to_sign` its own sign.

Short stretches are left out of the sequence: two long stretches of one sign with a short stretch of the other
sign between them give no change. A run that begins with a long stretch of the sign opposite to the start has a
change at frame 0, time 0.01 ns. A run shorter than 1 ns has no long stretch and no lasting change.

Written to `n2_configuration` of `summary.json`: the counts and fractions above; `sign_first_frame` and
`sign_last_frame` (`n2_sign` of frame 0 and of frame N − 1 as written, not the assigned sign); `n_stretches` (all
stretches), `n_long_stretches` and the list `long_stretches` (sign, first frame, number of frames, the times
(first frame + 1) Δt and (last frame + 1) Δt); `n_lasting_sign_changes`, the list `lasting_sign_changes` and
`first_lasting_sign_change_ns` (the time of the first one, none without one).

Worked examples (run through the primary implementation, `unit_checks.json`; Δt = 10 ps):

- Example 1, sign of the start +1, 325 frames with V = 0.5 (frames 0 to 2), 0 (3, 4), −0.4 (5 to 124), 0 (125), 0.3 (126 to 224), −0.2 (225 to 324): 102 frames positive, 220 negative, 3 zero. The assigned signs of the frames 3, 4 and 125 are +1, +1 and −1. There are 4 stretches, 2 of them long: sign −1 from frame 5 with 121 frames and sign −1 from frame 225 with 100 frames. The stretch of 99 frames (126 to 224) is short. One lasting sign change: frame 5, 0.06 ns, from +1 to −1; the second long stretch has the sign of the one before it and gives none.
- Example 2, sign of the start −1, V = 0 in the frames 0 and 1 and 0.2 in the 121 frames after them: the assigned signs of the frames 0, 1, 2 are −1, −1, +1; one lasting sign change at frame 2, 0.03 ns, from −1 to +1; `sign_first_frame` is 0 (the sign as written, not the assigned one) and `count_zero` is 2.
- Example 3, sign of the start +1, 60 frames with V = −0.3: 0 long stretches, 0 lasting sign changes; the window holds 60 frames.

Nothing depends on this section: the states, the verdict and the rules use d1 and d2 only.

## 9B. The hydrogen bond N2–H···OG1(Thr199) before and after the first passage (reported)

Let p be the index of the first frame of the first long visit in the opposite state (section 4, "first passage";
`orientation.first_passage_frame`; `first_passage_ns` = (p + 1) Δt). With the flag `hb_N2H_T199OG1` of section 5:

- **before**: the frames 0 … p − 1. `n_before` = p, `count_before` the number of them with the flag 1,
  `fraction_before` = `count_before` / `n_before` (none if p = 0);
- **after**: the frames p … N − 1, the frame of the passage included. `n_after` = N − p, `count_after`,
  `fraction_after` = `count_after` / `n_after`.

The six numbers are written for every run that has a first passage, whatever its verdict, and the verdict is
copied beside them; they are what the plan asks for a run that converts or is mixed. A run without a first passage
(no visit of 1 ns in the opposite state; the verdict counts frames of the last 20 ns and does not need such a
visit, so a CONVERTED or MIXED run can be one) has `available` false and none for all six. Worked examples: start state C, 150 frames of state R and then 50 of state C, the flag 1 in every second frame beginning with frame 0: the first passage is frame 0 (0.01 ns), `n_before` 0, `count_before` 0, `fraction_before` none, `n_after` 200, `count_after` 100, `fraction_after` 0.5. With 200 frames of state C there is no first passage: `available` false and the six numbers none.

## 9C. The zinc to the carboxylate oxygens of Glu106 (reported)

Solute trajectory, every frame, after step 1 of section 2:

- `d_Zn_E106OE1_A` = |Zn − OE1(Glu106)|, `d_Zn_E106OE2_A` = |Zn − OE2(Glu106)| (atoms `glu106_OE1`, `glu106_OE2`),
  `d_Zn_E106OE_min_A` the smaller of the two;
- flags, 1 or 0: `fifth_E106OE1` = 1 if `d_Zn_E106OE1_A` ≤ 2.6 Å, `fifth_E106OE2` likewise for OE2,
  `fifth_E106OE` = 1 if `d_Zn_E106OE_min_A` ≤ 2.6 Å (one of the two oxygens).

The limit is the 2.6 Å of the fifth ligand (section 6), the limit included. Statistics and fractions as in
section 10; asked for by name are the mean, the smallest value and the fraction of frames within 2.6 Å of each of
the two distances. In 6RL9 the distances are 4.004 Å (OE1) and 5.470 Å (OE2).

**What does not change.** OE1 and OE2 are protein oxygens, so in the full trajectory (every 100 ps) they are
already inside `d_Zn_Oprotein_min_A`, `fifth_proteinO` and `fifth_any`, and through `fifth_any` inside the full
fraction of GV3 d (section 6). The solute fraction of GV3 d stays the fraction of `fifth_sulfonylO` alone. The
three flags of this section enter no criterion; the numbers are printed beside the gate
(`gate_GV.reported_zn_glu106`, section 11) and are no part of `gate_GV.passed`.

## 9D. Water in the place of the crystal ring (reported)

Full trajectory, every frame, of the production and of the equilibration (section 9E). Two counts are written.
The first is **the count that the driver acts on** in its hydration check (plan, changes of 6 Oct 2026, 1:03 AM
point 3, 3:02 AM point 8, 8:01 AM point 7 and 8:56 AM point 6; `md/protocol.json`, key
`equilibration.hydration_check`; appendix A.5): the water oxygens within 2.5 Å of the positions that the six ring
carbons C1 to C6 have in the crystal pose, "a water inside the place of the ring". The second is the wider count of
version 1.2, kept under another name.

**Settings of the first count** (`water_ring_place` of `summary.json`). The list P of site atoms and the cut-off c
are those of the run, so that the analysis counts what the driver of that run counted.

- If `manifest.json` has the key `hydration_check.settings` (a value of null counts as no key) **and** its value
  is an object for which both of the following hold, P and c are taken from it. (a) `site_atoms` is a list of at
  least one name, every name is, letter for letter, one of the eleven heavy-atom names of the ligand (section 0:
  C1 to C6, N1, N2, O1, O2, S), and no name occurs twice. (b) `water_cutoff_A` is a JSON number (a whole number
  counts; a text, true or false do not), finite and above zero. Then P = `site_atoms` in the order of the manifest, c = `water_cutoff_A` (Å), and `source` is
  `manifest hydration_check.settings`.
- Otherwise P = (C1, C2, C3, C4, C5, C6) and c = 2.5 Å, **the defaults of this specification**, and `source` is
  `specification (defaults)`. This is the case without a manifest, with a manifest of the plain keys, with the
  manifest of a run whose protocol had no hydration check (the driver then writes `hydration_check` without
  `settings`), and with settings that fail (a) or (b). Both values come from the same place: one is never taken
  from the manifest and the other from the defaults.
- `settings_are_the_defaults` is true if P holds exactly the six names C1 to C6, in any order, and
  |c − 2.5| ≤ 10⁻⁹ Å, wherever P and c came from; else false.

Two warnings belong to the settings; each is given at most once for a run (`n_warnings` counts them). (i) `hydration_check.settings` is there and fails (a) or (b): the defaults
are used. (ii) P and c come from the manifest and `settings_are_the_defaults` is false: the run was made with a
hydration check other than that of the plan (a test, or a run started with a protocol older than 6 Oct 2026,
8:56 AM IST). The count then follows the run, so that it can still be set beside the driver's record, and it is not
the quantity of the plan's rule.

**The counts of a full frame.**

1. **Site positions.** r_a for every a in P: the position of the one atom of that name in `reference.pdb` (the
   crystal pose; the file holds it for every start).
2. **Fit of the full frame.** The CA atoms of the full frame (the atoms `protein_CA`, found at their indices, in
   list order) are superposed on the CA atoms of `reference.pdb` exactly as in step 2 of section 2 (least squares,
   equal weights, rotation without reflection). No atom is moved by a lattice vector before this fit: the protein
   is whole in the full trajectory (section 1, "What the driver must keep"). If the fit carries a point x of the
   frame to F(x), the site positions in the frame are s_a = F⁻¹(r_a): the crystal positions carried back into the
   frame. Written out, with g_x the centroid of the CA atoms of the frame, g_r the centroid of the CA atoms of the
   reference and Q the rotation matrix of the fit for row vectors: F(x) = (x − g_x) Q + g_r and
   s_a = (r_a − g_r) Qᵀ + g_x.
3. **Distance of a water oxygen.** For every atom w of `water_O` (found at its index): d_w = the smallest of the
   minimum-image distances |w − s_a| over the a in P, in the unit cell of the frame (section 2; without a unit
   cell: the plain distances).
4. `n_water_ring_place` = the number of atoms w with d_w ≤ c. A water oxygen counts once, however many of the
   positions it is near. `d_water_ring_place_min_A` = the smallest d_w of the frame, whether or not it is within c.
5. **The wider count** (the count of version 1.2, there `n_water_ring_site_3p5A` and `d_water_ring_site_min_A`):
   steps 1 to 4 with the seven positions of C1, C2, C3, C4, C5, C6 and N1 and the cut-off 3.5 Å. Both are
   constants of this specification, whatever the manifest says. Its two numbers are `n_water_ring_wide_3p5A` and
   `d_water_ring_wide_min_A`.

The four numbers are written for every full frame, those of the production to `timeseries_full.csv` and those of
the equilibration to `timeseries_equil_full.csv`, and get the statistics of section 10 (frames of the production
only; mean over the run, blocks of 100 full frames). They are computed for every start; they are meant for
start R. The driver's own record of the first count, at the end of stage 03 and of every extension step, is set
beside it in section 9E.

**How to read them.** In the starts X and N the ring itself lies at the six positions, so the first count is zero there as long as the ring stays. In the built systems (section 16.1) `n_water_ring_place` is 0 in X, 0 in N and 0 in R, with the nearest water oxygen 3.03, 3.03 and 4.15 Å from a ring-carbon position. In X and N that nearest water is the crystal water in contact with the ring; in R nothing is nearer, so the place of the ring is empty when the run begins, which is what the rule of the driver acts on. A count of one or more in a full frame whose solute frame of the same time has the state R says that a water oxygen lies within the cut-off of a position that a ring carbon has in the crystal pose: water in the place of the ring. The wider count is 1, 1 and 1 in the three built systems (nearest 3.03, 3.03 and 2.84 Å): it also sees water beside a ring that is in place and, in R, a water of the solvent near the place of the aniline nitrogen, which points into the solvent. It is kept because version 1.2 defined it; no rule, no step of the driver and no comparison uses it.

## 9E. Equilibration frames (reported; a series of their own)

The driver writes `solute_eqNN.dcd` and `full_eqNN.dcd` for every equilibration stage (`md/RUN.md`), and files of
the same kind for every extension step of the held part (`md/run_md.py`, `md/protocol.json`). If such
files are in the run folder they are read as follows; if none is there, this section writes
`equilibration.available` false and nothing else changes.

**Steps and files.** An equilibration step has a tag: `eq` followed by digits only is stage NN (`eq01`, `eq03`);
that followed by `x` and a whole number M from 1 on, written without a leading zero, is extension step M of stage
NN (`eq03x1`, `eq03x2`). NN and M are the numbers the digits give; a stage has M = 0. The driver extends the held
part of the equilibration by the steps `eq03x1` to `eq03x4`, one at a time, as long as its hydration check applies
and finds no water oxygen in the place of the crystal ring (the first count of section 9D is zero; appendix A.5;
plan, 1:03 AM point 3, 8:01 AM point 7 and 8:56 AM point 6). The steps are
ordered by NN and, for one NN, by M: `eq01`, `eq02`, `eq03`, `eq03x1`, `eq03x2`, `eq04`, … A file belongs to a
step if its name is `solute_` or `full_`, then the tag of the step, then `.dcd`. Atoms, atom order and frame
intervals are those of the production files of the same kind (`solute.pdb`, `full.pdb`; Δt of section 1).

**Step table**: the start S and the end E of every step, in ps, counted from the start of the equilibration.

- If `manifest.json` holds a list under `equilibration.stages` (the manifest of `md/run_md.py`): every entry of
  it whose `tag` is the tag of a step as defined above and which has the numbers `time_start_ps` and `time_end_ps`
  is one step, S = `time_start_ps`, E = `time_end_ps`. Only these steps exist
  (`stage_table_source` = `manifest equilibration.stages`). Entries with another tag are passed over.
- Otherwise (no manifest, a manifest with the plain keys, the key missing or not a list): the six stages of
  `md/protocol.json`, the constants of this specification, and no extension step
  (`stage_table_source` = `protocol (constants of the specification)`):

| Stage NN | 01 | 02 | 03 | 04 | 05 | 06 |
|---|---|---|---|---|---|---|
| S (ps) | 0 | 100 | 400 | 700 | 1100 | 1500 |
| E (ps) | 100 | 400 | 700 | 1100 | 1500 | 2000 |

- T_eq is the largest E of the table used (2000 ps for the constants; with extension steps the manifest gives
  more).
- A file whose step is not in the table used is not read (warning); so the files of extension steps are not
  read without the list of the manifest, which alone says where they lie in time. If the table is empty nothing is
  read (warning).

**Frames and time.** The files of a kind are read in the order of the steps, complete frames only; `frame` counts
through them from 0. Frame k (k = 0, 1, …) of the file of a step has the time

t = (S + (k + 1) Δt − T_eq) / 1000 ns,

with S of that step and Δt the frame interval of its kind in ps. The first frame of a step is written one interval
after the start of the step; every time is zero or negative; time zero is the start of the production, and the
last frame of the last step lies at zero if the length of that step is a whole number of intervals. A step from S
to E is expected to hold m frames, m the largest whole number with m Δt ≤ (E − S) (1 + 10⁻⁹) (a step shorter than
one interval holds none: a file with a header only). Each line of the two CSV files carries `stage` (NN) and
`extension` (M, 0 for a stage).

**Quantities.**

- Solute frames: imaging and fit as in section 2, then every quantity of the solute trajectory (sections 3 to 9C).
  `timeseries_equil_solute.csv` has the columns of `timeseries_solute.csv` (section 13) except the two that depend
  on the frames before a frame, `state_assigned` and `n2_sign_assigned`.
- Full frames: the quantities of sections 6, 7 and 9D. `timeseries_equil_full.csv` has the columns of
  `timeseries_full.csv`. They are read only if `full.pdb` is in the run folder.

**Kept apart.** No equilibration frame is counted in N or in the time of a production frame, and none enters a
statistic, a fraction, a block, the RMSF, the states assigned, the visits, the first passage, the verdict, rule 1,
gate G-V or a rule of the campaign. With or without the equilibration files, the analysis of the production of a
run writes the same numbers (checked, section 16.4).

**What is reported** (`equilibration` of `summary.json`): whether frames were read; where the step table came
from and T_eq; the numbers of solute and full frames; for every step with a file, in the order of the steps, its
stage, extension number and tag, its S and E, the frames read and expected, the numbers of frames of state C, R, I
and of sign +1 and −1 at N2, the last solute frame (state, ligand RMSD, signed volume and sign at N2) and the last
full frame (the four numbers of section 9D); over all solute frames of the equilibration the same counts, and the
first frame, in order of time, whose `state` is not the state of the start, the first whose `state` is the
opposite state, and the first whose `n2_sign` is not the sign of the start (each with frame, stage, extension, time;
none if there is no such frame). The last solute frame and the last full frame of a step are the last frames
read from its two files (none if the file holds no frame). The water in the place of the ring at the end of the
held part of the equilibration (plan, 1:03 AM, point 3) is `n_water_ring_place` in `last_full_frame` of stage 03
and of each of its extension steps; that of the last of these steps is the count on which the driver went on.

**The driver's own count.** At the end of stage 03 and of every extension step the driver makes the first count
of section 9D itself, on the coordinates it saves there and with the settings of the run, and writes it into the
record of the step, which the manifest holds as the entry of that step in `equilibration.stages` (appendix A.4):
`hydration_check.n_water_O` (the count), `hydration_check.nearest_water_O.distance_A` (the smallest d_w, in Å with
three decimals), `hydration_check.applies` (true if the rule acts on the run: the ligand's own site atoms are on
average more than 2.0 Å away from the site positions) and `hydration_check.decision` (`extend` or `go on`).

For every step that has an entry in `stages` of `summary.json` (a step of the step table with a file) and whose
entry in the manifest holds `hydration_check.n_water_O` as a whole number (a JSON integer), `driver_hydration_check`
of the step is written; for every other step it is none. Its keys:

- `n_water_O`, `nearest_water_O_A`, `applies`, `decision`: the four values copied from the manifest (none where the
  manifest lacks one);
- `compared`: true if the full file of the step holds m ≥ 1 frames and the last of them is written at the end of
  the step, that is |S + m Δt − E| ≤ 10⁻⁶ ps with S and E of the step table and Δt the interval of the full frames;
  else false. With the protocol of the campaign every step is a whole number of 100 ps long and ends on a full
  frame; a shortened test run need not.
- `count_agrees`: if `compared` is true, whether `n_water_ring_place` of that last full frame equals `n_water_O`;
  none otherwise;
- `nearest_agrees`: if `compared` is true and the driver's distance is a number, whether `d_water_ring_place_min_A`
  of that frame differs from it by at most 0.002 Å; none otherwise.

A `count_agrees` or a `nearest_agrees` that is false is a warning, one for the step whether the count, the
distance or both disagree: the analysis and the driver have not counted the same thing in the same coordinates. The comparison is like for like because section 9D takes P and c from the
same manifest. The driver acts on its count only when `applies` is true; the analysis does not judge that and
reports its counts for every run. The wider count is compared with nothing.

**The first production frame against the start as built** (`start_of_production`; written for every run, with or
without equilibration files). With the state of the start of section 4 (C for X and N, R for R) and the sign of
the start of section 9A:

- `state_first_frame` = `state` of solute frame 0 of the production (time Δt); `state_differs` is true if it is
  not the state of the start (an intermediate first frame differs); `state_first_frame_is_opposite` is true if it
  is the opposite state;
- `n2_sign_first_frame` = `n2_sign` of that frame; `n2_sign_differs` is true if it is not the sign of the start
  (a sign of 0 differs);
- `n2_volume_first_frame`, `lig_rmsd_first_frame_A`: the two columns of that frame.

**Guards of the equilibration** (warnings; each text begins with `equilibration:` or names the file). They concern
the series of this section only; the numbers of the production do not depend on them. They are counted in
`n_warnings` like every warning.

| Guard | Rule |
|---|---|
| step unknown | a file of a stage or of an extension step whose step is not in the step table used |
| frames | the frames read of a file are not the m frames expected for its step |
| count of the driver | the driver's count of water in the place of the ring at the end of a step, or its nearest distance, does not agree with `n_water_ring_place` or `d_water_ring_place_min_A` of the full frame of that time (see above) |
| interval | the interval in the DCD header of a file with frames differs from Δt by more than 0.1 % |
| header | the DCD header (first step × step length; appendix A.2) puts the first frame of a file at a time t_h that is neither S + Δt nor Δt, each within 0.1 % of Δt; a header that counts from the start of its own file (t_h = Δt) tells nothing. The driver writes the first step as a whole number of time steps of the file, so t_h can lie up to half a time step (1 or 2 fs) from S + Δt when S is not a whole number of these steps. With the protocol of the campaign every S is one, and t_h is exact; a test run shortened by an odd factor can get this warning with a right time axis (seen on the driver's own tests, section 16.4) |
| unreadable | a file that cannot be read as a DCD file (an empty file) is taken as holding no frame |
| no topology | `full_eqNN.dcd` without `full.pdb`: the full frames are not read |

## 10. Statistics

For every distance, angle and count of sections 3 to 9D, and for the signed volume `n2_volume` (arithmetic mean,
like a distance), per run, over the frames of the production only:

- over the whole run: mean, standard deviation (with N in the denominator), smallest and largest value;
- over the **last half**, the last ⌊N/2⌋ frames: mean and standard deviation;
- **block means** over blocks of 10 ns (1000 solute frames or 100 full frames; block k holds the frames
  k B … (k + 1) B − 1; an incomplete last block is left out), and `block_sem`, the standard deviation of the
  block means (denominator: blocks − 1) divided by the square root of their number (none with fewer than two
  blocks).

For the two signed torsions the mean is the circular mean atan2(⟨sin⟩, ⟨cos⟩), the standard deviation the circular
one, sqrt(−2 ln R) in degrees with R = sqrt(⟨sin⟩² + ⟨cos⟩²); block means likewise; no `block_sem`.

For every flag (hydrogen bonds, fifth ligand, the three flags of Glu106 of section 9C): the number of frames with
the flag, the number of frames and their ratio, over the whole run, the last half and each block.

**Nothing is averaged over runs.** The campaign table has one line per run.

## 11. Gate G-V

Gate G-V of the plan is the gate of the **validation run**: the run of role `validation` (section 1), run `X0` of
its zinc model, 5 ns from the crystal pose. It is analysed on its own, `analyse_md.py <runs>/<model>/X0`, with no
option. The criteria below are evaluated in the same way for every other run and written to its `gate_GV`, with
`applies_to_this_run` false and `verdict` none: for a production run (X1 to X3, and N and R as well) these are
numbers without a verdict, and nothing in section 12 depends on them.

Every mean is the mean over all solute frames read of the run (`stats.<column>.mean`), every count is over all
frames read of its trajectory. The mean of `tau_fold_deg`, which lies in 0–180°, is the arithmetic mean, not a
circular one. The comparisons are made with the means as they are written to `summary.json` (six decimals).

| Criterion | Part | Rule |
|---|---|---|
| GV1 | | mean `lig_rmsd_A` ≤ 1.0 Å |
| GV2 | | `hb_N2H_T199OG1` in at least 80 % of the solute frames: 5 × count ≥ 4 × N |
| GV3, zinc geometry "as in G-Z" | a | 1.83 Å ≤ mean `d_Zn_N2_A` ≤ 2.03 Å |
| | b | 1.91 Å ≤ mean ≤ 2.12 Å for each of `d_Zn_H94NE2_A`, `d_Zn_H96NE2_A`, `d_Zn_H119ND1_A` |
| | c | each of the six mean N–Zn–N angles within 10.0° of its own reference, the value of that angle in 6RL9 (table below) |
| | d | no fifth ligand: the solute fraction and the full fraction of section 6 are each at most 5 % (20 c_s ≤ N_s and N_f ≥ 1 and 20 c_f ≤ N_f) |
| GV4, ligand internal coordinates | a | mean `ang_C4_S_N2_deg` within 3.0° of 108.8°: 105.8° ≤ mean ≤ 111.8° |
| | b | mean `tau_fold_deg` within 20.0° of 132.1°: 112.1° ≤ mean ≤ 152.1° |

GV3 is met if its part a, all three distances of part b, all six angles of part c and part d are met; GV4 if both
of its parts are met.
`gate_GV.passed` is true if GV1, GV2, GV3 and GV4 are all met. `gate_GV.verdict` is `PASSED` or `FAILED` after
`passed` for the validation run, and none for every other run. `gate_GV.provisional` is true if the run is not
complete (fewer frames than planned): the verdict of an unfinished validation run is provisional.

**The six references of GV3 c.** Each is the angle of the crystal start X, that is of 6RL9: Zn 301
(`md/system/zinc.pdb`), N2 of SAN 304 (`md/system/ligand_X.pdb`), NE2 of His94, NE2 of His96 and ND1 of His119
(`md/system/protein_md.pdb`, deposited coordinates), rounded to one decimal as the plan quotes them.

| Column | Angle | In 6RL9, from `md/system/` (°) | Reference (°) | The mean must lie in (°) |
|---|---|---|---|---|
| `ang_N2_Zn_H94NE2_deg` | N2 of SAN 304 – Zn 301 – NE2 of His94 | 108.920 | 108.9 | 98.9 to 118.9 |
| `ang_N2_Zn_H96NE2_deg` | N2 of SAN 304 – Zn 301 – NE2 of His96 | 114.280 | 114.3 | 104.3 to 124.3 |
| `ang_N2_Zn_H119ND1_deg` | N2 of SAN 304 – Zn 301 – ND1 of His119 | 116.833 | 116.8 | 106.8 to 126.8 |
| `ang_H94NE2_Zn_H96NE2_deg` | NE2 of His94 – Zn 301 – NE2 of His96 | 103.819 | 103.8 | 93.8 to 113.8 |
| `ang_H94NE2_Zn_H119ND1_deg` | NE2 of His94 – Zn 301 – ND1 of His119 | 113.591 | 113.6 | 103.6 to 123.6 |
| `ang_H96NE2_Zn_H119ND1_deg` | NE2 of His96 – Zn 301 – ND1 of His119 | 98.409 | 98.4 | 88.4 to 108.4 |

The last two columns are the constants of this specification. An implementation takes these numbers as printed; it
does not compute them from a structure file. The limits of GV3 a and b are the literal numbers of the plan
("within 0.10 Å of 1.93 Å", "within 0.10 Å of 2.01–2.02 Å"); those of GV4 are the plan's 108.8° ± 3° and
132.1° ± 20° (the crystal start has 108.795° and 132.056°).

GV3 is the zinc geometry of gate G-Z of the plan (stage S2, after minimisation and 100 ps), which names the same
four parts. This specification defines them once; which run of stage S2 they are applied to is not laid down here.

**Reported beside the gate, no criterion** (version 1.2; plan, 6 Oct 2026, 3:02 AM, point 4): for each of
`d_Zn_E106OE1_A` and `d_Zn_E106OE2_A` the mean, the smallest value, and the number and fraction of solute frames
within 2.6 Å (section 9C), under `gate_GV.reported_zn_glu106`. These numbers are written for every run. They are
no part of GV3, of `gate_GV.passed` or of the verdict: a run in which OE1 of Glu106 is within 2.6 Å of the zinc in
more than 5 % of the solute frames and in no full frame still meets GV3 d (tested, section 16.2, run glu/X0).

## 12. Campaign

`analyse_md.py --campaign DIR` analyses every folder below DIR that holds a `solute_segNNNN.dcd` and whose name
is `<start><replica>` with replica 1, 2, … (X, N or R followed by a number that does not begin with 0). Other
folders with a trajectory (the validation run `X0`, a test) are listed as left out (`skipped`) and are not analysed
by this call. A folder with a production name whose run has the role `validation` by its manifest is analysed and
then left out in the same way. The validation run is therefore in no campaign table, in no campaign figure and in
no rule below; its gate is read from its own `summary.json`. Runs are ordered by model, start (X, N, R) and
replica. The rules of section 5 of the plan are evaluated per zinc model:

- **Rule 1** fails if `rule1_left_crystal` is true for two or three runs of start X.
- **Outcome** (rule 3), in this order:
  1. `validation_failed` if rule 1 fails;
  2. `rotated_converts` if at least two runs of start R are CONVERTED and every run of start X and N is STAYED;
  3. `both_persist` if every run is STAYED;
  4. `mixed` otherwise: the runs are reported one by one.
- **Rule 4:** `rule4_sensitivity_run_required` is true exactly for `rotated_converts`.
- The outcome is **final** only if the model has three complete runs of each start; otherwise it is marked as
  not final.

A run that is MIXED has not "stayed"; a campaign without a converted run but with a mixed one is therefore `mixed`,
not `both_persist`.

## 13. Output

**Per run**, in `<run folder>/analysis/` (or `--out DIR`):

| File | Content |
|---|---|
| `timeseries_solute.csv` | one line per solute frame: `frame`, `segment`, `time_ns`, then the columns listed below |
| `timeseries_full.csv` | one line per full frame: `frame`, `segment`, `time_ns`, `d_Zn_Owater_min_A`, `water_O_index_nearest_Zn`, `d_Zn_Oprotein_min_A`, `protein_O_index_nearest_Zn`, `d_Zn_Osulfonyl_min_full_A`, `d_Zn_O_any_min_A`, `fifth_water`, `fifth_proteinO`, `fifth_any`, `n_water_N1_3p5A`, `d_Zn_N2_full_A`, since version 1.2 `n_water_ring_wide_3p5A` (whole number), `d_water_ring_wide_min_A` (four decimals; these two were named `n_water_ring_site_3p5A` and `d_water_ring_site_min_A` in version 1.2), and since version 1.3, at the end, `n_water_ring_place` (whole number), `d_water_ring_place_min_A` (four decimals) |
| `timeseries_equil_solute.csv` | only if a solute file of a step of the step table is in the run folder: one line per solute frame of the equilibration: `frame`, `stage`, `extension`, `time_ns` (six decimals, zero or negative), then the columns of `timeseries_solute.csv` without `state_assigned` and `n2_sign_assigned` (section 9E) |
| `timeseries_equil_full.csv` | only if a full file of a step of the step table and `full.pdb` are there: `frame`, `stage`, `extension`, `time_ns` (six decimals), then the columns of `timeseries_full.csv` |
| `rmsf_ca.csv` | `residue_number_deposited`, `residue_name`, `rmsf_A` |
| `summary.json` | every number of the run (keys below) |
| `fig_orientation.png` | d1 and d2 against time with the midpoints (and N1 to ND2 of Asn67 where it is reported); the state of each frame |
| `fig_ligand.png` | ligand RMSD against time, with the 2.0 Å and 1.0 Å lines |
| `fig_zinc.png` | the four Zn–N distances; the six angles; the nearest sulfonyl, water and protein oxygen |
| `fig_contacts.png` | the three donor–acceptor distances; the three hydrogen bonds as % of frames (running mean over 1 ns); water at N1 |
| `fig_internal.png` | τ_fold; the two signed torsions; C4–S–N2; S–N2 |
| `fig_protein.png` | CA RMSD against time; CA RMSF per residue |
| `fig_n2_site.png` | since version 1.2: the signed volume at N2; the zinc to OE1 and OE2 of Glu106; the two water counts of section 9D; the distance of the nearest water oxygen to the place of the ring, with the cut-off |
| `fig_equilibration.png` | since version 1.2, only with equilibration frames: d1 and d2, ligand RMSD, signed volume at N2 and the two water counts of section 9D against the negative time |

A file `timeseries_equil_*.csv` or `fig_equilibration.png` left in the output folder by an earlier analysis is
removed when the run folder holds no equilibration frames for it. A file with a header line only is written when
the stage files hold no frame.

Columns of `timeseries_solute.csv` after `frame`, `segment`, `time_ns`, in this order (time with three decimals,
distances with four, angles with three, flags 0 or 1):

`lig_rmsd_A`, `lig_rmsd_labelling`, `lig_rmsd_identity_A`, `d1_N1_F131CZ_A`, `d2_N1_N67OD1_A`, `state`,
`state_assigned`, `d_N2_T199OG1_A`, `ang_N2_HN2_T199OG1_deg`, `hb_N2H_T199OG1`, `d_Zn_N2_A`, `d_Zn_H94NE2_A`,
`d_Zn_H96NE2_A`, `d_Zn_H119ND1_A`, `ang_N2_Zn_H94NE2_deg`, `ang_N2_Zn_H96NE2_deg`, `ang_N2_Zn_H119ND1_deg`,
`ang_H94NE2_Zn_H96NE2_deg`, `ang_H94NE2_Zn_H119ND1_deg`, `ang_H96NE2_Zn_H119ND1_deg`, `ang_S_N2_Zn_deg`,
`d_Zn_O1_A`, `d_Zn_O2_A`, `d_Zn_Osulfonyl_min_A`, `fifth_sulfonylO`, `d_T199N_O1_A`, `d_T199N_O2_A`,
`ang_T199N_H_O1_deg`, `ang_T199N_H_O2_deg`, `d_T199N_Osulfonyl_min_A`, `hb_T199NH_Osulfonyl`,
`d_T199OG1_E106OE1_A`, `d_T199OG1_E106OE2_A`, `ang_T199OG1_HG1_E106OE1_deg`, `ang_T199OG1_HG1_E106OE2_deg`,
`d_T199OG1_E106OE_min_A`, `hb_T199OG1H_E106OE`, `tors_C3_C4_S_N2_deg`, `tors_C5_C4_S_N2_deg`, `tau_fold_deg`,
`ang_C4_S_N2_deg`, `d_S_N2_A`, `ca_rmsd_A`,

and since version 1.2, in this order: `n2_volume` (four decimals), `n2_sign` (1, -1 or 0), `n2_sign_assigned` (1 or
-1), `d_Zn_E106OE1_A`, `d_Zn_E106OE2_A`, `d_Zn_E106OE_min_A`, `fifth_E106OE1`, `fifth_E106OE2`, `fifth_E106OE`.

Only if `selections.json` has the key `asn67_ND2`, two more columns follow at the end, after `fifth_E106OE`:
`d_N1_N67ND2_A` (four decimals), `asn67_ND2_nearer_N1` (0 or 1). Without the key the file ends with `fifth_E106OE`.

`segment` is the NNNN of the file the frame came from. The indices in `timeseries_full.csv` are atom indices of
the topology (zero-based).

**Keys of `summary.json`** (numbers rounded to six decimals; "none" is JSON null):

```
spec, script{name, sha256}, software{python, numpy, MDAnalysis, matplotlib}, random_numbers, units
run{run_dir, label, model, start, replica, role, system_dir, planned_ns, seen_ns, complete, solute_interval_ps,
    full_interval_ps, n_frames_solute, n_frames_full, seen_full_ns, n_atoms_solute, n_CA,
    segments_solute[], segments_full[]   each {index, file, n_frames, bytes, header_first_step,
                                               header_step_interval, header_interval_ps, header_periodic}
    header_check_solute{interval_matches, contiguous, details[]}, header_check_full{...},
    segment_overlap_check{suspected, boundaries[]}, manifest_segment_check{available, consistent, details[]},
    sources{...}, system_files_match_manifest{selections.json, reference.pdb},
    imaging_solute{frames_with_box, frames_without_box, frames_zinc_moved, frames_ligand_moved},
    full_trajectory{n_water_O, n_protein_O, frames_without_box, name_mismatches[]},
    largest_neighbour_CA_distance_A, consistency_full_vs_solute{ratio_full_to_solute_interval, n_compared,
    max_abs_diff_d_Zn_N2_A}}
inputs_sha256{selections.json, reference.pdb, solute.pdb, full.pdb, manifest.json, every DCD file}
constants{...}                      the constants of this specification as used
stats{<column>{kind, n, mean, sd, min, max, n_last_half, mean_last_half, sd_last_half, block_frames,
               block_means[], block_sem, trajectory}}
      for every distance, angle and torsion column of the two CSV files, for n2_volume (kind "vol"), for
      n_water_N1_3p5A, n_water_ring_wide_3p5A and n_water_ring_place (kind "count"; d_N1_N67ND2_A only with
      the key asn67_ND2)
fractions{<flag>{n, count, fraction, n_last_half, count_last_half, fraction_last_half, block_frames,
                 block_fractions[]}}
      for hb_N2H_T199OG1, hb_T199NH_Osulfonyl, hb_T199OG1H_E106OE, fifth_sulfonylO, fifth_water,
      fifth_proteinO, fifth_any, fifth_E106OE1, fifth_E106OE2, fifth_E106OE (and asn67_ND2_nearer_N1 only with
      the key asn67_ND2)
lig_rmsd_labelling_counts{identity, ring_flip, oxygen_swap, ring_flip+oxygen_swap}
orientation{start_state, opposite_state, midpoint_d1_A, midpoint_d2_A, window_frames, window_ns,
            window_complete, count_C_window, count_R_window, count_I_window, fraction_start_state_window,
            fraction_opposite_state_window, fraction_intermediate_window, verdict, provisional, dwell_frames,
            first_passage_ns, first_passage_frame, n_transitions, n_visits, n_long_visits,
            long_visits[]{state, first_frame, n_frames, t_first_ns, t_last_ns, long, fraction_strict},
            state_first_frame, count_C, count_R, count_I, fraction_C, fraction_R, fraction_I, n_last_half,
            count_C_last_half, count_R_last_half, count_I_last_half}
asn67_ND2{available, note}          without the key asn67_ND2; with it also:
          {n, count_ND2_nearer, fraction_ND2_nearer,
           by_state{C{n, count_ND2_nearer}, R{n, count_ND2_nearer}, I{n, count_ND2_nearer}}}
          n: solute frames (in by_state: the frames whose `state`, not `state_assigned`, is C, R or I);
          count_ND2_nearer: those of them with the flag 1; fraction_ND2_nearer = count_ND2_nearer / n
n2_configuration{order, sign_of_start, n, count_positive, count_negative, count_zero, fraction_positive,
                 fraction_negative, window_frames, count_positive_window, count_negative_window,
                 count_zero_window, fraction_positive_window, fraction_negative_window,
                 by_state{C{n, count_positive, count_negative, count_zero, fraction_positive, fraction_negative},
                          R{...}, I{...}},
                 sign_first_frame, sign_last_frame, dwell_frames, n_stretches, n_long_stretches,
                 long_stretches[]{sign, first_frame, n_frames, t_first_ns, t_last_ns, long},
                 n_lasting_sign_changes, lasting_sign_changes[]{frame, time_ns, from_sign, to_sign},
                 first_lasting_sign_change_ns}                                                     section 9A
hbond_first_passage{available, flag, verdict, first_passage_frame, first_passage_ns, n_before, count_before,
                    fraction_before, n_after, count_after, fraction_after}                         section 9B
water_ring_place{note, site_atoms[], water_cutoff_A, source, settings_are_the_defaults}           section 9D
                 the site atoms P and the cut-off c used for n_water_ring_place and d_water_ring_place_min_A;
                 source: "manifest hydration_check.settings" or "specification (defaults)"
start_of_production{note, time_ns, state_of_start, state_first_frame, state_differs,
                    state_first_frame_is_opposite, n2_sign_of_start, n2_sign_first_frame, n2_sign_differs,
                    n2_volume_first_frame, lig_rmsd_first_frame_A}                                 section 9E
equilibration{available, note, files_solute[], files_full[]}      without equilibration files; with them also:
             {stage_table_source, total_ps, n_frames_solute, n_frames_full, files_written[],
              imaging_solute{...}                                 as run.imaging_solute; only with solute frames
              state_of_start, n2_sign_of_start,
              stages[]{stage, extension, tag, time_start_ps, time_end_ps, time_start_ns, time_end_ns,
                       solute{file, n_frames, bytes, header_first_step, header_step_interval, header_interval_ps,
                              header_periodic, tag, stage, extension, n_frames_expected, header_first_frame_ps},
                       full{...},                                 as solute; none without a file of that kind
                       count_C, count_R, count_I, count_n2_positive, count_n2_negative,
                       last_solute_frame{frame, time_ns, state, lig_rmsd_A, n2_volume, n2_sign},
                       last_full_frame{frame, time_ns, n_water_ring_place, d_water_ring_place_min_A,
                                       n_water_ring_wide_3p5A, d_water_ring_wide_min_A},
                       driver_hydration_check{n_water_O, nearest_water_O_A, applies, decision, compared,
                                              count_agrees, nearest_agrees}},      none without a count in the manifest
              count_C, count_R, count_I, count_n2_positive, count_n2_negative,
              first_frame_not_in_state_of_start{frame, stage, extension, time_ns, state, n2_sign},
              first_frame_in_opposite_state{...}, first_frame_n2_sign_not_of_start{...},
              state_last_frame, n2_sign_last_frame, time_last_solute_frame_ns}                     section 9E
              one entry of stages[] per step with a file, stages and extension steps alike, in the order of the
              steps; total_ps is T_eq; times in ns from the start of the production; header_first_frame_ps is
              first step x step length of the DCD header and is written only for a file with frames
validation{lig_rmsd_mean_A, lig_rmsd_sd_A, lig_rmsd_mean_last_half_A, rule1_threshold_A, rule1_left_crystal,
           rule1_applies_to_this_run, hb_N2H_T199OG1_fraction}
gate_GV{applies_to_this_run, note,
        GV1_lig_rmsd{mean_A, limit_max_A, passed},
        GV2_hbond{fraction, count, n, limit_min, passed},
        GV3_zinc{criteria{d_Zn_N2_A{mean_A, window_A[2], passed},
                          d_Zn_H94NE2_A{...}, d_Zn_H96NE2_A{...}, d_Zn_H119ND1_A{...}      as d_Zn_N2_A
                          <each of the six angle columns>{mean_deg, reference_6RL9_deg, tolerance_deg,
                                                          deviation_deg, window_deg[2], passed},
                          fifth_ligand{limit_max_fraction,
                                       solute_trajectory{flag, count, n, fraction, passed},
                                       full_trajectory{flag, available, count, n, fraction, passed},
                                       passed}},
                 passed},
        GV4_ligand_internal{criteria{ang_C4_S_N2_deg{mean_deg, sd_deg, reference_deg, tolerance_deg,
                                                     deviation_deg, window_deg[2], passed},
                                     tau_fold_deg{...}},                                   as ang_C4_S_N2_deg
                            passed},
        passed, verdict, provisional,
        reported_zn_glu106{note, distance_max_A,
                           d_Zn_E106OE1_A{mean_A, min_A, count_within, n, fraction_within},
                           d_Zn_E106OE2_A{...}}}                   reported, no criterion (section 9C)
        deviation_deg = mean minus reference; window = [reference - tolerance, reference + tolerance];
        "passed" of a part, of a criterion and of the gate: the rule of section 11 is met by the numbers of
        this run (computed for every run); verdict: "PASSED" or "FAILED" for the validation run, none otherwise;
        full_trajectory.count is none without a full trajectory
rmsf{file, n_residues, mean_A, max_A, max_residue_number_deposited}
warnings[]
table_row{...}                      the line of the run in the campaign table
```

**Campaign**, in `DIR/analysis_campaign/` (or `--out DIR`):

| File | Content |
|---|---|
| `campaign_table.csv` | one line per run, the keys of `table_row` as columns |
| `campaign_table.md` | the outcome per model and a shorter table for reading |
| `campaign.json` | `spec`, `script`, `software`, `campaign_dir`, `rules_per_model{<model>{n_runs, n_complete, runs_per_start, all_nine_complete, final, rule1{n_X_runs, n_X_left_crystal, failed}, n_converted, n_stayed, n_mixed, outcome, rule4_sensitivity_run_required}}`, `outcome_text`, `runs[]{table_row, summary_file, summary_sha256}`, `failed[]`, `skipped[]` |
| `campaign_orientation.png`, `campaign_ligand_rmsd.png` | one panel per run (rows: starts; columns: replicas); with several models the model name is appended to the file name |

Columns of the campaign table (`table_row`): `run`, `model`, `start`, `replica`, `seen_ns`, `planned_ns`,
`complete`, `verdict`, `provisional`, `frac_start_state_last20ns`, `frac_opposite_state_last20ns`,
`frac_intermediate_last20ns`, `first_passage_ns`, `n_transitions`, `d1_N1_F131CZ_mean_last_half_A`,
`d2_N1_N67OD1_mean_last_half_A`, `d_N1_N67ND2_mean_last_half_A`, `asn67_ND2_nearer_N1_fraction` (both empty
without the key `asn67_ND2`), `lig_rmsd_mean_A`, `lig_rmsd_sd_A`, `lig_rmsd_mean_last_half_A`,
`lig_rmsd_sd_last_half_A`, `rule1_left_crystal`, `hb_N2H_T199OG1_fraction`, `hb_T199NH_Osulfonyl_fraction`,
`hb_T199OG1H_E106OE_fraction`, `d_Zn_N2_mean_A`, `d_Zn_N2_sd_A`, `d_Zn_H94NE2_mean_A`, `d_Zn_H96NE2_mean_A`,
`d_Zn_H119ND1_mean_A`, `ang_N_Zn_N_smallest_mean_deg`, `ang_N_Zn_N_largest_mean_deg` (the smallest and the
largest of the six mean angles), `ang_N_Zn_N_largest_abs_deviation_deg` (the largest of the six values
|mean − reference| of section 11), `d_Zn_Osulfonyl_min_mean_A`, `fifth_sulfonylO_fraction`,
`d_Zn_Owater_min_mean_A`, `d_Zn_Owater_min_min_A`, `fifth_any_fraction`, `n_water_N1_3p5A_mean`,
`tau_fold_mean_deg`, `tau_fold_sd_deg`, `ang_C4_S_N2_mean_deg`, `ang_C4_S_N2_sd_deg`, `d_S_N2_mean_A`,
`d_S_N2_sd_A`, `ca_rmsd_mean_A`, `ca_rmsd_mean_last_half_A`, `rmsf_ca_max_A`, `rmsf_ca_max_residue`,
`gate_GV_passed`, `n_warnings`, and since version 1.2, in this order: `n2_fraction_positive`,
`n2_fraction_positive_last20ns`, `n2_n_lasting_sign_changes`, `n2_first_lasting_sign_change_ns`,
`hb_N2H_T199OG1_fraction_before_first_passage`, `hb_N2H_T199OG1_fraction_after_first_passage`,
`d_Zn_E106OE1_mean_A`, `d_Zn_E106OE1_min_A`, `fifth_E106OE1_fraction`, `d_Zn_E106OE2_mean_A`, `d_Zn_E106OE2_min_A`,
`fifth_E106OE2_fraction`, `n_water_ring_wide_3p5A_mean`, `d_water_ring_wide_min_min_A` (the smallest value of
`d_water_ring_wide_min_A` over the run; these two were named `n_water_ring_site_3p5A_mean` and
`d_water_ring_site_min_min_A` in version 1.2), `state_first_frame`, `state_first_frame_differs_from_start`,
`n2_sign_first_frame`, `equilibration_frames_solute` (0 without equilibration frames), and since version 1.3, at
the end: `n_water_ring_place_mean` (the mean of `n_water_ring_place` over the full frames of the run),
`d_water_ring_place_min_min_A` (the smallest value of `d_water_ring_place_min_A` over the run).
Truth values are written `True` and `False`; a missing value is an empty field.
`gate_GV_passed` is `gate_GV.passed` of the run; in the campaign table, which holds production runs only, it is a
number without a verdict (section 11). The `summary.json` of the validation run has a `table_row` with the same
keys; it enters no table.

**Command line** of the primary implementation:

```
analyse_md.py RUN_DIR [--system DIR] [--start X|N|R] [--model NAME] [--planned-ns NS] [--solute-interval-ps PS]
              [--full-interval-ps PS] [--out DIR] [--no-figures] [--no-hash-trajectories] [--skip-name-check]
analyse_md.py --campaign DIR [--out DIR] [--reuse]        (--reuse takes an existing analysis/summary.json of a run
                                                           if it was written after this version of the specification)
analyse_md.py --check-starts [--out DIR]
```

The validation run needs no option: `analyse_md.py /root/cqc-sfa-md/runs/<model>/X0`. The equilibration files need
none either: they are read when they are in the run folder.

The analysis uses no random numbers. A second run on the same input gives the same files, byte for byte
(checked, section 16.4). Trajectories and the per-run files stay in `/root/cqc-sfa-md/`; the campaign table, the
JSON files and the figures are small enough for git (a time series of 10,000 frames is about 3 MB).

## 14. Agreement of two implementations

Compared key by key in `summary.json` (`stats`, `fractions`, `lig_rmsd_labelling_counts`, `orientation`,
`asn67_ND2`, `n2_configuration`, `hbond_first_passage`, `water_ring_place`, `start_of_production`, `equilibration`
without the file and header entries of its stages, `validation`, `gate_GV`, `rmsf`, `table_row`, the frame counts, `seen_ns`,
`planned_ns`, `complete` and `role` of `run`), line by line in the three CSV files and in the two CSV files of the
equilibration, and in `campaign.json` (`rules_per_model`, and the folders in `skipped`):

| Kind | Must agree within |
|---|---|
| distances, RMSD, RMSF and their means, standard deviations, block means | 0.001 Å |
| angles and torsions (torsions modulo 360°) and their statistics | 0.01° |
| the signed volume at N2 and its statistics | 0.001 |
| times | 0.000001 ns |
| numbers of frames, flags, states, labellings, visits, transitions, verdicts, outcomes, roles, `passed` of every part of gate G-V | identical |
| signs at N2, stretches, lasting sign changes with their frames, the frame of the first passage, the counts before and after it, water counts, the stage, the extension step and the frame of every equilibration line, the truth values of `start_of_production` and of `driver_hydration_check`, the site atoms, the cut-off, `source` and `settings_are_the_defaults` of `water_ring_place` | identical |
| fractions | follow from identical counts |

Coordinates are stored in single precision, so two correct programs can differ by about 0.0001 Å. A frame whose
value lies within the tolerance of a threshold (3.5 Å, 135°, 2.6 Å, 6.7 Å, 5.2 Å; for `asn67_ND2_nearer_N1` a
frame whose two distances differ by less than 0.001 Å; for `n2_sign` a frame with |V| below 0.001; for the water
counts of sections 7 and 9D a water oxygen within 0.001 Å of the cut-off of the count: 3.5 Å, and for
`n_water_ring_place` the cut-off c of the run, 2.5 Å by default) may be classed differently; such frames are listed, counts
may differ by at most their number, and a verdict that changes with them is reported as it is with both readings.
The same holds for a mean that lies within the tolerance of a limit of gate G-V (section 11) or of rule 1 (2.0 Å),
and for a count that sits at the limit of GV2 or of a fifth-ligand fraction because of such frames: the part is
listed and the gate is reported with both readings. A stretch at N2 or a visit whose length changes with such
frames is listed in the same way, and a lasting sign change that depends on them is reported with both readings.
In `driver_hydration_check`, `count_agrees` may differ between two programs only through such a water oxygen in the
frame compared, and `nearest_agrees` only if the difference from the driver's distance lies within 0.0001 Å of
0.002 Å; both cases are listed.
Times are compared as numbers (a time of zero may be written `0.000000` or `-0.000000`).
Paths, checksums, software versions, the `note` texts and the wording of warnings are not compared.

## 15. What the analysis does not do

- The orientation coordinate separates two states only. A third orientation shows as intermediate frames, a
  MIXED verdict and no transition.
- A turn of the Asn67 amide (OD1 and ND2 exchanged) changes d2 by about 0.7 Å in the rotated pose and 1.2 Å in
  the crystal pose (`md/SYSTEM.md`, section 9). Neither crosses the midpoint of 5.2 Å, and the states do not see
  the turn. With the key `asn67_ND2` it can be read from the flag of section 4, state by state; without the key it
  is not detected.
- The hydrogen bonds are geometric criteria; no energies. Water is seen every 100 ps only.
- Fractions of time in a state are not populations: nine runs of 100 ns do not sample an equilibrium, and no
  free energy is derived from them (plan, rule 5).
- The state of His64 and the charge of Glu106 are fixed by the system, not examined here.
- Gate G-Z (stage S2) is not evaluated by this program as a gate of its own; its zinc criteria are parts a to d
  of GV3 (section 11).
- The numbers of gate G-V of a production run are not a second validation: the gate is decided on run X0 alone.
- The sign at N2 (section 9A) says on which side of the plane of its neighbours the nitrogen is in a frame. It is
  not a barrier and not a population, and with the hydrogen placed by the force field it is only as good as the
  angle terms at N2 (plan, gate G-R). No rule uses it.
- `n_water_ring_place` (section 9D) counts the water oxygens within the cut-off of six points, the crystal
  positions of the ring carbons. It is not an occupancy of the volume of the ring and says nothing of how such a
  water is held. Read it with the state of the solute frame of the same time and with the distance of the nearest
  water oxygen. The wider count also sees water that lies beside a ring that is in place, and water near the place
  of the aniline nitrogen. Water is seen every 100 ps only.
- The equilibration series (section 9E) follows the manifest: the extension steps of the held part (plan, 1:03 AM,
  point 3) are read where the manifest lists them, and without such a list only the six stages of the protocol are
  known. The analysis does not apply the driver's rule (extend or go on); it reports the counts, and the driver's
  count beside its own. The driver's count is read at the end of stage 03 and of the extension steps only; the
  count that the driver also notes at the end of every other step (appendix A.4, `monitor.crystal_ring_site`) is
  not read. A run folder with equilibration files and no production frame is not analysed.
- Glu106 (section 9C) is followed by its two carboxylate oxygens only; another protein oxygen near the zinc shows
  in `d_Zn_Oprotein_min_A` of the full trajectory, with its atom index.

## 16. Test of the primary implementation

All files of the test are under
`/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/mdgr2/spec13/`
(throwaway; `run_all.sh` there repeats everything, log `run_all.log`); the reports and the scripts are copied to
`md/analysis/test_reports/`. The tests were run on the two scripts as they stood in that folder
(`dev/md/analysis/`); the files in `md/analysis/` are these files, byte for byte (checksums below), put there
together with this text. To repeat by hand, with `PY=/root/cqc-sfa-tools/md/bin/python` and an empty folder `T`:

```
$PY md/analysis/analyse_md.py --check-starts --out T/check_starts
$PY md/analysis/make_test_trajectory.py --out T/test
$PY md/analysis/analyse_md.py --campaign T/test/runs
$PY md/analysis/analyse_md.py T/test/runs/test/X0
$PY md/analysis/analyse_md.py T/test/runs_overlap/test/X9
for d in T/test/runs_gv/*/*/ T/test/runs_v12/*/*/ T/test/runs_v13/*/*/; do $PY md/analysis/analyse_md.py $d; done
$PY md/analysis/make_test_trajectory.py --compare T/test
```

Scripts tested (sha256): `analyse_md.py` `f75e6155de7007000b8c118b5e9989af7fbaf40994b3f684d4c9b3c9bbcf8ea3`, `make_test_trajectory.py` `b98c4a8d18a2b8326ba3a25b0393335cd73dad572ade0f40fcfd6edb46b2a14a`.

### 16.1 The three starts against `md/system/starts.csv`

Each start (protein CA atoms, site atoms, zinc, ligand; `md/system/`) was rotated and shifted as a whole and then analysed with the functions of the run analysis. File: `check_starts/check_starts.json`.

| Start | Ligand RMSD here (Å) | Labelling | `starts.csv` (Å) | RMSD with labelling 0 / 1 / 2 / 3 (Å) | d1 (Å) | d2 (Å) | N1 to ND2 of Asn67 (Å) | ND2 nearer than OD1 | State | CA RMSD after the fit (Å) |
|---|---|---|---|---|---|---|---|---|---|---|
| X | 0.0000 | identity | 0.0 | 0.000 / 1.454 / 1.043 / 1.790 | 3.8966 | 7.7443 | 6.5511 | yes | C | below 0.000001 |
| N | 0.6104 | identity | 0.61 | 0.610 / 1.565 / 1.220 / 1.888 | 3.3561 | 7.7255 | 6.5574 | yes | C | below 0.000001 |
| R | 3.5732 | identity | 3.573 | 3.573 / 3.762 / 3.691 / 3.874 | 9.5849 | 2.5968 | 3.3102 | no | R | below 0.000001 |

All 51 of 51 values agree (17 columns of `starts.csv` for three starts: `C3_C4_S_N2`, `C4_S_N2`, `C5_C4_S_N2`, `N1_Asn67ND2`, `N1_Asn67OD1`, `N1_Phe131CZ`, `N2_H_OG1`, `N2_OG1`, `S_N2`, `S_N2_Zn`, `T199N_O_min`, `Zn_N2`, `Zn_O1`, `Zn_O2`, `Zn_O_min`, `rmsd_to_crystal_A`, `tau_fold`); largest difference 0.0005 Å and 0.037° (tolerance 0.0015 Å and 0.06°, because `starts.csv` and the coordinate files have three decimals). The signs of the torsions agree. The exact midpoints from the two starts are 6.741 and 5.171 Å.

Crystal zinc site by the same functions: Zn–N(His) 2.010, 2.021, 2.018 Å; the six angles 108.920, 114.280, 116.833, 103.819, 113.591, 98.409° (order of section 6). A separate short program with plain formulas (`crystal_site.py`, arccosine of the scalar product) gives the same six angles from `zinc.pdb`, `ligand_X.pdb` and `protein_md.pdb`, and the same again from `tleap/protein_H.pdb`. Each of the eight references of section 11 (six angles, 108.8° and 132.1°) is the value of the crystal start to one decimal (largest difference 0.044°), and the crystal start as built lies inside every window of GV3 a to c and of GV4. All three hydrogen bonds of sections 5 and 7 are present in the starts X and N; in R the Thr199 N–H···O(sulfonyl) bond is absent (angle 118.5°).

**Versions 1.2 and 1.3: the reported quantities on the starts.** Signed volume at N2 (section 9A), the zinc to Glu106 (9C) and the crystal water at the crystal positions of the ring (9D; `md/system/waters_<start>.pdb`, the water kept in each start before solvation; the count that the driver acts on with the defaults, and the wider count):

| Start | Signed volume at N2 | Review | Sign | Sign of the start in section 9A | Zn to OE1 of Glu106 (Å) | Zn to OE2 (Å) | Crystal water oxygens within 2.5 Å of the six ring-carbon positions | Nearest one (Å) | Wider count: within 3.5 Å of the seven positions | Nearest one (Å) |
|---|---|---|---|---|---|---|---|---|---|---|
| X | +0.7240 | +0.72 | +1 | +1 | 4.004 | 5.470 | 0 of 215 | 3.025 (water 553) | 1 of 215 | 3.025 (water 553) |
| N | +0.7855 | +0.79 | +1 | +1 | 4.004 | 5.470 | 0 of 215 | 3.025 (water 553) | 1 of 215 | 3.025 (water 553) |
| R | −0.7104 | −0.71 | −1 | −1 | 4.004 | 5.470 | 0 of 211 | 4.148 (water 582) | 0 of 211 | 3.802 (water 582) |

The three signed volumes agree with the numbers of the review within 0.0045 (the review gives two decimals), and their signs are the signs of the start of section 9A. The start X is the crystal structure, so the first line is 6RL9 itself.

**The built systems as runs of one frame** (`built_starts.py`, `built_starts.json`). No trajectory of hCA II exists yet, so the three built systems of the restrained-zinc model (`/root/cqc-sfa-md/systems/restrained/`, 40902, 40911 and 40857 atoms) were written as runs of one frame, `system.inpcrd` as built, with the DCD and PDB writers of `md/run_md.py`, and analysed without a warning. This is no simulation; it checks the reading of the real `selections.json` and `reference.pdb` and the reported quantities on real coordinates. The one-frame manifests name no settings of a hydration check, so the count of the driver is made with the defaults of section 9D. A separate program with plain formulas (own reading of the Amber files, minimum image by a search over the 27 neighbouring cells) gives the same 39 numbers: largest difference 0.000051 Å and 0.000018 in the signed volume, every count identical.

| Built system | Signed volume at N2 | Zn to OE1 / OE2 of Glu106 (Å) | `n_water_ring_place`: water oxygens within 2.5 Å of the six ring-carbon positions | Nearest water oxygen to them (Å) | Wider count: within 3.5 Å of the seven positions | Nearest (Å) | Water oxygens within 3.5 Å of the aniline N | State | Ligand RMSD (Å) |
|---|---|---|---|---|---|---|---|---|---|
| X | +0.7240 | 4.004 / 5.470 | 0 of 12247 | 3.025 | 1 | 3.025 | 0 | C | 0.000 |
| N | +0.7855 | 4.004 / 5.470 | 0 of 12250 | 3.025 | 1 | 3.025 | 0 | C | 0.610 |
| R | −0.7104 | 4.004 / 5.470 | 0 of 12232 | 4.148 | 1 | 2.845 | 1 | R | 3.573 |

**In every built system the place of the ring holds no water oxygen** (`n_water_ring_place` 0, 0 and 0). In X and N the ring itself is there, and the nearest water oxygen is the crystal water beside it (3.025 Å from the nearest ring-carbon position); in R the nearest one is 4.148 Å away, so the rule of the driver finds the place empty when the run begins. Within 3.5 Å of the same six positions the built systems hold 1, 1 and 0 water oxygens: the built X has one with the ring still in place, which is why the plan (8:56 AM, point 6) took 2.5 Å. The wider count is 1, 1 and 1; in R its one water oxygen lies 2.845 Å from the place of the aniline nitrogen and was placed by the solvation (the case that made the plan leave N1 out, 8:01 AM, point 7).

### 16.2 Synthetic trajectories with planted values

`make_test_trajectory.py` writes a system of 410 atoms in the conventions (288 in the solute: the 258 CA atoms of 6RL9, the named site atoms, zinc, ligand; 40 water molecules, 2 ions; five water molecules stand between protein and zinc in the topology, so that an index is not a position in the solute trajectory), with OpenMM's DCD and PDB writers. A second system has ND2 of Asn67 as well (411 atoms, 289 in the solute) and the key `asn67_ND2`; the first has neither. 23 runs. The six runs of version 1.0 are kept as they were (same frames, same seed):

| Run | Frames (solute / full) | What is planted |
|---|---|---|
| X1 | 3000 / 300, three segments | crystal pose with atom labels swapped between equivalent atoms; a 0.5 ns excursion to the rotated pose (too short to count); rotated pose from 9.01 ns (the crystal ligand moved rigidly onto the rotated start), with 20 intermediate frames inside (ligand translated by 3.000 Å); back for 2 ns; rotated again; exactly 90 % of the last 20 ns rotated |
| X2 | 3000 / 300, one segment | crystal pose throughout, labels swapped in three of four frames; plain manifest |
| X3 | 3000 / 300, two segments | as X1 with one rotated frame fewer in the last 20 ns (1799 of 2000); plain manifest; DCD headers as OpenMM's own reporter writes them |
| N1 | 1234 / 123, two segments | near-native start throughout; planned 100 ns, so unfinished; half a frame of rubbish at the end of the last file |
| R1 | 3000 / 300, two segments | rotated start for 8 ns, then the crystal pose; the system folder named in the manifest is not on this machine |
| X9 | 250 / 20 (outside the campaign) | the second segment repeats the last 50 frames of the first (a restart from an old checkpoint); the manifest carries the checksum of another `selections.json` |

In each of these: the sulfonamide hydrogen is turned away from Thr199 in every fifth frame; the zinc is moved 0.300 Å away from N2 in every tenth frame; the CA atoms are scaled about their centroid by 1 ± 0.01 (known CA RMSD and RMSF, fit unchanged); every frame is rotated at random and shifted; water is wrapped into a truncated-octahedral cell whose edge changes with time; the ligand (every 7th frame) and the zinc (every 11th) are moved by lattice vectors; a water sits 2.400 Å from the zinc in every tenth full frame and 3.690 Å otherwise; 0 to 3 water molecules sit 3.000 Å from the aniline N.

Added with version 1.1: nine validation runs (folders named `X0`, analysed without any option) and one production run. They keep the hydrogen, zinc, CA, rotation and lattice schedules above; the water at 2.400 Å sits in chosen full frames only. A pose is bent by a rigid turn: ND1 of His119 about the zinc (one N–Zn–N angle), the ring about an axis through S (C4–S–N2) or about the C4–S bond (τ_fold); in chosen frames the zinc is put 2.500 Å from the nearer sulfonyl oxygen. The deviations and fractions in this table are the planted ones, computed by the generator:

| Run | Frames (solute / full) | What is planted | Gate G-V by design |
|---|---|---|---|
| test/X0 | 500 / 50 | the crystal geometry, in the folder of the campaign beside X1 to R1; manifest of the driver with a target of 5 ns; a water at 2.400 Å in 2 of 50 (4.0 %) of the full frames; system without `asn67_ND2` | **passes**; the campaign must leave the run out |
| edge/X0 | 500 / 50, two segments | every new criterion just inside its limit: His96 NE2–Zn–His119 ND1 +9.15° from 6RL9, C4–S–N2 +2.49°, τ_fold +18.96°; a sulfonyl oxygen within 2.6 Å in 25 of 500 (5.0 %) of the solute frames, one of them a full frame, so that with one water frame 2 of 50 (4.0 %) of the full frames have an oxygen near; the Asn67 amide turned in 20 frames; manifest with model and start only | **passes** |
| angle/X0 | 500 / 50 | His96 NE2–Zn–His119 ND1 +11.85° from 6RL9, the other five within 7.08° | fails by GV3 c alone |
| sulfonyl/X0 | 500 / 50 | a sulfonyl oxygen within 2.6 Å in 26 of 500 (5.2 %) of the solute frames, none of them a full frame; full fraction 2 of 50 (4.0 %); plain manifest with replica 0 and 5 ns | fails by the solute fraction of GV3 d alone |
| water/X0 | 500 / 50, two segments | full frames with an oxygen near: 3 of 50 (6.0 %), two by water and one by a sulfonyl oxygen; solute fraction 5 of 500 (1.0 %) | fails by the full fraction of GV3 d alone |
| csn/X0 | 500 / 50 | C4–S–N2 +3.49° from 108.8°; manifest with model and start only | fails by GV4 a alone |
| tau/X0 | 500 / 50 | τ_fold +20.96° from 132.1° | fails by GV4 b alone |
| long/X0 | 600 / 60, two segments | crystal geometry for 6 ns, the manifest of the driver with a target of 6 ns; water near in 3 of 60 (5.0 %) of the full frames, the limit | **passes**; planned length stays 5 ns; two warnings |
| short/X0 | 300 / 30 | crystal geometry for 3 ns; no manifest at all | **passes**, provisional (unfinished) |
| gate/R1 | 400 / 40 | a production run of the system with ND2: rotated start for 2.5 ns, then the crystal pose; the amide turned in every fourth frame | no verdict; numbers only |

Added with version 1.2: five runs with the new quantities planted (folders `runs_v12/`), four of them with equilibration files. In these runs the hydrogen on N2 stays in place unless a frame is chosen for one of two moves: the hydrogen is mirrored through the plane Zn–N2–S (the signed volume keeps its size and changes its sign), or OG1 of Thr199 is moved 2.000 Å away from N2 (no hydrogen bond, the signed volume unchanged). OE1 or OE2 of Glu106 is put at a chosen distance from the zinc; three water molecules are planted at the crystal ring: A at the centre of the ring, B 3.400 Å from the crystal N1, C 3.600 Å from the crystal C3. The old runs carry the new quantities as well: there the hydrogen that is turned away in every fifth frame has the other sign at N2, so no stretch of one sign is longer than four frames and no sign change lasts.

| Run | Frames (solute / full) | Equilibration files | What is planted |
|---|---|---|---|
| conv/R1 | 2200 / 220, two segments | six stages of the protocol, 200 / 20 frames; headers and manifest as the driver writes them | rotated start for 6 ns, then the crystal pose: first passage at 6.01 ns, verdict MIXED (80 % of the last 20 ns in the opposite state). Hydrogen mirrored for 99 frames (too short to count) and for exactly 100 frames: lasting sign changes at 4.01, 5.01, 6.01 ns; 1799 frames positive, 401 negative, in the last 20 ns 1799 and 201. Hydrogen bond in 301 of 600 frames before the first passage and in 800 of 1600 after it. OE1 of Glu106 at 2.500 Å in 44 frames (20 of them full frames) and at 2.700 Å in 22, OE2 at 2.550 Å in 22. Water at the ring: none, then B, then A and B (20 full frames each), then the water at the aniline N of the returned ligand; mean of the wider count 1.3636. In the equilibration the wider count rises from 0 to 1 to 2; it is 1 at the end of stage 3. (Version 1.3: the count of the driver sees water A alone, in 20 full frames of the production, mean 0.0909; it is 0 at the end of stage 3.) |
| flip/X1 | 300 / 30 | six stages of the protocol, no manifest; headers as OpenMM's own reporter writes them | the pose turns during the equilibration (first frame in the other state and with the other sign: frame 120, stage 5, 0 ns). The production begins with five intermediate frames with the other sign at N2 (first frame differs from the start in state and in sign), then the rotated pose, then the crystal pose: first passage at 0.06 ns, lasting sign changes at 0.01, 1.51 ns; hydrogen bond in 0 of 5 frames before the first passage, 245 of 295 after it; unfinished (100 ns planned by default) |
| scaled/N1 | 200 / 20 | a shortened equilibration as the driver records it: stages of 5, 15, 15, 20, 20 and 25 ps, so the files hold 0, 1, 1, 2, 2, 2 solute frames and no full frame | near-native start throughout; the times of the eight equilibration frames are -0.085, -0.07, -0.055, -0.045, -0.035, -0.025, -0.015, -0.005 ns: the last frame lies 5 ps before the production |
| ext/R1 | 200 / 20 | an equilibration whose held part was extended twice, as the driver of 6 Oct 2026 records it: the steps eq01, eq02, eq03, eq03x1, eq03x2, eq04, eq05, eq06, 3000 ps, 300 / 30 frames; the manifest carries the driver's count at the end of eq03, eq03x1 and eq03x2 | rotated start throughout (verdict STAYED). Water at the crystal ring positions (seven positions, 3.5 Å: C, then B and C) at the end of eq03, eq03x1 and eq03x2: 0, 0, 1, as the driver's records in the manifest say; solute frame 70 is the first of eq03x1 (-2.29 ns), frame 170 the first of eq04 (-1.29 ns) |
| glu/X0 | 500 / 50 | none | a validation run in the crystal geometry with OE1 of Glu106 at 2.500 Å in 30 of 500 solute frames (6 %) and in no full frame; hydrogen bond in exactly 80 % of the frames; gate G-V **passes**, because Glu106 is reported and is no criterion |

Added with version 1.3: the count that the driver acts on is planted in every run beside the wider count, by second formulas in the generator (crystal positions of the six ring carbons, 2.5 Å; for a run whose manifest names other settings, those). Three more water molecules can be planted at the crystal ring: D 2.400 Å from the crystal C5 along the normal of the ring (inside the place of the ring), E 2.600 Å from the crystal C2 on the other side (outside it, inside the wider count), F 2.300 Å from the crystal N1, away from the ring (near the aniline nitrogen only, so in the wider count alone). Water A, at the centre of the ring, is within 2.5 Å of all six positions and must count once. The 21 runs of version 1.2 keep their frames (section 16.4); one of them gets another manifest, and two runs are new (folder `runs_v13/`):

| Run | Frames (solute / full) | Equilibration files | What is planted |
|---|---|---|---|
| ext/R1 (of version 1.2) | 200 / 20 | as before | the same frames; the manifest now carries `hydration_check.settings` as the driver wrote them before 8:01 AM of 6 Oct 2026: 7 site atoms (with N1) and 3.5 Å. The analysis must take these, so that `n_water_ring_place` equals the wider count in this run, agree with the driver's counts 0, 0, 1 and give 1 warning (settings other than those of the plan) |
| place/R1 | 200 / 20 | an equilibration whose held part was extended three times: the steps eq01, eq02, eq03, eq03x1, eq03x2, eq03x3, eq04, eq05, eq06, 3500 ps, 350 / 35 frames; the manifest with the settings and the step records of the hydration check as the driver writes them now (six ring carbons, 2.5 Å) | rotated start throughout. Production, five full frames each: E alone, D, D with E and F, A with D: `n_water_ring_place` 0, 1, 1, 2 with the nearest water oxygen at 2.600, 2.400, 2.400, 1.382 Å; wider count 1, 1, 3, 2. Equilibration: E alone up to the last full frame but one of eq03x3, then D and F; the driver's records at the end of eq03, eq03x1, eq03x2, eq03x3 say 0, 0, 0, 1 (nearest 2.6, 2.6, 2.6, 2.4 Å), and the wider count there is 1, 1, 1, 2. Water C, 3.600 Å from the crystal C3, is there throughout and belongs to neither count |
| default/R1 | 200 / 20 | none | the production frames of place/R1 with a manifest of the plain keys, which says nothing of a hydration check: the analysis must use the defaults (`source` `specification (defaults)`) and return the same counts |

A planted value comes from one of three sources: **construction** (exact by design, typed into the generator or counted from its schedule), **starts.csv** (the numbers of the system builder, written with its own code), or **second formulas** (other formulas for distance, angle, torsion, RMSD and signed volume inside the generator, applied in double precision before the frame is rotated, wrapped and stored; the generator never imports `analyse_md.py`). Gate G-V is evaluated a second time inside the generator, from its own values of every frame and with the limits typed from section 11, and the verdicts typed by hand are checked against that when the runs are written. For the quantities of versions 1.2 and 1.3 the generator has its own code for the stretches of one sign, the lasting sign changes and the first passage; the design of the seven runs of these versions is typed by hand and checked against the frames as built before a file is written.

**Result: 7339 of 7339 planted values are returned within tolerance, 0 fail** (2098 by construction, 270 against `starts.csv`, 4971 by second formulas). Of these, 6517 belong to the 21 runs and the campaign of version 1.2 (6232 then; the others are the count of the driver, its nearest distance and its settings on those runs), 822 to the two new runs. This includes every column of every frame of 22 runs and of the equilibration series of 5 (1822 column comparisons): largest difference 0.00006 Å, 0.0013° and 0.000041 in the signed volume, and no state, flag, sign, labelling or count differs in any frame. Full list: `test_report.md` and `.json`; planted values per frame: `test/planted/` of the scratch folder. A selection:

| Run | Planted quantity | Planted | Returned | Source |
|---|---|---|---|---|
| X2 | ligand RMSD, mean (labels swapped in 3 of 4 frames) (Å) | 0 | 1e-06 | construction |
| X2 | ligand RMSD, largest value (Å) | 0 | 1.1e-05 | construction |
| X2 | ligand RMSD without the symmetry correction, largest value (must be large) | yes | yes | construction |
| X1 | frames 1602 and 1603 (labels swapped in one): same ligand RMSD (Å) | 0 | 0 | construction |
| X1 | frames 1602 and 1603: labelling differs | yes | yes | construction |
| X1 | frame 960 (ligand translated by 3.000 Å): ligand RMSD (Å) | 3 | 3 | construction |
| X1 | frame 960: d1 (starts.csv 3.897 + 3.000) (Å) | 6.897 | 6.8966 | starts.csv |
| X1 | frame 960: state | I | I | construction |
| X1 | frame 960: assigned state (carried from the visit it lies in) | R | R | construction |
| R1 | frame 1, lig_rmsd_A against starts.csv (R, rmsd_to_crystal_A) (Å) | 3.573 | 3.5732 | starts.csv |
| R1 | frame 1, d1_N1_F131CZ_A against starts.csv (R, N1_Phe131CZ) (Å) | 9.585 | 9.5849 | starts.csv |
| R1 | frame 1, d2_N1_N67OD1_A against starts.csv (R, N1_Asn67OD1) (Å) | 2.597 | 2.5968 | starts.csv |
| R1 | frame 1, tau_fold_deg against starts.csv (R, tau_fold) (°) | 142.744 | 142.744 | starts.csv |
| N1 | ligand RMSD, mean (starts.csv, near-native pose) (Å) | 0.61 | 0.610426 | starts.csv |
| X1 | lig_rmsd_A: mean (Å) | 2.31128 | 2.31128 | second formulas |
| X1 | lig_rmsd_A: mean_last_half (Å) | 3.32413 | 3.32413 | second formulas |
| X1 | rule 1: left the crystal pose | yes | yes | construction |
| X1 | verdict | CONVERTED | CONVERTED | construction |
| X1 | fraction of the last 20 ns in the opposite state (boundary of the 90 % rule) | 0.9 | 0.9 | construction |
| X1 | count_R_window | 1800 | 1800 | construction |
| X1 | first passage into the opposite state (ns) | 9.01 | 9.01 | construction |
| X1 | transitions | 3 | 3 | construction |
| X1 | long visits (at least 1 ns) | 5 | 5 | construction |
| X1 | count_C | 1050 | 1050 | construction |
| X1 | count_R | 1930 | 1930 | construction |
| X1 | count_I | 20 | 20 | construction |
| X3 | verdict | MIXED | MIXED | construction |
| X3 | count_R_window | 1799 | 1799 | construction |
| X2 | verdict | STAYED | STAYED | construction |
| X2 | transitions | 0 | 0 | construction |
| X2 | first passage into the opposite state (ns) | none | none | construction |
| R1 | verdict | CONVERTED | CONVERTED | construction |
| R1 | first passage into the opposite state (ns) | 8.01 | 8.01 | construction |
| R1 | transitions | 1 | 1 | construction |
| N1 | verdict | STAYED | STAYED | construction |
| N1 | verdict provisional (run not finished) | yes | yes | construction |
| N1 | solute frames seen | 1234 | 1234 | construction |
| N1 | time seen (ns) | 12.34 | 12.34 | construction |
| N1 | verdict window complete (20 ns) | no | no | construction |
| X2 | hydrogen bond N2–H···OG1, fraction of frames | 0.8 | 0.8 | construction |
| X2 | hydrogen bond N2–H···OG1, fraction in the last half | 0.8 | 0.8 | construction |
| X2 | hydrogen bond N2–H···OG1, the three block fractions | 0.8, 0.8, 0.8 | 0.8, 0.8, 0.8 | construction |
| X2 | Zn–N2, mean (starts.csv 1.931 + 0.300 in every tenth frame) (Å) | 1.961 | 1.96144 | starts.csv |
| X2 | Zn–N2, largest value (starts.csv 1.931 + 0.300) (Å) | 2.231 | 2.23145 | starts.csv |
| X2 | nearest water to the zinc, mean (0.1 × 2.400 + 0.9 × 3.690) (Å) | 3.561 | 3.561 | construction |
| X2 | nearest water to the zinc, smallest (Å) | 2.4 | 2.39999 | construction |
| X2 | fraction of full frames with a fifth ligand (water) | 0.1 | 0.1 | construction |
| X2 | water molecules within 3.5 Å of the aniline N, mean | 1.5 | 1.5 | construction |
| X2 | tau_fold, mean (starts.csv, unchanged by the label swaps) (°) | 132.056 | 132.056 | starts.csv |
| X2 | tau_fold, standard deviation (°) | 0 | 6.9e-05 | construction |
| X1 | tau_fold, mean (rigid moves and label swaps keep the crystal value of starts.csv) (°) | 132.056 | 132.056 | starts.csv |
| X1 | S–N2, standard deviation (rigid moves) (Å) | 0 | 1e-06 | construction |
| X2 | CA RMSD, mean (0.01 × radius of gyration of the CA atoms) (Å) | 0.171798 | 0.171798 | construction |
| X2 | CA RMSF of 258 residues, largest difference (Å) | 0 | 5e-05 | construction |
| X2 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | construction |
| X1 | frames in which the ligand was in another periodic image | 429 | 429 | construction |
| X1 | frames in which the zinc was in another periodic image | 273 | 273 | construction |
| X1 | DCD headers: segments contiguous (none = the headers tell nothing) | yes | yes | construction |
| X3 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | construction |
| X1 | manifest and files agree (none = no list) | yes | yes | construction |
| R1 | where the system folder was found | convention (the system_dir of the manifest is not on this machine) | convention (the system_dir of the manifest is not on this machine) | construction |
| X1 | selections.json is the file of the manifest (none = no checksum there) | yes | yes | construction |
| X2 | gate G-V, criterion 1 (mean ligand RMSD at most 1.0 Å) | yes | yes | construction |
| X2 | gate G-V, criterion 2 (hydrogen bond in at least 80 % of frames; here exactly 80 %) | yes | yes | construction |
| X2 | gate G-V, criterion 3, Zn–N2 window | yes | yes | construction |
| X2 | gate G-V, criterion 3, angle His96 NE2–Zn–His119 ND1 within 10° of 98.4° (version 1.0: in 100–120°, not met) | yes | yes | construction |
| X2 | gate G-V, criterion 3, fifth ligand in the solute trajectory alone (no sulfonyl oxygen near) | yes | yes | construction |
| X2 | gate G-V, criterion 3, full frames with an oxygen within 2.6 Å, of 300 | 30 | 30 | construction |
| X2 | gate G-V, criterion 3, no fifth ligand (planted: water at 2.400 Å in 10 % of the full frames) | no | no | construction |
| X2 | gate G-V, criterion 4 (C4–S–N2 within 3° of 108.8°, τ_fold within 20° of 132.1°) | yes | yes | construction |
| X2 | gate G-V passed | no | no | construction |
| X2 | gate G-V verdict (a production run has none) | none | none | construction |
| X2 | role of the run (X0 = validation) | production | production | construction |
| X2 | gate G-V applies to this run (only to X0) | no | no | construction |
| X2 | rule 1 applies to this run (production runs of start X) | yes | yes | construction |
| X9 (overlap) | solute frames seen | 250 | 250 | construction |
| X9 (overlap) | repeated frames found at the segment boundary | yes | yes | construction |
| X9 (overlap) | frame of the first segment that is repeated | 150 | 150 | construction |
| X9 (overlap) | DCD headers: segments contiguous (none = the headers tell nothing) | no | no | construction |
| X9 (overlap) | manifest and files agree (none = no list) | no | no | construction |
| X9 (overlap) | selections.json is the file of the manifest (none = no checksum there) | no | no | construction |
| X9 (overlap) | warnings written (headers, repeated frame, manifest, system file) | 4 | 4 | construction |
| campaign | rule1.n_X_left_crystal | 2 | 2 | construction |
| campaign | rule1.failed | yes | yes | construction |
| campaign | n_converted.X | 1 | 1 | construction |
| campaign | n_converted.R | 1 | 1 | construction |
| campaign | n_stayed.N | 1 | 1 | construction |
| campaign | n_mixed.X | 1 | 1 | construction |
| campaign | outcome | validation_failed | validation_failed | construction |
| campaign | final | no | no | construction |
| campaign | order of the runs | test/X1 test/X2 test/X3 test/N1 test/R1 | test/X1 test/X2 test/X3 test/N1 test/R1 | construction |
| campaign | n_runs | 5 | 5 | construction |
| campaign | rule1.n_X_runs | 3 | 3 | construction |
| campaign | folders left out of the campaign (the validation run X0) | test/X0 | test/X0 | construction |
| campaign | lines of campaign_table.csv | 5 | 5 | construction |
| test/X0 | role of the run (X0 = validation) | validation | validation | construction |
| test/X0 | planned length (X0: 5 ns by the specification) (ns) | 5 | 5 | construction |
| test/X0 | run complete | yes | yes | construction |
| test/X0 | gate G-V applies to this run (only to X0) | yes | yes | construction |
| test/X0 | rule 1 applies to this run (production runs of start X) | no | no | construction |
| test/X0 | gate G-V by design: verdict | PASSED | PASSED | construction |
| test/X0 | gate G-V by design: the criteria that are not met | none | none | construction |
| test/X0 | by design: full frames with any oxygen within 2.6 Å | 2 | 2 | construction |
| test/X0 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | no | no | construction |
| test/X0 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (°) | 0.807205 | 0.807202 | second formulas |
| edge/X0 | gate G-V by design: verdict | PASSED | PASSED | construction |
| edge/X0 | gate G-V by design: the criteria that are not met | none | none | construction |
| edge/X0 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (°) | 9.15181 | 9.15181 | second formulas |
| edge/X0 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (°) | 2.4946 | 2.4946 | second formulas |
| edge/X0 | gate G-V: tau_fold_deg, mean minus the reference (°) | 18.9561 | 18.9561 | second formulas |
| edge/X0 | by design: solute frames with a sulfonyl oxygen within 2.6 Å | 25 | 25 | construction |
| edge/X0 | gate G-V: solute frames counted for that fraction | 500 | 500 | construction |
| edge/X0 | by design: full frames with any oxygen within 2.6 Å | 2 | 2 | construction |
| edge/X0 | gate G-V: full frames counted for that fraction | 50 | 50 | construction |
| edge/X0 | where the system folder was found | convention | convention | construction |
| edge/X0 | by design: ND2 nearer in state C, frames of state C, ND2 nearer in state R, frames of state R | 480, 500, 0, 0 | 480, 500, 0, 0 | construction |
| angle/X0 | gate G-V by design: verdict | FAILED | FAILED | construction |
| angle/X0 | gate G-V by design: the criteria that are not met | GV3 ang_H96NE2_Zn_H119ND1_deg | GV3 ang_H96NE2_Zn_H119ND1_deg | construction |
| angle/X0 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (°) | 11.8451 | 11.8451 | second formulas |
| sulfonyl/X0 | gate G-V by design: verdict | FAILED | FAILED | construction |
| sulfonyl/X0 | gate G-V by design: the criteria that are not met | GV3 fifth ligand, solute trajectory | GV3 fifth ligand, solute trajectory | construction |
| sulfonyl/X0 | by design: solute frames with a sulfonyl oxygen within 2.6 Å | 26 | 26 | construction |
| sulfonyl/X0 | by design: full frames with any oxygen within 2.6 Å | 2 | 2 | construction |
| water/X0 | gate G-V by design: verdict | FAILED | FAILED | construction |
| water/X0 | gate G-V by design: the criteria that are not met | GV3 fifth ligand, full trajectory | GV3 fifth ligand, full trajectory | construction |
| water/X0 | by design: solute frames with a sulfonyl oxygen within 2.6 Å | 5 | 5 | construction |
| water/X0 | by design: full frames with any oxygen within 2.6 Å | 3 | 3 | construction |
| csn/X0 | gate G-V by design: verdict | FAILED | FAILED | construction |
| csn/X0 | gate G-V by design: the criteria that are not met | GV4 ang_C4_S_N2_deg | GV4 ang_C4_S_N2_deg | construction |
| csn/X0 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (°) | 3.4946 | 3.4946 | second formulas |
| tau/X0 | gate G-V by design: verdict | FAILED | FAILED | construction |
| tau/X0 | gate G-V by design: the criteria that are not met | GV4 tau_fold_deg | GV4 tau_fold_deg | construction |
| tau/X0 | gate G-V: tau_fold_deg, mean minus the reference (°) | 20.9561 | 20.9561 | second formulas |
| long/X0 | solute frames seen | 600 | 600 | construction |
| long/X0 | planned length (X0: 5 ns by the specification) (ns) | 5 | 5 | construction |
| long/X0 | warnings | 2 | 2 | construction |
| long/X0 | gate G-V by design: verdict | PASSED | PASSED | construction |
| long/X0 | by design: full frames with any oxygen within 2.6 Å | 3 | 3 | construction |
| long/X0 | gate G-V: full frames counted for that fraction | 60 | 60 | construction |
| short/X0 | role of the run (X0 = validation) | validation | validation | construction |
| short/X0 | planned length (X0: 5 ns by the specification) (ns) | 5 | 5 | construction |
| short/X0 | run complete | no | no | construction |
| short/X0 | gate G-V by design: verdict | PASSED | PASSED | construction |
| short/X0 | gate G-V by design: provisional | yes | yes | construction |
| short/X0 | where the system folder was found | convention | convention | construction |
| gate/R1 | role of the run (X0 = validation) | production | production | construction |
| gate/R1 | gate G-V: verdict (none for a run that is not X0) | none | none | second formulas |
| gate/R1 | gate G-V: the criteria that are not met | GV1 ligand RMSD, GV3 Zn-N2 | GV1 ligand RMSD, GV3 Zn-N2 | second formulas |
| gate/R1 | verdict | MIXED | MIXED | construction |
| gate/R1 | first passage into the opposite state (ns) | 2.51 | 2.51 | construction |
| gate/R1 | frame 2, d_N1_N67ND2_A against starts.csv (R, N1_Asn67ND2) (Å) | 3.31 | 3.3102 | starts.csv |
| gate/R1 | frames with ND2 nearer to the aniline N than OD1 | 176 | 176 | second formulas |
| gate/R1 | by design: ND2 nearer in state C, frames of state C, ND2 nearer in state R, frames of state R | 113, 150, 63, 250 | 113, 150, 63, 250 | construction |
| X2 | N2: frames with positive, negative and zero signed volume | 2400, 600, 0 | 2400, 600, 0 | second formulas |
| X2 | N2: stretches of one sign, and those of at least 1 ns | 1200, 0 | 1200, 0 | second formulas |
| X2 | N2: number of sign changes that last at least 1 ns | 0 | 0 | second formulas |
| X2 | n2_volume: mean | 0.451372 | 0.451372 | second formulas |
| X1 | N2: positive, negative and zero in the last 20 ns | 520, 1480, 0 | 520, 1480, 0 | second formulas |
| X1 | N2, frames of state I: number, positive, negative, zero | 20, 16, 4, 0 | 20, 16, 4, 0 | second formulas |
| X1 | hydrogen bond N2–H···OG1 before the first passage: frames, frames with the bond | 900, 720 | 900, 720 | second formulas |
| X1 | hydrogen bond N2–H···OG1 from the first passage on: frames, frames with the bond | 2100, 1664 | 2100, 1664 | second formulas |
| X2 | hydrogen bond around the first passage: the six numbers (none without a first passage) | none, none, none, none, none, none | none, none, none, none, none, none | construction |
| X2 | Glu106: smallest distance zinc to OE1 (Å) | 4.00373 | 4.00372 | second formulas |
| X2 | water at the crystal ring positions: mean number within 3.5 Å | 1.5 | 1.5 | second formulas |
| R1 | first production frame: sign at N2, and whether it differs from the sign of the start | -1, no | -1, no | second formulas |
| X2 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | construction |
| conv/R1 | verdict | MIXED | MIXED | construction |
| conv/R1 | first passage into the opposite state (ns) | 6.01 | 6.01 | construction |
| conv/R1 | by design: frames with positive and with negative signed volume at N2 | 1799, 401 | 1799, 401 | construction |
| conv/R1 | by design: positive and negative in the last 20 ns | 1799, 201 | 1799, 201 | construction |
| conv/R1 | by design, frames of state C: number, positive, negative | 1600, 1600, 0 | 1600, 1600, 0 | construction |
| conv/R1 | by design, frames of state R: number, positive, negative | 600, 199, 401 | 600, 199, 401 | construction |
| conv/R1 | by design: stretches of one sign, and those of at least 1 ns | 6, 5 | 6, 5 | construction |
| conv/R1 | by design: times of the sign changes that last at least 1 ns (ns) | 4.01, 5.01, 6.01 | 4.01, 5.01, 6.01 | construction |
| conv/R1 | by design: sign after each of these changes | 1, -1, 1 | 1, -1, 1 | construction |
| conv/R1 | by design: hydrogen bond before the first passage (frames, with the bond) and after it (frames, with the bond) | 600, 301, 1600, 800 | 600, 301, 1600, 800 | construction |
| conv/R1 | hydrogen bond N2–H···OG1, fraction before the first passage | 0.501667 | 0.501667 | second formulas |
| conv/R1 | hydrogen bond N2–H···OG1, fraction from the first passage on | 0.5 | 0.5 | second formulas |
| conv/R1 | by design: solute frames with OE1, with OE2, with either oxygen of Glu106 within 2.6 Å | 44, 22, 66 | 44, 22, 66 | construction |
| conv/R1 | by design: smallest distance zinc to OE1 of Glu106 (Å) | 2.5 | 2.5 | construction |
| conv/R1 | by design: smallest distance zinc to OE2 of Glu106 (Å) | 2.55 | 2.55 | construction |
| conv/R1 | Glu106: full frames with a protein oxygen within 2.6 Å of the zinc | 20 | 20 | second formulas |
| conv/R1 | by design: water oxygens at the crystal ring positions, full frames that differ | 0 | 0 | construction |
| conv/R1 | by design: mean number of water oxygens at the crystal ring positions | 1.36364 | 1.36364 | construction |
| conv/R1 | by design: smallest distance of a water oxygen to the crystal ring positions (Å) | 1.382 | 1.38236 | construction |
| conv/R1 | by design: state of the first production frame, differs from the start; sign at N2, differs from the start | R, no, -1, no | R, no, -1, no | construction |
| conv/R1 | by design: equilibration, solute frames, full frames, length (ps) | 200, 20, 2000 | 200, 20, 2000 | construction |
| conv/R1 | by design: equilibration, where the stage table came from | manifest equilibration.stages | manifest equilibration.stages | construction |
| conv/R1 | by design: equilibration, solute and full frames of each stage | 10, 1; 30, 3; 30, 3; 40, 4; 40, 4; 50, 5 | 10, 1; 30, 3; 30, 3; 40, 4; 40, 4; 50, 5 | construction |
| conv/R1 | by design: equilibration, time of the first and of the last solute frame (ns) | -1.99, 0 | -1.99, 0 | construction |
| conv/R1 | by design: equilibration, time of solute frame 10 (the first of stage 2) (ns) | -1.89, 2 | -1.89, 2 | construction |
| conv/R1 | by design: equilibration, water oxygens at the crystal ring positions in every full frame | 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 | 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 | construction |
| conv/R1 | by design: equilibration, water oxygens at the crystal ring positions at the end of stage 3 (the held part) | 1 | 1 | construction |
| conv/R1 | equilibration: frames in which the ligand and the zinc were in another periodic image | 29, 18 | 29, 18 | construction |
| flip/X1 | verdict | MIXED | MIXED | construction |
| flip/X1 | first passage into the opposite state (ns) | 0.06 | 0.06 | construction |
| flip/X1 | by design: state of the first production frame, differs from the start; sign at N2, differs from the start | I, yes, -1, yes | I, yes, -1, yes | construction |
| flip/X1 | first production frame: in the opposite state | no | no | second formulas |
| flip/X1 | by design: times of the sign changes that last at least 1 ns (ns) | 0.01, 1.51 | 0.01, 1.51 | construction |
| flip/X1 | by design: hydrogen bond before the first passage (frames, with the bond) and after it (frames, with the bond) | 5, 0, 295, 245 | 5, 0, 295, 245 | construction |
| flip/X1 | by design, frames of state I: number, positive, negative | 5, 0, 5 | 5, 0, 5 | construction |
| flip/X1 | by design: equilibration, where the stage table came from | protocol (constants of the specification) | protocol (constants of the specification) | construction |
| flip/X1 | by design: equilibration, first_frame_not_in_state_of_start (frame, stage, extension step, time) | 120, 5, 0, -0.79 | 120, 5, 0, -0.79 | construction |
| flip/X1 | by design: equilibration, first_frame_n2_sign_not_of_start (frame, stage, extension step, time) | 120, 5, 0, -0.79 | 120, 5, 0, -0.79 | construction |
| flip/X1 | by design: equilibration, state and sign at N2 of the last frame | R, -1 | R, -1 | construction |
| flip/X1 | equilibration, stage 5: frames of state C, R, I; frames with positive and negative sign at N2 | 10, 30, 0, 10, 30 | 10, 30, 0, 10, 30 | second formulas |
| scaled/N1 | by design: equilibration, solute frames, full frames, length (ps) | 8, 0, 100 | 8, 0, 100 | construction |
| scaled/N1 | by design: equilibration, solute and full frames of each stage | 0, 0; 1, 0; 1, 0; 2, 0; 2, 0; 2, 0 | 0, 0; 1, 0; 1, 0; 2, 0; 2, 0; 2, 0 | construction |
| scaled/N1 | by design: equilibration, the times of all solute frames (ns) | -0.085, -0.07, -0.055, -0.045, -0.035, -0.025, -0.015, -0.005 | -0.085, -0.07, -0.055, -0.045, -0.035, -0.025, -0.015, -0.005 | construction |
| scaled/N1 | equilibration: lines of timeseries_equil_full.csv | 0 | 0 | construction |
| scaled/N1 | warnings | 0 | 0 | construction |
| ext/R1 | verdict | STAYED | STAYED | construction |
| ext/R1 | warnings | 1 | 1 | construction |
| ext/R1 | by design: equilibration, the steps in the order of the run | eq01, eq02, eq03, eq03x1, eq03x2, eq04, eq05, eq06 | eq01, eq02, eq03, eq03x1, eq03x2, eq04, eq05, eq06 | construction |
| ext/R1 | by design: equilibration, solute frames, full frames, length (ps) | 300, 30, 3000 | 300, 30, 3000 | construction |
| ext/R1 | by design: equilibration, solute and full frames of each stage | 10, 1; 30, 3; 30, 3; 50, 5; 50, 5; 40, 4; 40, 4; 50, 5 | 10, 1; 30, 3; 30, 3; 50, 5; 50, 5; 40, 4; 40, 4; 50, 5 | construction |
| ext/R1 | by design: equilibration, time of the first and of the last solute frame (ns) | -2.99, 0 | -2.99, 0 | construction |
| ext/R1 | by design: equilibration, stage, extension step and time of solute frame 70 (the first of eq03x1) | 3, 1, -2.29 | 3, 1, -2.29 | construction |
| ext/R1 | by design: equilibration, stage, extension step and time of solute frame 170 (the first of eq04) | 4, 0, -1.29 | 4, 0, -1.29 | construction |
| ext/R1 | by design: site atoms and cut-off of the count that the driver acts on (Å) | C1, C2, C3, C4, C5, C6, N1, 3.5 | C1, C2, C3, C4, C5, C6, N1, 3.5 | construction |
| ext/R1 | by design: where these settings came from, and whether they are the defaults of the specification | manifest hydration_check.settings, no | manifest hydration_check.settings, no | construction |
| ext/R1 | by design: mean number of water oxygens in the place of the ring | 1 | 1 | construction |
| ext/R1 | by design: equilibration, water oxygens in the place of the ring (the count the driver acts on) at the end of eq03, eq03x1, eq03x2 | 0, 0, 1 | 0, 0, 1 | construction |
| ext/R1 | by design: equilibration, the count of the driver at the end of these steps agrees with the analysis | yes, yes, yes | yes, yes, yes | construction |
| ext/R1 | by design: equilibration, the nearest water oxygen of the driver at the end of these steps agrees with the analysis | yes, yes, yes | yes, yes, yes | construction |
| glu/X0 | gate G-V by design: verdict | PASSED | PASSED | construction |
| glu/X0 | gate G-V by design: the criteria that are not met | none | none | construction |
| glu/X0 | by design: solute frames with OE1, with OE2, with either oxygen of Glu106 within 2.6 Å | 30, 0, 30 | 30, 0, 30 | construction |
| glu/X0 | by design: solute frames with a sulfonyl oxygen within 2.6 Å | 0 | 0 | construction |
| glu/X0 | by design: full frames with any oxygen within 2.6 Å | 2 | 2 | construction |
| glu/X0 | Glu106 beside gate G-V, OE1: frames within 2.6 Å, frames counted | 30, 500 | 30, 500 | second formulas |
| glu/X0 | Glu106 beside gate G-V, OE1: smallest distance (Å) | 2.5 | 2.5 | second formulas |
| X2 | water in the place of the ring: site atoms and cut-off used (Å) | C1, C2, C3, C4, C5, C6, 2.5 | C1, C2, C3, C4, C5, C6, 2.5 | construction |
| X2 | water in the place of the ring: where the settings came from, and whether they are the defaults | specification (defaults), yes | specification (defaults), yes | construction |
| X2 | water in the place of the ring: mean number | 0 | 0 | second formulas |
| X2 | water in the place of the ring: smallest distance of a water oxygen (Å) | 4.02013 | 4.02012 | second formulas |
| conv/R1 | by design: mean number of water oxygens in the place of the ring | 0.090909 | 0.090909 | construction |
| conv/R1 | by design: smallest distance of a water oxygen to the place of the ring (Å) | 1.382 | 1.38236 | construction |
| conv/R1 | by design: equilibration, water oxygens in the place of the ring in every full frame | 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 | 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 | construction |
| conv/R1 | by design: equilibration, water oxygens in the place of the ring at the end of stage 3 | 0 | 0 | construction |
| place/R1 | verdict | STAYED | STAYED | construction |
| place/R1 | warnings | 0 | 0 | construction |
| place/R1 | by design: site atoms and cut-off of the count that the driver acts on (Å) | C1, C2, C3, C4, C5, C6, 2.5 | C1, C2, C3, C4, C5, C6, 2.5 | construction |
| place/R1 | by design: where these settings came from, and whether they are the defaults of the specification | manifest hydration_check.settings, yes | manifest hydration_check.settings, yes | construction |
| place/R1 | by design: water oxygens in the place of the ring in every full frame | 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2 | 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2 | construction |
| place/R1 | by design: mean number of water oxygens in the place of the ring | 1 | 1 | construction |
| place/R1 | by design: smallest distance of a water oxygen to the place of the ring (Å) | 1.382 | 1.38237 | construction |
| place/R1 | by design: nearest water oxygen to the place of the ring in every full frame, largest difference (Å) | 0 | 0.0004 | construction |
| place/R1 | by design: water oxygens at the crystal ring positions, full frames that differ | 0 | 0 | construction |
| place/R1 | by design: mean number of water oxygens at the crystal ring positions | 1.75 | 1.75 | construction |
| place/R1 | by design: equilibration, the steps in the order of the run | eq01, eq02, eq03, eq03x1, eq03x2, eq03x3, eq04, eq05, eq06 | eq01, eq02, eq03, eq03x1, eq03x2, eq03x3, eq04, eq05, eq06 | construction |
| place/R1 | by design: equilibration, solute frames, full frames, length (ps) | 350, 35, 3500 | 350, 35, 3500 | construction |
| place/R1 | by design: equilibration, time of the first and of the last solute frame (ns) | -3.49, 0 | -3.49, 0 | construction |
| place/R1 | by design: equilibration, stage, extension step and time of solute frame 170 (the first of eq03x3) | 3, 3, -1.79 | 3, 3, -1.79 | construction |
| place/R1 | by design: equilibration, water oxygens in the place of the ring in every full frame | 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 | 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 | construction |
| place/R1 | by design: equilibration, nearest water oxygen to the place of the ring in every full frame, largest difference (Å) | 0 | 0 | construction |
| place/R1 | by design: equilibration, water oxygens in the place of the ring (the count the driver acts on) at the end of eq03, eq03x1, eq03x2, eq03x3 | 0, 0, 0, 1 | 0, 0, 0, 1 | construction |
| place/R1 | by design: equilibration, these steps were compared with the driver (their last full frame is at the end of the step) | yes, yes, yes, yes | yes, yes, yes, yes | construction |
| place/R1 | by design: equilibration, the count of the driver at the end of these steps agrees with the analysis | yes, yes, yes, yes | yes, yes, yes, yes | construction |
| place/R1 | by design: equilibration, the nearest water oxygen of the driver at the end of these steps agrees with the analysis | yes, yes, yes, yes | yes, yes, yes, yes | construction |
| default/R1 | warnings | 0 | 0 | construction |
| default/R1 | by design: site atoms and cut-off of the count that the driver acts on (Å) | C1, C2, C3, C4, C5, C6, 2.5 | C1, C2, C3, C4, C5, C6, 2.5 | construction |
| default/R1 | by design: where these settings came from, and whether they are the defaults of the specification | specification (defaults), yes | specification (defaults), yes | construction |
| default/R1 | by design: water oxygens in the place of the ring in every full frame | 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2 | 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2 | construction |
| default/R1 | by design: mean number of water oxygens at the crystal ring positions | 1.75 | 1.75 | construction |

The rigidly moved crystal ligand has, by the second formulas, RMSD 3.562 Å, d1 9.556 Å, d2 2.634 Å (state R) and keeps τ_fold 132.056° and S–N2 1.5963 Å; the translated ligand has RMSD 3.000 Å, d1 6.897 Å, d2 6.732 Å (state I). The numbers of gate G-V of X2 are those of a production run and carry no verdict; its full fraction of 10 % would not meet GV3 d, its crystal angle of 98.4° now meets GV3 c. In the crystal-like frames of the runs with `asn67_ND2` the flag is 1 except where the amide was turned, in the rotated-like frames of gate/R1 it is 1 only where the amide was turned (section 4).

### 16.3 Does the test notice a wrong analysis?

83 copies of `analyse_md.py`, each with one rule of this specification broken, were run on the same synthetic runs and compared with the planted values (`mutate.py`, `mutation_results.json`): the 14 of version 1.0, the 26 of version 1.1, the 34 for the reported quantities of version 1.2 (sections 9A to 9E; those of section 9D follow the new code) and 9 new ones for the count that the driver acts on, its settings and its comparison with the driver's record. One copy of version 1.2 is reversed: there "the six ring carbons only" was the broken rule; now the count of the driver on seven positions is, and the wider count without the aniline N is a broken copy of its own. 

The machine was heavily loaded while this version was tested (one pass of the unbroken script over all runs and the comparison then took 10 minutes for 101 s of processor time; `timing_one_copy.txt`), and the copies were run in two ways. The first 16 of the list (M01 to M16; the measure with N − 1 described below is one of them) were run on all 23 runs and compared with all 7339 planted values, as in version 1.2 (`mutate.py`). The other 68 were run in stages by a rule fixed in `mutate_staged.py` before any of them was run. Stage A: the broken copy analyses only the runs named for it, every other run keeps the output of the unbroken script, and the planted values of these runs are compared (with the values typed for single runs and for the campaign, which the comparison always reads). The runs named are, for a copy that version 1.2 had, the two smallest runs in which that copy showed a difference in the test of version 1.2 (the campaign table alone if it showed only there; ext/R1 and glu/X0 if it stopped the analysis), and for the copies of section 9D and of the comparison with the driver's record that are new or follow new code (the eight from "crystal ring positions not carried into the frame" to "the count of the driver on seven positions" and the nine new ones) the three runs default/R1, place/R1 and ext/R1. Stage B, only if stage A returns every planted value: the copy on all runs, as for the first copies. A copy is noticed if a planted value is no longer returned or if the analysis stops. A copy noticed in stage A would be noticed on all runs as well, because the comparison on all runs contains the comparisons of stage A; the number of values given for it is that of the runs analysed and is a lower bound. Before the stages the unbroken script was checked in the same way: on all runs it returns 7339 of 7339 planted values, and through stage A on default/R1, place/R1, ext/R1 with the campaign table 1379 of 1379. 68 copies were noticed in stage A, 0 went on to stage B. All 83 copies were noticed:

| Broken rule | The copy analysed | Planted values no longer returned (of those compared) | Runs in which a value differs |
|---|---|---|---|
| no symmetry correction of the ligand RMSD | all 23 runs | 224 (of 7339) | 17 runs |
| no fit on the CA atoms | all 23 runs | 953 (of 7339) | 23 runs |
| no imaging of zinc and ligand | all 23 runs | 2433 (of 7339) | 24 runs |
| water distances without minimum image | all 23 runs | 678 (of 7339) | 23 runs |
| 90 % rule strict (more than 90 %) | all 23 runs | 3 (of 7339) | campaign, X1 |
| hydrogen-bond angle at the donor, not at the hydrogen | all 23 runs | 166 (of 7339) | 22 runs |
| time axis from zero | all 23 runs | 6 (of 7339) | conv/R1, flip/X1, gate/R1, R1, X1, X3 |
| dwell of 0.4 ns instead of 1 ns | all 23 runs | 14 (of 7339) | X1, X3 |
| τ_fold from one torsion only | all 23 runs | 165 (of 7339) | 17 runs |
| intermediate frames not carried | all 23 runs | 13 (of 7339) | X1, X3 |
| last half taken as the first half | all 23 runs | 134 (of 7339) | 11 runs |
| midpoint of d1 at 7.0 instead of 6.7 Å | all 23 runs | 43 (of 7339) | flip/X1, X1, X3 |
| fifth ligand at 2.3 instead of 2.6 Å | all 23 runs | 142 (of 7339) | 22 runs |
| solute index taken as the index of the full system | all 23 runs | the analysis stops with an error | all analysed |
| N–Zn–N angles in the window 100–120° (the rule of version 1.0) | all 23 runs | 53 (of 7339) | 21 runs |
| N–Zn–N angles within 12 instead of 10° | angle/X0 | 6 (of 334) | angle/X0 |
| N–Zn–N angles within 8 instead of 10° | edge/X0 | 6 (of 334) | edge/X0 |
| references of two N–Zn–N angles exchanged (108.9 and 98.4) | ext/R1, scaled/N1 | 8 (of 987) | ext/R1, scaled/N1 |
| fifth ligand: no frame allowed (the rule of version 1.0) | ext/R1, scaled/N1 | 4 (of 987) | ext/R1, scaled/N1 |
| fifth ligand: up to 10 % instead of 5 % of the frames | sulfonyl/X0, water/X0 | 14 (of 642) | sulfonyl/X0, water/X0 |
| fifth ligand, solute trajectory: fewer than 5 % (limit not included) | edge/X0 | 7 (of 334) | edge/X0 |
| fifth ligand, full trajectory: fewer than 5 % (limit not included) | ext/R1, scaled/N1 | 4 (of 987) | ext/R1, scaled/N1 |
| fifth ligand: the solute trajectory is not counted | sulfonyl/X0 | 5 (of 349) | sulfonyl/X0 |
| fifth ligand: the full trajectory is not counted | water/X0, test/N1 | 6 (of 629) | N1, water/X0 |
| fifth ligand: one fraction over the frames of both trajectories together | sulfonyl/X0, water/X0 | 9 (of 642) | sulfonyl/X0, water/X0 |
| fifth ligand, full trajectory: water and protein oxygen only, no sulfonyl oxygen | edge/X0, water/X0 | 13 (of 627) | edge/X0, water/X0 |
| fifth ligand: fractions over the last half of the run | ext/R1, scaled/N1 | 6 (of 987) | ext/R1, scaled/N1 |
| C4–S–N2 within 4 instead of 3° | csn/X0 | 6 (of 334) | csn/X0 |
| C4–S–N2 within 2 instead of 3° | edge/X0 | 6 (of 334) | edge/X0 |
| τ_fold within 25 instead of 20° | tau/X0 | 6 (of 334) | tau/X0 |
| τ_fold within 15 instead of 20° | scaled/N1, edge/X0 | 8 (of 764) | edge/X0, scaled/N1 |
| C4–S–N2 and τ_fold reported only, not part of the verdict (version 1.0) | csn/X0, tau/X0 | 6 (of 612) | csn/X0, tau/X0 |
| run X0 not treated as the validation run | short/X0, angle/X0 | 13 (of 627) | angle/X0, short/X0 |
| validation run: planned length from the manifest, then 100 ns | short/X0, csn/X0 | 6 (of 627) | csn/X0, short/X0 |
| run X0 entered in the campaign table | campaign table (from the summaries) | 7 (of 56) | campaign |
| gate G-V with a verdict for every run of start X (version 1.0) | flip/X1, test/X1 | 4 (of 802) | flip/X1, X1 |
| rule 1 applied to the validation run as well | short/X0, angle/X0 | 2 (of 627) | angle/X0, short/X0 |
| Asn67 flag reversed (set where OD1 is nearer) | ext/R1, scaled/N1 | 10 (of 987) | ext/R1, scaled/N1 |
| distance to Asn67 ND2 taken to OD1 | ext/R1, scaled/N1 | 11 (of 987) | ext/R1, scaled/N1 |
| orientation distance d2 to the nearer of OD1 and ND2 | scaled/N1, flip/X1 | 22 (of 949) | flip/X1, scaled/N1 |
| signed volume at N2 with HN21 and S exchanged (the sign reversed) | ext/R1, scaled/N1 | 109 (of 987) | ext/R1, scaled/N1 |
| signed volume at N2 from vectors that are not normalised | ext/R1, scaled/N1 | 27 (of 987) | ext/R1, scaled/N1 |
| sign at N2 of the rotated start taken as +1 | ext/R1, gate/R1 | 14 (of 851) | ext/R1, gate/R1 |
| a sign change at N2 counts after 0.5 ns instead of 1 ns | conv/R1 | 10 (of 531) | conv/R1 |
| a stretch of exactly 1 ns is not long (more than 100 frames asked) | conv/R1 | 10 (of 531) | conv/R1 |
| every change of sign at N2 counted (short stretches not left out) | short/X0, gate/R1 | 14 (of 643) | gate/R1, short/X0 |
| time of a sign change at N2 counted from zero | flip/X1, conv/R1 | 6 (of 994) | conv/R1, flip/X1 |
| sign at N2 counted per assigned state, not per state of the frame | flip/X1, test/X1 | 10 (of 802) | flip/X1, X1 |
| sign at N2 of the last 20 ns counted over the whole run | conv/R1, test/R1 | 4 (of 814) | conv/R1, R1 |
| hydrogen bond split at the first frame in the opposite state (not at the first long visit) | test/X1, test/X3 | 8 (of 622) | X1, X3 |
| the frame of the first passage counted with the frames before it | flip/X1, gate/R1 | 9 (of 813) | flip/X1, gate/R1 |
| hydrogen bond of Thr199 N–H to the sulfonyl oxygens split instead of N2–H···OG1 | flip/X1, gate/R1 | 7 (of 813) | flip/X1, gate/R1 |
| OE1 and OE2 of Glu106 exchanged | ext/R1, scaled/N1 | 36 (of 987) | ext/R1, scaled/N1 |
| Glu106 counted within 2.8 instead of 2.6 Å | conv/R1 | 6 (of 531) | conv/R1 |
| Glu106 counted within 2.52 instead of 2.6 Å | conv/R1 | 5 (of 531) | conv/R1 |
| Glu106 made a criterion: counted in the solute fraction of the fifth ligand of gate G-V | glu/X0, conv/R1 | 9 (of 840) | conv/R1, glu/X0 |
| crystal ring positions not carried into the frame (no fit of the full frame) | default/R1, place/R1, ext/R1 | 134 (of 1379) | default/R1, ext/R1, place/R1 |
| the fit applied the wrong way round (rotation not inverted) | default/R1, place/R1, ext/R1 | 128 (of 1379) | default/R1, ext/R1, place/R1 |
| ring positions of the ligand in the frame instead of the crystal positions | default/R1, place/R1, ext/R1 | 149 (of 1379) | default/R1, ext/R1, place/R1 |
| the wider count within 3.3 instead of 3.5 Å | default/R1, place/R1, ext/R1 | 11 (of 1379) | ext/R1 |
| the wider count within 3.7 instead of 3.5 Å | default/R1, place/R1, ext/R1 | 38 (of 1379) | default/R1, ext/R1, place/R1 |
| pairs of water and position counted instead of water molecules | default/R1, place/R1, ext/R1 | 32 (of 1379) | default/R1, place/R1 |
| water at the ring positions without the minimum image | default/R1, place/R1, ext/R1 | 115 (of 1379) | default/R1, ext/R1, place/R1 |
| the count of the driver on seven positions (the aniline N added, as version 1.2 counted) | default/R1, place/R1, ext/R1 | 12 (of 1379) | default/R1, place/R1 |
| equilibration time counted from the start of the equilibration (positive times) | ext/R1, scaled/N1 | 34 (of 987) | ext/R1, scaled/N1 |
| first frame of a stage put at the start of the stage | ext/R1, scaled/N1 | 32 (of 987) | ext/R1, scaled/N1 |
| stage table of the protocol used although the manifest lists the stages | ext/R1, scaled/N1 | 65 (of 893) | ext/R1, scaled/N1 |
| equilibration frames read as frames of the production | ext/R1, scaled/N1 | 76 (of 987) | ext/R1, scaled/N1 |
| first production frame said to differ only when it is in the opposite state | flip/X1 | 3 (of 519) | flip/X1 |
| sign at N2 of the first production frame compared with +1 for every start | ext/R1, gate/R1 | 3 (of 851) | ext/R1, gate/R1 |
| no imaging of zinc and ligand in the equilibration frames | ext/R1, scaled/N1 | 94 (of 987) | ext/R1, scaled/N1 |
| extension steps of the held part left out (only the files eqNN read) | ext/R1 | 24 (of 463) | ext/R1 |
| extension steps read after all stages, not after their own stage | ext/R1, glu/X0 | the analysis stops with an error | all analysed |
| the last full frame of a step, and with it the count set beside the count of the driver, taken as its first full frame | ext/R1, flip/X1 | 22 (of 1020) | ext/R1, flip/X1 |
| the wider count without the aniline N (six positions instead of seven) | default/R1, place/R1, ext/R1 | 45 (of 1379) | default/R1, ext/R1, place/R1 |
| the count of the driver within 2.7 instead of 2.5 Å (default cut-off) | default/R1, place/R1, ext/R1 | 9 (of 1379) | default/R1, place/R1 |
| the count of the driver within 2.3 instead of 2.5 Å (default cut-off) | default/R1, place/R1, ext/R1 | 10 (of 1379) | default/R1, place/R1 |
| site atoms and cut-off of the manifest not read (the defaults for every run) | default/R1, place/R1, ext/R1 | 30 (of 1379) | ext/R1, place/R1 |
| the cut-off of the manifest read, its site atoms not | default/R1, place/R1, ext/R1 | 27 (of 1379) | ext/R1 |
| the site atoms of the manifest read, its cut-off not | default/R1, place/R1, ext/R1 | 17 (of 1379) | ext/R1 |
| the count of the driver compared with the wider count (the comparison of version 1.2) | default/R1, place/R1, ext/R1 | 6 (of 1379) | place/R1 |
| the nearest distance of the driver compared with that of the wider count | default/R1, place/R1, ext/R1 | 2 (of 1379) | place/R1 |
| no warning when the settings of a run are not those of the plan | default/R1, place/R1, ext/R1 | 1 (of 1379) | ext/R1 |

One further change is run as a measure of what the tolerance of section 14 lets through: a standard deviation with N − 1 instead of N in the denominator. It changes 7 of 7339 planted values by more than the tolerance, all of them standard deviations, in gate/R1 (400 frames): 3; flip/X1 (300 frames): 4. The largest change is in `d1_N1_F131CZ_A` of flip/X1 (2.8056 planted, 2.8103 returned). All these runs are short and hold two poses or two values that lie far apart; the change falls with the number of frames, so for a production run of 10,000 frames the comparison of section 14 cannot tell the two denominators apart: the two programs must agree on N from this text.

### 16.4 Other checks

- **Version 1.3 against version 1.2 on the 21 runs of version 1.2** (`v12_against_v13.json`; the test folder of version 1.2 was kept). The generator writes the same runs as before: the coordinates of all 24625 frames in the 113 trajectory files, those of the equilibration included, are identical to those of the test of version 1.2 (largest difference 0.0 Å), and the planted values agree as text in every column that version 1.2 had. On these runs the analysis of version 1.3 writes the same text in every column of the five CSV files that version 1.2 had, in the same order, two of them under their new names, and the same content under every key of `summary.json` that version 1.2 had (`stats`, `fractions`, `orientation`, `validation`, `gate_GV`, `rmsf`, `asn67_ND2`, `n2_configuration`, `hbond_first_passage`, `start_of_production`, `equilibration`, the labelling counts, the table row, the constants, the checks of `run`), and the same campaign rules and the same folders left out. Renamed: `n_water_ring_site_3p5A` to `n_water_ring_wide_3p5A`, `d_water_ring_site_min_A` to `d_water_ring_wide_min_A`, `n_water_ring_site_3p5A_mean` to `n_water_ring_wide_3p5A_mean`, `d_water_ring_site_min_min_A` to `d_water_ring_wide_min_min_A`. New are 2 columns of `timeseries_full.csv` and of `timeseries_equil_full.csv` (`n_water_ring_place`, `d_water_ring_place_min_A`; none in `timeseries_solute.csv`), the key `water_ring_place`, 2 entries of the table row (`n_water_ring_place_mean`, `d_water_ring_place_min_min_A`) and 2 numbers in `last_full_frame` of the equilibration steps. The number of warnings differs in one run, as the change of its manifest asks: ext/R1 0 then, 1 now.
- **The equilibration changes nothing else** (`run_all.log`, section "edge cases"). Run conv/R1 was analysed a second time without its equilibration files: 3 of 3 CSV files of the production are byte-identical, and of `summary.json` only `equilibration`, `inputs_sha256` (the checksums of the files read), the folder name in `run` and the entry `equilibration_frames_solute` of the table row differ; no equilibration series is written and there is no warning.
- **Same files again.** Three analyses of run X1, of the validation run edge/X0, of conv/R1 and of place/R1 gave 11, 11, 14 and 14 byte-identical output files (CSV, JSON, PNG; `repro/`). The synthetic DCD and PDB files themselves carry the time of writing in a header line, so their checksums change from one generation to the next; the coordinates do not (seed 20261005).
- **Runs that are not as expected** (`run_all.log`, section "edge cases"): no full trajectory (water quantities missing, the two water counts of section 9D included, one warning, the full fraction does not exist and GV3 d is not met); a last segment that holds only a header, and an empty file (left out, the rest analysed); a folder whose name is not a start and has no manifest (stops with a message); a folder without a trajectory (stops with a message); a start given on the command line against the folder name (warning; with start R the sign of the start is −1 and the first frame of the crystal pose differs in state and in sign).
- **Equilibration files that are not as expected** (same log): a file `solute_eq07.dcd` beside the six stages of the manifest is not read (1 warning, the 8 frames of the six stages are read); the files of the shortened equilibration read without the stage list of the manifest fall back on the table of the protocol and are flagged (17 warnings, 12 on the number of frames and 5 on the time axis by the DCD headers); with one stage file of conv/R1 cut down to its header and one emptied, 120 of 200 frames are read with 3 warnings, the times and values of the other stages are unchanged (first and last time -1.99 and -0.5 ns) and the production series is byte-identical; the run with extension steps read without the stage list of the manifest loses the two extension files (2 warnings on files that belong to no step, 200 of 300 frames read, 9 warnings in all, 6 of them on the time axis of the stages after the extension and 1 on the settings of its manifest); and when the count of the driver recorded in the manifest of ext/R1 for eq03x2 is changed from 1 to 2, the analysis gives 1 warning on that count (beside the 1 on the settings) and writes `count_agrees` true, true, false for the three steps.
- **The settings of the count of the driver** (same log; run place/R1 with other entries under `hydration_check` of its manifest, the frames and the step records unchanged). Without `settings` (a protocol without the hydration check) and without the key at all: the defaults, no warning, the same counts, and the driver's four records still agree. With a site atom that the ligand does not have (C9), or a cut-off of zero: the defaults, one warning each, the same counts. With seven site atoms and 3.5 Å in place of what the records were counted with: the analysis follows the manifest (`n_water_ring_place` becomes the wider count), warns that the settings are not those of the plan, and warns 4 times that the driver's records do not agree (`count_agrees` false at all four steps, `nearest_agrees` false at the last). With the six ring carbons and 3.5 Å (the protocol between 8:01 and 8:56 AM): 3 disagreements. A manifest whose settings do not belong to its records is therefore noticed.
- **The role of a run, further cases** (same log): a folder named `X0` whose manifest says replica 2 is analysed as a production run, with the planned length of its manifest and a warning that name and manifest disagree; the validation run with `--planned-ns 2` takes 2 ns from the command line and warns that it holds more frames; in a campaign folder with `X0`, `X2` and a folder `X1` whose manifest says `X0`, the table holds X2 alone and both others are listed as left out; `--reuse` gives the same campaign table from the five existing summaries, and analyses a run again when its summary names another version of this specification.
- **Rules that cannot be planted in stored frames** (`unit_checks.py`, `unit_checks.json`): 32 of 32 checks of the functions of sections 9A, 9B, 9D and 9E on values typed by hand pass: a signed volume of exactly zero (which single-precision coordinates never give), a first passage at the first frame, a run shorter than 1 ns, the frames expected in a short stage, the stage table with and without a manifest; and for section 9D the settings taken from a manifest or from the defaults in 13 cases (no manifest; no hydration check in the protocol; the settings of today; the ring carbons in another order; the settings of 8:01 AM and those before; and 7 kinds of settings that cannot be used). The first are the worked examples of sections 9A and 9B.
- **A trajectory of the real driver, with its equilibration.** `md/run_md.py` (sha256 `45e4cf2872beddbf…`) had been run by another agent on a stand-in system (4600 atoms, 111 in the solute; stages of 2 fs and of 4 fs steps); a copy of its finished test run (`driver_snapshot/`, 20 solute and 10 full frames of production in two segments, 40 and 20 of equilibration in six stages) is read by version 1.3 without a warning: the manifest keys are found, the DCD headers are contiguous, the segment list agrees with the files, the checksums of the system files agree, Zn–N2 of the two trajectories agrees within 0.000001 Å. For the equilibration the stage table comes from the manifest (4 ps), every stage file holds the frames expected and recorded, the DCD headers put the first frame of each stage one interval after the start of the stage, and the last frame of every stage has the time of the end of that stage (first frame -0.0039 ns, last 0 ns). The 56 distances that the driver itself monitors at the end of each production segment and of each equilibration stage (14 and 42; Zn–N2, N2···OG1, d1, d2, the three Zn–N(His)) agree with the analysis within 0.0005 Å (the driver gives three decimals). File: `driver_crosscheck.json`. The stand-in system is not hCA II; this checks the reading of the files, not the chemistry. That run was made by the driver before its hydration check (no extension step, no count in the manifest).
- **Two runs of the driver as it stands, with its hydration check and its extension steps** (`driver_now/`, `run_driver.sh`, `driver_crosscheck.json`, key `hydration_runs_of_the_driver_of_today`). Version 1.2 could test the files of extension steps on a synthetic run only. For this version `md/run_md.py` (sha256 `6d797f080ff75679…`) was run twice on 6 Oct 2026, 10:05 AM IST, on the stand-in systems `dry` and `wet` of its maintainer's tests (`/root/cqc-sfa-md/test_driver/hydration/systems/`, read only; 3120 atoms, 161 in the solute, 10 CA atoms; not hCA II), one thread, CPU platform, with a test protocol that is `md/protocol.json` as it stands (sha256 `decda42f5340d06f…`: site atoms C1, C2, C3, C4, C5, C6, cut-off 2.5 Å) with output every 0.1 ps, and with `--equil-scale 0.001`, so that every step ends on a full frame. In `dry` the driver found no water at the site and made 4 extension steps; in `wet` it counted 9 water oxygens at the end of eq03 and went on. The analysis reads both runs without a warning: the settings come from the manifest and are the defaults; the steps are read in the order of the manifest (eq01, eq02, eq03, eq03x1, eq03x2, eq03x3, eq03x4, eq04, eq05, eq06 in `dry`), with the start, the end and the frames of every step as recorded, 4 extension steps among them; the last full frame of every step lies at the end of the step. All 6 hydration records of the driver are compared and agree in count and in nearest distance (`compared`, `count_agrees`, `nearest_agrees` true). Beyond what the analysis itself compares, the test sets the count that the driver notes at the end of every step (`monitor.crystal_ring_site`, also for the production segments) beside the analysis: 18 of 18 counts are identical (in `wet`: 6, 6, 9, 8, 7, 6, 6; in `dry` all zero, with the nearest water oxygen 7.89 to 13.75 Å away), and the nearest distances agree within 0.0005 Å (the driver gives three decimals). The DCD headers of these runs put the first frame of every step at S + Δt exactly.
- **The test runs of the driver's maintainer** (`/root/cqc-sfa-md/test_driver/hydration/runs/`, made on 6 Oct 2026 between 5:35 and 6:32 AM IST by the driver as version 1.2 named it (sha256 `b0f3e2b23b3eedf5…`) with the hydration block of that hour: 7 site atoms, 3.5 Å; read only, output in the scratch folder; `run_all.log`, `driver_0528.json`). The runs `dry`, `arrive`, `wet`, `skip` and `late` are read with exit code 0: the steps and their frames are those of the manifest (4 extension steps in `dry`, 3 in `arrive`), the settings are taken from each manifest (seven atoms and 3.5 Å; in `arrive` 11.2 Å, a device of that test) with the warning that they are not those of the plan, and the driver's counts are copied (0, 0, 0, 0, 0 in `dry`; 0, 0, 0, 1 in `arrive`; 15 in `wet`). None of them is compared: in these runs no step ends on a full frame (stage 03 lasts 0.15 ps and an extension step 0.25 ps, with a full frame every 0.2 ps), which is why the two runs above were made. Four of the five runs (`dry`, `wet`, `skip`, `late`) get 6 warnings of the guard "header" of section 9E each, with a right time axis: their stage 04 starts at 0.35 ps (87.5 steps) or at 1.35 ps (337.5 steps), not on a whole number of its 4 fs steps, and the driver rounds the first step of the header to a whole number (the header of `solute_eq04.dcd` of `dry` says 1.452 ps for a first frame at 1.45 ps). In `arrive` stage 04 starts at 1.1 ps, a whole number of steps, and there is no such warning. With the protocol of the campaign every step starts on a whole number of its time steps (section 9E, guard "header").
- **Appendix A read by a second reader** (`appendix_check.py`, `appendix_check.json`). A small reader was written from the text of appendix A and of sections 2, 9D and 9E, with the Python standard library and NumPy only (its own reading of the DCD and PDB files and of the manifest, its own fit and its own minimum image; no MDAnalysis and nothing of `analyse_md.py` or `run_md.py`). On the two runs of the driver above and on the synthetic run place/R1 it makes both water counts in all 123 full frames of the equilibration and of the production: every count is the count of `analyse_md.py`, the nearest distances agree within 0.00005 Å, the stage, the extension number and the time of every equilibration frame are the same, the DCD headers put the first frame of every file one interval after the start of its step, and the number of frames found from the file size is the NSET of the header. So the layout, the lattice vectors, the formula of the site positions and the time axis of the appendix are enough to make the count. This is no substitute for the second program of stage S6.
- **Time.** Not measured again. Version 1.3 adds, per full frame, one more set of distances from six positions to the water oxygens (the fit on the CA atoms is shared with the wider count); the figure of version 1.0 stands for the size of the first sketch (10,000 solute frames of 4,100 atoms, 1,000 full frames of 25,600 atoms: about 26 s per run with one thread, `probe_speed.json` of the test of version 1.0; random coordinates of that size, not a real run). The built systems hold 40,900 atoms; by the atom count alone the full frames would take about 1.6 times as long to read (an estimate, not measured).
- **Figures** were drawn for all synthetic runs. The code that draws those of version 1.1 is unchanged; in the two of version 1.2 (`fig_n2_site.png`, `fig_equilibration.png`) the panels of the water count now show both counts and the cut-off; they were looked at for place/R1 and conv/R1.

## 17. Changes to this specification

- 5 Oct 2026, about 9:30 PM IST: version 1.0, before any run.
- 5 Oct 2026, 10:55 PM IST: version 1.1, before any validation or production run (`/root/cqc-sfa-md/runs/` does
  not exist). The changes are listed at the top of this file ("Changes against version 1.0"): N–Zn–N angles
  against their values in 6RL9; the fifth ligand as two fractions of at most 5 %; C4–S–N2 and τ_fold in the verdict
  of gate G-V; the validation run as run X0; the optional key `asn67_ND2`; the tests extended.
- 6 Oct 2026, 7:18 AM IST: version 1.2, before any validation or production run (`/root/cqc-sfa-md/runs/` does
  not exist; 6 Oct 2026, 5:47 AM and 7:17 AM IST). The changes are listed at the top of this file ("Changes of
  version 1.2 against version 1.1"):
  reported quantities only. The configuration at N2 (signed volume, sign, lasting sign changes); the hydrogen bond
  N2–H···OG1 before and after the first passage; the zinc to OE1 and OE2 of Glu106 by name; water at the crystal
  positions of the ring; the equilibration frames, extension steps of the held part included, in a series of their
  own and the first production frame against the start as built; the tests extended. No verdict, gate, rule or
  tolerance was changed.
- 6 Oct 2026, 11:26 AM IST: version 1.3, before any validation or production run (`/root/cqc-sfa-md/runs/` holds 3 folders of tests of the rented machine, `bench`, `preflight_a` and `preflight_b`, and no run of the campaign; 6 Oct 2026, 9:56 AM and 11:22 AM IST). The changes
  are listed at the top of this file ("Changes against version 1.2"): the water count of section 9D is now the
  count that the driver acts on (the water oxygens within 2.5 Å of the crystal positions of the six ring carbons;
  site atoms and cut-off from the manifest of the run where it gives them), in every full frame of the
  equilibration and of the production, with the distance of the nearest one, and it is the count set beside the
  driver's record at the end of stage 03 and of every extension step; the count of version 1.2 stays under the
  names `n_water_ring_wide_3p5A` and `d_water_ring_wide_min_A`; appendix A on the files of the driver; wordings
  made exact; the tests extended. This version follows `md/run_md.py` with sha256 `6d797f080ff75679…` and
  `md/protocol.json` with sha256 `decda42f5340d06f…`. No verdict, gate, rule or tolerance was changed.

## Appendix A. Reading the driver's files

What a program that reads a run folder has to know of `md/run_md.py` (sha256 `6d797f080ff75679…`), so that it need
not open that file or `analyse_md.py`. Everything here is a statement about files; the definitions are in
sections 1 to 14.

### A.1 The files of a run folder

`<step>` is the name of a step of the run: `min` (minimisation), `eq01` to `eq06` (the equilibration stages),
`eq03x1` to `eq03x4` (extension steps of stage 03, present only if the hydration check asked for them; they run
between `eq03` and `eq04`), `seg0001`, `seg0002`, … (production segments, 5 ns each in the campaign; the validation
run X0 has one segment of 5 ns).

| File | Read | What it is |
|---|---|---|
| `manifest.json` | the keys of A.4 | the record of the whole run, written anew after every step |
| `solute.pdb`, `full.pdb` | names only | one `ATOM` line per atom, in the atom order of the DCD files of the same kind. Fixed columns: atom name 13 to 16, residue name 18 to 20, coordinates 31 to 54 (those of the built system; not used). Serial and residue numbers are cut to the width of their columns and repeat in the full system, so an atom is known by the position of its line (first atom line = atom 0), never by its serial number. The driver writes no chain letter and no `TER` or `HETATM` line; the files of the synthetic test (section 16.2) are written by OpenMM's PDB writer and have all three, so a reader takes every line that begins with `ATOM` or `HETATM` and passes over the others |
| `solute_<step>.dcd` | yes | the atoms of `selections.json` `solute`, in list order, one frame every 10 ps. Coordinates as integrated: the molecules are whole and are never moved into the cell, so the zinc or the ligand can stand in another periodic image than the protein (step 1 of section 2) |
| `full_<step>.dcd` | yes | all atoms in the order of the topology (index = position), one frame every 100 ps. The solute as in the solute file; every water molecule and every ion moved as a whole, by lattice vectors, to the image nearest to the centre of the solute. Distances to water therefore need the minimum image (section 2) |
| `rec_<step>.json` | no | the record of one step; the manifest holds the same record |
| `log_<step>.csv`, `state_<step>.xml`, `checkpoint_<step>.chk`, `final_solute.pdb`, `run_md.log`, `TEST_ONLY`, `benchmark.json` | no | state log, saved states, the console log, the mark of a test run, the speed test |

There is no DCD file for `min`. The files of a kind have the same atoms in every step. A step is complete when its
record exists; the driver writes a DCD file frame by frame, so the last file of a run in progress can hold the
header only (a step that has just begun) or end in the middle of a frame, and the manifest lists complete steps
only.

### A.2 The DCD files

CHARMM layout as OpenMM writes it: little-endian, coordinates in Å as 4-byte floats, a unit cell in every frame.
Every record stands between two 4-byte integers that give its length in bytes.

| Record | Bytes of the record | Content |
|---|---|---|
| 1 | 84 | the four characters `CORD`; nine 4-byte integers: NSET (number of frames; see below), ISTART (step number of the first frame), NSAVC (steps between two frames), a step number near the last frame that a reader does not need, five zeros; DELTA as a 4-byte float (the time step in AKMA units: time step in ps = DELTA × 0.04888821); ten 4-byte integers: 1 (a unit-cell record in every frame), eight zeros, 24 |
| 2 | 164 | the integer 2 and two title lines of 80 characters |
| 3 | 4 | NATOM, the number of atoms |
| every frame | 48 | the unit cell: six 8-byte floats in the order a, cos γ, b, cos β, cos α, c (lengths in Å; cosines of the angles, not degrees) |
| | 3 × 4 NATOM | three records: all x, then all y, then all z |

The header is 276 bytes long, a frame 80 + 12 NATOM bytes. The number of frames is taken from the size of the
file, ⌊(size − 276) / (80 + 12 NATOM)⌋, not from NSET, so that a frame cut off at the end is left out; a file of
276 bytes holds no frame, and an empty file is no DCD file (section 9E, guard "unreadable").

- **Lattice vectors** of a frame, from its unit-cell record: a = (a, 0, 0), b = (b cos γ, b sin γ, 0),
  c = (c_x, c_y, c_z) with c_x = c cos β, c_y = c (cos α − cos β cos γ) / sin γ, c_z = √(c² − c_x² − c_y²). For
  the truncated octahedron of the campaign a = b = c and α = β = γ = 109.47°.
- **Times in the header.** The time of the first frame is ISTART × time step, the interval NSAVC × time step.
  The driver counts ISTART from time zero of the kind of step, in time steps of that file: for a production
  segment from the start of the production, for an equilibration step from the start of `eq01`, so that
  ISTART × time step = S + Δt with S the start of the step (for an equilibration step S rounded to a whole number
  of its time steps; section 9E, guard "header"). The time step is 2 fs in `eq01` to `eq03` and in the extension
  steps and 4 fs in `eq04` to `eq06` and in the production. A header with ISTART = NSAVC in every file (OpenMM's
  own reporter, where each file begins a new count) tells nothing of the place of the file in the run. The
  analysis takes its time axis from A.3, not from the headers; the headers are guards (sections 1 and 9E).

### A.3 The time axis

- **Production.** Time zero is the end of the equilibration. The driver writes a frame at the end of every
  interval, never at time zero: solute frame n, counted from 0 through all segments in ascending NNNN, is at
  (n + 1) × 10 ps, full frame j at (j + 1) × 100 ps. A segment of 5 ns holds 500 solute and 50 full frames; segment
  NNNN begins (NNNN − 1) × 5 ns after time zero (`production.segments[].time_start_ps`, in ps).
- **Equilibration.** Its own clock starts at the beginning of `eq01` and runs through all steps, extension steps
  included; `equilibration.stages[].time_start_ps` and `time_end_ps` are on that clock, in ps. Frame k (from 0) of
  the file of a step that starts at S is at S + (k + 1) Δt on that clock; section 9E moves it to the clock of the
  production by subtracting T_eq, the end of the last step.
- **The steps of the campaign's protocol**: `eq01` 0 to 100 ps, `eq02` 100 to 400 ps, `eq03` 400 to 700 ps, then
  each extension step 500 ps, then `eq04` 400 ps, `eq05` 400 ps, `eq06` 500 ps. Without an extension step T_eq is
  2000 ps; with one to four it is 2500, 3000, 3500 or 4000 ps. Every step ends on a solute frame and on a full
  frame. A test run made with `--equil-scale` has shorter steps, which the manifest gives and the constants of
  section 9E do not.

### A.4 The manifest keys used

All are optional for the analysis: section 1 says what replaces a missing one.

| Key | Kind | Meaning | Used in |
|---|---|---|---|
| `run_id` | text, for example `X1` | the name the run was started with | start and replica (section 1) |
| `system_dir` | text | the folder of the built system on the machine that made the run (absolute path) | system folder (section 1) |
| `build_json.start`, `build_json.model` | text | from `build.json` of the built system, if it has them | start, model |
| `production.target_ns` | number | the length asked for | planned length of a production run |
| `protocol.production.length_ns` | number | the protocol as used, copied in full under `protocol` | planned length, after the key before |
| `protocol.output.solute_frame_interval_ps`, `protocol.output.full_frame_interval_ps` | numbers | the frame intervals, 10 and 100 in the campaign | Δt |
| `production.segments[]` | list of step records, complete segments only, in order | of each: `tag` (`seg0001`), `time_start_ps`, `frames.solute`, `frames.full` | guard "manifest" (section 1) |
| `inputs.<file>.sha256` for `selections.json` and `reference.pdb` | text | checksums of the files of the built system as the run read them | guard "system files" |
| `equilibration.stages[]` | list of step records, complete steps only, in the order of the run, extension steps among them | of each: `tag` (`eq03`, `eq03x1`), `time_start_ps`, `time_end_ps` | step table (section 9E) |
| `equilibration.stages[].hydration_check` | object; only in the records of `eq03` and of the extension steps, and only if the protocol holds the hydration check | `n_water_O` (whole number), `nearest_water_O.distance_A` (number, three decimals), `applies` (true or false), `decision` (`extend` or `go on`) | the driver's own count (section 9E) |
| `hydration_check.settings` | object; only if the protocol holds the hydration check | `site_atoms` (list of atom names), `water_cutoff_A` (number, Å) | P and c of section 9D |

Keys that are there and are not read: `hydration_check.requested`, `status`, `counts[]` (the counts of the steps
once more), `n_extensions`, `final_decision`; in every step record `monitor.crystal_ring_site` (the driver's count
at the end of that step, also of the production segments: `n_water_O_within_<cut-off>_A`, `nearest_water_O_A`,
`ligand_mean_distance_from_site_A`; used for the test of section 16.4 only); `test_only`, `equil_scale`,
`script.sha256`, and everything else. The part of a manifest that sections 9D and 9E read looks like this (other
keys left out):

```
{"run_id": "R1",
 "hydration_check": {"requested": true,
                     "settings": {"after_step": "eq03", "site_atoms": ["C1", "C2", "C3", "C4", "C5", "C6"],
                                  "water_cutoff_A": 2.5, "applies_if_ligand_mean_distance_above_A": 2.0,
                                  "max_extensions": 4}},
 "equilibration": {"stages": [
   {"tag": "eq01", "time_start_ps": 0.0, "time_end_ps": 100.0},
   {"tag": "eq02", "time_start_ps": 100.0, "time_end_ps": 400.0},
   {"tag": "eq03", "time_start_ps": 400.0, "time_end_ps": 700.0,
    "hydration_check": {"n_water_O": 0, "nearest_water_O": {"atom_index": 12345, "distance_A": 4.148},
                        "applies": true, "decision": "extend"}},
   {"tag": "eq03x1", "time_start_ps": 700.0, "time_end_ps": 1200.0,
    "hydration_check": {"n_water_O": 1, "nearest_water_O": {"atom_index": 23456, "distance_A": 2.214},
                        "applies": true, "decision": "go on"}},
   {"tag": "eq04", "time_start_ps": 1200.0, "time_end_ps": 1600.0},
   {"tag": "eq05", "time_start_ps": 1600.0, "time_end_ps": 2000.0},
   {"tag": "eq06", "time_start_ps": 2000.0, "time_end_ps": 2500.0}]}}
```

(The numbers of this example are made up to show the layout; they are no result.)

### A.5 The hydration check of the driver, in brief

At the end of `eq03`, the last stage in which the heavy atoms of the solute are held, the driver carries the
positions that the site atoms have in `reference.pdb` into the frame of its coordinates by the least-squares fit
of the CA atoms (as step 2 of section 9D) and counts the atoms of `water_O` within the cut-off of any of them
(minimum image; the limit included). The check **applies** if the ligand's own site atoms are on average more than
2.0 Å away from these positions, that is for a ligand that does not lie in the crystal pose (start R); for the
starts X and N it is recorded with its count and does not act. If it applies and the count is zero, an extension
step of 500 ps with the settings of `eq03` follows and the count is made again at its end, at most four times;
then the run goes on in any case, and no water is placed by hand. The site atoms (C1 to C6) and the cut-off
(0.25 nm) are in `md/protocol.json`, key `equilibration.hydration_check`, and the driver repeats them in Å under
`hydration_check.settings` of the manifest. The analysis makes the same count in every full frame (section 9D)
and sets it beside the driver's record (section 9E); it does not apply the rule.
