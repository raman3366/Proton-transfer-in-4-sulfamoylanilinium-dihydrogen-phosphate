# Running the MD: one trajectory, the campaign, resuming

Stage S5 of `md/PLAN.md` (and the validation run of S4). Written on 5 Oct 2026, about 10:15 PM IST; revised on
6 Oct 2026, 5:05 to 6:50 AM IST: the hydration check of the empty pocket (section 1A), the validation run as
run X0, sizes for the built systems of 40,900 atoms, the order of work on the rented machine, the open points
that `md/TOPOLOGY.md` and the review of the set-up had left with the driver.
Revised on 7 Oct 2026, 12:15 AM IST, to what was run: the hydration check as the plan changed it on 6 Oct 2026
(section 1A; section 5, row "Hydration check"; section 7, point 10), the GPU and the files that ran on it
(section 5, row "GPU precision"; the first lines of section 6.2; section 7, points 1 and 2), the specification of
the analysis and the extension steps (section 7, point 6) and the state below. That revision was checked and
corrected at about 1:05 AM IST; what is said of the state of the runs is as of 12:15 AM IST.
Files: `md/run_md.py` (the driver), `md/protocol.json` (every setting, with units), this file.
Nothing here rents or starts anything; the commands were written for the day a GPU is there, which came on
6 Oct 2026 (state below).

**State.** The driver is tested on the Raspberry Pi with small stand-in systems, CPU platform only (section 6:
the tests of 5 Oct all passed; of the 81 tests of 6 Oct 79 pass, and the two that fail compare a segment
repeated on 6 Oct with its file of 5 Oct and fail with the old driver in the same way). The tests of 6 Oct were
made with driver sha256 `b0f3e2b23b3eedf5582a25572ef80b0c3c214edb9414e94273f13be1b7a74fcc` and protocol sha256
`d886fe5b65cecfedb8515fdfd8a89f93ce8e54c54e6fba6326c64a1f38d60d39` (section 6.2), before the plan changed the
hydration check. **Since 6 Oct 2026 the driver runs on the rented GPU** (CUDA platform, NVIDIA L40S, mixed
precision; `minimisation.platform` in the manifest of run X1): first the speed test and the pre-flight test of
section 2 (`/root/cqc-sfa-md/runs/bench`, `preflight_a`, `preflight_b`; restrained model, start X; 6 Oct 2026,
8:55:57 to 9:10:30 AM IST by their `run_md.log`), then the rehearsal (`md/test/REHEARSAL_2e.md`), the validation
runs X0 and the campaign (`md/PRODUCTION.md`). The speed test gave 866.39 ns per day, and the pre-flight test
ended `DIFFERENT` (25 files compared, 18 differ; `speed.json` and `preflight_interrupt.json` in
`/root/cqc-sfa-md/status/vm/status/`): a run on this GPU cannot be regenerated bit for bit (section 7, point 2).
These three test runs were made with the driver of the tests (`b0f3e2b2…`) and a third protocol file, sha256
`83695bea2cf9703f03ba0055429c7d0931b0e4cc53f508e62e02fa6aaca27856` (the file of commit `7a7a988`: the six ring
carbons, still 3.5 Å; `script.sha256`, `inputs.protocol.sha256` and `hydration_check.settings` in their
manifests). Section 7 keeps what was open before that. **The campaign runs with** driver sha256
`6d797f080ff7567923469dabd86af71ca34183f5c482c14d55bf8038eb9d1822` and protocol sha256
`decda42f5340d06f5f695b4c56b59a7bb796e5872b3a8af1c1a5b426f9c09d9a` (`script.sha256` and `inputs.protocol.sha256` in
`/root/cqc-sfa-md/runs/bonded_2e/X1/manifest.json`, the one run of the campaign begun by 7 Oct 2026, 12:15 AM IST;
the two files of this folder have these checksums). Against the driver of the tests the campaign's driver has
one more line (the header of a new trajectory file is written to disk at once; commit `2fd9cdc`), and the
protocol has the hydration check of section 1A (six ring carbons, 2.5 Å) where the protocol of the tests had
seven positions and 3.5 Å. The records of the tests of section 6.2 (`test_results.json` and `unit_checks.json`,
written before 7 AM IST on 6 Oct 2026) name the earlier two files, and their folder holds no later results.

## 1. What a run is

A run is a chain of steps: `min`, then `eq01` to `eq06` (the equilibration stages of the protocol; after `eq03`
there may be up to four extension steps `eq03x1` to `eq03x4`, section 1A), then `seg0001`, `seg0002`, …
(production, 5 ns each; 20 segments for 100 ns).

- Every step starts from the saved state of the step before it (`state_<step>.xml`: positions, velocities,
  box), in a new OpenMM Context, with seeds that depend only on the run seed and the name of the step.
- A step is complete when its record `rec_<step>.json` exists. The record is written last and holds the
  sha256 of every file of the step.
- `--resume` deletes the files of a step that has no record and repeats that step from the state before it.
  An interruption can therefore leave neither duplicated nor missing frames. On a platform that repeats
  itself bit for bit, the repeated step is the step an uninterrupted run would have made (section 6).
- The run seed is `--seed` or, without it, derived from the run id (rule in `protocol.json`, value in the
  manifest). X1, X2 and X3 start from the same built structure (each run minimises it itself) and differ in
  the starting velocities and in the random sequences of thermostat and barostat.
- The run id is a name and nothing else. The driver reads no start and no replica number from it; X0, X1, R3
  and `bench` are treated alike (section 6.2, group "run id").

Needs: Python 3, NumPy, OpenMM 8.1 or newer (tested with 8.6.1; use one version for all runs).
`python -m openmm.testInstallation` lists the platforms of an installation.

## 1A. The hydration check

Plan, entries of 6 Oct 2026, 1:03 AM (point 3) and 3:02 AM (point 8): the rotated start has no water where the
ring of the crystal pose lies, and an empty place beside the ligand could draw the ring back. The plan changed
the count twice on 6 Oct 2026, before the rehearsal, the validation runs and the campaign: at 8:01 AM (point 7)
to the positions of the six ring carbons C1 to C6, without the aniline nitrogen N1, and at 8:56 AM (point 6) to
2.5 Å in place of 3.5 Å. Only the two pre-flight test runs on the GPU (restrained model, start X; 6 Oct 2026,
9:07:22 and 9:09:35 AM IST) still counted within 3.5 Å of the six ring carbons (1 and 2 water oxygens; check
skipped; their `run_md.log`). This section describes the check as it is run. Settings: `protocol.json`, key
`equilibration.hydration_check` (`site_atoms` C1 to C6, `water_cutoff_nm` 0.25,
`applies_if_ligand_mean_distance_above_nm` 0.2, `extension_length_ps` 500, `max_extensions` 4); the manifest of
every run repeats them (`hydration_check.settings`).

**What is counted.** At the end of `eq03`, the last stage in which the heavy atoms of the solute are held
(700 ps into the equilibration), on the coordinates that the step saves:

1. The positions of the six ring carbons C1 to C6 of the ligand in `reference.pdb` of the built system (the
   crystal pose, the same file content for every start) are carried into the frame of the current coordinates
   by the least-squares fit of the protein CA atoms (`protein_CA` of `selections.json` against the CA atoms of
   `reference.pdb`, no weights). The protein is whole in the coordinates of the driver (never wrapped).
2. Count: the water oxygens (`water_O`) at most 2.5 Å from any of the six positions, each counted once;
   minimum-image distances in the cell of that moment.
3. The check **applies** only if the ligand's own six ring carbons are on average more than 2.0 Å from these
   positions (atom by atom, same names). Otherwise the ligand itself is at the site: the check is recorded as
   skipped, with its count.

**The rule.** If the check applies and the count is zero, the held part is extended by one step of 500 ps with
the settings of `eq03` (NPT, 2 fs, 300 K, k = 250 on the heavy atoms of protein, zinc and ligand), and the count
is made again at its end; at most four times (2 ns). After that the run goes on in any case, and the manifest
says that the site was without water when the solute was released. No water is placed by hand.

**The extension steps** are ordinary steps: `eq03x1` … `eq03x4`, each with its own record, state log,
trajectories, saved state and seeds (the seed rule takes the step name), restartable like every step.
`--equil-scale` shortens them like the stages. The equilibration time counts on through them (`eq04` then
starts at 1,200 ps instead of 700 ps after one extension); production time is zero at the end of `eq06` as
before.

**Where it is recorded.**

- `rec_eq03.json` and the record of every extension step, key `hydration_check`: the count, the waters counted
  (atom index, name, distance), the nearest water oxygen, the nearest one for each of the six positions, the
  distance of each ligand atom from its position and their mean, `applies`, `skipped`, the r.m.s. deviation of
  the CA fit, the number of extension steps made, `decision` (`extend` or `go on`), the reason in words and the
  name of the next step. The chain of steps is read from these records, so `--resume` finds it again.
- `manifest.json`, key `hydration_check`: the settings, all counts in order, `n_extensions`, `extension_ps`,
  `final_decision`, `site_without_water_when_the_solute_was_released`; `equilibration.n_extension_steps`.
- The `monitor` block of **every** record (minimisation, stages, segments) carries the count, the nearest water
  oxygen and the mean distance of the ligand from the site (`crystal_ring_site`), reported only.
- The console and `run_md.log`: one line per count, for example (run X1 of the campaign, `run_md.log`, 6 Oct 2026,
  10:16:50 PM IST)
  `eq03: hydration check: 0 water oxygen(s) within 2.5 A of the crystal ring site (nearest 3.26 A); ligand 0.26 A from the site on average; check skipped: the ligand itself is at the site; next step eq04`.
  Where the check applies and the count is zero, the last two parts read
  `no water oxygen at the site after 0 of at most 4 extension steps: the held part is extended; next step eq03x1`
  (the words of `hydration_decision` in the driver; the example that stood here before 7 Oct 2026 showed the
  3.5 Å count).

**On the real system, as of 7 Oct 2026, 12:15 AM IST** (count within 2.5 Å of the six ring-carbon positions, at
the end of `eq03`): in the rehearsal the check applied to the run from start R and found 3 water oxygens (nearest
1.22 Å; ligand 3.94 Å from the site), so no extension step was made, and it was skipped in the runs from X and N
(`md/test/REHEARSAL_2e.md`); it was skipped in X0 of both zinc sets (`md/PRODUCTION.md`, section 5.1) and in X1
(count 0, nearest 3.263 Å, ligand 0.259 Å from the site; `hydration_check.counts` in the manifest of X1). No run
from start R of the campaign had begun by then, and the extension branch had not run on the real system.

**The built systems at their start coordinates** (restrained model, `system.inpcrd`; the bonded systems will
have the same heavy atoms) [`unit_checks.json`]. **This table and the two remarks after it show the rule as it
was worded before 8:01 AM on 6 Oct 2026** (seven positions, C1 to C6 and N1, and 3.5 Å): the numbers were made
with that rule (driver `b0f3e2b2…`, protocol `d886fe5b…`) and are left as they were made. They are the finding that
the plan answered at 8:01 AM (point 7) and 8:56 AM (point 6); section 7, point 10, has the outcome. The rule in
force counts within 2.5 Å of the six ring carbons; the nearest water oxygen to a ring-carbon position in these
start coordinates is 3.025, 3.025 and 4.148 Å (fifth row), farther than that in all three starts.

| | X | N | R |
|---|---:|---:|---:|
| Ligand, mean distance of C1–C6, N1 from the crystal positions (Å) | 0.00 | 0.67 | 4.04 |
| The check | skipped | skipped | applies |
| Water oxygens within 3.5 Å of the seven positions | 1 | 1 | 1 |
| Nearest water oxygen (Å) | 3.025 | 3.025 | 2.845 |
| Water oxygens within 3.5 Å of the six ring carbons alone; nearest (Å) | 1; 3.025 | 1; 3.025 | 0; 4.148 |
| CA fit of `system.inpcrd` on `reference.pdb`, 258 atoms, r.m.s. (Å) | 0.0000 | 0.0000 | 0.0000 |

Two things follow, both for the owner of the plan (section 7, point 10):

- The count is **not** zero in X and N: a crystal water lies 3.025 Å from a ring position, beside the ring. The
  check is skipped there because of the 2.0 Å rule, not because the count is zero.
- **In start R the count is 1 as built**: one water oxygen put there by the solvation lies 2.845 Å from the
  crystal position of N1, the outer end of the ring, while no water oxygen is within 3.5 Å of any of the six
  ring carbons (nearest 4.148 Å). As the rule is worded (seven positions, count zero), R is not extended
  unless that water leaves during the 700 ps. The driver applies the rule as written; `site_atoms` in the
  protocol is the list to change if the plan is changed (for example to the six carbons).

## 2. Commands

```bash
MD=/path/to/md          # run_md.py and protocol.json
SYS=/path/to/systems    # built systems: $SYS/<model>/<start>/ with system.prmtop, system.inpcrd,
                        # selections.json, extra_forces.json, reference.pdb (and build.json)
RUNS=/path/to/runs      # runs: $RUNS/<model>/<start><replica>/
MODEL=bonded_2e         # the primary zinc set (since 6 Oct 2026, md/PRODUCTION.md); bonded_2s is the zinc
                        # sensitivity set, bonded_2e_n2fit the one for the terms at N2; "restrained": tests only
```

On the Pi these are `/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md`,
`/root/cqc-sfa-md/systems` and `/root/cqc-sfa-md/runs`, and python is `/root/cqc-sfa-tools/md/bin/python`.

`reference.pdb` is now needed (the hydration check reads it); the driver stops with exit code 2 if it is
missing. The bonded systems are built by `md/build_topology.py --model bonded --set 2e --start X` (also N, R;
`--set 2s` for the Seminario set) into `/root/cqc-sfa-md/systems/bonded_2e/<start>/` (`bonded_2s`); they were
built on 6 Oct 2026 with the commands and checksums of `md/PRODUCTION.md`, section 3 (there is no folder
`bonded`). Their
`extra_forces.json` holds an empty list. Copy the six files of each start, not the folder `build/`, and compare
with the `SHA256SUMS` the builder writes. Until they exist, the speed test and the pre-flight test can be made
with `MODEL=restrained` (same size); **the restrained model is not simulated for the paper** (plan, 6 Oct 2026,
3:02 AM, point 4).

### Order of work on the rented machine

Plan, entry of 6 Oct 2026, 5:00 AM. The driver's part is steps 1, 2, 4 and 5; a step that fails stops the work.

| # | Step | Command below | Goes on if |
|---|---|---|---|
| 1 | Speed test, 10 minutes | "Speed test" | a speed is printed; it sets the time and the cost of the rest |
| 2 | Pre-flight test of an interruption, about 5 minutes | "Pre-flight test" | both runs end with `DONE`; `IDENTICAL` or not is noted (section 7, point 2) |
| 3 | Not the driver's: quantum chemistry, charges, gate G-R, the bonded systems with gate G-Z (`md/ZINC.md`, `md/GR.md`, `md/TOPOLOGY.md`) | | the gates pass |
| 4 | Validation run X0, 5 ns from the crystal pose, and its analysis | "Validation run X0" | gate G-V `PASSED` |
| 5 | The nine runs X1 N1 R1, X2 N2 R2, X3 N3 R3 | "The campaign" | every run ends with `DONE`; then `CAMPAIGN DONE` |
| 6 | If the outcome is a conversion: the repeat of plan rule 4 with the second bonded parameter set | "The repeat of rule 4" | |
| 7 | Copy back and check (section 4), then the machine is deleted | section 4 | `checksums OK` for every run on the Pi |

Steps 1 and 2 need no gate and can run while step 3 is computed on the CPUs. A machine that was stopped and
started again needs nothing but the same command (`--resume`).

**Speed test, 10 minutes (first thing on the rented machine).**

```bash
python $MD/run_md.py --system $SYS/$MODEL/X --out $RUNS/bench --run-id bench --platform CUDA --benchmark-seconds 600
```

In an empty folder it minimises, warms up for 1 ps at 1 fs steps (not timed) and then times NPT at 4 fs; in
a folder that holds an equilibrated run it starts from that. It prints
`BENCHMARK bench: <ns/day> … 100 ns would take <h> h` and writes `benchmark.json`. No trajectory is written
during the timing, so production is slightly slower.

**Pre-flight test of an interruption on the GPU, about 5 minutes.** The same short run twice, the second
one killed and resumed; the two must agree.

```bash
P="--system $SYS/$MODEL/X --run-id preflight --platform CUDA --equil-scale 0.05 --prod-ns 0.4 --segment-ns 0.2"
python $MD/run_md.py $P --out $RUNS/preflight_a
timeout -s KILL 45 python $MD/run_md.py $P --out $RUNS/preflight_b    # choose the time so that it falls inside the run
python $MD/run_md.py $P --out $RUNS/preflight_b --resume
(cd $RUNS/preflight_a && sha256sum *.dcd state_*.xml) > a.sum
(cd $RUNS/preflight_b && sha256sum *.dcd state_*.xml) > b.sum
diff a.sum b.sum && echo IDENTICAL
```

`IDENTICAL` means a repeated segment on this GPU is the segment that was lost. If the files differ, the
campaign is still sound (no frame is lost or doubled, section 1); the manifest then simply records runs
that cannot be regenerated bit for bit. Both folders carry a file `TEST_ONLY` and are not production data.
With start X the hydration check is skipped, so this test makes no extension step; those are ordinary steps
and were interrupted and resumed on the Pi (section 6.2).

**One trajectory.**

```bash
python $MD/run_md.py --system $SYS/$MODEL/X --out $RUNS/$MODEL/X1 --run-id X1 \
    --platform CUDA --stage all --prod-ns 100 --segment-ns 5 --resume
```

`--resume` on an empty folder starts the run, so this one command starts, continues and confirms a run.
The last line on success is `DONE X1`. The stages can be run one at a time with `--stage min`,
`--stage equil --resume`, `--stage prod --resume`. On the Pi: `--platform CPU --threads 4`, with
`OPENMM_CPU_THREADS=4` in the environment (section 7, point 2: `--threads` alone does not bind every thread).

**Validation run X0 (stage S4), 5 ns from the crystal pose, and gate G-V.**

```bash
python $MD/run_md.py --system $SYS/$MODEL/X --out $RUNS/$MODEL/X0 --run-id X0 \
    --platform CUDA --stage all --prod-ns 5 --segment-ns 5 --resume
python $MD/analysis/analyse_md.py $RUNS/$MODEL/X0        # needs MDAnalysis and Matplotlib; on the Pi after copying back, or there
```

Run X0 is the validation run of its zinc model (plan, 5 Oct 2026, 9:41 PM, point 4; `md/analysis/SPEC.md`,
section 1): the analysis knows it by its name, plans 5 ns for it without an option, gives it alone a verdict of
gate G-V (section 11 of the specification; `gate_GV.verdict` in its `summary.json`, `PASSED` or `FAILED`) and
leaves it out of the campaign table. **The nine runs start only after `PASSED`.** A failed gate stops the work
and is reported (plan, section 4). X0 is not one of the nine runs and is not lengthened into one: X1 starts
again from the built structure with its own seed. To the driver X0 is a run like any other (2 ns of
equilibration, one segment of 5 ns). A second zinc model has its own X0 (`$RUNS/<its model>/X0`).
(A run can be lengthened: a larger `--prod-ns` with `--resume` adds segments; tests 9 and 9b. Not for X0.)

**The campaign: nine runs, one after another, on one GPU.** Save as `campaign.sh` beside the runs; it can
be started again at any time (a finished run answers `DONE` in a second, an interrupted one continues).

```bash
#!/usr/bin/env bash
# Nine runs, one after another, on one GPU. Safe to start again at any time: a finished run
# answers "DONE" in a second, an interrupted one continues after its last complete step.
MD=${MD:?folder with run_md.py and protocol.json}
SYS=${SYS:?folder with the built systems}
RUNS=${RUNS:?folder for the runs}
MODEL=${MODEL:-bonded_2e}
PY=${PY:-python}
OPTS=${OPTS:---platform CUDA --prod-ns 100 --segment-ns 5}
mkdir -p "$RUNS/$MODEL"
for rep in 1 2 3; do
  for start in X N R; do
    id=$start$rep
    $PY "$MD/run_md.py" --system "$SYS/$MODEL/$start" --out "$RUNS/$MODEL/$id" --run-id "$id" \
        $OPTS --resume >> "$RUNS/$MODEL/$id.console" 2>&1
    code=$?
    if [ $code -ne 0 ]; then
      echo "$id stopped with exit code $code; see $RUNS/$MODEL/$id.console"
      exit $code
    fi
    tail -n 1 "$RUNS/$MODEL/$id.console"
  done
done
echo "CAMPAIGN DONE"
```

```bash
MD=… SYS=… RUNS=… nohup bash campaign.sh > campaign.out 2>&1 &
```

**On the rented machine (6 Oct 2026) the nine runs are not started by this loop** but as tasks of the queue, one
waiting for the other in the same order, with the same command for each run
(`md/gpu/tasks_production.json`; `md/gpu/GPU_RUN.md`, section 5; `md/PRODUCTION.md`, section 6). There are then no
`*.console` files: the driver's lines are in each run's `run_md.log` and in `/root/cqc-sfa-md/logs/<task>.log`.

The order is X1 N1 R1, X2 N2 R2, X3 N3 R3: if the machine or the budget ends early, each start has the
same number of finished runs. The loop stops at the first run that fails and names its console file.
The replicas start at 1; the loop never touches X0. An R run can take up to 2 ns of equilibration more than
the others (section 1A); `grep "hydration check" $RUNS/$MODEL/R1.console` shows what was counted.

**The repeat of rule 4** (only if the outcome is a conversion; plan, 6 Oct 2026, 3:02 AM, point 3): the rotated
start once more with the second bonded parameter set, built as a model of its own. The same command with the
system folder and the run folder of that model, for example

```bash
python $MD/run_md.py --system $SYS/<second model>/R --out $RUNS/<second model>/R1 --run-id R1 \
    --platform CUDA --stage all --prod-ns 100 --segment-ns 5 --resume
```

**Fixed on 6 Oct 2026** (plan, entry of 9:15 PM, point 3; `md/PRODUCTION.md`): the rotated start is run once
with each of two other parameter sets whatever the outcome, after the nine runs: with the other zinc set
(`bonded_2s`, run folder `runs/bonded_2s/R1`) and with the candidate terms at N2 that gate G-R2 did not accept
(`bonded_2e_n2fit`, run folder `runs/bonded_2e_n2fit/R1`); both are in the queue (`md/gpu/tasks_sensitivity.json`).
Before that date: the name of that model folder is not fixed anywhere yet (section 7, point 11); the analysis takes the model
from the name of the folder above the run.

**Resuming.** Run the same command (or the loop) again. What the driver does then:

- verifies that the inputs, the protocol, the run id, the seed, `--equil-scale` and `--segment-ns` are those
  of the run in the folder (otherwise exit code 3, and the run is left as it is);
- reads the chain of steps from the records (with the extension steps that the counts asked for), finds the
  last complete step, verifies the checksum of its saved state, removes the files of the interrupted step and
  continues;
- refuses to work in a folder in which another `run_md.py` is alive (`run_md.lock`).

One exception to the first point, made for runs that began before 6 Oct 2026, 5:28 AM IST: if the protocol
file has gained the block `equilibration.hydration_check` and is otherwise the protocol recorded in the
manifest of the run, the run keeps its recorded protocol and continues without the check; the console, the
invocation and `hydration_check.status` of the manifest say so. Any other change of the protocol is refused
as before.

To continue on **another machine**, copy the run folder (at least `manifest.json`, all `rec_*.json` and the
last `state_*.xml`) and the system folder; finished trajectory files may stay in the archive (test 9b).
A SIGTERM or Ctrl-C ends the run within 500 steps with exit code 143 or 130; a SIGKILL or a lost machine
needs nothing special.

Exit codes: 0 done; 2 wrong command line or input files; 3 the folder does not fit the request; 4 the
simulation failed (no usable platform, energy not finite); 130/143 stopped by a signal; 1 anything else.

## 3. What each output file is

In `$RUNS/<model>/<run>/`; `<step>` is `min`, `eq01`…`eq06`, `eq03x1`…`eq03x4` (only if the hydration check
extended the held part) or `seg0001`…

| File | Content |
|---|---|
| `manifest.json` | the whole record: run id, seed and seed rules, `test_only`, sha256 of every input and of `run_md.py`, the protocol in full, the system (atoms, masses, constraints, forces, box), the extra forces as applied, every invocation (command line, start and end in UTC and IST, exit code, host, versions of OpenMM, Python and NumPy, platform with its properties, GPU name, `nvidia-smi` line), the hydration check (`hydration_check`, section 1A), and the records of all steps with ns per day for each (`equilibration.stages` holds the records of the stages and of the extension steps in the order they were run) |
| `rec_<step>.json` | record of one step, written last: settings, seeds, start and end time, speed, final and mean temperature and density, sha256 and size of the step's files, sha256 of the state it started from, for `eq03` and the extension steps the hydration check, and a few distances for a quick look (Zn–N, the two orientation distances, water near the ligand, water at the crystal ring site; never used to steer the run, the hydration check apart) |
| `state_<step>.xml` | saved state at the end of the step (OpenMM XML; nm, nm ps⁻¹; coordinates as integrated, not wrapped). What `--resume` starts from |
| `checkpoint_seg<NNNN>.chk` | OpenMM binary checkpoint at the end of the newest segment only; tied to the machine; not needed by `--resume` |
| `log_<step>.csv` | state log every 10 ps: `step, time_ps, potential_energy_kj_per_mol, kinetic_energy_kj_per_mol, temperature_K, volume_nm3, density_g_per_cm3, position_restraint_energy_kj_per_mol, extra_bond_energy_kj_per_mol, thermostat_temperature_K, speed_ns_per_day` |
| `solute_<step>.dcd` | the atoms of `selections.json` `solute` (protein, zinc, ligand, in that order) every 10 ps; Å; as integrated, never wrapped, so the complex is whole in every frame |
| `full_<step>.dcd` | all atoms every 100 ps; Å; the solute exactly as in the solute file, every water and ion moved by cell vectors to the image nearest to the centre of the solute |
| `solute.pdb`, `full.pdb` | topology for the two kinds of trajectory, with the atom and residue names of the prmtop (HID, HIE, WAT, …), the start coordinates, serial number = prmtop atom index + 1, residue number = prmtop residue index + 1 (`md/system/residue_map.csv` leads back to the deposited numbers) |
| `final_solute.pdb` | the solute at the end of the last segment (for git) |
| `benchmark.json` | results of `--benchmark-seconds` |
| `run_md.log` | the messages of all invocations, with IST times |
| `TEST_ONLY` | present only if the run is not a production run (`--equil-scale` below 1, or a protocol whose purpose is not "production"); the manifest says the same in `test_only` |
| `run_md.lock` | present while a process works in the folder |

Time: production time is zero at the end of the equilibration; frame *k* of a production run (counted
through the segments) is at *k* × 10 ps, so a segment of 5 ns holds 500 solute frames (the last at the
segment's end) and 50 full frames. Equilibration files count from the start of `eq01`, through the extension
steps if there are any. The equilibration is part of the record: with the ligand free from `eq04` on, a
reorientation there would be seen in `solute_eq04.dcd` to `solute_eq06.dcd`.

**Size**, for the built systems (40,902 atoms in X, 40,911 in N, 40,857 in R; 4,094 solute atoms), estimated
from the file formats and from the sizes of the files the driver wrote for system X in the pipeline test on the
Pi, not measured on a production run [`sizes.json`; GB = 10⁹ bytes]:

| | Size |
|---|---:|
| One solute frame; one full frame | 49,208 bytes; 490,904 bytes |
| One saved state | 6.7 MB |
| A run of 100 ns: solute frames, full frames, saved states of the 20 segments, equilibration with minimisation and the rest | 0.49 + 0.49 + 0.13 + 0.07 = 1.19 GB |
| One extension step of 500 ps; four | 0.012 GB; 0.046 GB |
| The validation run X0 (5 ns) | 0.12 GB |
| The nine runs | 10.7 GB |
| The nine runs, X0 and one repeat of 100 ns, four extension steps in every R run | 12.2 GB |

The full frames, 0.31 GB per run in the estimate for 25,600 atoms, are now as large as the solute frames. The
plan says "a few GB"; it is about 12 GB. On the rented machine add one checkpoint file per run and the folders
of the speed test and the pre-flight test; 30 GB of free disk is ample.

## 4. What to copy back, and the check

Everything in `$RUNS/<model>/` (the run folders, X0 among them, and the `*.console` files) except `*.chk`,
`*.tmp` and `run_md.lock`; `campaign.out`; of the speed test `bench/benchmark.json`, `bench/manifest.json` and
`bench/run_md.log`; of the pre-flight test `a.sum`, `b.sum` and the two `manifest.json` (its trajectories are
test data and need not come back); the output folder of the analysis of X0 if it was made there. Files of
finished steps never change, so pulling while the campaign runs is safe, for example hourly from the Pi:

```bash
rsync -a --exclude '*.chk' --exclude '*.tmp' --exclude run_md.lock <gpu>:<RUNS>/ /root/cqc-sfa-md/runs/
```

Of this, the trajectories go to the Pi archive and R2, not into git (plan, section 6). `final_solute.pdb` of
every run is the final frame the plan wants in git; `manifest.json`, the records, the state logs and
`run_md.log` are small and can go with it; the saved states are 6.7 MB each and belong with the
trajectories.

Check on the Pi, against the checksums the driver recorded while it wrote the files:

```bash
/root/cqc-sfa-tools/md/bin/python - /root/cqc-sfa-md/runs/bonded_2e/* <<'EOF'
import hashlib, json, os, sys
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()
for d in sys.argv[1:]:
    if not os.path.isfile(os.path.join(d, "manifest.json")):
        continue
    m = json.load(open(os.path.join(d, "manifest.json")))
    steps = [s for s in [m["minimisation"]] + m["equilibration"]["stages"] + m["production"]["segments"] if s]
    missing, changed = [], []
    for s in steps:
        for name, info in s["files"].items():
            p = os.path.join(d, name)
            if not os.path.isfile(p):
                missing.append(name)
            elif sha256(p) != info["sha256"]:
                changed.append(name)
    h = m.get("hydration_check", {})
    print(m["run_id"], m["status"], m["production"]["completed_ns"], "ns", "TEST" if m["test_only"] else "production",
          "checksums OK" if not missing + changed else f"MISSING {len(missing)} CHANGED {len(changed)}: " + " ".join(missing + changed),
          "| hydration check:", h.get("final_decision") or h.get("status"), "| extension steps:", h.get("n_extensions", 0))
EOF
```

A complete production run prints, for example,
`X1 complete 100.0 ns production checksums OK | hydration check: go on: check skipped: the ligand itself is at the site | extension steps: 0`.

## 5. The protocol and the reasons

All values are in `md/protocol.json`; the manifest of every run carries a copy. Harmonic terms are
U = ½ k d². No source is cited that was not read: where a default is named, it was read in the software
installed here (OpenMM 8.6.1, ParmEd 4.3.1).

| Setting | Value | Reason |
|---|---|---|
| Force field | as in `system.prmtop`; the driver adds only the harmonic bonds of `extra_forces.json`, the restraints and the barostat | the driver simulates what the builder built; the choice of force field is the plan's |
| Hydrogen mass | 3.024 amu for every hydrogen bonded to a heavy atom, taken from that atom; water hydrogens stay at 1.008 amu; applied when the System is created | three hydrogen masses, the default of ParmEd's `HMassRepartition`; it slows the motions of hydrogens that remain once their bonds are fixed, which is what makes a 4 fs step possible; the total mass is unchanged (test 10) |
| Time step | 4 fs in production and in the last three equilibration stages; 2 fs in the first three and in the extension steps | 4 fs halves the GPU time of 900 ns; 2 fs while the system heats and finds its density, where a long step fails first |
| Constraints | bonds to hydrogen; rigid water; tolerance 10⁻⁶ | removes the fastest vibrations; TIP3P is a rigid model; ten times tighter than OpenMM's default (10⁻⁵), because the error of the constraints is paid at every 4 fs step |
| Electrostatics | PME, real-space cut-off 1.0 nm, Ewald error tolerance 5 × 10⁻⁴ | the cut-off is the project's; the tolerance is OpenMM's default. Measured on the stand-in: the PME forces differ from a 10⁻⁶ reference by 0.36 kJ mol⁻¹ nm⁻¹ r.m.s., 0.06 % of the r.m.s. force; 10⁻⁴ would give 0.013 % with a grid 3.5 times larger (`ewald_tolerance_check.json`) |
| Lennard-Jones | cut at 1.0 nm, no switching function, long-range dispersion correction on | no switch is the default for an Amber topology in OpenMM, and a switch would alter the force field as built; without the correction the attraction beyond the cut-off is missing from the pressure and the density in NPT comes out too low |
| Thermostat | LangevinMiddle integrator, 300 K, friction 1 ps⁻¹ | 300 K is the plan's; a stochastic thermostat makes the three runs of a start independent through their seeds; 1 ps⁻¹ equilibrates within picoseconds and is far from the nanosecond motions in question |
| Barostat | Monte Carlo, isotropic, 1 bar, an attempt every 25 steps | 1 bar is the plan's; isotropic for a protein in water; 25 steps is OpenMM's default |
| Centre-of-mass motion | removed every step | default of the topology reader; stops the drift of the whole box that a stochastic thermostat allows |
| Minimisation | L-BFGS to 10 kJ mol⁻¹ nm⁻¹ (r.m.s. force), at most 5,000 iterations per pass; restraints of 1,000 kJ mol⁻¹ nm⁻² on the heavy atoms of protein, zinc and ligand and on the oxygens of the crystal water (2,289 atoms in X and N, 2,285 in R) | relaxes what tleap leaves (hydrogens of water and hydroxyl groups placed without regard to neighbours, `SYSTEM.md` section 4) while the crystal and docked heavy atoms stay; the tolerance is OpenMM's default; the limit bounds the time and is per pass (section 7, point 8) |
| Restraint reference | the coordinates of `system.inpcrd` | the crystal or docked positions, not the minimised ones, so that the minimiser cannot move the reference |
| Stage 1 `heat` | NVT, 100 ps, 2 fs, 50 → 300 K (raised once per ps), k = 1,000 on the heavy atoms of protein, zinc and ligand | starts cold, so the energy left by the minimisation is taken up gradually; constant volume while the temperature is not yet the barostat's; k = 1,000 allows 0.50 Å r.m.s. per coordinate at 300 K, enough for bonds and angles to relax, too little to lose the pose |
| Stage 2 `density` | NPT, 300 ps, 2 fs, 300 K, same restraints | the box reaches its density and the water arranges itself around a solute that is held; this is the time for water to enter the space left empty in the rotated start (`SYSTEM.md` section 9, point 3) before anything is released |
| Stage 3 `loosen` | NPT, 300 ps, 2 fs, k = 250, same atoms | 1.0 Å r.m.s. per coordinate: side chains and ligand adjust to the water without leaving their places |
| Hydration check | at the end of stage 3: water oxygens within 2.5 Å of the crystal positions of the six ring carbons C1 to C6; applies if the ligand's own ring carbons are on average more than 2.0 Å from them; count zero: extension steps of 500 ps with the settings of stage 3, at most four; then the run goes on (section 1A) | the plan's rule (6 Oct 2026, 1:03 AM, point 3; 3:02 AM, point 8; changed at 8:01 AM, point 7, from seven positions with N1 to the six ring carbons, and at 8:56 AM, point 6, from 3.5 Å to 2.5 Å). The 2.0 Å separate the starts: as built, on the seven positions of the first wording, 0.00, 0.67 and 4.04 Å, with room for a ring turned by 180° in place (1.4 Å by the names of its atoms on the seven positions; 1.6 Å on the six ring carbons, that is 1.378 Å × 7 / 6, since C1, C4 and N1 keep their names in that calculation; `unit_checks.json` and `scripts/unit_checks.py`); at the end of stage 3, on the six ring carbons, 0.23 to 0.63 Å in the five runs from X and N up to 7 Oct 2026, 12:15 AM IST, and 3.94 Å in the one run from R (rehearsal runs TX0, TN0 and TR0, `md/test/REHEARSAL_2e.md`; X0 of both zinc sets, `md/PRODUCTION.md`, section 5.1; X1, its manifest) |
| Stage 4 `release_side_chains_zinc_ligand` | NPT, 400 ps, 4 fs, k = 250 on the backbone (N, CA, C, O) only | ligand and zinc are released together with the side chains, as fixed for the project, so the site relaxes as a whole while the fold is held; first stage at 4 fs |
| Stage 5 `loosen_backbone` | NPT, 400 ps, 4 fs, k = 50 on the backbone | 2.2 Å r.m.s. per coordinate: nearly free, so the last release is not a jump |
| Stage 6 `free` | NPT, 500 ps, 4 fs, no restraint | the conditions of the production; its end is time zero. Total 2.0 ns, up to 4.0 ns with extension steps |
| Crystal water | its oxygens are restrained in the minimisation only (`crystal_water_O` of `selections.json`: 215 in X and N, 211 in R) | fixed for the project; the builder lists them since 6 Oct 2026 (`md/TOPOLOGY.md`, section 3) |
| Production | NPT, 300 K, 1 bar, 4 fs, 100 ns, no restraint, no bias of any kind | the plan's; nothing acts on the orientation of the ring |
| Segments | 5 ns | an interruption costs at most one segment; 20 saved states per run |
| Output | solute every 10 ps, all atoms every 100 ps, state log every 10 ps | the shared conventions: 10,000 solute frames per run for the two orientation distances, 1,000 full frames for the water in the site |
| GPU precision | mixed | OpenMM's middle setting: forces in single precision, integration in double; not exercised in the tests of section 6 (no GPU on the Pi); in force in the runs on the rented machine (`minimisation.platform` in the manifest of X1: Precision mixed) |
| Deterministic forces | switched on where the platform offers it | so that a repeated segment can be the lost one |
| Seeds | run seed from the run id; step seeds from run seed, step name and purpose (velocities, integrator, barostat) | a run can be repeated from its name alone; no two steps share a random sequence |

The restraint constants in the Amber convention (U = k d², kcal mol⁻¹ Å⁻²) are 1.195, 0.299 and 0.060.

## 6. Tests on the Raspberry Pi

### 6.1 The first tests (5 Oct 2026)

Made with the driver as it stood before the hydration check (sha256 `45e4cf28…f3bc29`). They were not repeated
on 6 Oct 2026; that the driver of today makes the same files when the protocol holds no hydration block, and
continues a run of the old driver, is shown in section 6.2 (groups "without the block" and "run of the old
driver").

Stand-in system, **not hCA II**: a capped pentapeptide, one Zn²⁺, the sulfonamidate anion, 2 Na⁺, 2 Cl⁻
and 1,495 TIP3P waters in a truncated octahedron (4,600 atoms), with one harmonic bond in
`extra_forces.json` (Zn to NE2 of the histidine, r₀ = 3.0 Å, k = 200,000 kJ mol⁻¹ nm⁻²). CPU platform, at
most 2 threads, on a machine shared with other jobs. Test runs use `protocol_test.json` (the protocol
with output every 0.1 and 0.2 ps and at most 500 minimiser iterations per pass), `--equil-scale 0.002`
(six stages, 4 ps in all) and two production segments of 1 ps; every one of them is marked `TEST_ONLY`.
Everything is in `/root/cqc-sfa-md/test_driver/`: the scripts in `scripts/` (`build_standin.py`,
`make_test_protocols.py`, `test_driver.py run | run1 | check`, `campaign.sh`, `campaign_test.sh`,
`preflight_test.sh`, `ewald_check.py`, `long_check.py`), the runs in `runs/` and `campaign/`, the numbers below in
`test_results.json`, `long_check.json` and `ewald_tolerance_check.json`. (`make_test_protocols.py` there reads
`md/protocol.json`; run again today it would write test protocols with the hydration block.)

| # | Test | Result |
|---|---|---|
| 0 | Does the CPU platform repeat itself? (single-point forces in three new Contexts; 50 steps twice) | 1 thread with DeterministicForces: yes, differences 0.0. 2 threads with DeterministicForces: **no**, forces differ by up to 8.5 × 10⁻⁴ kJ mol⁻¹ nm⁻¹ and positions by 1.1 × 10⁻⁶ nm after 50 steps; two uninterrupted two-thread minimisations end 0.094 nm apart at most |
| 1 | Minimisation lowers the energy | **passed**: 17,574.5 → −73,957.3 kJ mol⁻¹ (1,601 L-BFGS iterations in four passes, three of them stopped by the test's limit of 500); restrained heavy atoms moved 0.57 Å r.m.s. With the production protocol (limit 5,000 per pass): → −76,884.7 kJ mol⁻¹, passes of 991, 3,356, 10 and 6 iterations, none stopped by the limit (`runs/realproto`) |
| 2 | Shortened equilibration and two segments run | **passed**: stages of 100, 300, 300 steps at 2 fs and 200, 200, 250 steps at 4 fs; temperature 94 → 291 K, density 0.834 → 0.957 g cm⁻³ over the 4 ps; two segments of 250 steps with 10 solute frames, 5 full frames and 10 log rows each; mean temperature 296.2 and 301.5 K at 4 fs |
| 3a | Interrupted and resumed = uninterrupted, **1 thread** | **passed**. The run was made in pieces, stopped by SIGTERM in the middle of `eq03` (exit code 143; 3 solute frames, 1 full frame and 3 log rows left without a record) and killed by SIGKILL in the middle of `seg0002` (4 solute frames, 2 full frames, 4 log rows left), each time continued with `--resume`. All 36 files (16 trajectories, 9 saved states, 3 PDB files by sha256; 8 state logs without the speed column) are identical to those of the uninterrupted run with the same seed; largest coordinate difference 0.0 Å; final positions, velocities and box differ by 0.0 |
| 3b | The same with **2 threads** | not bit for bit (test 0), so only the structure is compared, and that **passed**: same numbers of frames, frame times and log rows; every trajectory file holds whole frames only; every step started from the recorded state of the step before. The coordinates themselves differ (minimised positions by up to 0.15 nm), as two uninterrupted two-thread runs do |
| 4 | The extra harmonic bond acts | **passed**: Zn···NE2 starts at 5.00 Å; 2.99 Å after the minimisation; 2.98 ± 0.03 Å over the 40 equilibration frames; 2.99 ± 0.04 Å (2.92–3.05) over the 20 production frames; expected thermal width 0.035 Å. Control without the bond: 5.95 ± 0.18 Å (5.65–6.26) in production. The energy column of the state logs equals ½ k (r − r₀)² from the frames within 1.2 × 10⁻⁴ kJ mol⁻¹ (60 rows) |
| 5 | The manifest is complete | **passed**: no key missing (inputs, script, protocol in full, seeds, system, extra forces, invocations with versions, platform and properties, times, speed per step); recorded checksums of inputs and outputs match the files |
| 6 | Trajectory files | **passed**, read with MDAnalysis' DCD reader: frame times 1.1 … 2.0 ps for the second segment; the solute atoms of the full frames equal the solute frames exactly; last solute frame = saved state within 9.5 × 10⁻⁷ Å; every water oxygen is at its minimum-image distance from the centre of the solute (MDAnalysis, from the cell stored in the frame; largest difference 4.0 × 10⁻⁶ Å); water molecules intact |
| 7 | Exit codes | **passed**, 13 calls: folder in use without `--resume` 3; other seed 3; other system 3; `--stage prod` before the earlier stages 3; altered saved state 3; second process in a folder in use 3; system folder missing 2; `--prod-ns` not a multiple of `--segment-ns` 2; `--platform CUDA` on this machine 4; `--resume` on a complete run 0 with nothing changed; three continuations 0 |
| 8 | `--benchmark-seconds` | **passed**: 10.05 ns per day from the equilibrated state and 10.05 from the minimised state (4,600 atoms, 2 threads, busy machine; 292 steps in 10.0 s each) |
| 9 | Continuation to a longer run | **passed**: a copy of the finished run continued to three segments; earlier files untouched; only the newest checkpoint kept |
| 9b | Continuation after the finished files were moved away | **passed** (1 thread): with only `manifest.json`, the records and the last saved state in the folder, the third segment is identical to that of the full copy |
| 10 | Hydrogen masses | **passed**: total mass 27,988.132 amu in the prmtop and in the System; all 49 hydrogens outside water at 3.024 amu, water hydrogens at 1.008 amu, lightest heavy atom 5.962 amu; a saved state read and written again is the same text |
| 11 | Cost of the output at production size | 0.013 s per full frame of 25,600 atoms (wrapping and writing, on the Pi); measured again on 6 Oct 2026 for the built system X: section 6.2, group "cost" |
| 12 | Twenty-two segments: a copy of the run continued by 20 ps at 4 fs | **passed**: 22 of 22 segments, 220 solute frames and 220 log rows, times rising by 0.1 ps without a gap, segment headers contiguous, no frame repeated across a boundary; temperature 300.8 ± 4.5 K (4.4 K expected for 9,263 degrees of freedom); density from 0.986 to 1.003 g cm⁻³; the bond at 2.994 ± 0.034 Å over all 220 frames (`long_check.json`) |
| – | The launcher loop of section 2, nine tiny runs | **passed**: started, the loop and the running `run_md.py` killed by PID during N1, started again: X1 answered `DONE`, N1 continued with `eq03`, all nine ended `DONE` and the loop `CAMPAIGN DONE`; the check of section 4 prints `checksums OK` nine times (`campaign/`) |
| – | The pre-flight test of section 2, with the stand-in on one thread | **passed**: the second run was killed by `timeout -s KILL` during `eq05`, resumed, and `diff` of the 25 checksums (16 trajectories, 9 saved states) printed `IDENTICAL` (`logs/preflight.out`) |
| – | The production protocol file itself | `md/protocol.json` (as it then was) with `--equil-scale 0.002`, minimisation and six stages: ran, exit code 0, marked `TEST_ONLY` because of the scale; `--platform auto` fell back to the CPU platform with a warning (`runs/realproto`) |

Wall time of the final round: 485 s (two-thread runs), 442 s (one-thread runs), 255 s (launcher loop),
175 s (production protocol file), 162 s (twenty-two segments), 71 s (pre-flight test).

### 6.2 The hydration check, run X0 and the driver against its earlier self (6 Oct 2026)

Driver sha256 `b0f3e2b2…a74fcc`, protocol sha256 `d886fe5b…d60d39` (the files of this folder as they stood when
the tests were made, with the count of the hydration check on seven positions within 3.5 Å; the campaign's two
files are named in the state paragraph at the top). CPU platform,
**one thread** (`--threads 1` with `OPENMM_CPU_THREADS=1` and `OMP_NUM_THREADS=1`), nice 10, on the shared machine.
Everything is in `/root/cqc-sfa-md/test_driver/hydration/`: the scripts in `scripts/` (`build_standin2.py`,
`make_test_protocols.py`, `run_tests.py <group>`, `check_tests.py`, `unit_checks.py`, `sizes.py`,
`pi_run_numbers.py`, `check_copy.py`), the systems in `systems/`, the runs in `runs/`, the
results in `test_results.json` (81 tests in 9 groups), `unit_checks.json` (9 checks), `check_copy.out`,
`sizes.json`, `pi_run_numbers.json` and `systems.json`. A second call of `build_standin2.py` gave the same system files
(sha256).

**Second stand-in, not hCA II**: a capped ten-residue helix (ACE-ALA-HID-ALA-THR-GLU-ALA-PHE-ASN-ALA-ALA-NME),
one Zn²⁺ tethered as in the first stand-in, the sulfonamidate anion, 2 Na⁺, 2 Cl⁻ and 985 TIP3P waters in a
truncated octahedron (3,120 atoms, 10 CA atoms). It was built because the five CA atoms of the first stand-in
lie almost on a line (singular values of their positions 11.52, 1.19 and 0.00 Å; 14.90, 4.90 and 4.69 Å for
the helix): the fit that carries the ring site is loose about that line, and the ligand, held in place, was
2.70 Å from its own fitted site at the end of `eq03` of the old run `ref1` [`unit_checks.json`]. Four systems
share its topology and coordinates and differ in `reference.pdb` and `selections.json` only:

| System | What is changed | Meant to give |
|---|---|---|
| `skip` | nothing: the ring site of `reference.pdb` is where the ligand is built | the skipped case |
| `wet` | the ring site moved by 9 Å into the water, 6.1 Å clear of every solute atom | the check applies, water at the site |
| `dry` | as `wet`, and `water_O` lists only the 623 of 985 water oxygens that are at least 14 Å from the moved site as built (a device of the test: a stand-in without a pocket has water everywhere) | count zero at every count |
| `late` | as `dry` with a limit of 5 Å (951 listed) | count zero at first, water a little later |

Test runs use `protocol_fast_hyd.json` (the protocol of this folder with purpose "test", output every 0.1 and
0.2 ps, at most 30 minimiser iterations per pass), `--equil-scale 0.0005` (stages of 25, 75, 75 steps at 2 fs
and 50, 50, 62 steps at 4 fs; an extension step is 125 steps, 0.25 ps), one production segment of 100 steps,
run id `T1`. Every run is marked `TEST_ONLY`. They are short because the machine was busy: 2.9 steps per second
in the median over the 85 MD steps made for these tests (1.7 to 15.5); the nine groups took 62 minutes of wall
time.

| Group | Test | Result |
|---|---|---|
| cases | The skipped case (`skip`) | **passed**: the ligand is 1.48 Å from the site on average (limit 2.0 Å), so the check is recorded as skipped, with its count (3) and `go on: check skipped: the ligand itself is at the site`; no extension step; chain `min`, `eq01`–`eq06`, `seg0001` |
| cases | A case that does not extend (`wet`) | **passed**: the check applies (ligand 10.05 Å from the site), 15 water oxygens within 3.5 Å (nearest 0.88 Å), `eq04` is next |
| cases | Count zero that extends (`dry`) | **passed**: count 0 at the end of `eq03` and of each of the four extension steps (nearest listed water oxygen 11.0 to 11.7 Å); the chain is `eq03`, `eq03x1`, `eq03x2`, `eq03x3`, `eq03x4`, `eq04`, …; the four steps have the settings of `eq03` and 125 steps each; then the run goes on to `DONE`. Manifest: `n_extensions` 4, 6 stages and 4 extension steps, `final_decision` "go on: still no water oxygen at the site after 4 extension steps, the most allowed: the run goes on as built, with an empty site", `site_without_water_when_the_solute_was_released` true |
| arrive | Water found after some extension steps | **passed**, with a second device: the `dry` system and a test protocol whose cut-off (11.2 Å) lies between the nearest listed water of the first three counts of `dry` and that of its fourth (11.16 Å). Counts 0, 0, 0, 1: three extension steps, `eq04` starts from `state_eq03x3.xml`; the steps up to `eq03x3` are those of `dry`, file by file (25 files) |
| cases | `late`, as found | the device did not do what it was meant to: one listed water oxygen was already within 3.5 Å (3.23 Å) at the first count, so the run went on without an extension step, like `wet`. The group "arrive" covers the branch it was meant for |
| cases, arrive | The extension steps are ordinary steps | **passed** for all five runs: every step has its own seeds (14 per run, 22 in `dry`, 20 in the run of "arrive", all distinct, equal to the rule of the manifest applied to the step name); every step started from the recorded state of the step before; the equilibration time has no gap; frames and log rows as the length of each step asks; sha256 of every file equals the manifest; the check of section 4 prints `checksums OK` with the decision and the number of extension steps for each of them (`check_copy.out`) |
| cases, arrive | The count, made again | **passed**: for the 12 counts of the five runs a second calculation (SciPy's `Rotation.align_vectors` for the fit, MDAnalysis' `distance_array` with the cell) on the saved states gives the same number and the same waters; nearest distance and ligand distance within 0.0004 Å |
| cases | The check does not touch the dynamics | **passed**: all 32 files of `wet` (7 solute and 7 full trajectories, 8 saved states, 3 PDB files by sha256; 7 state logs without the speed column) equal those of `skip`; `dry` equals `skip` up to `eq03` |
| interrupt | Interrupted and resumed across extension steps = uninterrupted | **passed**. The `dry` run again in pieces: killed (SIGKILL) at the moment `rec_eq03.json` existed (the count made, `eq03x1` next, nothing of it written yet); stopped by SIGTERM inside `eq03x2` (exit code 143; 1 solute frame, a full file with its header only and 1 log row left without a record); killed inside `eq03x4` (1 solute frame, an empty full file, 1 log row left); each time continued with `--resume`. All 48 files (11 + 11 trajectories, 12 saved states, 3 PDB files by sha256; 11 state logs without the speed column) are identical to those of the uninterrupted `dry` run; the five check records are equal; the manifest lists the five invocations (two without an end, one with exit code 143) |
| run id | Run X0 is a run like any other | **passed**: runs `X0` and `X1` of the `skip` system with the same `--seed` give the same 8 saved states (sha256), the same 7 state logs and the same frames in all 14 DCD files (their headers carry the run id); the two manifests have the same keys throughout, the same steps and the same hydration check. In the source no line reads a start or a replica number; the run id enters the seed rule, the manifest, the titles and the messages (18 lines). Without `--seed` the rule gives 1793718340 for X0 and 154924513 for X1 |
| without the block | The driver of today = the driver before the change | **passed**: with a protocol without the hydration block (`protocol_fast_noblock.json`) the two drivers (sha256 `45e4cf28…` and `b0f3e2b2…`) write the same 32 files; no check is made, no record names one, and the manifest says "not made: the protocol has no block equilibration.hydration_check". With the block (`skip`) the same 32 files again; the record of `eq03` differs by the key `hydration_check` only |
| run of the old driver | A run begun by the old driver, continued with the protocol file of today | **passed**: the old driver, protocol without the block, killed inside `seg0001`; the driver of today with the protocol file that holds the block continues it (exit code 0, a line `NOTE: … the run keeps the protocol of its manifest (no hydration check)`, the note also in the invocation and in `hydration_check.status`); all 32 files equal the uninterrupted run that the old driver made today |
| compat_ref1 | The same on run `ref1` of 5 Oct 2026 | the continuation itself works as above (its three tests passed), but **2 tests failed**: the repeated `seg0002` is not the `seg0002` of 5 Oct. 5 of 36 files differ (the two trajectories, the log, the saved state, `final_solute.pdb`), from the first log row on (step 275: −64,527.2982 against −64,527.2906 kJ mol⁻¹, 1.2 × 10⁻⁷ of the value). The old driver, repeating the segment today, differs from 5 Oct in the same way, and its 36 files equal those of the driver of today. So the difference is between the two days, not between the two drivers. Suspected, **not verified** (it needs more than one thread): the first tests did not set `OPENMM_CPU_THREADS`, and OpenMM's PME library reads that variable, not the platform's `Threads`; today's runs set it to 1. See section 7, point 2 |
| errors | Refusals | **passed**: the protocol asks for the check and the system has no `reference.pdb`: exit code 2; `--resume` with another cut-off in the block: 3; a run of the old protocol with a file that has the block **and** another friction: 3; a protocol with the block and a record of `eq03` without a check (a run of the old driver): 3, "Start this run again in an empty folder"; a complete `eq03x1` in a run whose `eq03` says go on: 3, "the folder … is damaged" |
| errors | No-ops and continuation | **passed**: files of `eq03x1` without a record in such a run are removed and the run answers `DONE`; `--resume` on the complete `dry` run changes no file; the `dry` run continued to two segments keeps its four extension steps and every earlier file (`final_solute.pdb` is written anew) |
| bench | `--benchmark-seconds` in a folder with extension steps | **passed**: in a copy of the `dry` run (manifest, records, `state_eq06.xml`) the speed test starts from `state_eq06.xml`, times 5 s and prints its line (0.66 ns per day for 3,120 atoms on one busy thread; no warm-up needed). The speed test in an empty folder (minimisation, then 1 ps of warm-up) was not run again |
| unit | The rule alone | **passed**: seven sequences of counts give 0, 0, 0, 1, 2, 4 and 4 extension steps (skipped; skipped with water; water at once; water after one, after two, after four steps; never) |
| unit | The built hCA II systems, start coordinates | **passed** for X, N and R: the driver's count equals the second calculation (count, nearest distance, ligand distance); it does not change when every water molecule and the ligand are moved by whole cell vectors; the numbers are those of the table in section 1A |
| unit | Run X1 of the pipeline test on the Pi | **passed**, on a copy of its manifest and records: the driver of today, with the protocol file of today, accepts the run, keeps its protocol and plans `min`, `eq01`–`eq06`, `seg0001`, `seg0002` without a check |
| cost | Cost at the size of the built systems (40,902 atoms, one thread, shared machine) | 0.50 s for one count (one per step, at most 31 per run); 0.09 s to wrap and write one full frame (0.013 s for 25,600 atoms in the first tests): about 90 s for the 1,000 full frames of 100 ns at this speed |

Not repeated today: tests 0, 1, 3b, 4, 5, 6, 8, 9b, 10 and 12 of section 6.1, the launcher loop and the pre-flight
script (the change does not touch what they test, and several need two threads). Tests 2, 3a, 7 and 9 have
counterparts above: the shortened equilibration with a segment, the interruption with one thread, the refusals,
the continuation to a longer run.

## 7. Not tested, and what to watch

1. **No GPU here.** Not executed: the CUDA or OpenCL platform with mixed precision and
   DeterministicForces, the reading of the GPU name (the `DeviceName` property; the `nvidia-smi` line is
   recorded as well), and any speed. The property names are set only if the platform lists them, so a
   missing one cannot stop a run; the manifest says what was set. The speed test and the pre-flight test
   of section 2 come first on the rented machine.
   *Superseded on 6 Oct 2026:* the CUDA platform has run with mixed precision and DeterministicForces on an
   NVIDIA L40S, and the name of the GPU was read (`minimisation.platform` in the manifest of run X1); the speeds
   are in `md/PRODUCTION.md`, sections 5.1 and 7. The text above is the state before.
2. **Bit-for-bit repetition is not guaranteed.** On the Pi it holds with one thread only (test 0). With
   two or more threads, or on a GPU where the pre-flight test shows differences, a segment repeated after
   an interruption is another valid continuation from the same saved state: no frame is lost or doubled,
   but the run cannot be regenerated bit for bit. The manifest records every interruption (an invocation
   that ended with exit code 143 or 130, or has no end).
   **It also did not hold from 5 Oct to 6 Oct 2026**: a segment of the one-thread run `ref1` of the first
   tests, repeated on 6 Oct from the same saved state with the same seeds, differs from the first log row on
   (1.2 × 10⁻⁷ of the potential energy), with the old driver and with the new one alike, while every run made
   on 6 Oct repeats itself exactly (section 6.2, group compat_ref1; reported there as two failed tests).
   Suspected and not verified: the runs of 6 Oct set `OPENMM_CPU_THREADS=1`, those of 5 Oct did not, and the
   PME library of OpenMM takes its number of threads from that variable (the name is in
   `libOpenMMPME.so`), not from the platform property `Threads` that `--threads` sets. If so, a "one-thread"
   run without the variable computes the reciprocal sum on all cores, and `--threads N` on the Pi needs
   `OPENMM_CPU_THREADS=N` beside it to be what it says. The driver does not set the variable itself and
   does not record it; that would be a change of one line and a new round of tests. To verify: repeat
   `run_tests.py compat_ref1` without the two variables in `ENV` (it then uses all cores); the files of
   5 Oct should come back. None of this touches a GPU run.
   *On the rented GPU (6 Oct 2026, 9:10:30 AM IST):* the pre-flight test ended `DIFFERENT`: 18 of 25 files
   differ, `state_min.xml` among them, that is, the two runs differ from the minimisation on
   (`/root/cqc-sfa-md/status/vm/status/preflight_interrupt.json`). A segment repeated after an interruption
   there is therefore to be read as another valid continuation, not as the lost one.
3. **Crystal water in the minimisation: settled.** The builder now lists the crystal-water oxygens
   (`crystal_water_O`: 215 in X and N, 211 in R), and the driver holds them in the minimisation with the heavy
   atoms of the solute (2,289 atoms in X and N, 2,285 in R [`unit_checks.json`]; the manifest has the number in
   `system.restraint_atom_sets_n_atoms`). From `eq01` on they are free. A system without the list would be
   minimised with free crystal water; `system.crystal_water_O_listed_in_selections` shows which.
4. **The real systems on a GPU have not been run, and the bonded systems not at all** (they do not exist yet).
   *Superseded on 6 Oct 2026:* the bonded systems exist and ran on the GPU (rehearsal `md/test/REHEARSAL_2e.md`;
   validation runs and the start of production `md/PRODUCTION.md`). The rest of this point is the state before.
   The restrained system X is being run on the Pi in another agent's pipeline test (`md/test/`,
   `/root/cqc-sfa-md/runs_test/`), with the driver and the protocol as they stood before the hydration check;
   it is a test of the pipeline, not a result. Whether 2 ns settle temperature and density (gate of S3) and
   whether the zinc model is stable at 4 fs will show in run X0; the state logs and the means in each record
   are there for that.
5. **Restraints under the barostat.** The reference positions are fixed coordinates while the barostat
   rescales the box (stages 2 to 5 and the extension steps): the solute stays where it was built and the box
   settles around it. The barostat's moves are then somewhat smaller than without restraints; the density of
   the 300 ps stage in the log shows whether that mattered.
6. **Files beyond the conventions, and who reads them.** Equilibration trajectories (`solute_eqNN.dcd`,
   `full_eqNN.dcd`), the names `log_<step>.csv`, `state_<step>.xml`, `rec_<step>.json`, and since 6 Oct 2026 the
   files of the extension steps (`…_eq03x1.…` to `…_eq03x4.…`) and the manifest keys `hydration_check` and
   `equilibration.n_extension_steps`. **The analysis does not read the extension steps**:
   `md/analysis/SPEC.md` (version 1.2, section 9E) takes as a stage a file and a manifest entry whose name is
   `eq` followed by digits only. It places `eq04` to `eq06` correctly, because it takes each stage's start and
   end from the manifest, but the frames of an extension step are missing from its equilibration series, and
   its "count at the end of the held part" is the last full frame of `eq03`, not of the last extension step.
   The driver's own counts are in the manifest in any case. The checking scripts of the pipeline test
   (`md/test/pi_test_common.py`, `expected_steps`) build the chain from the six stages and would report a
   broken chain for a run with extension steps. Both are for their owners; nothing of them was changed here.
   *Superseded on 6 Oct 2026 for the specification:* `md/analysis/SPEC.md`, version 1.3 (fixed at 11:26 AM
   IST), section 9E, reads the files of the extension steps (`eq03x1` to `eq03x4`) as steps of their own, placed
   in time by `equilibration.stages` of the manifest, and sets the driver's count at the end of stage 3 and of
   every extension step beside its own (sections 9D and 9E). Version 1.2 as it was fixed at 7:18 AM IST (the
   file of commit `7a7a988`, section 9E) read them already; the sentences above that say the specification
   does not read the extension steps were written before that (by the header of this file, in its revision of
   5:05 to 6:50 AM IST) and hold for neither of the two fixed versions. Looked at again on 7 Oct 2026:
   `expected_steps` in `md/test/pi_test_common.py` (file of 6 Oct 2026, 1:07 AM IST) still builds the chain from
   the stages of the protocol and the segments, without extension steps.
7. **Sizes** in section 3 are estimates for the built systems, not measurements of a production run.
8. **The minimiser's limit is per pass: settled as documented.** OpenMM repeats L-BFGS with stiffer terms
   for the constrained bonds, and `max_iterations` bounds each pass. The passes and their iteration counts are
   in `rec_min.json` and now in the console line; if the last pass ends at the limit, the tolerance was not
   reached and the driver says so in a line that begins with `WARNING` (the run goes on: heating starts at
   50 K). For the built system X with the production limit of 5,000 the passes took 1,005, 3,619, 6 and 6
   iterations, none stopped by the limit, 1,653 s on the Pi with two threads (pipeline test of another agent,
   restrained model [`pi_run_numbers.json`]); `md/TOPOLOGY.md`, section 6, had not run it to its end.
9. The ligand files of the first stand-in are `md/ligand/sfa_anion.mol2` and `sfa_anion_gaff2.frcmod` as they
   stood at 8:45 PM on 5 Oct 2026 (checksums in the stand-in's `build.json`); those of the second stand-in as
   they stand on 6 Oct 2026 (its `build.json`). The tests do not depend on what they hold.
10. **The rule of the hydration check against the built start R** (section 1A): the count is 1 before any
    dynamics, from one solvent water at the place of N1, and no water is at the place of the ring carbons. The
    rule as worded may therefore never extend an R run although the place of the ring is empty. Whether the
    count should be made on the six ring carbons, or with another criterion, is the plan's to say **before the
    first R run**; in the driver it is the list `site_atoms` (and `water_cutoff_nm`) of the protocol. The
    analysis specification reports the same finding (its section 16.1).
    *Solved on 6 Oct 2026:* the plan changed the rule before the rehearsal, the validation runs and the
    campaign, in its entries of 8:01 AM (point 7: the count is taken on the positions of the six ring carbons C1
    to C6, without N1) and 8:56 AM (point 6: within 2.5 Å, because a water oxygen within 3.5 Å of a ring-carbon
    position is found in the crystal with the ring in place). `md/protocol.json` (`site_atoms` C1 to C6,
    `water_cutoff_nm` 0.25), the driver and the manifests of the runs follow the plan (section 1A); the rest of
    the rule stands.
11. **The second bonded parameter set has no folder name yet.** `md/build_topology.py --model bonded` writes
    to `<out-root>/bonded/<start>/`; a second set needs another `--out-root` or a second model name. The
    driver takes any folder.
    *Solved on 6 Oct 2026:* the builder takes `--set 2e|2s` and `--variant NAME` and writes to
    `<out-root>/bonded_<set>[_<variant>]/<start>/`; the folders are `bonded_2e`, `bonded_2s`, `bonded_2e_n2fit`
    and `bonded_2s_n2fit` (`md/PRODUCTION.md`, section 3).
12. **The fit of the CA atoms needs CA atoms that do not lie on a line.** With the first stand-in (five CA atoms
    of an extended pentapeptide) the fitted site swings about the chain by several Å, and a ligand held in
    place was 2.7 Å from "its own" site at the end of `eq03`. hCA II has 258 CA atoms; the r.m.s. deviation of
    the fit is in every check record (`ca_fit.rmsd_A`) and should stay well below 1 Å during the held part.
    The driver also records, and does not enforce, that the residue names of the CA atoms in `reference.pdb`
    equal those of the topology (`hydration_check.reference_guard`; 0 of 258 differ in the three restrained
    systems). For the bonded systems, whose zinc residues are renamed, look at this number in the first
    manifest.
13. **Runs and tools on the Pi that use the driver and the protocol of this folder.** (a) Run X1 of the
    pipeline test began before the change and will be resumed after it: by the exception of section 2 it
    continues without the check (verified on a copy of its manifest, section 6.2). (b) Any run that **starts**
    after 6 Oct 2026, 5:28 AM IST gets the check: run R1 of that test would make up to four extension steps of
    15 ps (`--equil-scale 0.03`), about an hour each at the speed of that test, and its checking scripts do
    not know them (point 6). (c) `md/build_topology.py --min-check --short-md` runs the driver with a copy of
    `md/protocol.json`: for start R it may now add up to four extension steps of 1 ps to its 4 ps test, and it
    reads the six stages by name, which still works.
