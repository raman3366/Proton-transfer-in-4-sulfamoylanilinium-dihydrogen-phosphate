# The rented GPU machine: what runs on it, in which order, and how to follow it

Stage S5 of `md/PLAN.md`, and the quantum chemistry of stage S2 that the Pi cannot do in time. Written on
6 Oct 2026, 5:00 to 7:00 AM IST, on the Pi. The user's decision of that day: "Go with spot on ibm01. Start, run
and report every 2 hours on status and eta. Monitor remaining credits."

**State when this was written: nothing has been created, started or spent.** These files only prepare the side
of the machine. Creating, starting and stopping the machine (`ops/ibm/cqc_gpu_up.py`, `ops/ibm/cqc_gpu.py`,
`ops/ibm/GPU.md`), the watchdog on the Pi and the credit check (`ops/pi/cqc_gpu_watchdog.py`,
`ops/pi/cqc_gpu_status.py`) are other files, written at the same time by others; this text names what they
can rely on, and where it says what they do, it says when I read them.

| File | Runs on | What it is |
|---|---|---|
| `md/gpu/push_to_vm.sh` | Pi | copies what the machine needs to the same absolute paths; checksum list; dry run by default |
| `ops/vm/setup_gpu.sh` | machine, once | NVIDIA driver, micromamba, the two environments, checks, `status/setup.json` |
| `md/gpu/queue.py` | machine | the task runner: one GPU task at a time, CPU tasks beside it, restartable, heartbeat |
| `md/gpu/cqc-md-queue.service` | machine | starts the runner at boot, so that work goes on after a pre-emption |
| `md/gpu/deadman.py`, `cqc-md-deadman.service`, `cqc-md-deadman.timer` | machine | the switch: shuts the machine down when the Pi has not been seen for 4 hours (section 4; added 6 Oct 2026 after the safety review) |
| `md/gpu/tasks_initial.json` | machine | the first six tasks (section 3) |
| `md/gpu/GPU_RUN.md` | | this text |

Machine assumed: IBM Cloud VPC profile `gx3-24x120x1l40s` (24 vCPU, 120 GB, one L40S with 48 GB, **no instance
disk**: everything lives on the boot volume), image `ibm-ubuntu-24-04-5-minimal-amd64-1`, login user `ubuntu`
with sudo. Everything on the machine runs as **root** and uses the paths of the Pi:
`/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md`, `/root/cqc-sfa-md/`, `/root/cqc-sfa-tools/`.
The boot volume should leave at least 40 GB free after the set-up (section 7).

## 1. In short

1. Push, set up, add the first tasks, enable the service. Four commands from the Pi (section 2).
2. From then on the machine works alone. After a pre-emption and restart the service starts the runner, the
   runner starts every unfinished task again, and every task of this project continues from its own records.
3. The Pi follows through four small files in `/root/cqc-sfa-md/status/` on the machine (section 6) and copies
   results back (section 7).
4. A task that fails three times stops, and so does everything that depends on it. The machine then idles and
   **keeps costing money**: `heartbeat.json` says `idle: true` with the reason under `attention`.
5. **The watchdog stops the machine when the runner reports nothing left to do.** As
   `ops/pi/cqc_gpu_watchdog.py` stood on 6 Oct 2026, 6:38 AM IST: every listed task done for 20 minutes means a
   last copy, the verification and a stop; idle without being done for 3 hours means a stop. So the next task
   file has to be added within 20 minutes after the last task ends, or the machine has to be started again when
   the file is ready (a new spot request).

## 2. Order of work (from the Pi; `<ip>` is the address of the machine)

```bash
SSH="ssh -i /root/.ssh/aitf_pi -o IdentitiesOnly=yes -o BatchMode=yes -o ConnectTimeout=15 -o StrictHostKeyChecking=accept-new ubuntu@<ip>"
REPO=/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate
Q="sudo python3 $REPO/md/gpu/queue.py"

# 1. copy (first a look at what would go: 762 files, 86.1 MB on 6 Oct 2026, 9:20 AM IST; it grows with the Pi's
#    own quantum-chemistry jobs: 459 files, 57.3 MB at 6:27 AM, 693 files, 72.4 MB at 8:11 AM)
$REPO/md/gpu/push_to_vm.sh --dry-run
$REPO/md/gpu/push_to_vm.sh --preview                # ON A MACHINE THAT ALREADY WORKS, first: asks the machine and lists what
                                                    # --yes would REPLACE there; changes nothing (see the notes below)
$REPO/md/gpu/push_to_vm.sh --yes                    # ends with VERIFIED and PUSH COMPLETE; the address comes
                                                    # from /root/.config/cqc-sfa/vm_gpu.json, or from --ip <ip>

# 2. software (detached, so that a lost connection does not end it); the LAST line of the log is the result
$SSH "sudo setsid nohup bash $REPO/ops/vm/setup_gpu.sh < /dev/null > /dev/null 2>&1 &"
$SSH "sudo tail -n 3 /root/setup_gpu.log"
#    SETUP_GPU RESULT: REBOOT_NEEDED  ->  $SSH "sudo reboot" ; wait until SSH answers ; start the script again
#    SETUP_GPU RESULT: COMPLETE       ->  go on.   FAILED -> read the log; nothing below is done.
$SSH "sudo cat /root/cqc-sfa-md/status/setup.json"   # versions, checks, what differs from the Pi's environment
#    The set-up also installs AND starts the switch that shuts the machine down when the Pi is silent (section 4).
#    On a machine that was set up before the switch existed (the one of 6 Oct 2026, 8:36 AM IST) and that is
#    WORKING: do not push everything for this (the notes below). Copy the four files of the switch, then run its step:
#    (these two rsync lines are the second reviewer's proposal; they have NOT been run from here)
rsync -a -e "ssh -i /root/.ssh/aitf_pi -o IdentitiesOnly=yes -o BatchMode=yes -o ConnectTimeout=15" --rsync-path="sudo -n rsync" \
  $REPO/md/gpu/deadman.py $REPO/md/gpu/cqc-md-deadman.service $REPO/md/gpu/cqc-md-deadman.timer ubuntu@<ip>:$REPO/md/gpu/
rsync -a -e "ssh -i /root/.ssh/aitf_pi -o IdentitiesOnly=yes -o BatchMode=yes -o ConnectTimeout=15" --rsync-path="sudo -n rsync" \
  $REPO/ops/vm/setup_gpu.sh ubuntu@<ip>:$REPO/ops/vm/
$SSH "sudo bash $REPO/ops/vm/setup_gpu.sh --deadman-only"        # seconds; nothing else is touched. It runs the switch's
#    dry run FIRST and refuses to enable the timer while that says the Pi has not been seen (a stale status/PI_SEEN
#    would shut the machine down at the first run): then look at `ls -l /root/cqc-sfa-md/status/PI_SEEN` and at the
#    watchdog on the Pi, and run the step again.

# 3. the first tasks, then the runner (this is the moment the machine starts to work)
$SSH "$Q add $REPO/md/gpu/tasks_initial.json"
$SSH "sudo systemctl enable --now cqc-md-queue.service"
$SSH "$Q status"
```

- The machine itself is created by `ops/ibm/cqc_gpu_up.py` (see `ops/ibm/GPU.md`); nothing here does that.
- `push_to_vm.sh` installs `rsync` on the machine if the image lacks it, and can be run again at any time; it
  copies only what changed on the Pi and never deletes on the machine. It leaves out the files that the
  machine's own quantum-chemistry tasks write (`qm/prod/status.json`, `qm/prod/opt*`, `qm/prod/hess*` and the
  like), so a later call cannot overwrite them with the Pi's. **A file that is newer on the machine than on the
  Pi is left as it is there** (`rsync --update`, since 6 Oct 2026): the machine's `qm_write` task rewrites files
  inside `md/` (`md/qm/results/CA2_qm_summary.json`, `md/qm/mcpb/CA2_large_mk.log`), and a later push must not
  put the Pi's older ones over them. Such files are named in `/root/cqc-sfa-md/status/push_kept_newer.txt`
  (on both sides) and in the log, and are left out of the checksum list, so the push still ends with VERIFIED.
  To force the Pi's version of such a file onto the machine, `touch` it on the Pi and push again.
- **A full push is not a neutral act while tasks run.** Every file that is newer on the Pi and of other content
  is put over the machine's: scripts (`md/run_md.py`, `md/qm/*.py`), `md/protocol.json`, task files, inputs. A
  task that reads such a file at its next step, or at its next start after a pre-emption, then runs with
  another one than it began with. On 6 Oct 2026, 10:16 AM IST, with `gr_points` and `qm_opt` running, 16 files
  of other content would have been replaced that way. **`push_to_vm.sh --preview` asks the machine first**: rsync
  dry runs with checksums on both sides, per part the files that would be new, those that would REPLACE a file
  of other content (each named), those that would only get the Pi's time stamp, and those that stay because
  they are newer on the machine. It copies, writes, installs and logs nothing. Do the full push at a moment the
  main session chooses (nothing running that reads the files named, or the queue drained), and read
  `push_kept_newer.txt` afterwards. To bring single files to the machine, copy just those (as for the switch,
  above).
- `setup_gpu.sh` can be run again at any time; `--checks-only` repeats only the checks.
- `systemctl enable` is what makes the runner start at boot. `queue.json` and `heartbeat.json` carry
  `starts_at_boot: true` when it is set; **without it a pre-empted machine comes back and does nothing.**
- To stop everything cleanly: `$SSH "$Q stop"` (section 4). Stopping or deleting the machine is the
  provisioning script's business.

## 3. The first tasks (`tasks_initial.json`)

They need nothing that did not exist on 6 Oct 2026.

| Task | Resource | After | What it does | Done when |
|---|---|---|---|---|
| `preflight_gpu` | gpu | | `nvidia-smi`; `python -m openmm.testInstallation`; the CUDA platform must compute forces and the platforms must agree | `status/preflight_gpu.json` |
| `preflight_speed` | gpu | `preflight_gpu` | the speed test of `md/RUN.md`: `run_md.py --benchmark-seconds 600` on `systems/restrained/X` (40,902 atoms, the size of production) in `runs/bench` | `status/speed.json` (ns per day, hours for 100 ns and for 900 ns, GPU name, which tasks ran beside it, load) |
| `preflight_interrupt` | gpu | `preflight_gpu` | the pre-flight interruption test of `md/RUN.md`: the same short run twice (`--equil-scale 0.05 --prod-ns 0.4 --segment-ns 0.2`), the second killed with SIGKILL after 60 % of the wall time of the first, then resumed; checksums of all trajectories and saved states compared | `status/preflight_interrupt.json` with `result` IDENTICAL or DIFFERENT |
| `qm_opt` | cpu:16 | | `run_qm.py config_prod.json --jobs opt`: optimisation of the small zinc model (64 atoms, B3LYP/6-31G*, `md/ZINC.md` section 7.3) | `qm/prod/status.json` says `jobs.opt.state = done`, `optimised = true`, and `opt_final.npy` exists |
| `qm_hess` | cpu:16 | `qm_opt` | `--jobs hess`: 37 gradients for the force constants of the zinc, its four nitrogens and the sulfur | `jobs.hess.state = done` and `hess_columns.npz` exists |
| `qm_write` | cpu:1 | `qm_hess` | `--jobs write`: `md/qm/mcpb/CA2_small_opt.fchk` (what MCPB.py step 2s reads) and the summary with Seminario's constants; seconds | the fchk file exists and `jobs.write.fchk.partial = true` |

- **Order.** The three GPU tasks run one after another in this order; the quantum chemistry starts at once and
  runs beside them on 16 of the 24 vCPUs (the runner sets `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and
  `MKL_NUM_THREADS` to 16, the variables `md/qm/run_qm.sh` sets to 4 on the Pi; `config_prod.json` is used
  unchanged). CPU tasks run at nice 10, GPU tasks at nice 0, so the one CPU thread that feeds the GPU goes first.
  **The speed test therefore runs with the quantum chemistry beside it**; `speed.json` says so
  (`tasks_running_beside_at_the_end`, `load_average_at_the_end`). The speed of every production segment is in
  its own record in any case.
- **DIFFERENT is not a failure** of the interruption test (`md/RUN.md`, section 2): no frame is lost or
  doubled, but a run cannot be regenerated bit for bit. The task fails only if the test itself is not sound: the
  second run was not killed inside the run, the resumed run is not complete, or the two runs do not have the
  same files and trajectory sizes. It then writes `status/preflight_interrupt_FAILED.json`, and no MD task that
  depends on it starts.
- **The pre-flight tests use the restrained model**, the only built system. `md/TOPOLOGY.md` (section 5.4)
  shows that this model does not keep the zinc site; that does not matter for speed and repeatability. It ran
  through the shortened equilibration on the Pi's CPU (`/root/cqc-sfa-md/runs_test/restrained/X1`, 6 Oct 2026).
  If a pre-flight run ends with exit code 4 (energy not finite), repeat the test with the bonded system when it
  exists (a new task with another `PF_SYSTEM`); do not start production on an unsound test.
  (The bonded systems exist since 6 Oct 2026 as `bonded_2e` and `bonded_2s`; both pre-flight runs ended without
  that exit code, so the test was not repeated.)
- `qm_write` was added to the two jobs named for this step because it costs seconds and gives the file the next
  step needs. On the machine it also rewrites `md/qm/results/CA2_qm_summary.json`; see section 7 for the copy back.
- All three QM jobs are restartable from their own files (`opt_last.npy`, `hess_disp/*.npy`); the driver skips a
  job whose result exists (read in `md/qm/run_qm.py`).
- **`qm_opt` and `qm_hess` have a time limit: `"timeout_s": 43200`, 12 hours for one attempt** (added 6 Oct 2026
  after the safety review). The runner then ends the attempt (SIGTERM, SIGKILL after 90 s) and counts it as
  failed; the next attempt goes on from the job's own files; after three attempts (36 hours at most) the task
  is `failed` and everything after it waits for a person. Without a limit a job that hangs, or an optimisation
  that does not converge in its 300 steps, would hold 16 CPUs for as long as the machine lives. **12 hours is a
  guess**, made before either job had run on this machine: when the first optimisation steps and the first
  gradients have been timed (`qm/prod/status.json`: seconds per step; 37 gradients for `qm_hess`), replace it
  by a number with a reason, and put the projected hours of the quantum chemistry into the two-hourly report
  beside the speed test (the watchdog's projection covers the MD runs only). A changed limit reaches a task
  that is already running only at its next start (the runner says so in `queue.log`); for the machine of
  6 Oct 2026, on which `tasks_initial.json` was read at 8:55 AM IST without the limits: push, then
  `$SSH "$Q add --replace $REPO/md/gpu/tasks_initial.json"`.

## 4. The runner (`queue.py`)

- **Tasks** are JSON files in `/root/cqc-sfa-md/queue/tasks.d/`; the runner reads the folder again every
  10 seconds. Fields of a task: `id`, `command` (text or a list of lines, run by `bash -eo pipefail -c`), `cwd`,
  `resources` (`"gpu"`, `"cpu:N"`, `"gpu+cpu:N"`), `depends_on`, and `done_marker` and/or `success_test`. The
  full description is the header of `queue.py`.
- **Done means: the success test holds.** A task whose test holds already is marked done without being run. A
  task that ran is done when its command ended with exit code 0 **and** the test holds. At every start of the
  runner the test of every task recorded as done is made again.
- **One GPU task at a time**; CPU tasks beside it while the sum of their CPUs fits into the 24 vCPUs.
- **Interruptions.** A task that was running when the runner or the machine went down is started again
  (`interruptions` counts them; after 30 the task is failed). Processes of a dead runner that are still alive
  are ended first, so a task never runs twice at once.
- **Failures.** A command that fails, or ends with 0 without its test holding, is tried again after 60 s, three
  times in all; then the task is `failed` and everything after it waits. `queue.py retry ID` gives it new
  attempts; a changed definition of the task does the same.
- **Stopping.** `queue.py stop` writes `/root/cqc-sfa-md/queue/STOP`: every running task gets SIGTERM
  (`run_md.py` then ends within 500 steps and leaves a consistent folder), SIGKILL after its `stop_grace_s`; the
  runner ends with exit code 0 and systemd leaves it stopped. **While the STOP file exists nothing starts, also
  not after a reboot.** `queue.py resume` removes the file; `systemctl start cqc-md-queue` starts the runner.
  `queue.py drain` lets running tasks finish and starts nothing new. `queue.py cancel ID` ends one task.
- **The switch on the machine (`md/gpu/deadman.py`).** Everything that enforces the cost ceiling runs on the
  Pi. For the case that the Pi, its network or its cron stops, the Pi's watchdog touches
  `/root/cqc-sfa-md/status/PI_SEEN` at every contact (every 5 minutes), and a systemd timer on the machine
  (`cqc-md-deadman.timer`: first 45 minutes after boot, then every 10 minutes) runs `deadman.py --yes`, which
  looks at the age of that file. **Older than 240 minutes:** `queue.py stop` (the clean stop above), up to
  10 minutes of waiting for the runner to end, the STOP file it made removed again (so that the queue goes on at
  the next start; a STOP file of a person stays), then `shutdown -h now`. The boot volume keeps everything.
  If `PI_SEEN` was never touched, the 240 minutes count from the boot, or from the first run of the switch if
  that was later (from the boot alone if the switch can neither read nor write its own `status/deadman.json`:
  it would otherwise begin its count again at every run). After a start of the machine the Pi has 45 minutes
  to be seen, else it goes down again. If the Pi reaches the machine but its touch of `PI_SEEN` fails, the
  watchdog raises `ALERT_PI_NOT_SEEN`.
  - Installed **and started** by `setup_gpu.sh` (step 6b; alone: `setup_gpu.sh --deadman-only`). The step runs
    the switch's dry run first (`deadman.py --state`: `armed`, `silent` or `off`) and **refuses to enable the
    timer while the answer is `silent`**: on a machine that is up for more than 45 minutes the timer runs the
    switch the moment it is enabled, and with a stale `PI_SEEN` that run would shut the machine down. Look:
    `$SSH "systemctl list-timers cqc-md-deadman.timer; sudo cat /root/cqc-sfa-md/status/deadman.json; sudo tail -n 5 /root/cqc-sfa-md/status/deadman.log"`.
    On the Pi, the row "Guard" of `python3 ops/pi/cqc_gpu_status.py` says whether the switch is armed, and
    `ALERT_DEADMAN_OFF` is up while it is not.
  - **Off:** `$SSH "sudo touch /root/cqc-sfa-md/status/NO_DEADMAN"` (for work by hand with the watchdog switched
    off on the Pi); remove the file to arm it again. `deadman.py` without `--yes` is a dry run.
  - **Order on a machine that already runs:** first the new watchdog on the Pi (it is what touches `PI_SEEN`:
    `$SSH "ls -l /root/cqc-sfa-md/status/PI_SEEN"` shows a time of the last 5 minutes), then the switch.
  - **Test it once, in the first hour, when nothing of value runs** (it was only ever run against a faked clock
    and faked commands on the Pi, and whether IBM reports an instance that was shut down from inside as
    `stopped` is not verified): `touch /root/cqc-sfa-md/status/HOLD` on the Pi, then
    `$SSH "sudo python3 $REPO/md/gpu/deadman.py --yes --limit-min 0"`; within a few minutes
    `python3 ops/ibm/cqc_gpu.py status` on the Pi must show `stopped`, and the balance must stop falling at the
    running rate. Then remove `HOLD`: the watchdog takes the stop for a pre-emption and starts the machine again
    (or `python3 ops/ibm/cqc_gpu.py start --yes`), and the queue goes on by itself. If IBM does **not** report
    `stopped`, the switch does not end the billing: say so at once, and stop the machine with
    `python3 ops/ibm/cqc_gpu.py stop --yes`.
- **Files it writes**, all under `/root/cqc-sfa-md/`:

| File | Content |
|---|---|
| `status/queue.json` | everything: state of the runner, and for every task its definition, state, attempts, interruptions, last test and its history (each change of state with UTC and IST time) |
| `status/heartbeat.json` | rewritten every minute and at every change: time, `queue_state`, `running` (id, pid, since, seconds), `counts`, `idle`, `idle_since_utc`, `all_done`, `attention`, CPUs in use, load, free disk, uptime, one line of `nvidia-smi`, `starts_at_boot` |
| `status/queue.log` | one line per event, UTC and IST |
| `logs/<id>.log` | output of the task, appended, with a line at the start and at the end of every attempt |
| `queue/done/<id>.done` | only for tasks with `"done_marker": "auto"` |

States of a task: `pending`, `running`, `done`, `retry` (failed once, waits for the next attempt), `interrupted`
(only between a stop and the next start), `failed`, `cancelled`.

## 5. Adding tasks

1. Write a JSON file on the Pi (any name ending in `.json`; the name decides the order among files).
2. Check it there: `python3 md/gpu/queue.py check FILE --root /tmp/none --cpus 24` (syntax and resources), copy it
   to the machine (`scp`, or put it under `md/gpu/` and run `push_to_vm.sh --yes`, which leaves alone every file
   that is newer on the machine: section 2).
3. On the machine: `sudo python3 $REPO/md/gpu/queue.py add FILE`. This checks the file against the tasks that
   are there (double ids, tasks that need more CPUs than the machine has, circles) and copies it into
   `tasks.d/` in one step. The running runner picks it up within 10 seconds. A file that fails the check is not
   added; a file that is not JSON and is put into the folder by hand is reported in `queue.json` and ignored.

Rules: ids are unique over all files. Every task must be safe to start again from the beginning of its
command. A task must have a test that is true only when its result is really there. Results that a later
task reads should be named in `depends_on`, not assumed from the order.

### Ready-made: validation run X0 (stage S4)

**As used on 6 Oct 2026** (`md/PRODUCTION.md`): the bonded systems are in one folder for each zinc set,
`bonded_2e` (empirical, the primary set) and `bonded_2s` (Seminario, the sensitivity set); there is no folder
`bonded`. The files that were added to the queue are `md/gpu/tasks_x0.json` (X0 of both sets, ids `X0_2s` and
`X0_2e`), `md/gpu/tasks_production.json` (the nine runs with `bonded_2e`; X1 waits for `X0_2e`) and
`md/gpu/tasks_sensitivity.json` (`R_zn` with `bonded_2s`, `R_n2` with `bonded_2e_n2fit`). They follow the entries
below, with the model name written out and a time limit of one attempt (`timeout_s`: one hour for X0, six hours
for a run of 100 ns). The entries below stay as the pattern.

Fill in `MODEL` (the folder under `/root/cqc-sfa-md/systems/` that holds `X`, `N`, `R`: `bonded_2e` for the primary
set). The
command is that of `md/RUN.md` with `--platform CUDA`; run id X0 and 5 ns as the plan's line of 5 Oct, 9:41 PM,
says. If the system is built by a task on the machine, add that task's id to `depends_on`.

```json
{
 "comment": "Validation run X0 (md/PLAN.md stage S4, gate G-V): 5 ns from the crystal pose. Fill in MODEL.",
 "env": {
  "MD": "/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md",
  "PY": "/root/cqc-sfa-tools/md/bin/python",
  "SYS": "/root/cqc-sfa-md/systems",
  "RUNS": "/root/cqc-sfa-md/runs",
  "MODEL": "FILL_IN_THE_MODEL_NAME"
 },
 "tasks": [
  {"id": "X0",
   "title": "validation run X0: 5 ns from the crystal pose (gate G-V)",
   "cwd": "$RUNS",
   "resources": "gpu",
   "depends_on": ["preflight_speed", "preflight_interrupt"],
   "stop_grace_s": 120,
   "command": "\"$PY\" \"$MD/run_md.py\" --system \"$SYS/$MODEL/X\" --out \"$RUNS/$MODEL/X0\" --run-id X0 --platform CUDA --stage all --prod-ns 5 --segment-ns 5 --resume",
   "success_test": "python3 \"$MD/gpu/queue.py\" json \"$RUNS/$MODEL/X0/manifest.json\" run_id=X0 status=complete production.complete=true production.completed_ns=5 test_only=false"}
 ]
}
```

### Ready-made: the nine production runs (stage S5)

**Add this file only after gate G-V has been judged on the Pi**: the runner cannot judge a gate, and X1 would
start as soon as X0 is complete. Order X1 N1 R1, X2 N2 R2, X3 N3 R3 as in `md/RUN.md`; each run waits for the
one before it, so the chain stops at the first run that fails, like the loop of `md/RUN.md`.

```json
{
 "comment": "The nine production runs (md/PLAN.md stage S5), one after another in the order of md/RUN.md: X1 N1 R1, X2 N2 R2, X3 N3 R3. Each waits for the one before it, so the chain stops at the first run that fails. Fill in MODEL. Add only after gate G-V has been judged.",
 "env": {
  "MD": "/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md",
  "PY": "/root/cqc-sfa-tools/md/bin/python",
  "SYS": "/root/cqc-sfa-md/systems",
  "RUNS": "/root/cqc-sfa-md/runs",
  "MODEL": "FILL_IN_THE_MODEL_NAME"
 },
 "tasks": [
  {"id": "X1",
   "title": "production run X1: 100 ns from start X, replica 1",
   "cwd": "$RUNS",
   "resources": "gpu",
   "depends_on": ["X0", "preflight_speed", "preflight_interrupt"],
   "stop_grace_s": 120,
   "command": "\"$PY\" \"$MD/run_md.py\" --system \"$SYS/$MODEL/X\" --out \"$RUNS/$MODEL/X1\" --run-id X1 --platform CUDA --stage all --prod-ns 100 --segment-ns 5 --resume",
   "success_test": "python3 \"$MD/gpu/queue.py\" json \"$RUNS/$MODEL/X1/manifest.json\" run_id=X1 status=complete production.complete=true production.completed_ns=100 test_only=false"},
  {"id": "N1",
   "title": "production run N1: 100 ns from start N, replica 1",
   "cwd": "$RUNS",
   "resources": "gpu",
   "depends_on": ["X1"],
   "stop_grace_s": 120,
   "command": "\"$PY\" \"$MD/run_md.py\" --system \"$SYS/$MODEL/N\" --out \"$RUNS/$MODEL/N1\" --run-id N1 --platform CUDA --stage all --prod-ns 100 --segment-ns 5 --resume",
   "success_test": "python3 \"$MD/gpu/queue.py\" json \"$RUNS/$MODEL/N1/manifest.json\" run_id=N1 status=complete production.complete=true production.completed_ns=100 test_only=false"},
  {"id": "R1",
   "title": "production run R1: 100 ns from start R, replica 1",
   "cwd": "$RUNS",
   "resources": "gpu",
   "depends_on": ["N1"],
   "stop_grace_s": 120,
   "command": "\"$PY\" \"$MD/run_md.py\" --system \"$SYS/$MODEL/R\" --out \"$RUNS/$MODEL/R1\" --run-id R1 --platform CUDA --stage all --prod-ns 100 --segment-ns 5 --resume",
   "success_test": "python3 \"$MD/gpu/queue.py\" json \"$RUNS/$MODEL/R1/manifest.json\" run_id=R1 status=complete production.complete=true production.completed_ns=100 test_only=false"},
  {"id": "X2",
   "title": "production run X2: 100 ns from start X, replica 2",
   "cwd": "$RUNS",
   "resources": "gpu",
   "depends_on": ["R1"],
   "stop_grace_s": 120,
   "command": "\"$PY\" \"$MD/run_md.py\" --system \"$SYS/$MODEL/X\" --out \"$RUNS/$MODEL/X2\" --run-id X2 --platform CUDA --stage all --prod-ns 100 --segment-ns 5 --resume",
   "success_test": "python3 \"$MD/gpu/queue.py\" json \"$RUNS/$MODEL/X2/manifest.json\" run_id=X2 status=complete production.complete=true production.completed_ns=100 test_only=false"},
  {"id": "N2",
   "title": "production run N2: 100 ns from start N, replica 2",
   "cwd": "$RUNS",
   "resources": "gpu",
   "depends_on": ["X2"],
   "stop_grace_s": 120,
   "command": "\"$PY\" \"$MD/run_md.py\" --system \"$SYS/$MODEL/N\" --out \"$RUNS/$MODEL/N2\" --run-id N2 --platform CUDA --stage all --prod-ns 100 --segment-ns 5 --resume",
   "success_test": "python3 \"$MD/gpu/queue.py\" json \"$RUNS/$MODEL/N2/manifest.json\" run_id=N2 status=complete production.complete=true production.completed_ns=100 test_only=false"},
  {"id": "R2",
   "title": "production run R2: 100 ns from start R, replica 2",
   "cwd": "$RUNS",
   "resources": "gpu",
   "depends_on": ["N2"],
   "stop_grace_s": 120,
   "command": "\"$PY\" \"$MD/run_md.py\" --system \"$SYS/$MODEL/R\" --out \"$RUNS/$MODEL/R2\" --run-id R2 --platform CUDA --stage all --prod-ns 100 --segment-ns 5 --resume",
   "success_test": "python3 \"$MD/gpu/queue.py\" json \"$RUNS/$MODEL/R2/manifest.json\" run_id=R2 status=complete production.complete=true production.completed_ns=100 test_only=false"},
  {"id": "X3",
   "title": "production run X3: 100 ns from start X, replica 3",
   "cwd": "$RUNS",
   "resources": "gpu",
   "depends_on": ["R2"],
   "stop_grace_s": 120,
   "command": "\"$PY\" \"$MD/run_md.py\" --system \"$SYS/$MODEL/X\" --out \"$RUNS/$MODEL/X3\" --run-id X3 --platform CUDA --stage all --prod-ns 100 --segment-ns 5 --resume",
   "success_test": "python3 \"$MD/gpu/queue.py\" json \"$RUNS/$MODEL/X3/manifest.json\" run_id=X3 status=complete production.complete=true production.completed_ns=100 test_only=false"},
  {"id": "N3",
   "title": "production run N3: 100 ns from start N, replica 3",
   "cwd": "$RUNS",
   "resources": "gpu",
   "depends_on": ["X3"],
   "stop_grace_s": 120,
   "command": "\"$PY\" \"$MD/run_md.py\" --system \"$SYS/$MODEL/N\" --out \"$RUNS/$MODEL/N3\" --run-id N3 --platform CUDA --stage all --prod-ns 100 --segment-ns 5 --resume",
   "success_test": "python3 \"$MD/gpu/queue.py\" json \"$RUNS/$MODEL/N3/manifest.json\" run_id=N3 status=complete production.complete=true production.completed_ns=100 test_only=false"},
  {"id": "R3",
   "title": "production run R3: 100 ns from start R, replica 3",
   "cwd": "$RUNS",
   "resources": "gpu",
   "depends_on": ["N3"],
   "stop_grace_s": 120,
   "command": "\"$PY\" \"$MD/run_md.py\" --system \"$SYS/$MODEL/R\" --out \"$RUNS/$MODEL/R3\" --run-id R3 --platform CUDA --stage all --prod-ns 100 --segment-ns 5 --resume",
   "success_test": "python3 \"$MD/gpu/queue.py\" json \"$RUNS/$MODEL/R3/manifest.json\" run_id=R3 status=complete production.complete=true production.completed_ns=100 test_only=false"}
 ]
}
```

Notes on these entries:

- The test asks for `test_only=false`: a run made with `--equil-scale` below 1 or with a test protocol can never
  count as a production run.
- No `*.console` files are written (`md/RUN.md`, section 2, uses them in its loop): the output of a run is in
  `/root/cqc-sfa-md/logs/<id>.log`, and every message of the driver is also in the run's own `run_md.log`.
- A second model (rule 4 of the plan: the sensitivity set) needs other ids, for example `R1_set2`, because an id
  can exist only once; `--run-id` should then stay `R1` only if a different seed is not wanted. Decide and
  write it down before adding.
  **Decided, 6 Oct 2026** (`md/gpu/tasks_sensitivity.json`, `md/PRODUCTION.md`, section 6): the task ids are
  `R_zn` and `R_n2`; both runs have `--run-id R1` in a folder `R1` of their own model (`runs/bonded_2s/R1`,
  `runs/bonded_2e_n2fit/R1`), so their seeds are those of run R1 of the primary set and the analysis reads them
  as start R, replica 1, of another model. Each task writes `SENSITIVITY_RUN.txt` into its run folder first.
- **Open, not solved here:** the plan's line of 6 Oct, 1:03 AM, point 3 (if no water is found where the crystal
  ring lies at the end of the held part of the equilibration of an R run, the held part is extended in steps of
  500 ps) is **not implemented in `md/run_md.py`** as it stands (searched 6 Oct 2026). The entries R1 to R3 run
  the protocol as it is. Someone has to decide before they are added.
  **No longer open, 6 Oct 2026:** the check is in the driver and the protocol (`md/RUN.md`, section 1A) and ran
  in the rehearsal (`md/test/REHEARSAL_2e.md`) and in both validation runs.

Later tasks named in the brief (charges, gate G-R, building the bonded systems) are not written here: their
commands and inputs did not exist when this was written.

## 6. Following the work

From the Pi (`$SSH` and `$Q` as in section 2):

```bash
$SSH "$Q status"                                              # one table: every task, state, attempts, note
$SSH "sudo cat /root/cqc-sfa-md/status/heartbeat.json"        # the small file for the 2-hourly report
$SSH "sudo tail -n 20 /root/cqc-sfa-md/status/queue.log"      # events
$SSH "sudo tail -n 5 /root/cqc-sfa-md/logs/qm_opt.log"        # output of one task
$SSH "sudo tail -n 3 /root/cqc-sfa-md/runs/<model>/X1/run_md.log"                        # a run: last finished step, ns per day
$SSH "sudo tail -n 3 /root/cqc-sfa-md/qm/prod/run_qm.log"     # quantum chemistry: step, energy, gradient
$SSH "systemctl is-active cqc-md-queue; systemctl is-enabled cqc-md-queue"
```

What a report can be computed from:

| Question | Where |
|---|---|
| Is the runner alive? | `heartbeat.json`: `epoch` not older than about 3 minutes and `queue_state` `running`. An old heartbeat with `queue_state` `running` means the runner or the machine is down; `stopped` means it was stopped on purpose |
| Was the machine restarted? | `boot_id` and `uptime_s` in the heartbeat; `interruptions` of the tasks in `queue.json` |
| Speed of the GPU | `speed.json`: `ns_per_day`, `hours_for_100_ns` |
| Progress of a run | `runs/<model>/<run>/manifest.json`: `production.completed_ns`, `production.ns_per_day`; remaining time = remaining ns / ns per day |
| Progress of the optimisation | `qm/prod/status.json`: `jobs.opt.steps`, `jobs.opt.last` (energy, gradient, seconds of the step). The number of steps it will need is not known beforehand |
| Progress of the force constants | `jobs.hess.gradients_done` of `jobs.hess.n_gradients` |
| Anything wrong? | `heartbeat.json`: `attention` (failed and blocked tasks, unreadable task files) |
| Nothing left to do? | `idle: true` with `idle_since_utc`; `all_done: true` when every listed task is done |

`idle` is true whenever nothing runs and nothing can start, also while the machine waits for the next task
file. The machine costs the same idle or busy.

## 7. What is copied back to the Pi

The copying is the watchdog's (`ops/pi/cqc_gpu_watchdog.py`); that file is the authority. As it stood on
6 Oct 2026, 6:38 AM IST, it copies every 15 minutes, with rsync, never deleting on either side:

| On the machine | To the Pi | What it holds, and what to know |
|---|---|---|
| `/root/cqc-sfa-md/status/`, `logs/`, `queue/` | `/root/cqc-sfa-md/status/vm/` | the runner's files (section 4), `setup.json`, `speed.json`, `preflight_*`, the two `*_env_linux64_explicit.txt`, the output of every task, the task files |
| `/root/cqc-sfa-md/runs/` | the same path | only the files of steps that have their `rec_<step>.json`; no `checkpoint_*.chk`. Each finished step is checked against the sha256 that `run_md.py` recorded. `bench` and `preflight_*` are tests |
| `/root/cqc-sfa-md/qm/`, `systems/` | the same path | a file the Pi made itself is never overwritten: the machine's version then goes to `from_gpu/conflicts/` |
| `<repository>/md/` on the machine | **`/root/cqc-sfa-md/from_gpu/repo_md/`**, with the same path below it | only files that differ in size or time from the Pi's repository. **The watchdog never writes into the repository**: what is wanted from here is copied into `md/` by hand |

What that rule means for the first tasks, so that nobody looks in the wrong place (read again in the watchdog
on 6 Oct 2026, 9:00 AM IST, after the safety review had found this passage wrong):

- The machine's `qm/prod/status.json` and `run_qm.log` (optimisation and force constants) **meet the Pi's own
  files of the same names** (the potential job and the `timing` job, which was running on the Pi at 6:07 AM
  IST). The Pi's stay; the machine's are under `from_gpu/conflicts/`.
- New files under `/root/cqc-sfa-md/` arrive at their own paths: `qm/prod/opt_steps.jsonl`, `opt_last.npy`,
  `opt_final.npy`, `hess_disp/`, `hess_columns.npz`.
- **Files of the repository folder `md/` do not arrive in the Pi's `md/`.** What the machine's tasks write
  there lands under `/root/cqc-sfa-md/from_gpu/repo_md/`: `qm/results/CA2_small_opt.xyz`,
  `qm/results/CA2_qm_summary.json` (which `qm_write` rewrites on the machine; the Pi's own stays where it is),
  and **`qm/mcpb/CA2_small_opt.fchk`**, the file `md/qm/after_qm.sh 2s` reads on the Pi. That one has to be
  copied by hand before the next step:
  `cp /root/cqc-sfa-md/from_gpu/repo_md/qm/mcpb/CA2_small_opt.fchk md/qm/mcpb/`.
- `md/qm/mcpb/CA2_large_mk.log` is written again on the machine by `qm_write` from the same `esp.npz`; it should
  come out identical to the Pi's (not checked here). If it differs in size or time, it too is under
  `from_gpu/repo_md/`.

Software (`/root/cqc-sfa-tools/`) is not copied; `setup.json` and the two `*_explicit.txt` files record what was
installed.

Size on the machine: a 100 ns run of 40,900 atoms is about 1.1 GB (10,000 solute frames of 4,094 atoms and 1,000
full frames, 12 bytes per atom and frame, and 20 saved states of 6.7 MB, the size measured in
`runs_test/restrained/X1`), so about 10 GB for the nine runs. The small model's fitting integrals fit into
memory; the large model's would be a 17 GB scratch file if that job were run there.

When the work is over: a last copy, the check of `md/RUN.md` section 4 for every run, and only then the end of
the machine (stage S5: "checksums verified; the machine deleted").

## 8. Facts used, and where they were read (6 Oct 2026)

| Fact | Source |
|---|---|
| Profile `gx3-24x120x1l40s`: 24 vCPU, 120 GB, one L40S (48 GB), `disks: []`, classes standard and spot; image `ibm-ubuntu-24-04-5-minimal-amd64-1` (created 23 Sep 2026) | the listings of `ibmcloud is instance-profile` and `ibmcloud is images` that the provisioning step saved at 5:00 AM IST (`scratchpad/gpuops/provision/`); I ran no ibmcloud command |
| Ubuntu 24.04 (amd64) has `nvidia-headless-no-dkms-580-server` and `nvidia-utils-580-server` 580.178.04-0ubuntu0.24.04.1, `linux-modules-nvidia-580-server-generic` 6.8.0-146.146+1 and `-ibm` 6.8.0-1066.67, and kernel-specific module packages for generic 6.8.0-101 to -146 and ibm 6.8.0-1045 to -1066. The 590-server user space depends on a 595 package that these indexes do not hold | `Packages.xz` of noble, noble-updates and noble-security (main, restricted) from archive.ubuntu.com; index dated 5 Oct 2026, 20:04 UTC |
| `ubuntu-drivers install --gpgpu` installs `nvidia-headless-no-dkms-<n>-server` and `linux-modules-nvidia-<n>-server-<flavour>` with `--no-install-recommends`, and does not pass on apt's exit code | source of `ubuntu-drivers-common` 1:0.9.7.6ubuntu3.7 (the package, unpacked in scratch) |
| OpenMM 8.6.1 for linux-64, Python 3.12: build `py312h5a97af1_0` needs `cuda-version >=12.9,<13`, build `py312hdfcc665_0` needs `>=13.0,<14`. `cuda-version` exists as 12.9 and 13.0 to 13.4. AmberTools 26.0 `cuda_None_nompi_py312hd652fd9_100`. PySCF 2.14.0 for Python 3.11 in builds for x86-64 levels 1, 3 and 4. No `gpu4pyscf` on conda-forge | conda-forge, with `micromamba search` and `micromamba create --dry-run --platform linux-64` (micromamba 2.9.0 on the Pi) |
| With the 85 pins of `setup_gpu.sh` the MD environment resolves for linux-64 with `cuda-version` 13.0, 13.4 and 12.9: 246 packages, 0.9 to 1.0 GB to download. Of the 242 packages of `md_env.yml`, 236 come out at the same version; other versions: `cuda-version`, `cuda-nvrtc`, `libcufft` (tied to the driver) and `libntlm` (1.8 for 1.4); two exist for aarch64 only; six are new (`mkl`, `llvm-openmp`, `tbb`, `libhwloc`, `rocm-core`, `ld_impl_linux-64`). Pinning all 242 does not resolve. The QM environment resolves: 70 packages, 185 MB | the same dry runs, and `setup_gpu.sh --solve-only` |
| `gpu4pyscf-cuda13x` and `-cuda12x` 1.8.1 (10 Aug 2026) need `pyscf>=2.8.0`, `cupy-cuda13x`/`12x`, `gpu4pyscf-libxc-…==0.8.1`; `cuda-toolkit` exists as 13.0.1 to 13.4.2 | pypi.org, JSON of each project |
| `https://micro.mamba.pm/api/micromamba/linux-64/2.9.0` delivers `micromamba-2.9.0-0.tar.bz2` (6,988,090 bytes) | fetched once from the Pi |
| Flags of micromamba (`create --dry-run --json --platform --override-channels -p -c -y -q`, `list --explicit --sha256`) | `micromamba --help` of 2.9.0 |
| The lines `N CUDA - Successfully computed forces` and `All differences are within tolerance.` | `openmm/testInstallation.py` of the installed OpenMM 8.6.1 |
| Job names `opt`, `hess`, `write`; what each skips and resumes; that the thread number comes from `OMP_NUM_THREADS` and not from the key `threads` of the configuration | `md/qm/run_qm.py`, `md/qm/run_qm.sh` |

Not read anywhere, and therefore decided by the script at run time or left open: which CUDA version the 580
driver reports (the script reads it from `nvidia-smi`); which kernel flavour the image runs (the script reads
`uname -r`); whether the reboot is needed. That a CUDA compiler library newer than the driver's CUDA version can
fail is my reason for the pin `cuda-version=<what the driver reports>`; I did not verify it here, the pin is the
careful choice either way.

## 9. Tested on the Pi, and what could not be

There is no GPU and no x86-64 processor here. Everything below ran on the Pi (aarch64, Ubuntu 24.04, Python
3.12.3, systemd 255), at nice 10, beside other jobs.

| What | How | Result |
|---|---|---|
| `queue.py`, the runner | `queue.py selftest`: dummy tasks (sleep, echo) against a runner started as its own process | 31 of 31 checks: dependencies; never two GPU tasks at once; CPU tasks beside a GPU task within the CPUs given; the runner killed with SIGKILL in the middle of a task, started again, the orphaned task process ended and the task run again once; a task file added while it ran; a broken file reported; double ids and oversized tasks refused; heartbeat rewritten; STOP (SIGTERM, then SIGKILL after the grace time for a task that ignores SIGTERM; nothing started while the file exists); DRAIN; retry, cancel, time limit; a task done without its test holding is failed; a task whose test holds is not run; a done task whose test no longer holds at the next start is run again; a second runner on the same folder refuses; a state file with unusable content does not keep the runner from starting |
| Does the self-test notice faults? | six copies of the final `queue.py` with one fault each | every fault makes checks fail: GPU limit ignored (3 checks fail), CPU limit ignored (1), dependencies ignored (5), orphans not ended (1), success test ignored (2), done not tested again at the start (1) |
| The runner under systemd | a transient unit (`systemd-run`) with the settings of the unit file (`KillMode=mixed`, `Restart=on-failure`), dummy tasks | runner killed with SIGKILL: systemd ended its task and started the runner again, which ran the task again; STOP file: runner ended with 0 and was **not** started again; `systemctl stop`: the task got SIGTERM; `systemd-analyze verify` of the unit file: no complaint |
| `tasks_initial.json` | `queue.py check` | 6 tasks accepted for 24 CPUs and 1 GPU |
| `preflight_gpu` without a GPU | run through the runner on the Pi as it is | fails closed: `nvidia-smi` not found, three attempts, then failed; the two tasks after it never start. With a stand-in `nvidia-smi` the OpenMM test runs (Reference and CPU platform) and the task fails on "the CUDA platform of OpenMM did not compute forces". Its JSON part was run on a text with a CUDA line added by hand |
| `preflight_speed`, `preflight_interrupt` | the commands unchanged, through the runner, with only the environment block changed: the 4,600-atom stand-in of `md/RUN.md` section 6, CPU platform with one thread, a test protocol, 15 s of timing | both done; `speed.json` written (2.69 ns per day, which says nothing about a GPU); the second run was killed during `eq05` after 54.3 s, 60 % of the first run's 90.5 s, resumed, and all 25 files were IDENTICAL |
| The ready-made run entries | X0, X1, N1, R1 with the stand-in; in the commands `--platform CUDA` and the lengths replaced by the CPU test options | with `test_only=false` X0 is refused after three attempts (it is a test run) and nothing after it starts; with that one condition changed for the test, X0, X1, N1 and R1 ran in this order |
| `setup_gpu.sh`: resolving the environments | `--solve-only` for linux-64 with a driver assumed to report CUDA 13.0, 13.4, 12.9 and 12.8 | the first three resolve (section 8); 12.8 stops with the message that a newer driver series is needed |
| `setup_gpu.sh`: the checks | `--checks-only` against the Pi's own two environments, without GPU | all pass: versions, tleap, MCPB.py, `run_md.py --help`, and the SCF of [Zn(NH₃)₄]²⁺ (−2004.9070316727 Hartree, the value the machine must reproduce within 10⁻⁵); the comparison with `md_env.yml` gives 242 of 242. With a GPU expected it fails closed, last line `SETUP_GPU RESULT: FAILED` |
| `push_to_vm.sh --yes`, after the review (6 Oct 2026, 9:20 AM IST) | a copy of the script in which only the four fixed paths at its top point into scratch, with a stand-in `ssh` that runs the remote commands on the Pi inside a second scratch folder, so rsync, mkdir and `sha256sum -c` are real: a first push, then a push after the "machine" had rewritten one file, the "Pi" another, and a third existed only on the "machine" | 13 of 13 checks: the file that is newer on the machine stays, is named in `push_kept_newer.txt` and left out of the list; the Pi's newer file arrives; an older file on the machine is replaced; nothing is deleted; VERIFIED and PUSH COMPLETE each time. The script as it was before fails 4 of the 13 (it puts the Pi's older file over the machine's) |
| `push_to_vm.sh --preview`, `deadman.py --state` and the switch step of `setup_gpu.sh`, after the second review (6 Oct 2026, 11:42 AM IST) | the same stand-ins, in `gpuops/fix2/tests/`: the preview on a "machine" that holds a file the Pi changed, one it rewrote itself, one of equal content with an older time and one that is older and different; the step with a `PI_SEEN` of 300 and of 3 minutes, with `NO_DEADMAN`, and with an older `deadman.py`; the switch with a status folder that cannot be written | `t_push.sh`: 18 of 18 checks, five of them for the preview (it names the two files that would be replaced, the one that stays and the one that would only get a time stamp, and changes nothing on either side; 15 of 18 before); `t_setup_deadman.sh`: 17 of 17, four of them new (with a `PI_SEEN` of 300 minutes nothing is installed and nothing enabled; 13 of 17 before); the switch: 27 of 27 as before, and 21 of 21 in the reviewer's set as extended (never seen by the Pi and with a status folder that cannot be written, it shuts down once the machine has been up for 240 minutes; before, it never did). **None of it has run on the machine or over a real connection** |
| `md/gpu/deadman.py`, after the review | 27 checks with a faked clock and faked commands (`gpuops/fix/tests/t_deadman.py`); the step of `setup_gpu.sh` that installs it, with stand-ins for `install` and `systemctl` (13 checks); `systemd-analyze verify` of the two unit files | all pass: nothing below 240 minutes; then `queue.py stop`, the wait for the runner (the shutdown comes when it has ended, or after 10 minutes), the STOP file removed, `shutdown -h now`; `NO_DEADMAN`; the dry run; and on the Pi (aarch64) `--yes` refuses, with a stand-in `shutdown` on the PATH that stayed untouched. **Never run on the machine** |
| `push_to_vm.sh` | the dry run; `--rehearse` (the checksum list is made and checked against the Pi's own files, no machine); the part that drops files changed during the copy, on made-up input | 459 files, 57.3 MB (md 8.8, systems 28.2, qm 10.5, ligand 0.3, system 9.6 MB), measured again at 6:55 AM IST while a 2 GB scratch file lay in the Pi's quantum-chemistry folder (left out); the list of 459 files verifies; a changed and a vanished file are taken out and named; `--yes` without an address, or with something that is not an address, stops with exit code 2 and writes nothing. An end of the dry run with exit code 141 (a pipe closed early), seen once, was removed; six runs in a row then ended with 0. The numbers are those of 6:27 and 6:55 AM IST; at 8:11 AM the dry run gave 693 files, 72.4 MB, at 9:20 AM 762 files, 86.1 MB (md 10.4, systems 28.2, qm 37.6, ligand 0.3, system 9.6 MB) |
| Do the watchdog and the task files fit? | `md_plans_of_task` of `ops/pi/cqc_gpu_watchdog.py` (as of 6:38 AM IST) on `tasks_initial.json` and on the ready-made entries | every `run_md.py` call is read: folder, run id, length, test or benchmark; none unresolved. The heartbeat has every key that file reads (`epoch`, `queue_state`, `idle`, `idle_since_utc`, `all_done`, `n_tasks`, `counts`, `running`, `attention`, `starts_at_boot`, `gpu`, `load_average`) |

**Not tested, because it cannot be tested here:**

1. Everything in `setup_gpu.sh` that changes a machine: apt, the NVIDIA packages, the module, the reboot path,
   the download of micromamba inside the script, the two `micromamba create` calls for real, gpu4pyscf.
2. The CUDA platform: `preflight_gpu` passing, any speed, bit-for-bit repetition on the GPU, the 60 % kill time
   falling inside the run on a fast machine (if the run was not killed, the task fails with a message that says so).
3. `push_to_vm.sh --yes` and `--verify` were not tested when this was written: no machine existed. Since then
   the script as it was ran once for real (6 Oct 2026, 8:45 AM IST, to the new machine: 723 files in the list,
   VERIFIED, PUSH COMPLETE; `/root/cqc-sfa-md/status/push_to_vm.log`). The version with `--update` has run only
   against the stand-in of the table above: its two extra dry runs per part over a real connection, and
   `sudo -n rsync` with them, have not.
4. The service started at boot through `systemctl enable` (the unit was not installed on the Pi), and what a
   spot pre-emption does to a running machine (a clean shutdown or a hard stop; the runner is built for both).
5. The quantum-chemistry tasks themselves. Their commands are the driver's own (`md/ZINC.md`, section 9); the
   optimisation and force-constant jobs of the driver have run only on the 16-atom complex (`md/ZINC.md`,
   section 10, point 9). How long they take on 16 threads is not known.
6. 16 PySCF threads beside a GPU run on 24 vCPUs: whether the MD loses speed.

## 10. Open points for whoever starts the machine

1. **Start R and the empty pocket**: section 5, last note. A decision, not a detail.
2. **The numerical cross-check of `md/ZINC.md` section 6.1** (finite differences against PySCF's analytic
   Hessian) is unresolved, and that section asks whoever runs the small model to repeat it first. `qm_hess`
   computes the finite-difference force constants as planned and does not repeat the check.
3. **`config_prod.json` allows PySCF 8,000 MB** (the Pi's value). The small model fits; the machine has 120 GB,
   and the value was left as it is so that the settings are those of the Pi.
4. **The machine idles when a task has failed** (section 1, points 4 and 5). The watchdog, as read, raises an
   alert after 10 minutes of that and stops the machine after 3 hours; until then it costs the full rate.
5. **Which model the pre-flight tests ran on** is in their result files (`system`); they are tests, their run
   folders carry `TEST_ONLY` (the interruption test) or are named `bench`.
6. **Automatic security upgrades are switched off on the machine** and driver and kernel are held, on purpose
   (`setup_gpu.sh`, header). The machine is meant to live for days and to accept SSH from the Pi only.
7. **PySCF's scratch files do not all begin with `tmp`.** On the Pi, at 6:51 AM IST, the gradient of the `timing`
   job had a scratch file `t71ug9bl` of 1.98 GB in `/root/cqc-sfa-md/qm/prod/`. `push_to_vm.sh` therefore takes
   from the quantum-chemistry folder only files that have an extension. **The watchdog now has the matching
   rule** (read on 6 Oct 2026, 9:00 AM IST: `SCRATCH_QM_RE` in `ops/pi/cqc_gpu_watchdog.py`; an earlier version
   of this point said it had not): under `qm/` on the machine, and under `md/qm/` there, a file whose name is
   eight lower-case letters, digits or underscores without an extension, with or without `tmp` before them, is
   taken for scratch and not pulled; one above 100 MB that is left out for its name is listed in
   `gpu_status.json`, in the status table and in the alert of the end. So **a result under `qm/` must have an
   extension**, or it stays on the machine.
8. **The switch on the machine has never run on a machine** (section 4), and the time limits of `qm_opt` and
   `qm_hess` are a guess (section 3).
