#!/usr/bin/env bash
# pi_test_after_stop.sh: what could still be done after run X1 was stopped from outside.
#
# TEST OF THE PIPELINE with the restrained-zinc model. NOT a result for the paper. See md/test/PI_TEST.md.
#
# On 6 Oct 2026, 8:04 AM IST, the main session killed run X1 inside its first production segment and
# stopped the chain (pi_test_chain.sh): /root/cqc-sfa-md/runs_test/STOPPED_README.txt. The run was not
# continued. This script works on what is on disk, with one thread and in a few minutes:
#
#   1. the state of the run folder as the kill left it (pi_test_check.py snapshot);
#   2. the recovery of the driver on COPIES of the folder: "--stage min --resume" makes the driver take the
#      lock of the dead process, verify the last saved state, remove the files of the cut-off segment and
#      stop before it would repeat the segment (which would take about 4 hours here). Once with the driver
#      and protocol the run was made with (git a98bd3f) and once with the driver as it is now;
#   3. frame counts, checksums and time axis of the complete steps, and what the cut-off step holds
#      (pi_test_check.py check --incomplete);
#   4. md/analysis/analyse_md.py, as it is now, on the unfinished run;
#   5. every snapshot from the built start to the last frame written (pi_test_look.py), also for the
#      minimised start R of md/TOPOLOGY.md (its test protocol, 300 iterations per pass; no run R1 was made);
#   6. speed and the processor time received (pi_test_report.py).
# Running it again gives the same files (the texts carry the time they were written).
set -u
REPO=${REPO:-/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate}
MD=$REPO/md
SYS=${SYS:-/root/cqc-sfa-md/systems}
RUNS=${RUNS:-/root/cqc-sfa-md/runs_test}
PY=${PY:-/root/cqc-sfa-tools/md/bin/python}
RES=$MD/test/results
run=$RUNS/restrained/X1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENMM_CPU_THREADS=1
mkdir -p "$RES"
status=0

# 1. the folder as the kill left it
nice -n 10 "$PY" "$MD/test/pi_test_check.py" snapshot "$run" --killed-pid "$(cat "$RUNS/restrained/X1.pid")" \
    --signal "SIGKILL, sent by the main session (STOPPED_README.txt), not by this test" \
    --out "$RES/X1_interrupt_snapshot.json" > "$RES/X1_interrupt_snapshot.console" 2>&1 || status=1
cat "$RES/X1_interrupt_snapshot.console"

# 2. recovery on copies
old=$RUNS/driver_of_the_run
mkdir -p "$old"
git -C "$REPO" show a98bd3f:md/run_md.py > "$old/run_md.py"
git -C "$REPO" show a98bd3f:md/protocol.json > "$old/protocol.json"
for which in driver_of_the_run driver_now; do
  copy=$RUNS/resume_check/$which/X1
  rm -rf "$RUNS/resume_check/$which"
  mkdir -p "$copy"
  cp -p "$run"/manifest.json "$run"/rec_*.json "$run"/state_eq06.xml "$run"/*_seg0001.* "$run"/run_md.lock \
        "$run"/solute.pdb "$run"/full.pdb "$run"/TEST_ONLY "$copy"/ 2>/dev/null
  if [ "$which" = driver_of_the_run ]; then drv="$old/run_md.py"; proto="--protocol $old/protocol.json"; else drv="$MD/run_md.py"; proto=""; fi
  { echo "driver: $drv  sha256 $(sha256sum "$drv" | cut -d' ' -f1)"
    echo "files of the cut-off segment before: $(cd "$copy" && ls -1 *_seg0001.* 2>/dev/null | tr '\n' ' ')"
    nice -n 10 "$PY" "$drv" --system "$SYS/restrained/X" --out "$copy" --run-id X1 --platform CPU --threads 1 \
        --stage min --equil-scale 0.05 --prod-ns 0.2 --segment-ns 0.1 --resume $proto
    echo "exit code $?"
    echo "files of the cut-off segment after: $(cd "$copy" && ls -1 *_seg0001.* 2>/dev/null | tr '\n' ' ')"
    echo "lock file after: $([ -f "$copy/run_md.lock" ] && echo present || echo absent)"
  } > "$RES/X1_resume_check_$which.console" 2>&1
  cat "$RES/X1_resume_check_$which.console"
done
nice -n 10 "$PY" "$MD/test/pi_test_check.py" recovery "$RUNS/resume_check" --snapshot "$RES/X1_interrupt_snapshot.json" \
    --out "$RES/X1_resume_check.json" || status=1

# 3. the files of the run
nice -n 10 "$PY" "$MD/test/pi_test_check.py" check "$run" --incomplete --out "$RES/X1_check.json" > "$RES/X1_check.console" 2>&1 || status=1
tail -n 1 "$RES/X1_check.console"

# 4. the analysis program on the unfinished run
t0=$(date +%s)
nice -n 10 "$PY" "$MD/analysis/analyse_md.py" "$run" > "$RES/X1_analyse_md.console" 2>&1; code=$?
echo "analyse_md.py X1: exit code $code, $(( $(date +%s) - t0 )) s; script sha256 $(sha256sum "$MD/analysis/analyse_md.py" | cut -d' ' -f1)" | tee -a "$RES/X1_analyse_md.console"
for f in summary.json timeseries_solute.csv timeseries_full.csv fig_zinc.png fig_ligand.png fig_orientation.png; do
  [ -f "$run/analysis/$f" ] && cp "$run/analysis/$f" "$RES/X1_analysis_$f"
done

# 5. every snapshot
nice -n 10 "$PY" "$MD/test/pi_test_look.py" "$run" --out-prefix "$RES/X1" || status=1
rmin=$SYS/_checks/restrained/R/min_300
[ -f "$rmin/manifest.json" ] && { nice -n 10 "$PY" "$MD/test/pi_test_look.py" "$rmin" --out-prefix "$RES/R_min300_of_TOPOLOGY" --no-figure || status=1; }
# the dense-output test of start N of md/TOPOLOGY.md: its state log has a row at every barostat attempt (used in step 6)
nmin=$SYS/_checks/restrained/N/min_300
[ -f "$nmin/manifest.json" ] && { nice -n 10 "$PY" "$MD/test/pi_test_look.py" "$nmin" --out-prefix "$RES/N_min300_of_TOPOLOGY" --no-figure || status=1; }

# 6. speed
nice -n 10 "$PY" "$MD/test/pi_test_report.py" --runs "$RUNS" --out "$RES/pi_test_speed.json" > "$RES/pi_test_speed.console" 2>&1 || status=1
tail -n 1 "$RES/pi_test_speed.console"
exit $status
