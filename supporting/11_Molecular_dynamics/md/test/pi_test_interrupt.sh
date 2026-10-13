#!/usr/bin/env bash
# pi_test_interrupt.sh: interrupt run X1 of the pipeline test once, for real, and resume it.
#
# TEST OF THE PIPELINE with the restrained-zinc model. NOT a result for the paper. See md/test/PI_TEST.md.
#
# Waits until the second production segment of X1 has written at least two state-log rows
# (20 ps into the segment), kills run_md.py by its process id with SIGKILL (as a lost machine
# would), records what the kill left behind (pi_test_check.py snapshot; the files of the killed
# attempt are copied to <runs>/restrained/X1_killed_attempt/), and continues the run with the same
# command and --resume (pi_test_run.sh X).
#
#   bash md/test/pi_test_interrupt.sh          # start it while pi_test_run.sh X is running
#
# It does nothing if the segment is already complete or if it has already been run
# (results/X1_interrupt_snapshot.json exists).
set -u
MD=${MD:-/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md}
RUNS=${RUNS:-/root/cqc-sfa-md/runs_test}
PY=${PY:-/root/cqc-sfa-tools/md/bin/python}
ROWS=${ROWS:-2}
run=$RUNS/restrained/X1
snap=$MD/test/results/X1_interrupt_snapshot.json
now() { TZ=Asia/Kolkata date '+%-d %b %Y, %-I:%M:%S %p IST'; }
mkdir -p "$MD/test/results"
if [ -f "$snap" ]; then echo "$(now) the interruption was already made ($snap); nothing done"; exit 0; fi
echo "$(now) waiting for $ROWS state-log rows of seg0002 in $run"
while true; do
  if [ -f "$run/rec_seg0002.json" ]; then echo "$(now) seg0002 is already complete; no interruption made"; exit 1; fi
  if [ -f "$run/log_seg0002.csv" ]; then
    rows=$(( $(wc -l < "$run/log_seg0002.csv") - 1 ))
    [ "$rows" -ge "$ROWS" ] && break
  fi
  if ! kill -0 "$(cat "$RUNS/restrained/X1.pid")" 2>/dev/null; then
    echo "$(now) run_md.py of X1 is not running and seg0002 has not got that far; no interruption made"; exit 2
  fi
  sleep 15
done
pid=$(cat "$RUNS/restrained/X1.pid")
if ! tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null | grep -q "run_md.py.*--run-id X1"; then
  echo "$(now) process $pid is not the run_md.py of X1; nothing killed"; exit 1
fi
echo "$(now) seg0002 has $rows log rows; sending SIGKILL to run_md.py, pid $pid; load $(cut -d' ' -f1-3 /proc/loadavg)"
kill -KILL "$pid"
while kill -0 "$pid" 2>/dev/null; do sleep 1; done
sleep 3
echo "$(now) pid $pid is gone; recording what it left"
# keep the files of the killed attempt (a few frames): --resume removes them from the run folder
kept=$RUNS/restrained/X1_killed_attempt
mkdir -p "$kept"
cp -p "$run"/*_seg0002.* "$kept"/ 2>/dev/null
"$PY" "$MD/test/pi_test_check.py" snapshot "$run" --killed-pid "$pid" --signal SIGKILL --out "$snap" || exit 1
echo "$(now) resuming with pi_test_run.sh X"
bash "$MD/test/pi_test_run.sh" X
code=$?
echo "$(now) the resumed run ended with exit code $code"
exit $code
