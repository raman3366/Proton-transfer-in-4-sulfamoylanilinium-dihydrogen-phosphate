#!/usr/bin/env bash
# pi_test_chain.sh: the whole pipeline test, one thing after another, never more than 2 threads.
#
# TEST OF THE PIPELINE with the restrained-zinc model. NOT a result for the paper. See md/test/PI_TEST.md.
#
#   1. run X1 (pi_test_run.sh X; started here unless it is already running)
#   2. the interruption inside its second segment and the resumption (pi_test_interrupt.sh)
#   3. run R1 (pi_test_run.sh R)
#   4. speed test with 1 and with 2 threads from the equilibrated X1 (pi_test_bench.sh)
#   5. analysis and checks of both runs (pi_test_analyse.sh)
# Safe to start again: every part continues or answers that it is done.
set -u
MD=${MD:-/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md}
RUNS=${RUNS:-/root/cqc-sfa-md/runs_test}
now() { TZ=Asia/Kolkata date '+%-d %b %Y, %-I:%M:%S %p IST'; }
alive() { [ -f "$RUNS/restrained/$1.pid" ] && tr '\0' ' ' < "/proc/$(cat "$RUNS/restrained/$1.pid")/cmdline" 2>/dev/null | grep -q "run_md.py.*--run-id $1"; }
if [ ! -f "$MD/test/results/X1_interrupt_snapshot.json" ]; then
  if ! alive X1; then
    echo "$(now) chain: starting X1"
    setsid nohup bash "$MD/test/pi_test_run.sh" X > "$RUNS/restrained/X1.launcher.out" 2>&1 < /dev/null &
    sleep 20
  fi
  echo "$(now) chain: waiting to interrupt X1"
  bash "$MD/test/pi_test_interrupt.sh"; code=$?
else
  echo "$(now) chain: X1 was interrupted before; continuing it"
  while alive X1; do sleep 30; done
  bash "$MD/test/pi_test_run.sh" X; code=$?
fi
echo "$(now) chain: X1 ended with exit code $code"
[ $code -eq 0 ] || exit $code
echo "$(now) chain: starting R1"
bash "$MD/test/pi_test_run.sh" R; code=$?
echo "$(now) chain: R1 ended with exit code $code"
[ $code -eq 0 ] || exit $code
bash "$MD/test/pi_test_bench.sh"; echo "$(now) chain: speed test ended with exit code $?"
bash "$MD/test/pi_test_analyse.sh"; code=$?
echo "$(now) chain: analysis ended with exit code $code"
exit $code
