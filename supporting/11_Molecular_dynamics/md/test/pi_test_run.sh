#!/usr/bin/env bash
# pi_test_run.sh: start or continue one run of the pipeline test on the Raspberry Pi.
#
# TEST OF THE PIPELINE with the restrained-zinc model (nonbonded Zn2+, four Zn-N distance
# restraints). NOT a result for the paper. See md/test/PI_TEST.md.
#
#   bash md/test/pi_test_run.sh X     # equilibration x 0.05 (100 ps), production 200 ps in two segments of 100 ps
#   bash md/test/pi_test_run.sh R     # equilibration x 0.03 ( 60 ps), production 100 ps in one segment (a segment must be a
#                                     # whole multiple of the 100 ps interval of the full frames, so 2 x 50 ps is refused)
#
# The production protocol md/protocol.json is used unchanged (minimisation with its own limit,
# solute frames every 10 ps, full frames every 100 ps); only the lengths differ.
# Safe to run again: the driver continues after its last complete step (--resume).
# CPU platform, nice 10, THREADS threads (default 2, the limit of this test on the shared machine);
# each step record names its thread count.
# The process id is written to <run>.pid; the console to <run>.console.
set -u
start=${1:?start: X or R}
MD=${MD:-/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md}
SYS=${SYS:-/root/cqc-sfa-md/systems}
RUNS=${RUNS:-/root/cqc-sfa-md/runs_test}
PY=${PY:-/root/cqc-sfa-tools/md/bin/python}
THREADS=${THREADS:-2}
MODEL=restrained
case "$start" in
  X) scale=0.05; prod=0.2; seg=0.1 ;;
  R) scale=0.03; prod=0.1; seg=0.1 ;;
  *) echo "start must be X or R"; exit 2 ;;
esac
id=${start}1
mkdir -p "$RUNS/$MODEL"
export OMP_NUM_THREADS=$THREADS OPENMM_CPU_THREADS=$THREADS
echo "$(TZ=Asia/Kolkata date '+%-d %b %Y, %-I:%M:%S %p IST') pi_test_run.sh $start: $THREADS thread(s); load $(cut -d' ' -f1-3 /proc/loadavg)" >> "$RUNS/$MODEL/$id.console"
nice -n 10 "$PY" "$MD/run_md.py" --system "$SYS/$MODEL/$start" --out "$RUNS/$MODEL/$id" --run-id "$id" \
    --platform CPU --threads "$THREADS" --stage all --equil-scale "$scale" --prod-ns "$prod" --segment-ns "$seg" \
    --resume >> "$RUNS/$MODEL/$id.console" 2>&1 &
pid=$!
echo "$pid" > "$RUNS/$MODEL/$id.pid"
wait "$pid"
code=$?
echo "$(TZ=Asia/Kolkata date '+%-d %b %Y, %-I:%M:%S %p IST') pi_test_run.sh $start: run_md.py (pid $pid) ended with exit code $code; load $(cut -d' ' -f1-3 /proc/loadavg)" >> "$RUNS/$MODEL/$id.console"
exit $code
