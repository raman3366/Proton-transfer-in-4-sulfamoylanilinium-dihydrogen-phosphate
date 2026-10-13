#!/usr/bin/env bash
# pi_test_bench.sh: speed of the equilibrated X system with 1 and with 2 threads, one after the other.
#
# TEST OF THE PIPELINE with the restrained-zinc model. NOT a result for the paper. See md/test/PI_TEST.md.
#
# The driver's own speed test (run_md.py --benchmark-seconds: NPT at 4 fs from the saved state of the last
# equilibration stage, no trajectory written) in copies of the X1 folder that hold only what the speed test
# needs (manifest, records, the last equilibration state). Beside the ns per day it records the processor
# time the process used, so that the two thread counts can be compared on a shared machine.
#
#   bash md/test/pi_test_bench.sh        # needs the finished equilibration of X1; about 8 minutes
#
# Each call adds one entry per thread count to <runs>/bench_<n>thr/benchmark.json and bench_cpu.json.
set -u
MD=${MD:-/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md}
SYS=${SYS:-/root/cqc-sfa-md/systems}
RUNS=${RUNS:-/root/cqc-sfa-md/runs_test}
PY=${PY:-/root/cqc-sfa-tools/md/bin/python}
SECONDS_EACH=${SECONDS_EACH:-180}
src=$RUNS/restrained/X1
[ -f "$src/rec_eq06.json" ] || { echo "X1 has no finished equilibration"; exit 1; }
for n in 1 2; do
  dst=$RUNS/bench_${n}thr
  if [ ! -f "$dst/manifest.json" ]; then
    mkdir -p "$dst"
    cp -p "$src"/manifest.json "$src"/rec_min.json "$src"/rec_eq0?.json "$src"/state_eq06.xml "$dst"/
  fi
  export OMP_NUM_THREADS=$n OPENMM_CPU_THREADS=$n
  load0=$(cut -d' ' -f1-3 /proc/loadavg)
  read -r _ u ni s i w q sq _ < /proc/stat; busy0=$((u+ni+s+q+sq)); idle0=$((i+w))
  t0=$(date +%s.%N)
  out=$( { /usr/bin/time -f "TIME user %U system %S elapsed %e" nice -n 10 "$PY" "$MD/run_md.py" --system "$SYS/restrained/X" \
        --out "$dst" --run-id X1 --platform CPU --threads "$n" --equil-scale 0.05 --prod-ns 0.2 --segment-ns 0.1 \
        --benchmark-seconds "$SECONDS_EACH"; } 2>&1 )
  code=$?
  t1=$(date +%s.%N)
  read -r _ u ni s i w q sq _ < /proc/stat; busy1=$((u+ni+s+q+sq)); idle1=$((i+w))
  echo "$out" | tail -n 4
  echo "$out" | "$PY" -c "
import json, re, sys
text = sys.stdin.read()
m = re.search(r'TIME user ([\d.]+) system ([\d.]+) elapsed ([\d.]+)', text)
b = re.search(r'BENCHMARK \S+ ([\d.]+) ns/day .*?(\d+) steps in ([\d.]+) s', text)
tck = 100.0
out = {'exit_code': $code, 'threads_requested': $n, 'seconds_requested': $SECONDS_EACH,
       'load_average_before_1_5_15_min': '$load0'.split(),
       'whole_process_user_s': float(m.group(1)) if m else None, 'whole_process_system_s': float(m.group(2)) if m else None,
       'whole_process_elapsed_s': float(m.group(3)) if m else None,
       'whole_process_cores': round((float(m.group(1)) + float(m.group(2))) / float(m.group(3)), 3) if m else None,
       'machine_busy_cores_during_the_call': round(($busy1 - $busy0) / tck / ($t1 - $t0), 3),
       'machine_idle_cores_during_the_call': round(($idle1 - $idle0) / tck / ($t1 - $t0), 3),
       'note': 'whole process = start-up, building the System and the timed part; the cores of the timed part alone are not separated'}
json.dump(out, open('$dst/bench_cpu.json', 'w'), indent=1)
print(out)
"
done
