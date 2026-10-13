#!/usr/bin/env bash
# pi_test_cpu_sampler.sh: once a minute, how much processor time the test runs really got.
#
# TEST OF THE PIPELINE with the restrained-zinc model. NOT a result for the paper. See md/test/PI_TEST.md.
#
# The Pi is shared, so "2 threads" does not mean two cores. This writes, every 60 s, the processor time
# used so far by each run_md.py named in $RUNS/restrained/*.pid and by the whole machine, and the load
# average. pi_test_report.py turns the samples into "cores received" for every step of a run.
#
#   bash md/test/pi_test_cpu_sampler.sh &        # ends by itself when $RUNS/cpu_sampler.stop exists
#
# Columns: epoch_s, run, pid, process_cpu_s (user + system, all threads, since the process began),
# machine_busy_s and machine_idle_s (all cores, since boot; busy = user + nice + system + irq + softirq,
# idle = idle + iowait), load_1min, load_5min, load_15min, n_cores.
set -u
RUNS=${RUNS:-/root/cqc-sfa-md/runs_test}
out=$RUNS/cpu_samples.csv
tck=$(getconf CLK_TCK)
ncores=$(nproc)
[ -f "$out" ] || echo "epoch_s,run,pid,process_cpu_s,machine_busy_s,machine_idle_s,load_1min,load_5min,load_15min,n_cores" > "$out"
while [ ! -f "$RUNS/cpu_sampler.stop" ]; do
  read -r _ u n s i w q sq _ < /proc/stat
  busy=$(awk -v a="$u" -v b="$n" -v c="$s" -v d="$q" -v e="$sq" -v t="$tck" 'BEGIN{printf "%.2f",(a+b+c+d+e)/t}')
  idle=$(awk -v a="$i" -v b="$w" -v t="$tck" 'BEGIN{printf "%.2f",(a+b)/t}')
  read -r l1 l5 l15 _ < /proc/loadavg
  now=$(date +%s)
  for f in "$RUNS"/restrained/*.pid; do
    [ -f "$f" ] || continue
    pid=$(cat "$f"); run=$(basename "$f" .pid)
    if [ -r "/proc/$pid/stat" ] && tr '\0' ' ' < "/proc/$pid/cmdline" | grep -q "run_md.py"; then
      # fields 14 and 15 of /proc/<pid>/stat are utime and stime in clock ticks (the name in field 2 has no blank here)
      cpu=$(awk -v t="$tck" '{printf "%.2f",($14+$15)/t}' "/proc/$pid/stat")
      echo "$now,$run,$pid,$cpu,$busy,$idle,$l1,$l5,$l15,$ncores" >> "$out"
    fi
  done
  echo "$now,machine,0,0,$busy,$idle,$l1,$l5,$l15,$ncores" >> "$out"
  sleep 60 &
  wait $!
done
