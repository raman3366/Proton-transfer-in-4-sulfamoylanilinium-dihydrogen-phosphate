#!/bin/bash
# one process at a time, niced; the eleven 100 ns runs, then the two validation runs for comparison
cd "$(dirname "$0")"
PY=/root/cqc-sfa-tools/md/bin/python
for mr in "bonded_2e R1" "bonded_2e R2" "bonded_2e R3" "bonded_2s R1" "bonded_2e_n2fit R1" "bonded_2e X1" "bonded_2e X2" "bonded_2e X3" "bonded_2e N1" "bonded_2e N2" "bonded_2e N3" "bonded_2s X0"; do
  set -- $mr
  nice -n 10 $PY extract.py $1 $2 2>/dev/null >> extract.log
done
echo ALL_DONE >> extract.log
