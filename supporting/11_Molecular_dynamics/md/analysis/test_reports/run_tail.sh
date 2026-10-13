#!/bin/bash
# throwaway (version 1.3): the end of run_all.sh after its loop over the broken copies was stopped (the machine was too
# loaded: 16 copies in 44 minutes): the remaining copies in stages (mutation/mutate_staged.py), then the last lines of
# run_all.sh as they stand there.  Appends to run_all.log.  One thread, niced.
set -u
export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
S=/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/mdgr2/spec13; A=${A:-/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md/analysis}; PY=/root/cqc-sfa-tools/md/bin/python
cd $S
{
echo
echo "== mutations, continued in stages $(TZ=Asia/Kolkata date '+%d %b %Y, %I:%M %p IST'): the loop above was stopped after M16 (the machine was too loaded to run every copy on all runs); run_tail.sh, mutation/mutate_staged.py"
nice -n 10 $PY mutation/mutate_staged.py $S $A $PY
echo "== scripts again"; sha256sum $A/analyse_md.py $A/make_test_trajectory.py; echo "driver and protocol:"; sha256sum $A/../run_md.py $A/../protocol.json
echo "== runs of the campaign on this machine (at the end)"; $PY runs_state.py; $PY runs_dates.py $S/runs_dates.json > /dev/null
echo "== finished $(TZ=Asia/Kolkata date '+%d %b %Y, %I:%M %p IST')"
echo ALLDONE
} >> run_all.log 2>&1
