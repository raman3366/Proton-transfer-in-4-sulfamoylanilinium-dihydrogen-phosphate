#!/bin/bash
# Gate of stage S6: the one documented command (md/analysis/RECHECK.md, section 1; md/PRODUCTION.md, "After the runs").
# Environment as in the gate checks of 6 Oct (one thread per process, no display). Output and exit status are kept.
set -u
L=/root/cqc-sfa-md/analysis_logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg TZ=Asia/Kolkata
cd /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate || exit 99
echo "start $(date '+%d %b %Y, %I:%M:%S %p IST') epoch $(date +%s)" > $L/gate_S6.times
nice -n 10 /root/cqc-sfa-tools/md/bin/python md/analysis/recheck_md.py --compare-campaign /root/cqc-sfa-md/runs --run-primary > $L/gate_S6.log 2>&1
st=$?
echo "end $(date '+%d %b %Y, %I:%M:%S %p IST') epoch $(date +%s)" >> $L/gate_S6.times
echo "$st" > $L/gate_S6.exit_status
echo "exit status $st" >> $L/gate_S6.times
