#!/bin/bash
# Start md/qm/run_qm.py detached, with 4 threads and low priority.
#
#   md/qm/run_qm.sh [jobs]          jobs: comma-separated, default as in md/qm/config_prod.json
#
#   md/qm/run_qm.sh                 what was started on 6 Oct 2026: esp,write,timing
#   md/qm/run_qm.sh opt,hess,write  the small model (1 to 3.5 days on this Pi; md/ZINC.md, section 7)
#
# The driver is restartable: run the same command again after a stop; finished jobs are skipped, the
# optimisation continues from its last geometry, the force constants from the last displaced gradient.
# Progress:  cat /root/cqc-sfa-md/qm/prod/status.json ;  tail /root/cqc-sfa-md/qm/prod/run_qm.log
# Stop:      kill $(cat /root/cqc-sfa-md/qm/prod/driver.pid)        (by PID; never pkill -f)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
BULK=/root/cqc-sfa-md/qm/prod
mkdir -p "$BULK"
if [ -s "$BULK/driver.pid" ] && kill -0 "$(cat "$BULK/driver.pid")" 2>/dev/null; then
  echo "a driver is running already, PID $(cat "$BULK/driver.pid")" >&2; exit 1
fi
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
ARGS=("$HERE/config_prod.json")
[ $# -ge 1 ] && ARGS+=(--jobs "$1")
cd "$HERE/../.."
rm -f "$BULK/driver.pid"
setsid nohup nice -n 10 /root/cqc-sfa-tools/env/bin/python "$HERE/run_qm.py" "${ARGS[@]}" >> "$BULK/run_qm.stdout" 2>&1 &
# the driver writes its own PID into driver.pid and into status.json when it starts
until [ -s "$BULK/driver.pid" ]; do sleep 1; done
echo "started, PID $(cat "$BULK/driver.pid"), $(date '+%Y-%m-%d %I:%M %p %Z')"
