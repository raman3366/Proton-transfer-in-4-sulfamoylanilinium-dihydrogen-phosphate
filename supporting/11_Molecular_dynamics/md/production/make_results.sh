#!/bin/bash
# Everything of stage S6 that follows the gate, in its order. The gate itself is the one command of
# md/analysis/RECHECK.md (wrapper: /root/cqc-sfa-md/analysis_logs/run_gate_S6.sh, copy in md/results/gate/);
# it must have ended with exit status 0 before this is run. Nothing here changes a file under /root/cqc-sfa-md/runs.
#   1. the test of the conversion-time program on made-up series (must pass)
#   2. md/results/conversion_times.json
#   3. md/results/fig_orientation.(png|pdf), fig_ligand_rmsd.(png|pdf), fig_sources.json
#   4. the copies of the result files, run_health.json, run_events.json, SOURCES.json
#   5. the checks on the raw frames that were added after the audits of 8 Oct 2026 (md/production/site_checks.py):
#      md/results/site_checks_solute.json, site_checks_full.json (one process at a time)
#   6. md/RESULTS.md and md/results/tables.md
# md/results/independent/ (the physical review and the blind third reading) is a copy made once, by hand, from
# the scratch folder of the session of 8 Oct 2026; its checksums are in md/results/independent/SOURCES.json.
set -eu
set -o pipefail      # every step is piped to tee: without this a failed step would not stop the script
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg TZ=Asia/Kolkata
PY=/root/cqc-sfa-tools/md/bin/python
L=/root/cqc-sfa-md/analysis_logs
cd /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate
[ "$(cat $L/gate_S6.exit_status)" = "0" ] || { echo "the gate of stage S6 has not passed"; exit 1; }
mkdir -p md/results
nice -n 10 $PY md/production/test_conversion_time.py   2>&1 | tee $L/test_conversion_time.log
cp $L/test_conversion_time.log md/results/conversion_time_test.log
nice -n 10 $PY md/production/conversion_time.py        2>&1 | tee $L/conversion_time.log
nice -n 10 $PY md/production/make_result_figures.py    2>&1 | tee $L/make_result_figures.log
nice -n 10 $PY md/production/collect_results.py        2>&1 | tee $L/collect_results.log
nice -n 10 $PY -W ignore md/production/site_checks.py --part solute 2>&1 | tee $L/site_checks_solute.log
nice -n 10 $PY -W ignore md/production/site_checks.py --part full   2>&1 | tee $L/site_checks_full.log
nice -n 10 $PY md/production/assemble/make_results_md.py 2>&1 | tee $L/make_results_md.log
