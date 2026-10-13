#!/bin/bash
# throwaway driver of the synthetic test (version 1.3): generate, analyse (the campaign, the validation runs one
# by one, the other single runs, the runs of versions 1.2 and 1.3 with their equilibration files), compare.
# One thread, niced.   usage: run_test.sh [--keep]     A (folder of the two scripts) and T (test folder) from the environment
set -u
export OMP_NUM_THREADS=1 AMBERHOME=/root/cqc-sfa-tools/md
PY=/root/cqc-sfa-tools/md/bin/python
A=${A:-/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md/analysis}
T=${T:-/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/mdgr2/spec13/test}
if [ "${1:-}" != "--keep" ]; then
  rm -rf $T
  nice -n 10 $PY $A/make_test_trajectory.py --out $T || exit 1
fi
nice -n 10 $PY $A/analyse_md.py --campaign $T/runs || exit 1
nice -n 10 $PY $A/analyse_md.py $T/runs/test/X0 || exit 1          # the validation run: on its own, no option
nice -n 10 $PY $A/analyse_md.py $T/runs_overlap/test/X9
for d in $T/runs_gv/*/*/ $T/runs_v12/*/*/ $T/runs_v13/*/*/; do nice -n 10 $PY $A/analyse_md.py $d || exit 1; done
nice -n 10 $PY $A/make_test_trajectory.py --compare $T
echo "exit code of the comparison: $?"
