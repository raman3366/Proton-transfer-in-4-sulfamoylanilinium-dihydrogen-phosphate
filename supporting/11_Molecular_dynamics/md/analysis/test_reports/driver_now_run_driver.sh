#!/bin/bash
# throwaway (version 1.3): two short runs of the current md/run_md.py on the stand-in systems of its own tests
# (/root/cqc-sfa-md/test_driver/hydration/systems/dry and wet, read only), with the hydration block of md/protocol.json
# as it stands.  One thread, nice 10.  Output only below this folder.
set -u
export OMP_NUM_THREADS=1 OPENMM_CPU_THREADS=1
cd /tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/mdgr2/spec13/driver_now
echo "== started $(TZ=Asia/Kolkata date '+%d %b %Y, %I:%M:%S %p IST')"
sha256sum /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md/run_md.py /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md/protocol.json protocol_test_now.json
for s in dry wet; do
  rm -rf runs/$s/R1
  nice -n 10 /root/cqc-sfa-tools/md/bin/python /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md/run_md.py --system /root/cqc-sfa-md/test_driver/hydration/systems/$s --out /tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/mdgr2/spec13/driver_now/runs/$s/R1 --run-id R1 \
      --platform CPU --threads 1 --equil-scale 0.001 --segment-ns 0.0004 --prod-ns 0.0004 --protocol /tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/mdgr2/spec13/driver_now/protocol_test_now.json > run_$s.console 2>&1
  echo "$s: exit code $?; $(tail -n 1 run_$s.console)"
  grep "hydration check" run_$s.console | sed 's/^/   /'
done
echo "== finished $(TZ=Asia/Kolkata date '+%d %b %Y, %I:%M:%S %p IST')"
echo DRIVERDONE
