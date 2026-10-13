#!/usr/bin/env bash
# AD4Zn redocking and cross-docking for the two hCA II targets.
# The maps must already exist (dock/run_ad4zn.py --target T --stage maps --vm); they are receptor
# properties and are reused by every ligand.
set -u
cd "$(dirname "$0")/.."
PY=/root/cqc-sfa-tools/dock/bin/python
export PATH=/root/cqc-sfa-tools/dock/bin:$PATH

run() {
  echo "=== $1 $2 ad4zn"
  nice -n 10 "$PY" dock/run_ad4zn.py --target "$1" --ligand "$2" --stage dock \
      --seeds 10 --exhaustiveness 32 --cpu 3 2>&1 | grep -v 'tagged as 2D'
}

# redocking: the species that actually binds zinc is the sulfonamidate anion, but AD4Zn was
# calibrated with neutral -NH2 (it drops the Zn-H term), so both are run and both are reported.
run T3_6RL9 sfa_neutral
run T3_6RL9 sfa_anion
run T4_3HS4 azm_neutral
run T4_3HS4 azm_anion
# cross-docking
run T3_6RL9 azm_neutral
run T4_3HS4 sfa_neutral
echo "AD4ZN DONE"
