#!/bin/bash
# SANITISED COPY (SAPH hand-over package): a network address was masked.
# The original is in the git repository; the key itself never was.
# pull_run.sh RUN: copy one finished rehearsal run from the machine into the scratch folder (read-only on the machine)
set -euo pipefail
run=$1
rsync -az --partial --stats -e "ssh -i /root/.ssh/aitf_pi -o IdentitiesOnly=yes -o BatchMode=yes -o ConnectTimeout=15 -o ServerAliveInterval=30" --rsync-path="sudo -n rsync" --exclude='*.chk' --exclude='*.tmp' --exclude=run_md.lock \
  ubuntu@<address of the rented GPU machine, which was deleted on 8 Oct 2026>:/root/cqc-sfa-md/runs_test/bonded_2e/$run/ /tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/rehearsal/runs_test/bonded_2e/$run/ | grep -E 'Number of regular files transferred|Total transferred file size|Total bytes received'
