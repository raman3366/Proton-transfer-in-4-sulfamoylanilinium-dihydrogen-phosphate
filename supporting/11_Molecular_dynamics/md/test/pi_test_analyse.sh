#!/usr/bin/env bash
# pi_test_analyse.sh: analysis and checks of the two test runs. One thread.
#
# TEST OF THE PIPELINE with the restrained-zinc model. NOT a result for the paper. See md/test/PI_TEST.md.
#
#   1. md/analysis/analyse_md.py on X1 and R1 (writes <run>/analysis/; summary and tables are copied to md/test/results/)
#   2. pi_test_check.py: frame counts, time axis, checksums, the interruption of X1
#   3. pi_test_look.py: every snapshot from the built start to the last frame, and the comparison with analyse_md.py
#   4. analyse_md.py --campaign on the folder of the two runs (exercises that path; two runs decide nothing)
#   5. pi_test_report.py: speed, cores received, estimate for a validation run on the Pi
# Running it again gives the same files (the texts carry the time they were written).
set -u
MD=${MD:-/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md}
RUNS=${RUNS:-/root/cqc-sfa-md/runs_test}
PY=${PY:-/root/cqc-sfa-tools/md/bin/python}
RES=$MD/test/results
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$RES"
status=0
for id in X1 R1; do
  run=$RUNS/restrained/$id
  [ -f "$run/manifest.json" ] || { echo "$id: no run folder"; status=1; continue; }
  t0=$(date +%s)
  nice -n 10 "$PY" "$MD/analysis/analyse_md.py" "$run" > "$RES/${id}_analyse_md.console" 2>&1; code=$?
  echo "analyse_md.py $id: exit code $code, $(( $(date +%s) - t0 )) s" | tee -a "$RES/${id}_analyse_md.console"
  [ $code -eq 0 ] || status=1
  for f in summary.json timeseries_solute.csv timeseries_full.csv rmsf_ca.csv; do
    [ -f "$run/analysis/$f" ] && cp "$run/analysis/$f" "$RES/${id}_analysis_$f"
  done
  for f in fig_zinc.png fig_ligand.png fig_orientation.png fig_contacts.png; do
    [ -f "$run/analysis/$f" ] && cp "$run/analysis/$f" "$RES/${id}_analysis_$f"
  done
  snap=""
  [ "$id" = X1 ] && [ -f "$RES/X1_interrupt_snapshot.json" ] && snap="--snapshot $RES/X1_interrupt_snapshot.json"
  nice -n 10 "$PY" "$MD/test/pi_test_check.py" check "$run" $snap --out "$RES/${id}_check.json" > "$RES/${id}_check.console" 2>&1 || status=1
  tail -n 1 "$RES/${id}_check.console"
  nice -n 10 "$PY" "$MD/test/pi_test_look.py" "$run" --out-prefix "$RES/$id" || status=1
done
# the campaign path of the analysis on the two test runs (a table of two runs; no rule of the plan can be read from it)
nice -n 10 "$PY" "$MD/analysis/analyse_md.py" --campaign "$RUNS/restrained" --reuse --out "$RUNS/analysis_campaign" \
    > "$RES/campaign_analyse_md.console" 2>&1; code=$?
echo "analyse_md.py --campaign: exit code $code" | tee -a "$RES/campaign_analyse_md.console"
[ $code -eq 0 ] || status=1
for f in campaign_table.csv campaign_table.md campaign.json; do
  [ -f "$RUNS/analysis_campaign/$f" ] && cp "$RUNS/analysis_campaign/$f" "$RES/$f"
done
nice -n 10 "$PY" "$MD/test/pi_test_report.py" --runs "$RUNS" --out "$RES/pi_test_speed.json" > "$RES/pi_test_speed.console" 2>&1 || status=1
tail -n 1 "$RES/pi_test_speed.console"
exit $status
