#!/usr/bin/env bash
# Input-geometry diagnostic for the AutoDock4Zn ring orientation of sulfanilamide in hCA II
# (production 2026-09-28).
#
# What it tests.  Docking the DFT geometry instead of the ETKDG/MMFF94s conformer moved the AD4Zn pose
# of sulfanilamide: in 6RL9 every seed now lands on the near-native pose (0.61 A, 20/20, where the
# validation with ETKDG conformers was bimodal at 0.5), and in 3HS4 every seed lands on the pose
# rotated about the Zn-N/S axis (3.2 A, 20/20, where the validation with ETKDG conformers gave 0.62 A).
# The zinc anchor is unchanged in both.  AD4Zn separates the two orientations by ~0.02 kcal/mol
# (dock/VALIDATION.md section 3), and the two inputs differ most in the C4-S1-N2 angle -- MMFF94s
# 101.4-102.0 deg against 108.1-109.2 deg in four crystal structures and 109.5 deg at
# B3LYP-D3(BJ)/PCM.  Vina holds bond angles rigid, so that angle is the candidate, and this is the
# 2 x 2 test of it:
#
#              C4-S1-N2 as in the source        C4-S1-N2 swapped
#   DFT        the production run (__prod)      __csn102   (DFT geometry, angle -> 102.0)
#   ETKDG      __etkdg (fixed seed-1 conformer) __csn109   (ETKDG seed 1, angle -> 109.5)
#
# __etkdg also separates "one fixed input geometry for every seed" (production) from "a fresh ETKDG
# conformer per seed" (validation), which is the other difference between the two protocols.
#
# Geometries: dock/results/diagnostics/*.xyz, written by dock/geometry_variant.py (the DFT one) and by
# the same routine for the ETKDG conformer; each changes exactly one angle and nothing else.
#
# 6 jobs, 10 seeds, exhaustiveness 32, --cpu 3, nice -- about 35 min on the Pi.  Run it only when no
# other docking is running (one process tree at a time).
set -u
cd "$(dirname "$0")/.."
PY=${PY:-/root/cqc-sfa-tools/dock/bin/python}
export PATH=$(dirname "$PY"):$PATH
CPU=${CPU:-3}
SEEDS=${SEEDS:-10}
EXH=${EXH:-32}
LOG=dock/results/diagnostics/diagnostic.log
mkdir -p dock/results/diagnostics
ist () { TZ=Asia/Kolkata date '+%d %b %I:%M %p IST'; }

for t in ${TARGETS:-T3_6RL9 T4_3HS4}; do
  for pair in "sfa_neutral_csn102.xyz:__csn102" "sfa_neutral_etkdg_csn109.xyz:__csn109" \
              "sfa_neutral_etkdg_seed1.xyz:__etkdg"; do
    xyz=dock/results/diagnostics/${pair%%:*}
    suf=${pair##*:}
    [ -f "$xyz" ] || { echo "[$(ist)] missing $xyz" | tee -a "$LOG"; exit 1; }
    echo "[$(ist)] START $t sfa_neutral__ad4zn$suf  ($xyz)" | tee -a "$LOG"
    nice -n 10 "$PY" dock/run_ad4zn.py --target "$t" --ligand sfa_neutral --stage dock \
      --seeds "$SEEDS" --exhaustiveness "$EXH" --cpu "$CPU" --geometry-xyz "$xyz" \
      --tag-suffix "$suf" >> "$LOG" 2>&1 \
      && echo "[$(ist)] DONE  $t sfa_neutral__ad4zn$suf  | $(tail -1 "$LOG")" | tee -a "$LOG" \
      || echo "[$(ist)] FAIL  $t sfa_neutral__ad4zn$suf" | tee -a "$LOG"
  done
done
echo "[$(ist)] DIAGNOSTIC DONE" | tee -a "$LOG"
