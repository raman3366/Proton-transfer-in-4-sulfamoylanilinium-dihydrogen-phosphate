#!/usr/bin/env bash
# Production docking of the DFT ligand geometries.  Run this only after gate G2 has been reviewed:
# redocking validates the protocol, and a protocol that fails G2 for a target must not be used to
# make claims about that target (dock/VALIDATION.md sections 10-11 say which target supports what).
#
# Usage
#   bash dock/run_production.sh                    # Pi: sequential, --cpu 3, nice 10 (~5-6 h)
#   RUN_ANALYSIS=0 bash dock/run_production.sh     # docking only
#   PARALLEL=12 CPU=4 bash dock/run_production.sh  # 48-vCPU VM (not used: the VM is busy with DFT)
#
# Changes for the production run (2026-09-28), each for a reason found while preparing it:
#   1. The DFT geometries were never docked.  prep_ligands.py --dft wrote them to dock/ligands/, but
#      run_vina.py and run_ad4zn.py rebuilt an ETKDG conformer from the SMILES for every seed and never
#      read that file.  Both runners now take --geometry dft, and the three sulfanilamide species
#      are docked with it (the controls keep ETKDG: they were never optimised).
#   2. The production run directories were the validation run directories.  A production job
#      tagged <ligand>__<scoring> would have overwritten, e.g., the 6RL9 sulfanilamide redock that
#      gate G2 rests on.  Every production job now carries --tag-suffix __prod.
#   3. AD4Zn only for the species with a neutral sulfonamide (sulfanilamide, sulfamethoxazole,
#      acetazolamide): its TZ pseudo-atom models the Zn-N(H) bond of the neutral form, and the anion
#      gave Zn-N 4.86 A for acetazolamide in validation.  20 seeds for sulfanilamide (bimodal in 6RL9).
#   4. Each target's co-crystal ligand is docked as its reference score.  For 3TZF (08D =
#      sulfamethoxazole), 6RL9 (SAN = sulfanilamide) and 3HS4 (AZM = acetazolamide) that ligand is
#      already a species of the matrix; 7S2K's dihydropteroate (78H) is added for Vina and Vinardo.
#   5. Jobs run one at a time in priority order (primary ligand first), resume after an interruption
#      (a job whose results.json already holds every seed with the same settings is skipped), and
#      write a progress log with an ETA (IST) to dock/results/production.log.
#   T2 1AJ0 stays out of the matrix: it failed validation (11.6 A decoy cleft in every seed).
#
# Reproducibility (dock/VALIDATION.md section 7): at fixed --seed a rerun is bit-identical, and
# --cpu does not change the result on this build, only the wall time.
set -u
cd "$(dirname "$0")/.."
PY=${PY:-/root/cqc-sfa-tools/dock/bin/python}
export PATH=$(dirname "$PY"):$PATH
PARALLEL=${PARALLEL:-1}
CPU=${CPU:-3}
SEEDS=${SEEDS:-10}
POSE_SEEDS=${POSE_SEEDS:-20}
EXH=${EXH:-32}
SUFFIX=${SUFFIX:-__prod}
RUN_ANALYSIS=${RUN_ANALYSIS:-1}
TARGETS=${TARGETS:-"T5_7S2K T4_3HS4 T3_6RL9 T1_3TZF"}
RES=dock/results
LOG=$RES/production.log
mkdir -p "$RES/logs"
echo $$ > "$RES/production.pid"
ist () { TZ=Asia/Kolkata date '+%d %b %I:%M %p IST'; }
say () { echo "[$(ist)] $*" | tee -a "$LOG"; }

# ---------------------------------------------------------------------------------------------
# 1. Ligands from the DFT geometries (B3LYP-D3(BJ)/6-311++G(d,p), IEF-PCM water).
# ---------------------------------------------------------------------------------------------
DFT_ARGS=()
add_dft () { [ -f "$2" ] && DFT_ARGS+=("$1=$2") || { say "ABORT: $2 missing"; exit 1; }; }
add_dft sfa_neutral molecules/sfa/results/sulfanilamide_pcm_d3/opt.xyz
add_dft sfa_anion   molecules/sfa/results/sulfonamidate_pcm_d3/opt.xyz
add_dft sfa_cation  molecules/sfa/results/cation_pcm_d3/opt.xyz
say "=== ligand preparation (DFT geometries for sfa_neutral, sfa_anion, sfa_cation)"
nice -n 10 "$PY" dock/prep_ligands.py --dft "${DFT_ARGS[@]}" 2>&1 | grep -v 'tagged as 2D\|matching pattern\|^\[RDKit\] WARNING:$\|^$' | tee -a "$LOG"

# ---------------------------------------------------------------------------------------------
# 2. The production matrix, in priority order.  Each line: runner target ligand scoring seeds.
#    Expected seconds per seed on the Pi at --cpu 3 (validation runs, dock/T*/runs/*/results.json)
#    are used for the ETA only.
# ---------------------------------------------------------------------------------------------
cost () {  # scoring ligand -> expected s/seed
  case "$1:$2" in
    vina:sfa_*) echo 14;;  vina:paba_*) echo 12;;  vina:azm_*) echo 16;;  vina:smx_*) echo 30;;
    vina:dhp_*) echo 57;;
    vinardo:sfa_*) echo 22;; vinardo:paba_*) echo 18;; vinardo:azm_*) echo 22;;
    vinardo:smx_*) echo 31;; vinardo:dhp_*) echo 60;;
    ad4zn:sfa_*) echo 34;;  ad4zn:azm_*) echo 36;;  ad4zn:*) echo 60;;
    *) echo 30;;
  esac
}
JOBS=()
vv () { for t in $TARGETS; do for s in vina vinardo; do JOBS+=("vina $t $1 $s $SEEDS"); done; done; }
zn () { for t in T4_3HS4 T3_6RL9; do case " $TARGETS " in *" $t "*) JOBS+=("ad4zn $t $1 ad4zn $2");; esac; done; }
# A. the compound's primary species
vv sfa_neutral;  zn sfa_neutral "$POSE_SEEDS"
# B. reference ligands and the key controls
case " $TARGETS " in *" T5_7S2K "*) JOBS+=("vina T5_7S2K dhp_anion vina $SEEDS" "vina T5_7S2K dhp_anion vinardo $SEEDS");; esac
vv smx_anion;  vv azm_neutral;  zn azm_neutral "$SEEDS";  vv sfa_anion
# C. the remaining controls and the cation sensitivity case
vv smx_neutral;  zn smx_neutral "$SEEDS";  vv azm_anion;  vv paba_anion;  vv sfa_cation

geom_of () { case "$1" in sfa_*) echo dft;; *) echo etkdg;; esac; }
total_cost=0
for j in "${JOBS[@]}"; do set -- $j; total_cost=$((total_cost + $(cost "$4" "$3") * $5)); done
say "=== ${#JOBS[@]} jobs, PARALLEL=$PARALLEL CPU=$CPU EXH=$EXH SEEDS=$SEEDS (AD4Zn sulfanilamide $POSE_SEEDS) suffix=$SUFFIX"
say "    expected docking time $((total_cost / 3600)) h $(((total_cost % 3600) / 60)) min (validation rates)"

if [ "${DRY_RUN:-0}" = 1 ]; then printf '%s\n' "${JOBS[@]}"; rm -f "$RES/production.pid"; exit 0; fi

done_already () {  # rundir seeds geometry -> 0 if the job is complete with these settings
  "$PY" - "$1" "$2" "$3" "$EXH" <<'EOF'
import json, sys, pathlib
p = pathlib.Path(sys.argv[1]) / 'results.json'
if not p.exists():
    sys.exit(1)
r = json.loads(p.read_text())
ok = [s for s in r.get('seeds', []) if s.get('ok')]
sys.exit(0 if (len(ok) == int(sys.argv[2]) and r.get('geometry', 'etkdg') == sys.argv[3]
               and int(r.get('exhaustiveness', -1)) == int(sys.argv[4])) else 1)
EOF
}

run_one () {  # runner target ligand scoring seeds
  local runner=$1 t=$2 l=$3 s=$4 n=$5 g
  g=$(geom_of "$l")
  if [ "$runner" = ad4zn ]; then
    nice -n 10 "$PY" dock/run_ad4zn.py --target "$t" --ligand "$l" --stage dock --seeds "$n" \
      --exhaustiveness "$EXH" --cpu "$CPU" --geometry "$g" --tag-suffix "$SUFFIX"
  else
    nice -n 10 "$PY" dock/run_vina.py --target "$t" --ligand "$l" --scoring "$s" --seeds "$n" \
      --exhaustiveness "$EXH" --cpu "$CPU" --geometry "$g" --tag-suffix "$SUFFIX"
  fi
}

if [ "$PARALLEL" -gt 1 ]; then
  # VM path: no progress accounting, one line per job
  export -f run_one geom_of; export PY EXH CPU SUFFIX
  printf '%s\n' "${JOBS[@]}" | xargs -P "$PARALLEL" -I{} bash -c 'run_one {} > /dev/null 2>&1 || echo "FAILED: {}"'
else
  i=0; done_cost=0; t_start=$(date +%s); spent_cost=0
  for j in "${JOBS[@]}"; do
    i=$((i + 1)); set -- $j
    runner=$1; t=$2; l=$3; s=$4; n=$5; c=$(( $(cost "$s" "$l") * n ))
    rundir="dock/$t/runs/${l}__${s}${SUFFIX}"
    if done_already "$rundir" "$n" "$(geom_of "$l")"; then
      say "SKIP  $i/${#JOBS[@]} $t ${l}__${s}${SUFFIX} (complete)"; done_cost=$((done_cost + c)); continue
    fi
    say "START $i/${#JOBS[@]} $t ${l}__${s}${SUFFIX} seeds=$n geometry=$(geom_of "$l")"
    t0=$(date +%s)
    run_one "$runner" "$t" "$l" "$s" "$n" > "$RES/logs/${t}__${l}__${s}${SUFFIX}.log" 2>&1
    rc=$?
    dt=$(( $(date +%s) - t0 ))
    done_cost=$((done_cost + c)); spent_cost=$((spent_cost + c))
    elapsed=$(( $(date +%s) - t_start ))
    # ETA from the measured speed of this session relative to the validation rates
    rem=$(( (total_cost - done_cost) * elapsed / (spent_cost > 0 ? spent_cost : 1) ))
    pct=$(( 100 * done_cost / total_cost ))
    tail1=$(grep -E 'kcal/mol \(n=|ad4zn: ' "$RES/logs/${t}__${l}__${s}${SUFFIX}.log" | tail -1)
    if [ $rc -ne 0 ]; then
      say "FAIL  $i/${#JOBS[@]} $t ${l}__${s}${SUFFIX} rc=$rc after ${dt}s (see $RES/logs/)"
    else
      say "DONE  $i/${#JOBS[@]} $t ${l}__${s}${SUFFIX} ${dt}s  | ${tail1##*: }"
    fi
    say "      progress ${pct}% (by expected cost)  elapsed $((elapsed / 60)) min  ETA $(TZ=Asia/Kolkata date -d "@$(( $(date +%s) + rem ))" '+%I:%M %p IST') (~$((rem / 60)) min left)"
  done
fi
say "=== docking finished"

# ---------------------------------------------------------------------------------------------
# 3. Analysis: analyse.py --production (scores, clusters, RMSD to the co-crystal ligand and to its
#    shared core, ligand efficiency, rank order, torsions and MMFF strain, poses) and
#    contacts.py --production (PLIP + ProLIF + geometric contacts), then the report.
# ---------------------------------------------------------------------------------------------
if [ "$RUN_ANALYSIS" = 1 ]; then
  say "=== analysis"
  nice -n 10 "$PY" dock/analyse.py --production 2>&1 | grep -v 'tagged as 2D' | tee -a "$LOG"
  nice -n 10 "$PY" dock/contacts.py --production 2>&1 | grep -v 'tagged as 2D' | tee -a "$LOG"
  nice -n 10 "$PY" dock/production_report.py 2>&1 | tee -a "$LOG"
fi
say "PRODUCTION DONE"
rm -f "$RES/production.pid"
