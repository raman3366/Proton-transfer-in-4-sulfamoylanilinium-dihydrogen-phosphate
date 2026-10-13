#!/usr/bin/env bash
# Four ordered crystal waters in the 6RL9 receptor: does keeping them remove the ring-rotated pose of
# sulfanilamide from the AutoDock4Zn results?  (6 Oct 2026; report: dock/WATER4.md)
#
# The ring-rotated pose lies on waters 423, 516, 553 and 568 of 6RL9 chain A, which are in the crystal
# together with the ligand in its crystal pose (md/SYSTEM.md section 4).  The reported receptors hold no
# water.  This script docks into the variants written by dock/water_variant.py, in which only those
# four waters are added, with the settings of the reported runs, for both states of Glu106
# (dock/GLU106.md), and compares each run with the run of the same name in the dry receptor.
#
#   SCRATCH=<dir> bash dock/run_water4_sensitivity.sh <stage> [<stage> ...]
#
#   variants   the water variants of 6RL9 (hydrogen rules hbond, rev, bare, fixed; hbond for both Glu106 states)
#   maps       AutoDock4Zn maps of every variant that has none: TZ pseudo-atom by vendor/zinc_pseudo.py,
#              the STORED grid parameter file of the reported maps, the autogrid4 validated in dock/GLU106.md
#   rescore    stored modes of the dry runs scored without a search in the maps of each variant
#   key        sulfanilamide, DFT geometry, 20 seeds, AutoDock4Zn, both Glu106 states (the decisive runs)
#   val        the ETKDG validation run (10 seeds), both Glu106 states: the gates
#   rule2      the key run with the other hydrogen rules (rev, bare), Glu106 as reported
#   anion      the sulfonamidate (10 seeds), both Glu106 states
#   vina       Vina and Vinardo for sulfanilamide, validation and production, both Glu106 states
#   hs4        the same four sites in 3HS4 (variant T4_3HS4_wat4 from its own waters), maps, the two runs
#   compare    dock/water_variant.py --compare for every variant -> <variant>/comparison_water.json
#
# Settings are those of the stored results.json of the reported runs (exhaustiveness 32, seeds 1..10 or
# 1..20, the 20 A box, 20 modes within 3 kcal/mol), with --cpu 2 instead of --cpu 3: dock/GLU106.md
# section 3 shows byte-identical results for the two.  Vina runs at nice 15 (run_ad4zn.py adds its own
# nice -n 10, so its python is started at nice 5).  Before each job the 5-minute load is read; above
# $LOAD_MAX the job waits, at most $LOAD_WAIT_MAX minutes.  While $SCRATCH/PAUSE exists no job starts.
set -u
cd "$(dirname "$0")/.."
PY=${PY:-/root/cqc-sfa-tools/dock/bin/python}
export PATH=$(dirname "$PY"):$PATH
: "${SCRATCH:?set SCRATCH to a scratch directory}"
CPU=${CPU:-2}
EXH=32
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
AUTOGRID_PREFIX=${AUTOGRID_PREFIX:-/root/cqc-sfa-tools/autogrid}
LOAD_MAX=${LOAD_MAX:-8}
LOAD_WAIT_MAX=${LOAD_WAIT_MAX:-20}
WATERS="423 516 553 568"
LOG=$SCRATCH/run_water4.log
mkdir -p "$SCRATCH/logs"
ist () { TZ=Asia/Kolkata date '+%d %b %I:%M:%S %p IST'; }
say () { echo "[$(ist)] $*" | tee -a "$LOG"; }
wait_turn () {
  while [ -e "$SCRATCH/PAUSE" ]; do sleep 15; done
  local waited=0
  while awk -v m="$LOAD_MAX" '{exit !($2 > m)}' /proc/loadavg && [ "$waited" -lt "$LOAD_WAIT_MAX" ]; do
    [ "$waited" = 0 ] && say "WAIT  5-minute load $(cut -d' ' -f2 /proc/loadavg) above $LOAD_MAX"
    sleep 60; waited=$((waited + 1))
  done
  [ "$waited" -gt 0 ] && say "      waited $waited min, 5-minute load now $(cut -d' ' -f2 /proc/loadavg)"
  return 0
}

done_already () {
  "$PY" - "$1" "$2" "$EXH" <<'EOF'
import json, sys, pathlib
p = pathlib.Path(sys.argv[1]) / 'results.json'
if not p.exists():
    sys.exit(1)
r = json.loads(p.read_text())
ok = [s for s in r.get('seeds', []) if s.get('ok')]
sys.exit(0 if len(ok) == int(sys.argv[2]) and int(r.get('exhaustiveness', -1)) == int(sys.argv[3]) else 1)
EOF
}

# one job: runner target ligand scoring seeds geometry suffix
job () {
  local runner=$1 t=$2 l=$3 s=$4 n=$5 g=$6 suf=$7
  local tag="${l}__${s}${suf}" rundir="dock/$t/runs/${l}__${s}${suf}"
  [ -f "dock/$t/receptor.pdbqt" ] || { say "SKIP  $t $tag (no such variant)"; return 0; }
  if done_already "$rundir" "$n"; then say "SKIP  $t $tag (complete)"; return 0; fi
  wait_turn
  say "START $t $tag seeds=$n geometry=$g cpu=$CPU"
  local t0; t0=$(date +%s)
  if [ "$runner" = ad4zn ]; then
    nice -n 5 "$PY" dock/run_ad4zn.py --target "$t" --ligand "$l" --stage dock --seeds "$n" \
      --exhaustiveness "$EXH" --cpu "$CPU" --geometry "$g" --tag-suffix "$suf" > "$SCRATCH/logs/${t}__${tag}.log" 2>&1
  else
    nice -n 15 "$PY" dock/run_vina.py --target "$t" --ligand "$l" --scoring "$s" --seeds "$n" \
      --exhaustiveness "$EXH" --cpu "$CPU" --geometry "$g" --tag-suffix "$suf" > "$SCRATCH/logs/${t}__${tag}.log" 2>&1
  fi
  say "DONE  $t $tag rc=$? $(( $(date +%s) - t0 )) s | $(grep -E 'kcal/mol \(n=|ad4zn: ' "$SCRATCH/logs/${t}__${tag}.log" | tail -1)"
}

stage_variants () {
  for spec in "T3_6RL9 hbond" "T3_6RL9_glu106GLU hbond" "T3_6RL9 rev" "T3_6RL9 bare" "T3_6RL9 fixed"; do
    set -- $spec
    nice -n 15 "$PY" dock/water_variant.py --target "$1" --waters $WATERS --rule "$2" || exit 1
  done
}

stage_maps () {
  export PATH=$AUTOGRID_PREFIX/bin:$PATH
  local want; want=$("$PY" -c "import json; print(json.load(open('dock/T3_6RL9_glu106GLU/ad4zn/autogrid_validation.json'))['autogrid4_sha256'])")
  [ "$(sha256sum "$(command -v autogrid4)" | cut -d' ' -f1)" = "$want" ] || { say "autogrid4 is not the validated binary: stop"; exit 1; }
  for d in dock/T3_6RL9*wat* dock/T4_3HS4*wat*; do
    [ -f "$d/receptor.pdbqt" ] || continue
    local v; v=$(basename "$d"); local W=$d/ad4zn
    ls "$W"/receptor_tz.*.map > /dev/null 2>&1 && continue
    local src; src=$("$PY" -c "import json; print(json.load(open('$d/variant.json'))['source'])")
    local base=${v:0:7}
    mkdir -p "$W"
    wait_turn
    say "START maps $v (dry: $src; gpf of dock/$base)"
    "$PY" dock/vendor/zinc_pseudo.py -r "$d/receptor.pdbqt" -o "$W/receptor_tz.pdbqt" > "$SCRATCH/logs/maps__$v.log" 2>&1 || { say "zinc_pseudo.py failed"; exit 1; }
    cp "dock/$base/ad4zn/receptor_tz.gpf" "$W/receptor_tz.gpf"
    cp dock/vendor/AD4Zn.dat "$W/"
    ( cd "$W" && nice -n 15 autogrid4 -p receptor_tz.gpf -l receptor_tz.glg ) >> "$SCRATCH/logs/maps__$v.log" 2>&1
    say "DONE  maps $v rc=$? | $(tail -2 "$W/receptor_tz.glg" | tr '\n' ' ')"
    "$PY" - "$src" "$v" "$base" "$(command -v autogrid4)" <<'EOF' | tee -a "$LOG"
import hashlib, json, pathlib, sys
import numpy as np
src, v, base, exe = sys.argv[1:5]
B, W = pathlib.Path('dock') / src / 'ad4zn', pathlib.Path('dock') / v / 'ad4zn'
sha = lambda p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
atoms = lambda p: [l for l in pathlib.Path(p).read_text().splitlines() if l.startswith(('ATOM', 'HETATM'))]
xyz = lambda l: np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
ob, nv = atoms(B / 'receptor_tz.pdbqt'), atoms(W / 'receptor_tz.pdbqt')
key = lambda l: (l[12:26], l[30:54], l[66:])          # everything but the serial number
kb = {key(l) for l in ob}
new = [l for l in nv if key(l) not in kb]
gone = [l for l in ob if key(l) not in {key(x) for x in nv}]
tzb = [xyz(l) for l in ob if l[77:79].strip() == 'TZ']
tzv = [xyz(l) for l in nv if l[77:79].strip() == 'TZ']
glg = (W / 'receptor_tz.glg').read_text().splitlines()
pick = lambda k: next((l.strip() for l in glg if k in l), None)
rep = dict(target=v, dry=src, autogrid4=exe, autogrid4_sha256=sha(exe),
           gpf_identical_to_reported=sha(pathlib.Path('dock') / base / 'ad4zn' / 'receptor_tz.gpf') == sha(W / 'receptor_tz.gpf'),
           receptor_tz_atoms_added=len(new), receptor_tz_atoms_added_not_water=sum(1 for l in new if l[17:20] != 'HOH'),
           receptor_tz_atoms_lost=len(gone), tz_shift_vs_dry_A=float(np.linalg.norm(tzv[0] - tzb[0])),
           glg=dict(version=pick('AutoGrid 4'), total_charge=pick('Total charge'), completion=pick('Successful Completion'),
                    warnings=[l.strip() for l in glg if 'WARNING' in l or 'no bonded atoms' in l][:12]), maps={})
for m in sorted(W.glob('receptor_tz.*.map')):
    a = np.array([float(x) for x in (B / m.name).read_text().splitlines()[6:]])
    c = np.array([float(x) for x in m.read_text().splitlines()[6:]])
    rep['maps'][m.name] = dict(sha256=sha(m), n_values=int(len(c)), n_values_changed=int((a != c).sum()),
                               header_identical_to_dry=(B / m.name).read_text().splitlines()[:6] == m.read_text().splitlines()[:6])
(W / 'maps_report.json').write_text(json.dumps(rep, indent=2) + '\n')
ok = rep['gpf_identical_to_reported'] and rep['receptor_tz_atoms_added_not_water'] == 0 and rep['receptor_tz_atoms_lost'] == 0 and len(rep['maps']) == 9
print(f"  {v}: gpf identical {rep['gpf_identical_to_reported']}; atoms added {rep['receptor_tz_atoms_added']} (not water: "
      f"{rep['receptor_tz_atoms_added_not_water']}), lost {rep['receptor_tz_atoms_lost']}; TZ shift {rep['tz_shift_vs_dry_A']:.3f} A; "
      f"{rep['glg']['total_charge']}; maps {len(rep['maps'])}; warnings {len(rep['glg']['warnings'])}")
sys.exit(0 if ok else 1)
EOF
    [ "${PIPESTATUS[0]}" = 0 ] || { say "maps of $v failed their checks: stop"; exit 1; }
  done
}

stage_rescore () {
  mkdir -p "$SCRATCH/rescore"
  nice -n 15 "$PY" dock/water_variant.py --rescore --scratch "$SCRATCH/rescore" 2>&1 | grep -v 'tagged as 2D\|RDKit\] WARNING\|^$' | tee -a "$LOG"
}

stage_key () {
  job ad4zn T3_6RL9_wat4           sfa_neutral ad4zn 20 dft __prod
  job ad4zn T3_6RL9_glu106GLU_wat4 sfa_neutral ad4zn 20 dft __prod
  say "KEY STAGE DONE"
}
stage_val () {
  job ad4zn T3_6RL9_wat4           sfa_neutral ad4zn 10 etkdg ""
  job ad4zn T3_6RL9_glu106GLU_wat4 sfa_neutral ad4zn 10 etkdg ""
  say "VAL STAGE DONE"
}
stage_rule2 () {
  job ad4zn T3_6RL9_wat4rev  sfa_neutral ad4zn 20 dft __prod
  job ad4zn T3_6RL9_wat4bare sfa_neutral ad4zn 20 dft __prod
  say "RULE2 STAGE DONE"
}
stage_anion () {
  job ad4zn T3_6RL9_wat4           sfa_anion ad4zn 10 etkdg ""
  job ad4zn T3_6RL9_glu106GLU_wat4 sfa_anion ad4zn 10 etkdg ""
  say "ANION STAGE DONE"
}
stage_vina () {
  for t in T3_6RL9_wat4 T3_6RL9_glu106GLU_wat4; do
    job vina "$t" sfa_neutral vina    10 dft __prod
    job vina "$t" sfa_neutral vinardo 10 dft __prod
    job vina "$t" sfa_neutral vina    10 etkdg ""
    job vina "$t" sfa_neutral vinardo 10 etkdg ""
  done
  say "VINA STAGE DONE"
}
stage_hs4 () {
  [ -f dock/T4_3HS4_wat4/receptor.pdbqt ] || { say "T4_3HS4_wat4 has not been made (see dock/WATER4.md): nothing to run"; return 0; }
  stage_maps
  job ad4zn T4_3HS4_wat4 sfa_neutral ad4zn 20 dft __prod
  job ad4zn T4_3HS4_wat4 sfa_neutral ad4zn 10 etkdg ""
  say "3HS4 STAGE DONE"
}
stage_compare () {
  nice -n 15 "$PY" dock/water_variant.py --compare-all --scratch "$SCRATCH/compare_water" ${COMPARE_ARGS:-} 2>&1 \
    | grep -v 'tagged as 2D\|RDKit\] WARNING\|^$' | tee -a "$LOG"
}

for stage in "$@"; do
  case "$stage" in
    variants) stage_variants;; maps) stage_maps;; rescore) stage_rescore;; key) stage_key;; val) stage_val;;
    rule2) stage_rule2;; anion) stage_anion;; vina) stage_vina;; hs4) stage_hs4;; compare) stage_compare;;
    *) echo "unknown stage $stage"; exit 2;;
  esac
done
