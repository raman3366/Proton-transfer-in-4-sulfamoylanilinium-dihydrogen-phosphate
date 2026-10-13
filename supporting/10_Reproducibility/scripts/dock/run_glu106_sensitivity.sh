#!/usr/bin/env bash
# Glu106 protonation sensitivity test for the two hCA II receptors (5 Oct 2026); report: dock/GLU106.md.
#
# Both prepared hCA II receptors carry Glu106 as the neutral acid (PROPKA, run without the zinc).  This
# script docks into the two variants written by dock/glu_protonation_variant.py, in which only Glu106
# differs (carboxylate), with the settings of the reported runs, and compares each run with the stored
# run of the same name in the unmodified receptor.
#
#   SCRATCH=<dir> bash dock/run_glu106_sensitivity.sh <stage> [<stage> ...]
#
#   variants        write dock/T3_6RL9_glu106GLU and dock/T4_3HS4_glu106GLU
#   baseline        rerun two stored runs of the UNMODIFIED 6RL9 receptor (sulfanilamide, Vina and
#                   Vinardo) in a copy of the drivers under $SCRATCH and compare inputs and poses with
#                   the stored runs byte for byte -> dock/T3_6RL9_glu106GLU/baseline_check.json
#   vina            the Vina and Vinardo runs of sulfanilamide (neutral, anion) and acetazolamide
#   autogrid-check  regenerate the AutoDock4Zn maps of the UNMODIFIED receptors with the stored grid
#                   parameter files and compare them with the stored maps, number by number
#                   -> dock/<variant>/ad4zn/autogrid_validation.json
#   maps            maps of the two variants: TZ pseudo-atom by vendor/zinc_pseudo.py, the STORED grid
#                   parameter file, local autogrid4; refuses unless autogrid-check passed
#   ad4zn-baseline  rerun the stored AutoDock4Zn run of sulfanilamide in the UNMODIFIED 6RL9 receptor
#                   with the stored maps -> dock/T3_6RL9_glu106GLU/baseline_check_ad4zn.json
#   ad4zn           the AutoDock4Zn runs of sulfanilamide and acetazolamide
#   ad4zn-diag      the six input-geometry diagnostics of sulfanilamide (csn102, csn109, etkdg)
#   rescore         stored AutoDock4Zn modes of sulfanilamide scored, without a search, in the reported
#                   maps and in the maps of the variant -> dock/<variant>/rescore_fixed_poses.json
#   rescore-vina    the same for stored Vina and Vinardo top poses, in the reported receptor and in the
#                   variant -> dock/<variant>/rescore_fixed_poses_vina.json
#   compare         dock/glu_protonation_variant.py --compare for both variants
#
# Settings.  Everything is read from the stored results.json of the reported runs: exhaustiveness 32,
# seeds 1..10 (1..20 for the two production AutoDock4Zn jobs of sulfanilamide), the 20 A box of
# box.json, 20 modes within 3 kcal/mol, ETKDG conformers for the validation runs and the controls, the
# DFT geometries for the production runs of sulfanilamide.  One setting differs: --cpu 2 instead of
# --cpu 3, because the machine is shared.  --cpu changes the wall time and not the result on this
# build (dock/VALIDATION.md section 7); the `baseline` stage checks that again, for Vina and for
# Vinardo, before any variant is docked.
#
# Run names are those of the stored runs, so that every variant run has its baseline next door.  A
# production job of an ETKDG control repeats its validation run bit for bit in the stored data, so the
# controls are docked once and stand for both.
#
# The stage `maps` needs autogrid4 (AutoGrid 4.2.7.x) in $AUTOGRID_PREFIX/bin; how that binary was
# built is recorded in dock/GLU106.md.  Jobs run one at a time under nice; while the file
# $SCRATCH/PAUSE exists no new job starts.
set -u
cd "$(dirname "$0")/.."
PY=${PY:-/root/cqc-sfa-tools/dock/bin/python}
export PATH=$(dirname "$PY"):$PATH
: "${SCRATCH:?set SCRATCH to a scratch directory}"
CPU=${CPU:-2}
EXH=32
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}       # autogrid4 only; Vina takes --cpu
AUTOGRID_PREFIX=${AUTOGRID_PREFIX:-/root/cqc-sfa-tools/autogrid}
VARIANTS="T3_6RL9_glu106GLU T4_3HS4_glu106GLU"
LOG=$SCRATCH/run.log
mkdir -p "$SCRATCH/logs"
ist () { TZ=Asia/Kolkata date '+%d %b %I:%M:%S %p IST'; }
say () { echo "[$(ist)] $*" | tee -a "$LOG"; }
wait_pause () { while [ -e "$SCRATCH/PAUSE" ]; do sleep 15; done; }

done_already () {  # rundir seeds -> 0 if results.json holds that many good seeds at EXH
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

# one job: runner target ligand scoring seeds geometry suffix [xyz]
job () {
  local runner=$1 t=$2 l=$3 s=$4 n=$5 g=$6 suf=$7 xyz=${8:-}
  local tag="${l}__${s}${suf}" rundir="dock/$t/runs/${l}__${s}${suf}" geo=(--geometry "$g")
  [ -n "$xyz" ] && geo=(--geometry-xyz "$xyz")
  wait_pause
  if done_already "$rundir" "$n"; then say "SKIP  $t $tag (complete)"; return 0; fi
  say "START $t $tag seeds=$n geometry=${xyz:-$g} cpu=$CPU"
  local t0; t0=$(date +%s)
  if [ "$runner" = ad4zn ]; then
    # no nice here: run_ad4zn.py puts its own `nice -n 10` on vina, and a second one makes it nice 19,
    # which other nice-10 jobs starve (seen on 5 Oct 2026; the first AutoDock4Zn jobs ran that way)
    "$PY" dock/run_ad4zn.py --target "$t" --ligand "$l" --stage dock --seeds "$n" \
      --exhaustiveness "$EXH" --cpu "$CPU" "${geo[@]}" --tag-suffix "$suf" \
      > "$SCRATCH/logs/${t}__${tag}.log" 2>&1
  else
    nice -n 10 "$PY" dock/run_vina.py --target "$t" --ligand "$l" --scoring "$s" --seeds "$n" \
      --exhaustiveness "$EXH" --cpu "$CPU" "${geo[@]}" --tag-suffix "$suf" \
      > "$SCRATCH/logs/${t}__${tag}.log" 2>&1
  fi
  local rc=$?
  say "DONE  $t $tag rc=$rc $(( $(date +%s) - t0 )) s | $(grep -E 'kcal/mol \(n=|ad4zn: ' "$SCRATCH/logs/${t}__${tag}.log" | tail -1)"
}

stage_variants () {
  for t in T3_6RL9 T4_3HS4; do "$PY" dock/glu_protonation_variant.py --target "$t" --resnum 106 || exit 1; done
}

# Rerun of stored runs in the unmodified receptor.  The drivers resolve every path from their own
# location, so a copy of them under $SCRATCH, with the unmodified receptor and box linked in, docks
# the same receptor with the same code and writes nothing into dock/.
stage_baseline () {
  local M=$SCRATCH/mirror/dock
  mkdir -p "$M/T3_6RL9"
  cp dock/run_vina.py dock/run_ad4zn.py dock/prep_ligands.py dock/targets.py "$M/"
  ln -sfn "$PWD/dock/ccd" "$M/ccd"; ln -sfn "$PWD/dock/ligands" "$M/ligands"
  ln -sfn "$PWD/dock/T3_6RL9/receptor.pdbqt" "$M/T3_6RL9/receptor.pdbqt"
  ln -sfn "$PWD/dock/T3_6RL9/box.json" "$M/T3_6RL9/box.json"
  for s in vina vinardo; do
    wait_pause
    say "START baseline rerun T3_6RL9 sfa_neutral__$s (unmodified receptor, cpu=$CPU, stored run: cpu 3)"
    nice -n 10 "$PY" "$M/run_vina.py" --target T3_6RL9 --ligand sfa_neutral --scoring "$s" --seeds 10 \
      --exhaustiveness "$EXH" --cpu "$CPU" > "$SCRATCH/logs/baseline__sfa_neutral__$s.log" 2>&1
    say "DONE  baseline rerun sfa_neutral__$s | $(tail -1 "$SCRATCH/logs/baseline__sfa_neutral__$s.log")"
  done
  "$PY" - "$M" <<'EOF' | tee -a "$LOG"
import hashlib, json, pathlib, sys
M = pathlib.Path(sys.argv[1]); D = pathlib.Path('dock')
out = dict(purpose='rerun of stored runs in the UNMODIFIED receptor T3_6RL9, so that stored baseline '
                   'and new variant runs are comparable', runs={})
for tag in ('sfa_neutral__vina', 'sfa_neutral__vinardo'):
    new = json.loads((M / 'T3_6RL9/runs' / tag / 'results.json').read_text())
    old = json.loads((D / 'T3_6RL9/runs' / tag / 'results.json').read_text())
    rows = []
    for sn, so in zip(new['seeds'], old['seeds']):
        a, b = M / 'T3_6RL9/runs' / tag / f"seed{sn['seed']:02d}", D / 'T3_6RL9/runs' / tag / f"seed{so['seed']:02d}"
        rows.append(dict(seed=sn['seed'],
                         ligand_pdbqt_identical=(a / 'ligand.pdbqt').read_bytes() == (b / 'ligand.pdbqt').read_bytes(),
                         out_pdbqt_identical=(a / 'out.pdbqt').read_bytes() == (b / 'out.pdbqt').read_bytes(),
                         scores_identical=sn['scores'] == so['scores']))
    out['runs'][tag] = dict(
        stored=dict(cpu=old['cpu'], exhaustiveness=old['exhaustiveness'], vina_version=old['vina_version'],
                    best_mean=old['summary']['best_mean'], best_sd=old['summary']['best_sd'],
                    wall_s_total=old['summary']['wall_s_total']),
        rerun=dict(cpu=new['cpu'], exhaustiveness=new['exhaustiveness'], vina_version=new['vina_version'],
                   best_mean=new['summary']['best_mean'], best_sd=new['summary']['best_sd'],
                   wall_s_total=new['summary']['wall_s_total']),
        receptor_sha256=hashlib.sha256((D / 'T3_6RL9/receptor.pdbqt').read_bytes()).hexdigest(),
        n_seeds=len(rows),
        all_ligand_inputs_identical=all(r['ligand_pdbqt_identical'] for r in rows),
        all_out_pdbqt_identical=all(r['out_pdbqt_identical'] for r in rows),
        all_pose_scores_identical=all(r['scores_identical'] for r in rows), per_seed=rows)
    r = out['runs'][tag]
    print(f"baseline {tag}: inputs identical {r['all_ligand_inputs_identical']}, out.pdbqt identical "
          f"{r['all_out_pdbqt_identical']}, scores identical {r['all_pose_scores_identical']} "
          f"(stored cpu {old['cpu']}, rerun cpu {new['cpu']})")
out['reproduced'] = all(r['all_out_pdbqt_identical'] and r['all_ligand_inputs_identical'] for r in out['runs'].values())
(D / 'T3_6RL9_glu106GLU' / 'baseline_check.json').write_text(json.dumps(out, indent=2) + '\n')
EOF
}

stage_vina () {
  # sulfanilamide first, in both receptors; then acetazolamide (crystal ligand of 3HS4)
  for t in $VARIANTS; do
    job vina "$t" sfa_neutral vina 10 etkdg ""
    job vina "$t" sfa_anion   vina 10 etkdg ""
    job vina "$t" sfa_neutral vina    10 dft __prod
    job vina "$t" sfa_neutral vinardo 10 dft __prod
    job vina "$t" sfa_anion   vina    10 dft __prod
    job vina "$t" sfa_anion   vinardo 10 dft __prod
  done
  job vina T3_6RL9_glu106GLU sfa_neutral vinardo 10 etkdg ""
  for t in $VARIANTS; do
    job vina "$t" azm_neutral vina 10 etkdg ""
    job vina "$t" azm_anion   vina 10 etkdg ""
    job vina "$t" azm_anion   vinardo 10 etkdg __prod
  done
  job vina T4_3HS4_glu106GLU azm_neutral vinardo 10 etkdg ""
  job vina T3_6RL9_glu106GLU azm_neutral vinardo 10 etkdg __prod
  say "VINA STAGE DONE"
}

# AutoGrid must reproduce the stored maps of the unmodified receptors before it may make new ones.
stage_autogrid_check () {
  export PATH=$AUTOGRID_PREFIX/bin:$PATH
  command -v autogrid4 > /dev/null || { say "autogrid4 not found in $AUTOGRID_PREFIX/bin"; exit 1; }
  for t in T3_6RL9 T4_3HS4; do
    local W=$SCRATCH/autogrid_check/$t
    rm -rf "$W"; mkdir -p "$W"
    cp "dock/$t/ad4zn/receptor_tz.gpf" "dock/$t/ad4zn/receptor_tz.pdbqt" dock/vendor/AD4Zn.dat "$W/"
    wait_pause
    say "START autogrid-check $t ($(command -v autogrid4), OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset})"
    ( cd "$W" && nice -n 10 autogrid4 -p receptor_tz.gpf -l receptor_tz.glg ) > "$SCRATCH/logs/autogrid_check__$t.log" 2>&1
    say "DONE  autogrid-check $t rc=$? | $(tail -2 "$W/receptor_tz.glg" | tr '\n' ' ')"
    "$PY" - "$t" "$W" "$(command -v autogrid4)" <<'EOF' | tee -a "$LOG"
import hashlib, json, pathlib, sys
import numpy as np
t, W, exe = sys.argv[1], pathlib.Path(sys.argv[2]), sys.argv[3]
S = pathlib.Path('dock') / t / 'ad4zn'
def read(p):
    L = p.read_text().splitlines()
    return L[:6], np.array([float(x) for x in L[6:]])
def glg_facts(p):
    L = p.read_text().splitlines()
    pick = lambda key: next((l.strip() for l in L if key in l), None)
    return dict(version=pick('AutoGrid 4'), revision=pick('$Revision'), openmp=pick('_OPENMP'),
                threads=pick('CPU threads will be used'), host=pick('using:'),
                total_charge=pick('Total charge'), completion=pick('Successful Completion'))
rep = dict(target=t, autogrid4=exe, autogrid4_sha256=hashlib.sha256(pathlib.Path(exe).read_bytes()).hexdigest(),
           gpf_sha256=hashlib.sha256((S / 'receptor_tz.gpf').read_bytes()).hexdigest(),
           receptor_tz_sha256=hashlib.sha256((S / 'receptor_tz.pdbqt').read_bytes()).hexdigest(),
           parameter_file_sha256=hashlib.sha256(pathlib.Path('dock/vendor/AD4Zn.dat').read_bytes()).hexdigest(),
           stored_glg=glg_facts(S / 'receptor_tz.glg'), new_glg=glg_facts(W / 'receptor_tz.glg'), maps={})
ok = True
for sm in sorted(S.glob('receptor_tz.*.map')):
    nm = W / sm.name
    if not nm.exists():
        rep['maps'][sm.name] = dict(error='not produced'); ok = False; continue
    (h0, a), (h1, b) = read(sm), read(nm)
    same_shape = len(a) == len(b)
    d = np.abs(a - b) if same_shape else None
    # the maps are written with three decimals; the relative figure matters for the clamped walls
    rel = (d / np.maximum(1.0, np.abs(a))) if same_shape else None
    rep['maps'][sm.name] = dict(
        header_identical=h0 == h1, n_values=int(len(a)), byte_identical=sm.read_bytes() == nm.read_bytes(),
        max_abs_diff=float(d.max()) if same_shape else None,
        n_values_differing=int((d > 0).sum()) if same_shape else None,
        n_differing_by_more_than_0p001=int((d > 0.0011).sum()) if same_shape else None,
        max_rel_diff=float(rel.max()) if same_shape else None,
        max_abs_diff_where_stored_below_10=float(d[np.abs(a) < 10].max()) if same_shape else None,
        stored_min=float(a.min()), new_min=float(b.min()) if same_shape else None)
    m = rep['maps'][sm.name]
    ok = ok and same_shape and h0 == h1 and m['max_abs_diff_where_stored_below_10'] <= 0.0011 and m['max_rel_diff'] <= 1e-5
    print(f"  {t} {sm.name:22s} byte-identical {m['byte_identical']}  max|d| {m['max_abs_diff']}  "
          f"max|d| (|E|<10) {m['max_abs_diff_where_stored_below_10']}  differing values {m['n_values_differing']}/{m['n_values']}")
rep['criterion'] = ('headers identical; every value within 0.001 (the last printed decimal) where the stored '
                    'energy is below 10 kcal/mol, and within a relative 1e-5 everywhere')
rep['agrees_to_file_precision'] = bool(ok)
dst = pathlib.Path('dock') / f'{t}_glu106GLU' / 'ad4zn'
dst.mkdir(parents=True, exist_ok=True)
(dst / 'autogrid_validation.json').write_text(json.dumps(rep, indent=2) + '\n')
print(f"autogrid-check {t}: agrees to file precision = {ok}")
EOF
  done
}

stage_maps () {
  # The two steps of run_ad4zn.py --stage maps, done by hand for one reason: that stage writes a new
  # grid parameter file from the ligand types of every species in dock/ligands/, and since the
  # reported maps were made the control-receptor ligands were added there (ciprofloxacin brings type
  # F), so today it would write a different file with one more map.  Here the TZ pseudo-atom is placed
  # by the same vendor script (zinc_pseudo.py, as run_ad4zn.add_tz does) and the STORED grid parameter
  # file of the reported maps is used unchanged.
  export PATH=$AUTOGRID_PREFIX/bin:$PATH
  for v in $VARIANTS; do
    local b=${v%_glu106GLU} W=dock/$v/ad4zn
    "$PY" -c "import json,sys; sys.exit(0 if json.load(open('$W/autogrid_validation.json'))['agrees_to_file_precision'] else 1)" \
      || { say "REFUSED: autogrid-check has not passed for $v; no maps are made"; exit 1; }
    wait_pause
    say "START maps $v ($(command -v autogrid4), OMP_NUM_THREADS=$OMP_NUM_THREADS)"
    "$PY" dock/vendor/zinc_pseudo.py -r "dock/$v/receptor.pdbqt" -o "$W/receptor_tz.pdbqt" > "$SCRATCH/logs/maps__$v.log" 2>&1 || { say "zinc_pseudo.py failed"; exit 1; }
    cp "dock/$b/ad4zn/receptor_tz.gpf" "$W/receptor_tz.gpf"
    cp dock/vendor/AD4Zn.dat "$W/"                  # autogrid reads the parameter file from its cwd
    ( cd "$W" && nice -n 10 autogrid4 -p receptor_tz.gpf -l receptor_tz.glg ) >> "$SCRATCH/logs/maps__$v.log" 2>&1
    say "DONE  maps $v rc=$? | $(tail -2 "$W/receptor_tz.glg" | tr '\n' ' ')"
    "$PY" - "$b" "$v" "$(command -v autogrid4)" <<'EOF' | tee -a "$LOG"
import hashlib, json, pathlib, sys
import numpy as np
b, v, exe = sys.argv[1:4]
B, W = pathlib.Path('dock') / b / 'ad4zn', pathlib.Path('dock') / v / 'ad4zn'
sha = lambda p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
atoms = lambda p: [l for l in pathlib.Path(p).read_text().splitlines() if l.startswith(('ATOM', 'HETATM'))]
xyz = lambda l: np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
ob, nv = atoms(B / 'receptor_tz.pdbqt'), atoms(W / 'receptor_tz.pdbqt')
diff = [l for l in ob if l not in set(nv)] + [l for l in nv if l not in set(ob)]
tzb = [xyz(l) for l in ob if l[77:79].strip() == 'TZ']
tzv = [xyz(l) for l in nv if l[77:79].strip() == 'TZ']
glg = (W / 'receptor_tz.glg').read_text().splitlines()
pick = lambda key: next((l.strip() for l in glg if key in l), None)
rep = dict(target=v, made_from=b, autogrid4=exe, autogrid4_sha256=sha(exe), omp_threads=pick('CPU threads will be used'),
           gpf='copied from dock/%s/ad4zn/receptor_tz.gpf' % b,
           gpf_identical_to_reported=sha(B / 'receptor_tz.gpf') == sha(W / 'receptor_tz.gpf'),
           gpf_sha256=sha(W / 'receptor_tz.gpf'), parameter_file_sha256=sha(W / 'AD4Zn.dat'),
           receptor_tz_sha256=dict(reported=sha(B / 'receptor_tz.pdbqt'), variant=sha(W / 'receptor_tz.pdbqt')),
           receptor_tz_lines_differing=len(diff),
           receptor_tz_lines_differing_outside_Glu106=sum(
               1 for l in diff if not (l[17:20] == 'GLU' and int(l[22:26]) == 106)),
           n_tz_atoms=len(tzv),
           tz_shift_vs_reported_A=float(np.linalg.norm(tzv[0] - tzb[0])) if tzv and tzb else None,
           glg=dict(version=pick('AutoGrid 4'), revision=pick('$Revision'), total_charge=pick('Total charge'),
                    completion=pick('Successful Completion')),
           maps={})
for m in sorted(W.glob('receptor_tz.*.map')):
    a = np.array([float(x) for x in (B / m.name).read_text().splitlines()[6:]])
    c = np.array([float(x) for x in m.read_text().splitlines()[6:]])
    d = c - a
    lo = np.abs(a) < 10
    rep['maps'][m.name] = dict(
        sha256=sha(m),
        header_identical_to_reported=(B / m.name).read_text().splitlines()[:6] == m.read_text().splitlines()[:6],
        n_values=int(len(c)), n_values_changed=int((d != 0).sum()),
        max_abs_change=float(np.abs(d).max()),
        max_abs_change_where_reported_below_10=float(np.abs(d[lo]).max()),
        mean_change_where_reported_below_10=float(d[lo].mean()),
        reported_min=float(a.min()), variant_min=float(c.min()))
    r = rep['maps'][m.name]
    print(f"  {v} {m.name:20s} values changed {r['n_values_changed']:6d}/{r['n_values']}  max|change| (|E|<10) "
          f"{r['max_abs_change_where_reported_below_10']:.3f}  mean change {r['mean_change_where_reported_below_10']:+.4f}  "
          f"min {r['reported_min']} -> {r['variant_min']}")
rep['maps_ready'] = len(rep['maps']) == len(list(B.glob('receptor_tz.*.map')))
(W / 'maps_report.json').write_text(json.dumps(rep, indent=2) + '\n')
print(f"  {v}: gpf identical to the reported one {rep['gpf_identical_to_reported']}; receptor_tz lines differing "
      f"{rep['receptor_tz_lines_differing']} ({rep['receptor_tz_lines_differing_outside_Glu106']} outside Glu106); "
      f"TZ shift {rep['tz_shift_vs_reported_A']} A; {rep['glg']['total_charge']}")
sys.exit(0 if rep['maps_ready'] and rep['gpf_identical_to_reported']
         and rep['receptor_tz_lines_differing_outside_Glu106'] == 0 else 1)
EOF
    [ "${PIPESTATUS[0]}" = 0 ] || { say "maps of $v failed their checks: stop"; exit 1; }
  done
}

# The same check for AutoDock4Zn: the stored bimodal run of sulfanilamide in the UNMODIFIED 6RL9
# receptor, with the stored maps, at --cpu 2 (stored: --cpu 3).
stage_ad4zn_baseline () {
  local M=$SCRATCH/mirror/dock
  mkdir -p "$M/T3_6RL9"
  cp dock/run_vina.py dock/run_ad4zn.py dock/prep_ligands.py dock/targets.py "$M/"
  ln -sfn "$PWD/dock/ccd" "$M/ccd"; ln -sfn "$PWD/dock/ligands" "$M/ligands"
  ln -sfn "$PWD/dock/T3_6RL9/receptor.pdbqt" "$M/T3_6RL9/receptor.pdbqt"
  ln -sfn "$PWD/dock/T3_6RL9/box.json" "$M/T3_6RL9/box.json"
  ln -sfn "$PWD/dock/T3_6RL9/ad4zn" "$M/T3_6RL9/ad4zn"          # the stored maps, read only
  wait_pause
  say "START baseline rerun T3_6RL9 sfa_neutral__ad4zn (unmodified receptor, stored maps, cpu=$CPU, stored run: cpu 3)"
  "$PY" "$M/run_ad4zn.py" --target T3_6RL9 --ligand sfa_neutral --stage dock --seeds 10 \
    --exhaustiveness "$EXH" --cpu "$CPU" > "$SCRATCH/logs/baseline__sfa_neutral__ad4zn.log" 2>&1
  say "DONE  baseline rerun sfa_neutral__ad4zn | $(tail -1 "$SCRATCH/logs/baseline__sfa_neutral__ad4zn.log")"
  "$PY" - "$M" <<'EOF' | tee -a "$LOG"
import json, pathlib, sys
M = pathlib.Path(sys.argv[1]); D = pathlib.Path('dock')
tag = 'sfa_neutral__ad4zn'
new = json.loads((M / 'T3_6RL9/runs' / tag / 'results.json').read_text())
old = json.loads((D / 'T3_6RL9/runs' / tag / 'results.json').read_text())
rows = []
for sn, so in zip(new['seeds'], old['seeds']):
    a, b = M / 'T3_6RL9/runs' / tag / f"seed{sn['seed']:02d}", D / 'T3_6RL9/runs' / tag / f"seed{so['seed']:02d}"
    rows.append(dict(seed=sn['seed'],
                     ligand_pdbqt_identical=(a / 'ligand.pdbqt').read_bytes() == (b / 'ligand.pdbqt').read_bytes(),
                     out_pdbqt_identical=(a / 'out.pdbqt').read_bytes() == (b / 'out.pdbqt').read_bytes(),
                     scores_identical=sn['scores'] == so['scores']))
out = dict(purpose='rerun of the stored AutoDock4Zn run of sulfanilamide in the UNMODIFIED receptor T3_6RL9 '
                   'with the stored maps', run=tag,
           stored=dict(cpu=old['cpu'], exhaustiveness=old['exhaustiveness'], best_mean=old['summary']['best_mean'],
                       best_sd=old['summary']['best_sd'], wall_s_total=old['summary']['wall_s_total']),
           rerun=dict(cpu=new['cpu'], exhaustiveness=new['exhaustiveness'], best_mean=new['summary']['best_mean'],
                      best_sd=new['summary']['best_sd'], wall_s_total=new['summary']['wall_s_total']),
           n_seeds=len(rows),
           all_ligand_inputs_identical=all(r['ligand_pdbqt_identical'] for r in rows),
           all_out_pdbqt_identical=all(r['out_pdbqt_identical'] for r in rows),
           all_pose_scores_identical=all(r['scores_identical'] for r in rows), per_seed=rows)
out['reproduced'] = out['all_out_pdbqt_identical'] and out['all_ligand_inputs_identical']
(D / 'T3_6RL9_glu106GLU' / 'baseline_check_ad4zn.json').write_text(json.dumps(out, indent=2) + '\n')
print(f"baseline {tag}: inputs identical {out['all_ligand_inputs_identical']}, out.pdbqt identical "
      f"{out['all_out_pdbqt_identical']}, scores identical {out['all_pose_scores_identical']} "
      f"(stored cpu {old['cpu']}, rerun cpu {new['cpu']})")
EOF
}

stage_ad4zn () {
  # the runs of main-text Table 7 first
  job ad4zn T3_6RL9_glu106GLU sfa_neutral ad4zn 10 etkdg ""
  job ad4zn T4_3HS4_glu106GLU sfa_neutral ad4zn 10 etkdg ""
  job ad4zn T3_6RL9_glu106GLU sfa_neutral ad4zn 20 dft __prod
  job ad4zn T4_3HS4_glu106GLU sfa_neutral ad4zn 20 dft __prod
  job ad4zn T3_6RL9_glu106GLU sfa_anion   ad4zn 10 etkdg ""
  job ad4zn T4_3HS4_glu106GLU azm_neutral ad4zn 10 etkdg ""
  job ad4zn T3_6RL9_glu106GLU azm_neutral ad4zn 10 etkdg ""
  job ad4zn T4_3HS4_glu106GLU azm_anion   ad4zn 10 etkdg ""
  say "AD4ZN STAGE DONE"
}

stage_ad4zn_diag () {
  for t in $VARIANTS; do
    job ad4zn "$t" sfa_neutral ad4zn 10 dft __csn102 dock/results/diagnostics/sfa_neutral_csn102.xyz
    job ad4zn "$t" sfa_neutral ad4zn 10 dft __csn109 dock/results/diagnostics/sfa_neutral_etkdg_csn109.xyz
    job ad4zn "$t" sfa_neutral ad4zn 10 dft __etkdg  dock/results/diagnostics/sfa_neutral_etkdg_seed1.xyz
  done
  say "AD4ZN DIAGNOSTIC STAGE DONE"
}

# Fixed-pose rescoring (no search): the best near-native and the best ring-rotated mode of every seed
# of the STORED AutoDock4Zn runs of sulfanilamide, scored once in the reported maps and once in the
# maps of the variant (vina --score_only --scoring ad4).  It separates what the changed maps do to a
# given pose from what the search then does with them.  Not a number of the reported protocol.
stage_rescore () {
  mkdir -p "$SCRATCH/rescore"
  nice -n 10 "$PY" - "$SCRATCH/rescore" <<'EOF' 2>&1 | grep -v 'tagged as 2D' | tee -a "$LOG"
import json, pathlib, re, subprocess, sys
import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.warning')
sys.path.insert(0, 'dock')
import analyse as A
S, D = pathlib.Path(sys.argv[1]), pathlib.Path('dock').resolve()

def models(p):
    out, cur = [], None
    for ln in p.read_text().splitlines():
        if ln.startswith('MODEL'):
            cur = []
        elif ln.startswith('ENDMDL'):
            out.append(cur); cur = None
        elif cur is not None:
            cur.append(ln)
    return out

def score(lines, maps):
    f = S / 'pose.pdbqt'
    f.write_text('\n'.join(lines) + '\n')
    r = subprocess.run(['vina', '--ligand', str(f), '--maps', str(maps), '--scoring', 'ad4', '--score_only',
                        '--cpu', '1'], capture_output=True, text=True)
    return float(re.search(r'Estimated Free Energy of Binding\s*:\s*(-?[\d.]+)', r.stdout).group(1))

REF = {'T3_6RL9': 'ref_SAN.sdf', 'T4_3HS4': 'ref_cross_T3_6RL9_SAN.sdf'}   # stored references, read only
for t in ('T3_6RL9', 'T4_3HS4'):
    ref = Chem.MolFromMolFile(str(D / t / REF[t]))
    out = dict(target=t, reference=REF[t], what='vina --score_only --scoring ad4 on stored modes; '
               'near-native: RMSD <= 1.0 A, ring-rotated: > 2.0 A; best mode of each kind per seed',
               maps=dict(reported=f'{t}/ad4zn/receptor_tz', variant=f'{t}_glu106GLU/ad4zn/receptor_tz'), runs={})
    for tag in ('sfa_neutral__ad4zn', 'sfa_neutral__ad4zn__prod'):
        rows = []
        for sd in sorted((D / t / 'runs' / tag).glob('seed*')):
            pm = A.load_poses(sd / 'out.pdbqt')
            rr = A.symm_rmsd(ref, pm, [c.GetId() for c in pm.GetConformers()])
            ms = models(sd / 'out.pdbqt')
            for kind, idx in (('near-native', [i for i, v in enumerate(rr) if v <= 1.0][:1]),
                              ('ring-rotated', [i for i, v in enumerate(rr) if v > 2.0][:1])):
                for i in idx:
                    stored = float(next(l for l in ms[i] if l.startswith('REMARK VINA RESULT')).split()[3])
                    b = score(ms[i], D / t / 'ad4zn' / 'receptor_tz')
                    v = score(ms[i], D / f'{t}_glu106GLU' / 'ad4zn' / 'receptor_tz')
                    rows.append(dict(seed=sd.name, kind=kind, mode=i + 1, rmsd=round(rr[i], 3), stored_score=stored,
                                     rescored_reported_maps=b, rescored_variant_maps=v, change=round(v - b, 3)))
        summ = {}
        for kind in ('near-native', 'ring-rotated'):
            x = [r for r in rows if r['kind'] == kind]
            if x:
                d = np.array([r['change'] for r in x])
                summ[kind] = dict(n=len(x), change_mean=round(float(d.mean()), 3), change_min=float(d.min()),
                                  change_max=float(d.max()),
                                  max_abs_rescored_minus_stored=round(max(abs(r['rescored_reported_maps'] - r['stored_score']) for r in x), 3))
        if len(summ) == 2:
            summ['ring_rotated_minus_near_native_change'] = round(
                summ['ring-rotated']['change_mean'] - summ['near-native']['change_mean'], 3)
        out['runs'][tag] = dict(summary=summ, modes=rows)
        print(f"rescore {t} {tag}: " + '; '.join(
            f"{k} n={v['n']} change {v['change_mean']:+.3f} ({v['change_min']:+.3f} to {v['change_max']:+.3f})"
            for k, v in summ.items() if isinstance(v, dict))
            + (f"; rotated minus near-native {summ['ring_rotated_minus_near_native_change']:+.3f}"
               if 'ring_rotated_minus_near_native_change' in summ else ''))
    (D / f'{t}_glu106GLU' / 'rescore_fixed_poses.json').write_text(json.dumps(out, indent=2) + '\n')
EOF
}

# The same for Vina and Vinardo: the top pose of every seed of the STORED sulfanilamide runs, scored
# without a search in the reported receptor and in the variant (vina --score_only).  If a fixed pose
# scores the same in both, any difference between the docking runs comes from the search.
stage_rescore_vina () {
  mkdir -p "$SCRATCH/rescore"
  nice -n 10 "$PY" - "$SCRATCH/rescore" <<'EOF' 2>&1 | tee -a "$LOG"
import json, pathlib, re, subprocess, sys
S, D = pathlib.Path(sys.argv[1]), pathlib.Path('dock').resolve()

def top_model(p):
    out = []
    for ln in p.read_text().splitlines():
        if ln.startswith('ENDMDL'):
            break
        if not ln.startswith('MODEL'):
            out.append(ln)
    return out

def score(lines, receptor, scoring):
    f = S / 'pose_vina.pdbqt'
    f.write_text('\n'.join(lines) + '\n')
    cmd = ['vina', '--receptor', str(receptor), '--ligand', str(f), '--score_only', '--autobox', '--cpu', '1']
    if scoring != 'vina':
        cmd += ['--scoring', scoring]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return float(re.search(r'Estimated Free Energy of Binding\s*:\s*(-?[\d.]+)', r.stdout).group(1))

for t in ('T3_6RL9', 'T4_3HS4'):
    out = dict(target=t, what='vina --score_only on the top pose of every seed of the stored runs, in the '
               'reported receptor.pdbqt and in the variant receptor.pdbqt (--autobox: grid around the pose)', runs={})
    for tag in ('sfa_neutral__vina', 'sfa_neutral__vina__prod', 'sfa_neutral__vinardo__prod', 'sfa_anion__vina__prod'):
        scoring = tag.split('__')[1]
        rows = []
        for sd in sorted((D / t / 'runs' / tag).glob('seed*')):
            pose = top_model(sd / 'out.pdbqt')
            stored = float(next(l for l in pose if l.startswith('REMARK VINA RESULT')).split()[3])
            b = score(pose, D / t / 'receptor.pdbqt', scoring)
            v = score(pose, D / f'{t}_glu106GLU' / 'receptor.pdbqt', scoring)
            rows.append(dict(seed=sd.name, stored_score=stored, rescored_reported_receptor=b,
                             rescored_variant_receptor=v, change=round(v - b, 3)))
        ch = [r['change'] for r in rows]
        out['runs'][tag] = dict(n=len(rows), change_min=min(ch), change_max=max(ch),
                                max_abs_rescored_minus_stored=round(max(abs(r['rescored_reported_receptor'] - r['stored_score']) for r in rows), 3),
                                poses=rows)
        print(f"rescore {t} {tag}: n={len(rows)} change {min(ch):+.3f} to {max(ch):+.3f} "
              f"(rescored minus stored, reported receptor: max {out['runs'][tag]['max_abs_rescored_minus_stored']})")
    (D / f'{t}_glu106GLU' / 'rescore_fixed_poses_vina.json').write_text(json.dumps(out, indent=2) + '\n')
EOF
}

stage_compare () {
  for t in T3_6RL9 T4_3HS4; do
    nice -n 10 "$PY" dock/glu_protonation_variant.py --target "$t" --resnum 106 --compare \
      --scratch "$SCRATCH/compare" ${COMPARE_ARGS:-} 2>&1 | grep -v 'tagged as 2D' | tee -a "$LOG"
  done
}

for stage in "$@"; do
  case "$stage" in
    variants) stage_variants;; baseline) stage_baseline;; vina) stage_vina;;
    autogrid-check) stage_autogrid_check;; maps) stage_maps;; ad4zn-baseline) stage_ad4zn_baseline;;
    ad4zn) stage_ad4zn;;
    ad4zn-diag) stage_ad4zn_diag;; rescore) stage_rescore;; rescore-vina) stage_rescore_vina;;
    compare) stage_compare;;
    *) echo "unknown stage $stage"; exit 2;;
  esac
done
