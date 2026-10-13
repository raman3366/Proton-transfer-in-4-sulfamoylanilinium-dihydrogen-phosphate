#!/bin/bash
# throwaway: every check of the analysis (version 1.3) in one go: synthetic test, start check, crystal site, unit
# checks, the built systems as one-frame runs, version 1.2 against 1.3, repeat run, edge cases, runs of the driver
# (its first test run with equilibration files; two runs of the driver of today with its hydration check), broken
# copies.  One thread, niced.   A: the folder of analyse_md.py and make_test_trajectory.py (default: md/analysis).
set -u
export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1      # no bytecode file is left in md/analysis
S=/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/mdgr2/spec13; A=${A:-/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md/analysis}; PY=/root/cqc-sfa-tools/md/bin/python
V12=/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/mdprep3/spec12/test
HYD=/root/cqc-sfa-md/test_driver/hydration
export A T=$S/test
cd $S
echo "== started $(TZ=Asia/Kolkata date '+%d %b %Y, %I:%M %p IST')"
echo "== runs of the campaign on this machine (none may exist)"; $PY runs_state.py
echo "== scripts"; sha256sum $A/analyse_md.py $A/make_test_trajectory.py; echo "driver and protocol:"; sha256sum $A/../run_md.py $A/../protocol.json
echo "== the driver as it stands against the driver that version 1.2 named (commit 7a7a988)"
REPO=/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate
echo "md/run_md.py of commit 7a7a988: sha256 $(git -C $REPO show 7a7a988:md/run_md.py | sha256sum | cut -d' ' -f1)"
git -C $REPO diff 7a7a988 -- md/run_md.py > $S/driver_diff.txt
echo "lines added: $(grep -c '^+[^+]' $S/driver_diff.txt); lines removed: $(grep -c '^-[^-]' $S/driver_diff.txt)"; grep '^[+-][^+-]' $S/driver_diff.txt | cut -c1-200
echo "== synthetic test"; ./run_test.sh
echo "== start check"; nice -n 10 $PY $A/analyse_md.py --check-starts --out $S/check_starts
echo "== crystal site by plain formulas"; nice -n 10 $PY crystal_site.py $A/../system $S/crystal_site.json
echo "== unit checks"; nice -n 10 $PY unit_checks.py $A $S/unit_checks.json
echo "== the built systems as one-frame runs"; rm -rf built; nice -n 10 $PY built_starts.py $S $A
echo "== version 1.2 against version 1.3 on the runs of version 1.2"; nice -n 10 $PY compare_v12_v13.py $V12 $S/test $S/v12_against_v13.json | tail -n 24
echo "== repeat run"; rm -rf repro
for r in runs/test/X1 runs_gv/edge/X0 runs_v12/conv/R1 runs_v13/place/R1; do
  k=$(basename $(dirname $r))_$(basename $r)
  nice -n 10 $PY $A/analyse_md.py test/$r --out repro/$k/a > /dev/null
  nice -n 10 $PY $A/analyse_md.py test/$r --out repro/$k/b > /dev/null
  (cd repro/$k/a && sha256sum * > ../a.sha); (cd repro/$k/b && sha256sum * > ../b.sha); (cd test/$r/analysis && sha256sum * > $S/repro/$k/c.sha)
  diff repro/$k/a.sha repro/$k/b.sha && diff repro/$k/a.sha repro/$k/c.sha && echo "analysis outputs of $r byte-identical in three runs: $(wc -l < repro/$k/a.sha) files"
done
echo "== edge cases"; rm -rf edge; mkdir -p edge/runs/test/X5 edge/runs/test/X6 edge/runs/test/Q7 edge/camp/test/X1 edge/camp/test/X2 edge/camp/test/X0 edge/runs/test/X0; cp -r test/systems edge/
for f in solute.pdb full.pdb manifest.json solute_seg0001.dcd; do cp test/runs/test/X2/$f edge/runs/test/X5/; cp test/runs/test/X2/$f edge/runs/test/X6/; done
cp test/runs/test/X2/full_seg0001.dcd edge/runs/test/X6/
head -c 276 test/runs/test/X2/solute_seg0001.dcd > edge/runs/test/X6/solute_seg0002.dcd
head -c 276 test/runs/test/X2/full_seg0001.dcd > edge/runs/test/X6/full_seg0002.dcd
: > edge/runs/test/X6/solute_seg0003.dcd
cp test/runs/test/X2/solute.pdb edge/runs/test/Q7/
echo "-- X5: no full trajectory, planned 5 ns"; nice -n 10 $PY $A/analyse_md.py edge/runs/test/X5 --planned-ns 5 --no-figures; echo "exit $?"
$PY -c "import json; s=json.load(open('edge/runs/test/X5/analysis/summary.json')); g=s['gate_GV']; f=g['GV3_zinc']['criteria']['fifth_ligand']; print('   X5: full fraction available', f['full_trajectory']['available'], '| fifth-ligand criterion met', f['passed'], '| GV3 met', g['GV3_zinc']['passed'], '| verdict', g['verdict'], '| water in the place of the ring, frames', s['stats']['n_water_ring_place']['n'], 'mean', s['stats']['n_water_ring_place']['mean'], '| wider count, frames', s['stats']['n_water_ring_wide_3p5A']['n'])"
echo "-- X6: last segments hold a header only or nothing"; nice -n 10 $PY $A/analyse_md.py edge/runs/test/X6 --no-figures; echo "exit $?"
echo "-- Q7: folder name is not a start, no manifest"; nice -n 10 $PY $A/analyse_md.py edge/runs/test/Q7 --no-figures; echo "exit $?"
echo "-- Q7 with --start X: no trajectory"; nice -n 10 $PY $A/analyse_md.py edge/runs/test/Q7 --start X --system edge/systems/test/X --no-figures; echo "exit $?"
echo "-- X5 with start R on the command line (folder says X)"; nice -n 10 $PY $A/analyse_md.py edge/runs/test/X5 --start R --no-figures --out edge/out_R; echo "exit $?"
$PY -c "import json; s=json.load(open('edge/out_R/summary.json')); print('   start R on the command line: sign of the start', s['n2_configuration']['sign_of_start'], '| first frame differs in state', s['start_of_production']['state_differs'], 'and in sign', s['start_of_production']['n2_sign_differs'])"
echo "-- version 1.1: a folder X0 without a full trajectory (files of X2: 30 ns, plain manifest with replica 2 and 30 ns)"
for f in solute.pdb manifest.json solute_seg0001.dcd; do cp test/runs/test/X2/$f edge/runs/test/X0/; done
nice -n 10 $PY $A/analyse_md.py edge/runs/test/X0 --no-figures; echo "exit $?"
$PY -c "import json; s=json.load(open('edge/runs/test/X0/analysis/summary.json')); print('   folder X0 with manifest replica 2: role', s['run']['role'], '| replica', s['run']['replica'], 'from', s['run']['sources']['replica'], '| planned', s['run']['planned_ns'], 'ns')"
echo "-- version 1.1: the validation run with --planned-ns 2 on the command line"; nice -n 10 $PY $A/analyse_md.py test/runs/test/X0 --planned-ns 2 --no-figures --out edge/out_X0_2ns; echo "exit $?"
$PY -c "import json; s=json.load(open('edge/out_X0_2ns/summary.json')); print('   role', s['run']['role'], '| planned', s['run']['planned_ns'], 'ns from', s['run']['sources']['planned_ns'], '| warnings', len(s['warnings']), '| verdict', s['gate_GV']['verdict'])"
echo "-- version 1.1: a campaign folder with X1 (really the validation run: its manifest says X0), X2 and X0"
for f in solute.pdb full.pdb manifest.json solute_seg0001.dcd full_seg0001.dcd; do cp test/runs/test/X0/$f edge/camp/test/X1/; cp test/runs/test/X0/$f edge/camp/test/X0/; cp test/runs/test/X2/$f edge/camp/test/X2/; done
nice -n 10 $PY $A/analyse_md.py --campaign edge/camp --no-figures; echo "exit $?"
$PY -c "import json; c=json.load(open('edge/camp/analysis_campaign/campaign.json')); print('   table:', [r['table_row']['run'] for r in c['runs']], '| left out:', ['/'.join(x.split('/')[-2:]) for x in c['skipped']], '| failed:', len(c['failed']))"
echo "-- version 1.1: the campaign of the test again with --reuse"; nice -n 10 $PY $A/analyse_md.py --campaign test/runs --reuse --no-figures --out edge/reuse | cut -c1-60
cmp edge/reuse/campaign_table.csv test/runs/analysis_campaign/campaign_table.csv && echo "   campaign table with --reuse identical"
echo "-- --reuse does not take a summary of another version (a summary that names version 1.2)"
mkdir -p edge/old/test/X2/analysis; for f in solute.pdb full.pdb manifest.json solute_seg0001.dcd full_seg0001.dcd; do cp test/runs/test/X2/$f edge/old/test/X2/; done
sed 's/version 1.3 (6 Oct 2026)/version 1.2 (6 Oct 2026)/' test/runs/test/X2/analysis/summary.json > edge/old/test/X2/analysis/summary.json
nice -n 10 $PY $A/analyse_md.py --campaign edge/old --reuse --no-figures | cut -c1-60
echo "-- version 1.2: conv/R1 without its equilibration files (the manifest still lists the stages)"
mkdir -p edge/noeq/conv/R1; for f in solute.pdb full.pdb manifest.json; do cp test/runs_v12/conv/R1/$f edge/noeq/conv/R1/; done; cp test/runs_v12/conv/R1/*_seg*.dcd edge/noeq/conv/R1/
nice -n 10 $PY $A/analyse_md.py edge/noeq/conv/R1 --no-figures; echo "exit $?"
$PY - <<'PYEOF'
import json, filecmp, os
a, b = 'test/runs_v12/conv/R1/analysis', 'edge/noeq/conv/R1/analysis'
same = [f for f in ('timeseries_solute.csv', 'timeseries_full.csv', 'rmsf_ca.csv') if filecmp.cmp(f'{a}/{f}', f'{b}/{f}', shallow=False)]
sa, sb = json.load(open(f'{a}/summary.json')), json.load(open(f'{b}/summary.json'))
diff = sorted(k for k in sa if sa[k] != sb[k])
row = sorted(k for k in sa['table_row'] if sa['table_row'][k] != sb['table_row'][k])
run = sorted(k for k in sa['run'] if sa['run'][k] != sb['run'][k])
print('   without the equilibration files:', len(same), 'of 3 CSV files identical | keys of summary.json that differ:', diff, '| of table_row:', row, '| of run:', run,
      '| equilibration available', sb['equilibration']['available'], '| equilibration series written', os.path.exists(f'{b}/timeseries_equil_solute.csv'), '| warnings', len(sb['warnings']))
PYEOF
echo "-- version 1.2: a stage file that the stage table does not have (solute_eq07.dcd beside six stages)"
mkdir -p edge/eq07/scaled; cp -r test/runs_v12/scaled/N1 edge/eq07/scaled/N1; rm -rf edge/eq07/scaled/N1/analysis; cp edge/eq07/scaled/N1/solute_eq06.dcd edge/eq07/scaled/N1/solute_eq07.dcd
nice -n 10 $PY $A/analyse_md.py edge/eq07/scaled/N1 --no-figures; echo "exit $?"
$PY -c "import json; s=json.load(open('edge/eq07/scaled/N1/analysis/summary.json')); e=s['equilibration']; print('   solute_eq07.dcd: warnings', len(s['warnings']), '| equilibration frames read', e['n_frames_solute'], '| stages', [a['stage'] for a in e['stages']])"
echo "-- version 1.2: the files of a shortened equilibration read without the stage list of the manifest"
mkdir -p edge/nolist/scaled; cp -r test/runs_v12/scaled/N1 edge/nolist/scaled/N1; rm -rf edge/nolist/scaled/N1/analysis
$PY -c "import json; p='edge/nolist/scaled/N1/manifest.json'; m=json.load(open(p)); del m['equilibration']; json.dump(m, open(p,'w'), indent=1)"
nice -n 10 $PY $A/analyse_md.py edge/nolist/scaled/N1 --no-figures; echo "exit $?"
$PY -c "import json; s=json.load(open('edge/nolist/scaled/N1/analysis/summary.json')); e=s['equilibration']; w=s['warnings']; print('   no stage list: table from', e['stage_table_source'], '| warnings', len(w), '| of them on the number of frames', sum('expected for a stage' in x for x in w), '| on the time axis', sum('time axis of the equilibration is wrong' in x for x in w))"
echo "-- version 1.2: one stage file cut down to its header, one emptied (conv/R1)"
mkdir -p edge/cut/conv; cp -r test/runs_v12/conv/R1 edge/cut/conv/R1; rm -rf edge/cut/conv/R1/analysis
head -c 276 test/runs_v12/conv/R1/solute_eq03.dcd > edge/cut/conv/R1/solute_eq03.dcd; : > edge/cut/conv/R1/solute_eq06.dcd
nice -n 10 $PY $A/analyse_md.py edge/cut/conv/R1 --no-figures; echo "exit $?"
$PY - <<'PYEOF'
import csv, json, filecmp
a, b = 'test/runs_v12/conv/R1/analysis', 'edge/cut/conv/R1/analysis'
s = json.load(open(f'{b}/summary.json')); e = s['equilibration']
t = [r for r in csv.DictReader(open(f'{b}/timeseries_equil_solute.csv'))]
o = {r['frame']: r for r in csv.DictReader(open(f'{a}/timeseries_equil_solute.csv'))}
kept = [r for r in csv.DictReader(open(f'{a}/timeseries_equil_solute.csv')) if r['stage'] not in ('3', '6')]
print('   cut files: warnings', len(s['warnings']), '| equilibration frames read', e['n_frames_solute'], '| first and last time', t[0]['time_ns'], t[-1]['time_ns'],
      '| times and values of the other stages unchanged', [(x['stage'], x['time_ns'], x['lig_rmsd_A'], x['n2_volume']) for x in t] == [(x['stage'], x['time_ns'], x['lig_rmsd_A'], x['n2_volume']) for x in kept],
      '| production series identical', filecmp.cmp(f'{a}/timeseries_solute.csv', f'{b}/timeseries_solute.csv', shallow=False))
PYEOF
echo "-- version 1.2: extension steps read without the stage list of the manifest (ext/R1)"
mkdir -p edge/extnolist/ext; cp -r test/runs_v12/ext/R1 edge/extnolist/ext/R1; rm -rf edge/extnolist/ext/R1/analysis
$PY -c "import json; p='edge/extnolist/ext/R1/manifest.json'; m=json.load(open(p)); del m['equilibration']; json.dump(m, open(p,'w'), indent=1)"
nice -n 10 $PY $A/analyse_md.py edge/extnolist/ext/R1 --no-figures 2>&1 | cut -c1-230; echo "exit $?"
$PY -c "import json; s=json.load(open('edge/extnolist/ext/R1/analysis/summary.json')); e=s['equilibration']; w=s['warnings']; print('   extension steps without a stage list: table from', e['stage_table_source'], '| warnings', len(w), '| on files that belong to no step', sum('belong to no step' in x for x in w), '| on the time axis', sum('time axis of the equilibration is wrong' in x for x in w), '| on the settings of the hydration check', sum('hydration check of this run' in x for x in w), '| solute frames read', e['n_frames_solute'], '| steps read', [a['tag'] for a in e['stages']])"
echo "-- version 1.2: the count of the driver in the manifest changed (ext/R1, eq03x2: 2 water oxygens instead of 1)"
mkdir -p edge/extcount/ext; cp -r test/runs_v12/ext/R1 edge/extcount/ext/R1; rm -rf edge/extcount/ext/R1/analysis
$PY -c "import json; p='edge/extcount/ext/R1/manifest.json'; m=json.load(open(p)); r=[x for x in m['equilibration']['stages'] if x['tag']=='eq03x2'][0]; r['hydration_check']['n_water_O']=2; json.dump(m, open(p,'w'), indent=1)"
nice -n 10 $PY $A/analyse_md.py edge/extcount/ext/R1 --no-figures 2>&1 | cut -c1-330; echo "exit $?"
$PY -c "import json; s=json.load(open('edge/extcount/ext/R1/analysis/summary.json')); e=s['equilibration']; w=s['warnings']; print('   changed count of the driver: warnings', len(w), '| on the count', sum('the driver counted' in x for x in w), '| on the settings', sum('hydration check of this run' in x for x in w), '| count agrees at the end of eq03, eq03x1, eq03x2:', [a['driver_hydration_check']['count_agrees'] for a in e['stages'] if a['driver_hydration_check']], '| nearest agrees:', [a['driver_hydration_check']['nearest_agrees'] for a in e['stages'] if a['driver_hydration_check']])"
echo "-- version 1.3: place/R1 with other entries under hydration_check of its manifest (the frames and the records of the steps unchanged)"
$PY - <<'PYEOF'
import json, os, shutil, subprocess, sys
A, PY = os.environ['A'], sys.executable
base = json.load(open('test/runs_v13/place/R1/analysis/summary.json'))
cases = (('nosettings', 'the protocol of the run had no hydration block (requested false, no settings)',
          lambda m: m.__setitem__('hydration_check', {'requested': False, 'status': 'not made: the protocol has no block equilibration.hydration_check'})),
         ('nokey', 'no key hydration_check at all', lambda m: m.pop('hydration_check')),
         ('badatom', 'a site atom that the ligand does not have (C9)', lambda m: m['hydration_check']['settings'].__setitem__('site_atoms', ['C1', 'C2', 'C3', 'C4', 'C5', 'C9'])),
         ('badcut', 'a cut-off of zero', lambda m: m['hydration_check']['settings'].__setitem__('water_cutoff_A', 0.0)),
         ('seven', 'seven site atoms and 3.5 A, which is not what the records of the steps were counted with',
          lambda m: m['hydration_check']['settings'].update(site_atoms=['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1'], water_cutoff_A=3.5)),
         ('sixwide', 'the six ring carbons with 3.5 A (the protocol between 8:01 and 8:56 AM of 6 Oct 2026)',
          lambda m: m['hydration_check']['settings'].update(water_cutoff_A=3.5)))
for name, what, change in cases:
    d = f'edge/place_{name}/place/R1'
    shutil.rmtree(f'edge/place_{name}', ignore_errors=True)
    shutil.copytree('test/runs_v13/place/R1', d, ignore=shutil.ignore_patterns('analysis'))
    m = json.load(open(f'{d}/manifest.json'))
    change(m)
    json.dump(m, open(f'{d}/manifest.json', 'w'), indent=1)
    r = subprocess.run(['nice', '-n', '10', PY, f'{A}/analyse_md.py', d, '--no-figures', '--system', 'test/systems/place/R'], capture_output=True, text=True)
    s = json.load(open(f'{d}/analysis/summary.json'))
    w, p, e = s['warnings'], s['water_ring_place'], s['equilibration']
    full = [x.split(',') for x in open(f'{d}/analysis/timeseries_full.csv').read().splitlines()]
    col = full[0].index('n_water_ring_place')
    print(f'   {name} ({what}): exit {r.returncode} | site atoms {",".join(p["site_atoms"])}, cut-off {p["water_cutoff_A"]:g} A, from: {p["source"]}, defaults: {p["settings_are_the_defaults"]}'
          f' | warnings {len(w)}: settings not usable {sum("without a usable" in x for x in w)}, settings not those of the plan {sum("hydration check of this run" in x for x in w)},'
          f' count of the driver {sum("the driver counted" in x for x in w)}'
          f' | count in the 20 production frames {"".join(x[col] for x in full[1:])}'
          f' | same production counts as with the manifest as written: {s["stats"]["n_water_ring_place"] == base["stats"]["n_water_ring_place"]}'
          f' | count agrees at the end of eq03, eq03x1, eq03x2, eq03x3: {[a["driver_hydration_check"]["count_agrees"] for a in e["stages"] if a["driver_hydration_check"]]}'
          f' | nearest agrees: {[a["driver_hydration_check"]["nearest_agrees"] for a in e["stages"] if a["driver_hydration_check"]]}')
PYEOF
echo "== the copy of a finished test run of the driver (5 Oct 2026, before its hydration check)"; rm -rf driver_snapshot/analysis; (cd driver_snapshot && sha256sum -c snapshot.sha256 | grep -vc ': OK$' | sed 's/^/files of the snapshot that changed: /'); nice -n 10 $PY $A/analyse_md.py driver_snapshot/runs/standin/ref --start X --system driver_snapshot/systems/standin --out driver_snapshot/analysis
echo "== two runs of the driver as it stands, with its hydration check (driver_now/run_driver.sh, made now; its log follows)"
(cd driver_now && $PY make_protocol.py $A/.. protocol_test_now.json > /dev/null && ./run_driver.sh > run_driver.log 2>&1); cut -c1-330 driver_now/run_driver.log
rm -rf driver_now/analysis
for s in dry wet; do nice -n 10 $PY $A/analyse_md.py driver_now/runs/$s/R1 --out driver_now/analysis/$s; done
$PY driver_crosscheck.py $S $A/..
echo "== the test runs of the driver's maintainer (6 Oct 2026, early morning: the driver as version 1.2 named it, hydration block of that hour: seven site atoms, 3.5 A), read only"
for s in dry arrive wet skip late; do
  nice -n 10 $PY $A/analyse_md.py $HYD/runs/$s --start R --out $S/driver_0528/$s --no-figures --no-hash-trajectories > $S/driver_0528_$s.console 2>&1; code=$?
  $PY - $s $code <<'PYEOF'
import json, sys
s, code = sys.argv[1], sys.argv[2]
S = json.load(open(f'driver_0528/{s}/summary.json'))
man = json.load(open(f'/root/cqc-sfa-md/test_driver/hydration/runs/{s}/manifest.json'))
w, p, e = S['warnings'], S['water_ring_place'], S['equilibration']
hs = man['hydration_check']['settings']
print(f'   {s}: exit {code} | driver {man["script"]["sha256"][:16]} | steps read {" ".join(a["tag"] for a in e["stages"])} | as in the manifest {[a["tag"] for a in e["stages"]] == [x["tag"] for x in man["equilibration"]["stages"]]}'
      f' | frames read as recorded {all([a["solute"]["n_frames"], a["full"]["n_frames"]] == [x["frames"]["solute"], x["frames"]["full"]] for a, x in zip(e["stages"], man["equilibration"]["stages"]))}'
      f' | settings used {",".join(p["site_atoms"])} {p["water_cutoff_A"]:g} A ({p["source"]}), as in the manifest {p["site_atoms"] == hs["site_atoms"] and p["water_cutoff_A"] == hs["water_cutoff_A"]}'
      f' | counts of the driver {[a["driver_hydration_check"]["n_water_O"] for a in e["stages"] if a["driver_hydration_check"]]}, compared {[a["driver_hydration_check"]["compared"] for a in e["stages"] if a["driver_hydration_check"]]}'
      f' | warnings {len(w)}: settings not those of the plan {sum("hydration check of this run" in x for x in w)}, DCD header against the stage table {sum("time axis of the equilibration is wrong" in x for x in w)}, other {sum(("hydration check of this run" not in x) and ("time axis of the equilibration is wrong" not in x) for x in w)}')
PYEOF
done
$PY driver_0528_facts.py $S $S/driver_0528.json | sed 's/^/   facts: /' | cut -c1-260
if [ "${SKIP_MUT:-0}" = 1 ]; then echo "== mutations: left out in this trial"; else
echo "== mutations"; rm -rf test_mut; cp -r test test_mut
nice -n 10 $PY mutation/mutate.py $S $A $PY
fi
echo "== scripts again"; sha256sum $A/analyse_md.py $A/make_test_trajectory.py; echo "driver and protocol:"; sha256sum $A/../run_md.py $A/../protocol.json
echo "== runs of the campaign on this machine (at the end)"; $PY runs_state.py; $PY runs_dates.py $S/runs_dates.json > /dev/null
echo "== finished $(TZ=Asia/Kolkata date '+%d %b %Y, %I:%M %p IST')"
echo ALLDONE
