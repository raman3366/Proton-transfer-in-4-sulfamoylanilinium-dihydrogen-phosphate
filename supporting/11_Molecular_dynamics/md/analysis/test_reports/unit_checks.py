# throwaway (versions 1.2 and 1.3): rules of sections 9A, 9B, 9D and 9E that cannot be planted in stored frames (a
# signed volume of exactly zero) or that sit at an edge (first passage at the first frame, a run shorter than 1 ns,
# frames of a short stage; version 1.3: the settings of the water count of section 9D from a manifest).  The functions of analyse_md.py are called on arrays typed here; the expected values are typed by hand
# and are the worked examples of SPEC.md.   usage: unit_checks.py ANALYSIS_DIR OUT.json
import hashlib, json, sys
from pathlib import Path
A = Path(sys.argv[1])
sys.path.insert(0, str(A))
import analyse_md as m
res, ok = [], True


def check(name, got, want):
    global ok
    good = got == want
    ok &= good
    res.append({'check': name, 'returned': got, 'expected': want, 'passed': good})


# 1. zero volume keeps the sign before it; the stretch of 99 frames is short, the one of 100 frames long
vol = [0.5] * 3 + [0.0] * 2 + [-0.4] * 120 + [0.0] + [0.3] * 99 + [-0.2] * 100
state = ['C'] * len(vol)
assigned, d = m.n2_configuration(vol, state, 1, 0.01)
check('example 1 (325 frames, start sign +1): frames positive, negative, zero', [d['count_positive'], d['count_negative'], d['count_zero']], [102, 220, 3])
check('example 1: assigned sign of frames 3, 4 (zero after +) and 125 (zero after -)', [int(assigned[3]), int(assigned[4]), int(assigned[125])], [1, 1, -1])
check('example 1: stretches, long stretches', [d['n_stretches'], d['n_long_stretches']], [4, 2])
check('example 1: long stretches (sign, first frame, frames)', [[x['sign'], x['first_frame'], x['n_frames']] for x in d['long_stretches']], [[-1, 5, 121], [-1, 225, 100]])
check('example 1: lasting sign changes (frame, time, from, to)', [[c['frame'], c['time_ns'], c['from_sign'], c['to_sign']] for c in d['lasting_sign_changes']], [[5, 0.06, 1, -1]])
# 2. zero volume in the first frames takes the sign of the start
assigned, d = m.n2_configuration([0.0, 0.0, 0.2] + [0.2] * 120, ['R'] * 123, -1, 0.01)
check('example 2 (start sign -1, two frames of zero volume first): assigned signs of frames 0, 1, 2', [int(v) for v in assigned[:3]], [-1, -1, 1])
check('example 2: lasting sign changes (frame, time, from, to)', [[c['frame'], c['time_ns'], c['from_sign'], c['to_sign']] for c in d['lasting_sign_changes']], [[2, 0.03, -1, 1]])
check('example 2: sign of the first frame as written (0), count_zero', [d['sign_first_frame'], d['count_zero']], [0, 2])
# 3. a run shorter than 1 ns has no long stretch and no lasting change
assigned, d = m.n2_configuration([-0.3] * 60, ['R'] * 60, 1, 0.01)
check('example 3 (60 frames of the other sign): long stretches, lasting changes, window frames', [d['n_long_stretches'], d['n_lasting_sign_changes'], d['window_frames']], [0, 0, 60])
# 4. first passage at the first frame: nothing before it
import numpy as np
strict = np.array(['R'] * 150 + ['C'] * 50)
_, o = m.orientation(strict, 'C', 0.01, 200)
h = m.hbond_first_passage([1, 0] * 100, o)
check('example 4 (start state C, 150 frames R then 50 frames C): first passage frame and time', [o['first_passage_frame'], o['first_passage_ns']], [0, 0.01])
check('example 4: before (frames, count, fraction) and after (frames, count, fraction)',
      [h['n_before'], h['count_before'], h['fraction_before'], h['n_after'], h['count_after'], h['fraction_after']], [0, 0, None, 200, 100, 0.5])
# 5. no first passage
_, o = m.orientation(np.array(['C'] * 200), 'C', 0.01, 200)
h = m.hbond_first_passage([1] * 200, o)
check('example 5 (no first passage): available, and the six numbers', [h['available'], h['n_before'], h['count_before'], h['fraction_before'], h['n_after'], h['count_after'], h['fraction_after']],
      [False, None, None, None, None, None, None])
# 6. frames of an equilibration stage
check('frames expected in a stage: 0-100 ps at 10 ps; 5-20 ps at 10 ps; 0-5 ps at 10 ps; 75-100 ps at 10 ps; 0-0.6 ps at 0.1 ps; 100-400 ps at 100 ps',
      [m.expected_frames(0, 100, 10), m.expected_frames(5, 20, 10), m.expected_frames(0, 5, 10), m.expected_frames(75, 100, 10), m.expected_frames(0, 0.6, 0.1), m.expected_frames(100, 400, 100)],
      [10, 1, 0, 2, 6, 3])
tab, src = m.equilibration_stage_table({})
check('stage table without a manifest: source, steps, end of the last stage', [src, [m.step_tag(c) for c in sorted(tab)], max(e for _, e in tab.values())],
      ['protocol (constants of the specification)', ['eq01', 'eq02', 'eq03', 'eq04', 'eq05', 'eq06'], 2000.0])
tab, src = m.equilibration_stage_table({'equilibration': {'stages': [{'tag': 'eq01', 'time_start_ps': 0.0, 'time_end_ps': 5.0}, {'tag': 'eq02', 'time_start_ps': 5.0, 'time_end_ps': 20.0}]}})
check('stage table from a manifest with two stages: source and table', [src, {m.step_tag(k): list(v) for k, v in tab.items()}], ['manifest equilibration.stages', {'eq01': [0.0, 5.0], 'eq02': [5.0, 20.0]}])
# 7. extension steps: their place in the order of the run, and tags that are not steps
man = {'equilibration': {'stages': [{'tag': t, 'time_start_ps': a, 'time_end_ps': b} for t, a, b in
                                    (('eq04', 1700.0, 2100.0), ('eq03x2', 1200.0, 1700.0), ('eq03', 400.0, 700.0), ('eq03x1', 700.0, 1200.0),
                                     ('eq03x0', 0.0, 1.0), ('eq3b', 0.0, 1.0), ('min', 0.0, 0.0), ('eq10', 2100.0, 2200.0))]}}
tab, src = m.equilibration_stage_table(man)
check('stage table from a manifest with extension steps: the steps in the order of the run (eq03x0, eq3b and min are no steps)',
      [m.step_tag(c) for c in sorted(tab)], ['eq03', 'eq03x1', 'eq03x2', 'eq04', 'eq10'])
check('stage and extension number of eq03x2, of eq04 and of eq10', [[c // 1000, c % 1000] for c in (m.step_code(3, 2), m.step_code(4), m.step_code(10))], [[3, 2], [4, 0], [10, 0]])
import tempfile, os
with tempfile.TemporaryDirectory() as d:
    for name in ('solute_eq03.dcd', 'solute_eq03x2.dcd', 'solute_eq03x1.dcd', 'solute_eq04.dcd', 'solute_eq03x0.dcd', 'solute_eq3b.dcd', 'solute_seg0001.dcd', 'full_eq03x1.dcd', 'solute_eq03x1.dcd.tmp'):
        open(os.path.join(d, name), 'w').close()
    check('files taken as equilibration steps of the solute, in the order of the run', [p.name for _, p in m.equilibration_files(d, 'solute')],
          ['solute_eq03.dcd', 'solute_eq03x1.dcd', 'solute_eq03x2.dcd', 'solute_eq04.dcd'])
hc = m.driver_hydration_checks({'equilibration': {'stages': [{'tag': 'eq03', 'hydration_check': {'n_water_O': 0, 'nearest_water_O': {'distance_A': 3.61}, 'applies': True, 'decision': 'extend'}},
                                                               {'tag': 'eq03x1', 'hydration_check': {'n_water_O': 2}}, {'tag': 'eq04'}]}})
check('count of the driver read from a manifest: steps and counts', {m.step_tag(k): v['n_water_O'] for k, v in hc.items()}, {'eq03': 0, 'eq03x1': 2})
# 8. version 1.3, section 9D: site atoms and cut-off of the count that the driver acts on
ref_lig = {a: None for a in m.LIG_HEAVY}
SIX, SEVEN = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1']


def settings(man):
    w = []
    s = m.ring_place_settings(man, ref_lig, w.append)
    return [s['site_atoms'], s['water_cutoff_A'], s['source'], s['settings_are_the_defaults'], len(w)]


def hc(**kw):
    return {'hydration_check': {'requested': True, 'settings': kw}}


check('settings of section 9D without a manifest: six ring carbons, 2.5 A, defaults, no warning', settings({}),
      [SIX, 2.5, 'specification (defaults)', True, 0])
check('a manifest whose protocol had no hydration block (no settings): defaults, no warning',
      settings({'hydration_check': {'requested': False, 'status': 'not made: the protocol has no block equilibration.hydration_check'}}),
      [SIX, 2.5, 'specification (defaults)', True, 0])
check('settings of the manifest equal to the defaults (as md/run_md.py writes them now): taken from the manifest, no warning',
      settings(hc(after_step='eq03', site_atoms=SIX, water_cutoff_A=2.5)), [SIX, 2.5, 'manifest hydration_check.settings', True, 0])
check('the six ring carbons in another order: still the defaults, the order of the manifest kept',
      settings(hc(site_atoms=SIX[::-1], water_cutoff_A=2.5)), [SIX[::-1], 2.5, 'manifest hydration_check.settings', True, 0])
check('settings of the driver before 8:01 AM of 6 Oct 2026 (seven atoms, 3.5 A): followed, one warning',
      settings(hc(site_atoms=SEVEN, water_cutoff_A=3.5)), [SEVEN, 3.5, 'manifest hydration_check.settings', False, 1])
check('settings of the protocol between 8:01 and 8:56 AM (six ring carbons, 3.5 A): followed, one warning',
      settings(hc(site_atoms=SIX, water_cutoff_A=3.5)), [SIX, 3.5, 'manifest hydration_check.settings', False, 1])
for what, kw in (('a site atom that is no heavy atom of the ligand', dict(site_atoms=SIX + ['C9'], water_cutoff_A=2.5)),
                 ('a site atom named twice', dict(site_atoms=SIX + ['C1'], water_cutoff_A=2.5)),
                 ('an empty list of site atoms', dict(site_atoms=[], water_cutoff_A=2.5)),
                 ('a cut-off of zero', dict(site_atoms=SIX, water_cutoff_A=0.0)),
                 ('a cut-off that is text', dict(site_atoms=SIX, water_cutoff_A='2.5')),
                 ('no cut-off', dict(site_atoms=SIX)), ('no site atoms', dict(water_cutoff_A=2.5))):
    check(f'settings of the manifest that cannot be used ({what}): both defaults, one warning', settings(hc(**kw)),
          [SIX, 2.5, 'specification (defaults)', True, 1])
out = {'what': 'rules of sections 9A, 9B, 9D and 9E on arrays typed by hand (unit_checks.py)', 'spec': m.SPEC_VERSION,
       'analysis_script_sha256': hashlib.sha256((A / 'analyse_md.py').read_bytes()).hexdigest(),
       'n_checks': len(res), 'n_passed': sum(r['passed'] for r in res), 'all_passed': bool(ok), 'checks': res}
Path(sys.argv[2]).write_text(json.dumps(out, indent=1) + '\n')
for r in res:
    print('passed' if r['passed'] else 'FAILED', '|', r['check'], '|', r['returned'])
print(f"{out['n_passed']} of {out['n_checks']} unit checks passed")
sys.exit(0 if ok else 1)
