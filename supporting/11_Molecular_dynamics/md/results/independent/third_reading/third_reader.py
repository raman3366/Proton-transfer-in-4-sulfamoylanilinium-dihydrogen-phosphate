#!/usr/bin/env python
"""Blind third reading of the MD campaign: the eleven 100 ns runs.

One process.  Reads solute.pdb, manifest.json and solute_*.dcd of each run
folder, nothing of any analysis folder.  Writes per-run series (npz) and
runs_raw.json into the scratch folder; third_reading.json is put together by
make_reading.py from runs_raw.json.
"""
import json
import os
import re
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import third_lib as T

HERE = os.path.dirname(os.path.abspath(__file__))
AKMA_FS = 48.88821            # one AKMA time unit in fs (the unit of the step length in a DCD header)
KEYS = ['d1', 'd2', 'd1_mi', 'd2_mi', 'zn_n2', 'zn_h94', 'zn_h96', 'zn_h119', 'zn_o1', 'zn_o2',
        'zn_omin', 'n2_og1', 'ang_n2_h_og1', 'hb', 'n2_vol', 'ca_rmsd', 'lig_rmsd',
        'lig_rmsd_identity', 'lig_rmsd_labelling', 'zn_moved', 'lig_moved']


def read_files(folder, names, sel, ref):
    """Read DCD files in the given order; return joined series and per-file facts."""
    series = {k: [] for k in KEYS}
    facts = []
    for name in names:
        xyz, lat, head = T.read_dcd(os.path.join(folder, name), sel['n_atoms'])
        head['file'] = name
        if xyz.shape[0]:
            q = T.frame_quantities(xyz, lat, sel, ref)
            for k in KEYS:
                series[k].append(q[k])
        facts.append(head)
        del xyz
    out = {}
    for k in KEYS:
        out[k] = np.concatenate(series[k]) if series[k] else np.zeros(0)
    return out, facts


def r(x, n=4):
    return None if x is None else round(float(x), n)


def stats(a):
    if len(a) == 0:
        return None
    return {'mean': r(a.mean()), 'sd': r(a.std()), 'min': r(a.min()), 'max': r(a.max())}


def visit_rec(v, t_of_first):
    return {'state': v[0], 'first_frame': int(v[1]), 'n_frames': int(v[2]), 't_first': t_of_first(v[1])}


def one_run(zset, run, start, role, ref):
    folder = os.path.join(T.RUNS, zset, run)
    sel = T.select(os.path.join(folder, 'solute.pdb'))
    man = T.load_manifest(folder)
    start_state = 'R' if start == 'R' else 'C'
    opp = 'C' if start_state == 'R' else 'R'
    start_sign = -1 if start == 'R' else 1

    res = {'set': zset, 'run': run, 'start': start, 'role': role, 'folder': folder,
           'state_of_start': start_state, 'opposite_state': opp}

    # ---------------- production
    segs = T.segment_files(folder)
    P, facts = read_files(folder, segs, sel, ref)
    N = len(P['d1'])
    hdr_ok = True
    for i, h in enumerate(facts):
        # first frame of segment i+1 expected at i*5000 + 10 ps; the header counts steps of 4 fs
        step_ps = h['delta_akma'] * AKMA_FS / 1000.0      # step length of the header (4 fs in the production)
        t_first_ps = h['first_step'] * step_ps
        if h['n_frames_read'] != h['nset_header'] or h['bytes_left_over'] != 0:
            hdr_ok = False
        if abs(h['interval_steps'] * step_ps - T.DT_PS) > 1e-4:
            hdr_ok = False
        exp = sum(f['n_frames_read'] for f in facts[:i]) * T.DT_PS + T.DT_PS
        if abs(t_first_ps - exp) > 1e-3 * max(1.0, exp / 1000.0):
            hdr_ok = False
    mp = man.get('production', {})
    res['frames'] = {
        'n_production_frames': int(N),
        'n_segment_files': len(segs),
        'frames_per_file': sorted(set(h['n_frames_read'] for h in facts)),
        'manifest_frames_solute': mp.get('frames_solute'),
        'manifest_completed_ns': mp.get('completed_ns'),
        'manifest_n_segments_complete': mp.get('n_segments_complete'),
        'manifest_status': man.get('status'),
        'dcd_headers_consistent_with_10ps_and_contiguous_time': bool(hdr_ok),
        'time_of_frame_n_ns': '(n + 1) * 0.01, n from 0',
    }
    st = T.states_of(P['d1'], P['d2'])
    half0 = N // 2                       # last half: frames N//2 ... N-1
    w = min(N, T.WINDOW)

    res['imaging_and_fit'] = {
        'frames_zinc_moved_by_a_lattice_vector': int(P['zn_moved'].sum()),
        'frames_ligand_moved_by_a_lattice_vector': int(P['lig_moved'].sum()),
        'max_abs_d1_minus_minimum_image_d1_A': r(np.abs(P['d1'] - P['d1_mi']).max(), 6),
        'max_abs_d2_minus_minimum_image_d2_A': r(np.abs(P['d2'] - P['d2_mi']).max(), 6),
        'ca_rmsd_to_reference_A': stats(P['ca_rmsd']),
    }
    res['d1_A'] = {'whole': stats(P['d1']), 'last_half': stats(P['d1'][half0:]), 'last_20ns': stats(P['d1'][N - w:])}
    res['d2_A'] = {'whole': stats(P['d2']), 'last_half': stats(P['d2'][half0:]), 'last_20ns': stats(P['d2'][N - w:])}
    res['frames_exactly_on_a_midpoint'] = int(np.sum((P['d1'] == T.M1) | (P['d2'] == T.M2)))
    res['states'] = {'whole': T.counts(st), 'last_half': T.counts(st[half0:]), 'last_20ns': T.counts(st[N - w:])}
    v, w_, c_start, c_opp = T.verdict(st, start_state)
    res['verdict'] = {'verdict': v, 'W': int(w_), 'c_start': c_start, 'c_opposite': c_opp,
                      'c_intermediate': int(w_ - c_start - c_opp),
                      'test_converted': '10*%d >= 9*%d : %s' % (c_opp, w_, 10 * c_opp >= 9 * w_),
                      'test_stayed': '10*%d >= 9*%d : %s' % (c_start, w_, 10 * c_start >= 9 * w_),
                      'window_complete': bool(N >= T.WINDOW)}

    # ---------------- visits in the production
    a, vis, longv, first, ntr = T.passage(st, start_state, opp, 'I')
    tprod = lambda n: r((n + 1) * 0.01, 2)
    a_arr = np.array(a)
    res['visits'] = {
        'assigned': {s: int(np.sum(a_arr == s)) for s in 'CR'},
        'n_visits': len(vis),
        'n_visits_by_state': {s: sum(1 for x in vis if x[0] == s) for s in 'CR'},
        'n_long_visits': len(longv),
        'long_visits': [dict(visit_rec(x, tprod),
                             fraction_strict=r(np.mean(st[x[1]:x[1] + x[2]] == x[0])))
                        for x in longv],
        'longest_visit_in_opposite_state_frames': max([x[2] for x in vis if x[0] == opp] or [0]),
        'longest_visit_in_start_state_frames': max([x[2] for x in vis if x[0] == start_state] or [0]),
        'first_passage_ns': (tprod(first[1]) if first else None),
        'first_passage_frame': (int(first[1]) if first else None),
        'n_transitions': int(ntr),
        'first_frame_state': str(st[0]),
        'first_frame_not_in_start_state': (lambda ix: {'frame': int(ix[0]), 't_ns': tprod(ix[0]), 'state': str(st[ix[0]])}
                                           if len(ix) else None)(np.nonzero(st != start_state)[0]),
        'first_frame_in_opposite_state': (lambda ix: {'frame': int(ix[0]), 't_ns': tprod(ix[0])}
                                          if len(ix) else None)(np.nonzero(st == opp)[0]),
    }

    # ---------------- equilibration, all steps with a solute file (for the record)
    eqfiles = sorted(f for f in os.listdir(folder) if re.fullmatch(r'solute_eq\d+(x[1-9]\d*)?\.dcd', f))

    def eq_key(f):
        m = re.fullmatch(r'solute_eq(\d+)(?:x(\d+))?\.dcd', f)
        return (int(m.group(1)), int(m.group(2) or 0))
    eqfiles.sort(key=eq_key)
    mstages = {s['tag']: s for s in man.get('equilibration', {}).get('stages', [])}
    eq_steps = []
    eq_series = {}
    for f in eqfiles:
        tag = f[len('solute_'):-len('.dcd')]
        E, fa = read_files(folder, [f], sel, ref)
        ste = T.states_of(E['d1'], E['d2'])
        ms = mstages.get(tag, {})
        S = ms.get('time_start_ps')
        hd = fa[0]
        hdr_step_ps = hd['delta_akma'] * AKMA_FS / 1000.0
        hdr_first_ps = hd['first_step'] * hdr_step_ps
        eq_steps.append({
            'tag': tag, 'frames': int(len(ste)), 'time_start_ps': S, 'time_end_ps': ms.get('time_end_ps'),
            'header_step_fs': r(hdr_step_ps * 1000.0, 3),
            'header_first_frame_ps': r(hdr_first_ps, 3),
            'header_agrees_with_S_plus_10ps': (None if S is None else bool(abs(hdr_first_ps - (S + T.DT_PS)) < 1e-2)),
            'restraint_atoms': ms.get('settings', {}).get('restraint_atoms'),
            'states': {s: int(np.sum(ste == s)) for s in 'CRI'},
            'n2_sign': {'plus': int(np.sum(E['n2_vol'] > 0)), 'minus': int(np.sum(E['n2_vol'] < 0))},
            'last_frame': ({'state': str(ste[-1]), 'd1': r(E['d1'][-1], 3), 'd2': r(E['d2'][-1], 3),
                            'lig_rmsd': r(E['lig_rmsd'][-1], 3), 'n2_vol': r(E['n2_vol'][-1], 3)} if len(ste) else None),
        })
        eq_series[tag] = (E, ste)
    res['equilibration_steps'] = eq_steps

    # ---------------- joined series: stages 4, 5, 6 then the production (plan, 7 Oct 12:10 AM, point 1)
    rel_tags = ['eq04', 'eq05', 'eq06']
    missing = [t for t in rel_tags if t not in eq_series]
    if missing:
        res['conversion'] = {'error': 'missing equilibration files: %s' % missing}
        n_eq = 0
    else:
        st_eq = np.concatenate([eq_series[t][1] for t in rel_tags])
        vol_eq = np.concatenate([eq_series[t][0]['n2_vol'] for t in rel_tags])
        n_eq = len(st_eq)
        joined = np.concatenate([st_eq, st])
        aj, visj, longj, firstj, ntrj = T.passage(joined, start_state, opp, 'I')
        conv = {
            'n_frames_stage4_5_6': [int(len(eq_series[t][1])) for t in rel_tags],
            'n_equilibration_frames_joined': int(n_eq),
            'expected_130_frames_1300_ps': bool(n_eq == 130),
            'time_of_joined_frame_k_from_release_ps': '(k + 1) * 10, k from 0',
            'states_in_stage4_5_6': {s: int(np.sum(st_eq == s)) for s in 'CRI'},
            'first_frame_in_opposite_state_joined': (lambda ix: {'joined_frame': int(ix[0]), 't_from_release_ps': r((ix[0] + 1) * 10.0, 1)}
                                                     if len(ix) else None)(np.nonzero(joined == opp)[0]),
            'last_frame_in_start_state_joined': (lambda ix: {'joined_frame': int(ix[-1]), 't_from_release_ps': r((ix[-1] + 1) * 10.0, 1)}
                                                 if len(ix) else None)(np.nonzero(joined == start_state)[0]),
            'n_long_visits_joined': len(longj),
            'long_visits_joined': [{'state': x[0], 'first_joined_frame': int(x[1]), 'n_frames': int(x[2]),
                                    't_first_from_release_ps': r((x[1] + 1) * 10.0, 1)} for x in longj],
            'n_transitions_joined': int(ntrj),
            'production_frames_in_start_state': int(np.sum(st == start_state)),
        }
        if firstj is None:
            conv.update({'has_conversion_time': False, 'where': None, 't_from_release_ps': None,
                         't_from_release_ns': None, 't_in_production_ns': None, 'report': 'no long visit in the opposite state: no time'})
        else:
            k = firstj[1]
            t_rel = (k + 1) * 10.0
            if k < n_eq:
                conv.update({'has_conversion_time': True, 'where': 'equilibration',
                             't_from_release_ps': r(t_rel, 1), 't_from_release_ns': r(t_rel / 1000.0, 3),
                             't_in_production_ns': None,
                             'stage_of_that_frame': ('eq04' if k < len(eq_series['eq04'][1]) else
                                                     'eq05' if k < len(eq_series['eq04'][1]) + len(eq_series['eq05'][1]) else 'eq06'),
                             'report': 'converted during the equilibration, %d ps after release' % round(t_rel)})
            else:
                tp = (k - n_eq + 1) * 0.01
                conv.update({'has_conversion_time': True, 'where': 'production',
                             't_from_release_ps': r(t_rel, 1), 't_from_release_ns': r(t_rel / 1000.0, 3),
                             't_in_production_ns': r(tp, 2),
                             'report': 'first long visit in the opposite state begins %.2f ns into the production, %.2f ns after release'
                                       % (tp, t_rel / 1000.0)})
        res['conversion'] = conv

    # ---------------- ligand RMSD
    lab = P['lig_rmsd_labelling']
    res['lig_rmsd_A'] = {
        'mean_whole': r(P['lig_rmsd'].mean()), 'mean_last_half': r(P['lig_rmsd'][half0:].mean()),
        'mean_first_half': r(P['lig_rmsd'][:half0].mean()), 'mean_last_20ns': r(P['lig_rmsd'][N - w:].mean()),
        'sd_whole': r(P['lig_rmsd'].std()), 'min': r(P['lig_rmsd'].min()), 'max': r(P['lig_rmsd'].max()),
        'n_frames_last_half': int(N - half0),
        'identity_labelling_mean_whole': r(P['lig_rmsd_identity'].mean()),
        'identity_labelling_mean_last_half': r(P['lig_rmsd_identity'][half0:].mean()),
        'frames_by_best_labelling': {T.LABELLINGS[j][0]: int(np.sum(lab == j)) for j in range(4)},
        'frames_above_2A_whole': int(np.sum(P['lig_rmsd'] > 2.0)),
        'frames_above_2A_last_half': int(np.sum(P['lig_rmsd'][half0:] > 2.0)),
        'rule1_mean_last_half_above_2A': bool(P['lig_rmsd'][half0:].mean() > 2.0),
    }
    # ---------------- hydrogen bond
    res['hbond_N2H_T199OG1'] = {
        'count': int(P['hb'].sum()), 'n': int(N), 'fraction': r(P['hb'].mean()),
        'fraction_last_half': r(P['hb'][half0:].mean()),
        'mean_N2_OG1_A': r(P['n2_og1'].mean()), 'mean_angle_deg': r(P['ang_n2_h_og1'].mean(), 2),
        'fraction_distance_ok': r(np.mean(P['n2_og1'] <= 3.5)), 'fraction_angle_ok': r(np.mean(P['ang_n2_h_og1'] >= 135.0)),
        'by_state': {s: ({'n': int(np.sum(st == s)), 'fraction': r(P['hb'][st == s].mean())} if np.any(st == s) else None) for s in 'CRI'},
    }
    # ---------------- zinc site
    res['zinc'] = {
        'mean_Zn_N2_A': r(P['zn_n2'].mean()), 'mean_Zn_H94NE2_A': r(P['zn_h94'].mean()),
        'mean_Zn_H96NE2_A': r(P['zn_h96'].mean()), 'mean_Zn_H119ND1_A': r(P['zn_h119'].mean()),
        'sd_Zn_N2_A': r(P['zn_n2'].std()), 'sd_Zn_H94NE2_A': r(P['zn_h94'].std()),
        'sd_Zn_H96NE2_A': r(P['zn_h96'].std()), 'sd_Zn_H119ND1_A': r(P['zn_h119'].std()),
        'max_of_the_four_A': r(max(P[k].max() for k in ('zn_n2', 'zn_h94', 'zn_h96', 'zn_h119'))),
        'sulfonyl_O_within_2.6A': {'count': int(np.sum(P['zn_omin'] <= 2.6)), 'n': int(N),
                                   'fraction': r(np.mean(P['zn_omin'] <= 2.6), 5),
                                   'at_most_5_percent_integers': bool(20 * int(np.sum(P['zn_omin'] <= 2.6)) <= N)},
        'mean_Zn_Osulfonyl_min_A': r(P['zn_omin'].mean()), 'min_Zn_Osulfonyl_min_A': r(P['zn_omin'].min()),
        'mean_Zn_O1_A': r(P['zn_o1'].mean()), 'mean_Zn_O2_A': r(P['zn_o2'].mean()),
    }
    # ---------------- configuration at N2
    vol = P['n2_vol']
    sg = np.sign(vol).astype(int)
    sa, svis, slong, sfirst, sntr = T.passage(list(sg), start_sign, -start_sign, 0)
    n2 = {
        'sign_of_start': start_sign,
        'frames_plus': int(np.sum(sg > 0)), 'frames_minus': int(np.sum(sg < 0)), 'frames_zero': int(np.sum(sg == 0)),
        'fraction_plus': r(np.mean(sg > 0), 5), 'fraction_minus': r(np.mean(sg < 0), 5),
        'mean_volume': r(vol.mean()), 'mean_volume_last_half': r(vol[half0:].mean()),
        'min_volume': r(vol.min()), 'max_volume': r(vol.max()),
        'frames_abs_volume_below_0.1': int(np.sum(np.abs(vol) < 0.1)),
        'n_sign_visits': len(svis), 'n_long_sign_visits': len(slong),
        'long_sign_visits': [{'sign': int(x[0]), 'first_frame': int(x[1]), 'n_frames': int(x[2]), 't_first_ns': tprod(x[1])} for x in slong[:20]],
        'n_lasting_sign_changes': int(sntr),
        'first_frame_with_other_sign_production': (lambda ix: {'frame': int(ix[0]), 't_ns': tprod(ix[0])} if len(ix) else None)(np.nonzero(sg != start_sign)[0]),
        'first_lasting_change_production_ns': (tprod(sfirst[1]) if sfirst else None),
        'sign_by_state': {s: {'plus': int(np.sum((st == s) & (sg > 0))), 'minus': int(np.sum((st == s) & (sg < 0)))} for s in 'CRI'},
    }
    if not missing:
        sgj = np.concatenate([np.sign(vol_eq).astype(int), sg])
        sja, sjvis, sjlong, sjfirst, sjntr = T.passage(list(sgj), start_sign, -start_sign, 0)
        ixo = np.nonzero(sgj != start_sign)[0]
        n2['joined_release_series'] = {
            'sign_in_stage4_5_6': {'plus': int(np.sum(vol_eq > 0)), 'minus': int(np.sum(vol_eq < 0))},
            'first_frame_with_other_sign': ({'joined_frame': int(ixo[0]), 't_from_release_ps': r((ixo[0] + 1) * 10.0, 1)} if len(ixo) else None),
            'first_lasting_change': (None if sjfirst is None else {
                'joined_frame': int(sjfirst[1]), 't_from_release_ps': r((sjfirst[1] + 1) * 10.0, 1),
                'where': 'equilibration' if sjfirst[1] < n_eq else 'production',
                't_in_production_ns': (None if sjfirst[1] < n_eq else r((sjfirst[1] - n_eq + 1) * 0.01, 2)),
                'length_of_that_visit_frames': int(sjfirst[2])}),
            'n_lasting_sign_changes': int(sjntr),
        }
    # sign in the held steps before the release (eq01..eq03 and extensions)
    held = [t for t in eq_series if t not in rel_tags]
    n2['sign_in_held_steps_before_release'] = {t: {'plus': int(np.sum(eq_series[t][0]['n2_vol'] > 0)),
                                                    'minus': int(np.sum(eq_series[t][0]['n2_vol'] < 0))} for t in held}
    res['n2_configuration'] = n2

    # ---------------- series kept in the scratch folder
    np.savez_compressed(os.path.join(HERE, 'series_%s_%s.npz' % (zset, run)),
                        state=st, assigned=a_arr, **{k: P[k] for k in KEYS},
                        **{'eq_%s_%s' % (t, k): eq_series[t][0][k] for t in eq_series for k in ('d1', 'd2', 'n2_vol', 'lig_rmsd', 'hb')})
    return res


def main():
    ref = T.build_reference()
    only = sys.argv[1:]
    out = {}
    raw_fn = os.path.join(HERE, 'runs_raw.json')
    if os.path.exists(raw_fn) and only:
        out = json.load(open(raw_fn))
    for zset, run, start, role in T.RUN_LIST:
        key = zset + '/' + run
        if only and key not in only:
            continue
        t0 = time.time()
        out[key] = one_run(zset, run, start, role, ref)
        out[key]['wall_s'] = round(time.time() - t0, 1)
        v = out[key]
        print(key, v['frames']['n_production_frames'], v['verdict']['verdict'], v['states']['last_20ns'],
              'rmsd', v['lig_rmsd_A']['mean_whole'], v['lig_rmsd_A']['mean_last_half'], '%.0fs' % (time.time() - t0), flush=True)
        with open(raw_fn, 'w') as f:
            json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
