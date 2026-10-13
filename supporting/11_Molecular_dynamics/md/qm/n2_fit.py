"""Gate G-R2, step 3: the fit of the ladder of models on the FITTING set (md/N2.md, sections 2.4 and 6).

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/n2_fit.py [--out DIR]
        [--points-per-parameter 8]        (8 is the rule; another value marks the output as a test of the script)
    N2_FF=ff2 in the environment: the force-field copy of the declared second step (charges from two placements,
    md/qm/n2_charges.py); its output then goes to md/qm/n2/fit_ff2 unless --out is given.

It cannot read a test energy: its only door to the energies is n2_common.load_points('fit') (through
n2_common.fit_rows), which never opens the result file of a point flagged "test" in md/qm/n2/split.json.

DATA.  The 36 rigid points of gate G-R and the new points flagged "fit" whose single points are finished.
y = E_QM - E(rest): B3LYP/6-31G* + D3(BJ) minus the force field without the terms that may change
(md/qm/n2/points.json).  Equal weights, one constant.  The 25 points of part A of gate G-R are among the 36 old
points and are fitted like every other point; the number of part A that this program prints is therefore
IN-SAMPLE (the leave-one-out errors of the same points are printed beside it).

CANDIDATES.  The five models of n2_common.LADDER times ridge 0, 0.01, 0.1; each fitted with the bounds of the
plan as constraints (n2_common.fit_model) and judged by its leave-one-out errors (n2_common.leave_one_out).
ELIGIBLE: at least 8 fitting points per parameter and a converged solver.

ORDER OF PREFERENCE (md/N2.md, section 2.4, with the addition of the design check of 6 Oct 2026, finding 1, both
fixed before any energy of a new point was read by a fit).  For every candidate two numbers from its leave-one-out
errors e_i: (a) their root mean square; (b) the same errors read against ONE reference, as condition 1 of the gate
reads the test set: r = the fitting row of lowest E_QM, root mean square of e_i - e_r over the other rows.  The
ranking number is the LARGER of (a) and (b).  Smallest first; candidates within 5 % of the smallest are ordered
among themselves by fewer parameters, then by the larger ridge; the rest by the ranking number.

Writes <out>/candidates.json (every candidate: numbers, parameters with standard errors, mismatch by family
before and after, part A in-sample, the R-like configuration against the crystal-like one, and which conditions
of the preview of the gate it is expected to meet), <out>/fit_points.csv (every fitting point with its residual
and leave-one-out error for the first candidate), <out>/MODEL_<rank>_<name>.json for the eligible ones (the files
for `n2_gate.py preview --model ...`) and <out>/n2_terms_<rank>_<name>.frcmod.  <out> is md/qm/n2/fit.
"""
import argparse
import json
import pathlib

import numpy as np

import n2_common as nc

RIDGES = (0.0, 0.01, 0.1)
FAMILY_X = ['X', 'a2', 'a3', 'a4', 'a5', 'a6', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'e_p015', 'e_p030', 'e_m120',
            'e_p120', 'f_p08', 'f_m08']                                    # as md/qm/n2_gate.py, md/GR.md 6.2
FAMILY_RX = ['Rx', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6']
PRE = dict(loo_rms=1.0, one_reference_rms=1.0, one_reference_largest=2.5, part_A_rms=1.0)   # limits of the gate


def rms(v):
    v = np.asarray(list(v), float)
    return float(np.sqrt(np.mean(v ** 2))) if len(v) else None


def r4(x):
    return None if x is None else round(float(x), 4)


def part_a(rows, dev):
    """Part A of gate G-R, two families, each against its point of lowest E_QM; dev = MM - QM up to a constant."""
    by = {r['id']: k for k, r in enumerate(rows)}
    if not all(i in by for i in FAMILY_X + FAMILY_RX):
        return None
    out, d_all = {}, []
    for name, ids in (('X', FAMILY_X), ('Rx', FAMILY_RX)):
        k0 = min((by[i] for i in ids), key=lambda k: rows[k]['e_qm'])
        dd = {i: float(dev[by[i]] - dev[k0]) for i in ids if by[i] != k0}
        out[name] = dict(reference=rows[k0]['id'], deviations={k: r4(v) for k, v in dd.items()})
        d_all += list(dd.values())
    k1 = min((by[i] for i in FAMILY_X + FAMILY_RX), key=lambda k: rows[k]['e_qm'])
    single = [dev[by[i]] - dev[k1] for i in FAMILY_X + FAMILY_RX if by[i] != k1]
    return dict(rms=r4(rms(d_all)), n=len(d_all), largest=r4(max(d_all, key=abs)), families=out,
                rms_single_reference=r4(rms(single)), single_reference=rows[k1]['id'])


def by_family(rows, dev):
    """Old points by family, each against the point it was made from (X; Rx for g and c2); new points by region,
    kind and energy, about the mean of all rows.  dev = MM - QM up to a constant."""
    by = {r['id']: k for k, r in enumerate(rows)}
    out = dict(old={}, new={})
    for fam in ('a', 'd', 'g', 'e', 'f', 'c1', 'c2'):
        base = by.get('Rx' if fam in ('g', 'c2') else 'X')
        ks = [k for k, r in enumerate(rows) if r['family'] == fam]
        if ks and base is not None:
            d = np.array([dev[k] - dev[base] for k in ks])
            out['old'][fam] = dict(n=len(ks), rms=r4(rms(d)), mean=r4(d.mean()), largest=r4(d[np.argmax(np.abs(d))]))
    for a, b in (('Rx', 'X'), ('Nx', 'X')):
        if a in by and b in by:
            out['old'][f'{a}_against_{b}'] = r4(dev[by[a]] - dev[by[b]])
    c = dev - dev.mean()
    e_low = min(r['e_qm'] for r in rows)
    groups = dict(
        all_rows=list(range(len(rows))),
        old=[k for k, r in enumerate(rows) if r['set'] == 'old'],
        new=[k for k, r in enumerate(rows) if r['set'] != 'old'],
        new_X=[k for k, r in enumerate(rows) if r['set'] != 'old' and r['region'] == 'X'],
        new_between=[k for k, r in enumerate(rows) if r['set'] != 'old' and r['region'] == 'between'],
        new_R=[k for k, r in enumerate(rows) if r['set'] != 'old' and r['region'] == 'R'],
        new_design_core=[k for k, r in enumerate(rows) if r['kind'] == 'design' and r['layer'] == 'core'],
        new_design_margin=[k for k, r in enumerate(rows) if r['kind'] == 'design' and r['layer'] == 'margin'],
        new_frames=[k for k, r in enumerate(rows) if r['kind'] == 'frame'],
        within_10_of_the_lowest_E_QM=[k for k, r in enumerate(rows) if r['e_qm'] - e_low <= 10.0],
        above_10=[k for k, r in enumerate(rows) if r['e_qm'] - e_low > 10.0],
        not_flagged_protein=[k for k, r in enumerate(rows) if r['set'] != 'old' and not r['protein_contact']])
    for g, ks in groups.items():
        if ks:
            out['new'][g] = dict(n=len(ks), rms=r4(rms(c[ks])), mean=r4(c[ks].mean()),
                                 largest=r4(c[ks][np.argmax(np.abs(c[ks]))]))
    xr = [k for k, r in enumerate(rows) if r['set'] != 'old' and r['region'] == 'R']
    xx = [k for k, r in enumerate(rows) if r['set'] != 'old' and r['region'] == 'X']
    out['mean_of_R_like_minus_mean_of_X_like_new_points'] = r4(c[xr].mean() - c[xx].mean()) if xr and xx else None
    return out


def candidate(rows, name, groups, ridge, th_opt, per_parameter):
    f = nc.fit_model(rows, groups, ridge, True, th_opt)
    lo = nc.leave_one_out(rows, groups, ridge, True, th_opt)
    res, err = np.array(f['residuals']), np.array(lo['errors'])
    e_qm = np.array([r['e_qm'] for r in rows])
    one = nc.one_reference(err, e_qm)
    one_fit = nc.one_reference(res, e_qm)
    se = nc.standard_errors(f)
    # the AMBER numbers as the frcmod will carry them (three and four decimals), and their bounds
    amber = nc.amber_terms(f['parameters'])
    rounded = json.loads(json.dumps(amber))
    for v in rounded['angles'].values():
        v['K'], v['theta0'] = round(v['K'], 3), round(v['theta0'], 3)
    for terms in rounded['torsions'].values():
        for t in terms:
            t['PK'] = round(t['PK'], 4)
    for v in (rounded['improper'] or {}).values():
        v['PK'] = round(v['PK'], 4)
    outside = nc.within_bounds(rounded, th_opt)
    pa, pa_loo = part_a(rows, res), part_a(rows, err)
    eligible = bool(len(rows) >= per_parameter * f['n_parameters'] and f['solver']['success'])
    pre = dict(
        c3_part_A_in_sample=None if pa is None else bool(pa['rms'] <= PRE['part_A_rms']),
        c4_inside_the_bounds=not outside,
        c5_leave_one_out=bool(lo['rms'] <= PRE['loo_rms'] and eligible),
        pre_test_one_reference=bool(one['rms'] <= PRE['one_reference_rms'] and abs(one['largest']) <= PRE['one_reference_largest']),
        c6_minimum='not computed here: md/qm/n2_behaviour.py')
    by = {r['id']: k for k, r in enumerate(rows)}
    return dict(
        name=f'{name}_ridge{ridge:g}', ladder=name, groups=list(groups), ridge=ridge, n_parameters=f['n_parameters'],
        n_points=len(rows), points_per_parameter=round(len(rows) / f['n_parameters'], 2), eligible=eligible,
        residual_rms=r4(f['rms']), residual_largest=r4(res[np.argmax(np.abs(res))]),
        residual_largest_at=rows[int(np.argmax(np.abs(res)))]['id'],
        leave_one_out_rms=r4(lo['rms']), leave_one_out_largest=r4(err[np.argmax(np.abs(err))]),
        leave_one_out_largest_at=rows[int(np.argmax(np.abs(err)))]['id'],
        one_reference=dict(reference=rows[one['reference_row']]['id'], leave_one_out_rms=r4(one['rms']),
                           leave_one_out_largest=r4(one['largest']), leave_one_out_mean=r4(one['mean']),
                           leave_one_out_error_of_the_reference=r4(one['error_of_the_reference']),
                           residual_rms=r4(one_fit['rms']), residual_largest=r4(one_fit['largest'])),
        ranking_number=r4(max(lo['rms'], one['rms'])),
        expected_in_the_preview_of_the_gate=pre,
        expected_to_pass_the_preview_without_condition_6=bool(all(v is True for k, v in pre.items() if k != 'c6_minimum')),
        part_A_in_sample=pa, part_A_from_leave_one_out_errors=pa_loo,
        by_family=by_family(rows, res), by_family_leave_one_out=by_family(rows, err),
        Rx_against_X=(dict(QM=r4(e_qm[by['Rx']] - e_qm[by['X']]),
                           force_field_as_built=r4(rows[by['Rx']]['mm_total'] - rows[by['X']]['mm_total']),
                           fitted=r4(e_qm[by['Rx']] - e_qm[by['X']] + res[by['Rx']] - res[by['X']]),
                           fitted_minus_QM=r4(res[by['Rx']] - res[by['X']]),
                           leave_one_out_minus_QM=r4(err[by['Rx']] - err[by['X']])) if 'Rx' in by and 'X' in by else None),
        solver=f['solver'], parameters=f['parameters'], amber=amber, amber_as_written=rounded,
        standard_errors=se, outside_bounds=outside,
        _residuals=[float(v) for v in res], _loo=[float(v) for v in err])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=None)
    ap.add_argument('--points-per-parameter', type=float, default=8.0)
    ap.add_argument('--only', default=None, help='comma-separated ladder names (a test of the script)')
    ap.add_argument('--sensitivity', choices=['hydrogen-contacts', 'visited-ranges', 'both'], default=None,
                    help='NOT the fit of the rule; subsets chosen by GEOMETRY alone. hydrogen-contacts: without the '
                         'fitting points that have a contact with a hydrogen shorter than the limits gate G-R used for '
                         'its turns (heavy atom to hydrogen 2.1 A, hydrogen to hydrogen 1.9 A). visited-ranges: without '
                         'the designed points of the widened ranges (layer "margin"). both: without either. Writes no '
                         'model file; --out is required')
    args = ap.parse_args()
    out = pathlib.Path(args.out) if args.out else nc.N2 / ('fit' if nc.FF_NAME == 'ff' else f'fit_{nc.FF_NAME}')
    out.mkdir(parents=True, exist_ok=True)
    rule = args.points_per_parameter == 8.0 and not args.only and not args.sensitivity
    if args.sensitivity and not args.out:
        ap.error('--sensitivity needs --out (it must not be written where the fit of the rule is)')
    _, th_opt = nc.qm_optimum()
    ff_what = 'as built' if nc.FF_NAME == 'ff' else json.loads((nc.FFDIR / 'mm_points.json').read_text())['what']
    rows = nc.fit_rows('fit')
    assert not any(r['set'] == 'test' for r in rows)
    left_out = []
    if args.sensitivity:
        pj = json.loads((nc.N2 / 'points.json').read_text())
        con = {p['id']: p['contacts']['model'] for p in pj['points']}
        con.update({p['id']: p['contacts']['model'] for p in json.loads((nc.GR / 'points.json').read_text())['points']})
        if args.sensitivity in ('hydrogen-contacts', 'both'):
            left_out += [r['id'] for r in rows if con[r['id']]['heavy_h']['d'] < 2.1 or con[r['id']]['h_h']['d'] < 1.9]
        if args.sensitivity in ('visited-ranges', 'both'):
            left_out += [r['id'] for r in rows if r['layer'] == 'margin' and r['id'] not in left_out]
        rows = [r for r in rows if r['id'] not in left_out]
    d0 = np.array([r['mm_total'] - r['e_qm'] for r in rows])
    cands = []
    for name, groups in nc.LADDER:
        if args.only and name not in args.only.split(','):
            continue
        for ridge in RIDGES:
            cands.append(candidate(rows, name, groups, ridge, th_opt, args.points_per_parameter))
            c = cands[-1]
            print(f"  fitted {c['name']:14s} rms {c['residual_rms']:.3f}  LOO {c['leave_one_out_rms']:.3f}  "
                  f"one reference {c['one_reference']['leave_one_out_rms']:.3f}", flush=True)
    el = [c for c in cands if c['eligible']]
    order = []
    if el:
        best = min(c['ranking_number'] for c in el)
        near = sorted([c for c in el if c['ranking_number'] <= 1.05 * best],
                      key=lambda c: (c['n_parameters'], -c['ridge'], c['ranking_number']))
        far = sorted([c for c in el if c['ranking_number'] > 1.05 * best], key=lambda c: c['ranking_number'])
        order = near + far
    for old in list(out.glob('MODEL_*.json')) + list(out.glob('n2_terms_*.frcmod')):
        old.unlink()
    mark = '' if rule else (f' (SENSITIVITY "{args.sensitivity}": fitting points left out by their geometry; NOT the fit of the rule)'
                            if args.sensitivity else ' (TEST OF THE SCRIPT, not the rule of 8 points per parameter)')
    for k, c in enumerate(order, 1):
        c['rank'] = k
        if args.sensitivity:
            continue                                                       # no model file: it cannot go to the gate
        nc.write_json(out / f"MODEL_{k:02d}_{c['name']}.json", dict(
            name=c['name'] + mark, groups=c['groups'], ridge=c['ridge'], parameters=c['parameters'], rank=k,
            leave_one_out_rms=c['leave_one_out_rms'], ranking_number=c['ranking_number'], n_points=c['n_points'],
            force_field=nc.FF_NAME, force_field_what=ff_what, written=nc.now()))
        nc.write_frcmod(c['parameters'], out / f"n2_terms_{k:02d}_{c['name']}.frcmod",
                        f"terms at N2, candidate {c['name']}{mark} of gate G-R2 (md/N2.md); written by md/qm/n2_fit.py; "
                        f"NOT ACCEPTED by the gate")
    first = order[0] if order else None
    lines = ['id,set,region,family,kind,layer,E_QM_above_the_lowest,MM_as_built_minus_QM_about_the_mean,'
             'residual_first_candidate,leave_one_out_error_first_candidate']
    e_low = min(r['e_qm'] for r in rows)
    for k, r in enumerate(rows):
        lines.append(','.join(str(v) for v in (
            r['id'], r['set'], r['region'], r['family'] or '', r['kind'], r['layer'] or '', round(r['e_qm'] - e_low, 4),
            round(float(d0[k] - d0.mean()), 4), round(first['_residuals'][k], 4) if first else '',
            round(first['_loo'][k], 4) if first else '')))
    (out / 'fit_points.csv').write_text('\n'.join(lines) + '\n')
    for c in cands:
        c.pop('_residuals'), c.pop('_loo')
    doc = dict(what='fit of the ladder of gate G-R2 on the fitting set (md/N2.md, sections 2.4 and 6); written by md/qm/n2_fit.py'
                    + ('' if rule else '; ' + mark.strip(' ()')),
               sensitivity=args.sensitivity, fitting_points_left_out=left_out,
               written=nc.now(), force_field=nc.FF_NAME, force_field_what=ff_what,
               points=dict(old=sum(r['set'] == 'old' for r in rows), fit=sum(r['set'] == 'fit' for r in rows),
                           test_points_read=0, ids=[r['id'] for r in rows]),
               lowest_E_QM_of_the_fitting_set=min(rows, key=lambda r: r['e_qm'])['id'],
               points_per_parameter_required=args.points_per_parameter,
               as_built=dict(rms_about_the_mean=r4(rms(d0 - d0.mean())), by_family=by_family(rows, d0), part_A=part_a(rows, d0),
                             one_reference_rms=r4(nc.one_reference(d0, [r['e_qm'] for r in rows])['rms'])),
               order_of_preference=[c['name'] for c in order],
               rule='ranking number = the larger of the leave-one-out r.m.s. and the leave-one-out r.m.s. against one '
                    'reference (the fitting row of lowest E_QM); smallest first; within 5 % of the smallest: fewer '
                    'parameters, then the larger ridge',
               limits_of_the_gate_used_for_the_expectation=PRE,
               any_candidate_expected_to_pass=bool(any(c['expected_to_pass_the_preview_without_condition_6'] for c in el)),
               candidates=sorted(cands, key=lambda c: c.get('rank', 999)),
               split_sha256=nc.sha(nc.N2 / 'split.json'), points_sha256=nc.sha(nc.N2 / 'points.json'),
               n2_common_sha256=nc.sha(nc.HERE / 'n2_common.py'), versions=dict(numpy=np.__version__))
    nc.write_json(out / 'candidates.json', doc)
    print(doc['what'])
    print('points:', {k: v for k, v in doc['points'].items() if k != 'ids'}, '| as built, r.m.s. about the mean:',
          doc['as_built']['rms_about_the_mean'], '| lowest E_QM:', doc['lowest_E_QM_of_the_fitting_set'])
    for c in doc['candidates']:
        pa = c['part_A_in_sample'] or {}
        print(f"{c.get('rank', '-'):>3} {c['name']:14s} P {c['n_parameters']:2d} ({c['points_per_parameter']:5.1f}/par) "
              f"{'eligible' if c['eligible'] else 'not eligible':12s} rms {c['residual_rms']:6.3f} LOO {c['leave_one_out_rms']:6.3f} "
              f"one-ref {c['one_reference']['leave_one_out_rms']:6.3f} (largest {c['one_reference']['leave_one_out_largest']:+.2f}) "
              f"part A {pa.get('rms')} Rx-X {c['Rx_against_X']['fitted_minus_QM'] if c['Rx_against_X'] else None} "
              f"at a bound {len(c['standard_errors']['at_a_bound'])} expected to pass {c['expected_to_pass_the_preview_without_condition_6']}")
    print('order of preference:', doc['order_of_preference'] or 'no eligible candidate')


if __name__ == '__main__':
    main()
