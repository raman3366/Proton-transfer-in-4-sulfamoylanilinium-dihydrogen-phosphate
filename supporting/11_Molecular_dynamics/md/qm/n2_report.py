"""Gate G-R2: the tables of md/N2.md, sections 5 and 6, printed as Markdown from the files of the fit.

    /root/cqc-sfa-tools/md/bin/python md/qm/n2_report.py [--fit md/qm/n2/fit] [--rank 1]

Reads <fit>/candidates.json, <fit>/fit_points.csv and md/qm/n2/points.json; computes nothing new except sums over
the fitting points (root mean squares by the closest contact of a point).  It reads no result file of a point.
"""
import argparse
import csv
import json
import pathlib

import numpy as np

import n2_common as nc


def f2(x, d=2, sign=False):
    if x is None:
        return '–'
    s = f'{x:+.{d}f}' if sign else f'{x:.{d}f}'
    return s.replace('-', '−')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fit', default=str(nc.N2 / 'fit'))
    ap.add_argument('--rank', type=int, default=1)
    args = ap.parse_args()
    fit = pathlib.Path(args.fit)
    doc = json.loads((fit / 'candidates.json').read_text())
    cands = doc['candidates']
    print(f"<!-- {doc['what']}; written {doc['written']}; force field {doc.get('force_field')} -->")
    print(f"Fitting rows: {doc['points']['old']} old + {doc['points']['fit']} new = {doc['points']['old'] + doc['points']['fit']}; "
          f"lowest E_QM of the fitting set: {doc['lowest_E_QM_of_the_fitting_set']}; as built, r.m.s. about the mean "
          f"{doc['as_built']['rms_about_the_mean']}, against one reference {doc['as_built']['one_reference_rms']}; part A as built "
          f"{doc['as_built']['part_A']['rms']} ({doc['as_built']['part_A']['rms_single_reference']} single reference)\n")
    print('| Rank | Candidate | Free parameters | Points per parameter | Residual r.m.s. | Leave-one-out r.m.s. | Leave-one-out against one '
          'reference: r.m.s., largest | Ranking number | Part A in-sample (single reference) | Rx against X, fitted − QM | At a bound | '
          'Expected to pass the preview (without condition 6) |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|')
    for c in cands:
        pa = c['part_A_in_sample'] or {}
        o = c['one_reference']
        print(f"| {c.get('rank', '–')} | {c['name']} | {c['n_parameters']} | {c['points_per_parameter']} | {f2(c['residual_rms'])} | "
              f"{f2(c['leave_one_out_rms'])} | {f2(o['leave_one_out_rms'])}, {f2(o['leave_one_out_largest'], 2, True)} | "
              f"{f2(c['ranking_number'])} | {f2(pa.get('rms'))} ({f2(pa.get('rms_single_reference'))}) | "
              f"{f2(c['Rx_against_X']['fitted_minus_QM'], 2, True)} | {len(c['standard_errors']['at_a_bound'])} | "
              f"{'yes' if c['expected_to_pass_the_preview_without_condition_6'] else 'no'}{'' if c['eligible'] else ' (not eligible)'} |")
    c = [x for x in cands if x.get('rank') == args.rank][0]
    se, am = c['standard_errors'], c['amber_as_written']
    print(f"\n**Candidate of rank {args.rank}: {c['name']}** ({c['n_parameters']} free parameters, {c['n_points']} fitting rows; sigma "
          f"{f2(se['sigma'])}; at a bound: {se['at_a_bound'] or 'none'})\n")
    print('| Term | As built | Fitted | Standard error | Window or bound |')
    print('|---|---|---|---|---|')
    built = {'ZnN2S': ('Zn–N2–S', 'M1-Y4-sy', 50.0, 117.55), 'ZnN2H': ('Zn–N2–H', 'M1-Y4-hn', 50.0, 107.41), 'SN2H': ('S–N2–H', 'hn-Y4-sy', 61.55, 118.06)}
    _, opt = nc.qm_optimum()
    for a, (lab, typ, k0, t0) in built.items():
        v = se['angles'][a]
        print(f"| {lab}, K (kcal mol⁻¹ rad⁻²) | {k0} | {f2(am['angles'][typ]['K'], 3)} | {f2(v['K_se'])} | 10 to 300 |")
        print(f"| {lab}, θ₀ (°) | {t0} | {f2(am['angles'][typ]['theta0'], 3)} | {f2(v['theta0_se'])} | {f2(opt[a] - 15)} to {f2(opt[a] + 15)} "
              f"(optimum {f2(opt[a])}) |")
    lab = dict(H_NS_C4='H–N2–S–C4', Zn_NS_C4='Zn–N2–S–C4')
    built_t = {'H_NS_C4_1': -0.2624, 'H_NS_C4_2': 2.703}
    for k, v in se['coefficients'].items():
        if k == 'improper':
            print(f"| improper at N2, c on 1 − cos 2ω | none | {f2(v['value'], 4, True)} | {f2(v['se'], 3)} | \\|c\\| ≤ 5 |")
        else:
            stem, n = k.rsplit('_', 1)
            print(f"| {lab[stem]}, coefficient of cos {n}φ | {f2(built_t.get(k, 0.0), 4, True)} | {f2(v['value'], 4, True)} | {f2(v['se'], 3)} | "
                  f"\\|c\\| ≤ 5 |")
    for k, v in se['pairs'].items():
        end = 'S' if k.endswith('_S') else 'H'
        print(f"| N(His)–Zn–N2–{end}, threefold: amplitude; phase (°) | 0 | {f2(v['PK'], 4)}; {f2(v['phase'], 1)} | "
              f"{f2(v['PK_se'], 3)}; {f2(v['phase_se'], 1)} | amplitude ≤ 5 |")
    print(f"\nResidual r.m.s. {f2(c['residual_rms'])} (largest {f2(c['residual_largest'], 2, True)} at {c['residual_largest_at']}); "
          f"leave-one-out r.m.s. {f2(c['leave_one_out_rms'])} (largest {f2(c['leave_one_out_largest'], 2, True)} at "
          f"{c['leave_one_out_largest_at']}); against one reference ({c['one_reference']['reference']}): "
          f"{f2(c['one_reference']['leave_one_out_rms'])}, largest {f2(c['one_reference']['leave_one_out_largest'], 2, True)}.\n")
    b0, b1, b2 = doc['as_built']['by_family'], c['by_family'], c['by_family_leave_one_out']
    print('| Family | n | As built: r.m.s. (mean) | Fitted: r.m.s. (mean) | Leave-one-out: r.m.s. (mean) |')
    print('|---|---|---|---|---|')
    names = dict(a='a, inversion path (against X)', d='d, hydrogen moved on X', g='g, hydrogen moved on Rx (against Rx)',
                 e='e, turns about Zn–N2', f='f, bend of Zn–N2–S', c1='c1, turns about N2–S from X', c2='c2, turns about N2–S from Rx (against Rx)')
    for k, lab_ in names.items():
        if k in b0['old']:
            print(f"| {lab_} | {b0['old'][k]['n']} | {f2(b0['old'][k]['rms'])} ({f2(b0['old'][k]['mean'], 2, True)}) | "
                  f"{f2(b1['old'][k]['rms'])} ({f2(b1['old'][k]['mean'], 2, True)}) | {f2(b2['old'][k]['rms'])} ({f2(b2['old'][k]['mean'], 2, True)}) |")
    for k in ('Rx_against_X', 'Nx_against_X'):
        print(f"| {k.replace('_', ' ')} | 1 | {f2(b0['old'][k], 2, True)} | {f2(b1['old'][k], 2, True)} | {f2(b2['old'][k], 2, True)} |")
    names = dict(all_rows='all fitting rows, about their mean', old='the 36 old points', new='the new fitting points', new_X='new, X-like',
                 new_between='new, between', new_R='new, R-like', new_design_core='new, designed, visited ranges',
                 new_design_margin='new, designed, widened ranges', new_frames='new, frames of the rehearsal',
                 within_10_of_the_lowest_E_QM='E_QM within 10 of the lowest', above_10='E_QM more than 10 above the lowest',
                 not_flagged_protein='new, not flagged "protein"')
    for k, lab_ in names.items():
        if k in b0['new']:
            print(f"| {lab_} | {b0['new'][k]['n']} | {f2(b0['new'][k]['rms'])} ({f2(b0['new'][k]['mean'], 2, True)}) | "
                  f"{f2(b1['new'][k]['rms'])} ({f2(b1['new'][k]['mean'], 2, True)}) | {f2(b2['new'][k]['rms'])} ({f2(b2['new'][k]['mean'], 2, True)}) |")
    k = 'mean_of_R_like_minus_mean_of_X_like_new_points'
    print(f"| mean of the R-like new points minus mean of the X-like | | {f2(b0[k], 2, True)} | {f2(b1[k], 2, True)} | {f2(b2[k], 2, True)} |")
    r = c['Rx_against_X']
    print(f"\nRx against X: quantum chemistry {f2(r['QM'], 2, True)}, force field as built {f2(r['force_field_as_built'], 2, True)}, "
          f"fitted {f2(r['fitted'], 2, True)} (fitted − QM {f2(r['fitted_minus_QM'], 2, True)}; from the leave-one-out errors "
          f"{f2(r['leave_one_out_minus_QM'], 2, True)}).")
    pa, pl = c['part_A_in_sample'], c['part_A_from_leave_one_out_errors']
    if pa:
        print(f"Part A (23 comparisons, two families; references {pa['families']['X']['reference']} and {pa['families']['Rx']['reference']}): "
              f"as built {doc['as_built']['part_A']['rms']}, fitted in-sample {pa['rms']} (single reference {pa['rms_single_reference']}, largest "
              f"{pa['largest']}); from the leave-one-out errors {pl['rms']} ({pl['rms_single_reference']}).")
    print('Expected in the preview of the gate:', c['expected_in_the_preview_of_the_gate'])
    # residuals by the closest contacts of the point (fitting points only)
    pj = json.loads((nc.N2 / 'points.json').read_text())
    con = {p['id']: p['contacts']['model'] for p in pj['points']}
    old = {p['id']: p['contacts']['model'] for p in json.loads((nc.GR / 'points.json').read_text())['points']}
    con.update(old)
    rows = list(csv.DictReader(open(fit / 'fit_points.csv')))
    if args.rank == 1:
        print('\n| Fitting points by their closest heavy atom to hydrogen contact (ligand to imidazoles) | n | As built about the mean, r.m.s. | '
              'Residual r.m.s. | Leave-one-out r.m.s. | Largest leave-one-out error |')
        print('|---|---|---|---|---|---|')
        bins = [(0, 2.1, 'below 2.1 Å (the limit gate G-R used for its own turns)'), (2.1, 2.3, '2.1 to 2.3 Å'), (2.3, 2.6, '2.3 to 2.6 Å'), (2.6, 99, 'above 2.6 Å')]
        for lo, hi, lab_ in bins:
            sel = [r_ for r_ in rows if lo <= con[r_['id']]['heavy_h']['d'] < hi]
            if not sel:
                continue
            a = np.array([float(r_['MM_as_built_minus_QM_about_the_mean']) for r_ in sel])
            b = np.array([float(r_['residual_first_candidate']) for r_ in sel])
            l = np.array([float(r_['leave_one_out_error_first_candidate']) for r_ in sel])
            print(f"| {lab_} | {len(sel)} | {f2(float(np.sqrt((a ** 2).mean())))} | {f2(float(np.sqrt((b ** 2).mean())))} | "
                  f"{f2(float(np.sqrt((l ** 2).mean())))} | {f2(float(l[np.argmax(np.abs(l))]), 2, True)} |")
        worst = sorted(rows, key=lambda r_: -abs(float(r_['leave_one_out_error_first_candidate'])))[:8]
        print('\nLargest leave-one-out errors: ' + '; '.join(
            f"{r_['id']} {f2(float(r_['leave_one_out_error_first_candidate']), 2, True)} ({r_['region']}, heavy–H "
            f"{con[r_['id']]['heavy_h']['d']} {con[r_['id']]['heavy_h']['pair']}, H–H {con[r_['id']]['h_h']['d']})" for r_ in worst))


if __name__ == '__main__':
    main()
