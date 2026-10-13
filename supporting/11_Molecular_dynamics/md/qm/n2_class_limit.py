"""Gate G-R2, a DIAGNOSIS on the fitting set: how far could the terms that the plan allows go at all?

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/n2_class_limit.py

NOT a candidate, NOT the rule of the design, and nothing here can go to the gate (no model file is written).
The ladder of md/N2.md, section 2.4, was fixed before any energy and restricts the plan's point (2): cosines only
about N2-S, and about Zn-N2 one threefold term common to the three histidines.  After the ladder failed on the
fitting set (section 6.2) the question for the user's decision is whether that restriction is what stops it.
This script answers with plain least squares WITHOUT bounds and without ridge (so every number is a LOWER limit
of the residual that a bounded fit of the same terms could reach), on the fitting rows only
(n2_common.fit_rows('fit'): no test result is opened), for both charge sets:

    L5            the largest model of the ladder (17 parameters): must return the numbers of n2_fit.py, L5_ridge0
    L5+sin        + sin n*phi, n = 1..3, on H-N2-S-C4 and on Zn-N2-S-C4 (what the O-ended torsions or free phases
                  would add)                                                              23 parameters
    all allowed   the angles, cos and sin n = 1..3 on both N2-S series, the improper, and cos and sin n = 1..3 on
                  N(His)-Zn-N2-S and N(His)-Zn-N2-H SEPARATELY for each histidine       55 parameters
                  (3.2 rows per parameter: far below the 8 the design requires; shown as the limit of the class)

Reported for each: residual r.m.s., the exact leave-one-out r.m.s. (hat matrix), rows per parameter; the same
without the fitting points with a hydrogen contact below the limits of gate G-R (the sensitivity of section 6.4).
Writes md/qm/n2/class_limit.json.
"""
import importlib
import json
import os

import numpy as np


def design(rows, kind, nc):
    names = [k for k, _ in nc.FEATURES]
    cols, labels = [], []
    common = ('His_ZnN_S_3c', 'His_ZnN_S_3s', 'His_ZnN_H_3c', 'His_ZnN_H_3s')
    for k in names:
        if kind == 'all allowed' and k in common:
            continue                                                       # replaced by the terms per histidine
        cols.append([r['features'][k] for r in rows])
        labels.append(k)
    rad = np.radians
    if kind in ('L5+sin', 'all allowed'):
        for key, lab in (('H_N2_S_C4', 'H_NS_C4'), ('tau', 'Zn_NS_C4')):
            for n in (1, 2, 3):
                cols.append([np.sin(n * rad(r['coords'][key])) for r in rows])
                labels.append(f'{lab}_{n}s')
    if kind == 'all allowed':
        for h in ('N94', 'N96', 'N119'):
            for end in ('S', 'H'):
                for n in (1, 2, 3):
                    for fn, tag in ((np.cos, 'c'), (np.sin, 's')):
                        cols.append([fn(n * rad(r['coords'][f'{h}_Zn_N2_{end}'])) for r in rows])
                        labels.append(f'{h}_ZnN_{end}_{n}{tag}')
    cols.append([1.0] * len(rows))
    labels.append('const')
    return np.array(cols).T, labels


def fit(rows, kind, nc):
    x, labels = design(rows, kind, nc)
    y = np.array([r['y'] for r in rows])
    pinv = np.linalg.pinv(x)
    res = y - x @ (pinv @ y)
    lev = np.clip(np.einsum('ij,ji->i', x, pinv), 0.0, 1.0 - 1e-9)
    loo = res / (1.0 - lev)
    rank = int(np.linalg.matrix_rank(x))
    npar = len(labels) - 1
    return dict(parameters=npar, rank_with_the_constant=rank, rows=len(rows), rows_per_parameter=round(len(rows) / npar, 2),
                residual_rms=round(float(np.sqrt(np.mean(res ** 2))), 4),
                leave_one_out_rms=round(float(np.sqrt(np.mean(loo ** 2))), 4),
                largest_leverage=round(float(lev.max()), 3),
                residual_largest=round(float(res[np.argmax(np.abs(res))]), 3),
                residual_largest_at=rows[int(np.argmax(np.abs(res)))]['id'])


def main():
    out = dict(what='DIAGNOSIS, not a candidate and not the rule: unbounded least squares of three sets of terms '
                    'that the plan allows, on the fitting rows only; written by md/qm/n2_class_limit.py',
               note='no bounds, no ridge: the residual is a lower limit for a bounded fit of the same terms; the '
                    'leave-one-out error is exact for the unbounded fit (hat matrix)', charge_sets={})
    for ffn in ('ff', 'ff2'):
        os.environ['N2_FF'] = ffn
        import n2_common as nc
        importlib.reload(nc)
        rows = nc.fit_rows('fit')
        assert not any(r['set'] == 'test' for r in rows)
        pj = json.loads((nc.N2 / 'points.json').read_text())
        con = {p['id']: p['contacts']['model'] for p in pj['points']}
        con.update({p['id']: p['contacts']['model'] for p in json.loads((nc.GR / 'points.json').read_text())['points']})
        short = [r['id'] for r in rows if con[r['id']]['heavy_h']['d'] < 2.1 or con[r['id']]['h_h']['d'] < 1.9]
        sub = [r for r in rows if r['id'] not in short]
        rec = dict(all_fitting_rows={}, without_the_short_hydrogen_contacts={}, fitting_points_left_out=short)
        for kind in ('L5', 'L5+sin', 'all allowed'):
            rec['all_fitting_rows'][kind] = fit(rows, kind, nc)
            rec['without_the_short_hydrogen_contacts'][kind] = fit(sub, kind, nc)
        prog = json.loads((nc.N2 / ('fit' if ffn == 'ff' else f'fit_{ffn}') / 'candidates.json').read_text())
        c0 = next(c for c in prog['candidates'] if c['name'] == 'L5_ridge0')
        rec['check_against_n2_fit_L5_ridge0'] = dict(
            residual_rms=c0['residual_rms'], leave_one_out_rms=c0['leave_one_out_rms'],
            agrees=bool(abs(c0['residual_rms'] - rec['all_fitting_rows']['L5']['residual_rms']) < 2e-4
                        and abs(c0['leave_one_out_rms'] - rec['all_fitting_rows']['L5']['leave_one_out_rms']) < 2e-3))
        out['charge_sets'][ffn] = rec
        for scope in ('all_fitting_rows', 'without_the_short_hydrogen_contacts'):
            for kind, v in rec[scope].items():
                print(f"{ffn:4s} {scope:36s} {kind:12s} P {v['parameters']:2d} rows {v['rows']:3d} ({v['rows_per_parameter']:5.2f}/par) "
                      f"residual {v['residual_rms']:.3f}  leave-one-out {v['leave_one_out_rms']:.3f}")
        print('     check against n2_fit.py:', rec['check_against_n2_fit_L5_ridge0'])
    import n2_common as nc
    out['written'] = nc.now()
    out['split_sha256'] = nc.sha(nc.N2 / 'split.json')
    nc.write_json(nc.N2 / 'class_limit.json', out)


if __name__ == '__main__':
    main()
