#!/usr/bin/env python
"""Gate G-R2: the examiner's numbers (conditions_fitting_data.json) beside those of the designer's own gate script,
run by the examiner in PREVIEW mode (md/qm/n2_gate.py preview; fitting data only; outputs in
md/qm/n2/final_test/designer_preview_<candidate>/), and one sensitivity that the builder reported, computed again
with the examiner's code (NOT the rule, NOT a candidate: the limits were chosen after fitting residuals were seen).

    nice -n 10 /root/cqc-sfa-tools/md/bin/python md/qm/n2/final_test/n2_exam_compare.py

Writes md/qm/n2/final_test/comparison_with_the_designer.json.  Fitting data only.
"""
import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import n2_exam_conditions as ex                                            # noqa: E402


def main():
    mine = json.loads((ex.OUT / 'conditions_fitting_data.json').read_text())
    out = dict(what="the examiner's numbers beside the designer's gate script in preview mode (fitting data only)", written=ex.now(), candidates={})
    for name in ex.CAND:
        g = json.loads((ex.OUT / f'designer_preview_{name}' / 'n2_gate.json').read_text())
        m = mine['candidates'][name]
        c = g['conditions']
        rows = [
            ('condition 3, part A r.m.s.', m['condition_3_part_A']['rms'], c['c3_part_A']['rms']),
            ('condition 3, single reference', m['condition_3_part_A']['single_reference']['rms'], c['c3_part_A']['not_deciding']['rms_single_reference']),
            ('condition 3, largest', m['condition_3_part_A']['largest'], c['c3_part_A']['largest']),
            ('condition 5, leave-one-out r.m.s.', m['condition_5_leave_one_out']['leave_one_out_rms'], c['c5_leave_one_out']['leave_one_out_rms']),
            ('condition 5, residual r.m.s.', m['condition_5_leave_one_out']['residual_rms'], c['c5_leave_one_out']['residual_rms_of_the_fit']),
            ('pre-test, leave-one-out against one reference', m['condition_5_leave_one_out']['against_one_reference']['rms'], c['pre_test_one_reference']['rms']),
            ('pre-test, largest', m['condition_5_leave_one_out']['against_one_reference']['largest'], c['pre_test_one_reference']['largest']),
        ]
        for tag in ('from_X', 'from_the_optimum'):
            for a in ('Zn_N2_S', 'Zn_N2_H', 'S_N2_H'):
                rows.append((f'condition 6, {tag}, {a} minus the optimum', m['condition_6_minimum']['minima'][tag]['difference'][a],
                             c['c6_minimum']['minima'][tag]['difference'][a]))
        out['candidates'][name] = dict(
            designer_file=str(ex.OUT / f'designer_preview_{name}' / 'n2_gate.json'), designer_verdict=g['verdict'], designer_failed=g['failed'],
            designer_passed={k: v.get('passed') for k, v in c.items()},
            mine_holds=dict(c3=m['condition_3_part_A']['holds'], c4=m['condition_4_bounds']['holds'], c5=m['condition_5_leave_one_out']['holds'],
                            c6=m['condition_6_minimum']['holds']),
            numbers=[dict(quantity=q, examiner=a, designer=b, difference=round(a - b, 4)) for q, a, b in rows],
            largest_difference=max(abs(a - b) for _, a, b in rows),
            designer_linear_model_against_topology=g['linear_model_against_topology'], designer_model_sha256=g['model']['sha256'])
        print(name, '| designer:', g['verdict'], g['failed'], '| largest difference of a number:', round(out['candidates'][name]['largest_difference'], 4))
        for r in out['candidates'][name]['numbers']:
            print('    ', r)
    # ---- the builder's sensitivity, again (NOT the rule): fitting points with short hydrogen contacts left out
    rows, notes, split, mm = ex.load(final=False)
    x_opt, _ = ex.optimum()
    r0 = {k: np.radians(v) for k, v in ex.optimum_angles().items()}
    short = []
    for r in rows:
        f = (ex.GR if r['set'] == 'old' else ex.N2) / 'geom' / f"{r['id']}.xyz"
        l = f.read_text().splitlines()
        sym = [t.split()[0].capitalize() for t in l[2:66]]
        x = np.array([[float(v) for v in t.split()[1:4]] for t in l[2:66]])
        fix, lig = np.arange(46), np.arange(46, 64)
        d = np.linalg.norm(x[lig][:, None, :] - x[fix][None, :, :], axis=2)
        hl, hf = np.array([sym[i] == 'H' for i in lig]), np.array([sym[i] == 'H' for i in fix])
        heavy_h = min(d[~hl][:, hf].min(), d[hl][:, ~hf].min())
        h_h = d[hl][:, hf].min()
        r['heavy_H'], r['H_H'] = float(heavy_h), float(h_h)
        if heavy_h < 2.1 or h_h < 1.9:
            short.append(r['id'])
    sens = {}
    for name, c in ex.CAND.items():
        keep = [r for r in rows if r['id'] not in short]
        f = ex.fit_and_loo(keep, ex.LADDER[c['ladder']], c['ridge'], c['rest'], r0)
        ref = min(keep, key=lambda r: r['e_qm'])
        ir = keep.index(ref)
        one = np.delete(f['loo'] - f['loo'][ir], ir)
        sens[name] = dict(rows=len(keep), left_out=len(short), residual_rms=round(ex.rms(f['residuals']), 4), leave_one_out_rms=round(ex.rms(f['loo']), 4),
                          against_one_reference=dict(rms=round(ex.rms(one), 4), largest=round(float(one[np.argmax(np.abs(one))]), 4)))
    out['sensitivity_short_hydrogen_contacts'] = dict(
        note='NOT the rule and NOT a candidate: fitting points with a contact heavy atom to hydrogen below 2.1 A or hydrogen to hydrogen below '
             '1.9 A (ligand against zinc site) left out; limits chosen by the builder after fitting residuals had been seen; fitting data only',
        left_out=short, old_points_among_them=[p for p in short if p in split['old_rigid_points']], same_model_L5_ridge_0_01=sens)
    print('short hydrogen contacts:', len(short), short)
    print('sensitivity:', json.dumps(sens))
    # ---- the builder's diagnosis of the class, again (NOT candidates): plain least squares without bounds and ridge
    def hat_loo(A, y):
        p_, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = A @ p_ - y
        h = np.einsum('ij,ji->i', A, np.linalg.pinv(A))
        return round(ex.rms(r), 4), round(ex.rms(r / (1.0 - h)), 4)

    def wide(rws, kind, rest):
        A, names, _ = ex.design(rws, ex.LADDER['L5'], r0)
        it = [r['mm']['internal'] for r in rws]
        cols = [A[:, k] for k in range(A.shape[1])]
        if kind in ('sines', 'all'):
            for key in ('phi_HNSC', 'phi_ZnNSC'):
                for n in (1, 2, 3):
                    cols.append(np.array([np.sin(n * i[key]) for i in it]))
        if kind == 'all':
            keep = [k for k, n in enumerate(names) if not n.startswith('tor_ZnN_')]
            cols = [cols[k] for k in keep] + cols[len(names):]
            for key in ('phi_His_S', 'phi_His_H'):
                for hsd in range(3):
                    for n in (1, 2, 3):
                        cols.append(np.array([np.cos(n * i[key][hsd]) for i in it]))
                        cols.append(np.array([np.sin(n * i[key][hsd]) for i in it]))
        A2 = np.array(cols).T
        y = np.array([r['e_qm'] - r['mm'][rest] for r in rws])
        res, loo = hat_loo(A2, y)
        return dict(parameters=A2.shape[1] - 1, rows=len(rws), residual_rms=res, leave_one_out_rms=loo)
    limit = {}
    for name, c in ex.CAND.items():
        limit[name] = {f'{sub}, {kind}': wide(rws, kind, c['rest'])
                       for sub, rws in (('all fitting rows', rows), ('without the short hydrogen contacts', [r for r in rows if r['id'] not in short]))
                       for kind in ('ladder L5', 'sines', 'all')}
    out['class_limit'] = dict(
        note='NOT candidates: plain least squares without bounds and without ridge, leave-one-out from the hat matrix; "sines" = L5 with '
             'sin n phi, n = 1 to 3, on both series about N2-S (23 parameters); "all" = the same with separate terms n = 1 to 3 for each '
             'histidine instead of the common threefold term (55 parameters); fitting data only', by_charge_set=limit)
    for name, v in limit.items():
        print('class limit', name, json.dumps(v))
    (ex.OUT / 'comparison_with_the_designer.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
