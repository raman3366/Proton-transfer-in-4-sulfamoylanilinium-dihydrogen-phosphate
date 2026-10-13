#!/usr/bin/env python
"""Gate G-R2: conditions 1 and 2, the only two that need the TEST SET (md/PLAN.md, 6 Oct 2026, 1:34 PM IST, point 4).

    nice -n 10 /root/cqc-sfa-tools/md/bin/python md/qm/n2/final_test/n2_exam_final.py --dry-run
        a test of this program on FITTING data: every fourth new fitting point plays the test set; no file of a
        test point is opened; output in the examiner's scratch folder only, marked as a test of the program.
    N2_FINAL_TEST=1 nice -n 10 /root/cqc-sfa-tools/md/bin/python md/qm/n2/final_test/n2_exam_final.py --final
        THE READ OF THE TEST SET.  ONCE: it refuses when md/qm/n2/final_test/test_set_read.json exists.

STATE ON 6 OCT 2026: --final HAS NOT BEEN RUN.  The examiner left the test set unread, because the verdict (NOT
ACCEPTED) is fixed by conditions 3, 5 and 6 on fitting data and the read cannot be taken back (md/N2.md, section 9).
Whoever runs --final spends the test set; that is a decision for the user, in a dated line of md/PLAN.md.

What --final does, exactly as the plan words it: E_QM = B3LYP/6-31G* + D3(BJ) two-body; r = the point of lowest
E_QM among ALL points of the gate (36 old, new fitting, new test); d_i = [E_MM(i) - E_MM(r)] - [E_QM(i) - E_QM(r)];
condition 1: root mean square of d_i over the test points (r left out if it is one) at most 1.0, with at least 20
test points; condition 2: the largest |d_i| at most 2.5.  For both candidates and, beside, for the terms as built.
The test result files are taken from the watchdog's copy of the machine's folder and are checked before use: set
flag, checksum of the geometry, converged, no test level, the settings of the gate.
"""
import argparse
import json
import os
import pathlib
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import n2_exam_conditions as ex                                            # noqa: E402

DONE = ex.OUT / 'test_set_read.json'


def evaluate(rows, pseudo=None):
    """rows: every point of the gate with its energies; pseudo: identifiers that play the test set in a dry run."""
    is_test = (lambda r: r['id'] in pseudo) if pseudo is not None else (lambda r: r['set'] == 'test')
    ref = min(rows, key=lambda r: r['e_qm'])
    out = dict(reference=ref['id'], reference_set=ref['set'], n_rows=len(rows), candidates={})
    for name, c in ex.CAND.items():
        res = {}
        for label, key in (('new', c['new']), ('as_built', c['built'])):
            t = [r for r in rows if is_test(r) and r is not ref]
            d = np.array([(r['mm'][key] - ref['mm'][key]) - (r['e_qm'] - ref['e_qm']) for r in t])
            dn = np.array([(r['mm'][key] - ref['mm'][key]) - (r['e_qm_no_d3'] - ref['e_qm_no_d3']) for r in t])
            w = int(np.argmax(np.abs(d)))
            res[label] = dict(n=len(t), rms=round(ex.rms(d), 4), largest=round(float(d[w]), 4), largest_at=t[w]['id'],
                              mean=round(float(d.mean()), 4), rms_about_the_mean=round(ex.rms(d - d.mean()), 4),
                              rms_without_D3=round(ex.rms(dn), 4),
                              rms_by_region={g: (round(ex.rms([x for x, r in zip(d, t) if r['region'] == g]), 4), sum(r['region'] == g for r in t))
                                             for g in ('X', 'between', 'R')},
                              deviations={r['id']: round(float(x), 4) for r, x in zip(t, d)})
        n = res['new']
        out['candidates'][name] = dict(
            condition_1=dict(holds=bool(n['n'] >= 20 and n['rms'] <= ex.LIMITS['test_rms']), rms=n['rms'], n=n['n'], limit=ex.LIMITS['test_rms']),
            condition_2=dict(holds=bool(abs(n['largest']) <= ex.LIMITS['test_largest']), largest=n['largest'], at=n['largest_at'],
                             limit=ex.LIMITS['test_largest']), new=n, as_built=res['as_built'])
    return out


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--dry-run', action='store_true')
    g.add_argument('--final', action='store_true')
    args = ap.parse_args()
    if args.dry_run:
        rows, notes, split, _ = ex.load(final=False)
        new = [p for p in split['order'] if split['set'][p] == 'fit' and any(r['id'] == p for r in rows)]
        pseudo = set(new[3::4])
        res = evaluate(rows, pseudo)
        res.update(what='TEST OF THE PROGRAM on fitting data: every fourth new fitting point plays the test set (in-sample for the '
                        'candidates); NOT a result of the gate; no file of a test point was opened', written=ex.now(), pseudo_test_points=len(pseudo))
        out = ex.SCR / 'dry_run_final.json'
        out.write_text(json.dumps(res, indent=1) + '\n')
        for k, v in res['candidates'].items():
            print(k, 'pseudo-test:', {a: v[a] for a in ('condition_1', 'condition_2')}, '| as built rms', v['as_built']['rms'])
        print('reference', res['reference'], '| written to', out)
        return
    if os.environ.get('N2_FINAL_TEST') != '1':
        raise SystemExit('the test set is read once, with N2_FINAL_TEST=1 in the environment, on the decision of the user')
    if DONE.exists():
        raise SystemExit(f'{DONE} exists: the test set has been read; it is read once')
    cond = json.loads((ex.OUT / 'conditions_fitting_data.json').read_text())
    with open(ex.OUT / 'test_set_access.log', 'a') as fh:
        fh.write(f"{ex.now()}  pid {os.getpid()}  {' '.join(sys.argv)}\n")
    gate_settings = json.loads((ex.GR / 'qm' / 'X.json').read_text())['settings']
    split = json.loads((ex.N2 / 'split.json').read_text())
    bad = []
    for pid in split['order']:
        f = ex.WATCHDOG / f'{pid}.json'
        if split['set'][pid] == 'test' and f.exists():
            r = json.loads(f.read_text())
            if r.get('settings') != gate_settings or not r.get('converged') or r.get('test_level') or r.get('second_order_solver') is None:
                bad.append(pid)
    rows, notes, split, _ = ex.load(final=True)
    res = evaluate(rows)
    res.update(what='gate G-R2, conditions 1 and 2 on the test set, read ONCE', written=ex.now(), data_notes=notes,
               test_files_with_other_settings_or_not_converged=bad,
               rows=dict(old=sum(r['set'] == 'old' for r in rows), fit=sum(r['set'] == 'fit' for r in rows), test=sum(r['set'] == 'test' for r in rows)),
               conditions_on_fitting_data=str(ex.OUT / 'conditions_fitting_data.json'), sha256_of_that_file=ex.sha(ex.OUT / 'conditions_fitting_data.json'))
    DONE.write_text(json.dumps(res, indent=1) + '\n')
    for k, v in res['candidates'].items():
        print(k, {a: v[a] for a in ('condition_1', 'condition_2')}, '| as built rms', v['as_built']['rms'], 'largest', v['as_built']['largest'])
    print('reference', res['reference'], res['reference_set'], '| rows', res['rows'], '| notes', notes, '| bad test files', bad)
    print('written:', DONE, '(the verdict file md/qm/n2/n2_gate_result.json has to be brought up to date with these two conditions)')


if __name__ == '__main__':
    main()
