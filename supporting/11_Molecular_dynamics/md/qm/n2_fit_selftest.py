"""Gate G-R2: self-test of the fit (md/qm/n2_fit.py, n2_common.fit_model) on SYNTHETIC energies.

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/n2_fit_selftest.py SCRATCH [--draws 40]
        [--proper-test-energies]      (for the test of md/qm/n2_gate.py in the same scratch tree: the test points
                                       get energies of the known model instead of poison)

No quantum-chemical energy of a new point is read: the script makes a COPY of the programs and of the geometry
records in SCRATCH/md/qm, writes synthetic result files there for every point (old, fitting AND test), and runs
the copy.  Nothing is written under md/qm/n2/ except the record md/qm/n2/selftest/fit_selftest.json.

Synthetic energies: E(i) = E(rest, i) + sum_j p*_j f_j(i) + constant + noise, on the real geometries (their
features and E(rest) from md/qm/n2/points.json), with a known parameter set p* inside the bounds of the plan.
The result files of the TEST points hold a poison energy (+1000 Hartree): any program that read one would give
nonsense, and the test also records every file the fit opens.

  T1  exact data, every model of the ladder with its own part of p*: the parameters come back.
  T2  noise 0.3 kcal/mol, `draws` draws (seed 20261006), models L2 and L5: mean and scatter of every recovered
      parameter against p*, and the standard errors of n2_common.standard_errors against that scatter.
  T3  md/qm/n2_fit.py itself (the copy), once, on data of the L2 part of p* with noise 0.3: it runs, reads no
      test file, chooses a model; its leave-one-out error is near the noise.
  T4  the door to the energies: load_points('fit') gives no test point although their files exist;
      load_points('gate') stops without N2_FINAL_TEST=1; a result file that contradicts the split stops the
      loader; a result of a test level stops it.
  T5  data from a parameter set OUTSIDE the bounds (an equilibrium angle 25 deg from the optimum, a coefficient of
      7): the fit stays inside, says which bounds it rests on, and the numbers as written to the frcmod (three
      decimals) are inside the window as the gate reads it.
"""
import argparse
import builtins
import importlib
import json
import os
import pathlib
import shutil
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
SEED = 20261006
NOISE = 0.3
P_STAR = dict(  # AMBER numbers: K, theta0 (deg); coefficients (kcal/mol)
    angles=dict(ZnN2S=(40.0, 112.0), ZnN2H=(30.0, 118.0), SN2H=(45.0, 113.0)),
    coef=dict(H_NS_C4_1=-1.2, H_NS_C4_2=2.0, H_NS_C4_3=0.5, Zn_NS_C4_1=-2.0, Zn_NS_C4_2=0.8, Zn_NS_C4_3=0.6,
              His_ZnN_S_3c=0.5, His_ZnN_S_3s=-0.4, improper=1.0, His_ZnN_H_3c=-0.3, His_ZnN_H_3s=0.2))
P_OUT = dict(angles=dict(ZnN2S=(40.0, 112.0), ZnN2H=(30.0, 122.7565 + 25.0), SN2H=(45.0, 113.0)),
             coef=dict(H_NS_C4_1=-1.2, H_NS_C4_2=7.0, H_NS_C4_3=0.5, Zn_NS_C4_1=-2.0, Zn_NS_C4_2=0.8, Zn_NS_C4_3=0.6))


def parameters(spec, names):
    p = {}
    for a, (k, th) in spec['angles'].items():
        p[f'K_{a}'] = k
        p[f'B_{a}'] = -2.0 * k * np.radians(th)
    p.update(spec['coef'])
    return {n: p.get(n, 0.0) for n in names}


def make_tree(scr, proper_test):
    q = scr / 'md' / 'qm'
    if q.exists():
        shutil.rmtree(q)
    (q / 'gr' / 'qm').mkdir(parents=True)
    for f in ('n2_common.py', 'n2_fit.py', 'n2_gate.py', 'gr_points.py'):
        shutil.copyfile(HERE / f, q / f)
    shutil.copytree(HERE / 'gr' / 'geom', q / 'gr' / 'geom')
    for f in ('points.json', 'gr_gate.json'):
        shutil.copyfile(HERE / 'gr' / f, q / 'gr' / f)
    (q / 'n2').mkdir()
    for f in ('points.json', 'split.json', 'split.json.sha256', 'qm_optimum.json'):
        shutil.copyfile(HERE / 'n2' / f, q / 'n2' / f)
    shutil.copytree(HERE / 'n2' / 'geom', q / 'n2' / 'geom')
    shutil.copytree(HERE / 'n2' / 'ff', q / 'n2' / 'ff')
    (q / 'n2' / 'qm').mkdir()
    return q


def write_energies(q, nc, par, rng, noise, proper_test, poison=1000.0, lists=('main',)):
    """Synthetic result files in the scratch tree.  Returns the constant used."""
    pj = json.loads((q / 'n2' / 'points.json').read_text())
    split = json.loads((q / 'n2' / 'split.json').read_text())
    const = -2249000.0

    def e_of(rec):
        return rec['mm_as_built']['rest'] + sum(par[k] * rec['features'][k] for k in par) + const
    d3 = {}
    for pid, rec in pj['old_rigid_points'].items():
        e = e_of(rec) + noise * rng.standard_normal()
        (q / 'gr' / 'qm' / f'{pid}.json').write_text(json.dumps(dict(id=pid, converged=True, e_tot=e / nc.HARTREE, synthetic=True)))
        d3[pid] = dict(e_d3bj=0.0)
    (q / 'gr' / 'd3.json').write_text(json.dumps(dict(points=d3, synthetic=True)))
    for f in (q / 'n2' / 'qm').glob('*.json'):
        f.unlink()
    for rec in pj['points']:
        if rec['list'] not in lists:
            continue
        pid, flag = rec['id'], split['set'][rec['id']]
        e = (e_of(rec) + noise * rng.standard_normal()) / nc.HARTREE
        if flag == 'test' and not proper_test:
            e = poison
        (q / 'n2' / 'qm' / f'{pid}.json').write_text(json.dumps(dict(
            id=pid, set=flag, converged=True, e_tot=e, d3_hartree=0.0, xyz_sha256=split['xyz_sha256'][pid], synthetic=True)))
    return const


def amber_flat(nc, par):
    a = nc.amber_terms(par)
    out = {}
    for t, v in a['angles'].items():
        if v['K']:
            out[f'K {t}'], out[f'theta0 {t}'] = v['K'], v['theta0']
    for k, v in par.items():
        if not k.startswith(('K_', 'B_')):
            out[k] = v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('scratch')
    ap.add_argument('--draws', type=int, default=40)
    ap.add_argument('--proper-test-energies', action='store_true')
    ap.add_argument('--skip-program', action='store_true', help='leave out T3 (the slow part)')
    ap.add_argument('--record', default=None, help='where the record goes (default md/qm/n2/selftest/fit_selftest.json)')
    args = ap.parse_args()
    scr = pathlib.Path(args.scratch).resolve()
    assert HERE not in scr.parents and scr != HERE, 'the scratch tree must not be inside md/qm'
    q = make_tree(scr, args.proper_test_energies)
    sys.path[:] = [p for p in sys.path if pathlib.Path(p or '.').resolve() != HERE]
    sys.path.insert(0, str(q))
    for m in ('n2_common', 'n2_fit', 'gr_points'):
        sys.modules.pop(m, None)
    nc = importlib.import_module('n2_common')
    assert pathlib.Path(nc.__file__).resolve().parent == q, 'the copy was not imported'
    _, th_opt = nc.qm_optimum()
    rng = np.random.default_rng(SEED)
    rec = dict(what='self-test of the fit of gate G-R2 on SYNTHETIC energies (md/qm/n2_fit_selftest.py); no '
                    'quantum-chemical energy of a new point was read', written=nc.now(), seed=SEED, noise_kcal=NOISE,
               scratch=str(scr), known_parameters=P_STAR, tests={})
    ok = {}

    # ---- T1: exact data
    t1 = {}
    for name, groups in nc.LADDER:
        names = nc.model_names(groups)
        par = parameters(P_STAR, names)
        write_energies(q, nc, par, rng, 0.0, False)
        rows = nc.fit_rows('fit')
        f = nc.fit_model(rows, groups, 0.0, True, th_opt)
        want, got = amber_flat(nc, par), amber_flat(nc, f['parameters'])
        t1[name] = dict(n_rows=len(rows), n_parameters=len(names), residual_rms=f['rms'],
                        largest_difference=max(abs(got[k] - want[k]) for k in want), solver=f['solver']['success'])
    ok['T1'] = all(v['largest_difference'] < 1e-3 and v['residual_rms'] < 1e-4 for v in t1.values())
    rec['tests']['T1_exact_data'] = dict(passed=ok['T1'], tolerance=1e-3, models=t1)

    # ---- T2: noise
    t2 = {}
    for name in ('L2', 'L5'):
        groups = dict(nc.LADDER)[name]
        names = nc.model_names(groups)
        par = parameters(P_STAR, names)
        want = amber_flat(nc, par)
        got, ses, sig, zs = {k: [] for k in want}, {k: [] for k in want}, [], []
        for _ in range(args.draws):
            write_energies(q, nc, par, rng, NOISE, False)
            rows = nc.fit_rows('fit')
            f = nc.fit_model(rows, groups, 0.0, True, th_opt)
            se = nc.standard_errors(f)
            sig.append(se['sigma'])
            g = amber_flat(nc, f['parameters'])
            e = {}
            for a, t in (('ZnN2S', 'M1-Y4-sy'), ('ZnN2H', 'M1-Y4-hn'), ('SN2H', 'hn-Y4-sy')):
                e[f'K {t}'], e[f'theta0 {t}'] = se['angles'][a]['K_se'], se['angles'][a]['theta0_se']
            for k, v in se['coefficients'].items():
                e[k] = v['se']
            for k, v in se['pairs'].items():
                e[f'{k}_3c'], e[f'{k}_3s'] = v['cos_se'], v['sin_se']
            for k in want:
                got[k].append(g[k])
                ses[k].append(e[k])
                if e[k]:
                    zs.append((g[k] - want[k]) / e[k])
        tab = {k: dict(known=want[k], mean=round(float(np.mean(got[k])), 4), sd_over_draws=round(float(np.std(got[k], ddof=1)), 4),
                       mean_standard_error=round(float(np.mean([s for s in ses[k] if s])), 4)) for k in want}
        ratio = [v['mean_standard_error'] / v['sd_over_draws'] for v in tab.values()]
        bias = max(abs(v['mean'] - v['known']) / (v['sd_over_draws'] / np.sqrt(args.draws)) for v in tab.values())
        zs = np.array(zs)
        t2[name] = dict(draws=args.draws, n_rows=len(rows), sigma_mean=round(float(np.mean(sig)), 4),
                        fraction_within_two_standard_errors=round(float(np.mean(np.abs(zs) < 2.0)), 4),
                        standard_error_over_scatter_min_max=[round(min(ratio), 3), round(max(ratio), 3)],
                        largest_bias_in_standard_errors_of_the_mean=round(float(bias), 2), parameters=tab)
    ok['T2'] = all(0.90 <= v['fraction_within_two_standard_errors'] and 0.6 <= v['standard_error_over_scatter_min_max'][0]
                   and v['standard_error_over_scatter_min_max'][1] <= 1.6 and abs(v['sigma_mean'] - NOISE) < 0.05
                   and v['largest_bias_in_standard_errors_of_the_mean'] < 4.0 for v in t2.values())
    rec['tests']['T2_noise'] = dict(passed=ok['T2'], criteria='at least 90 % of the recovered parameters within two '
                                    'standard errors of the known value; mean standard error within 0.6 to 1.6 of the '
                                    'scatter over the draws; sigma within 0.05 of the noise; no mean off by 4 standard '
                                    'errors of the mean', models=t2)

    # ---- T3: the program itself, with a record of every file it opens
    names2 = nc.model_names(dict(nc.LADDER)['L2'])
    par2 = parameters(P_STAR, names2)
    write_energies(q, nc, par2, rng, NOISE, args.proper_test_energies)
    split = json.loads((q / 'n2' / 'split.json').read_text())
    test_files = {str((q / 'n2' / 'qm' / f'{p}.json').resolve()) for p, s in split['set'].items() if s == 'test'}
    n_test_files = sum(pathlib.Path(f).exists() for f in test_files)
    if not args.skip_program:
        opened = []
        real_open, real_read = builtins.open, pathlib.Path.read_text

        def spy_open(file, *a, **k):
            opened.append(str(pathlib.Path(file).resolve()) if isinstance(file, (str, os.PathLike)) else str(file))
            return real_open(file, *a, **k)

        def spy_read(self, *a, **k):
            opened.append(str(self.resolve()))
            return real_read(self, *a, **k)
        builtins.open, pathlib.Path.read_text = spy_open, spy_read
        try:
            nf = importlib.import_module('n2_fit')
            argv, sys.argv = sys.argv, ['n2_fit.py', '--out', str(scr / 'fit')]
            nf.main()
            sys.argv = argv
        finally:
            builtins.open, pathlib.Path.read_text = real_open, real_read
        doc = json.loads((scr / 'fit' / 'candidates.json').read_text())
        first = doc['candidates'][0]
        want, got = amber_flat(nc, par2), amber_flat(nc, {k: first['parameters'].get(k, 0.0) for k in first['parameters']})
        diff = {k: round(got.get(k, 0.0) - want[k], 4) for k in want}
        read_test = sorted(set(opened) & test_files)
        t3 = dict(order_of_preference=doc['order_of_preference'][:5], first=first['name'], rows=doc['points'],
                  residual_rms=first['residual_rms'], leave_one_out_rms=first['leave_one_out_rms'],
                  one_reference_leave_one_out_rms=first['one_reference']['leave_one_out_rms'],
                  first_minus_known=diff, test_result_files_present=n_test_files, test_result_files_opened=read_test,
                  files_opened=len(opened), model_files=sorted(p.name for p in (scr / 'fit').glob('MODEL_*.json'))[:3],
                  expected_to_pass=first['expected_to_pass_the_preview_without_condition_6'])
        extra = {k: round(v, 4) for k, v in first['parameters'].items() if k not in par2 and not k.startswith(('K_', 'B_'))}
        t3['terms_beyond_the_known_model'] = extra
        t3['note'] = ('the data are of the L2 part of p*; which model is first depends on the noise, because the ranking '
                      'number contains the error of one reference row; a larger model is right when its extra terms are small')
        ok['T3'] = bool(not read_test and n_test_files >= 30 and doc['points']['test_points_read'] == 0
                        and first['ladder'] != 'L1' and abs(first['leave_one_out_rms'] - NOISE) < 0.1
                        and max(abs(v) for k, v in diff.items() if not k.startswith(('K ', 'theta0'))) < 0.5
                        and max(abs(v) for k, v in diff.items() if k.startswith('theta0')) < 3.0
                        and max(abs(v) for k, v in diff.items() if k.startswith('K ')) < 8.0
                        and all(abs(v) < 0.5 for v in extra.values()))
        rec['tests']['T3_the_program'] = dict(passed=ok['T3'], **t3)

    # ---- T4: the door to the energies
    t4 = {}
    data = nc.load_points('fit')
    t4['fit_returns_test_points'] = sum(p['set'] == 'test' or split['set'].get(p['id']) == 'test' for p in data)
    t4['fit_returns_points'] = len(data)
    t4['largest_abs_energy_hartree_returned'] = max(abs(p['e_b3lyp_hartree']) for p in data)
    os.environ.pop('N2_FINAL_TEST', None)
    try:
        nc.load_points('gate')
        t4['gate_without_N2_FINAL_TEST'] = 'RETURNED DATA'
    except SystemExit as e:
        t4['gate_without_N2_FINAL_TEST'] = f'stopped: {e}'
    fit_id = next(p for p, s in split['set'].items() if s == 'fit' and (q / 'n2' / 'qm' / f'{p}.json').exists())
    f = q / 'n2' / 'qm' / f'{fit_id}.json'
    keep = f.read_text()
    f.write_text(json.dumps(dict(json.loads(keep), set='test')))
    try:
        nc.load_points('fit')
        t4['file_that_contradicts_the_split'] = 'ACCEPTED'
    except SystemExit as e:
        t4['file_that_contradicts_the_split'] = f'stopped: {e}'
    f.write_text(json.dumps(dict(json.loads(keep), test_level=True)))
    try:
        nc.load_points('fit')
        t4['file_of_a_test_level'] = 'ACCEPTED'
    except SystemExit as e:
        t4['file_of_a_test_level'] = f'stopped: {e}'
    f.write_text(keep)
    (q / 'n2' / 'split.json').write_text((q / 'n2' / 'split.json').read_text().replace('"test"', '"fit"', 1))
    try:
        nc.load_points('fit')
        t4['split_file_changed'] = 'ACCEPTED'
    except SystemExit as e:
        t4['split_file_changed'] = f'stopped: {e}'
    shutil.copyfile(HERE / 'n2' / 'split.json', q / 'n2' / 'split.json')
    ok['T4'] = bool(t4['fit_returns_test_points'] == 0 and all(str(t4[k]).startswith('stopped') for k in (
        'gate_without_N2_FINAL_TEST', 'file_that_contradicts_the_split', 'file_of_a_test_level', 'split_file_changed'))
        and (args.proper_test_energies or t4['largest_abs_energy_hartree_returned'] > 3000))
    rec['tests']['T4_the_door_to_the_energies'] = dict(passed=ok['T4'], **t4)

    # ---- T5: data from outside the bounds
    groups = dict(nc.LADDER)['L2']
    names = nc.model_names(groups)
    write_energies(q, nc, parameters(P_OUT, names), rng, 0.0, False)
    rows = nc.fit_rows('fit')
    f = nc.fit_model(rows, groups, 0.0, True, th_opt)
    a = json.loads(json.dumps(f['amber']))
    for v in a['angles'].values():
        v['K'], v['theta0'] = round(v['K'], 3), round(v['theta0'], 3)
    for terms in a['torsions'].values():
        for t in terms:
            t['PK'] = round(t['PK'], 4)
    gate_reads = abs(a['angles']['M1-Y4-hn']['theta0'] - th_opt['ZnN2H']) <= nc.BOUNDS['theta0_window_deg']
    t5 = dict(data_from=P_OUT, fitted_as_written=a, at_a_bound=[list(x) for x in nc.active_bounds(f)], residual_rms=f['rms'],
              outside_bounds_as_written=nc.within_bounds(a, th_opt), theta0_ZnN2H_as_the_gate_reads_it_is_inside=bool(gate_reads),
              theta0_ZnN2H_minus_optimum=round(a['angles']['M1-Y4-hn']['theta0'] - th_opt['ZnN2H'], 4))
    ok['T5'] = bool(not t5['outside_bounds_as_written'] and gate_reads and len(t5['at_a_bound']) >= 2 and f['rms'] > 0.05)
    rec['tests']['T5_outside_the_bounds'] = dict(passed=ok['T5'], **t5)

    rec['passed'] = bool(all(ok.values()))
    rec['copies_sha256'] = {f: nc.sha(q / f) for f in ('n2_common.py', 'n2_fit.py')}
    # leave the scratch tree with proper data for whoever tests the gate there
    if args.proper_test_energies:
        write_energies(q, nc, par2, np.random.default_rng(SEED + 1), NOISE, True)
    target = pathlib.Path(args.record) if args.record else HERE / 'n2' / 'selftest' / 'fit_selftest.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(rec, indent=1, default=float) + '\n')
    print(json.dumps({k: v['passed'] for k, v in rec['tests'].items()}), '| passed:', rec['passed'])
    return 0 if rec['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
