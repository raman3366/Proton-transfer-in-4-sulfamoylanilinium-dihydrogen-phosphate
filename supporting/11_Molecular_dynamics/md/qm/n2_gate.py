"""Gate G-R2: the acceptance test of a model of the sulfonamide nitrogen (md/PLAN.md, entry of 6 Oct 2026,
1:34 PM IST, point 4; md/N2.md, section 4).  Written on 6 Oct 2026 before any energy of a new point existed.

    AMBERHOME=/root/cqc-sfa-tools/md MD=/root/cqc-sfa-tools/md/bin/python
    $MD md/qm/n2_gate.py preview --model MODEL.json [--out DIR]     conditions 3 to 6 only; the test set is NOT read
    N2_FINAL_TEST=1 $MD md/qm/n2_gate.py final --model MODEL.json   all six conditions; ONCE (it refuses a second
                                                                    run while md/qm/n2/gate/n2_gate.json exists)

MODEL.json (written by the fitting stage): {"name": ..., "groups": [...], "ridge": ..., "parameters": {...}} in
the names of n2_common.FEATURES (K_*, B_* = -2 K theta0 in kcal/mol/rad, Fourier and improper coefficients in
kcal/mol).  The gate does not take a frcmod from anybody: it writes the frcmod from these parameters itself,
builds the small model with tleap, and reads every number it judges from that topology and from OpenMM.

THE SIX CONDITIONS, all of which must hold (E_QM = B3LYP/6-31G* + D3(BJ); force field = the topology built here):
 1  test set, one reference: r = the point of lowest E_QM among ALL points of the gate (the 36 rigid points of
    gate G-R, the new fitting points, the new test points); d_i = [E_MM(i) - E_MM(r)] - [E_QM(i) - E_QM(r)];
    root mean square of d_i over the test points (r itself left out if it is a test point) at most 1.0 kcal/mol.
 2  the largest |d_i| over the test points at most 2.5 kcal/mol.
 3  part A of gate G-R as read on 6 Oct 2026: family X = X, a2..a6, d1..d6, e_p015, e_p030, e_m120, e_p120,
    f_p08, f_m08; family Rx = Rx, g1..g6; in each family the reference is its point of lowest E_QM; root mean
    square over the 23 other points at most 1.0.  Beside it, not deciding: the same against one reference.
 4  every changed term within the bounds of the plan, read from the topology: the three angle force constants
    at N2 in 10..300 kcal/mol/rad^2; their equilibrium angles within 15 deg of the quantum-chemical optimum of
    the small model (md/qm/n2/qm_optimum.json); every torsion term about N2-S and about Zn-N2 with periodicity
    1, 2 or 3 and a coefficient of at most 5 kcal/mol; at most one improper at N2, at most 5 kcal/mol; AND
    nothing else in the topology differs from the topology as built (charges, Lennard-Jones types, bonds, all
    other angles, torsions and impropers, the list of 1-4 pairs and its scaling).
 5  leave-one-out error on the fitting set (36 old points and the new fitting points) at most 1.0 kcal/mol,
    recomputed here with the fit of n2_common.fit_model for the groups and the ridge of MODEL.json; the
    parameters of MODEL.json must be the result of that fit (largest difference of a predicted energy over the
    fitting points at most 0.01 kcal/mol), otherwise the model is not the declared one and is not accepted.
 6  the force-field minimum of the small model, all atoms free, within 5 deg of the quantum-chemical optimum
    in Zn-N2-S, Zn-N2-H and S-N2-H; minimised from start X and from the optimum itself, and BOTH must hold.
A model that misses any condition is not used.  Nothing here changes a condition after a number has been seen:
the tolerances are the constants below.

ADDED AFTER THE DESIGN CHECK (6 Oct 2026, 2:15 to 2:29 PM IST), before any energy of a new point was read by a fit:
 P  PREVIEW only, a pre-test on FITTING data: the leave-one-out errors of condition 5 read against ONE reference
    (r = the fitting row of lowest E_QM; e_i - e_r over the other rows), root mean square at most 1.0 and
    largest at most 2.5.  It is what condition 1 will ask of the test set, asked beforehand of data the fit has
    seen only through the other points; a candidate that misses it does not go to the final test.  It makes the
    choice stricter and changes no condition of the final test, where it is reported and does not decide.
 F  FINAL only: the model must have been fitted to exactly the fitting rows the gate sees (their number against
    `n_points` of MODEL.json); otherwise the run stops BEFORE the test set is read and writes no n2_gate.json, so
    that a late result file cannot burn the single final run.
 N  reported, not deciding: the test r.m.s. over the points within 10 kcal/mol of the lowest E_QM.
A preview writes into its own folder, md/qm/n2/gate_preview/<name of the model file>/, unless --out is given.

Writes <out>/n2_gate.json (verdict, every condition with its numbers), n2_gate_points.csv (every point: set,
E_QM, E_MM, deviation), n2_terms_candidate.frcmod, and, only when a FINAL run accepts, n2_terms_accepted.frcmod
(the file for md/build_topology.py --n2-frcmod).  <out> is md/qm/n2/gate for `final`.
"""
import argparse
import json
import os
import pathlib
import shutil
import sys

import numpy as np

import n2_common as nc

TOL = dict(test_rms=1.0, test_largest=2.5, part_A_rms=1.0, loo_rms=1.0, minimum_deg=5.0, declared_fit_kcal=0.01,
           linear_against_topology_kcal=0.01, pre_test_rms=1.0, pre_test_largest=2.5)
FAMILY_X = ['X', 'a2', 'a3', 'a4', 'a5', 'a6', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'e_p015', 'e_p030', 'e_m120',
            'e_p120', 'f_p08', 'f_m08']
FAMILY_RX = ['Rx', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6']


def topology_terms(prmtop):
    """Everything of an AMBER topology that the comparison needs, keyed by atom names."""
    import parmed
    top = parmed.load_file(str(prmtop))
    nm = lambda a: f'{a.residue.name}:{a.name}'                            # noqa: E731
    key = lambda *a: min(tuple(nm(x) for x in a), tuple(nm(x) for x in a)[::-1])   # noqa: E731
    out = dict(atoms={nm(a): (round(a.charge, 6), a.type, round(a.rmin, 5), round(a.epsilon, 6)) for a in top.atoms},
               bonds={key(b.atom1, b.atom2): (round(b.type.k, 4), round(b.type.req, 5)) for b in top.bonds},
               angles={key(a.atom1, a.atom2, a.atom3): (round(a.type.k, 4), round(a.type.theteq, 4)) for a in top.angles},
               torsions={}, impropers={}, pairs14=set())
    for d in top.dihedrals:
        at = (d.atom1, d.atom2, d.atom3, d.atom4)
        term = (int(d.type.per), round(d.type.phi_k, 5), round(d.type.phase % 360.0, 3))
        if d.improper:
            out['impropers'].setdefault(tuple(sorted(nm(x) for x in at)) + (nm(d.atom3),), []).append(term)
        else:
            out['torsions'].setdefault(key(*at), []).append(term)
            if not d.ignore_end:
                out['pairs14'].add((tuple(sorted((nm(d.atom1), nm(d.atom4)))), round(d.type.scee, 4), round(d.type.scnb, 4)))
    return out


def compare_topologies(built, new):
    """Differences between the topology as built and the candidate; which of them the plan allows."""
    n2, zn, s, h = 'SN1:N2', 'ZN1:ZN', 'SN1:S', 'SN1:HN21'
    his = {'HD1:NE2', 'HD2:NE2', 'HE1:ND1'}
    other, allowed = [], dict(angles={}, torsions={}, improper={})
    for sec in ('atoms', 'bonds'):
        if built[sec] != new[sec]:
            other.append(f'{sec} differ')
    if built['pairs14'] != new['pairs14']:
        other.append('the 1-4 pairs or their scaling differ')
    for k in set(built['angles']) | set(new['angles']):
        if built['angles'].get(k) != new['angles'].get(k):
            if k[1] == n2 and {k[0], k[2]} in ({zn, s}, {zn, h}, {s, h}) and k in built['angles'] and k in new['angles']:
                allowed['angles']['-'.join(k)] = dict(as_built=built['angles'][k], K=new['angles'][k][0], theta0=new['angles'][k][1])
            else:
                other.append(f'angle {k}')
    for k in set(built['angles']):
        if k[1] == n2 and '-'.join(k) not in allowed['angles']:
            allowed['angles']['-'.join(k)] = dict(as_built=built['angles'][k], K=new['angles'][k][0], theta0=new['angles'][k][1])
    for k in set(built['torsions']) | set(new['torsions']):
        a, b = sorted(built['torsions'].get(k, [])), sorted(new['torsions'].get(k, []))
        mid, ends = {k[1], k[2]}, {k[0], k[3]}
        ok_class = (mid == {n2, s} and (ends & {h, zn}) and (ends - {h, zn}) <= {'SN1:C4', 'SN1:O1', 'SN1:O2'}) or \
                   (mid == {zn, n2} and len(ends & his) == 1 and len(ends & {s, h}) == 1)
        if ok_class and k in new['torsions']:
            allowed['torsions']['-'.join(k)] = dict(as_built=a, terms=b, changed=a != b)
        elif a != b:
            other.append(f'torsion {k}')
    for k in set(built['impropers']) | set(new['impropers']):
        a, b = sorted(built['impropers'].get(k, [])), sorted(new['impropers'].get(k, []))
        if a != b:
            if k[-1] == n2 and set(k[:4]) == {zn, n2, s, h} and not a:
                allowed['improper']['-'.join(k[:4])] = b
            else:
                other.append(f'improper {k}')
    return allowed, other


def bounds_from_topology(allowed, theta_opt):
    bad = []
    opt = {frozenset(('ZN1:ZN', 'SN1:S')): theta_opt['ZnN2S'], frozenset(('ZN1:ZN', 'SN1:HN21')): theta_opt['ZnN2H'],
           frozenset(('SN1:S', 'SN1:HN21')): theta_opt['SN2H']}
    for k, v in allowed['angles'].items():
        a = k.split('-')
        o = opt[frozenset((a[0], a[2]))]
        v['optimum'] = round(o, 3)
        v['theta0_minus_optimum'] = round(v['theta0'] - o, 3)
        if not (nc.BOUNDS['k_angle'][0] <= v['K'] <= nc.BOUNDS['k_angle'][1]):
            bad.append(f"angle {k}: K {v['K']}")
        if abs(v['theta0'] - o) > nc.BOUNDS['theta0_window_deg']:
            bad.append(f"angle {k}: theta0 {v['theta0']} is {v['theta0'] - o:+.2f} deg from the optimum")
    if len(allowed['angles']) != 3:
        bad.append(f"{len(allowed['angles'])} angles at N2 found, 3 expected")
    for k, v in allowed['torsions'].items():
        for per, pk, phase in v['terms']:
            if per not in (1, 2, 3) or abs(pk) > nc.BOUNDS['torsion_abs_max'] + 1e-9:
                bad.append(f'torsion {k}: periodicity {per}, coefficient {pk}')
    if len(allowed['improper']) > 1:
        bad.append('more than one improper at N2')
    for k, terms in allowed['improper'].items():
        if len(terms) != 1 or abs(terms[0][1]) > nc.BOUNDS['improper_abs_max'] + 1e-9:
            bad.append(f'improper {k}: {terms}')
    return bad


def minimum(ff, xyz0):
    mm, u = ff.mm, ff.u
    ctx = ff.ctx['total']
    ctx.setPositions(np.asarray(xyz0) * 0.1)
    mm.LocalEnergyMinimizer.minimize(ctx, 1e-4, 50000)
    st = ctx.getState(getEnergy=True, getPositions=True)
    x = np.array(st.getPositions(asNumpy=True).value_in_unit(u.angstrom))
    ix = {n: k for k, n in enumerate(ff.names)}
    zn, n2, s, h = (x[ix[k]] for k in (('ZN1', 'ZN'), ('SN1', 'N2'), ('SN1', 'S'), ('SN1', 'HN21')))
    z, hh, uu = (nc.unit(v - n2) for v in (zn, h, s))
    return dict(Zn_N2_S=round(nc.angle(zn, n2, s), 3), Zn_N2_H=round(nc.angle(zn, n2, h), 3), S_N2_H=round(nc.angle(s, n2, h), 3),
                signed_volume=round(float(np.dot(z, np.cross(hh, uu))), 4), Zn_N2=round(float(np.linalg.norm(zn - n2)), 4),
                energy=round(st.getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole), 4)), x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['preview', 'final'])
    ap.add_argument('--model', required=True)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    final = args.mode == 'final'
    sfx = '' if nc.FF_NAME == 'ff' else f'_{nc.FF_NAME}'
    out = (pathlib.Path(args.out) if (args.out and not final) else
           nc.N2 / 'gate' if final else nc.N2 / f'gate_preview{sfx}' / pathlib.Path(args.model).stem)
    if final and (out / 'n2_gate.json').exists():
        raise SystemExit(f'{out}/n2_gate.json exists: the final test has been run; it is run once')
    model = json.loads(pathlib.Path(args.model).read_text())
    if final:
        # F: the fitting rows of now against those the model was fitted to; checked before anything is written
        # and before the test set is read (load_points('fit') opens no test file)
        n_now = len(nc.load_points('fit'))
        if n_now != model.get('n_points'):
            raise SystemExit(f"the model was fitted to {model.get('n_points')} fitting rows, the gate now finds {n_now}: "
                             'fit again (md/qm/n2_fit.py) on the rows as they are; the test set was NOT read and '
                             'nothing was written')
    out.mkdir(parents=True, exist_ok=True)
    groups, ridge, par = tuple(model['groups']), float(model.get('ridge', 0.0)), model['parameters']
    assert set(par) == set(nc.model_names(groups)), 'the parameters of MODEL.json are not those of its groups'
    rec_opt, th_opt = nc.qm_optimum()
    pj = json.loads((nc.N2 / 'points.json').read_text())
    feat = {pid: v for pid, v in pj['old_rigid_points'].items()}
    feat.update({p['id']: p for p in pj['points']})
    mm0 = nc.mm_energies()                                                 # of the force-field copy in use (N2_FF)
    if model.get('force_field', 'ff') != nc.FF_NAME:
        raise SystemExit(f"the model was fitted with the force-field copy {model.get('force_field')!r}, this run uses "
                         f'{nc.FF_NAME!r} (environment variable N2_FF); nothing was read or written')

    # ---- the topology of the candidate, and what differs from the topology as built
    frc = out / 'n2_terms_candidate.frcmod'
    nc.write_frcmod(par, frc, f"terms at N2, candidate {model.get('name', '')} of gate G-R2 (md/N2.md); written by md/qm/n2_gate.py")
    prmtop = nc.build_topology(frc, out / 'work')
    built_top = nc.build_topology(None, out / 'work_as_built')
    allowed, other = compare_topologies(topology_terms(built_top), topology_terms(prmtop))
    kept = topology_terms(nc.FFDIR / 'small_as_built.prmtop') == topology_terms(built_top)
    bad_bounds = bounds_from_topology(allowed, th_opt)
    c4 = dict(passed=bool(not bad_bounds and not other and kept), outside_bounds=bad_bounds,
              differences_the_plan_does_not_allow=other, rebuilt_as_built_equals_the_kept_topology=kept,
              angles=allowed['angles'], torsions=allowed['torsions'], improper=allowed['improper'],
              bounds=nc.BOUNDS, quantum_chemical_optimum=rec_opt['angles_at_N2'])
    ff = nc.FF(prmtop)

    # ---- the data
    data = nc.load_points('gate' if final else 'fit')
    rows = []
    for p in data:
        _, xyz = nc.read_xyz(p['xyz'])
        f = feat[p['id']]
        e_qm = (p['e_b3lyp_hartree'] + p['d3_hartree']) * nc.HARTREE
        e_mm = ff.energy(xyz)['total']
        lin = mm0[p['id']]['rest'] + sum(par[k] * f['features'][k] for k in par)
        rows.append(dict(id=p['id'], set=p['set'], e_qm=e_qm, e_qm_no_d3=p['e_b3lyp_hartree'] * nc.HARTREE, e_mm=e_mm,
                         linear=lin, features=f['features'], y=e_qm - mm0[p['id']]['rest'],
                         region=f.get('region', ''), protein_contact=f.get('protein_contact', '')))
    off = np.array([r['e_mm'] - r['linear'] for r in rows])
    lin_spread = float(off.max() - off.min())
    by = {r['id']: r for r in rows}
    ref = min(rows, key=lambda r: r['e_qm'])
    for r in rows:
        r['dev'] = (r['e_mm'] - ref['e_mm']) - (r['e_qm'] - ref['e_qm'])
        r['dev_no_d3'] = (r['e_mm'] - ref['e_mm']) - (r['e_qm_no_d3'] - ref['e_qm_no_d3'])
    rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else None   # noqa: E731

    # ---- conditions 1 and 2 (final only)
    if final:
        t = [r for r in rows if r['set'] == 'test' and r is not ref]
        d = np.array([r['dev'] for r in t])
        worst = max(t, key=lambda r: abs(r['dev']))
        c1 = dict(passed=bool(len(t) >= 20 and rms(d) <= TOL['test_rms']), rms=round(rms(d), 4), n=len(t), criterion=TOL['test_rms'],
                  reference=ref['id'], reference_set=ref['set'], at_least_20_test_points=bool(len(t) >= 20),
                  not_deciding=dict(rms_about_the_mean_deviation=round(rms(d - d.mean()), 4), mean_deviation=round(float(d.mean()), 4),
                                    rms_without_D3=round(rms([r['dev_no_d3'] for r in t]), 4),
                                    rms_by_region={g: (round(rms([r['dev'] for r in t if r['region'] == g]), 4)
                                                       if any(r['region'] == g for r in t) else None) for g in ('X', 'between', 'R')},
                                    rms_without_protein_contact=(round(rms([r['dev'] for r in t if not r['protein_contact']]), 4)
                                                                 if any(not r['protein_contact'] for r in t) else None),
                                    rms_within_10_kcal_of_the_lowest_E_QM=(
                                        round(rms([r['dev'] for r in t if r['e_qm'] - ref['e_qm'] <= 10.0]), 4)
                                        if any(r['e_qm'] - ref['e_qm'] <= 10.0 for r in t) else None),
                                    n_within_10_kcal_of_the_lowest_E_QM=sum(r['e_qm'] - ref['e_qm'] <= 10.0 for r in t)))
        c2 = dict(passed=bool(abs(worst['dev']) <= TOL['test_largest']), largest=round(worst['dev'], 4), at=worst['id'],
                  criterion=TOL['test_largest'])
    else:
        c1 = c2 = dict(passed=None, note='preview: the test set was not read')

    # ---- condition 3: part A of gate G-R, two families
    dev_a, fam = [], {}
    for name, ids in (('X', FAMILY_X), ('Rx', FAMILY_RX)):
        pts = [by[i] for i in ids]
        r0 = min(pts, key=lambda r: r['e_qm'])
        dd = {r['id']: (r['e_mm'] - r0['e_mm']) - (r['e_qm'] - r0['e_qm']) for r in pts if r is not r0}
        fam[name] = dict(reference=r0['id'], deviations={k: round(v, 4) for k, v in dd.items()})
        dev_a += list(dd.values())
    all_a = [by[i] for i in FAMILY_X + FAMILY_RX]
    r1 = min(all_a, key=lambda r: r['e_qm'])
    single = [(r['e_mm'] - r1['e_mm']) - (r['e_qm'] - r1['e_qm']) for r in all_a if r is not r1]
    gate_old = json.loads((nc.GR / 'gr_gate.json').read_text())
    c3 = dict(passed=bool(len(dev_a) == 23 and rms(dev_a) <= TOL['part_A_rms']), rms=round(rms(dev_a), 4), n=len(dev_a),
              criterion=TOL['part_A_rms'], families=fam, largest=round(float(max(dev_a, key=abs)), 4),
              references_as_on_6_Oct=bool(fam['X']['reference'] == 'e_p120' and fam['Rx']['reference'] == 'g3'),
              not_deciding=dict(rms_single_reference=round(rms(single), 4), single_reference=r1['id'],
                                as_built_on_6_Oct=gate_old.get('with_d3', {}).get('part_A', {}).get('rms')))

    # ---- condition 5: leave-one-out on the fitting set, and that the model is the declared fit
    fit_rows = [r for r in rows if r['set'] in ('old', 'fit')]
    refit = nc.fit_model(fit_rows, groups, ridge, True, th_opt)
    given = np.array([r['linear'] for r in fit_rows])
    again = np.array([nc._predict_one(refit, r, th_opt) + r['e_qm'] - r['y'] for r in fit_rows])
    dd = (given - again) - (given - again).mean()
    loo = nc.leave_one_out(fit_rows, groups, ridge, True, th_opt)
    n_par = len(par)
    c5 = dict(passed=bool(loo['rms'] <= TOL['loo_rms'] and np.abs(dd).max() <= TOL['declared_fit_kcal']
                          and len(fit_rows) >= 8 * n_par and refit['solver']['success']),
              leave_one_out_rms=round(loo['rms'], 4), leave_one_out_largest=round(loo['largest'], 4), criterion=TOL['loo_rms'],
              n_fitting_points=len(fit_rows), n_parameters=n_par, points_per_parameter=round(len(fit_rows) / n_par, 2),
              at_least_8_points_per_parameter=bool(len(fit_rows) >= 8 * n_par),
              residual_rms_of_the_fit=round(refit['rms'], 4), ridge=ridge, groups=list(groups),
              model_is_the_declared_fit=bool(np.abs(dd).max() <= TOL['declared_fit_kcal']),
              largest_difference_given_against_refitted_kcal=round(float(np.abs(dd).max()), 5), solver=refit['solver'])
    # ---- P: the leave-one-out errors against one reference (fitting data only)
    one = nc.one_reference(loo['errors'], [r['e_qm'] for r in fit_rows])
    pre = dict(passed=bool(one['rms'] <= TOL['pre_test_rms'] and abs(one['largest']) <= TOL['pre_test_largest']),
               rms=round(one['rms'], 4), largest=round(one['largest'], 4), reference=fit_rows[one['reference_row']]['id'],
               mean=round(one['mean'], 4), leave_one_out_error_of_the_reference=round(one['error_of_the_reference'], 4),
               criterion_rms=TOL['pre_test_rms'], criterion_largest=TOL['pre_test_largest'],
               note='fitting data only; decides in the preview, reported in the final test')

    # ---- condition 6: the minimum of the small model against the quantum-chemical optimum
    _, x_start = nc.read_xyz(nc.GR / 'geom' / 'X.xyz')
    mins = {}
    for tag, x0 in (('from_X', x_start), ('from_the_optimum', np.array(rec_opt['xyz_angstrom']))):
        mres, _ = minimum(ff, x0)
        mres['difference'] = {k: round(mres[k] - rec_opt['angles_at_N2'][k], 3) for k in ('Zn_N2_S', 'Zn_N2_H', 'S_N2_H')}
        mres['within'] = bool(all(abs(v) <= TOL['minimum_deg'] for v in mres['difference'].values()))
        mins[tag] = mres
    c6 = dict(passed=bool(all(m['within'] for m in mins.values())), criterion_deg=TOL['minimum_deg'],
              optimum=rec_opt['angles_at_N2'], optimum_source=rec_opt['source'], optimum_sha256=rec_opt['source_sha256'],
              minima=mins)

    consistent = bool(lin_spread <= TOL['linear_against_topology_kcal'])
    conds = dict(c1_test_rms=c1, c2_test_largest=c2, c3_part_A=c3, c4_bounds_and_nothing_else=c4, c5_leave_one_out=c5,
                 c6_minimum=c6)
    need = list(conds) if final else ['c3_part_A', 'c4_bounds_and_nothing_else', 'c5_leave_one_out', 'c6_minimum',
                                      'pre_test_one_reference']
    conds['pre_test_one_reference'] = pre
    accepted = bool(consistent and all(conds[k]['passed'] is True for k in need))
    doc = dict(
        what='gate G-R2 (md/PLAN.md, 6 Oct 2026, 1:34 PM IST, point 4), ' + ('FINAL test' if final else
             'PREVIEW without the test set: not a verdict'), mode=args.mode, written=nc.now(),
        model=dict(file=str(args.model), sha256=nc.sha(args.model), name=model.get('name'), groups=list(groups), ridge=ridge,
                   parameters=par, amber=nc.amber_terms(par)),
        verdict=('ACCEPTED' if accepted else 'NOT ACCEPTED') if final else
                ('conditions 3 to 6 and the pre-test hold' if accepted else 'conditions 3 to 6 and the pre-test do not all hold'),
        accepted=accepted if final else None, conditions=conds,
        failed=[k for k in need if conds[k]['passed'] is not True] + ([] if consistent else ['linear model against topology']),
        linear_model_against_topology=dict(spread_kcal=round(lin_spread, 6), tolerance=TOL['linear_against_topology_kcal'],
                                           consistent=consistent),
        points=dict(old=sum(r['set'] == 'old' for r in rows), fit=sum(r['set'] == 'fit' for r in rows),
                    test=sum(r['set'] == 'test' for r in rows), lowest_E_QM=ref['id']),
        not_deciding=dict(rms_fitting_points_one_reference=round(rms([r['dev'] for r in fit_rows if r is not ref]), 4),
                          Rx_against_X=dict(QM=round(by['Rx']['e_qm'] - by['X']['e_qm'], 3), MM=round(by['Rx']['e_mm'] - by['X']['e_mm'], 3),
                                            as_built_MM=round(mm0['Rx']['total'] - mm0['X']['total'], 3)),
                          force_field_copy=nc.FF_NAME),
        tolerances=TOL, files=dict(frcmod=str(frc), frcmod_sha256=nc.sha(frc), prmtop=str(prmtop), split_sha256=nc.sha(nc.N2 / 'split.json'),
                                   points_sha256=nc.sha(nc.N2 / 'points.json')),
        versions=dict(numpy=np.__version__, openmm=ff.mm.__version__, python=sys.version.split()[0]))
    nc.write_json(out / 'n2_gate.json', doc)
    cols = ['id', 'set', 'region', 'protein_contact', 'E_QM_rel', 'E_MM_rel', 'MM_minus_QM', 'MM_minus_QM_no_D3', 'linear_minus_topology']
    lines = [','.join(cols)]
    for r, o in zip(rows, off):
        lines.append(','.join(str(v) for v in (r['id'], r['set'], r['region'], r['protein_contact'], round(r['e_qm'] - ref['e_qm'], 4),
                                               round(r['e_mm'] - ref['e_mm'], 4), round(r['dev'], 4), round(r['dev_no_d3'], 4),
                                               round(float(o - off.mean()), 6))))
    (out / 'n2_gate_points.csv').write_text('\n'.join(lines) + '\n')
    if final and accepted:
        shutil.copyfile(frc, out / 'n2_terms_accepted.frcmod')
    print(doc['what'])
    print('verdict:', doc['verdict'], '| failed:', doc['failed'])
    for k, c in conds.items():
        brief = {kk: vv for kk, vv in c.items() if kk in ('passed', 'rms', 'n', 'largest', 'at', 'leave_one_out_rms',
                                                          'n_fitting_points', 'n_parameters', 'outside_bounds',
                                                          'differences_the_plan_does_not_allow', 'note')}
        print(' ', k, brief)
    print('  minima:', {k: (v['difference'], v['within']) for k, v in mins.items()})
    return 0 if accepted else 1


if __name__ == '__main__':
    sys.exit(main())
