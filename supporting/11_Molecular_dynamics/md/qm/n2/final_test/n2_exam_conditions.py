#!/usr/bin/env python
"""Gate G-R2, the examiner's evaluation of the acceptance conditions (md/PLAN.md, 6 Oct 2026, 1:34 PM IST, point 4;
md/N2.md, section 9).  Own route throughout: nothing is imported from md/qm/n2_common.py, n2_fit.py or n2_gate.py.

    AMBERHOME=/root/cqc-sfa-tools/md nice -n 10 /root/cqc-sfa-tools/md/bin/python md/qm/n2/final_test/n2_exam_conditions.py
        conditions 3, 4, 5 and 6 and everything else that needs FITTING data only (36 old points and the new points
        flagged "fit"); the result file of a test point is not opened.
    ... n2_exam_conditions.py --final          (needs N2_FINAL_TEST=1; reads the test set ONCE; see final_test())

Inputs: md/qm/n2/final_test/mm_energies.json and d3_check.json (n2_exam_mm.py: force-field energies of four
topologies built by the examiner, D3(BJ) by the gate's routine; geometry only), the B3LYP energies of the result
files (key e_tot), md/qm/n2/split.json, /root/cqc-sfa-md/qm/prod/opt_final.npy (the quantum-chemical optimum).

Two candidates are examined, each with its own charge set:
    'charges_as_built'   T1 against T0   : md/qm/n2/n2_terms_new_charges_as_built.frcmod (L5, ridge 0.01)
    'second_step'        T2 against T0b  : md/qm/n2/n2_terms_new.frcmod (L5, ridge 0.01; charges of two placements)

E_QM = B3LYP/6-31G* + D3(BJ) two-body (the gate's rule), kcal/mol.  Deviation = force field minus quantum chemistry.
"""
import argparse
import csv
import hashlib
import json
import os
import pathlib
import sys
import time

import numpy as np

REPO = pathlib.Path('/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate')
QM = REPO / 'md' / 'qm'
N2 = QM / 'n2'
GR = QM / 'gr'
OUT = N2 / 'final_test'
SCR = pathlib.Path('/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/'
                   '82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/n2/final')
WATCHDOG = pathlib.Path('/root/cqc-sfa-md/from_gpu/repo_md/qm/n2/qm')
HARTREE = 627.5094740631
LIMITS = dict(test_rms=1.0, test_largest=2.5, part_A_rms=1.0, loo_rms=1.0, minimum_deg=5.0,
              k_angle=(10.0, 300.0), theta0_window_deg=15.0, torsion_abs_max=5.0, improper_abs_max=5.0)
FAMILY_X = ['X', 'a2', 'a3', 'a4', 'a5', 'a6', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'e_p015', 'e_p030', 'e_m120',
            'e_p120', 'f_p08', 'f_m08']
FAMILY_RX = ['Rx', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6']
LADDER = dict(L1=('angle', 'tor_H'), L2=('angle', 'tor_H', 'tor_Zn'), L3=('angle', 'tor_H', 'tor_Zn', 'tor_ZnN_S'),
              L4=('angle', 'tor_H', 'tor_Zn', 'tor_ZnN_S', 'improper'),
              L5=('angle', 'tor_H', 'tor_Zn', 'tor_ZnN_S', 'improper', 'tor_ZnN_H'))
CAND = dict(charges_as_built=dict(new='T1', built='T0', rest='rest_T0', frcmod=N2 / 'n2_terms_new_charges_as_built.frcmod',
                                  builder_fit=N2 / 'fit', ladder='L5', ridge=0.01),
            second_step=dict(new='T2', built='T0b', rest='rest_T0b', frcmod=N2 / 'n2_terms_new.frcmod',
                             builder_fit=N2 / 'fit_ff2', ladder='L5', ridge=0.01))


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def now():
    return time.strftime('%Y-%m-%d %I:%M:%S %p IST', time.gmtime(time.time() + 19800))


def write_json(path, obj):
    path = pathlib.Path(path)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(obj, indent=1) + '\n')
    os.replace(tmp, path)


def rms(v):
    v = np.asarray(list(v), float)
    return float(np.sqrt(np.mean(v ** 2))) if len(v) else None


def ang(a, b, c):
    u, v = a - b, c - b
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1.0, 1.0))))


# ---------------------------------------------------------------------------------------------------- the data
def load(final=False):
    """Rows of the gate.  Fitting data: the 36 rigid points of gate G-R (md/qm/gr/qm) and the new points flagged
    "fit" (md/qm/n2/qm).  With final=True also the test points, from the watchdog's copy of the machine's folder."""
    split = json.loads((N2 / 'split.json').read_text())
    assert sha(N2 / 'split.json') == (N2 / 'split.json.sha256').read_text().split()[0], 'split.json changed'
    mm = json.loads((OUT / 'mm_energies.json').read_text())
    assert mm['sha256']['split'] == sha(N2 / 'split.json')
    d3 = json.loads((OUT / 'd3_check.json').read_text())['points']
    old_d3 = json.loads((GR / 'd3.json').read_text())['points']
    rows, notes = [], dict(d3_stored_minus_recomputed_largest_hartree=0.0, not_converged=[], missing=[])
    for pid in split['old_rigid_points']:
        r = json.loads((GR / 'qm' / f'{pid}.json').read_text())
        if not r.get('converged'):
            notes['not_converged'].append(pid)
            continue
        notes['d3_stored_minus_recomputed_largest_hartree'] = max(notes['d3_stored_minus_recomputed_largest_hartree'],
                                                                 abs(old_d3[pid]['e_d3bj'] - d3[pid]['gate_routine']))
        rows.append(dict(id=pid, set='old', e_b3lyp=r['e_tot'], d3=d3[pid]['gate_routine']))
    for pid in split['order']:
        flag = split['set'][pid]
        if flag == 'test' and not final:
            continue                                                       # the file of a test point is not opened
        f = (N2 / 'qm' / f'{pid}.json') if flag == 'fit' else (WATCHDOG / f'{pid}.json')
        if not f.exists():
            notes['missing'].append(pid)
            continue
        r = json.loads(f.read_text())
        assert r.get('set') == flag and not r.get('test_level') and r.get('xyz_sha256') == split['xyz_sha256'][pid], pid
        if not r.get('converged'):
            notes['not_converged'].append(pid)
            continue
        notes['d3_stored_minus_recomputed_largest_hartree'] = max(notes['d3_stored_minus_recomputed_largest_hartree'],
                                                                 abs(r['d3_hartree'] - d3[pid]['gate_routine']))
        rows.append(dict(id=pid, set=flag, e_b3lyp=r['e_tot'], d3=d3[pid]['gate_routine']))
    for r in rows:
        r['e_qm'] = (r['e_b3lyp'] + r['d3']) * HARTREE
        r['e_qm_no_d3'] = r['e_b3lyp'] * HARTREE
        r['mm'] = mm['points'][r['id']]
        sv = r['mm']['internal']['signed_volume']
        r['region'] = 'R' if sv < -0.3 else 'X' if sv > 0.3 else 'between'
    return rows, notes, split, mm


# ------------------------------------------------------------------------------------- the fit, written again
def optimum():
    """The quantum-chemical optimum of the small model from the file of task qm_opt; the three angles at N2."""
    f = pathlib.Path('/root/cqc-sfa-md/qm/prod/opt_final.npy')
    x = np.load(f)
    pinned = json.loads((N2 / 'qm_optimum.json').read_text())
    xp = np.array(pinned['xyz_angstrom'])
    if np.abs(x - xp).max() > 1e-6 and np.abs(x * 0.52917721067 - xp).max() < 1e-6:
        x = x * 0.52917721067
    names = json.loads((OUT / 'mm_energies.json').read_text())
    return x, dict(file=str(f), sha256=sha(f), largest_difference_from_the_pinned_geometry=float(np.abs(x - xp).max()),
                   pinned_sha256=pinned['source_sha256'], pinned_angles=pinned['angles_at_N2'])


def design(rows, groups, r0):
    """Columns of the linear model (my own internal coordinates, n2_exam_mm.py), the constant last."""
    cols, names, pen = [], [], []
    it = [r['mm']['internal'] for r in rows]
    if 'angle' in groups:
        for a in ('ZnN2S', 'ZnN2H', 'SN2H'):
            th = np.array([i[f'th_{a}'] for i in it]) - r0[a]
            cols += [th ** 2, th]
            names += [f'K_{a}', f'b_{a}']
            pen += [0.0, 0.0]
    for g, key in (('tor_H', 'phi_HNSC'), ('tor_Zn', 'phi_ZnNSC')):
        if g in groups:
            for n in (1, 2, 3):
                cols.append(np.array([np.cos(n * i[key]) for i in it]))
                names.append(f'{g}_{n}')
                pen.append(1.0)
    for g, key in (('tor_ZnN_S', 'phi_His_S'), ('tor_ZnN_H', 'phi_His_H')):
        if g in groups:
            cols.append(np.array([sum(np.cos(3 * p) for p in i[key]) for i in it]))
            cols.append(np.array([sum(np.sin(3 * p) for p in i[key]) for i in it]))
            names += [f'{g}_c', f'{g}_s']
            pen += [1.0, 1.0]
    if 'improper' in groups:
        cols.append(np.array([1.0 - np.cos(2 * i['omega']) for i in it]))
        names.append('improper')
        pen.append(1.0)
    cols.append(np.ones(len(rows)))
    names.append('const')
    pen.append(0.0)
    return np.array(cols).T, names, np.array(pen)


def outside(p, names):
    """Which bounds of the plan a solution violates (point 2 of the plan entry)."""
    bad = []
    dth = np.radians(LIMITS['theta0_window_deg'])
    for k, n in enumerate(names):
        if n.startswith('K_'):
            if not (LIMITS['k_angle'][0] <= p[k] <= LIMITS['k_angle'][1]):
                bad.append(f'{n} = {p[k]:.3f}')
            elif abs(p[k + 1]) > 2.0 * p[k] * dth:
                bad.append(f'theta0 of {n[2:]} {np.degrees(-p[k + 1] / (2 * p[k])):+.2f} deg from the optimum')
        elif n.endswith('_c'):
            if np.hypot(p[k], p[k + 1]) > LIMITS['torsion_abs_max']:
                bad.append(f'amplitude {n[:-2]}')
        elif n not in ('const',) and not n.startswith('b_') and not n.endswith('_s'):
            if abs(p[k]) > LIMITS['torsion_abs_max']:
                bad.append(f'{n} = {p[k]:.3f}')
    return bad


def solve(A, y, pen, ridge, names):
    """Ridge least squares, penalty ridge * N * c^2 on the Fourier and improper coefficients.  The unconstrained
    solution is taken when it lies inside the bounds of the plan (the problem is convex, so it is then also the
    bounded solution); otherwise a bounded solution by SLSQP, written here."""
    lam = ridge * len(y)
    p = np.linalg.solve(A.T @ A + lam * np.diag(pen), A.T @ y)
    bad = outside(p, names)
    if not bad:
        return p, False
    from scipy.optimize import minimize
    dth = np.radians(LIMITS['theta0_window_deg'] - 0.001)
    cmax = LIMITS['torsion_abs_max']
    bnd, cons = [], []
    for k, n in enumerate(names):
        if n.startswith('K_'):
            bnd.append(LIMITS['k_angle'])
            cons.append(dict(type='ineq', fun=lambda q, k=k: 2.0 * q[k] * dth - q[k + 1]))
            cons.append(dict(type='ineq', fun=lambda q, k=k: 2.0 * q[k] * dth + q[k + 1]))
        elif n.startswith('b_') or n == 'const':
            bnd.append((None, None))
        else:
            bnd.append((-cmax, cmax))
            if n.endswith('_c'):
                cons.append(dict(type='ineq', fun=lambda q, k=k: cmax ** 2 - q[k] ** 2 - q[k + 1] ** 2))
    p0 = np.zeros(len(names))
    for k, n in enumerate(names):
        if n.startswith('K_'):
            p0[k] = 50.0
    p0[-1] = float(np.mean(y - A @ p0))
    f = lambda q: float((A @ q - y) @ (A @ q - y) + lam * (pen * q) @ q)                    # noqa: E731
    g = lambda q: 2.0 * A.T @ (A @ q - y) + 2.0 * lam * pen * q                             # noqa: E731
    res = minimize(f, p0, jac=g, method='SLSQP', bounds=bnd, constraints=cons, options=dict(maxiter=3000, ftol=1e-13))
    return res.x, True


def fit_and_loo(rows, groups, ridge, rest_key, r0):
    A, names, pen = design(rows, groups, r0)
    y = np.array([r['e_qm'] - r['mm'][rest_key] for r in rows])
    p, bounded = solve(A, y, pen, ridge, names)
    resid = A @ p - y
    loo, n_bounded = np.zeros(len(y)), 0
    for i in range(len(y)):
        keep = np.arange(len(y)) != i
        q, b = solve(A[keep], y[keep], pen, ridge, names)
        n_bounded += b
        loo[i] = A[i] @ q - y[i]
    return dict(names=names, p=p, residuals=resid, loo=loo, bounded=bounded, loo_refits_on_a_bound=int(n_bounded),
                n_parameters=len(names) - 1, prediction=A @ p - p[-1] + np.array([r['mm'][rest_key] for r in rows]))


def amber_numbers(names, p, r0):
    out = {}
    for k, n in enumerate(names):
        if n.startswith('K_'):
            out[n[2:]] = dict(K=float(p[k]), theta0=float(np.degrees(r0[n[2:]] - p[k + 1] / (2.0 * p[k]))))
        elif n.endswith('_c'):
            out[n[:-2]] = dict(PK=float(np.hypot(p[k], p[k + 1])), phase=float(np.degrees(np.arctan2(p[k + 1], p[k])) % 360.0))
        elif n != 'const' and not n.startswith('b_') and not n.endswith('_s'):
            out[n] = float(p[k])
    return out


# -------------------------------------------------------------------------------------------- the conditions
def part_A(by, key):
    """The 23 comparisons of part A of gate G-R as read on 6 Oct: two families, each against its own point of
    lowest E_QM; beside it the number against ONE reference (the lowest E_QM of the 25 points)."""
    dev, fam = [], {}
    for name, ids in (('X', FAMILY_X), ('Rx', FAMILY_RX)):
        pts = [by[i] for i in ids]
        r0 = min(pts, key=lambda r: r['e_qm'])
        d = {r['id']: (r['mm'][key] - r0['mm'][key]) - (r['e_qm'] - r0['e_qm']) for r in pts if r is not r0}
        fam[name] = dict(reference=r0['id'], deviations={k: round(v, 4) for k, v in d.items()})
        dev += list(d.values())
    al = [by[i] for i in FAMILY_X + FAMILY_RX]
    r1 = min(al, key=lambda r: r['e_qm'])
    single = [(r['mm'][key] - r1['mm'][key]) - (r['e_qm'] - r1['e_qm']) for r in al if r is not r1]
    return dict(rms=round(rms(dev), 4), n=len(dev), largest=round(float(max(dev, key=abs)), 4), families=fam,
                single_reference=dict(reference=r1['id'], rms=round(rms(single), 4), n=len(single),
                                      largest=round(float(max(single, key=abs)), 4)))


def bounds_of_topology(tag, built):
    """Condition 4 from the topology the examiner built: the terms at N2 inside the bounds of point 2, and nothing
    else different from the topology as built (md/qm/n2/final_test/topology_differences.json)."""
    td = json.loads((OUT / 'topology_differences.json').read_text())
    x_opt, _ = optimum()
    mmj = json.loads((OUT / 'mm_energies.json').read_text())
    terms = td['terms_at_N2'][tag]
    diff = td['differences'][f'{tag}_against_{built}']
    opt = optimum_angles()
    bad, ang_, tor, imp = [], [], [], []
    for t in terms:
        if t['kind'] == 'angle':
            a = t['atoms'].split('-')
            key = {frozenset(('ZN', 'S')): 'ZnN2S', frozenset(('ZN', 'HN21')): 'ZnN2H', frozenset(('S', 'HN21')): 'SN2H'}[frozenset((a[0], a[2]))]
            k, th0 = t['numbers']
            ang_.append(dict(angle=key, K=k, theta0=th0, optimum=round(opt[key], 3), theta0_minus_optimum=round(th0 - opt[key], 3)))
            if not (LIMITS['k_angle'][0] <= k <= LIMITS['k_angle'][1]):
                bad.append(f'{key}: K {k}')
            if abs(th0 - opt[key]) > LIMITS['theta0_window_deg']:
                bad.append(f'{key}: theta0 {th0}')
        else:
            per, k, ph = t['numbers']
            (imp if t['kind'] == 'improper' else tor).append(dict(atoms=t['atoms'], n=per, PK=k, phase=ph))
            if per not in (1, 2, 3) or abs(k) > LIMITS['torsion_abs_max']:
                bad.append(f"{t['kind']} {t['atoms']}: n {per}, PK {k}")
    if len(imp) > 1:
        bad.append('more than one improper term at N2')
    if len(ang_) != 3:
        bad.append(f'{len(ang_)} angles at N2')
    # every difference must be one of the terms at N2
    allowed_angles = {('ZN1:ZN', 'SN1:N2', 'SN1:S'), ('ZN1:ZN', 'SN1:N2', 'SN1:HN21'), ('SN1:S', 'SN1:N2', 'SN1:HN21')}
    other = []
    for sec in ('bonds', 'charges', 'lj', 'exceptions'):
        if diff[sec]:
            other.append(f'{len(diff[sec])} {sec} differ')
    for d in diff['angles']:
        a = tuple(d['atoms'].split('-'))
        if a not in allowed_angles and a[::-1] not in allowed_angles:
            other.append('angle ' + d['atoms'])
    for d in diff['torsions']:
        a = d['atoms'].split('-')
        mid, ends, al = {a[1], a[2]}, {a[0], a[3]}, set(a)
        proper = (mid == {'SN1:N2', 'SN1:S'} and ends & {'SN1:HN21', 'ZN1:ZN'} and (ends - {'SN1:HN21', 'ZN1:ZN'}) <= {'SN1:C4', 'SN1:O1', 'SN1:O2'}) or \
                 (mid == {'ZN1:ZN', 'SN1:N2'} and len(ends & {'HD1:NE2', 'HD2:NE2', 'HE1:ND1'}) == 1 and len(ends & {'SN1:S', 'SN1:HN21'}) == 1)
        improper = al == {'ZN1:ZN', 'SN1:N2', 'SN1:S', 'SN1:HN21'}
        if not (proper or improper):
            other.append('torsion ' + d['atoms'])
    return dict(holds=bool(not bad and not other), outside_bounds=bad, differences_the_plan_does_not_allow=other,
                angles=ang_, torsion_terms=tor, improper=imp,
                counts_of_differences={k: len(v) for k, v in diff.items()})


_OPT = {}


def optimum_angles():
    if not _OPT:
        x, _ = optimum()
        mmj = json.loads((OUT / 'mm_energies.json').read_text())
        import openmm.app as app
        top = app.AmberPrmtopFile(mmj['topologies']['T0']['prmtop'])
        ix = {(a.residue.name, a.name): k for k, a in enumerate(top.topology.atoms())}
        zn, n2, s, h = (x[ix[k]] for k in (('ZN1', 'ZN'), ('SN1', 'N2'), ('SN1', 'S'), ('SN1', 'HN21')))
        _OPT.update(ZnN2S=ang(zn, n2, s), ZnN2H=ang(zn, n2, h), SN2H=ang(s, n2, h))
    return _OPT


def minimum(prmtop, x0):
    """Minimum of the 64-atom model, all atoms free: L-BFGS (SciPy) on the OpenMM energy and forces, then OpenMM's
    own minimiser from that point as a cross-check.  Returns the three angles at N2 and the largest force."""
    import openmm
    from openmm import app, unit as u
    from scipy.optimize import minimize
    top = app.AmberPrmtopFile(str(prmtop))
    ix = {(a.residue.name, a.name): k for k, a in enumerate(top.topology.atoms())}
    system = top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False, removeCMMotion=False)
    ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName('Reference'))

    def fg(v):
        ctx.setPositions(v.reshape(-1, 3) * 0.1)
        st = ctx.getState(getEnergy=True, getForces=True)
        e = st.getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
        g = -st.getForces(asNumpy=True).value_in_unit(u.kilocalorie_per_mole / u.angstrom)
        return e, np.asarray(g).ravel()
    res = minimize(fg, np.asarray(x0, float).ravel(), jac=True, method='L-BFGS-B',
                   options=dict(maxiter=200000, maxfun=400000, ftol=1e-15, gtol=1e-6, maxcor=30))
    x = res.x.reshape(-1, 3)
    e, g = fg(res.x)

    def angles(x):
        zn, n2, s, h, c4 = (x[ix[k]] for k in (('ZN1', 'ZN'), ('SN1', 'N2'), ('SN1', 'S'), ('SN1', 'HN21'), ('SN1', 'C4')))
        z, hh, uu = [(v - n2) / np.linalg.norm(v - n2) for v in (zn, h, s)]
        return dict(Zn_N2_S=round(ang(zn, n2, s), 3), Zn_N2_H=round(ang(zn, n2, h), 3), S_N2_H=round(ang(s, n2, h), 3),
                    signed_volume=round(float(np.dot(z, np.cross(hh, uu))), 4), Zn_N2=round(float(np.linalg.norm(zn - n2)), 4),
                    N2_S=round(float(np.linalg.norm(s - n2)), 4), C4_S_N2=round(ang(c4, s, n2), 3))
    out = angles(x)
    out.update(energy=round(e, 4), largest_force_kcal_mol_A=float(np.abs(g).max()), lbfgs_iterations=int(res.nit))
    ctx.setPositions(x * 0.1)
    openmm.LocalEnergyMinimizer.minimize(ctx, 1e-4, 50000)
    x2 = np.array(ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(u.angstrom))
    out['after_openmm_minimiser'] = angles(x2)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--final', action='store_true', help='read the test set (once; needs N2_FINAL_TEST=1)')
    args = ap.parse_args()
    if args.final:
        raise SystemExit('the final test is a separate program (n2_exam_final.py); this one reads fitting data only')
    rows, notes, split, mm = load(final=False)
    assert not any(r['set'] == 'test' for r in rows)
    by = {r['id']: r for r in rows}
    x_opt, opt_rec = optimum()
    oa = optimum_angles()
    r0 = {k: np.radians(v) for k, v in oa.items()}
    ref = min(rows, key=lambda r: r['e_qm'])
    doc = dict(what='gate G-R2: the conditions that need fitting data only, by the examiner (own force-field route, own fit, '
                    'own leave-one-out, own minimiser); no result file of a test point was opened by this program',
               written=now(), rows=dict(old=sum(r['set'] == 'old' for r in rows), fit=sum(r['set'] == 'fit' for r in rows),
                                        test=sum(r['set'] == 'test' for r in rows), ids=[r['id'] for r in rows]),
               data_notes=notes, lowest_E_QM_of_the_fitting_rows=ref['id'], limits=LIMITS,
               quantum_chemical_optimum=dict(opt_rec, angles_recomputed={k: round(v, 4) for k, v in oa.items()}),
               candidates={}, all_candidates_of_the_ladder={})
    xX = np.array([[float(v) for v in t.split()[1:4]] for t in (GR / 'geom' / 'X.xyz').read_text().splitlines()[2:66]])
    for name, c in CAND.items():
        new, built, rest = c['new'], c['built'], c['rest']
        groups = LADDER[c['ladder']]
        f = fit_and_loo(rows, groups, c['ridge'], rest, r0)
        # is the delivered file that fit?  my refit against the topology built from the file
        d = np.array([r['mm'][new] for r in rows]) - f['prediction']
        declared = float(np.abs(d - d.mean()).max())
        # deviations of the built candidate topology
        dev = np.array([(r['mm'][new] - ref['mm'][new]) - (r['e_qm'] - ref['e_qm']) for r in rows])
        dev0 = np.array([(r['mm'][built] - ref['mm'][built]) - (r['e_qm'] - ref['e_qm']) for r in rows])
        notref = np.array([r is not ref for r in rows])
        loo = f['loo']
        ir = rows.index(ref)
        loo_one = np.delete(loo - loo[ir], ir)
        c3 = part_A(by, new)
        c3b = part_A(by, built)
        # part A from my own leave-one-out errors (what condition 3 would measure out of sample)
        e_loo = {r['id']: v for r, v in zip(rows, loo)}
        la = []
        for ids, r_id in ((FAMILY_X, c3['families']['X']['reference']), (FAMILY_RX, c3['families']['Rx']['reference'])):
            la += [e_loo[i] - e_loo[r_id] for i in ids if i != r_id]
        c4 = bounds_of_topology(new, built)
        mins = {}
        prm = mm['topologies'][new]['prmtop']
        for tag, x0 in (('from_X', xX), ('from_the_optimum', x_opt)):
            m = minimum(prm, x0)
            m['difference'] = dict(Zn_N2_S=round(m['Zn_N2_S'] - oa['ZnN2S'], 3), Zn_N2_H=round(m['Zn_N2_H'] - oa['ZnN2H'], 3),
                                   S_N2_H=round(m['S_N2_H'] - oa['SN2H'], 3))
            m['within_5_deg'] = bool(all(abs(v) <= LIMITS['minimum_deg'] for v in m['difference'].values()))
            mins[tag] = m
        mins_built = {}
        for tag, x0 in (('from_X', xX), ('from_the_optimum', x_opt)):
            m = minimum(mm['topologies'][built]['prmtop'], x0)
            m['difference'] = dict(Zn_N2_S=round(m['Zn_N2_S'] - oa['ZnN2S'], 3), Zn_N2_H=round(m['Zn_N2_H'] - oa['ZnN2H'], 3),
                                   S_N2_H=round(m['S_N2_H'] - oa['SN2H'], 3))
            mins_built[tag] = m
        # the builder's numbers for the same candidate
        bc = json.loads((c['builder_fit'] / 'candidates.json').read_text())
        b1 = [k for k in bc['candidates'] if k['name'] == f"{c['ladder']}_ridge{c['ridge']:g}"][0]
        bp = {r['id']: r for r in csv.DictReader((c['builder_fit'] / 'fit_points.csv').open())}
        d_res = np.array([f['residuals'][k] - float(bp[r['id']]['residual_first_candidate']) for k, r in enumerate(rows)])
        d_loo = np.array([loo[k] - float(bp[r['id']]['leave_one_out_error_first_candidate']) for k, r in enumerate(rows)])
        worst = int(np.argmax(np.abs(loo)))
        w1 = int(np.argmax(np.abs(dev * notref)))
        doc['candidates'][name] = dict(
            frcmod=str(c['frcmod']), frcmod_sha256=sha(c['frcmod']), topology=new, topology_as_built=built,
            model=f"{c['ladder']}, ridge {c['ridge']}", n_parameters=f['n_parameters'], n_fitting_rows=len(rows),
            condition_1_test_rms=dict(evaluated=False, note='needs the test set'),
            condition_2_test_largest=dict(evaluated=False, note='needs the test set'),
            condition_3_part_A=dict(holds=bool(c3['rms'] <= LIMITS['part_A_rms']), limit=LIMITS['part_A_rms'], **c3,
                                    in_sample=True, as_built=dict(rms=c3b['rms'], single_reference_rms=c3b['single_reference']['rms']),
                                    from_my_leave_one_out_errors=round(rms(la), 4)),
            condition_4_bounds=c4,
            condition_5_leave_one_out=dict(
                holds=bool(rms(loo) <= LIMITS['loo_rms']), limit=LIMITS['loo_rms'], leave_one_out_rms=round(rms(loo), 4),
                largest=round(float(loo[worst]), 4), largest_at=rows[worst]['id'], residual_rms=round(rms(f['residuals']), 4),
                refit_solution_on_a_bound=bool(f['bounded']), leave_one_out_refits_on_a_bound=f['loo_refits_on_a_bound'],
                points_per_parameter=round(len(rows) / f['n_parameters'], 2),
                against_one_reference=dict(reference=ref['id'], rms=round(rms(loo_one), 4),
                                           largest=round(float(loo_one[np.argmax(np.abs(loo_one))]), 4)),
                the_delivered_file_is_this_fit=dict(largest_difference_kcal=round(declared, 5), holds=bool(declared <= 0.01)),
                refit_in_amber_numbers=amber_numbers(f['names'], f['p'], r0),
                builder=dict(leave_one_out_rms=b1['leave_one_out_rms'], residual_rms=b1['residual_rms'],
                             largest_difference_of_a_residual=float(np.abs(d_res).max()),
                             largest_difference_of_a_leave_one_out_error=float(np.abs(d_loo).max()))),
            condition_6_minimum=dict(holds=bool(all(m['within_5_deg'] for m in mins.values())), limit_deg=LIMITS['minimum_deg'],
                                     optimum={k: round(v, 3) for k, v in oa.items()}, minima=mins, as_built=mins_built),
            fitting_rows_built_topology=dict(
                rms_against_one_reference=round(rms(dev[notref]), 4), largest=round(float(dev[w1]), 4), largest_at=rows[w1]['id'],
                rms_about_the_mean=round(rms(dev - dev.mean()), 4),
                as_built_rms_against_one_reference=round(rms(dev0[notref]), 4), as_built_rms_about_the_mean=round(rms(dev0 - dev0.mean()), 4),
                note='in-sample for the candidate; NOT condition 1'),
            for_the_record=record(by, new, built))
        # every candidate of the ladder with this charge set: is there one whose leave-one-out error reaches 1.0?
        table = []
        for lad, groups in LADDER.items():
            for ridge in (0.0, 0.01, 0.1):
                g = fit_and_loo(rows, groups, ridge, rest, r0)
                bname = f'{lad}_ridge{ridge:g}'
                bb = [k for k in bc['candidates'] if k['name'] == bname]
                table.append(dict(candidate=bname, n_parameters=g['n_parameters'], residual_rms=round(rms(g['residuals']), 4),
                                  leave_one_out_rms=round(rms(g['loo']), 4), on_a_bound=bool(g['bounded']),
                                  leave_one_out_refits_on_a_bound=g['loo_refits_on_a_bound'],
                                  builder_leave_one_out_rms=bb[0]['leave_one_out_rms'] if bb else None,
                                  builder_residual_rms=bb[0]['residual_rms'] if bb else None))
        doc['all_candidates_of_the_ladder'][name] = dict(
            smallest_leave_one_out_rms=min(t['leave_one_out_rms'] for t in table),
            any_at_most_1=bool(any(t['leave_one_out_rms'] <= 1.0 for t in table)),
            largest_difference_from_the_builder=max(abs(t['leave_one_out_rms'] - t['builder_leave_one_out_rms']) for t in table
                                                    if t['builder_leave_one_out_rms'] is not None), candidates=table)
    doc['sha256'] = dict(this_program=sha(__file__), mm_energies=sha(OUT / 'mm_energies.json'), d3_check=sha(OUT / 'd3_check.json'),
                         split=sha(N2 / 'split.json'))
    import scipy
    doc['versions'] = dict(numpy=np.__version__, scipy=scipy.__version__, python=sys.version.split()[0])
    write_json(OUT / 'conditions_fitting_data.json', doc)
    for name, c in doc['candidates'].items():
        print('==', name)
        print('  c3 part A:', c['condition_3_part_A']['rms'], 'single ref', c['condition_3_part_A']['single_reference']['rms'],
              'as built', c['condition_3_part_A']['as_built'], 'from LOO', c['condition_3_part_A']['from_my_leave_one_out_errors'],
              'refs', c['condition_3_part_A']['families']['X']['reference'], c['condition_3_part_A']['families']['Rx']['reference'])
        print('  c4:', c['condition_4_bounds']['holds'], c['condition_4_bounds']['outside_bounds'], c['condition_4_bounds']['differences_the_plan_does_not_allow'],
              c['condition_4_bounds']['counts_of_differences'])
        c5 = c['condition_5_leave_one_out']
        print('  c5 LOO:', c5['leave_one_out_rms'], 'largest', c5['largest'], c5['largest_at'], 'residual', c5['residual_rms'], 'one ref', c5['against_one_reference'],
              'file is the fit:', c5['the_delivered_file_is_this_fit'], 'builder:', c5['builder'], 'bounded', c5['refit_solution_on_a_bound'], c5['leave_one_out_refits_on_a_bound'])
        print('  c6:', c['condition_6_minimum']['holds'], {k: (v['difference'], round(v['largest_force_kcal_mol_A'], 6), v['after_openmm_minimiser']['Zn_N2_S']) for k, v in c['condition_6_minimum']['minima'].items()})
        print('     as built:', {k: v['difference'] for k, v in c['condition_6_minimum']['as_built'].items()})
        print('  fitting rows, built topology:', c['fitting_rows_built_topology'])
        print('  record:', json.dumps(c['for_the_record']))
    for name, t in doc['all_candidates_of_the_ladder'].items():
        print('ladder', name, 'smallest LOO', t['smallest_leave_one_out_rms'], 'any <= 1:', t['any_at_most_1'], 'largest difference from the builder', t['largest_difference_from_the_builder'])
        print('   ', [(c['candidate'], c['leave_one_out_rms'], c['builder_leave_one_out_rms']) for c in t['candidates']])
    print('data notes:', notes, '| rows', doc['rows']['old'], doc['rows']['fit'], doc['rows']['test'], '| lowest E_QM', ref['id'])


def record(by, new, built):
    """For the record, without weight in the verdict (old points of gate G-R; kcal/mol)."""
    def rel(a, b, key):
        return by[a]['mm'][key] - by[b]['mm'][key]

    def qm(a, b):
        return by[a]['e_qm'] - by[b]['e_qm']
    out = dict(Rx_against_X=dict(QM=round(qm('Rx', 'X'), 3), as_built=round(rel('Rx', 'X', built), 3), new=round(rel('Rx', 'X', new), 3),
                                 mismatch_before=round(rel('Rx', 'X', built) - qm('Rx', 'X'), 3),
                                 mismatch_after=round(rel('Rx', 'X', new) - qm('Rx', 'X'), 3)))
    path = ['a2', 'a3', 'a4', 'a5', 'a6']
    sv = {p: by[p]['mm']['internal']['signed_volume'] for p in path}
    flat = min(path, key=lambda p: abs(sv[p]))
    top = max(path, key=lambda p: qm(p, 'X'))
    out['inversion_path_against_X'] = dict(
        points={p: dict(signed_volume=round(sv[p], 3), QM=round(qm(p, 'X'), 3), as_built=round(rel(p, 'X', built), 3),
                        new=round(rel(p, 'X', new), 3)) for p in path},
        flattest_point=flat, highest_QM_point=top,
        cost_of_the_planar_nitrogen=dict(point=flat, QM=round(qm(flat, 'X'), 3), as_built=round(rel(flat, 'X', built), 3),
                                         new=round(rel(flat, 'X', new), 3)))
    dpts = ['d1', 'd2', 'd3', 'd4', 'd5', 'd6']
    dflat = min(dpts, key=lambda p: abs(by[p]['mm']['internal']['signed_volume']))
    out['hydrogen_moved_into_the_plane'] = dict(point=dflat, signed_volume=round(by[dflat]['mm']['internal']['signed_volume'], 3),
                                                QM=round(qm(dflat, 'X'), 3), as_built=round(rel(dflat, 'X', built), 3),
                                                new=round(rel(dflat, 'X', new), 3))
    out['turn_about_Zn_N2_against_X'] = {p: dict(QM=round(qm(p, 'X'), 3), as_built=round(rel(p, 'X', built), 3),
                                                new=round(rel(p, 'X', new), 3)) for p in ('e_p015', 'e_p030', 'e_m120', 'e_p120')}
    return out


if __name__ == '__main__':
    main()
