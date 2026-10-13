#!/usr/bin/env python
"""Gate G-R2, for the record (no weight in the verdict): the force field as built and the two candidates at RELAXED
geometries of the 64-atom model.  None of these geometries is a fitting or a test point of the gate.

    AMBERHOME=/root/cqc-sfa-tools/md nice -n 10 /root/cqc-sfa-tools/md/bin/python md/qm/n2/final_test/n2_exam_behaviour.py

 1  Along the two quantum-chemical relaxations of gate G-R (ligand relaxed, zinc and imidazoles fixed; files of the
    gate read from the Pi's pulled copy /root/cqc-sfa-md/qm/gr/relax/, not changed): E_QM with D3(BJ) and the
    force-field energy of each of the four topologies of n2_exam_mm.py, all against ONE reference, the last step of
    the run from X.  The cost of a planar nitrogen at a relaxed geometry: the step of the run from Rx with the
    smallest |signed volume|.
 2  The examiner's own relaxation with zinc and imidazoles fixed, from X and from Rx, for each topology: a descent
    that only goes downhill, in boxes of +/- 0.02 A per coordinate (L-BFGS-B inside each box), until the largest
    force on a ligand atom is below 0.001 kcal/mol/A.  Does the R-like placement end in a minimum of its own (signed volume at N2 still
    negative)?  At every end point the Hessian of the 18 ligand atoms (central differences of forces) says
    whether it is a minimum.

Writes md/qm/n2/final_test/behaviour.json.
"""
import hashlib
import json
import pathlib
import time

import numpy as np

REPO = pathlib.Path('/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate')
N2 = REPO / 'md' / 'qm' / 'n2'
GR = REPO / 'md' / 'qm' / 'gr'
OUT = N2 / 'final_test'
RELAX = pathlib.Path('/root/cqc-sfa-md/qm/gr/relax')
HARTREE = 627.5094740631
LABEL = dict(T0='terms and charges as built', T1='candidate, charges as built', T0b='charges of two placements, terms as built',
             T2='candidate of the second step')


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def now():
    return time.strftime('%Y-%m-%d %I:%M:%S %p IST', time.gmtime(time.time() + 19800))


def ang(a, b, c):
    u, v = a - b, c - b
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))


def dih(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1n = b1 / np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1n) * b1n, b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


class Top:
    def __init__(self, prmtop):
        import openmm
        from openmm import app, unit as u
        self.u = u
        top = app.AmberPrmtopFile(str(prmtop))
        self.ix = {(a.residue.name, a.name): k for k, a in enumerate(top.topology.atoms())}
        self.lig = np.array([k for k, a in enumerate(top.topology.atoms()) if a.residue.name == 'SN1'])
        system = top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False, removeCMMotion=False)
        self.ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName('Reference'))

    def ef(self, x):
        self.ctx.setPositions(np.asarray(x) * 0.1)
        st = self.ctx.getState(getEnergy=True, getForces=True)
        return (st.getPotentialEnergy().value_in_unit(self.u.kilocalorie_per_mole),
                np.asarray(st.getForces(asNumpy=True).value_in_unit(self.u.kilocalorie_per_mole / self.u.angstrom)))

    def geo(self, x):
        zn, n2, s, h, c4 = (x[self.ix[k]] for k in (('ZN1', 'ZN'), ('SN1', 'N2'), ('SN1', 'S'), ('SN1', 'HN21'), ('SN1', 'C4')))
        z, hh, uu = [(v - n2) / np.linalg.norm(v - n2) for v in (zn, h, s)]
        return dict(signed_volume=round(float(np.dot(z, np.cross(hh, uu))), 4), Zn_N2_S=round(ang(zn, n2, s), 2), Zn_N2_H=round(ang(zn, n2, h), 2),
                    S_N2_H=round(ang(s, n2, h), 2), Zn_N2_S_C4=round(dih(zn, n2, s, c4), 1))

    def descend(self, x0, box=0.02, fmax=1e-3, limit=20000):
        """Only downhill, ligand atoms only, in small boxes: from the current point the energy is minimised by
        L-BFGS-B (a descent method: every accepted point is lower than the one before) with every ligand coordinate
        confined to +/- box A about the current point; the result is the next point.  So no coordinate moves more than
        `box` between two points of the path and the energy never rises between them.  Ends when the largest force on
        a ligand atom is below fmax (kcal/mol/A)."""
        from scipy.optimize import minimize
        x = np.array(x0, float)
        lig = self.lig

        def fg(v):
            y = x.copy()
            y[lig] = v.reshape(-1, 3)
            e, f = self.ef(y)
            return e, -f[lig].ravel()
        e, f = self.ef(x)
        n, crossed_at, e0 = 0, None, e
        sv0 = self.geo(x)['signed_volume']
        evals = 0
        while n < limit and float(np.abs(f[lig]).max()) >= fmax:
            v0 = x[lig].ravel()
            res = minimize(fg, v0, jac=True, method='L-BFGS-B', bounds=list(zip(v0 - box, v0 + box)),
                           options=dict(maxiter=200, ftol=0.0, gtol=fmax * 0.3, maxcor=20))
            evals += res.nfev
            if not res.fun < e - 1e-12:
                break
            x[lig] = res.x.reshape(-1, 3)
            e, f = self.ef(x)
            n += 1
            if crossed_at is None and np.sign(self.geo(x)['signed_volume']) != np.sign(sv0):
                crossed_at = dict(step=n, energy_below_the_start=round(e0 - e, 3))
        return x, dict(steps=n, energy_evaluations=int(evals), energy=e, energy_below_the_start=round(e0 - e, 3),
                       largest_force=float(np.abs(f[lig]).max()), converged=bool(np.abs(f[lig]).max() < fmax), sign_changed=crossed_at,
                       **self.geo(x))

    def hessian_lowest(self, x, h=1e-3):
        n = len(self.lig) * 3
        H = np.zeros((n, n))
        for k in range(n):
            a, c = divmod(k, 3)
            xp, xm = x.copy(), x.copy()
            xp[self.lig[a], c] += h
            xm[self.lig[a], c] -= h
            H[k] = -(self.ef(xp)[1][self.lig].ravel() - self.ef(xm)[1][self.lig].ravel()) / (2 * h)
        w = np.linalg.eigvalsh(0.5 * (H + H.T))
        return [round(float(v), 4) for v in w[:3]]


def read_xyz(path):
    l = pathlib.Path(path).read_text().splitlines()
    return np.array([[float(v) for v in t.split()[1:4]] for t in l[2:2 + int(l[0])]])


def main():
    mm = json.loads((OUT / 'mm_energies.json').read_text())
    tops = {t: Top(mm['topologies'][t]['prmtop']) for t in ('T0', 'T1', 'T0b', 'T2')}
    doc = dict(what='force field as built and candidates at relaxed geometries of the 64-atom model; for the record of gate G-R2, no '
                    'weight in the verdict; none of these geometries is a fitting or a test point', written=now(), topologies=LABEL,
               files={f.name: sha(f) for f in sorted(RELAX.glob('*')) if f.suffix in ('.npz', '.jsonl')})
    # ---- 1. along the quantum-chemical relaxations
    paths = {}
    for start in ('X', 'Rx'):
        xyz = np.load(RELAX / f'{start}_traj.npz')['xyz']
        steps = [json.loads(l) for l in (RELAX / f'{start}_steps.jsonl').read_text().splitlines() if l.strip()]
        n = min(len(xyz), len(steps))
        assert np.abs(xyz[0] - read_xyz(GR / 'geom' / f'{start}.xyz')).max() < 2e-3
        paths[start] = (xyz[:n], steps[:n])
    xr, sr = paths['X'][0][-1], paths['X'][1][-1]
    e_ref = {t: tops[t].ef(xr)[0] for t in tops}
    rel = {}
    for start, (xyz, steps) in paths.items():
        rows = []
        for x, s in zip(xyz, steps):
            r = dict(step=s['step'], signed_volume=s['signed_volume'], QM=round((s['e_total'] - sr['e_total']) * HARTREE, 3))
            for t in tops:
                r[t] = round(tops[t].ef(x)[0] - e_ref[t], 3)
            rows.append(r)
        d = {t: np.array([r[t] - r['QM'] for r in rows]) for t in tops}
        rel[start] = dict(n_steps=len(rows), last_step=steps[-1]['step'], last_time=steps[-1]['time'],
                          last=dict(signed_volume=steps[-1]['signed_volume'], grad_rms_ligand=steps[-1]['grad_rms_ligand'],
                                    Zn_N2_S=steps[-1].get('Zn_N2_S'), Zn_N2_H=steps[-1].get('Zn_N2_H'), S_N2_H=steps[-1].get('S_N2_H'),
                                    rmsd_ligand_heavy_from_X_run=steps[-1].get('rmsd_ligand_heavy_from_X_run'),
                                    e_total_minus_X_run_kcal_mol=steps[-1].get('e_total_minus_X_run_kcal_mol')),
                          rms_force_field_minus_QM={t: round(float(np.sqrt(np.mean(v ** 2))), 3) for t, v in d.items()},
                          rms_about_the_mean_of_the_path={t: round(float(np.sqrt(np.mean((v - v.mean()) ** 2))), 3) for t, v in d.items()},
                          rows=rows)
    rows = rel['Rx']['rows']
    flat = min(rows, key=lambda r: abs(r['signed_volume']))
    doc['along_the_quantum_chemical_relaxations'] = dict(
        reference='the last step of the run from X', runs=rel,
        flattest_nitrogen_of_the_run_from_Rx=dict(step=flat['step'], signed_volume=flat['signed_volume'], QM=flat['QM'],
                                                  force_field={t: flat[t] for t in tops},
                                                  force_field_minus_QM={t: round(flat[t] - flat['QM'], 3) for t in tops}),
        first_step_of_the_run_from_Rx=dict(QM=rows[0]['QM'], force_field={t: rows[0][t] for t in tops}),
        last_step_of_the_run_from_Rx=dict(QM=rows[-1]['QM'], force_field={t: rows[-1][t] for t in tops}))
    # ---- 2. own relaxations with zinc and imidazoles fixed
    xX, xR = read_xyz(GR / 'geom' / 'X.xyz'), read_xyz(GR / 'geom' / 'Rx.xyz')
    own = {}
    for t in tops:
        ex, rx_ = tops[t].descend(xX)
        er, rr = tops[t].descend(xR)
        d = np.sqrt(np.mean(np.sum((ex[tops[t].lig] - er[tops[t].lig]) ** 2, axis=1)))
        own[t] = dict(from_X=dict(rx_, lowest_eigenvalues_of_the_ligand_hessian=tops[t].hessian_lowest(ex)),
                      from_Rx=dict(rr, lowest_eigenvalues_of_the_ligand_hessian=tops[t].hessian_lowest(er)),
                      end_from_Rx_minus_end_from_X_kcal=round(rr['energy'] - rx_['energy'], 3),
                      rms_distance_of_the_two_end_points_A=round(float(d), 3),
                      R_like_placement_has_a_minimum_of_its_own=bool(rr['converged'] and rr['signed_volume'] < 0))
        print(t, LABEL[t], '| from X:', {k: rx_[k] for k in ('steps', 'signed_volume', 'Zn_N2_S', 'Zn_N2_H', 'S_N2_H', 'converged')},
              '| from Rx:', {k: rr[k] for k in ('steps', 'signed_volume', 'Zn_N2_S', 'converged', 'sign_changed')},
              '| Rx end - X end', own[t]['end_from_Rx_minus_end_from_X_kcal'], '| own minimum:', own[t]['R_like_placement_has_a_minimum_of_its_own'],
              '| hessian', own[t]['from_Rx']['lowest_eigenvalues_of_the_ligand_hessian'])
    doc['own_relaxations_zinc_and_imidazoles_fixed'] = own
    doc['this_program'] = sha(__file__)
    (OUT / 'behaviour.json').write_text(json.dumps(doc, indent=1) + '\n')
    a = doc['along_the_quantum_chemical_relaxations']
    print('QM run from Rx: last', a['runs']['Rx']['last_step'], a['runs']['Rx']['last_time'], a['runs']['Rx']['last'])
    print('flattest nitrogen:', a['flattest_nitrogen_of_the_run_from_Rx'])
    print('rms along X:', a['runs']['X']['rms_force_field_minus_QM'], '| along Rx:', a['runs']['Rx']['rms_force_field_minus_QM'])


if __name__ == '__main__':
    main()
