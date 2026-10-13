"""Gate G-R2: a FIRST LOOK at how the small model behaves with a candidate set of terms at N2.  A TEST, cheap, on
the 64-atom model in vacuum: not the gate, not a result for the paper (md/N2.md, section 6.7).

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/n2_behaviour.py --model MODEL.json
        [--out DIR] [--descent-steps 300000]

For the topology as built and for the topology with the terms of MODEL.json (built here with tleap from the
inputs of gate G-R, as md/qm/n2_gate.py does):

 1  MINIMUM, ALL ATOMS FREE, from start X and from the quantum-chemical optimum of the small model
    (md/qm/n2/qm_optimum.json): the three angles at N2 against that optimum.  This is what condition 6 of the gate
    will compute; here it is a look, the gate decides.
 2  RELAXATION WITH ZINC AND IMIDAZOLES FIXED (atoms 1 to 46), from X and from the R-like placement Rx:
    (a) OpenMM's minimiser (L-BFGS), which may step over a barrier;
    (b) from Rx a descent that only accepts a lower energy and moves no atom by more than 0.02 A in a step, with
        the signed volume at N2 recorded: does the R-like placement have a minimum of its own, or does it invert?
        The descent ends when the largest force is below 0.001 kcal/mol/A, when the nitrogen has gone over to the
        configuration of X (signed volume above +0.3), or at the step limit; the minimiser is then started at its
        end to say where that leads.
 3  THE FORCE FIELD AT THE GEOMETRIES OF THE QUANTUM-CHEMICAL RELAXATIONS as far as they have got (files of gate
    G-R, read and not changed: md/qm/gr/relax/qm_<start>.json and the geometry beside it, or the newer copies the
    watchdog brought from the machine): E(from Rx) - E(from X) in both descriptions.
 4  THE CONSTRUCTED INVERSION PATH a1 (= X) to a7 (= Rx) of gate G-R: quantum chemistry, as built, candidate.

Writes <out>/behaviour.json and the end geometries as .xyz; <out> is md/qm/n2/behaviour/<name of the model file>.
"""
import argparse
import json
import pathlib
import sys

import numpy as np

import n2_common as nc

N_FIXED = 46
WATCHDOG_GR = pathlib.Path('/root/cqc-sfa-md/from_gpu/repo_md/qm/gr')


class Small:
    """The small model in OpenMM: all atoms free, and with the first 46 atoms fixed (mass zero)."""

    def __init__(self, prmtop):
        import openmm
        from openmm import app, unit as u
        self.mm, self.u = openmm, u
        top = app.AmberPrmtopFile(str(prmtop))
        self.names = [(a.residue.name, a.name) for a in top.topology.atoms()]
        self.el = [a.element.symbol for a in top.topology.atoms()]
        ix = {n: k for k, n in enumerate(self.names)}
        self.ix = {k: ix[v] for k, v in dict(zn=('ZN1', 'ZN'), n2=('SN1', 'N2'), s=('SN1', 'S'), h=('SN1', 'HN21'),
                                             c4=('SN1', 'C4'), n94=('HD1', 'NE2')).items()}
        self.ctx = {}
        for which in ('free', 'fixed'):
            sysm = top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False, removeCMMotion=False)
            if which == 'fixed':
                for k in range(N_FIXED):
                    sysm.setParticleMass(k, 0.0)
            self.ctx[which] = openmm.Context(sysm, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName('Reference'))

    def energy(self, x, which='free', forces=False):
        c = self.ctx[which]
        c.setPositions(np.asarray(x) * 0.1)
        st = c.getState(getEnergy=True, getForces=forces)
        e = st.getPotentialEnergy().value_in_unit(self.u.kilocalorie_per_mole)
        if not forces:
            return e
        f = np.array(st.getForces(asNumpy=True).value_in_unit(self.u.kilocalorie_per_mole / self.u.angstrom))
        return e, f

    def minimise(self, x, which):
        c = self.ctx[which]
        c.setPositions(np.asarray(x) * 0.1)
        self.mm.LocalEnergyMinimizer.minimize(c, 1e-4, 50000)
        st = c.getState(getEnergy=True, getPositions=True)
        return (np.array(st.getPositions(asNumpy=True).value_in_unit(self.u.angstrom)),
                st.getPotentialEnergy().value_in_unit(self.u.kilocalorie_per_mole))

    def describe(self, x):
        i = self.ix
        zn, n2, s, h, c4, n94 = (x[i[k]] for k in ('zn', 'n2', 's', 'h', 'c4', 'n94'))
        z, hh, uu = (nc.unit(v - n2) for v in (zn, h, s))
        a = dict(Zn_N2_S=nc.angle(zn, n2, s), Zn_N2_H=nc.angle(zn, n2, h), S_N2_H=nc.angle(s, n2, h))
        return dict({k: round(v, 2) for k, v in a.items()}, angle_sum=round(sum(a.values()), 2), C4_S_N2=round(nc.angle(c4, s, n2), 2),
                    signed_volume=round(float(np.dot(z, np.cross(hh, uu))), 4),
                    Zn_N2=round(float(np.linalg.norm(zn - n2)), 4), N2_S=round(float(np.linalg.norm(s - n2)), 4),
                    turn_about_Zn_N2=round(nc.dihedral(n94, zn, n2, s), 2), Zn_N2_S_C4=round(nc.dihedral(zn, n2, s, c4), 2))

    def descent(self, x0, max_steps, step=0.02):
        """Downhill only, no atom moved by more than `step` A; fixed atoms stay.  Returns the end and the record."""
        x = np.array(x0, float)
        e, f = self.energy(x, 'fixed', True)
        f[:N_FIXED] = 0.0
        e0, s, sv0 = e, step, self.describe(x)['signed_volume']
        crossed_at, n = None, 0
        trace = [dict(step=0, e=round(e - e0, 4), signed_volume=sv0)]
        sv = sv0
        while n < max_steps and s > 1e-6:
            fmax = np.linalg.norm(f, axis=1).max()
            if fmax < 1e-3 or (sv0 < 0 and sv > 0.3):
                break
            xn = x + f * (s / fmax)
            en, fn = self.energy(xn, 'fixed', True)
            if en < e:
                x, e, f = xn, en, fn
                f[:N_FIXED] = 0.0
                n += 1
                s = min(step, s * 1.2)
                sv = self.describe(x)['signed_volume']
                if crossed_at is None and sv * sv0 < 0:
                    crossed_at = dict(step=n, energy_below_the_start=round(e0 - e, 3))
                if n % 2000 == 0:
                    trace.append(dict(step=n, e=round(e - e0, 4), signed_volume=sv))
            else:
                s *= 0.5
        fmax = float(np.linalg.norm(f, axis=1).max())
        trace.append(dict(step=n, e=round(e - e0, 4), signed_volume=self.describe(x)['signed_volume']))
        return x, dict(steps=n, energy_start=round(e0, 3), energy_end=round(e, 3), fell_by=round(e0 - e, 3),
                       largest_force_at_the_end_kcal_per_A=round(fmax, 4),
                       ended=('largest force below 0.001 kcal/mol/A' if fmax < 1e-3 else
                              'went over to the configuration of X (signed volume above +0.3)' if (sv0 < 0 and sv > 0.3) else
                              'step limit' if n >= max_steps else 'no lower energy within 1e-6 A'),
                       sign_changed=crossed_at, trace=trace)


def heavy_rmsd(el, a, b):
    k = [i for i in range(N_FIXED, len(el)) if el[i] != 'H']
    return round(float(np.sqrt(((a[k] - b[k]) ** 2).sum(1).mean())), 3)


def write_xyz(path, el, x, comment):
    pathlib.Path(path).write_text('\n'.join([str(len(el)), comment] + [f'{e:2s} {v[0]:14.8f} {v[1]:14.8f} {v[2]:14.8f}'
                                                                      for e, v in zip(el, x)]) + '\n')


def qm_relaxations():
    """State of the quantum-chemical relaxations of gate G-R: the newer of the repository copy and the watchdog copy."""
    out = {}
    for start in ('X', 'Rx'):
        best = None
        for root in (nc.GR, WATCHDOG_GR):
            f = root / 'relax' / f'qm_{start}.json'
            if not f.exists():
                continue
            r = json.loads(f.read_text())
            steps = r.get('steps') or []
            if best is None or len(steps) > len(best[1].get('steps') or []):
                best = (root, r)
        if best is None:
            continue
        root, r = best
        steps = r.get('steps') or []
        last = steps[-1] if steps else {}
        geo = None
        for name in (f'qm_{start}.xyz', f'qm_{start}_last.xyz'):
            if (root / 'relax' / name).exists():
                geo = root / 'relax' / name
                break
        out[start] = dict(file=str(root / 'relax' / f'qm_{start}.json'), sha256=nc.sha(root / 'relax' / f'qm_{start}.json'),
                          state=r.get('state'), converged=r.get('converged'), crossed=r.get('crossed'),
                          n_steps=len(steps), last_step={k: last.get(k) for k in last if not isinstance(last.get(k), (list, dict))},
                          geometry=str(geo) if geo else None, geometry_sha256=nc.sha(geo) if geo else None)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', required=True)
    ap.add_argument('--out', default=None)
    ap.add_argument('--descent-steps', type=int, default=300000)
    args = ap.parse_args()
    model = json.loads(pathlib.Path(args.model).read_text())
    sfx = '' if nc.FF_NAME == 'ff' else f'_{nc.FF_NAME}'
    out = pathlib.Path(args.out) if args.out else nc.N2 / f'behaviour{sfx}' / pathlib.Path(args.model).stem
    if model.get('force_field', 'ff') != nc.FF_NAME:
        raise SystemExit(f"the model was fitted with the force-field copy {model.get('force_field')!r}, this run uses {nc.FF_NAME!r} (N2_FF)")
    out.mkdir(parents=True, exist_ok=True)
    frc = out / 'n2_terms_candidate.frcmod'
    nc.write_frcmod(model['parameters'], frc, f"terms at N2, candidate {model.get('name', '')}; written by md/qm/n2_behaviour.py for a look, NOT ACCEPTED")
    tops = dict(as_built=nc.FFDIR / 'small_as_built.prmtop', candidate=nc.build_topology(frc, out / 'work'))
    rec_opt, _ = nc.qm_optimum()
    x_opt = np.array(rec_opt['xyz_angstrom'])
    _, x_X = nc.read_xyz(nc.GR / 'geom' / 'X.xyz')
    _, x_Rx = nc.read_xyz(nc.GR / 'geom' / 'Rx.xyz')
    qm = qm_relaxations()
    doc = dict(what='first look at the behaviour of the small model with a candidate set of terms at N2 (md/N2.md, section 6.7); '
                    'A TEST on 64 atoms in vacuum, not the gate and not a result; written by md/qm/n2_behaviour.py',
               written=nc.now(), model=dict(file=str(args.model), sha256=nc.sha(args.model), name=model.get('name')),
               force_field_copy=nc.FF_NAME,
               quantum_chemical_optimum=dict(angles_at_N2=rec_opt['angles_at_N2'], signed_volume=rec_opt['signed_volume'],
                                             C4_S_N2=rec_opt['C4_S_N2'], N2_S=rec_opt['N2_S'], Zn_N2=rec_opt['Zn_N2'],
                                             Zn_N2_S_C4=rec_opt['Zn_N2_S_C4'],
                                             source=rec_opt['source'], note='B3LYP/6-31G* without a dispersion term, all atoms free'),
               quantum_chemical_relaxations_zinc_and_imidazoles_fixed=qm, force_fields={})
    for tag, top in tops.items():
        sm = Small(top)
        r = dict(prmtop=str(top), prmtop_sha256=nc.sha(top))
        # 1 minimum with all atoms free
        r['minimum_all_atoms_free'] = {}
        for start, x0 in (('from_X', x_X), ('from_the_optimum', x_opt)):
            x, e = sm.minimise(x0, 'free')
            d = sm.describe(x)
            diff = {k: round(d[k] - rec_opt['angles_at_N2'][k], 2) for k in ('Zn_N2_S', 'Zn_N2_H', 'S_N2_H')}
            bc = dict(C4_S_N2=round(d['C4_S_N2'] - rec_opt['C4_S_N2'], 2), N2_S=round(d['N2_S'] - rec_opt['N2_S'], 4),
                      Zn_N2_S=diff['Zn_N2_S'])
            r['minimum_all_atoms_free'][start] = dict(d, energy=round(e, 3), minus_the_optimum=diff,
                                                      within_5_deg=bool(all(abs(v) <= 5.0 for v in diff.values())),
                                                      binding_check_quantities_minus_the_optimum=dict(
                                                          bc, tolerances=dict(C4_S_N2=2.0, N2_S=0.02, Zn_N2_S=3.0),
                                                          within=bool(abs(bc['C4_S_N2']) <= 2.0 and abs(bc['N2_S']) <= 0.02 and abs(bc['Zn_N2_S']) <= 3.0),
                                                          note='reported; by the plan (6 Oct, 1:09 PM, point 3) the ligand terms are not re-tuned to this optimum'))
            write_xyz(out / f'{tag}_free_{start}.xyz', sm.el, x, f'{tag}, all atoms free, minimised {start}; a test')
        # 2 relaxations with zinc and imidazoles fixed
        rel = {}
        ends = {}
        for start, x0 in (('X', x_X), ('Rx', x_Rx)):
            x, e = sm.minimise(x0, 'fixed')
            ends[start] = x
            rel[f'minimiser_from_{start}'] = dict(start=sm.describe(x0), end=sm.describe(x), energy_start=round(sm.energy(x0), 3),
                                                  energy_end=round(e, 3), ligand_heavy_atoms_moved_A=heavy_rmsd(sm.el, x, x0),
                                                  largest_shift_of_a_fixed_atom_A=round(float(np.abs(x[:N_FIXED] - x0[:N_FIXED]).max()), 6))
            write_xyz(out / f'{tag}_fixed_minimiser_from_{start}.xyz', sm.el, x, f'{tag}, atoms 1-46 fixed, minimiser from {start}; a test')
        xd, dr = sm.descent(x_Rx, args.descent_steps)
        xdm, edm = sm.minimise(xd, 'fixed')
        dr.update(start=sm.describe(x_Rx), end=sm.describe(xd), ligand_heavy_atoms_moved_A=heavy_rmsd(sm.el, xd, x_Rx),
                  ligand_heavy_atoms_from_the_minimum_reached_from_X_A=heavy_rmsd(sm.el, xd, ends['X']),
                  energy_above_the_minimum_reached_from_X=round(dr['energy_end'] - rel['minimiser_from_X']['energy_end'], 3),
                  minimiser_started_at_the_end_of_the_descent=dict(end=sm.describe(xdm), energy=round(edm, 3),
                                                                   ligand_heavy_atoms_moved_A=heavy_rmsd(sm.el, xdm, xd)))
        write_xyz(out / f'{tag}_fixed_descent_from_Rx.xyz', sm.el, xd, f'{tag}, atoms 1-46 fixed, downhill-only descent from Rx; a test')
        rel['descent_from_Rx'] = dr
        sv_x, sv_r = rel['minimiser_from_X']['end']['signed_volume'], dr['end']['signed_volume']
        sv_rm = dr['minimiser_started_at_the_end_of_the_descent']['end']['signed_volume']
        rel['reading'] = dict(
            R_like_placement_has_a_minimum_of_its_own=(None if dr['ended'] == 'step limit' else bool(
                dr['ended'].startswith('largest force') and dr['end']['signed_volume'] < -0.3)),
            descent_from_Rx_ends_with_a_nitrogen=('pyramidal, R-like' if sv_r < -0.3 else 'pyramidal, as in X' if sv_r > 0.3
                                                  else 'planar or nearly planar'),
            minimiser_from_the_end_of_the_descent_ends_with_a_nitrogen=(
                'pyramidal, R-like' if sv_rm < -0.3 else 'pyramidal, as in X' if sv_rm > 0.3 else 'planar or nearly planar'),
            descent_from_Rx_changed_sign=dr['sign_changed'] is not None,
            two_relaxed_configurations=bool(sv_x > 0.3 and sv_rm < -0.3),
            E_relaxed_from_Rx_minus_from_X=round(dr['minimiser_started_at_the_end_of_the_descent']['energy']
                                                 - rel['minimiser_from_X']['energy_end'], 3),
            note='signed volume at N2: above +0.3 the configuration of X, below -0.3 the R-like one, between: planar or nearly')
        r['zinc_and_imidazoles_fixed'] = rel
        # 3 the force field at the quantum-chemical geometries
        at = {}
        for start, q in qm.items():
            if q['geometry']:
                _, xq = nc.read_xyz(q['geometry'])
                at[start] = dict(energy=round(sm.energy(xq), 3), **sm.describe(xq),
                                 largest_shift_of_a_fixed_atom_A=round(float(np.abs(xq[:N_FIXED] - x_X[:N_FIXED]).max()), 4),
                                 ligand_heavy_atoms_from_the_force_field_minimum_reached_from_X_A=heavy_rmsd(sm.el, xq, ends['X']),
                                 ligand_heavy_atoms_from_the_force_field_end_reached_from_Rx_A=heavy_rmsd(sm.el, xq, xdm),
                                 energy_above_the_force_field_minimum_reached_from_X=round(sm.energy(xq) - rel['minimiser_from_X']['energy_end'], 3))
        if len(at) == 2:
            at['from_Rx_minus_from_X'] = round(at['Rx']['energy'] - at['X']['energy'], 3)
        r['at_the_geometries_of_the_quantum_chemical_relaxations'] = at
        # 4 the constructed inversion path
        path = {}
        for pid in ('X', 'a2', 'a3', 'a4', 'a5', 'a6', 'Rx', 'd1', 'g1'):
            _, xp = nc.read_xyz(nc.GR / 'geom' / f'{pid}.xyz')
            path[pid] = round(sm.energy(xp) - sm.energy(x_X), 3)
        r['constructed_points_against_X'] = path
        doc['force_fields'][tag] = r
    # quantum chemistry on the same constructed points (the 38 points of gate G-R are not test data)
    d3 = json.loads((nc.GR / 'd3.json').read_text())['points']
    eq = {}
    for pid in ('X', 'a2', 'a3', 'a4', 'a5', 'a6', 'Rx', 'd1', 'g1'):
        q = json.loads((nc.GR / 'qm' / f'{pid}.json').read_text())
        eq[pid] = (q['e_tot'] + d3[pid]['e_d3bj']) * nc.HARTREE
    doc['quantum_chemistry_constructed_points_against_X'] = {k: round(v - eq['X'], 3) for k, v in eq.items()}
    if len(qm) == 2 and all('e_total' in q['last_step'] for q in qm.values()):
        doc['quantum_chemistry_from_Rx_minus_from_X_at_the_last_steps'] = dict(
            kcal=round((qm['Rx']['last_step']['e_total'] - qm['X']['last_step']['e_total']) * nc.HARTREE, 3),
            with_D3=True, a_result=bool(all(q['state'] == 'done' for q in qm.values())),
            note='the distance between two runs on their way unless both have converged')
    nc.write_json(out / 'behaviour.json', doc)
    for tag, r in doc['force_fields'].items():
        print('==', tag)
        for k, v in r['minimum_all_atoms_free'].items():
            print('  free minimum', k, {a: v[a] for a in ('Zn_N2_S', 'Zn_N2_H', 'S_N2_H', 'signed_volume')}, 'minus optimum', v['minus_the_optimum'], v['within_5_deg'])
        z = r['zinc_and_imidazoles_fixed']
        for k in ('minimiser_from_X', 'minimiser_from_Rx'):
            print('  fixed', k, z[k]['energy_start'], '->', z[k]['energy_end'], 'signed volume', z[k]['start']['signed_volume'], '->', z[k]['end']['signed_volume'])
        d = z['descent_from_Rx']
        print('  descent from Rx:', d['steps'], 'steps,', d['ended'], ', fell by', d['fell_by'], ', sign changed:', d['sign_changed'],
              ', end', {a: d['end'][a] for a in ('Zn_N2_S', 'Zn_N2_H', 'S_N2_H', 'signed_volume', 'turn_about_Zn_N2', 'Zn_N2_S_C4')})
        print('  reading:', z['reading'])
        print('  at the QM relaxation geometries:', {k: (v if not isinstance(v, dict) else v['energy']) for k, v in r['at_the_geometries_of_the_quantum_chemical_relaxations'].items()})
        print('  constructed points against X:', r['constructed_points_against_X'])
    print('== quantum chemistry, constructed points against X:', doc['quantum_chemistry_constructed_points_against_X'])
    print('== quantum-chemical relaxations:', {k: (v['state'], v['n_steps'], v['last_step'].get('signed_volume'),
                                                   v['last_step'].get('grad_rms_ligand')) for k, v in qm.items()},
          doc.get('quantum_chemistry_from_Rx_minus_from_X_at_the_last_steps'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
