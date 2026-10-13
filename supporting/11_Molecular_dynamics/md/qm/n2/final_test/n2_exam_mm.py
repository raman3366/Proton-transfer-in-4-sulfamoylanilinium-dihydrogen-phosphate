#!/usr/bin/env python
"""Gate G-R2, the examiner's own force-field route (md/N2.md, section 9).  Geometry only: NO quantum-chemical energy is
read here, so the program may be run on fitting and test geometries alike.

    AMBERHOME=/root/cqc-sfa-tools/md nice -n 10 /root/cqc-sfa-tools/md/bin/python md/qm/n2/final_test/n2_exam_mm.py

Four topologies of the 64-atom model [Zn(ethyl-imidazole)3(sulfonamidate)]+, each built HERE with tleap from the
inputs kept in md/qm/n2/ff/inputs (charges as built) or md/qm/n2/ff2/inputs (charges of the declared second step):
    T0   charges as built,        terms at N2 as built
    T1   charges as built,        md/qm/n2/n2_terms_new_charges_as_built.frcmod loaded last
    T0b  charges of two placements, terms at N2 as built
    T2   charges of two placements, md/qm/n2/n2_terms_new.frcmod loaded last        ("the new terms")
Energies: OpenMM, Reference platform, no cut-off, at the 36 rigid geometries of gate G-R and the 192 new geometries.

Two routes to the energy of a candidate, compared with each other and with the builder's numbers:
    A  the topology built with the candidate's frcmod (T1, T2);
    B  the topology as built (T0, T0b) minus its own terms at N2 plus the candidate's terms, the terms computed here
       in NumPy from the numbers of the frcmod file (angles K (th - th0)^2, torsions PK [1 + cos(n phi - gamma)],
       improper PK [1 + cos(n w - gamma)] on the atom order tleap wrote).
Also: every parameter of T1 and T2 against T0 and T0b (what differs, and whether anything outside N2 differs), a
few energies by sander (another engine), the as-built energies against the table of gate G-R, and D3(BJ) of every
geometry by the gate's own routine (imported from md/qm/gr_queue.py when run with --d3 in the quantum-chemistry
environment) and by the library s-dftd3 (another implementation).

Writes md/qm/n2/final_test/mm_energies.json, topology_differences.json and (with --d3) d3_check.json.
"""
import argparse
import hashlib
import json
import os
import pathlib
import shutil
import subprocess
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
AMBERHOME = os.environ.get('AMBERHOME', '/root/cqc-sfa-tools/md')
KCAL = 4.184
HARTREE = 627.5094740631
TOPS = dict(T0=('ff', None), T1=('ff', N2 / 'n2_terms_new_charges_as_built.frcmod'),
            T0b=('ff2', None), T2=('ff2', N2 / 'n2_terms_new.frcmod'))


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def now():
    return time.strftime('%Y-%m-%d %I:%M:%S %p IST', time.gmtime(time.time() + 19800))


def read_xyz(path):
    l = pathlib.Path(path).read_text().splitlines()
    n = int(l[0])
    return [t.split()[0] for t in l[2:2 + n]], np.array([[float(v) for v in t.split()[1:4]] for t in l[2:2 + n]])


def write_json(path, obj):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(obj, indent=1) + '\n')
    os.replace(tmp, path)


# ------------------------------------------------------------------------------------------------ geometry
def ang(a, b, c):
    u, v = a - b, c - b
    return float(np.arccos(np.clip(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1.0, 1.0)))


def dih(a, b, c, d):
    """IUPAC dihedral in radians (the sign convention of AMBER)."""
    b0, b1, b2 = a - b, c - b, d - c
    b1n = b1 / np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1n) * b1n, b2 - np.dot(b2, b1n) * b1n
    return float(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


# ------------------------------------------------------------------------------------------------ topologies
def build(tag):
    """tleap in scratch: the inputs of the force-field copy, the candidate's frcmod loaded last."""
    copy, frc = TOPS[tag]
    work = SCR / 'tops' / tag
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    for f in (N2 / copy / 'inputs').iterdir():
        shutil.copyfile(f, work / f.name)
    lines = []
    for l in (work / 'leap.in').read_text().splitlines():
        if l.startswith('mol = loadpdb') and frc is not None:
            shutil.copyfile(frc, work / 'candidate.frcmod')
            lines.append('loadamberparams candidate.frcmod')
        lines.append(l)
    (work / 'leap_exam.in').write_text('\n'.join(lines) + '\n')
    env = dict(os.environ, AMBERHOME=AMBERHOME, PATH=f'{AMBERHOME}/bin:' + os.environ.get('PATH', ''))
    p = subprocess.run(['tleap', '-s', '-f', 'leap_exam.in'], cwd=work, env=env, capture_output=True, text=True)
    (work / 'leap_exam.log').write_text(p.stdout + p.stderr)
    log = p.stdout + p.stderr
    bad = [l for l in log.splitlines() if 'Could not find' in l or 'FATAL' in l or 'Could not open' in l]
    if bad or not (work / 'small.prmtop').exists():
        raise SystemExit(f'tleap failed for {tag}: {bad[:5]}')
    warn = [l for l in log.splitlines() if 'arning' in l]
    return work / 'small.prmtop', dict(inputs=str(N2 / copy / 'inputs'), frcmod=str(frc) if frc else None,
                                       frcmod_sha256=sha(frc) if frc else None, prmtop=str(work / 'small.prmtop'),
                                       prmtop_sha256=sha(work / 'small.prmtop'), tleap_warnings=len(warn))


class Top:
    def __init__(self, prmtop):
        import openmm
        from openmm import app, unit as u
        self.mm, self.u = openmm, u
        top = app.AmberPrmtopFile(str(prmtop))
        self.atoms = [(a.residue.name, a.name, a.element.symbol) for a in top.topology.atoms()]
        self.ix = {(r, n): k for k, (r, n, _) in enumerate(self.atoms)}
        self.system = top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False,
                                       removeCMMotion=False)
        self.ctx = openmm.Context(self.system, openmm.VerletIntegrator(0.001),
                                  openmm.Platform.getPlatformByName('Reference'))

    def energy(self, xyz):
        self.ctx.setPositions(np.asarray(xyz) * 0.1)
        return self.ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(self.u.kilocalorie_per_mole)

    def parameters(self):
        """Every parameter of the system, keyed by atom indices, in AMBER units."""
        u = self.u
        out = dict(bonds={}, angles={}, torsions={}, charges=[], lj=[], exceptions={})
        for f in self.system.getForces():
            t = type(f).__name__
            if t == 'HarmonicBondForce':
                for k in range(f.getNumBonds()):
                    a, b, r, kk = f.getBondParameters(k)
                    out['bonds'][tuple(sorted((a, b)))] = (round(r.value_in_unit(u.angstrom), 6),
                                                           round(kk.value_in_unit(u.kilocalorie_per_mole / u.angstrom ** 2) / 2, 5))
            elif t == 'HarmonicAngleForce':
                for k in range(f.getNumAngles()):
                    a, b, c, th, kk = f.getAngleParameters(k)
                    key = (a, b, c) if a < c else (c, b, a)
                    out['angles'][key] = (round(kk.value_in_unit(u.kilocalorie_per_mole / u.radian ** 2) / 2, 5),
                                          round(th.value_in_unit(u.degree), 5))
            elif t == 'PeriodicTorsionForce':
                for k in range(f.getNumTorsions()):
                    a, b, c, d, per, ph, kk = f.getTorsionParameters(k)
                    out['torsions'].setdefault((a, b, c, d), []).append(
                        (int(per), round(kk.value_in_unit(u.kilocalorie_per_mole), 6), round(ph.value_in_unit(u.degree) % 360.0, 4)))
            elif t == 'NonbondedForce':
                for k in range(f.getNumParticles()):
                    q, sg, ep = f.getParticleParameters(k)
                    out['charges'].append(round(q.value_in_unit(u.elementary_charge), 6))
                    out['lj'].append((round(sg.value_in_unit(u.angstrom), 6), round(ep.value_in_unit(u.kilocalorie_per_mole), 6)))
                for k in range(f.getNumExceptions()):
                    a, b, qq, sg, ep = f.getExceptionParameters(k)
                    out['exceptions'][tuple(sorted((a, b)))] = (round(qq.value_in_unit(u.elementary_charge ** 2), 6),
                                                                round(sg.value_in_unit(u.angstrom), 6),
                                                                round(ep.value_in_unit(u.kilocalorie_per_mole), 6))
        for v in out['torsions'].values():
            v.sort()
        return out


def n2_terms_of_topology(top):
    """The terms the plan allows to change, read from a built topology: angles with N2 in the middle; every torsion
    whose four atoms hold N2 and whose middle bond is N2-S or Zn-N2 (proper), or whose atoms are {Zn, N2, S, HN21}
    (the improper).  Returned as lists of (kind, atoms, numbers)."""
    ix = top.ix
    zn, n2, s, h = ix[('ZN1', 'ZN')], ix[('SN1', 'N2')], ix[('SN1', 'S')], ix[('SN1', 'HN21')]
    par = top.parameters()
    terms = []
    for (a, b, c), (k, th0) in par['angles'].items():
        if b == n2:
            terms.append(('angle', (a, b, c), (k, th0)))
    for at, lst in par['torsions'].items():
        a, b, c, d = at
        if {b, c} in ({n2, s}, {n2, zn}) or set(at) == {zn, n2, s, h}:
            for per, k, ph in lst:
                terms.append(('improper' if set(at) == {zn, n2, s, h} else 'torsion', at, (per, k, ph)))
    return terms


def terms_energy(terms, x):
    e = 0.0
    for kind, at, v in terms:
        if kind == 'angle':
            e += v[0] * (ang(x[at[0]], x[at[1]], x[at[2]]) - np.radians(v[1])) ** 2
        else:
            per, k, ph = v
            e += k * (1.0 + np.cos(per * dih(x[at[0]], x[at[1]], x[at[2]], x[at[3]]) - np.radians(ph)))
    return float(e)


def parse_frcmod(path, top_new):
    """The terms of a candidate frcmod as (kind, atoms, numbers) on the atoms of the 64-atom model.  Atom types of
    MCPB.py: M1 zinc, Y1 Y2 Y3 the histidine nitrogens (HD1 NE2, HD2 NE2, HE1 ND1), Y4 N2, sy S, hn HN21, ca C4.
    A proper torsion type applies to every quadruple of these types that is bonded in a row; the improper to the
    quadruple tleap wrote (its atom order is read from the topology `top_new`, its numbers from the file)."""
    ix = top_new.ix
    name = dict(M1=[ix[('ZN1', 'ZN')]], Y1=[ix[('HD1', 'NE2')]], Y2=[ix[('HD2', 'NE2')]], Y3=[ix[('HE1', 'ND1')]],
                Y4=[ix[('SN1', 'N2')]], sy=[ix[('SN1', 'S')]], hn=[ix[('SN1', 'HN21')]], ca=[ix[('SN1', 'C4')]])
    sec, terms = None, []
    for l in pathlib.Path(path).read_text().splitlines()[1:]:
        t = l.strip()
        if t in ('MASS', 'BOND', 'ANGLE', 'DIHE', 'IMPROPER', 'NONBON'):
            sec = t
            continue
        if not t:
            continue
        if sec == 'ANGLE':
            ty = [v.strip() for v in l[:8].split('-')]
            k, th0 = float(l[8:].split()[0]), float(l[8:].split()[1])
            terms.append(('angle', tuple(name[v][0] for v in ty), (k, th0)))
        elif sec == 'DIHE':
            ty = [v.strip() for v in l[:11].split('-')]
            f = l[11:].split()
            idivf, pk, ph, per = float(f[0]), float(f[1]), float(f[2]), abs(int(float(f[3])))
            terms.append(('torsion', tuple(name[v][0] for v in ty), (per, pk / idivf, ph % 360.0)))
        elif sec == 'IMPROPER':
            ty = [v.strip() for v in l[:11].split('-')]
            f = l[11:].split()
            pk, ph, per = float(f[0]), float(f[1]), abs(int(float(f[2])))
            want = {name[v][0] for v in ty}
            order = [at for kind, at, _ in n2_terms_of_topology(top_new) if kind == 'improper' and set(at) == want]
            assert len(order) == 1, 'the improper of the file is not once in the topology'
            terms.append(('improper', order[0], (per, pk, ph % 360.0)))
        elif sec in ('MASS', 'BOND', 'NONBON'):
            raise SystemExit(f'{path}: a line in section {sec}: {l!r} (the plan allows none)')
    return terms


def compare(par0, par1, top):
    """What differs between two parameter sets of the same atoms; names instead of indices."""
    nm = lambda k: '-'.join(f'{top.atoms[i][0]}:{top.atoms[i][1]}' for i in k)   # noqa: E731
    d = dict(bonds=[], angles=[], torsions=[], charges=[], lj=[], exceptions=[])
    for sec in ('bonds', 'angles', 'torsions', 'exceptions'):
        for k in sorted(set(par0[sec]) | set(par1[sec])):
            if par0[sec].get(k) != par1[sec].get(k):
                d[sec].append(dict(atoms=nm(k), before=par0[sec].get(k), after=par1[sec].get(k)))
    for sec in ('charges', 'lj'):
        for i, (a, b) in enumerate(zip(par0[sec], par1[sec])):
            if a != b:
                d[sec].append(dict(atom=nm((i,)), before=a, after=b))
    return d


def sander_energy(prmtop, xyz, work):
    """Single-point energy by sander (vacuum, no cut-off): another engine than OpenMM."""
    work.mkdir(parents=True, exist_ok=True)
    n = len(xyz)
    flat = [f'{v:12.7f}' for v in np.asarray(xyz).ravel()]
    rows = [''.join(flat[i:i + 6]) for i in range(0, len(flat), 6)]
    (work / 'in.rst7').write_text('exam\n' + f'{n:5d}\n' + '\n'.join(rows) + '\n')
    (work / 'sp.in').write_text('single point\n&cntrl\n imin=1, maxcyc=0, ncyc=0, ntb=0, igb=0, cut=999.0, ntpr=1,\n/\n')
    env = dict(os.environ, AMBERHOME=AMBERHOME, PATH=f'{AMBERHOME}/bin:' + os.environ.get('PATH', ''))
    subprocess.run(['sander', '-O', '-i', 'sp.in', '-p', str(prmtop), '-c', 'in.rst7', '-o', 'sp.out', '-r', 'sp.rst7',
                    '-inf', 'sp.info'], cwd=work, env=env, capture_output=True, text=True)
    txt = (work / 'sp.out').read_text().splitlines()
    for i, l in enumerate(txt):
        if l.strip().startswith('NSTEP') and 'ENERGY' in l:
            return float(txt[i + 1].split()[1])
    raise SystemExit(f'sander gave no energy in {work}')


def d3_check():
    """D3(BJ), two-body, of every geometry: the gate's routine (md/qm/gr_queue.py, class D3, imported) against the
    library s-dftd3 with the B3LYP parameters and s9 = 0.  Geometry only."""
    sys.path.insert(0, str(QM))
    import gr_queue                                                        # imported, not edited
    from dftd3.interface import DispersionModel, RationalDampingParam
    import dftd3
    split = json.loads((N2 / 'split.json').read_text())
    ids = [(p, GR / 'geom' / f'{p}.xyz') for p in split['old_rigid_points']] + \
          [(p, N2 / 'geom' / f'{p}.xyz') for p in split['order']]
    import gaucom                                                          # the bohr of the gate
    ref = json.loads((GR / 'd3_reference.json').read_text())
    how = 'gr_queue.D3(z, md/qm/gr/d3_reference.json).energy(xyz / gaucom.BOHR)'
    out = {}
    for pid, f in ids:
        sym, x = read_xyz(f)
        z = [gr_queue.D3_Z[s] for s in sym]
        e_gate = gr_queue.D3(z, ref).energy(x / gaucom.BOHR)
        model = DispersionModel(np.array(z), x / gaucom.BOHR)
        e_lib = float(model.get_dispersion(RationalDampingParam(s6=1.0, s8=1.9889, a1=0.3981, a2=4.4211, s9=0.0),
                                           grad=False)['energy'])
        out[pid] = dict(gate_routine=float(e_gate), s_dftd3=e_lib)
    dd = np.array([v['gate_routine'] - v['s_dftd3'] for v in out.values()])
    print('gate routine minus s-dftd3, largest (Hartree):', float(np.abs(dd).max()))
    write_json(OUT / 'd3_check.json', dict(
        what='D3(BJ) two-body of the 36 + 192 geometries of gate G-R2: the routine of gate G-R against the library s-dftd3 '
             '(B3LYP: s6 1.0, s8 1.9889, a1 0.3981, a2 4.4211, s9 0); Hartree; geometry only',
        written=now(), gr_queue_sha256=sha(QM / 'gr_queue.py'), gate_routine_call=how,
        dftd3_version=getattr(dftd3, '__version__', None), numpy=np.__version__,
        d3_reference_sha256=sha(GR / 'd3_reference.json'), b3lyp_bj_of_the_gate=ref['b3lyp_bj'],
        largest_abs_difference_hartree=float(np.abs(dd).max()), points=out))
    print('d3_check.json written:', len(out), 'geometries;', how)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--d3', action='store_true', help='only the D3 check (quantum-chemistry environment)')
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.d3:
        return d3_check()
    import openmm
    split = json.loads((N2 / 'split.json').read_text())
    assert sha(N2 / 'split.json') == (N2 / 'split.json.sha256').read_text().split()[0]
    geoms = {p: GR / 'geom' / f'{p}.xyz' for p in split['old_rigid_points']}
    geoms.update({p: N2 / 'geom' / f'{p}.xyz' for p in split['order']})
    for p in split['order']:
        assert sha(geoms[p]) == split['xyz_sha256'][p], p
    rec, tops = {}, {}
    for tag in TOPS:
        prm, rec[tag] = build(tag)
        tops[tag] = Top(prm)
    sym, _ = read_xyz(geoms['X'])
    assert [s.capitalize() for s in sym] == [a[2] for a in tops['T0'].atoms], 'atom order of the geometry files'
    # the kept topologies of the builder against the ones built here
    kept = {}
    for tag, copy in (('T0', 'ff'), ('T0b', 'ff2')):
        k = Top(N2 / copy / 'small_as_built.prmtop')
        kept[copy] = dict(parameters_equal=k.parameters() == tops[tag].parameters(),
                          kept_sha256=sha(N2 / copy / 'small_as_built.prmtop'))
    # what differs
    par = {t: tops[t].parameters() for t in tops}
    diff = {'T1_against_T0': compare(par['T0'], par['T1'], tops['T0']),
            'T2_against_T0b': compare(par['T0b'], par['T2'], tops['T0']),
            'T0b_against_T0': compare(par['T0'], par['T0b'], tops['T0']),
            'T2_against_T0': compare(par['T0'], par['T2'], tops['T0'])}
    counts = {k: {s: len(v) for s, v in d.items()} for k, d in diff.items()}
    write_json(OUT / 'topology_differences.json', dict(
        what='parameters of the 64-atom topologies built by the examiner, candidate against as built (OpenMM system: '
             'bonds, angles, torsions with impropers, charges, Lennard-Jones, 1-4 exceptions)', written=now(),
        topologies=rec, kept_topologies_of_the_builder=kept, counts=counts, differences=diff,
        terms_at_N2={t: [dict(kind=k, atoms='-'.join(tops[t].atoms[i][1] for i in at), numbers=v)
                         for k, at, v in n2_terms_of_topology(tops[t])] for t in tops}))
    # energies, route A and route B
    asb = {t: n2_terms_of_topology(tops[t]) for t in ('T0', 'T0b')}
    new = {'T1': parse_frcmod(TOPS['T1'][1], tops['T1']), 'T2': parse_frcmod(TOPS['T2'][1], tops['T2'])}
    pts = {}
    for pid, f in geoms.items():
        _, x = read_xyz(f)
        e = {t: tops[t].energy(x) for t in tops}
        e['T1_route_B'] = e['T0'] - terms_energy(asb['T0'], x) + terms_energy(new['T1'], x)
        e['T2_route_B'] = e['T0b'] - terms_energy(asb['T0b'], x) + terms_energy(new['T2'], x)
        e['rest_T0'] = e['T0'] - terms_energy(asb['T0'], x)
        e['rest_T0b'] = e['T0b'] - terms_energy(asb['T0b'], x)
        zn, n2, s, h, c4 = (x[tops['T0'].ix[k]] for k in (('ZN1', 'ZN'), ('SN1', 'N2'), ('SN1', 'S'), ('SN1', 'HN21'), ('SN1', 'C4')))
        his = [x[tops['T0'].ix[k]] for k in (('HD1', 'NE2'), ('HD2', 'NE2'), ('HE1', 'ND1'))]
        imp = [at for kind, at, _ in n2_terms_of_topology(tops['T1']) if kind == 'improper'][0]
        z, hh, uu = [(v - n2) / np.linalg.norm(v - n2) for v in (zn, h, s)]
        e['internal'] = dict(th_ZnN2S=ang(zn, n2, s), th_ZnN2H=ang(zn, n2, h), th_SN2H=ang(s, n2, h),
                             phi_HNSC=dih(h, n2, s, c4), phi_ZnNSC=dih(zn, n2, s, c4),
                             phi_His_S=[dih(p, zn, n2, s) for p in his], phi_His_H=[dih(p, zn, n2, h) for p in his],
                             omega=dih(*(x[i] for i in imp)), signed_volume=float(np.dot(z, np.cross(hh, uu))))
        pts[pid] = e
    spread = {}
    for a, b in (('T1', 'T1_route_B'), ('T2', 'T2_route_B')):
        d = np.array([pts[p][a] - pts[p][b] for p in pts])
        spread[f'{a}_minus_route_B'] = dict(mean=float(d.mean()), largest_about_the_mean=float(np.abs(d - d.mean()).max()),
                                            largest_abs=float(np.abs(d).max()))
    # sander on a few geometries (fitting and old points; geometry only)
    sand = {}
    for tag in ('T0', 'T1', 'T2'):
        for pid in ('X', 'Rx', 'a4', 'e_p030', 'n001', 'n050', 'f02'):
            es = sander_energy(SCR / 'tops' / tag / 'small.prmtop', read_xyz(geoms[pid])[1], SCR / 'sander' / tag)
            sand[f'{tag}:{pid}'] = dict(sander=es, openmm=pts[pid][tag], difference=es - pts[pid][tag])
    # the as-built energies against the table of gate G-R (another program, md/qm/gr_mm.py)
    table = {}
    tf = GR / 'gr_table.csv'
    if tf.exists():
        import csv
        rows = list(csv.DictReader(tf.open()))
        d = np.array([(pts[r['id']]['T0'] - pts['X']['T0']) - float(r['E_MM']) for r in rows if r['id'] in pts])
        table = dict(column='E_MM (relative to X, md/qm/gr_mm.py)', n=len(d), largest_abs_difference=float(np.abs(d).max()))
    write_json(OUT / 'mm_energies.json', dict(
        what='force-field energies (kcal/mol) of the 36 + 192 geometries of gate G-R2 by the examiner: OpenMM '
             f'{openmm.__version__}, Reference platform, no cut-off; four topologies built here with tleap; geometry only, no '
             'quantum-chemical energy was read', written=now(), topologies=rec,
        route_A_against_route_B=spread, sander_against_openmm=sand, as_built_against_the_table_of_gate_G_R=table,
        kept_topologies_of_the_builder=kept, new_terms_as_parsed={k: [dict(kind=a, atoms='-'.join(tops['T0'].atoms[i][1] for i in at), numbers=v)
                                                                    for a, at, v in t] for k, t in new.items()},
        sha256=dict(split=sha(N2 / 'split.json'), this_program=sha(__file__)),
        versions=dict(openmm=openmm.__version__, numpy=np.__version__, python=sys.version.split()[0]),
        points=pts))
    print('topologies:', {t: r['prmtop_sha256'][:12] for t, r in rec.items()})
    print('kept topologies equal to the ones built here:', kept)
    print('differences (counts):', json.dumps(counts))
    print('route A against route B:', json.dumps(spread))
    print('sander against OpenMM, largest:', max(abs(v['difference']) for v in sand.values()))
    print('gate G-R table:', table)


if __name__ == '__main__':
    main()
