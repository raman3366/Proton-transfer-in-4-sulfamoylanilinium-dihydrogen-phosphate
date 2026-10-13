"""Gate G-R2 (md/PLAN.md, entry of 6 Oct 2026, 1:34 PM IST): what the programs n2_*.py share.  md/N2.md.

Nothing here computes a quantum-chemical energy.  Four things:

 1. THE SIX COORDINATES of the ligand in the rigid small model (zinc and imidazoles fixed, ligand body B rigid,
    Zn-N2 1.931 A, N2-H 1.021 A), `coords()`, and the construction of a geometry from them, `construct()`:
      w     unit vector Zn -> N2 (the place of N2 on its sphere about the zinc); as a design coordinate it is
            the point `s` on the great circle through its directions in X (s = 0) and in Rx (s = 1), with `t`
            the elevation off that circle in degrees (0 in every designed point);
      th1   angle Zn-N2-S;                chi   torsion NE2(His94)-Zn-N2-S (the turn about Zn-N2);
      tau   torsion Zn-N2-S-C4 (the turn about N2-S);
      th2   angle Zn-N2-H;                eta   azimuth of N2->H about the axis N2->Zn, counted from the
            azimuth of N2->S, right-handed about N2->Zn, in [0, 360): 238 in X, 180 for a planar nitrogen,
            about 122 in the R-like configuration.  Signed volume at N2 = -sin(th1) sin(th2) sin(eta).
    `construct(coords(c))` gives c back for every rigid point of gate G-R (checked by n2_points.py).

 2. THE FORCE FIELD AS BUILT on the 64 atoms, `FF`: OpenMM, Reference platform, no cut-off, from the topology
    md/qm/n2/ff/small_as_built.prmtop (a copy of the topology of gate G-R, cap scheme "methyl", checksum in
    md/qm/n2/ff/SHA256SUMS).  `energy()` returns the total and `rest`: the total without the terms that the plan
    allows to change (the three angles at N2; every torsion about N2-S; every torsion about Zn-N2).

 3. THE MODEL CLASS, `FEATURES` and `features()`: the allowed terms as functions that are LINEAR in the free
    parameters (md/N2.md, section 2), and `write_frcmod()`, which turns a parameter vector into AMBER lines for
    the atom types of MCPB.py; `build_topology()` makes the small model with such a file loaded last.

 4. THE SPLIT, `load_points()`: the only door to the quantum-chemical energies of the new points.  With
    purpose='fit' it never opens the result file of a point flagged "test" in md/qm/n2/split.json and stops if
    asked for one; purpose='gate' is for n2_gate.py alone and needs the environment variable N2_FINAL_TEST=1,
    and it appends a line to md/qm/n2/test_set_access.log every time it is used.

Copied or imported from the code of gate G-R (none of it edited): md/qm/gr_points.py is IMPORTED (class Model
and the vector helpers; sha256 658c03ae9de62f39f6e8e1952da6b83a0cf80d6fd75ca97672aac0ae85b14847 when this was
written; `GR_POINTS_SHA` below is checked at import and a change is reported, not ignored).
"""
import hashlib
import json
import os
import pathlib
import subprocess
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[1]
N2 = HERE / 'n2'
GR = HERE / 'gr'
# The force-field copy of the small model: 'ff' = as built (the charges of gate G-R).  The declared second step of
# the plan (charges from two placements, md/qm/n2_charges.py) has its own folder, chosen with N2_FF=ff2 in the
# environment; every program then reads its tleap inputs, its topology and its energies from there.
FF_NAME = os.environ.get('N2_FF', 'ff')
FFDIR = N2 / FF_NAME
AMBERHOME = os.environ.get('AMBERHOME', '/root/cqc-sfa-tools/md')
HARTREE = 627.5094740631
KCAL = 4.184
GR_POINTS_SHA = '658c03ae9de62f39f6e8e1952da6b83a0cf80d6fd75ca97672aac0ae85b14847'
SEED = 20261006                                                            # every random choice of gate G-R2

sys.path.insert(0, str(HERE))
import gr_points as gp                                                     # noqa: E402

unit, angle, dihedral, rot, wrap = gp.unit, gp.angle, gp.dihedral, gp.rot, gp.wrap


def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def gr_points_unchanged():
    return sha(HERE / 'gr_points.py') == GR_POINTS_SHA


def now():
    """IST (UTC+5:30), 12-hour, whatever the time zone of the machine."""
    return time.strftime('%Y-%m-%d %I:%M:%S %p IST', time.gmtime(time.time() + 19800))


def write_json(path, obj):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f'.tmp{os.getpid()}')
    tmp.write_text(json.dumps(obj, indent=1) + '\n')
    os.replace(tmp, path)


def read_xyz(path):
    l = pathlib.Path(path).read_text().splitlines()
    n = int(l[0])
    return [t.split()[0] for t in l[2:2 + n]], np.array([[float(v) for v in t.split()[1:4]] for t in l[2:2 + n]])


# ------------------------------------------------------------------------------------- 1. the six coordinates
class Geometry:
    """The rigid small model of gate G-R (gr_points.Model) with the six coordinates of gate G-R2."""

    def __init__(self):
        self.m = m = gp.Model()
        pts = json.loads((GR / 'points.json').read_text())
        lig = {p['id']: np.array([p['ligand_coordinates'][n] for n in m.names]) for p in pts['points']}
        self.old = lig
        self.old_points = {p['id']: p for p in pts['points']}
        self.wx = unit(lig['X'][m.i['N2']] - m.zn)
        self.wr = unit(lig['Rx'][m.i['N2']] - m.zn)
        self.arc = float(np.arccos(np.clip(np.dot(self.wx, self.wr), -1, 1)))
        self.nrm = unit(np.cross(self.wx, self.wr))
        self.his = {}
        for res, at, key in ((93, 'NE2', 'N94'), (95, 'NE2', 'N96'), (118, 'ND1', 'N119')):
            self.his[key] = [a['xyz'] for a in m.atoms if a['resi'] == res and a['name'] == at][0]

    def w_of(self, s, t=0.0):
        """Direction Zn -> N2 at the place s on the great circle X -> Rx, elevation t (degrees) off it."""
        a = s * self.arc
        on = np.cos(a) * self.wx + np.sin(a) * np.cross(self.nrm, self.wx)
        return unit(np.cos(np.radians(t)) * on + np.sin(np.radians(t)) * self.nrm)

    def coords(self, c):
        m, i, zn = self.m, self.m.i, self.m.zn
        n2, s, hn, c4 = c[i['N2']], c[i['S']], c[i['HN21']], c[i['C4']]
        w = unit(n2 - zn)
        z, u, h = -w, unit(s - n2), unit(hn - n2)
        a = unit(u - np.dot(u, z) * z)
        b = np.cross(z, a)
        hp = h - np.dot(h, z) * z
        eta = float(np.degrees(np.arctan2(np.dot(hp, b), np.dot(hp, a)))) % 360.0
        t = float(np.degrees(np.arcsin(np.clip(np.dot(w, self.nrm), -1, 1))))
        wp = unit(w - np.dot(w, self.nrm) * self.nrm)
        spar = float(np.arctan2(np.dot(np.cross(self.wx, wp), self.nrm), np.dot(self.wx, wp)) / self.arc)
        out = dict(s=spar, t=t, th1=angle(zn, n2, s), chi=dihedral(m.n94, zn, n2, s), tau=dihedral(zn, n2, s, c4),
                   th2=angle(zn, n2, hn), eta=eta, th3=angle(s, n2, hn),
                   signed_volume=float(np.dot(z, np.cross(h, u))),
                   Zn_N2=float(np.linalg.norm(n2 - zn)), N2_H=float(np.linalg.norm(hn - n2)),
                   N2_S=float(np.linalg.norm(s - n2)),
                   H_N2_S_C4=dihedral(hn, n2, s, c4), omega_improper=dihedral(zn, hn, n2, s))
        out['angle_sum'] = out['th1'] + out['th2'] + out['th3']
        for key, p in self.his.items():
            out[f'{key}_Zn_N2'] = angle(p, zn, n2)
            out[f'{key}_Zn_N2_S'] = dihedral(p, zn, n2, s)
            out[f'{key}_Zn_N2_H'] = dihedral(p, zn, n2, hn)
        return out

    def construct(self, th1, chi, tau, th2, eta, s=0.0, t=0.0, w=None):
        """The 18 ligand atoms from the six coordinates (w given directly, or by s and t)."""
        m, i, zn = self.m, self.m.i, self.m.zn
        w = self.w_of(s, t) if w is None else unit(np.asarray(w, float))
        n2 = zn + m.d_zn_n * w
        z = -w
        r = m.n94 - zn
        p0 = unit(r - np.dot(r, w) * w)
        q0 = np.cross(w, p0)
        t1 = np.radians(th1)
        u = None
        for sg in (1.0, -1.0):
            al = np.radians(sg * chi)
            cand = np.cos(t1) * z + np.sin(t1) * (np.cos(al) * p0 + np.sin(al) * q0)
            if abs(wrap(dihedral(m.n94, zn, n2, n2 + cand) - chi)) < 1e-6:
                u = cand
                break
        assert u is not None, 'the turn about Zn-N2 could not be set'
        b0 = m.x0 - m.x0[i['N2']]
        u0 = unit(b0[i['S']])
        ax = np.cross(u0, u)
        r1 = np.eye(3) if np.linalg.norm(ax) < 1e-12 else rot(ax, np.degrees(np.arctan2(np.linalg.norm(ax), np.dot(u0, u))))
        c = (r1 @ b0.T).T + n2
        c = m.turn_ns(c, wrap(tau - dihedral(zn, n2, c[i['S']], c[i['C4']])))
        a = unit(u - np.dot(u, z) * z)
        b = np.cross(z, a)
        t2, e = np.radians(th2), np.radians(eta)
        c[i['HN21']] = n2 + gp.NH * (np.cos(t2) * z + np.sin(t2) * (np.cos(e) * a + np.sin(e) * b))
        return c

    def full(self, c):
        """64 x 3: the fixed atoms and the ligand, in the order of the quantum-chemical input."""
        return np.vstack([self.m.fxyz, c])

    def xyz_text(self, c, comment):
        sym = [a['el'] for a in self.m.atoms]
        rows = ['64', comment] + [f'{s:2s} {v[0]:14.8f} {v[1]:14.8f} {v[2]:14.8f}' for s, v in zip(sym, self.full(c))]
        return '\n'.join(rows) + '\n'


# --------------------------------------------------------------------------------- 2. the force field as built
class FF:
    """Energies of the 64 atoms from an AMBER topology; `rest` = without the terms the plan allows to change."""

    def __init__(self, prmtop=None):
        import openmm
        from openmm import app, unit as u
        self.mm, self.u = openmm, u
        self.prmtop = pathlib.Path(prmtop or FFDIR / 'small_as_built.prmtop')
        top = app.AmberPrmtopFile(str(self.prmtop))
        self.names = [(a.residue.name, a.name) for a in top.topology.atoms()]
        ix = {n: k for k, n in enumerate(self.names)}
        self.zn, self.n2, self.s = ix[('ZN1', 'ZN')], ix[('SN1', 'N2')], ix[('SN1', 'S')]
        self.h = ix[('SN1', 'HN21')]
        self.ctx, self.removed = {}, {}
        for which in ('total', 'rest'):
            sysm = top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False,
                                    removeCMMotion=False)
            rem = []
            if which == 'rest':
                for f in sysm.getForces():
                    if type(f).__name__ == 'HarmonicAngleForce':
                        for k in range(f.getNumAngles()):
                            a, b, c, th, kk = f.getAngleParameters(k)
                            if b == self.n2:
                                rem.append(('angle', self.names[a][1], self.names[b][1], self.names[c][1],
                                            kk.value_in_unit(u.kilojoule_per_mole / u.radian ** 2) / 2 / KCAL,
                                            th.value_in_unit(u.degree)))
                                f.setAngleParameters(k, a, b, c, th, 0.0)
                    if type(f).__name__ == 'PeriodicTorsionForce':
                        for k in range(f.getNumTorsions()):
                            a, b, c, d, per, ph, kk = f.getTorsionParameters(k)
                            if {b, c} in ({self.n2, self.s}, {self.n2, self.zn}) or \
                                    {a, b, c, d} == {self.zn, self.n2, self.s, self.h}:
                                kv = kk.value_in_unit(u.kilojoule_per_mole) / KCAL
                                if kv != 0.0:
                                    rem.append(('torsion', '-'.join(self.names[x][1] for x in (a, b, c, d)), int(per),
                                                ph.value_in_unit(u.degree), kv))
                                f.setTorsionParameters(k, a, b, c, d, per, ph, 0.0)
            self.removed[which] = rem
            self.ctx[which] = openmm.Context(sysm, openmm.VerletIntegrator(0.001),
                                             openmm.Platform.getPlatformByName('Reference'))

    def energy(self, xyz64):
        out = {}
        for which, ctx in self.ctx.items():
            ctx.setPositions(np.asarray(xyz64) * 0.1)
            out[which] = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(self.u.kilocalorie_per_mole)
        return out


# ---------------------------------------------------------------------------------------- 3. the model class
# One entry per free parameter: (name, group).  E_model = rest + sum_j p_j f_j + c0.
FEATURES = (
    ('K_ZnN2S', 'angle'), ('B_ZnN2S', 'angle'), ('K_ZnN2H', 'angle'), ('B_ZnN2H', 'angle'),
    ('K_SN2H', 'angle'), ('B_SN2H', 'angle'),
    ('H_NS_C4_1', 'tor_H'), ('H_NS_C4_2', 'tor_H'), ('H_NS_C4_3', 'tor_H'),
    ('Zn_NS_C4_1', 'tor_Zn'), ('Zn_NS_C4_2', 'tor_Zn'), ('Zn_NS_C4_3', 'tor_Zn'),
    ('His_ZnN_S_3c', 'tor_ZnN_S'), ('His_ZnN_S_3s', 'tor_ZnN_S'),
    ('His_ZnN_H_3c', 'tor_ZnN_H'), ('His_ZnN_H_3s', 'tor_ZnN_H'),
    ('improper', 'improper'),
)
# The ladder of models, smallest first (md/N2.md, section 2.4); a model is a list of groups.
LADDER = (
    ('L1', ('angle', 'tor_H')),
    ('L2', ('angle', 'tor_H', 'tor_Zn')),
    ('L3', ('angle', 'tor_H', 'tor_Zn', 'tor_ZnN_S')),
    ('L4', ('angle', 'tor_H', 'tor_Zn', 'tor_ZnN_S', 'improper')),
    ('L5', ('angle', 'tor_H', 'tor_Zn', 'tor_ZnN_S', 'improper', 'tor_ZnN_H')),
)
BOUNDS = dict(k_angle=(10.0, 300.0), theta0_window_deg=15.0, torsion_abs_max=5.0, improper_abs_max=5.0)


def features(cd):
    """The 17 functions at one geometry, from the dictionary of Geometry.coords(); angles in radians."""
    r = np.radians
    f = {}
    for key, name in (('th1', 'ZnN2S'), ('th2', 'ZnN2H'), ('th3', 'SN2H')):
        f[f'K_{name}'] = r(cd[key]) ** 2
        f[f'B_{name}'] = r(cd[key])
    for n in (1, 2, 3):
        f[f'H_NS_C4_{n}'] = np.cos(n * r(cd['H_N2_S_C4']))
        f[f'Zn_NS_C4_{n}'] = np.cos(n * r(cd['tau']))
    for end in ('S', 'H'):
        ph = [r(cd[f'{h}_Zn_N2_{end}']) for h in ('N94', 'N96', 'N119')]
        f[f'His_ZnN_{end}_3c'] = float(sum(np.cos(3 * p) for p in ph))
        f[f'His_ZnN_{end}_3s'] = float(sum(np.sin(3 * p) for p in ph))
    f['improper'] = 1.0 - np.cos(2 * r(cd['omega_improper']))
    return {k: float(f[k]) for k, _ in FEATURES}


def as_built_parameters(ff):
    """The terms of the topology `ff` that `rest` leaves out, in the parameters of the model class (for the
    check that rest + sum p f reproduces the total up to a constant).  As built: Zn-N2-S 50.0 / 117.55,
    Zn-N2-H 50.0 / 107.41, S-N2-H 61.55 / 118.06, H-N2-S-C4 0.2624 (n = 1, phase 180) and 2.703 (n = 2, phase 0)."""
    p = {k: 0.0 for k, _ in FEATURES}
    name = {frozenset(('ZN', 'S')): 'ZnN2S', frozenset(('ZN', 'HN21')): 'ZnN2H', frozenset(('S', 'HN21')): 'SN2H'}
    for r in ff.removed['rest']:
        if r[0] == 'angle':
            n = name[frozenset((r[1], r[3]))]
            p[f'K_{n}'] = r[4]
            p[f'B_{n}'] = -2.0 * r[4] * np.radians(r[5])
        else:
            atoms, per, phase, k = r[1].split('-'), r[2], r[3], r[4]
            assert set(atoms) == {'C4', 'S', 'N2', 'HN21'} and abs(phase % 180.0) < 1e-3, r
            p[f'H_NS_C4_{per}'] = k * np.cos(np.radians(phase))             # k (1 + cos(n phi - phase))
    return p


def amber_terms(p):
    """Parameter dictionary -> AMBER numbers: angles (K, theta0), torsions (PK, phase, n)."""
    out = dict(angles={}, torsions={}, improper=None)
    for name, types in (('ZnN2S', 'M1-Y4-sy'), ('ZnN2H', 'M1-Y4-hn'), ('SN2H', 'hn-Y4-sy')):
        k, b = p.get(f'K_{name}', 0.0), p.get(f'B_{name}', 0.0)
        out['angles'][types] = dict(K=k, theta0=float(np.degrees(-b / (2.0 * k))) if k else None)
    for stem, types in (('H_NS_C4', 'ca-sy-Y4-hn'), ('Zn_NS_C4', 'M1-Y4-sy-ca')):
        if any(f'{stem}_{n}' in p for n in (1, 2, 3)):
            out['torsions'][types] = [dict(n=n, PK=abs(p.get(f'{stem}_{n}', 0.0)),
                                           phase=0.0 if p.get(f'{stem}_{n}', 0.0) >= 0 else 180.0) for n in (1, 2, 3)]
    for end, t in (('S', 'sy'), ('H', 'hn')):
        if f'His_ZnN_{end}_3c' in p:
            a, b = p[f'His_ZnN_{end}_3c'], p.get(f'His_ZnN_{end}_3s', 0.0)
            term = [dict(n=3, PK=float(np.hypot(a, b)), phase=float(np.degrees(np.arctan2(b, a)) % 360.0))]
            for y in ('Y1', 'Y2', 'Y3'):
                out['torsions'][f'{y}-M1-Y4-{t}'] = term
    if 'improper' in p:
        c = p['improper']                                                 # c (1 - cos 2w) = |c| (1 + cos(2w - 180)) for c > 0
        out['improper'] = {'M1-hn-Y4-sy': dict(PK=abs(c), phase=180.0 if c >= 0 else 0.0, n=2)}
    return out


def write_frcmod(p, path, title):
    """The frcmod of a model: every periodicity 1 to 3 of a refitted torsion is written, zero included, so that
    no term of an earlier file survives beside it."""
    t = amber_terms(p)
    o = [title, 'MASS', '', 'BOND', '', 'ANGLE']
    for types, v in t['angles'].items():
        if v['K']:
            o.append(f"{types:<10s}{v['K']:9.3f}  {v['theta0']:9.3f}    gate G-R2, md/N2.md")
    o += ['', 'DIHE']
    for types, terms in t['torsions'].items():
        for k, term in enumerate(terms):
            last = k == len(terms) - 1
            o.append(f"{types:<12s}  1  {term['PK']:9.4f}  {term['phase']:8.3f}  {term['n'] if last else -term['n']:5.1f}"
                     f"    gate G-R2, md/N2.md")
    o += ['', 'IMPROPER']
    for types, v in (t['improper'] or {}).items():
        o.append(f"{types:<12s}  {v['PK']:9.4f}  {v['phase']:8.1f}  {v['n']:5.1f}    gate G-R2, md/N2.md")
    o += ['', 'NONBON', '', '']
    pathlib.Path(path).write_text('\n'.join(o))
    return t


def build_topology(extra_frcmod, work):
    """tleap: the small model from the inputs kept in md/qm/n2/ff/ (those of gate G-R, cap scheme "methyl"),
    with `extra_frcmod` loaded last (None: the terms as built).  Returns the prmtop path."""
    import shutil
    work = pathlib.Path(work)
    work.mkdir(parents=True, exist_ok=True)
    for f in (FFDIR / 'inputs').iterdir():
        shutil.copyfile(f, work / f.name)
    leap = (work / 'leap.in').read_text().splitlines()
    out = []
    for l in leap:
        if l.startswith('mol = loadpdb') and extra_frcmod:
            shutil.copyfile(extra_frcmod, work / 'n2_extra.frcmod')
            out.append('loadamberparams n2_extra.frcmod')
        out.append(l)
    (work / 'leap_n2.in').write_text('\n'.join(out) + '\n')
    if (work / 'small.prmtop').exists():
        (work / 'small.prmtop').unlink()
    env = dict(os.environ, AMBERHOME=AMBERHOME, PATH=f'{AMBERHOME}/bin:' + os.environ.get('PATH', ''))
    p = subprocess.run(['tleap', '-s', '-f', 'leap_n2.in'], cwd=work, env=env, capture_output=True, text=True)
    (work / 'leap_n2.log').write_text(p.stdout + p.stderr)
    log = p.stdout + p.stderr
    bad = [l for l in log.splitlines() if 'Could not find' in l or 'FATAL' in l]
    if not (work / 'small.prmtop').exists() or bad:
        raise SystemExit(f'tleap failed: {bad[:5]}; see {work}/leap_n2.log')
    return work / 'small.prmtop'


# --------------------------------------------------------------------------------------------- 4. the split
def split_record():
    f = N2 / 'split.json'
    rec = json.loads(f.read_text())
    want = (N2 / 'split.json.sha256').read_text().split()[0]
    if sha(f) != want:
        raise SystemExit('md/qm/n2/split.json does not have the checksum written when it was drawn: stop')
    return rec


def load_points(purpose):
    """The data of gate G-R2.  purpose = 'fit': the 36 rigid points of gate G-R and the new points flagged "fit";
    the result file of a "test" point is never opened.  purpose = 'gate': everything; only for n2_gate.py, only
    with N2_FINAL_TEST=1 in the environment, and every use is logged.
    Returns a list of dictionaries: id, set ('old', 'fit', 'test'), e_b3lyp_hartree, d3_hartree, xyz."""
    if purpose not in ('fit', 'gate'):
        raise ValueError(purpose)
    split = split_record()
    if purpose == 'gate':
        if os.environ.get('N2_FINAL_TEST') != '1':
            raise SystemExit('the test set is read once, by n2_gate.py, with N2_FINAL_TEST=1 (md/N2.md, section 3.3)')
        with open(N2 / 'test_set_access.log', 'a') as fh:
            fh.write(f"{now()}  pid {os.getpid()}  {' '.join(sys.argv)}\n")
    out = []
    d3 = json.loads((GR / 'd3.json').read_text())
    old_d3 = {k: v['e_d3bj'] for k, v in d3['points'].items()}
    for pid in split['old_rigid_points']:
        r = json.loads((GR / 'qm' / f'{pid}.json').read_text())
        if not r.get('converged'):
            continue
        out.append(dict(id=pid, set='old', e_b3lyp_hartree=r['e_tot'], d3_hartree=old_d3.get(pid),
                        xyz=str(GR / 'geom' / f'{pid}.xyz')))
    for pid, flag in split['set'].items():
        if flag == 'test' and purpose != 'gate':
            continue                                                       # never opened by a fit
        f = N2 / 'qm' / f'{pid}.json'
        if not f.exists():
            continue
        r = json.loads(f.read_text())
        if r.get('set') != flag:
            raise SystemExit(f'{pid}: the result file says set {r.get("set")!r}, the split says {flag!r}')
        if r.get('test_level'):
            raise SystemExit(f'{pid}: the result file is of a test level (rehearsal), not a result of the gate')
        if not r.get('converged') or r.get('xyz_sha256') != split['xyz_sha256'][pid]:
            continue
        out.append(dict(id=pid, set=flag, e_b3lyp_hartree=r['e_tot'], d3_hartree=r['d3_hartree'],
                        xyz=str(N2 / 'geom' / f'{pid}.xyz')))
    if purpose == 'fit':
        assert not any(p['set'] == 'test' for p in out)
    return out


# ------------------------------------------------------------------------------ 5. the fit (shared with the gate)
def qm_optimum():
    """The quantum-chemical optimum of the small model (task qm_opt, all atoms free, B3LYP/6-31G*), pinned in
    md/qm/n2/qm_optimum.json by md/qm/n2_diagnosis.py before any new energy: the three angles at N2 in degrees."""
    rec = json.loads((N2 / 'qm_optimum.json').read_text())
    return rec, {'ZnN2S': rec['angles_at_N2']['Zn_N2_S'], 'ZnN2H': rec['angles_at_N2']['Zn_N2_H'],
                 'SN2H': rec['angles_at_N2']['S_N2_H']}


def model_names(groups):
    return [k for k, g in FEATURES if g in groups]


def fit_model(rows, groups, ridge=0.0, bounded=True, theta_opt=None):
    """Least squares of  y_i = E_QM(i) - rest(i)  on the model terms of `groups` and one constant.

    rows: dictionaries with 'features' (n2_common.features), 'y' (kcal/mol; any common zero).
    ridge: penalty  ridge * N * sum c_j^2  on the Fourier and improper coefficients only (not on the angles,
    not on the constant).  bounded: the bounds of the plan as constraints (BOUNDS): K in 10..300, theta0 within
    15 deg of the quantum-chemical optimum, every Fourier coefficient at most 5 (for a cosine and sine pair:
    the amplitude), the improper at most 5.  Solved in the centred form  K (th - r)^2 + b (th - r), r = the
    optimum angle, in which the bounds on theta0 are the linear conditions |b| <= 2 K * 15 deg; SLSQP from the
    start K = 50, theta0 = r, all coefficients 0.  Returns the parameter dictionary in the form of FEATURES
    (K, B = -2 K theta0, coefficients), the constant, the residuals and the AMBER numbers."""
    from scipy.optimize import minimize
    names = model_names(groups)
    if theta_opt is None:
        theta_opt = qm_optimum()[1]
    ang = [n[2:] for n in names if n.startswith('K_')]
    r0 = {a: np.radians(theta_opt[a]) for a in ang}
    key = {'ZnN2S': 'K_ZnN2S', 'ZnN2H': 'K_ZnN2H', 'SN2H': 'K_SN2H'}
    cols, col_names = [], []
    for n in names:
        if n.startswith('K_'):
            a = n[2:]
            th = np.array([r['features'][f'B_{a}'] for r in rows])        # the angle itself, radians
            cols += [(th - r0[a]) ** 2, th - r0[a]]
            col_names += [f'K_{a}', f'b_{a}']
        elif not n.startswith('B_'):
            cols.append(np.array([r['features'][n] for r in rows]))
            col_names.append(n)
    cols.append(np.ones(len(rows)))
    col_names.append('const')
    A = np.array(cols).T
    y = np.array([r['y'] for r in rows])
    pen = np.array([0.0 if (c.startswith(('K_', 'b_')) or c == 'const') else 1.0 for c in col_names])
    lam = ridge * len(rows)

    def obj(p):
        d = A @ p - y
        return float(d @ d + lam * (pen * p) @ p), 2.0 * (A.T @ d) + 2.0 * lam * pen * p
    if not bounded:
        M = A.T @ A + lam * np.diag(pen)
        p = np.linalg.lstsq(M, A.T @ y, rcond=None)[0]
        info = dict(method='normal equations, no bounds', success=True)
    else:
        # 0.001 deg inside the window: write_frcmod rounds theta0 to three decimals, and a value resting on the
        # bound must still be inside it when the gate reads it back from the topology (design check, 6 Oct 2026)
        dth = np.radians(BOUNDS['theta0_window_deg'] - 0.001)
        cmax = BOUNDS['torsion_abs_max']
        bnd, cons = [], []
        for k, c in enumerate(col_names):
            if c.startswith('K_'):
                bnd.append(BOUNDS['k_angle'])
                cons.append(dict(type='ineq', fun=lambda p, k=k: 2.0 * p[k] * dth - p[k + 1],
                                 jac=lambda p, k=k: np.eye(len(p))[k] * 2.0 * dth - np.eye(len(p))[k + 1]))
                cons.append(dict(type='ineq', fun=lambda p, k=k: 2.0 * p[k] * dth + p[k + 1],
                                 jac=lambda p, k=k: np.eye(len(p))[k] * 2.0 * dth + np.eye(len(p))[k + 1]))
            elif c.startswith('b_') or c == 'const':
                bnd.append((None, None))
            else:
                bnd.append((-cmax, cmax))
                if c.endswith('_3c'):
                    cons.append(dict(type='ineq', fun=lambda p, k=k: cmax ** 2 - p[k] ** 2 - p[k + 1] ** 2,
                                     jac=lambda p, k=k: -2.0 * p[k] * np.eye(len(p))[k] - 2.0 * p[k + 1] * np.eye(len(p))[k + 1]))
        p0 = np.zeros(len(col_names))
        for k, c in enumerate(col_names):
            if c.startswith('K_'):
                p0[k] = 50.0
        p0[-1] = float(np.mean(y - A @ p0))
        res = minimize(lambda p: obj(p)[0], p0, jac=lambda p: obj(p)[1], method='SLSQP', bounds=bnd, constraints=cons,
                       options=dict(maxiter=2000, ftol=1e-12))
        p = res.x
        info = dict(method='SLSQP', success=bool(res.success), message=str(res.message), iterations=int(res.nit))
    par = {}
    for k, c in enumerate(col_names):
        if c.startswith('K_'):
            a = c[2:]
            kk, b = float(p[k]), float(p[k + 1])
            par[f'K_{a}'] = kk
            par[f'B_{a}'] = b - 2.0 * kk * r0[a]                           # K (th - r)^2 + b (th - r) = K th^2 + (b - 2 K r) th + const
        elif not c.startswith('b_') and c != 'const':
            par[c] = float(p[k])
    resid = A @ p - y
    return dict(groups=list(groups), names=names, n_parameters=len(names), n_points=len(rows), ridge=ridge,
                bounded=bounded, parameters=par, const_centred=float(p[-1]), residuals=[float(v) for v in resid],
                rms=float(np.sqrt(np.mean(resid ** 2))), largest=float(np.abs(resid).max()), solver=info,
                amber=amber_terms(par), _A=A, _p=p, _cols=col_names, _pen=pen, _lam=lam, _r0=r0)


def leave_one_out(rows, groups, ridge=0.0, bounded=True, theta_opt=None):
    """Every point predicted by the fit to all the others (the same fit, bounds and ridge)."""
    err = []
    for i in range(len(rows)):
        sub = rows[:i] + rows[i + 1:]
        f = fit_model(sub, groups, ridge, bounded, theta_opt)
        one = fit_model.__globals__['_predict_one'](f, rows[i], theta_opt)
        err.append(one - rows[i]['y'])
    err = np.array(err)
    return dict(rms=float(np.sqrt(np.mean(err ** 2))), largest=float(np.abs(err).max()), errors=[float(v) for v in err])


def _predict_one(f, row, theta_opt=None):
    if theta_opt is None:
        theta_opt = qm_optimum()[1]
    v = []
    for n in f['names']:
        if n.startswith('K_'):
            a = n[2:]
            th = row['features'][f'B_{a}'] - np.radians(theta_opt[a])
            v += [th ** 2, th]
        elif not n.startswith('B_'):
            v.append(row['features'][n])
    v.append(1.0)
    return float(np.dot(v, f['_p']))


def active_bounds(f, tol=1e-6):
    """Which bounds of the plan the solution of fit_model rests on (names of the centred columns)."""
    cols, p = f['_cols'], f['_p']
    dth = np.radians(BOUNDS['theta0_window_deg'] - 0.001)
    cmax = BOUNDS['torsion_abs_max']
    act = []
    for k, c in enumerate(cols):
        if c.startswith('K_'):
            if p[k] <= BOUNDS['k_angle'][0] + tol:
                act.append((c, 'K at its lower bound'))
            if p[k] >= BOUNDS['k_angle'][1] - tol:
                act.append((c, 'K at its upper bound'))
            if abs(p[k + 1]) >= 2.0 * p[k] * dth - tol * max(1.0, p[k]):
                act.append((cols[k + 1], 'theta0 at the upper edge of its window' if p[k + 1] < 0 else
                            'theta0 at the lower edge of its window'))
        elif not c.startswith('b_') and c != 'const':
            if abs(p[k]) >= cmax - tol:
                act.append((c, 'coefficient at its bound'))
            if c.endswith('_3c') and p[k] ** 2 + p[k + 1] ** 2 >= cmax ** 2 - tol:
                act.append((c, 'amplitude at its bound'))
    return act


def standard_errors(f):
    """Standard errors of a fit of fit_model, by the usual linear theory in the centred parameters
    (K, b, coefficients, constant): cov = s^2 M^-1 A'A M^-1 with M = A'A + ridge, s^2 = RSS / (N - trace of the
    hat matrix).  A parameter that rests on a bound has no standard error (it is not free): its column is taken
    out, and for an equilibrium angle on the edge of its window b is tied to K.  theta0 = r - b / (2 K) and the
    amplitude and phase of a cosine and sine pair by the delta method.  Returns AMBER-form numbers with errors:
    angles K (kcal/mol/rad^2) and theta0 (deg); coefficients (kcal/mol)."""
    A, p, cols, pen, lam, r0 = f['_A'], f['_p'], f['_cols'], f['_pen'], f['_lam'], f['_r0']
    n, m = A.shape
    act = active_bounds(f)
    dth = np.radians(BOUNDS['theta0_window_deg'] - 0.001)
    T = np.eye(m)
    drop = set()
    for c, why in act:
        k = cols.index(c)
        if c.startswith('b_'):                                            # b = -/+ 2 K dth: tied to K
            T[k, k - 1] = (-1.0 if p[k] < 0 else 1.0) * 2.0 * dth
            drop.add(k)
        elif 'amplitude' not in why:
            drop.add(k)
    keep = [k for k in range(m) if k not in drop]
    T = T[:, keep]
    At = A @ T
    M = At.T @ At + lam * (T.T @ np.diag(pen) @ T)
    Mi = np.linalg.pinv(M)
    hat = float(np.trace(At @ Mi @ At.T))
    resid = A @ p - np.array(f['_y']) if '_y' in f else np.array(f['residuals'])
    s2 = float(resid @ resid) / max(1.0, n - hat)
    cov = T @ (s2 * Mi @ (At.T @ At) @ Mi) @ T.T
    se = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    out = dict(sigma=float(np.sqrt(s2)), effective_parameters=round(hat, 3), at_a_bound=[list(a) for a in act],
               angles={}, coefficients={}, pairs={}, constant=dict(value=float(p[-1]), se=float(se[-1])))
    fixed = {c for c, _ in act}
    for k, c in enumerate(cols):
        if c.startswith('K_'):
            a = c[2:]
            K, b = p[k], p[k + 1]
            g = np.zeros(m)
            g[k], g[k + 1] = b / (2.0 * K * K), -1.0 / (2.0 * K)
            th_se = float(np.degrees(np.sqrt(max(0.0, g @ cov @ g))))
            out['angles'][a] = dict(K=float(K), K_se=None if c in fixed else float(se[k]),
                                    theta0=float(np.degrees(r0[a] - b / (2.0 * K))),
                                    theta0_se=None if f'b_{a}' in fixed else th_se)
        elif c.endswith('_3c'):
            a_, b_ = p[k], p[k + 1]
            amp = float(np.hypot(a_, b_))
            if amp > 1e-12:
                ga = np.zeros(m); ga[k], ga[k + 1] = a_ / amp, b_ / amp
                gp = np.zeros(m); gp[k], gp[k + 1] = -b_ / amp ** 2, a_ / amp ** 2
                out['pairs'][c[:-3]] = dict(PK=amp, PK_se=float(np.sqrt(max(0.0, ga @ cov @ ga))),
                                            phase=float(np.degrees(np.arctan2(b_, a_)) % 360.0),
                                            phase_se=float(np.degrees(np.sqrt(max(0.0, gp @ cov @ gp)))),
                                            cos=float(a_), cos_se=float(se[k]), sin=float(b_), sin_se=float(se[k + 1]))
            else:
                out['pairs'][c[:-3]] = dict(PK=0.0, PK_se=None, phase=0.0, phase_se=None, cos=float(a_), cos_se=float(se[k]),
                                            sin=float(b_), sin_se=float(se[k + 1]))
        elif not c.startswith('b_') and not c.endswith('_3s') and c != 'const':
            out['coefficients'][c] = dict(value=float(p[k]), se=None if c in fixed else float(se[k]))
    return out


def one_reference(errors, e_qm):
    """The leave-one-out errors read against ONE reference, as condition 1 of the gate reads the test set: r =
    the row of lowest E_QM; root mean square and largest of e_i - e_r over the other rows (design check of
    6 Oct 2026, finding 1).  Fitting data only."""
    e = np.asarray(errors, float)
    r = int(np.argmin(np.asarray(e_qm, float)))
    d = np.delete(e - e[r], r)
    return dict(reference_row=r, rms=float(np.sqrt(np.mean(d ** 2))), largest=float(d[np.argmax(np.abs(d))]),
                mean=float(d.mean()), error_of_the_reference=float(e[r]))


def old_family(pid):
    """The family of a rigid point of gate G-R by its name (md/GR.md, section 3): a path, d hydrogen on X, g
    hydrogen on Rx, e turns about Zn-N2, f bend of Zn-N2-S, c1 and c2 turns about N2-S, and X, Rx, Nx."""
    if pid in ('X', 'Rx', 'Nx'):
        return pid
    for pre in ('c1', 'c2', 'a', 'd', 'e', 'f', 'g'):
        if pid.startswith(pre):
            return pre
    return 'other'


def mm_energies():
    """{id: dict(total=, rest=)}: the force-field energy of every point (old and new) without a fitted term changed,
    and `rest`, the same without the terms that may change.  Copy 'ff': from md/qm/n2/points.json (written by
    n2_points.py).  Another copy (N2_FF): from <copy>/mm_points.json, written by md/qm/n2_charges.py ff."""
    if FF_NAME == 'ff':
        pj = json.loads((N2 / 'points.json').read_text())
        out = {k: v['mm_as_built'] for k, v in pj['old_rigid_points'].items()}
        out.update({p['id']: p['mm_as_built'] for p in pj['points']})
        return out
    return json.loads((FFDIR / 'mm_points.json').read_text())['points']


def fit_rows(purpose='fit', mm=None):
    """The rows of a fit: one per point of load_points(purpose) with its features, y = E_QM - E(rest) and labels.
    The force-field energies are those of the copy in use (mm_energies())."""
    pj = json.loads((N2 / 'points.json').read_text())
    feat = dict(pj['old_rigid_points'])
    feat.update({p['id']: p for p in pj['points']})
    if mm is None:
        mm = mm_energies()
    rows = []
    for p in load_points(purpose):
        f = feat[p['id']]
        e_qm = (p['e_b3lyp_hartree'] + p['d3_hartree']) * HARTREE
        e = mm[p['id']]
        sv = f['coords']['signed_volume']
        rows.append(dict(id=p['id'], set=p['set'], features=f['features'], y=e_qm - e['rest'], e_qm=e_qm,
                         e_qm_no_d3=p['e_b3lyp_hartree'] * HARTREE, mm_total=e['total'], mm_rest=e['rest'],
                         region=f.get('region') or ('R' if sv < -0.3 else 'X' if sv > 0.3 else 'between'),
                         family=old_family(p['id']) if p['set'] == 'old' else None,
                         kind=f.get('kind', 'old'), layer=f.get('layer'), protein_contact=f.get('protein_contact'),
                         coords=f['coords']))
    return rows


def within_bounds(amber, theta_opt):
    """The bounds of the plan on AMBER numbers (as read from a parameter set or from a topology)."""
    bad = []
    opt = {'M1-Y4-sy': theta_opt['ZnN2S'], 'M1-Y4-hn': theta_opt['ZnN2H'], 'hn-Y4-sy': theta_opt['SN2H']}
    for t, v in amber['angles'].items():
        if not v['K']:
            continue
        if not (BOUNDS['k_angle'][0] - 1e-6 <= v['K'] <= BOUNDS['k_angle'][1] + 1e-6):
            bad.append(f"{t}: K {v['K']:.3f} outside {BOUNDS['k_angle']}")
        if abs(v['theta0'] - opt[t]) > BOUNDS['theta0_window_deg'] + 1e-6:
            bad.append(f"{t}: theta0 {v['theta0']:.3f} more than {BOUNDS['theta0_window_deg']} deg from {opt[t]:.3f}")
    for t, terms in amber['torsions'].items():
        for term in terms:
            if term['n'] not in (1, 2, 3) or abs(term['PK']) > BOUNDS['torsion_abs_max'] + 1e-6:
                bad.append(f"{t}: n {term['n']}, PK {term['PK']:.4f}")
    for t, v in (amber['improper'] or {}).items():
        if abs(v['PK']) > BOUNDS['improper_abs_max'] + 1e-6:
            bad.append(f"{t}: improper PK {v['PK']:.4f}")
    return bad
