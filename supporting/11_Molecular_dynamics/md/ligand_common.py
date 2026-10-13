"""Helpers of stage S1 (ligand parameters) of md/PLAN.md.  Imported by ligand_params.py and ligand_fit.py.

Holds: paths, the atom naming (heavy atoms as in SAN of 6RL9), internal-coordinate functions, readers for the DFT
reference data of the paper (molecules/sfa/results/), the crystal ligand, gaff2.dat and frcmod files, and a small
molecular-mechanics engine (OpenMM, Reference platform, double precision) with a tight minimiser and an exact-enough
dihedral constraint for relaxed scans.

Run with the MD set-up environment:  /root/cqc-sfa-tools/md/bin/python
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '2')          # thread limit of this task (shared machine)
import gzip, hashlib, json, math, pathlib, re, shutil, subprocess, sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
MD = ROOT / 'md'
OUT = MD / 'ligand'
BULK = pathlib.Path('/root/cqc-sfa-md/ligand')
AMBERHOME = pathlib.Path('/root/cqc-sfa-tools/md')
GAFF2_DAT = AMBERHOME / 'dat' / 'leap' / 'parm' / 'gaff2.dat'
GAFF1_DAT = AMBERHOME / 'dat' / 'leap' / 'parm' / 'gaff.dat'
RESULTS = ROOT / 'molecules' / 'sfa' / 'results'
CIF_6RL9 = ROOT / 'dock' / 'pdb' / '6RL9.cif.gz'
KCAL = 4.184                      # kJ per kcal
HARTREE_KCAL = 627.5094740631

# Atom order of every DFT file of the paper (dock/results/poses_for_dft/manifest.json: "C1..C6 N1 N2 O1 O2 S1 H...").
# Heavy-atom names follow SAN of 6RL9 (the sulfur is "S" there, "S1" in the DFT files and in the paper).
HEAVY = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1', 'N2', 'O1', 'O2', 'S']
SPECIES = {
    'anion':   dict(case='sulfonamidate_pcm_d3', charge=-1, label='sulfonamidate anion',
                    hnames=['H2', 'H3', 'H5', 'H6', 'HN11', 'HN12', 'HN21']),
    'neutral': dict(case='sulfanilamide_pcm_d3', charge=0, label='sulfanilamide (neutral)',
                    hnames=['H2', 'H3', 'H5', 'H6', 'HN11', 'HN12', 'HN21', 'HN22']),
}
BONDS_HEAVY = [('C1', 'C2'), ('C2', 'C3'), ('C3', 'C4'), ('C4', 'C5'), ('C5', 'C6'), ('C6', 'C1'), ('C1', 'N1'),
               ('C4', 'S'), ('S', 'N2'), ('S', 'O1'), ('S', 'O2')]
H_PARENT = {'H2': 'C2', 'H3': 'C3', 'H5': 'C5', 'H6': 'C6', 'HN11': 'N1', 'HN12': 'N1', 'HN21': 'N2', 'HN22': 'N2'}


def env():
    e = dict(os.environ)
    e['AMBERHOME'] = str(AMBERHOME)
    e['PATH'] = f"{AMBERHOME / 'bin'}:{e.get('PATH', '')}"
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        e[v] = '2'
    return e


def run(cmd, cwd, log=None):
    """Run an Amber tool (niced); returns stdout+stderr, raises on a non-zero exit."""
    p = subprocess.run(['nice', '-n', '10'] + [str(c) for c in cmd], cwd=str(cwd), env=env(),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if log:
        pathlib.Path(log).write_text(p.stdout)
    if p.returncode != 0:
        raise RuntimeError(f"{' '.join(map(str, cmd))} failed in {cwd}:\n{p.stdout[-2000:]}")
    return p.stdout


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def versions():
    import openmm, parmed, scipy, matplotlib
    meta = sorted((AMBERHOME / 'conda-meta').glob('ambertools-*.json'))
    return dict(python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__,
                matplotlib=matplotlib.__version__, openmm=openmm.__version__, parmed=parmed.__version__,
                ambertools=meta[0].stem if meta else None, amberhome=str(AMBERHOME),
                gaff2_dat=dict(path=str(GAFF2_DAT), title=GAFF2_DAT.open().readline().strip(),
                               sha256=sha256(GAFF2_DAT)),
                openmm_platform='Reference (double precision, one thread)')


# ------------------------------------------------------------------------------------------- geometry --
def dist(x, i, j):
    return float(np.linalg.norm(x[i] - x[j]))


def angle(x, i, j, k):
    a, b = x[i] - x[j], x[k] - x[j]
    return math.degrees(math.acos(max(-1.0, min(1.0, float(a @ b) / np.linalg.norm(a) / np.linalg.norm(b)))))


def dihedral(x, i, j, k, l):
    """Dihedral i-j-k-l in degrees, (-180, 180], IUPAC sign."""
    b0, b1, b2 = x[i] - x[j], x[k] - x[j], x[l] - x[k]
    b1n = b1 / np.linalg.norm(b1)
    v, w = b0 - (b0 @ b1n) * b1n, b2 - (b2 @ b1n) * b1n
    return math.degrees(math.atan2(float(np.cross(b1n, v) @ w), float(v @ w)))


def wrap(d):
    return (d + 180.0) % 360.0 - 180.0


def measures(x, names):
    """Internal coordinates at sulfur and the two torsions, for a coordinate array in Angstrom and its atom names."""
    ix = {n: i for i, n in enumerate(names)}
    C3, C4, C5, S, N2, O1, O2, N1, C1 = (ix[n] for n in ('C3', 'C4', 'C5', 'S', 'N2', 'O1', 'O2', 'N1', 'C1'))
    t3, t5 = dihedral(x, C3, C4, S, N2), dihedral(x, C5, C4, S, N2)
    m = dict(a_C4_S_N2=angle(x, C4, S, N2), a_O1_S_O2=angle(x, O1, S, O2),
             a_O1_S_N2=angle(x, O1, S, N2), a_O2_S_N2=angle(x, O2, S, N2),
             a_C4_S_O1=angle(x, C4, S, O1), a_C4_S_O2=angle(x, C4, S, O2),
             d_S_N2=dist(x, S, N2), d_C4_S=dist(x, C4, S), d_S_O1=dist(x, S, O1), d_S_O2=dist(x, S, O2),
             d_C1_N1=dist(x, C1, N1),
             tau_C3=t3, tau_C5=t5, tau_fold=max(abs(t3), abs(t5)))
    # sulfur out of the ring plane: angle between C4->S and the plane C3,C4,C5 (sign: side of N2... not used)
    n = np.cross(x[C3] - x[C4], x[C5] - x[C4]); n /= np.linalg.norm(n)
    v = x[S] - x[C4]
    m['S_out_of_plane'] = math.degrees(math.asin(float(abs(v @ n)) / np.linalg.norm(v)))
    hs = [ix[h] for h in ('HN21', 'HN22') if h in ix]
    for h in hs:
        m[f'tau_C4_S_N2_{names[h]}'] = dihedral(x, C4, S, N2, h)
        m[f'a_S_N2_{names[h]}'] = angle(x, S, N2, h)
    if len(hs) == 2:
        m['N2_angle_sum'] = angle(x, S, N2, hs[0]) + angle(x, S, N2, hs[1]) + angle(x, hs[0], N2, hs[1])
    h1 = [ix[h] for h in ('HN11', 'HN12') if h in ix]
    if len(h1) == 2:
        m['N1_angle_sum'] = angle(x, C1, N1, h1[0]) + angle(x, C1, N1, h1[1]) + angle(x, h1[0], N1, h1[1])
    return m


GEOM_ROWS = [  # (key, label, kind) in the order of the report tables
    ('a_C4_S_N2', 'C4-S-N2', 'deg'), ('a_O1_S_O2', 'O1-S-O2', 'deg'), ('a_O_S_N2_mean', 'O-S-N2 (mean of two)', 'deg'),
    ('a_C4_S_O_mean', 'C4-S-O (mean of two)', 'deg'), ('d_S_N2', 'S-N2', 'A'), ('d_C4_S', 'C4-S', 'A'),
    ('d_S_O_mean', 'S-O (mean of two)', 'A'), ('tau_fold', 'tau_fold', 'deg'), ('tau_C3', 'C3-C4-S-N2', 'deg'),
    ('tau_C4_S_N2_HN21', 'C4-S-N2-HN21', 'deg'), ('tau_C4_S_N2_HN22', 'C4-S-N2-HN22', 'deg'),
    ('N2_angle_sum', 'angle sum at N2', 'deg'), ('S_out_of_plane', 'S out of ring plane', 'deg')]


def add_means(m):
    m = dict(m)
    m['a_O_S_N2_mean'] = 0.5 * (m['a_O1_S_N2'] + m['a_O2_S_N2'])
    m['a_C4_S_O_mean'] = 0.5 * (m['a_C4_S_O1'] + m['a_C4_S_O2'])
    m['d_S_O_mean'] = 0.5 * (m['d_S_O1'] + m['d_S_O2'])
    return m


# ------------------------------------------------------------------------------------------ DFT data --
def read_xyz(path):
    lines = pathlib.Path(path).read_text().split('\n')
    n = int(lines[0])
    el, xyz = [], []
    for l in lines[2:2 + n]:
        s = l.split()
        el.append(s[0]); xyz.append([float(v) for v in s[1:4]])
    return el, np.array(xyz), lines[1]


def names_for(species, el, x):
    """Atom names for a DFT geometry of the species; checks the documented atom order against the geometry."""
    sp = SPECIES[species]
    names = HEAVY + sp['hnames']
    assert el[:11] == ['C'] * 6 + ['N', 'N', 'O', 'O', 'S'] and all(e == 'H' for e in el[11:]), el
    assert len(el) == len(names), (len(el), names)
    ix = {n: i for i, n in enumerate(names)}
    for a, b in BONDS_HEAVY:
        assert dist(x, ix[a], ix[b]) < 1.95, (a, b, dist(x, ix[a], ix[b]))
    for h in sp['hnames']:                       # each hydrogen sits on the heavy atom its name says
        d = [dist(x, ix[h], j) for j in range(11)]
        assert int(np.argmin(d)) == ix[H_PARENT[h]] and min(d) < 1.15, (h, HEAVY[int(np.argmin(d))], min(d))
    return names


def dft_scan(species):
    """Relaxed DFT scan of the paper: list of points (dict), the unconstrained reference first.

    Each point: C3-C4-S1-N2 held at +tau (geomeTRIC constraint), every other coordinate relaxed, started from the
    unconstrained minimum with the sulfur side rotated rigidly (dft/torsion_scan.py).  B3LYP-D3(BJ)/6-311++G(d,p),
    IEF-PCM water.  Energies in kcal/mol relative to the unconstrained re-optimisation (ref.json)."""
    d = RESULTS / SPECIES[species]['case'] / 'torsion_scan'
    ref = json.loads((d / 'ref.json').read_text())
    pts = []
    for f in sorted(d.glob('*.json')):
        if f.name == 'torsion_scan_summary.json':
            continue
        j = json.loads(f.read_text())
        if 'metrics' not in j:
            continue
        el, x, _ = read_xyz(f.with_suffix('.xyz'))
        names = names_for(species, el, x)
        m = measures(x, names)
        # self-check of the atom mapping and of the internal-coordinate functions against the stored metrics
        assert abs(wrap(m['tau_C3'] - j['metrics']['tau_C3_deg'])) < 2e-3, (f.name, m['tau_C3'], j['metrics'])
        assert abs(m['a_C4_S_N2'] - j['metrics']['a_C4_S1_N2_deg']) < 2e-3
        assert abs(m['d_S_N2'] - j['metrics']['d_S1_N2_A']) < 2e-5
        pts.append(dict(file=str(f.relative_to(ROOT)), kind=j['kind'], branch=j.get('branch'),
                        target=j.get('target_tau_fold_deg'), converged=bool(j['converged']),
                        E_Ha=j['E_Ha'], E_rel=(j['E_Ha'] - ref['E_Ha']) * HARTREE_KCAL,
                        level=j['level'], x=x, names=names, m=m))
    pts.sort(key=lambda p: (p['kind'] != 'reference', p['branch'] != 'C3', p['target'] or 0.0))
    assert pts[0]['kind'] == 'reference'
    return pts


def crystal_ligand(seq=304):
    """Heavy atoms of SAN <seq> of the deposited 6RL9 (mmCIF, atom_site loop parsed here)."""
    lines = gzip.open(CIF_6RL9, 'rt').read().split('\n')
    i = next(k for k, l in enumerate(lines) if l.startswith('_atom_site.group_PDB'))
    cols = []
    while lines[i].startswith('_atom_site.'):
        cols.append(lines[i].split('.')[1].strip()); i += 1
    c = {n: k for k, n in enumerate(cols)}
    atoms = {}
    while i < len(lines) and lines[i] and not lines[i].startswith('#'):
        s = lines[i].split(); i += 1
        if s[c['label_comp_id']] == 'SAN' and int(s[c['auth_seq_id']]) == seq:
            assert s[c['label_alt_id']] == '.', 'alternate location'
            atoms[s[c['label_atom_id']]] = dict(xyz=[float(s[c[k]]) for k in ('Cartn_x', 'Cartn_y', 'Cartn_z')],
                                                occ=float(s[c['occupancy']]), b=float(s[c['B_iso_or_equiv']]))
    assert sorted(atoms) == sorted(HEAVY), sorted(atoms)
    x = np.array([atoms[n]['xyz'] for n in HEAVY])
    return HEAVY, x, atoms


# --------------------------------------------------------------------------- force-field parameter files --
def _key(types):
    return tuple(t.strip() for t in types)


def read_parm_dat(path):
    """Bond, angle, dihedral and improper sections of an Amber parm .dat file (gaff2.dat), with line numbers.

    Returns dict(bond={(a,b): rec}, angle={(a,b,c): rec}, dihedral={(a,b,c,d): [rec, ...]}, improper={...}).
    Keys are stored as written; look-ups try both directions.  Dihedral force constants are PK/IDIVF."""
    out = dict(bond={}, angle={}, dihedral={}, improper={})
    lines = pathlib.Path(path).read_text().split('\n')
    sec = 'atoms'
    k = 1
    while k < len(lines):
        l = lines[k]
        if sec == 'atoms':
            if not l.strip():
                sec = 'hydro'
        elif sec == 'hydro':
            sec = 'bond'            # one line of hydrophilic atom types
        elif sec == 'bond':
            if not l.strip():
                sec = 'angle'
            else:
                a, b = l[0:2], l[3:5]; v = l[5:].split()
                out['bond'][_key((a, b))] = dict(k=float(v[0]), eq=float(v[1]), line=k + 1, text=l.rstrip())
        elif sec == 'angle':
            if not l.strip():
                sec = 'dihedral'
            else:
                a, b, c = l[0:2], l[3:5], l[6:8]; v = l[8:].split()
                out['angle'][_key((a, b, c))] = dict(k=float(v[0]), eq=float(v[1]), line=k + 1, text=l.rstrip())
        elif sec == 'dihedral':
            if not l.strip():
                sec = 'improper'
            else:
                t = _key((l[0:2], l[3:5], l[6:8], l[9:11])); v = l[11:].split()
                rec = dict(idivf=int(float(v[0])), pk=float(v[1]), phase=float(v[2]), pn=abs(int(float(v[3]))),
                           line=k + 1, text=l.rstrip())
                rec['k'] = rec['pk'] / rec['idivf']
                out['dihedral'].setdefault(t, []).append(rec)
        elif sec == 'improper':
            if not l.strip():
                break
            t = _key((l[0:2], l[3:5], l[6:8], l[9:11])); v = l[11:].split()
            out['improper'][t] = dict(pk=float(v[0]), phase=float(v[1]), pn=abs(int(float(v[2]))), line=k + 1,
                                      text=l.rstrip())
        k += 1
    return out


def read_frcmod(path):
    """Sections of a frcmod file: same record layout as read_parm_dat, plus the free-text comment of each line."""
    out = dict(bond={}, angle={}, dihedral={}, improper={}, mass={}, nonbon={})
    sec = None
    for n, l in enumerate(pathlib.Path(path).read_text().split('\n')):
        s = l.strip()
        if s[:4] in ('MASS', 'BOND', 'ANGL', 'DIHE', 'IMPR', 'NONB'):
            sec = s[:4]; continue
        if not s or sec is None:
            continue
        if sec == 'BOND':
            v = l[5:].split()
            out['bond'][_key((l[0:2], l[3:5]))] = dict(k=float(v[0]), eq=float(v[1]), comment=' '.join(v[2:]),
                                                      line=n + 1)
        elif sec == 'ANGL':
            v = l[8:].split()
            out['angle'][_key((l[0:2], l[3:5], l[6:8]))] = dict(k=float(v[0]), eq=float(v[1]),
                                                                comment=' '.join(v[2:]), line=n + 1)
        elif sec == 'DIHE':
            t = _key((l[0:2], l[3:5], l[6:8], l[9:11])); v = l[11:].split()
            rec = dict(idivf=int(float(v[0])), pk=float(v[1]), phase=float(v[2]), pn=abs(int(float(v[3]))),
                       comment=' '.join(v[4:]), line=n + 1)
            rec['k'] = rec['pk'] / rec['idivf']
            out['dihedral'].setdefault(t, []).append(rec)
        elif sec == 'IMPR':
            t = _key((l[0:2], l[3:5], l[6:8], l[9:11])); v = l[11:].split()
            out['improper'][t] = dict(pk=float(v[0]), phase=float(v[1]), pn=abs(int(float(v[2]))),
                                      comment=' '.join(v[3:]), line=n + 1)
    return out


def lookup(table, types, wildcard=False):
    """Exact-type look-up in either direction; for dihedrals optionally the X -b-c-X form.  Returns (key, rec)."""
    t = tuple(types)
    for key in (t, t[::-1]):
        if key in table:
            return key, table[key]
    if wildcard and len(t) == 4:
        for key in (('X', t[1], t[2], 'X'), ('X', t[2], t[1], 'X')):
            if key in table:
                return key, table[key]
    return None, None


# ------------------------------------------------------------------------------------ molecular mechanics --
class MM:
    """Vacuum (or implicit-solvent) molecular mechanics of one small molecule from an Amber topology.

    OpenMM Reference platform (double precision).  Units at the interface: Angstrom and kcal/mol.
    A harmonic restraint on one dihedral with a very large force constant (K_RESTR) serves as the constraint of the
    relaxed scans; the reported energy never contains it, and the deviation from the target is returned."""
    K_RESTR = 2.0e6           # kJ mol-1 rad-2: a torque of 10 kJ/mol/rad moves the dihedral by 3e-4 degree

    def __init__(self, prmtop, scan_atoms=('C3', 'C4', 'S', 'N2'), implicit=None, dielectric=78.3553):
        import parmed, openmm
        from openmm import app, unit
        self.parm = parmed.load_file(str(prmtop))
        self.names = [a.name for a in self.parm.atoms]
        self.types = [a.type for a in self.parm.atoms]
        self.ix = {n: i for i, n in enumerate(self.names)}
        kw = dict(nonbondedMethod=app.NoCutoff, constraints=None, removeCMMotion=False)
        if implicit:
            kw.update(implicitSolvent=getattr(app, implicit), solventDielectric=dielectric, soluteDielectric=1.0)
        self.system = self.parm.createSystem(**kw)
        for f in self.system.getForces():
            f.setForceGroup(0)
        self.scan_idx = [self.ix[n] for n in scan_atoms]
        r = openmm.CustomTorsionForce('0.5*kr*min(dt, 6.283185307179586-dt)^2; dt=abs(theta-theta0)')
        r.addGlobalParameter('kr', 0.0)
        r.addPerTorsionParameter('theta0')
        r.addTorsion(*self.scan_idx, [0.0])
        r.setForceGroup(1)
        self.restraint = r
        self.system.addForce(r)
        self.torsion_force = next(f for f in self.system.getForces() if isinstance(f, openmm.PeriodicTorsionForce))
        self.bond_force = next(f for f in self.system.getForces() if isinstance(f, openmm.HarmonicBondForce))
        self.angle_force = next(f for f in self.system.getForces() if isinstance(f, openmm.HarmonicAngleForce))
        self._openmm, self._unit = openmm, unit
        self.context = None
        self._build()

    def _build(self):
        openmm = self._openmm
        integ = openmm.VerletIntegrator(0.001)
        self.context = openmm.Context(self.system, integ, openmm.Platform.getPlatformByName('Reference'))

    def rebuild(self):
        """After adding terms to a force (not only changing parameters)."""
        self._build()

    def _set(self, x):
        self.context.setPositions(np.asarray(x) * 0.1)

    def energy(self, x, with_forces=False, restrained=False):
        """Potential energy in kcal/mol (force group 0 only unless restrained=True) and, optionally, gradient."""
        unit = self._unit
        self._set(x)
        groups = {0, 1} if restrained else {0}
        st = self.context.getState(getEnergy=True, getForces=with_forces, groups=groups)
        e = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) / KCAL
        if not with_forces:
            return e
        f = st.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
        return e, -np.asarray(f) * 0.1 / KCAL            # gradient, kcal/mol/A

    def minimise(self, x0, tau=None, gtol=1e-6, maxiter=20000):
        """Minimise from x0 (Angstrom).  tau (deg): hold the scan dihedral there.  Returns dict(x, E, ...).

        L-BFGS-B on energy and analytic gradient, restarted until the gradient criterion (kcal/mol/A, largest
        component of the gradient of the minimised function) is met or nothing changes."""
        from scipy.optimize import minimize
        if tau is None:
            self.context.setParameter('kr', 0.0)
        else:
            self.context.setParameter('kr', self.K_RESTR)
            self.restraint.setTorsionParameters(0, *self.scan_idx, [math.radians(tau)])
            self.restraint.updateParametersInContext(self.context)
        n = len(self.names)

        def fun(v):
            e, g = self.energy(v.reshape(n, 3), with_forces=True, restrained=tau is not None)
            return e, g.ravel()
        v = np.asarray(x0, float).ravel().copy()
        nit = 0
        for _ in range(40):
            r = minimize(fun, v, jac=True, method='L-BFGS-B',
                         options=dict(maxiter=maxiter, maxfun=10 * maxiter, ftol=1e-15, gtol=gtol, maxcor=30))
            nit += r.nit
            moved = float(np.abs(r.x - v).max())
            v = r.x
            if float(np.abs(r.jac).max()) < gtol or moved < 1e-9:
                break
        x = v.reshape(n, 3)
        e, g = self.energy(x, with_forces=True, restrained=tau is not None)
        res = dict(x=x, E=self.energy(x), gmax=float(np.abs(g).max()), nit=nit)
        if tau is not None:
            res['tau_dev'] = wrap(dihedral(x, *self.scan_idx) - tau)
            res['E_restraint'] = e - res['E']
        self.context.setParameter('kr', 0.0)
        return res

    def hessian_min_eig(self, x, h=1e-3):
        """Lowest eigenvalues of the finite-difference Hessian (kcal/mol/A^2): six near zero and no negative one
        beyond them means a true minimum."""
        n = len(self.names) * 3
        H = np.zeros((n, n))
        x = np.asarray(x, float)
        for i in range(n):
            xp = x.ravel().copy(); xm = xp.copy(); xp[i] += h; xm[i] -= h
            gp = self.energy(xp.reshape(-1, 3), with_forces=True)[1].ravel()
            gm = self.energy(xm.reshape(-1, 3), with_forces=True)[1].ravel()
            H[i] = (gp - gm) / (2 * h)
        w = np.linalg.eigvalsh(0.5 * (H + H.T))
        return w

    # -- parameter access for the repair (types -> terms of the OpenMM system) -------------------------
    def bonds_of(self, t1, t2):
        out = []
        for k in range(self.bond_force.getNumBonds()):
            i, j, r0, kk = self.bond_force.getBondParameters(k)
            if {self.types[i], self.types[j]} == {t1, t2} and (t1 != t2 or self.types[i] == self.types[j]):
                out.append(k)
        return out

    def angles_of(self, t1, t2, t3):
        out = []
        for k in range(self.angle_force.getNumAngles()):
            i, j, l, th, kk = self.angle_force.getAngleParameters(k)
            ty = (self.types[i], self.types[j], self.types[l])
            if ty == (t1, t2, t3) or ty == (t3, t2, t1):
                out.append(k)
        return out

    def set_bond_eq(self, t1, t2, r0_A):
        unit = self._unit
        for k in self.bonds_of(t1, t2):
            i, j, r0, kk = self.bond_force.getBondParameters(k)
            self.bond_force.setBondParameters(k, i, j, r0_A * 0.1, kk)
        self.bond_force.updateParametersInContext(self.context)

    def set_angle_eq(self, t1, t2, t3, th_deg):
        for k in self.angles_of(t1, t2, t3):
            i, j, l, th, kk = self.angle_force.getAngleParameters(k)
            self.angle_force.setAngleParameters(k, i, j, l, math.radians(th_deg), kk)
        self.angle_force.updateParametersInContext(self.context)

    def proper_torsions(self):
        """Proper torsion terms of the topology as (index in the OpenMM force, i, j, k, l, n, phase_deg, k_kcal)."""
        unit = self._unit
        impropers = set()
        for d in self.parm.dihedrals:
            if d.improper:
                impropers.add((d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx))
        out = []
        for t in range(self.torsion_force.getNumTorsions()):
            i, j, k, l, n, ph, kk = self.torsion_force.getTorsionParameters(t)
            if (i, j, k, l) in impropers:
                continue
            out.append((t, i, j, k, l, n, ph.value_in_unit(unit.degree),
                        kk.value_in_unit(unit.kilojoule_per_mole) / KCAL))
        return out


def rotate_about(x, a, b, moving, deg):
    """Rotate the atoms `moving` about the axis a->b by deg (right-handed)."""
    x = np.array(x, float)
    ax = x[b] - x[a]; ax /= np.linalg.norm(ax)
    th = math.radians(deg)
    c, s = math.cos(th), math.sin(th)
    for m in moving:
        v = x[m] - x[b]
        x[m] = x[b] + v * c + np.cross(ax, v) * s + ax * (ax @ v) * (1 - c)
    return x


def set_scan_dihedral(x, names, tau):
    """Rigid rotation of the sulfur side (S, O1, O2, N2 and its hydrogens) about C4->S so that C3-C4-S-N2 = tau,
    as the DFT scan built its starts."""
    ix = {n: i for i, n in enumerate(names)}
    cur = dihedral(x, ix['C3'], ix['C4'], ix['S'], ix['N2'])
    moving = [ix[n] for n in ('O1', 'O2', 'N2', 'HN21', 'HN22') if n in ix]
    y = rotate_about(x, ix['C4'], ix['S'], moving, tau - cur)
    assert abs(wrap(dihedral(y, ix['C3'], ix['C4'], ix['S'], ix['N2']) - tau)) < 1e-6
    return y


def write_csv(path, rows, cols=None):
    cols = cols or list(rows[0].keys())
    with open(path, 'w') as f:
        f.write(','.join(cols) + '\n')
        for r in rows:
            f.write(','.join('' if r.get(c) is None else (f'{r[c]:.6g}' if isinstance(r[c], float) else
                                                           '"' + str(r[c]).replace('"', "'") + '"'
                                                           if (',' in str(r[c])) else str(r[c])) for c in cols) + '\n')


def now_ist():
    import datetime
    t = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=5, minutes=30)))
    return t.strftime('%-d %b %Y, %-I:%M %p IST')
