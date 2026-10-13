"""Stage S1 of md/PLAN.md, second half: the repair of the ligand terms at sulfur and gate G-L after it.

  /root/cqc-sfa-tools/md/bin/python md/ligand_fit.py            (after md/ligand_params.py; about 30 minutes)

Environment of the fit.  The DFT reference of the paper is in IEF-PCM water.  The molecular mechanics is therefore
minimised in implicit water as well (generalised Born, OBC2, dielectric 78.3553, no surface term): fitting a vacuum
model to energies in a dielectric would put the screening by the solvent into the torsion terms, and explicit water
would then count it a second time.  Vacuum results are computed and reported beside.

What is repaired, in this order, and looped until nothing changes (all terms carry the sulfonamide nitrogen, so the
shared sulfonyl terms o-sy-o, ca-sy-o, ca-sy and o-sy stay GAFF2 for both species):
 A. Rotation about S-N.  GAFF2 has no term for hn-N-sy-X of these nitrogen types; parmchk2 substitutes one by
    analogy (anion: penalty 160.5).  It is replaced by a two-term Fourier series on H-N-S-C4 (n = 1, 2) and an
    explicit zero on H-N-S-O, fitted by linear least squares to the DFT single points that the paper holds for
    rigidly rotated geometries (molecules/sfa/results/<case>/dock_strain.json: the DFT minimum with the C-S, S-N
    and C-N(aniline) rotors turned, 16 geometries for the anion, 20 for the neutral molecule).  The MM energies are
    single points on the same geometries.  One further coefficient absorbs the error of the aniline rotor in those
    geometries (cos 2 on ca-ca-nv-hn); it is reported and NOT put into the repaired force field.
 B. Equilibrium values of the S-N bond and of the angles C-S-N, O-S-N and S-N-H (and H-N-H), each moved until the
    MM minimum of the DFT conformer reproduces the DFT value.  Force constants stay GAFF2 (no Hessian is stored).
 C. Rotation about C4-S.  The relaxed scan is run with the protocol of the DFT scan; only if its largest deviation
    exceeds 0.5 kcal/mol is ca-ca-sy-N refitted (n = 2, then n = 2 and 4), together with A, against the relaxed
    scan (sigma 0.1 kcal/mol) and the rigid single points (sigma 0.5 kcal/mol), the scan being relaxed again after
    every fit.
The result is written as a frcmod (type order as in gaff2.dat, so that it replaces the native entries in tleap and
in MCPB.py alike), the topology is rebuilt from it with tleap, and every number of the "after" tables comes from
that rebuilt topology, not from the fitting model in memory.
"""
import json, math, pathlib, shutil
import numpy as np
import ligand_common as lc
from ligand_common import OUT, BULK, RESULTS, ROOT, SPECIES
import ligand_eval as le
import ligand_params as lp

ENV_FIT = 'OBC2'
SIGMA_RIGID, SIGMA_SCAN = 0.5, 0.1          # kcal/mol, fixed before fitting
NTYPE = {'anion': 'n2', 'neutral': 'n8'}


def rigid_points(sp):
    """DFT single points on rigidly rotated geometries (paper, burst 2).  Energies relative to the single point
    at opt.xyz, kcal/mol.  The first entry is opt.xyz itself."""
    case = SPECIES[sp]['case']
    d = json.loads((RESULTS / case / 'dock_strain.json').read_text())
    rel = d['relative']['b3lyp-d3bj']['vs_first_kcal']
    el, x0, _ = lc.read_xyz(RESULTS / case / 'opt.xyz')
    names = lc.names_for(sp, el, x0)
    pts = [dict(name='opt', x=x0, E=0.0, file=f'molecules/sfa/results/{case}/opt.xyz')]
    for g in d['geometries'][1:]:
        el, x, _ = lc.read_xyz(ROOT / g['file'])
        lc.names_for(sp, el, x)
        # same frame: every bond length of the heavy-atom skeleton as in opt.xyz
        ix = {n: i for i, n in enumerate(names)}
        for a, b in lc.BONDS_HEAVY:
            assert abs(lc.dist(x, ix[a], ix[b]) - lc.dist(x0, ix[a], ix[b])) < 2e-3, (g['name'], a, b)
        pts.append(dict(name=g['name'], x=x, E=rel[g['name']], file=g['file']))
    return names, pts, d['provenance']


class Model:
    """The MM model being repaired: an MM object with extra torsion slots whose coefficients can be set."""

    def __init__(self, sp, prmtop, implicit):
        self.sp, self.nt = sp, NTYPE[sp]
        self.mm = mm = lc.MM(prmtop, implicit=implicit)
        ix = mm.ix
        self.hs = [ix[h] for h in ('HN21', 'HN22') if h in ix]
        S, N2, C4 = ix['S'], ix['N2'], ix['C4']
        tf = mm.torsion_force
        self.orig_sn, self.orig_cs = [], []
        for (t, i, j, k, l, n, ph, kk) in mm.proper_torsions():
            if {j, k} == {S, N2}:
                self.orig_sn.append((t, i, j, k, l, n, ph, kk))
            if {j, k} == {C4, S} and N2 in (i, l):
                self.orig_cs.append((t, i, j, k, l, n, ph, kk))
        self.groups = {
            'rotor': dict(dihedrals=[(C4, S, N2, h) for h in self.hs], ns=(1, 2)),
            'cs': dict(dihedrals=[(ix['C3'], C4, S, N2), (ix['C5'], C4, S, N2)], ns=(2, 4)),
        }
        self.aniline = [(ix[c], ix['C1'], ix['N1'], ix[h]) for c in ('C2', 'C6') for h in ('HN11', 'HN12')]
        self.slots = {}
        for g, d in self.groups.items():
            for dih in d['dihedrals']:
                for n in d['ns']:
                    self.slots[(g, dih, n)] = tf.addTorsion(*dih, n, 0.0, 0.0)
        mm.rebuild()
        self.coef = {g: {n: 0.0 for n in d['ns']} for g, d in self.groups.items()}
        self.cs_replaced = False
        self.eq = {}                       # (kind, types) -> new equilibrium value

    def _apply(self):
        tf = self.mm.torsion_force
        for (g, dih, n), t in self.slots.items():
            c = self.coef[g][n]
            tf.setTorsionParameters(t, *dih, n, 0.0 if c >= 0 else math.pi, abs(c) * lc.KCAL)
        tf.updateParametersInContext(self.mm.context)

    def zero_original(self, which):
        """Remove the as-assigned terms that a fitted group replaces."""
        tf = self.mm.torsion_force
        for (t, i, j, k, l, n, ph, kk) in (self.orig_sn if which == 'rotor' else self.orig_cs):
            tf.setTorsionParameters(t, i, j, k, l, n, math.radians(ph), 0.0)
        tf.updateParametersInContext(self.mm.context)
        if which == 'cs':
            self.cs_replaced = True

    def set_coef(self, g, coef):
        self.coef[g].update(coef)
        self._apply()

    def basis(self, x, g, n):
        return sum(math.cos(n * math.radians(lc.dihedral(x, *d))) for d in self.groups[g]['dihedrals'])

    def aniline_basis(self, x):
        return sum(math.cos(2 * math.radians(lc.dihedral(x, *d))) for d in self.aniline)

    def e_rest(self, xs, fit_groups):
        """Energies with the coefficients of fit_groups set to zero."""
        keep = {g: dict(self.coef[g]) for g in fit_groups}
        for g in fit_groups:
            self.set_coef(g, {n: 0.0 for n in self.groups[g]['ns']})
        e = np.array([self.mm.energy(x) for x in xs])
        for g in fit_groups:
            self.set_coef(g, keep[g])
        return e


def fit_torsions(model, rigid, fit_cs, cs_ns, scan_rows=None):
    """Weighted linear least squares for the Fourier coefficients (cos n phi) of the rotor group and, if asked,
    of the C-S group.  Returns the coefficients, the nuisance aniline coefficient and fit statistics."""
    groups = ['rotor'] + (['cs'] if fit_cs else [])
    cols = [('rotor', n) for n in model.groups['rotor']['ns']] + ([('cs', n) for n in cs_ns] if fit_cs else [])
    xs = [p['x'] for p in rigid]
    er = model.e_rest(xs, groups)
    y_r = np.array([p['E'] for p in rigid])[1:] - (er[1:] - er[0])
    B = np.array([[model.basis(x, g, n) for (g, n) in cols] for x in xs])
    A_r = np.hstack([(B - B[0])[1:], np.array([[model.aniline_basis(x) - model.aniline_basis(xs[0])] for x in xs[1:]])])
    rows, rhs = [np.hstack([A_r, np.zeros((len(y_r), 1))]) / SIGMA_RIGID], [y_r / SIGMA_RIGID]
    if fit_cs:
        xs_s = [r['x'] for r in scan_rows]
        es = model.e_rest(xs_s, groups)
        y_s = np.array([r['E_dft_rel'] for r in scan_rows]) - es
        A_s = np.hstack([np.array([[model.basis(x, g, n) for (g, n) in cols] for x in xs_s]),
                         np.zeros((len(xs_s), 1)), np.ones((len(xs_s), 1))])
        rows.append(A_s / SIGMA_SCAN); rhs.append(y_s / SIGMA_SCAN)
    A, y = np.vstack(rows), np.concatenate(rhs)
    if not fit_cs:
        A = A[:, :-1]
    c, res, rank, sv = np.linalg.lstsq(A, y, rcond=None)
    nr = len(y_r)
    pred_r = A_r @ c[:A_r.shape[1]]
    # leave-one-out on the rigid points
    loo = []
    for k in range(nr):
        m = np.ones(len(y), bool); m[k] = False
        ck = np.linalg.lstsq(A[m], y[m], rcond=None)[0]
        loo.append(float(A_r[k] @ ck[:A_r.shape[1]] - y_r[k]))
    out = dict(columns=[f'{g}:cos{n}' for g, n in cols] + ['aniline:cos2'] + (['scan offset'] if fit_cs else []),
               coefficients=[float(v) for v in c], rank=int(rank), condition=float(sv[0] / sv[-1]),
               rigid_rms=float(np.sqrt(((pred_r - y_r) ** 2).mean())), rigid_max=float(np.abs(pred_r - y_r).max()),
               rigid_loo_rms=float(np.sqrt(np.mean(np.square(loo)))), n_rigid=nr,
               n_scan=len(scan_rows) if fit_cs else 0)
    coef = {'rotor': {}, 'cs': {}}
    for (g, n), v in zip(cols, c):
        coef[g][n] = float(v)
    out['aniline_cos2'] = float(c[len(cols)])
    return coef, out


def targets_for(sp, m_dft):
    nt = NTYPE[sp]
    t = [dict(kind='bond', types=(nt, 'sy'), key='d_S_N2', target=m_dft['d_S_N2'], tol=1e-4),
         dict(kind='angle', types=('ca', 'sy', nt), key='a_C4_S_N2', target=m_dft['a_C4_S_N2'], tol=0.01),
         dict(kind='angle', types=(nt, 'sy', 'o'), key='a_O_S_N2_mean', target=m_dft['a_O_S_N2_mean'], tol=0.01),
         dict(kind='angle', types=('hn', nt, 'sy'), key='a_S_N2_H_mean', target=m_dft['a_S_N2_H_mean'], tol=0.01)]
    if sp == 'neutral':
        t.append(dict(kind='angle', types=('hn', nt, 'hn'), key='a_H_N2_H', target=m_dft['a_H_N2_H'], tol=0.01))
    return t


def more_measures(x, names):
    m = lc.add_means(lc.measures(x, names))
    hs = [k for k in ('a_S_N2_HN21', 'a_S_N2_HN22') if k in m]
    m['a_S_N2_H_mean'] = float(np.mean([m[k] for k in hs]))
    ix = {n: i for i, n in enumerate(names)}
    if 'HN22' in ix:
        m['a_H_N2_H'] = lc.angle(x, ix['HN21'], ix['N2'], ix['HN22'])
    return m


def current_eq(mm, t):
    if t['kind'] == 'bond':
        k = mm.bonds_of(*t['types'])[0]
        return mm.bond_force.getBondParameters(k)[2].value_in_unit(mm._unit.angstrom)
    k = mm.angles_of(*t['types'])[0]
    return mm.angle_force.getAngleParameters(k)[3].value_in_unit(mm._unit.degree)


def set_eq(mm, t, v):
    if t['kind'] == 'bond':
        mm.set_bond_eq(*t['types'], v)
    else:
        mm.set_angle_eq(*t['types'], v)


def fit_equilibria(model, targets, x_start, names, log, maxit=60):
    """Move the equilibrium values until the minimum (from x_start) has the DFT values."""
    mm = model.mm
    x = x_start
    for it in range(maxit):
        r = mm.minimise(x)
        x = r['x']
        m = more_measures(x, names)
        worst = 0.0
        for t in targets:
            d = t['target'] - m[t['key']]
            worst = max(worst, abs(d) / t['tol'])
            if abs(d) > t['tol']:
                v = current_eq(mm, t) + d
                set_eq(mm, t, v)
                model.eq[(t['kind'], t['types'])] = v
        if worst <= 1.0:
            break
    log(f'    equilibrium values converged after {it + 1} minimisations (worst miss {worst:.2f} tolerances)')
    return x, m


def repair(sp, ref, log):
    names = ref[sp]['names']
    x_dft = ref[sp]['x']
    m_dft = more_measures(x_dft, names)
    rigid_names, rigid, prov = rigid_points(sp)
    assert rigid_names == names
    model = Model(sp, OUT / f'sfa_{sp}_gaff2.prmtop', ENV_FIT)
    mm = model.mm
    targets = targets_for(sp, m_dft)
    old_eq = {(t['kind'], t['types']): current_eq(mm, t) for t in targets}
    history = []
    # A (first pass): rotor against the rigid single points
    model.zero_original('rotor')
    fit_cs, cs_ns, scan_rows = False, (), None
    x_min = x_dft
    pruned = None
    for cycle in range(16):
        coef, stat = fit_torsions(model, rigid, fit_cs, cs_ns, scan_rows)
        change = max([abs(coef[g][n] - model.coef[g][n]) for g in coef for n in coef[g]] or [0.0])
        model.set_coef('rotor', coef['rotor'])
        if fit_cs:
            model.set_coef('cs', coef['cs'])
        x_min, m_min = fit_equilibria(model, targets, x_min, names, log)
        scan_rows = le.scan(mm, sp, x_min, names, global_search=False)
        ssum = le.scan_summary(scan_rows)
        history.append(dict(cycle=cycle + 1, fit_cs=fit_cs, cs_terms=list(cs_ns), fit=stat, coef_change=change,
                            scan_max_abs_dev=ssum['max_abs_dev'],
                            tau_H=[m_min.get('tau_C4_S_N2_HN21'), m_min.get('tau_C4_S_N2_HN22')]))
        log(f"    cycle {cycle + 1}: rotor {[round(coef['rotor'][n], 3) for n in (1, 2)]}"
            + (f", C-S {[round(coef['cs'][n], 3) for n in cs_ns]}" if fit_cs else '')
            + f"; rigid rms {stat['rigid_rms']:.2f} (LOO {stat['rigid_loo_rms']:.2f}); scan max|dev| "
              f"{ssum['max_abs_dev']:.3f}; coefficient change {change:.4f}")
        if ssum['max_abs_dev'] > 0.5 and (not fit_cs or (change < 2e-3 and cs_ns == (2,))):
            # C: the C-S torsion misses -> refit it (first n = 2, then n = 2 and 4)
            cs_ns = (2,) if not fit_cs else (2, 4)
            if not fit_cs:
                model.zero_original('cs')
                scan_rows = le.scan(mm, sp, x_min, names, global_search=False)
            fit_cs = True
            log(f'    scan outside 0.5 kcal/mol: refitting ca-ca-sy-{NTYPE[sp]} with n = {cs_ns}')
            continue
        if change < 2e-3 and cycle > 0:
            if pruned is None:
                # keep the number of changed terms small: an equilibrium value that moved by less than 1 degree
                # (0.005 A) goes back to GAFF2 if the DFT value is then still met within 0.5 degree (0.005 A)
                pruned = []
                for t in list(targets):
                    key = (t['kind'], t['types'])
                    lim = 1.0 if t['kind'] == 'angle' else 0.005
                    tol = 0.5 if t['kind'] == 'angle' else 0.005
                    v_new = current_eq(mm, t)
                    if abs(v_new - old_eq[key]) >= lim:
                        continue
                    set_eq(mm, t, old_eq[key])
                    rest = [u for u in targets if u is not t]
                    x2, m2 = fit_equilibria(model, rest, x_min, names, log)
                    miss = m2[t['key']] - t['target']
                    if abs(miss) <= tol:
                        targets, x_min = rest, x2
                        model.eq.pop(key, None)
                        pruned.append(dict(kind=t['kind'], types='-'.join(t['types']), fitted_value=v_new,
                                           gaff2_value=old_eq[key], miss_with_gaff2_value=miss))
                        log(f"    {'-'.join(t['types'])} goes back to GAFF2 ({old_eq[key]:.3f}; fitted {v_new:.3f}; "
                            f"miss {miss:+.3f})")
                    else:
                        set_eq(mm, t, v_new)
                        x_min, _ = fit_equilibria(model, targets, x_min, names, log)
                if pruned:
                    continue
            break
    new_eq = {(t['kind'], t['types']): current_eq(mm, t) for t in targets}
    return dict(model=model, targets=targets, old_eq=old_eq, new_eq=new_eq, history=history, fit=stat,
                fit_cs=fit_cs, cs_ns=cs_ns, x_min=x_min, rigid=rigid, rigid_provenance=prov, m_dft=m_dft,
                pruned=pruned or [], all_targets=targets_for(sp, m_dft))


def fourier_lines(types, coef):
    """frcmod DIHE lines of a Fourier series given as coefficients of cos(n phi): PK = |a|, phase 0 or 180."""
    terms = [(n, a) for n, a in sorted(coef.items())]
    lines = []
    for k, (n, a) in enumerate(terms):
        last = k == len(terms) - 1
        lines.append(f"{types}   1  {abs(a):8.4f}    {0.0 if a >= 0 else 180.0:8.3f}    {n if last else -n:6.1f}")
    return lines


def write_repair_frcmod(sp, rep, gaff2, path):
    nt = NTYPE[sp]
    model = rep['model']

    def order(kind, types):                 # the order in which gaff2.dat writes the native term
        key, rec = lc.lookup(gaff2[kind], types)
        return key if key else types

    def name(t):
        return '-'.join(f'{a:<2s}' for a in t)
    L = [f'Repair of GAFF2 terms at sulfur for the {SPECIES[sp]["label"]} of sulfanilamide (md/ligand_fit.py; '
         f'load after leaprc.gaff2 and sfa_{sp}_gaff2.frcmod)', 'MASS', '', 'BOND']
    for t in rep['targets']:
        if t['kind'] == 'bond':
            key = order('bond', t['types']); k0 = lc.lookup(gaff2['bond'], t['types'])[1]['k']
            L.append(f"{name(key)}  {k0:7.2f}   {rep['new_eq'][('bond', t['types'])]:7.4f}       "
                     f"r0 was {rep['old_eq'][('bond', t['types'])]:.4f}; force constant GAFF2")
    L += ['', 'ANGLE']
    for t in rep['targets']:
        if t['kind'] == 'angle':
            key = order('angle', t['types']); k0 = lc.lookup(gaff2['angle'], t['types'])[1]['k']
            L.append(f"{name(key)}  {k0:8.3f}   {rep['new_eq'][('angle', t['types'])]:9.3f}       "
                     f"theta0 was {rep['old_eq'][('angle', t['types'])]:.3f}; force constant GAFF2")
    L += ['', 'DIHE']
    c = model.coef['rotor']
    lines = fourier_lines(name(('hn', nt, 'sy', 'ca')), c)
    lines[0] += '      rotation about S-N, fitted to DFT single points; replaces the term by analogy'
    L += lines
    L.append(f"{name(('hn', nt, 'sy', 'o'))}   1  {0.0:8.4f}    {0.0:8.3f}    {2.0:6.1f}"
             '      set to zero (the series on H-N-S-C4 carries the rotor)')
    if rep['fit_cs']:
        cs = {n: model.coef['cs'][n] for n in rep['cs_ns']}
        lines = fourier_lines(name(('ca', 'ca', 'sy', nt)), cs)
        lines[0] += '      rotation about C4-S, fitted to the relaxed DFT scan; replaces X -ca-sy-X for this path'
        L += lines
    L += ['', 'IMPROPER', '', 'NONBON', '', '']
    pathlib.Path(path).write_text('\n'.join(L))


def rebuild(sp):
    work = BULK / f'repaired_{sp}'
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    for f in (f'sfa_{sp}.mol2', f'sfa_{sp}_gaff2.frcmod', f'sfa_{sp}_repair.frcmod'):
        shutil.copy(OUT / f, work / f)
    stem = f'sfa_{sp}_repaired'
    (work / f'{stem}.leap.in').write_text(
        f'source leaprc.gaff2\nloadamberparams sfa_{sp}_gaff2.frcmod\nloadamberparams sfa_{sp}_repair.frcmod\n'
        f'SAN = loadmol2 sfa_{sp}.mol2\ncheck SAN\nsaveamberparm SAN {stem}.prmtop {stem}.inpcrd\nquit\n')
    out = lc.run(['tleap', '-f', f'{stem}.leap.in'], work, log=work / f'{stem}.leap.out')
    import re
    m = re.search(r'Errors = (\d+)', out)
    assert m and int(m.group(1)) == 0 and (work / f'{stem}.prmtop').exists(), out[-1500:]
    for ext in ('prmtop', 'inpcrd'):
        shutil.copy(work / f'{stem}.{ext}', OUT / f'{stem}.{ext}')
    return OUT / f'{stem}.prmtop'


def charge_sensitivity(sp, ref):
    """AM1-BCC gives the two sulfonyl oxygens of the anion different charges, because antechamber draws one S=O
    and one S-O(-).  Check of the repaired terms with the two charges replaced by their mean: minimum from the DFT
    geometry and relaxed scan, implicit water.  The topology of this check is not delivered."""
    work = BULK / f'charge_check_{sp}'
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    lines = (OUT / f'sfa_{sp}.mol2').read_text().split('\n')
    a0, a1 = lines.index('@<TRIPOS>ATOM') + 1, lines.index('@<TRIPOS>BOND')
    q = {l.split()[1]: float(l.split()[8]) for l in lines[a0:a1]}
    qm = round(0.5 * (q['O1'] + q['O2']), 6)
    for k in range(a0, a1):
        f = lines[k].split()
        if f[1] in ('O1', 'O2'):
            lines[k] = lines[k][:lines[k].rindex(f[8])] + f'{qm:9.6f}'.strip().rjust(len(f[8]))
    (work / 'sym.mol2').write_text('\n'.join(lines))
    for f in (f'sfa_{sp}_gaff2.frcmod', f'sfa_{sp}_repair.frcmod'):
        shutil.copy(OUT / f, work / f)
    (work / 'leap.in').write_text(
        f'source leaprc.gaff2\nloadamberparams sfa_{sp}_gaff2.frcmod\nloadamberparams sfa_{sp}_repair.frcmod\n'
        'SAN = loadmol2 sym.mol2\nsaveamberparm SAN sym.prmtop sym.inpcrd\nquit\n')
    lc.run(['tleap', '-f', 'leap.in'], work, log=work / 'leap.out')
    mm = lc.MM(work / 'sym.prmtop', implicit=ENV_FIT)
    names = ref[sp]['names']
    r = mm.minimise(ref[sp]['x'])
    m = more_measures(r['x'], names)
    rows = le.scan(mm, sp, r['x'], names, global_search=False)
    ss = le.scan_summary(rows)
    return dict(what='O1 and O2 charges replaced by their mean; repaired terms unchanged; implicit water',
                charges_as_assigned=dict(O1=q['O1'], O2=q['O2']), charge_used=qm,
                min_from_dft=dict(a_C4_S_N2=m['a_C4_S_N2'], d_S_N2=m['d_S_N2'], a_O_S_N2_mean=m['a_O_S_N2_mean'],
                                  tau_fold=m['tau_fold'], tau_C4_S_N2_HN21=m.get('tau_C4_S_N2_HN21')),
                scan_max_abs_dev=ss['max_abs_dev'])


def rigid_table(sp, rigid, names, prmtops):
    """MM single points on the rigid DFT geometries, before and after, in both environments."""
    ix = {n: i for i, n in enumerate(names)}
    rows = []
    cols = {}
    for tag, prm in prmtops.items():
        for envn in ('vacuum', 'OBC2'):
            mm = lc.MM(prm, implicit=None if envn == 'vacuum' else envn)
            e = np.array([mm.energy(p['x']) for p in rigid])
            cols[f'mm_{tag}_{envn}'] = e - e[0]
    t0 = lc.dihedral(rigid[0]['x'], ix['C2'], ix['C1'], ix['N1'], ix['HN11'])
    for k, p in enumerate(rigid):
        m = lc.measures(p['x'], names)
        r = dict(species=sp, geometry=p['name'], tau_C3_C4_S_N2=m['tau_C3'],
                 tau_C4_S_N2_HN21=m.get('tau_C4_S_N2_HN21'), tau_C4_S_N2_HN22=m.get('tau_C4_S_N2_HN22'),
                 aniline_rotation=lc.wrap(lc.dihedral(p['x'], ix['C2'], ix['C1'], ix['N1'], ix['HN11']) - t0),
                 E_dft=p['E'])
        for c, v in cols.items():
            r[c] = float(v[k])
        rows.append(r)
    stats = {}
    d = np.array([p['E'] for p in rigid])
    clean = np.array([abs(r['aniline_rotation']) < 5 for r in rows])
    for c, v in cols.items():
        stats[c] = dict(rms_all=float(np.sqrt(((v - d)[1:] ** 2).mean())),
                        rms_aniline_unrotated=float(np.sqrt(((v - d)[clean][1:] ** 2).mean())),
                        max_abs_aniline_unrotated=float(np.abs((v - d)[clean][1:]).max()),
                        n_all=int(len(d) - 1), n_aniline_unrotated=int(clean.sum() - 1))
    return rows, stats


def main():
    print(f'ligand_fit.py started {lc.now_ist()}')
    ref = lp.reference_data()
    gaff2 = lc.read_parm_dat(lc.GAFF2_DAT)
    before = json.loads((OUT / 'gaff2_as_assigned.json').read_text())
    summary = dict(script='md/ligand_fit.py', started=lc.now_ist(), versions=lc.versions(), seed=le.SEED,
                   environment_of_fit='implicit water, OBC2 generalised Born, solvent dielectric 78.3553, solute '
                                      'dielectric 1, no surface-area term (OpenMM CustomGBForce built by ParmEd; '
                                      'radii of the topology, tleap default mbondi)',
                   sigma_rigid=SIGMA_RIGID, sigma_scan=SIGMA_SCAN, species={})
    geom_rows, min_rows, rigid_rows = [], [], []
    for sp in SPECIES:
        print(f'== {sp}')
        rep = repair(sp, ref, print)
        model = rep['model']
        write_repair_frcmod(sp, rep, gaff2, OUT / f'sfa_{sp}_repair.frcmod')
        prm = rebuild(sp)
        # the rebuilt topology must be the fitted model
        chk = lc.MM(prm, implicit=ENV_FIT)
        xs = [p['x'] for p in rep['rigid']] + [rep['x_min']]
        e1 = np.array([model.mm.energy(x) for x in xs]); e2 = np.array([chk.energy(x) for x in xs])
        d12 = float(np.abs((e1 - e1[0]) - (e2 - e2[0])).max())
        print(f'  rebuilt topology against the fitted model: largest difference of relative energies {d12:.4f} kcal/mol')
        assert d12 < 0.02, d12
        rows, overridden = lp.term_tables(prm, OUT / f'sfa_{sp}_gaff2.frcmod', gaff2,
                                          repair_path=OUT / f'sfa_{sp}_repair.frcmod')
        uniq = lp.unique_terms(rows)
        lc.write_csv(OUT / f'terms_repaired_{sp}.csv', uniq,
                     ['kind', 'types', 'k', 'eq', 'n', 'phase', 'instances', 'atoms', 'source'])
        res = le.evaluate(prm, sp, ref[sp]['x'], ref[sp]['names'], ref[sp]['m'], log=print)
        geom_rows += lp.geometry_table(sp, ref, res)
        min_rows += lp.minima_table(sp, res)
        lc.write_csv(OUT / f'scan_after_{sp}.csv', le.scan_table(res['OBC2']['scan']))
        lc.write_csv(OUT / f'scan_after_{sp}_vacuum.csv', le.scan_table(res['vacuum']['scan']))
        rrows, rstats = rigid_table(sp, rep['rigid'], ref[sp]['names'],
                                    dict(gaff2=OUT / f'sfa_{sp}_gaff2.prmtop', repaired=prm))
        rigid_rows += rrows
        changed = []
        for t in rep['targets']:
            key = (t['kind'], t['types'])
            native = lc.lookup(gaff2[t['kind']], t['types'])[1]
            changed.append(dict(kind=t['kind'], types='-'.join(t['types']), quantity=t['key'], dft_target=t['target'],
                                old=rep['old_eq'][key], new=rep['new_eq'][key], force_constant=native['k'],
                                involves_ligating_N=True))
        tors = [dict(kind='torsion', types=f'hn-{NTYPE[sp]}-sy-ca', path='H-N2-S-C4',
                     old='by analogy, see terms_gaff2 table',
                     new=[dict(n=n, coefficient_of_cos=model.coef['rotor'][n], PK=abs(model.coef['rotor'][n]),
                               phase=0.0 if model.coef['rotor'][n] >= 0 else 180.0) for n in (1, 2)],
                     involves_ligating_N=True),
                dict(kind='torsion', types=f'hn-{NTYPE[sp]}-sy-o', path='H-N2-S-O1, H-N2-S-O2',
                     old='by analogy, see terms_gaff2 table', new='zero', involves_ligating_N=True)]
        if rep['fit_cs']:
            tors.append(dict(kind='torsion', types=f'ca-ca-sy-{NTYPE[sp]}', path='C3-C4-S-N2, C5-C4-S-N2',
                             old='X -ca-sy-X: PK/IDIVF 1.3, phase 180, n 2',
                             new=[dict(n=n, coefficient_of_cos=model.coef['cs'][n], PK=abs(model.coef['cs'][n]),
                                       phase=0.0 if model.coef['cs'][n] >= 0 else 180.0) for n in rep['cs_ns']],
                             involves_ligating_N=True))
        summary['species'][sp] = dict(
            changed_equilibrium_values=changed, changed_torsions=tors, fit_history=rep['history'],
            final_fit=rep['fit'], aniline_rotor_note=dict(
                coefficient_of_cos2_per_dihedral=rep['fit']['aniline_cos2'],
                note='correction that the rigid single points ask for on the four ca-ca-nv-hn dihedrals '
                     '(GAFF2 by analogy: PK/IDIVF 1.05, phase 180, n 2); not part of the repaired force field'),
            rebuilt_vs_fitted_max_diff=d12, native_terms_differing_from_gaff2=overridden,
            rigid_single_points=rstats, rigid_provenance=rep['rigid_provenance'],
            results={e: dict(n_starts=r['n_starts'], n_minima=len(r['minima']), global_min=r['global_min']['m'],
                             global_min_xyz=np.round(r['global_min']['x'], 5).tolist(),
                             from_dft_min=dict(E_rel=r['from_dft_min']['E_rel'], m=r['from_dft_min']['m']),
                             hessian_lowest7=r['hessian_lowest7'], scan_summary=r['scan_summary'], gate=r['gate'],
                             gate_from_dft_min=r['gate_from_dft_min']) for e, r in res.items()},
            gate_before={e: before['species'][sp]['results'][e]['gate'] for e in ('vacuum', 'OBC2')},
            targets_back_to_gaff2=rep['pruned'],
            charge_sensitivity=charge_sensitivity(sp, ref) if sp == 'anion' else None)
        if sp == 'anion':
            cs_ = summary['species'][sp]['charge_sensitivity']
            print(f"  equal oxygen charges: C-S-N {cs_['min_from_dft']['a_C4_S_N2']:.2f}, S-N "
                  f"{cs_['min_from_dft']['d_S_N2']:.4f}, scan max|dev| {cs_['scan_max_abs_dev']:.3f}")
        # final frame of the repaired model (global minimum in implicit water), for the next stage
        g = res['OBC2']['global_min']
        with open(OUT / f'sfa_{sp}_repaired_min.xyz', 'w') as f:
            f.write(f"{len(ref[sp]['names'])}\n{SPECIES[sp]['label']}, repaired GAFF2, minimum in implicit water (OBC2)\n")
            for n, xyz in zip(ref[sp]['names'], g['x']):
                f.write(f"{n[0]:2s} {xyz[0]:12.6f} {xyz[1]:12.6f} {xyz[2]:12.6f}\n")
    lc.write_csv(OUT / 'geometry_after.csv', geom_rows)
    lc.write_csv(OUT / 'minima_after.csv', min_rows)
    lc.write_csv(OUT / 'rigid_single_points.csv', rigid_rows)
    summary['finished'] = lc.now_ist()
    (OUT / 'repair.json').write_text(json.dumps(le.strip(summary), indent=1) + '\n')
    print(f'ligand_fit.py finished {lc.now_ist()}')


if __name__ == '__main__':
    main()
