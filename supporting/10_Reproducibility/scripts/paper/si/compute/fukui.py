#!/usr/bin/env python3
"""
Condensed Fukui functions f+, f-, f0 (and the dual descriptor) by finite differences of atomic populations of the
N, N+1 and N-1 electron systems at one fixed geometry, for the SI (sulfanilamide_pcm_d3 and cation_pcm_d3).

Usage (from the repository root, pyscf env):
    nice -n 10 /root/cqc-sfa-tools/env/bin/python paper/si/compute/fukui.py <case> --threads 2
    ... --basis sto-3g --grid 1 --outdir <scratch>        (code test)
Writes molecules/sfa/results/<case>/fukui/{fukui.json, fukui.csv, fukui.md, LOG.md, N.chk, Nm1.chk, Np1.chk}.

Why this script lives under paper/si/compute/ and not dft/: the SI-completion task was allowed to write only under
paper/ and molecules/sfa/results/<case>/fukui/.  It goes through the same code path as dft/dscf.py (props_common:
thread pinning, scf_for -> common.make_mf, converge), so the level of theory is the case's own.

Method
  * Geometry: results/<case>/opt.xyz (the burst-1a IEF-PCM optimum).  B3LYP-D3(BJ)/6-311++G(d,p), density fitting,
    grid level 4, IEF-PCM water (eps 78.3553, Lebedev 29), conv_tol 1e-9: common.make_mf via props_common.scf_for.
  * N: RKS.  N-1 and N+1: UKS doublets at the same geometry, initial guess from the N orbitals with the beta HOMO
    emptied (N-1) or an alpha electron added to the LUMO (N+1), exactly as dft/dscf.py.  Each state is solvated
    at equilibrium (as in dscf.py).  The total energies are compared with dscf.json: the same state reached twice
    gives the same energy to ~1e-7 Ha, and dscf.py ran the UKS internal stability analysis on those states (no
    instability found), so an agreeing energy also carries over that check.  --stability repeats it here.
  * Populations: NPA from PySCF natural atomic orbitals (pyscf.lo.nao, Reed-Weinstock-Weinhold 1985), computed
    for each state from its own total density, as the NPA step of dft/mwfn_extra.py does; and Hirshfeld, from the
    molecular density on the SCF integration grid and a promolecule of spherically averaged neutral free atoms
    computed with the same functional and basis (pyscf.scf.atom_ks, fractional-occupation spherical average).
    The promolecule is the same for all three electron counts (the usual convention).
  * Condensed Fukui functions (q = atomic charge, Yang & Mortier 1986):
        f+_A = q_A(N) - q_A(N+1)     (site for nucleophilic attack, i.e. where an added electron goes)
        f-_A = q_A(N-1) - q_A(N)     (site for electrophilic attack)
        f0_A = [q_A(N-1) - q_A(N+1)] / 2
        dual descriptor Df_A = f+_A - f-_A   (Morell, Grand, Toro-Labbe 2005)
    Each of f+, f-, f0 sums to 1 over the atoms (checked).
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3] / 'dft'))
from props_common import *          # noqa: F401,F403 -- pins threads / scratch, so it must precede numpy and PySCF

ap = parser(__doc__)
ap.add_argument('--stability', action='store_true', help='run the UKS internal stability analysis of N+-1')
args = ap.parse_args()
t_start = time.time()
case, xyz, outdir0, tag = setup(args)
if case['spin'] != 0:
    raise SystemExit('[fukui] closed-shell N-electron cases only')
outdir = pathlib.Path(args.outdir).resolve() if args.outdir else case['dir'] / 'fukui'
outdir.mkdir(parents=True, exist_ok=True)
LOG = outdir / out_name('LOG', tag, '.md')


def log(msg):
    import datetime
    line = f"- {datetime.datetime.now():%Y-%m-%d %H:%M}  {msg}"
    print('[fukui] ' + msg, flush=True)
    with open(LOG, 'a') as f:
        f.write(line + '\n')


from pyscf import scf as _scf                      # noqa: E402
from pyscf.lo import nao as nao_mod                # noqa: E402
from pyscf.scf import atom_ks                      # noqa: E402
from pyscf.dft import numint                       # noqa: E402

syms, XYZ = read_xyz(xyz)
labels = case.get('cif_labels') or atom_labels(syms)


# ---------------------------------------------------------------------------------------------------- states
def run_n():
    mol = mol_for_case(case, xyz, verbose=args.verbose)
    mf = scf_for(case, mol, args, conv_tol=1e-9, chkfile=outdir / out_name('N', tag, '.chk'))
    t = time.time()
    converge(mf, what='N-electron SCF')
    log(f'N ({mol.nelectron} e): E = {mf.e_tot:.10f} Ha, {time.time() - t:.0f} s')
    return mol, mf


def run_ion(dq, mf_n):
    """Vertical N-dq electron state (dq = +1: N-1; dq = -1: N+1) as a UKS doublet (the dscf.py recipe)."""
    mol = mol_for_case(dict(case, charge=case['charge'] + dq, spin=1), xyz, verbose=args.verbose)
    stem = 'Nm1' if dq == 1 else 'Np1'
    mf = scf_for(case, mol, args, conv_tol=1e-9, chkfile=outdir / out_name(stem, tag, '.chk'))
    mf.max_cycle = 200
    C, occ = mf_n.mo_coeff, mf_n.mo_occ
    n = int((occ > 0).sum())
    oa, ob = np.zeros(C.shape[1]), np.zeros(C.shape[1])
    if dq == +1:
        oa[:n] = 1; ob[:n - 1] = 1
    else:
        oa[:n + 1] = 1; ob[:n] = 1
    dm0 = np.array([(C * oa) @ C.T, (C * ob) @ C.T])
    what = 'N-1' if dq == 1 else 'N+1'
    t = time.time()
    converge(mf, dm0=dm0, what=what)
    n_unstable = 0
    if args.stability:
        for _ in range(3):
            mo_i, _, stable, _ = mf.stability(return_status=True)
            if stable:
                break
            n_unstable += 1
            log(f'{what}: internal instability found, restarting from the lower solution')
            converge(mf, dm0=mf.make_rdm1(mo_i, mf.mo_occ), what=what)
    s2, _ = mf.spin_square()
    log(f'{what} ({mol.nelectron} e, charge {mol.charge:+d}): E = {mf.e_tot:.10f} Ha, <S^2> = {s2:.4f}, '
        f'{time.time() - t:.0f} s')
    return mol, mf, {'charge': mol.charge, 'multiplicity': 2, 'E_Ha': float(mf.e_tot), 'S2': float(s2),
                     'stability_checked_here': bool(args.stability), 'instabilities_followed': n_unstable,
                     'converged': bool(mf.converged), 'wall_s': round(time.time() - t, 1)}


# ---------------------------------------------------------------------------------------------------- populations
def total_dm(mf, mo=None):
    dm = mf.make_rdm1(mo if mo is not None else mf.mo_coeff, mf.mo_occ)
    return dm[0] + dm[1] if np.asarray(dm).ndim == 3 else dm


def npa_charges(mol, mf):
    """NPA exactly as dft/mwfn_extra.py npa(): NAOs of this state's total density, populations = NAO diagonal."""
    C = nao_mod.nao(mol, mf)
    if np.asarray(mf.mo_occ).ndim == 1:
        mo = np.linalg.solve(C, mf.mo_coeff)
    else:
        mo = [np.linalg.solve(C, c) for c in mf.mo_coeff]
    pop = np.diag(total_dm(mf, mo))
    return np.array([mol.atom_charge(a) - pop[s0:s1].sum() for a, (_, _, s0, s1) in enumerate(mol.aoslice_by_atom())])


FREE = {}


def free_atom_radial(mol):
    """Spherically averaged neutral free-atom densities rho_Z(r) (same functional and basis), tabulated on a radial
    grid for interpolation in log(rho)."""
    r = np.concatenate([[0.0], np.geomspace(1e-5, 30.0, 4000)])          # bohr
    # 26 directions (cube faces, edges, corners): the fractional-occupation average is spherical; averaging over
    # directions removes any residual anisotropy from the d/f polarisation shells.
    dirs = np.array([(x, y, z) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1) if (x, y, z) != (0, 0, 0)],
                    float)
    dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    xc = os.environ.get('DFT_XC', XC)
    for el in sorted(set(syms)):
        if el in FREE:
            continue
        atm = gto.M(atom=[[el, (0.0, 0.0, 0.0)]], basis=mol.basis if isinstance(mol.basis, str) else mol.basis,
                    charge=0, spin=None, verbose=0, cart=False)
        a = atom_ks.AtomSphericAverageRKS(atm)
        a.xc = xc
        a.grids.atom_grid = (150, 590)
        a.verbose = 0
        a.run()
        dm = (a.mo_coeff * a.mo_occ) @ a.mo_coeff.T
        pts = (r[:, None, None] * dirs[None]).reshape(-1, 3)
        ao = numint.eval_ao(atm, pts)
        rho = numint.eval_rho(atm, ao, dm).reshape(len(r), len(dirs)).mean(axis=1)
        aniso = float(np.max(np.abs(numint.eval_rho(atm, ao, dm).reshape(len(r), len(dirs)).std(axis=1)
                                    / np.maximum(rho, 1e-30))[r < 5]))
        # normalisation check: 4 pi int r^2 rho dr = Z
        nel = float(np.trapezoid(4 * np.pi * r ** 2 * rho, r))
        FREE[el] = {'r': r, 'logrho': np.log(np.maximum(rho, 1e-300)), 'nelec_int': nel, 'Z': int(atm.atom_charge(0)),
                    'E_Ha': float(a.e_tot), 'max_rel_anisotropy_r<5': aniso}
        log(f'free atom {el}: E = {a.e_tot:.8f} Ha, 4pi int r^2 rho = {nel:.5f} (Z = {atm.atom_charge(0)}), '
            f'anisotropy {aniso:.1e}')


def promolecule_weights(coords, mol):
    R = mol.atom_coords()                                   # bohr
    w = np.empty((mol.natm, coords.shape[0]))
    for A in range(mol.natm):
        f = FREE[syms[A]]
        d = np.linalg.norm(coords - R[A], axis=1)
        w[A] = np.exp(np.interp(d, f['r'], f['logrho'], right=-700.0))
    tot = w.sum(axis=0)
    return w / np.maximum(tot, 1e-300)


def hirshfeld_charges(mol, grids, dms):
    """Hirshfeld charges of several densities on one molecular grid (same geometry, same promolecule)."""
    ni = numint.NumInt()
    pops = np.zeros((len(dms), mol.natm)); nint = np.zeros(len(dms))
    for ao, mask, weight, coords in ni.block_loop(mol, grids, mol.nao, deriv=0, max_memory=2000):
        w = promolecule_weights(coords, mol)
        for k, dm in enumerate(dms):
            rho = ni.eval_rho(mol, ao, dm, mask, 'LDA')
            pops[k] += w @ (rho * weight)
            nint[k] += float(np.dot(rho, weight))
    Z = np.array([mol.atom_charge(a) for a in range(mol.natm)], float)
    return [Z - p for p in pops], nint


# ---------------------------------------------------------------------------------------------------- run
log(f'start: case {case["name"]}, geometry {xyz.relative_to(ROOT) if str(xyz).startswith(str(ROOT)) else xyz}, '
    f'threads {os.environ.get("DFT_THREADS")}, basis {basis_description(case)}, grid {grid_level(args)}, '
    f'{environment_label(case)}')
mol_n, mf_n = run_n()
mol_m, mf_m, st_m = run_ion(+1, mf_n)
mol_p, mf_p, st_p = run_ion(-1, mf_n)

t = time.time()
q = {'NPA': {}, 'Hirshfeld': {}}
for k, (mol, mf) in {'N': (mol_n, mf_n), 'N-1': (mol_m, mf_m), 'N+1': (mol_p, mf_p)}.items():
    q['NPA'][k] = npa_charges(mol, mf)
log(f'NPA done in {time.time() - t:.0f} s; sums N {q["NPA"]["N"].sum():+.4f}, N-1 {q["NPA"]["N-1"].sum():+.4f}, '
    f'N+1 {q["NPA"]["N+1"].sum():+.4f}')

t = time.time()
free_atom_radial(mol_n)
grids = mf_n.grids if hasattr(mf_n, 'grids') else mf_n._scf.grids
if grids.coords is None:
    grids.build()
hq, nint = hirshfeld_charges(mol_n, grids, [total_dm(mf_n), total_dm(mf_m), total_dm(mf_p)])
q['Hirshfeld'] = dict(zip(['N', 'N-1', 'N+1'], hq))
log(f'Hirshfeld done in {time.time() - t:.0f} s on {grids.coords.shape[0]} grid points; integrated electrons '
    f'N {nint[0]:.5f} ({mol_n.nelectron}), N-1 {nint[1]:.5f} ({mol_m.nelectron}), N+1 {nint[2]:.5f} ({mol_p.nelectron}); '
    f'charge sums N {hq[0].sum():+.4f}, N-1 {hq[1].sum():+.4f}, N+1 {hq[2].sum():+.4f}')

fk = {}
for s in ('NPA', 'Hirshfeld'):
    qn, qm, qp = q[s]['N'], q[s]['N-1'], q[s]['N+1']
    fp, fm = qn - qp, qm - qn
    fk[s] = {'f_plus': fp, 'f_minus': fm, 'f_zero': (qm - qp) / 2, 'dual_descriptor': fp - fm}
    sums = {k: float(v.sum()) for k, v in fk[s].items()}
    log(f'{s}: sums f+ {sums["f_plus"]:.4f}, f- {sums["f_minus"]:.4f}, f0 {sums["f_zero"]:.4f}, Df {sums["dual_descriptor"]:+.4f}; '
        f'max f+ {labels[int(np.argmax(fp))]} {fp.max():.3f}, max f- {labels[int(np.argmax(fm))]} {fm.max():.3f}')

# consistency with dft/dscf.py (same code path, same geometry): energies of the three states
cmp = None
dscf_p = case['dir'] / 'dscf.json'
if dscf_p.exists() and not args.basis and not args.grid:
    d = json.loads(dscf_p.read_text())
    same_geom = d['provenance'].get('geometry_sha256') == sha256(xyz)
    cmp = {'dscf_json': str(dscf_p.relative_to(ROOT)), 'same_geometry_sha256': same_geom,
           'dE_N_Ha': float(mf_n.e_tot) - d['E_N_Ha'], 'dE_Nm1_Ha': st_m['E_Ha'] - d['N_minus_1']['E_Ha'],
           'dE_Np1_Ha': st_p['E_Ha'] - d['N_plus_1']['E_Ha'],
           'dscf_stability_checked': [d['N_minus_1'].get('stability_checked'), d['N_plus_1'].get('stability_checked')],
           'dscf_instabilities_followed': [d['N_minus_1'].get('instabilities_followed'),
                                           d['N_plus_1'].get('instabilities_followed')]}
    cmp['agree_1e-6'] = bool(max(abs(cmp['dE_N_Ha']), abs(cmp['dE_Nm1_Ha']), abs(cmp['dE_Np1_Ha'])) < 1e-6)
    log(f'vs dscf.json: dE N {cmp["dE_N_Ha"]:+.2e}, N-1 {cmp["dE_Nm1_Ha"]:+.2e}, N+1 {cmp["dE_Np1_Ha"]:+.2e} Ha '
        f'(same geometry: {same_geom})')

# consistency with dft/mwfn_extra.py: N-state NPA (same recipe) and Hirshfeld (Multiwfn built-in free atoms)
cmp_q = None
chp = case['dir'] / 'mwfn_extra' / 'charges.json'
if chp.exists() and not args.basis and not args.grid:
    old = json.loads(chp.read_text())['atomic']
    cmp_q = {'charges_json': str(chp.relative_to(ROOT))}
    for s in ('NPA', 'Hirshfeld'):
        if s in old:
            cmp_q[f'max_abs_diff_{s}'] = float(np.max(np.abs(np.array(old[s]) - q[s]['N'])))
    log('vs mwfn_extra charges.json (N state): ' + ', '.join(f'{k} {v:.2e}' for k, v in cmp_q.items() if k.startswith('max')))

HA = HARTREE_EV
I = (st_m['E_Ha'] - mf_n.e_tot) * HA
A = (mf_n.e_tot - st_p['E_Ha']) * HA
res = {
    'provenance': provenance(case, args, xyz, t_start, extra={'script': 'paper/si/compute/fukui.py'}),
    'method': ('condensed Fukui functions by finite differences of atomic charges of the N (RKS), N-1 and N+1 '
               '(UKS doublets, vertical, same geometry, equilibrium IEF-PCM for every state) systems; NPA from PySCF '
               'natural atomic orbitals of each state; Hirshfeld on the SCF grid with a promolecule of spherically '
               'averaged neutral free atoms at the same functional and basis'),
    'formulas': 'f+ = q(N) - q(N+1); f- = q(N-1) - q(N); f0 = [q(N-1) - q(N+1)]/2; dual descriptor = f+ - f-',
    'atoms': labels, 'elements': syms,
    'states': {'N': {'charge': mol_n.charge, 'E_Ha': float(mf_n.e_tot), 'nelectron': mol_n.nelectron},
               'N-1': st_m | {'nelectron': mol_m.nelectron}, 'N+1': st_p | {'nelectron': mol_p.nelectron}},
    'I_eV': I, 'A_eV': A,
    'charges': {s: {k: v.tolist() for k, v in q[s].items()} for s in q},
    'fukui': {s: {k: v.tolist() for k, v in fk[s].items()} for s in fk},
    'sums': {s: {k: float(np.sum(v)) for k, v in fk[s].items()} for s in fk},
    'hirshfeld_grid': {'points': int(grids.coords.shape[0]), 'level': grid_level(args),
                       'integrated_electrons': dict(zip(['N', 'N-1', 'N+1'], map(float, nint)))},
    'free_atoms': {el: {k: v for k, v in f.items() if k not in ('r', 'logrho')} for el, f in FREE.items()},
    'check_vs_dscf': cmp, 'check_vs_mwfn_extra_N_charges': cmp_q,
}
save_json(outdir / out_name('fukui', tag), res)

rows = []
for i, lab in enumerate(labels):
    row = {'index': i + 1, 'atom': lab}
    for s in ('NPA', 'Hirshfeld'):
        for k in ('N', 'N-1', 'N+1'):
            row[f'q_{s}_{k}'] = f'{q[s][k][i]:.5f}'
        for k in ('f_plus', 'f_minus', 'f_zero', 'dual_descriptor'):
            row[f'{k}_{s}'] = f'{fk[s][k][i]:.5f}'
    rows.append(row)
hdr = list(rows[0])
(outdir / out_name('fukui', tag, '.csv')).write_text(','.join(hdr) + '\n' + '\n'.join(','.join(str(r[h]) for h in hdr) for r in rows) + '\n')
md = (f"Condensed Fukui functions, {case['name']} ({environment_label(case)}), {basis_description(case)}\n\n"
      + md_table([(lab, *[f'{fk[s][k][i]:.3f}' for s in ('NPA', 'Hirshfeld') for k in ('f_plus', 'f_minus', 'f_zero',
                                                                                       'dual_descriptor')])
                  for i, lab in enumerate(labels)],
                 ['atom'] + [f'{k} ({s})' for s in ('NPA', 'Hirshfeld') for k in ('f+', 'f-', 'f0', 'Df')]))
(outdir / out_name('fukui', tag, '.md')).write_text(md)
log(f'done in {(time.time() - t_start) / 60:.1f} min; I = {I:.4f} eV, A = {A:.4f} eV')
