"""
Functional check: single points with several methods on a list of geometries of one case, and the relative
energies between them.  Used on the proton-scan end points and the top of the profile (referee point 2.4:
B3LYP tends to underestimate proton-transfer barriers).

Usage: python dft/fsp.py <case> GEOM1.xyz [GEOM2.xyz ...] [--names n1,n2,...] [--label NAME]
                          [--methods b3lyp-d3bj,wb97x-d,m06-2x,pbe0-d3bj,df-mp2] [common switches]
       python dft/fsp.py <case> --scan results/<case>/scan [...]   # end points + top of the lower envelope
Writes results/<case>/fsp_<label>.json and .md.  Charge, spin, basis and environment are the case's
(6-311++G(d,p); IEF-PCM water for PCM cases), so the geometries must have the case's atoms in its order.

Methods (all with density fitting, the case grid, and the case solvent):
  b3lyp-d3bj  the reference level of the project, for a like-for-like comparison.
  wb97x-d     Chai & Head-Gordon 2008 (PCCP 10, 6615).  PySCF 2.14 refuses xc = 'wb97x-d' (it is on the
              black list of pyscf/scf/dispersion.py because its dispersion term is not implemented), so the
              exchange-correlation part is libxc's HYB_GGA_XC_WB97X_D and the -D term is added here: the
              Chai-Head-Gordon damped C6/R^6 sum, E = -sum_ij C6ij/R^6 / (1 + 6 (R/Rr)^-12), C6ij = sqrt(C6i C6j),
              Rr = R0i + R0j, with Grimme's D2 C6 and R0 (J. Comput. Chem. 27, 2006, 1787), s6 = 1.  This is how
              Psi4 defines "-CHG" (psi4/src/psi4/libdisp/dispersion.cc, d = 6.0, C6_D2, RvdW_D2, geometric mean),
              checked against that source on 27 Sep 2026.  The D-term is geometry-only and is reported separately.
  m06-2x      Zhao & Truhlar 2008; libxc 'm062x', no dispersion correction (as defined).  Meta-GGA: same grid as
              the rest; relative energies between nearby geometries are insensitive to it.
  pbe0-d3bj   PBE0 with Grimme D3(BJ) (simple-dftd3 parameters for pbe0).
  df-mp2      RI-MP2 on a DF-RHF reference, frozen core (1s for C/N/O, 1s2s2p for P/S, as Gaussian's FC), with
              the MP2-fitting auxiliary basis PySCF assigns (aug-cc-pVTZ-RI for 6-311++G(d,p)); the JK fit of the
              RHF is aug-cc-pVTZ-JKFIT.  In PCM the MP2 correlation is computed with the PCM-polarised RHF orbitals
              and the reaction field frozen at the RHF density (the "PTE" scheme of Cammi, Mennucci & Tomasi):
              PySCF's own solvent-MP2 wrapper would iterate the reaction field with relaxed MP2 densities, which
              its source notes "often causes convergence issue", and costs one MP2 per iteration.
Relative energies are given against the first geometry and against the lowest one, per method, in kcal/mol and
kJ/mol; electronic energies only (ZPE is dft/zpe_barrier.py's job).
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from props_common import *          # noqa: F401,F403 -- pins threads / scratch, so it must precede numpy and PySCF
import glob
from pyscf import scf as _scf, df as _df
from pyscf.mp import dfmp2
from pyscf.data.elements import chemcore

ALL_METHODS = ['b3lyp-d3bj', 'wb97x-d', 'm06-2x', 'pbe0-d3bj', 'df-mp2']
LABELS = {'b3lyp-d3bj': 'B3LYP-D3(BJ)', 'wb97x-d': 'wB97X-D', 'm06-2x': 'M06-2X', 'pbe0-d3bj': 'PBE0-D3(BJ)',
          'df-mp2': 'DF-MP2 (FC)'}
DFT = {'b3lyp-d3bj': ('b3lyp', 'd3bj'), 'wb97x-d': ('HYB_GGA_XC_WB97X_D', 'none'), 'm06-2x': ('m062x', 'none'),
       'pbe0-d3bj': ('pbe0', 'd3bj')}

# Grimme D2 parameters (J. Comput. Chem. 27, 2006, 1787, Table 1): C6 in J nm^6 mol^-1, R0 in Angstrom.
# Cross-checked against Psi4's C6_D2_ / RvdW_D2_ tables (dispersion_defines.h) for these elements.
D2 = {'H': (0.14, 1.001), 'C': (1.75, 1.452), 'N': (1.23, 1.397), 'O': (0.70, 1.342), 'P': (7.84, 1.705),
      'S': (5.57, 1.683)}
JNM6_TO_AU = 17.34527758          # 1 J nm^6 mol^-1 in Eh a0^6


def chg_dispersion(mol):
    """Chai-Head-Gordon -D energy of wB97X-D (Hartree) for the real atoms of mol."""
    sy = [mol.atom_pure_symbol(i) for i in range(mol.natm)]
    x = mol.atom_coords()
    e = 0.0
    for i in range(mol.natm):
        for j in range(i):
            c6 = np.sqrt(D2[sy[i]][0] * D2[sy[j]][0]) * JNM6_TO_AU
            rr = (D2[sy[i]][1] + D2[sy[j]][1]) / BOHR_A
            r = np.linalg.norm(x[i] - x[j])
            e -= c6 / r ** 6 / (1 + 6.0 * (r / rr) ** -12)
    return e


def from_scan(scan_dir):
    """End points and top(s) of a dft/scan.py lower envelope -> list of (name, xyz path)."""
    sd = pathlib.Path(scan_dir)
    summ = json.load(open(sd / 'scan_summary.json'))
    env = summ['lower_envelope']
    picks = [('end_N', env[0])] + [(f"top{k + 1}", min(env, key=lambda r: abs(r['d_N1_H1A_A'] - m['d_N1_H1A_A'])))
                                  for k, m in enumerate(summ.get('envelope_maxima', []))] + [('end_O', env[-1])]
    out = []
    for nm, r in picks:
        hits = [p for p in glob.glob(str(sd / f"{r['chain']}_{r['index']:02d}_*.xyz")) if not p.endswith('.partial.xyz')]
        if len(hits) != 1:
            raise SystemExit(f"[fsp] scan point {r['chain']}_{r['index']:02d}: expected one .xyz, found {hits}")
        out.append((f"{nm}_d{r['d_N1_H1A_A']:.3f}", pathlib.Path(hits[0])))
    return out


ap = parser(__doc__, geometry=False)
ap.add_argument('geometries', nargs='*', help='xyz files (same atoms and order as the case)')
ap.add_argument('--scan', help='a dft/scan.py output directory: use its end points and top(s)')
ap.add_argument('--names', help='comma-separated names for the geometries (default: file stems)')
ap.add_argument('--label', default='fsp', help='output name stem (fsp_<label>.json)')
ap.add_argument('--methods', default=','.join(ALL_METHODS))
args = ap.parse_args()
t_start = time.time()
case, _, outdir, tag = setup(args, need_geometry=False)
methods = [m.strip() for m in args.methods.split(',')]
for m in methods:
    if m not in ALL_METHODS:
        raise SystemExit(f'[fsp] unknown method {m}; choose from {ALL_METHODS}')
geoms = from_scan(args.scan) if args.scan else []
names = args.names.split(',') if args.names else [pathlib.Path(g).stem for g in args.geometries]
geoms += [(n, pathlib.Path(g).resolve()) for n, g in zip(names, args.geometries)]
if not geoms:
    raise SystemExit('[fsp] no geometries given')
syms0, _ = read_xyz(ROOT / case['xyz'])
for n, g in geoms:
    s, _ = read_xyz(g)
    if [a.capitalize() for a in s] != [a.capitalize() for a in syms0]:
        raise SystemExit(f'[fsp] {g}: atoms differ from the case input {case["xyz"]}')
solv = case.get('solvent')


def rhf_mp2(mol, dm0=None):
    mf = _scf.RHF(mol).density_fit()
    mf.conv_tol, mf.max_cycle = 1e-9, 150
    if solv:
        mf = mf.PCM()
        mf.with_solvent.method = 'IEF-PCM'
        mf.with_solvent.eps = SOLVENT_EPS[solv.lower()]
        mf.with_solvent.lebedev_order = int(os.environ.get('PCM_LEBEDEV_ORDER', '29'))
    converge(mf, dm0=dm0, what='RHF')
    # PTE: a plain DF-MP2 object on the PCM-polarised orbitals (not mf.MP2(), which would wrap it in PySCF's
    # iterative solvent-MP2 scheme); the RHF total energy already contains the HF-level solvation energy.
    pt = dfmp2.DFMP2(mf, frozen=chemcore(mol))
    pt.with_df = _df.DF(mol, auxbasis=_df.make_auxbasis(mol, mp2fit=True))
    # energy only: holding the t2 amplitudes in core needs more than DFT_MAX_MEMORY_MB for the 27-atom pair
    # (burst 1b failed here with "Insufficient memory for holding t2 incore", 28 Sep 2026)
    pt.kernel(with_t2=False)
    return mf, {'E_Ha': float(mf.e_tot + pt.e_corr), 'E_RHF_Ha': float(mf.e_tot), 'E_corr_Ha': float(pt.e_corr),
                'n_frozen': int(chemcore(mol)), 'aux_mp2': str(pt.with_df.auxbasis)}


E = {m: {} for m in methods}
for m in methods:
    dm_prev = None
    for n, g in geoms:
        t = time.time()
        mol = mol_for_case(case, g, verbose=args.verbose)
        if m == 'df-mp2':
            mf, rec = rhf_mp2(mol, dm_prev)
        else:
            xc, disp = DFT[m]
            mf = scf_for(case, mol, args, xc=xc, disp=disp, conv_tol=1e-9)
            converge(mf, dm0=dm_prev, what=f'{m} SCF')
            rec = {'E_Ha': float(mf.e_tot)}
            if m == 'wb97x-d':
                ed = chg_dispersion(mol)
                rec.update({'E_scf_Ha': float(mf.e_tot), 'E_disp_CHG_Ha': float(ed), 'E_Ha': float(mf.e_tot + ed)})
            elif disp != 'none':
                rec['E_disp_Ha'] = float(mf.scf_summary.get('dispersion', 0.0))
        dm_prev = mf.make_rdm1()        # guess for the next (nearby) geometry
        rec['wall_s'] = round(time.time() - t, 1)
        E[m][n] = rec
        print(f"[fsp] {LABELS[m]:14s} {n:24s} E = {rec['E_Ha']:.10f} Ha  ({rec['wall_s']:.0f} s)", flush=True)
        # keep what is finished if a later method fails (burst 1b lost four DFT methods to an MP2 crash)
        save_json(outdir / out_name(args.label + '_partial', tag), {'energies': E, 'geometries': [n for n, _ in geoms]})

rel = {}
for m in methods:
    e = np.array([E[m][n]['E_Ha'] for n, _ in geoms])
    rel[m] = {'vs_first_kcal': dict(zip([n for n, _ in geoms], ((e - e[0]) * KCAL_PER_HA).tolist())),
              'vs_first_kJ': dict(zip([n for n, _ in geoms], ((e - e[0]) * KJ_PER_HA).tolist())),
              'vs_lowest_kcal': dict(zip([n for n, _ in geoms], ((e - e.min()) * KCAL_PER_HA).tolist()))}
res = {'provenance': provenance(case, args, None, t_start),
       'geometries': [{'name': n, 'file': str(g.relative_to(ROOT)) if str(g).startswith(str(ROOT)) else str(g),
                       'sha256': sha256(g)} for n, g in geoms],
       'methods': {m: LABELS[m] for m in methods}, 'energies': E, 'relative': rel,
       'notes': ('Electronic energies; wB97X-D = libxc HYB_GGA_XC_WB97X_D + Chai-Head-Gordon D2-type dispersion '
                 '(added here); DF-MP2 frozen core' + (', PCM via PTE (reaction field frozen at RHF)' if solv else '')
                 + '. Relative energies per method.')}
name = out_name(args.label, tag)
save_json(outdir / name, res)
hdr = ['Method'] + [n for n, _ in geoms]
rows = [[LABELS[m]] + [f"{rel[m]['vs_first_kcal'][n]:.2f}" for n, _ in geoms] for m in methods]
(outdir / out_name(args.label, tag, '.md')).write_text(
    f"Relative energies (kcal/mol) vs {geoms[0][0]}, {case['name']} ({environment_label(case)}), "
    f"{basis_description(case)}\n\n" + md_table(rows, hdr))
if outdir == case['dir']:
    stage_note(case, f"functional check ({name}): " + '; '.join(
        f"{LABELS[m]} " + ', '.join(f"{rel[m]['vs_first_kcal'][n]:+.2f}" for n, _ in geoms[1:]) for m in methods)
        + ' kcal/mol vs ' + geoms[0][0])
print(open(outdir / out_name(args.label, tag, '.md')).read())
print(f"[fsp] {len(geoms)} geometries x {len(methods)} methods, {(time.time() - t_start) / 60:.1f} min -> {outdir / name}")
