"""
Vertical singlet excitations (UV-Vis) by linear-response TD-DFT: CAM-B3LYP primary, B3LYP for comparison, with
oscillator strengths, dominant orbital transitions, natural transition orbitals, fragment charge-transfer
analysis and a Gaussian-broadened spectrum.

Usage: python dft/tddft.py <case> [--xc camb3lyp[,b3lyp]] [--solvation noneq[,eq]] [--nstates 20] [--tda]
                            [--fwhm 0.3] [common switches]
Writes, per functional and solvation mode, results/<case>/tddft_<xc>_<mode>[_tda][_<tag>].json, a spectrum CSV
(spectrum_<xc>_<mode>...csv: energy, wavelength, molar absorptivity) and a short .md table; the leading NTO pairs
go to nto_<...>.npz for rendering.  <mode> is gas, noneq or eq.

Why CAM-B3LYP: the gas-phase B3LYP TDA benchmark put the first excitation of the X-ray ion pair at 2.44 eV
(508 nm) for a crystal the CIF records as colourless -- the classic charge-transfer underestimate of a global
hybrid (HOMO on H2PO4-, LUMO on the cation).  A range-separated hybrid (CAM-B3LYP: 19% -> 65% exact exchange,
mu = 0.33) removes most of it; B3LYP is kept to show the artefact.  The fragment analysis below measures it.

Method
  * Ground state: the case's geometry (B3LYP-D3(BJ) optimum), the chosen functional, DF, the case basis and grid.
    D3 is omitted in the TD step: it is a geometry-only energy term and does not touch the density, so it has
    no effect on excitation energies.
  * Full linear-response TD-DFT (RPA, "TDDFT"), 20 singlets by default; --tda switches to the Tamm-Dancoff
    approximation.  PCM TD-DFT is supported by PySCF 2.14 for both RPA and TDA (solvent/_attach_solvent.py,
    TDSCFWithSolvent; regression-tested against Q-Chem in pyscf/solvent/test/test_pcm_tdscf.py at tag v2.14.0,
    non-equilibrium only, which the Pi code test reproduces):
      noneq  non-equilibrium (default): the fast solvent response uses the optical dielectric constant,
             PySCF's fixed eps_inf = 1.78 (water: n^2 = 1.776), the slow part stays with the ground state --
             the appropriate model for a vertical absorption;
      eq     equilibrium: the full static eps responds to the transition density (comparison only; it
             over-stabilises charge-transfer states).
    Linear-response PCM only; no state-specific correction.
  * Per state: E (eV, nm), oscillator strength (length gauge), transition dipole, the dominant occupied ->
    virtual pairs with weights (X^2 - Y^2 normalised to 1, labelled H-n -> L+m), NTOs from the SVD of X+Y
    (first-pair weight), and, when the case defines fragments (M1: cation / H2PO4-), a charge-transfer analysis
    of the X+Y transition density following Plasser & Lischka (JCTC 8, 2012, 2777): Omega_AB with Mulliken
    partitioning, hole and electron populations per fragment, the CT number (share of Omega off the fragment
    diagonal) and the hole-electron centroid distance D_CT from atom-resolved Omega.
  * Spectrum: sum of Gaussians of FWHM 0.3 eV (in energy) weighted by f, converted to molar absorptivity with
    eps = f * g(E) * 2.8706e4 L mol-1 cm-1 per eV^-1 (from f = 4.319e-9 int eps dnu, 1 eV = 8065.544 cm-1).
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from props_common import *          # noqa: F401,F403 -- pins threads / scratch, so it must precede numpy and PySCF

HC_EV_NM = 1239.841984            # h c in eV nm
EPS_PER_F = 1.0 / (4.319e-9 * 8065.544)     # L mol-1 cm-1 per (unit f x unit eV^-1 lineshape)
XC_LABEL = {'camb3lyp': 'CAM-B3LYP', 'b3lyp': 'B3LYP', 'wb97x': 'wB97X', 'pbe0': 'PBE0'}

ap = parser(__doc__)
ap.add_argument('--xc', default='camb3lyp', help='comma-separated functionals (PySCF names)')
ap.add_argument('--solvation', default='noneq', help='PCM cases: comma-separated noneq, eq')
ap.add_argument('--nstates', type=int, default=20)
ap.add_argument('--tda', action='store_true', help='Tamm-Dancoff approximation instead of full TD-DFT')
ap.add_argument('--fwhm', type=float, default=0.3, help='Gaussian FWHM of the broadened spectrum, eV')
ap.add_argument('--conv', type=float, default=1e-5, help='Davidson residual threshold')
args = ap.parse_args()
t_start = time.time()
case, xyz, outdir, tag = setup(args)
if case['spin'] != 0:
    raise SystemExit('[tddft] closed-shell cases only (singlet excitations of an RKS ground state)')
modes = [m.strip() for m in args.solvation.split(',')] if case.get('solvent') else ['gas']
for m in modes:
    if m not in ('gas', 'noneq', 'eq'):
        raise SystemExit(f'[tddft] unknown solvation mode {m!r}')

mol = mol_for_case(case, xyz, verbose=args.verbose)
S = mol.intor_symmetric('int1e_ovlp')
aoslice = mol.aoslice_by_atom()
ao_atom = np.zeros(mol.nao, int)
for a, (_, _, p0, p1) in enumerate(aoslice):
    ao_atom[p0:p1] = a
R_A = mol.atom_coords() * BOHR_A
frags = case.get('fragments') or {}


def mulliken_atoms(c):
    """Mulliken population of a normalised orbital (AO coefficients c) per atom."""
    return np.bincount(ao_atom, weights=c * (S @ c), minlength=mol.natm)


def frag_sum(pop_atoms):
    return {k: float(pop_atoms[v].sum()) for k, v in frags.items()}


def omega_analysis(T):
    """Plasser-Lischka Omega matrix (atoms x atoms) of an AO transition density T (hole index first)."""
    TS, ST = T @ S, S @ T
    om = 0.5 * (TS * ST + T * (S @ T @ S))
    Om = np.zeros((mol.natm, mol.natm))
    np.add.at(Om, (ao_atom[:, None].repeat(mol.nao, 1), ao_atom[None, :].repeat(mol.nao, 0)), om)
    tot = Om.sum()
    Om /= tot
    hole, elec = Om.sum(1), Om.sum(0)
    rh, re = hole @ R_A, elec @ R_A
    out = {'omega_total_raw': float(tot), 'D_CT_A': float(np.linalg.norm(re - rh)),
           'hole_centroid_A': rh, 'electron_centroid_A': re}
    if frags:
        names = list(frags)
        F = np.array([[Om[np.ix_(frags[a], frags[b])].sum() for b in names] for a in names])
        out.update({'fragments': names, 'omega_frag': F,
                    'hole_on': {n: float(F[i].sum()) for i, n in enumerate(names)},
                    'electron_on': {n: float(F[:, i].sum()) for i, n in enumerate(names)},
                    'CT_number': float(F.sum() - np.trace(F))})
    return out


def label(i, a, nocc):
    h = 'H' if i == nocc - 1 else f'H-{nocc - 1 - i}'
    l = 'L' if a == 0 else f'L+{a}'
    return f'{h}->{l}'


results = []
for xc in [x.strip() for x in args.xc.split(',')]:
    t = time.time()
    mf = scf_for(case, mol, args, xc=xc, disp='none', conv_tol=1e-9)
    converge(mf, what=f'{xc} ground state')
    nocc = int((mf.mo_occ > 0).sum())
    Cocc, Cvir = mf.mo_coeff[:, :nocc], mf.mo_coeff[:, nocc:]
    e_mo = mf.mo_energy * HARTREE_EV
    print(f'[tddft] {xc}: ground state E = {mf.e_tot:.8f} Ha, HOMO {e_mo[nocc - 1]:.3f} eV, LUMO {e_mo[nocc]:.3f} eV '
          f'({time.time() - t:.0f} s)', flush=True)
    for mode in modes:
        t = time.time()
        if mode == 'gas':
            td = mf.TDA() if args.tda else mf.TDDFT()
        else:
            eq = (mode == 'eq')
            td = mf.TDA(equilibrium_solvation=eq) if args.tda else mf.TDDFT(equilibrium_solvation=eq)
        td.nstates = args.nstates
        td.conv_tol = args.conv
        td.max_cycle = 300
        td.kernel()
        conv = np.atleast_1d(td.converged)
        E = np.asarray(td.e) * HARTREE_EV
        f = np.asarray(td.oscillator_strength(gauge='length'))
        tdm = np.asarray(td.transition_dipole())
        states, ntos = [], {}
        for k in range(len(E)):
            X, Y = td.xy[k]
            X = np.asarray(X).reshape(nocc, -1)
            Y = np.zeros_like(X) if (args.tda or np.ndim(Y) == 0) else np.asarray(Y).reshape(nocc, -1)
            w = X ** 2 - Y ** 2
            w /= w.sum()
            order = np.argsort(-w.ravel())
            dom = []
            for idx in order[:6]:
                i, a = divmod(int(idx), X.shape[1])
                if w.ravel()[idx] < 0.02:
                    break
                dom.append({'from_mo': i + 1, 'to_mo': nocc + a + 1, 'label': label(i, a, nocc),
                            'weight': float(w.ravel()[idx])})
            Tmo = X + Y
            U, sv, Vt = np.linalg.svd(Tmo, full_matrices=False)
            nto_w = sv ** 2 / (sv ** 2).sum()
            hole1, part1 = Cocc @ U[:, 0], Cvir @ Vt[0]
            if k < 10:
                ntos[f's{k + 1}_hole'] = Cocc @ U[:, :2]; ntos[f's{k + 1}_particle'] = Cvir @ Vt[:2].T
                ntos[f's{k + 1}_weights'] = nto_w[:2]
            st = {'state': k + 1, 'E_eV': float(E[k]), 'lambda_nm': float(HC_EV_NM / E[k]), 'f': float(f[k]),
                  'transition_dipole_au': tdm[k], 'dominant': dom,
                  'nto_first_pair_weight': float(nto_w[0]), 'nto_pairs_for_90pct': int(np.searchsorted(np.cumsum(nto_w), 0.9) + 1),
                  'converged': bool(conv[k]) if k < len(conv) else None}
            T_ao = Cocc @ Tmo @ Cvir.T
            st['ct'] = omega_analysis(T_ao)
            if frags:
                st['nto1_hole_on'] = frag_sum(mulliken_atoms(hole1))
                st['nto1_particle_on'] = frag_sum(mulliken_atoms(part1))
            states.append(st)
        # broadened spectrum
        sig = args.fwhm / (2 * np.sqrt(2 * np.log(2)))
        grid = np.arange(max(0.5, E.min() - 1.5), E.max() + 1.5, 0.005)
        g = np.exp(-0.5 * ((grid[:, None] - E[None, :]) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
        spec_f = g @ f                                  # oscillator strength per eV
        eps = spec_f * EPS_PER_F
        mtag = f"{xc}_{mode}{'_tda' if args.tda else ''}"
        csv = outdir / out_name(f'spectrum_{mtag}', tag, '.csv')
        with open(csv, 'w') as fh:
            fh.write('energy_eV,wavelength_nm,f_per_eV,epsilon_L_mol-1_cm-1\n')
            for x, y, z in zip(grid, spec_f, eps):
                fh.write(f'{x:.4f},{HC_EV_NM / x:.3f},{y:.6e},{z:.3f}\n')
        bright = [s for s in states if s['f'] >= 0.01]
        ipk = int(np.argmax(eps))
        env = ('gas phase' if mode == 'gas' else environment_label(case)
               + ('; non-equilibrium, optical eps_inf = 1.78 (PySCF)' if mode == 'noneq' else '; equilibrium (static eps)'))
        res = {
            'provenance': provenance(case, args, xyz, t_start, {'td_wall_min': round((time.time() - t) / 60, 2)}),
            'functional': XC_LABEL.get(xc, xc), 'xc_pyscf': xc, 'td_method': 'TDA' if args.tda else 'TDDFT (RPA)',
            'solvation': mode, 'environment': env, 'nstates': len(E), 'all_converged': bool(conv.all()),
            'ground_state_E_Ha': float(mf.e_tot), 'HOMO_eV': float(e_mo[nocc - 1]), 'LUMO_eV': float(e_mo[nocc]),
            'n_occ': nocc, 'fragments': list(frags) or None,
            'S1_eV': float(E[0]), 'S1_nm': float(HC_EV_NM / E[0]), 'S1_f': float(f[0]),
            'first_bright_state': bright[0]['state'] if bright else None,
            'first_bright_eV': bright[0]['E_eV'] if bright else None,
            'first_bright_nm': bright[0]['lambda_nm'] if bright else None,
            'spectrum_peak_nm': float(HC_EV_NM / grid[ipk]), 'spectrum_peak_eps': float(eps[ipk]),
            'spectrum_csv': csv.name, 'fwhm_eV': args.fwhm,
            'states': states,
            'notes': ('D3 omitted in the TD step (no effect on excitation energies). Weights = (X^2 - Y^2) normalised; '
                      'CT analysis on X+Y with Mulliken partitioning (Plasser & Lischka 2012).'),
        }
        name = out_name(f'tddft_{mtag}', tag)
        save_json(outdir / name, res)
        np.savez_compressed(outdir / out_name(f'nto_{mtag}', tag, '.npz'), **ntos)
        rows = []
        for s in states[:min(len(states), 12)]:
            ct = s['ct']
            rows.append((s['state'], f"{s['E_eV']:.3f}", f"{s['lambda_nm']:.1f}", f"{s['f']:.4f}",
                         ', '.join(f"{d['label']} ({100 * d['weight']:.0f}%)" for d in s['dominant'][:3]),
                         f"{s['nto_first_pair_weight']:.2f}",
                         f"{ct['CT_number']:.2f}" if 'CT_number' in ct else '-', f"{ct['D_CT_A']:.2f}"))
        (outdir / out_name(f'tddft_{mtag}', tag, '.md')).write_text(
            f"{XC_LABEL.get(xc, xc)} {'TDA' if args.tda else 'TD-DFT'}, {case['name']}, {env}\n\n"
            + md_table(rows, ['S', 'E (eV)', 'lambda (nm)', 'f', 'dominant', 'NTO1 w', 'CT number', 'D_CT (A)']))
        if outdir == case['dir']:
            stage_note(case, f"TD-DFT ({name}): S1 {E[0]:.3f} eV ({HC_EV_NM / E[0]:.0f} nm, f {f[0]:.3f}); first bright "
                             f"{res['first_bright_nm'] and round(res['first_bright_nm'])} nm; "
                             f"{'all converged' if conv.all() else 'NOT all roots converged'}")
        s1 = states[0]
        print(f"[tddft] {xc}/{mode}: S1 = {E[0]:.3f} eV ({HC_EV_NM / E[0]:.1f} nm, f = {f[0]:.4f}, "
              f"{', '.join(d['label'] for d in s1['dominant'][:2])}"
              f"{', CT number %.2f' % s1['ct']['CT_number'] if 'CT_number' in s1['ct'] else ''}); first bright "
              f"{res['first_bright_eV'] and round(res['first_bright_eV'], 3)} eV; converged {int(conv.sum())}/{len(conv)}; "
              f"{(time.time() - t) / 60:.1f} min -> {outdir / name}", flush=True)
        results.append(name)
print(f"[tddft] done: {', '.join(results)} ({(time.time() - t_start) / 60:.1f} min)")
