"""
Vertical ionisation potential and electron affinity by Delta-SCF, and the conceptual-DFT descriptors built from
them, with the Koopmans (frontier-orbital) values alongside for comparison, clearly labelled.

Usage: python dft/dscf.py <case> [--geom-from CASE | --xyz FILE] [--no-stability] [common switches]
Writes results/<case>/dscf[_<tag>].json and dscf[_<tag>].md.  Intended for M2 (cation_*), M4 (sulfanilamide_*)
and M1 (sfa_pair_*), PCM and gas.

Why Delta-SCF: the M1 frontier orbitals sit on different ions (HOMO on H2PO4-, LUMO on the cation), so the
Koopmans gap of the ion pair is a charge-transfer gap and hardness / electrophilicity derived from it have no
clear meaning (referee critique 2.7).  Energy differences of the N-1, N and N+1 electron systems do not depend
on that orbital picture.

Method
  * N electrons: the case's own model (B3LYP-D3(BJ), DF, IEF-PCM water or gas) at the optimised geometry.
  * N-1 and N+1: UKS doublets at the SAME geometry (vertical), same functional, basis, grid and environment.
    Initial guess from the neutral orbitals (beta HOMO emptied / alpha LUMO filled), so the SCF lands on the
    Koopmans-like state; an internal UKS stability analysis follows (on by default) and, if an instability is
    found, the SCF is restarted from the lower solution (up to 3 times).  <S^2> is reported for both ions.
  * D3(BJ) depends on the geometry only, so it cancels exactly in I and A; it is kept in the total energies.
  * PCM: each charge state is solvated at equilibrium (the reaction field relaxes to the ion's density with the
    nuclei fixed), i.e. Gaussian's default SCRF on each state.  A strictly vertical process in solution would
    use non-equilibrium solvation; the equilibrium value is the standard choice for conceptual-DFT descriptors
    and is stated as such.  In gas phase an anion's extra electron may be unbound (A < 0 with diffuse
    functions); this is reported, not hidden.
  Descriptors (both from Delta-SCF and from Koopmans I = -e_HOMO, A = -e_LUMO of the N-electron system):
    chi = (I + A)/2,  mu = -chi,  eta = (I - A)/2,  S = 1/(2 eta),  omega = mu^2/(2 eta),  DeltaN_max = -mu/eta.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from props_common import *          # noqa: F401,F403 -- pins threads / scratch, so it must precede numpy and PySCF

ap = parser(__doc__)
ap.add_argument('--no-stability', action='store_true', help='skip the UKS internal stability analysis of N+-1')
args = ap.parse_args()
t_start = time.time()
case, xyz, outdir, tag = setup(args)
if case['spin'] != 0:
    raise SystemExit('[dscf] closed-shell N-electron cases only')


def descriptors(I, A):
    chi = (I + A) / 2; mu = -chi; eta = (I - A) / 2
    return {'I_eV': I, 'A_eV': A, 'electronegativity_chi_eV': chi, 'chemical_potential_mu_eV': mu,
            'hardness_eta_eV': eta, 'softness_S_eV-1': 1 / (2 * eta), 'electrophilicity_omega_eV': mu ** 2 / (2 * eta),
            'DeltaN_max': -mu / eta}


def run_state(dq, mf_n):
    """Vertical N-dq electron state (dq = +1: cation, N-1; dq = -1: anion, N+1) as a UKS doublet."""
    mol = mol_for_case(dict(case, charge=case['charge'] + dq, spin=1), xyz, verbose=args.verbose)
    mf = scf_for(case, mol, args, conv_tol=1e-9)
    mf.max_cycle = 200
    C, occ = mf_n.mo_coeff, mf_n.mo_occ
    n = int((occ > 0).sum())
    oa, ob = np.zeros(C.shape[1]), np.zeros(C.shape[1])
    if dq == +1:
        oa[:n] = 1; ob[:n - 1] = 1          # remove a beta electron from the HOMO
    else:
        oa[:n + 1] = 1; ob[:n] = 1          # add an alpha electron to the LUMO
    dm0 = np.array([(C * oa) @ C.T, (C * ob) @ C.T])
    what = 'N-1 (cation)' if dq == 1 else 'N+1 (anion)'
    t = time.time()
    converge(mf, dm0=dm0, what=what)
    n_unstable = 0
    if not args.no_stability:
        for _ in range(3):
            mo_i, _, stable, _ = mf.stability(return_status=True)
            if stable:
                break
            n_unstable += 1
            print(f'[dscf] {what}: internal instability found, restarting from the lower solution', flush=True)
            converge(mf, dm0=mf.make_rdm1(mo_i, mf.mo_occ), what=what)
    s2, mult = mf.spin_square()
    print(f'[dscf] {what}: E = {mf.e_tot:.10f} Ha, <S^2> = {s2:.4f}, {time.time() - t:.0f} s', flush=True)
    return {'charge': mol.charge, 'multiplicity': 2, 'E_Ha': float(mf.e_tot), 'S2': float(s2),
            'S2_ideal': 0.75, 'stability_checked': not args.no_stability, 'instabilities_followed': n_unstable,
            'converged': bool(mf.converged), 'wall_s': round(time.time() - t, 1)}


mol = mol_for_case(case, xyz, verbose=args.verbose)
mf = scf_for(case, mol, args, conv_tol=1e-9)
t = time.time()
converge(mf, what='N-electron SCF')
e_n = float(mf.e_tot)
nocc = int((mf.mo_occ > 0).sum())
homo, lumo = float(mf.mo_energy[nocc - 1] * HARTREE_EV), float(mf.mo_energy[nocc] * HARTREE_EV)
print(f'[dscf] N: E = {e_n:.10f} Ha, HOMO {homo:.4f} eV, LUMO {lumo:.4f} eV ({time.time() - t:.0f} s)', flush=True)
disp = float(mf.scf_summary.get('dispersion', 0.0)) if hasattr(mf, 'scf_summary') else 0.0

cat = run_state(+1, mf)
ani = run_state(-1, mf)
I = (cat['E_Ha'] - e_n) * HARTREE_EV
A = (e_n - ani['E_Ha']) * HARTREE_EV
d_scf = descriptors(I, A)
d_koop = descriptors(-homo, -lumo)
res = {
    'provenance': provenance(case, args, xyz, t_start),
    'method': ('vertical Delta-SCF: RKS N / UKS doublets N-1, N+1 at the N-electron geometry, same functional, '
               'basis, grid and environment' + ('; PCM at equilibrium for every charge state' if case.get('solvent') else '')),
    'E_N_Ha': e_n, 'dispersion_Ha': disp, 'HOMO_eV': homo, 'LUMO_eV': lumo, 'gap_KS_eV': lumo - homo,
    'N_minus_1': cat, 'N_plus_1': ani,
    'delta_scf': d_scf,
    'koopmans': d_koop,
    'koopmans_label': 'Koopmans-type estimate from the N-electron Kohn-Sham orbitals: I = -e(HOMO), A = -e(LUMO)',
    'fundamental_gap_I_minus_A_eV': I - A,
    'flags': {'electron_unbound_in_gas_phase': bool(A < 0 and not case.get('solvent')),
              'spin_contamination_over_0.05': [k for k, s in (('N-1', cat['S2']), ('N+1', ani['S2'])) if s - 0.75 > 0.05]},
    'formulas': 'chi=(I+A)/2, mu=-chi, eta=(I-A)/2, S=1/(2 eta), omega=mu^2/(2 eta), DeltaN_max=-mu/eta',
}
name = out_name('dscf', tag)
save_json(outdir / name, res)
labels = [('I', 'I_eV', 'eV'), ('A', 'A_eV', 'eV'), ('chi = (I+A)/2', 'electronegativity_chi_eV', 'eV'),
          ('mu = -chi', 'chemical_potential_mu_eV', 'eV'), ('eta = (I-A)/2', 'hardness_eta_eV', 'eV'),
          ('S = 1/(2 eta)', 'softness_S_eV-1', 'eV^-1'), ('omega = mu^2/(2 eta)', 'electrophilicity_omega_eV', 'eV'),
          ('DeltaN_max = -mu/eta', 'DeltaN_max', '')]
(outdir / out_name('dscf', tag, '.md')).write_text(
    f"Vertical Delta-SCF descriptors, {case['name']} ({environment_label(case)})\n\n"
    + md_table([(l, f'{d_scf[k]:.4f}', f'{d_koop[k]:.4f}', u) for l, k, u in labels],
               ['Descriptor', 'Delta-SCF', 'Koopmans (orbital)', 'unit']))
if outdir == case['dir']:
    stage_note(case, f"Delta-SCF ({name}): I = {I:.3f} eV, A = {A:.3f} eV, eta = {d_scf['hardness_eta_eV']:.3f} eV "
                     f"(Koopmans I = {-homo:.3f}, A = {-lumo:.3f})")
print(f"[dscf] I = {I:.4f} eV, A = {A:.4f} eV (Koopmans {-homo:.4f} / {-lumo:.4f}); eta = {d_scf['hardness_eta_eV']:.4f} eV, "
      f"omega = {d_scf['electrophilicity_omega_eV']:.4f} eV; <S^2> {cat['S2']:.3f} / {ani['S2']:.3f}; "
      f"{(time.time() - t_start) / 60:.1f} min -> {outdir / name}", flush=True)
