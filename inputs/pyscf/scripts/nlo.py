"""
Static electric properties by finite field: dipole moment mu, polarizability alpha (tensor, isotropic value,
anisotropy) and first hyperpolarizability beta (full tensor, beta_tot, beta along mu), in a.u. and esu.

Usage: python dft/nlo.py <case> [--geom-from CASE | --xyz FILE] [--field 0.001] [--no-analytic] [common switches]
Writes results/<case>/nlo[_<tag>].json and nlo[_<tag>].md.  Run urea_gas_d3 at the same level as the reference.

Why finite field: PySCF 2.14 has no analytic hyperpolarizability in its core (pyscf.prop, the separate
pyscf-properties package, is not installed on either machine), and a finite-field scheme works unchanged for
gas and PCM, with density fitting and D3.  It costs 37 SCFs, each started from the zero-field density.

Method
  * The field enters as h = h_core + F.r for the electrons and -sum_A Z_A F.R_A for the nuclei, with r measured
    from the centre of nuclear charge.  Then E(F) = E0 - mu_i F_i - 1/2 alpha_ij F_i F_j - 1/6 beta_ijk F_i F_j F_k
    - ... (the "Taylor" convention, as Gaussian's Polar keyword), and mu_i(F) = -dE/dF_i.
  * Derivatives are taken from the SCF dipole (Hellmann-Feynman, exact for a variational SCF), one order lower
    than from the energy:
        alpha_ij = [mu_i(+F e_j) - mu_i(-F e_j)] / 2F
        beta_ijj = [mu_i(+F e_j) + mu_i(-F e_j) - 2 mu_i(0)] / F^2
        beta_ijk = [mu_i(+F,+F) - mu_i(+F,-F) - mu_i(-F,+F) + mu_i(-F,-F)] / 4F^2     (fields along j and k)
    Each has a leading error O(F^2) (gamma and delta terms), removed by one Richardson (Romberg) step with the
    same formula at 2F:  X = [4 X(F) - X(2F)] / 3.  Fields: 0, +-F, +-2F on x, y, z, and (+-F, +-F),
    (+-2F, +-2F) on the xy, xz and yz planes -- 37 SCFs at the default F = 0.001 a.u. (2F = 0.002 a.u.).
  * The static beta is a third derivative of E, so it is fully symmetric; symmetry is NOT imposed on the raw
    tensor, its largest asymmetry is reported as a noise indicator, and the symmetrised tensor is used for the
    invariants.  beta_i = sum_j beta_ijj, beta_tot = |beta_vec|, beta along mu = beta_vec . mu/|mu|; also
    beta_par = (3/5) beta_vec . mu/|mu| (the EFISH-type convention), all clearly labelled in the JSON.  beta_tot is
    also given in the "B" convention (= T/2), because literature values often use it without saying so.
  * Checks written to the JSON: alpha_ii from the energy route -[E(+F)+E(-F)-2E0]/F^2 against the dipole route,
    and the analytic (coupled-perturbed Kohn-Sham) alpha at zero field.  CPKS is implemented below with the
    algorithm of pyscf-properties' prop/polarizability/rhf.py (orbital response via mf.gen_response), because
    pyscf.prop is not installed; if it ever is importable, its alpha is recorded as a second opinion.  The CPKS
    first-order density also seeds every field's SCF (dm0 + F.dm1).
  * The SCF is converged on the orbital gradient, NLO_CONV_GRAD (default 3e-9), with the energy threshold at
    1e-10 Ha: the dipole error is LINEAR in the density error, and beta comes from dipole differences divided
    by F^2 = 1e-6, so PySCF's usual 1e-7-ish gradient left ~0.1 a.u. of noise in beta (a 1% asymmetry of the
    urea tensor in the code test); at 3e-9 repeated runs agree to ~0.02 a.u. (0.1% of urea's beta_tot), which is
    also the size of the reported tensor asymmetry.  An energy threshold of 1e-12 Ha is below the numerical noise
    of E (the SCF then idles for ~10 cycles at |g| ~ 1e-9), hence 1e-10.

PCM (documented choice): the external field acts on the solute only, and the IEF-PCM surface charges relax
self-consistently with the field-perturbed density (equilibrium, i.e. static, solvation at eps = 78.36).  This is
the solute's screened response in solution, as Gaussian's SCRF + finite-field Polar computes it.  No cavity- or
local-field factors are applied, so PCM values are molecular properties in solution, not bulk susceptibilities.
The analytic alpha check in PCM uses the same equilibrium solvent response (with_solvent.equilibrium_solvation).
The primary NLO numbers for the paper are gas-phase (M4, and M1 at its PCM geometry via --geom-from), plus the
PCM value; molecular beta says nothing about the crystal's chi(2) by itself (see the crystal-frame tensor below).

Frames: the tensors are given in the frame of the geometry file (for the crystal-derived models the optimised
structure stays close to the CIF Cartesian frame), and also rotated onto the input structure (the CIF
coordinates) by a Kabsch fit of the heavy atoms, so an oriented-gas estimate for the Pc lattice can be built later.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from props_common import *          # noqa: F401,F403 -- pins threads / scratch, so it must precede numpy and PySCF
from pyscf import lib
from pyscf.scf import cphf

# ---- unit conversions (from CODATA-2018 constants in Gaussian units, computed rather than copied) ----
A0_CM = 0.529177210903e-8                   # Bohr radius, cm
E_ESU = 4.80320471e-10                      # elementary charge, statC
EH_ERG = 4.3597447222071e-11                # hartree, erg
AU_ALPHA_ESU = A0_CM ** 3                   # a0^3 = 1.4818e-25 cm^3 (= esu)
AU_BETA_ESU = E_ESU ** 3 * A0_CM ** 3 / EH_ERG ** 2      # e^3 a0^3 / Eh^2 = 8.6392e-33 esu
AU_ALPHA_A3 = (A0_CM * 1e8) ** 3            # 0.14818 A^3

ap = parser(__doc__)
ap.add_argument('--field', type=float, default=0.001, help='base field strength F in a.u. (2F is also used)')
ap.add_argument('--no-analytic', action='store_true', help='skip the pyscf.prop second opinion (own CPKS always runs)')
args = ap.parse_args()
t_start = time.time()
case, xyz, outdir, tag = setup(args)
if case['spin'] != 0:
    raise SystemExit('[nlo] closed-shell cases only')

mol = mol_for_case(case, xyz, verbose=args.verbose)
Z = mol.atom_charges().astype(float); R = mol.atom_coords()
origin = Z @ R / Z.sum()                     # centre of nuclear charge: alpha and beta do not depend on it,
with mol.with_common_orig(origin):           # mu does for an ion (reported with the origin)
    rint = mol.intor_symmetric('int1e_r', comp=3)
nuc_dip = Z @ (R - origin)

mf = scf_for(case, mol, args, conv_tol=1e-10)
mf.conv_tol_grad = float(os.environ.get('NLO_CONV_GRAD', '3e-9'))
mf.max_cycle = 200
h0 = mf.get_hcore()
F1 = args.field; F2 = 2 * F1
print(f"[nlo] {case['name']}: {mol.natm} atoms, {mol.nao} basis functions, {environment_label(case)}; "
      f"F = {F1} and {F2} a.u.; geometry {xyz}", flush=True)


def dipole(dm):
    dmt = dm[0] + dm[1] if getattr(dm, 'ndim', 2) == 3 else dm
    return nuc_dip - np.einsum('xij,ji->x', rint, dmt).real


records = {}
dm_zero = None          # zero-field density
dm_resp = None          # first-order density per unit field (CPKS), 3 x nao x nao


def at_field(F):
    """SCF in the static field F (3-vector, a.u.); returns (E_total incl. the nuclear-field term, mu)."""
    key = tuple(np.round(F, 8))
    if key in records:
        return records[key]
    F = np.asarray(F, float)
    h1 = h0 + np.einsum('x,xij->ij', F, rint)
    mf.get_hcore = lambda *a, **k: h1
    # Initial guess: zero-field density plus the linear (CPKS) response, dm0 + F.dm1, which is off by O(F^2)
    # only.  It saves a few cycles; most of the 10-20 cycles per field go into reaching |g| < 3e-9 from anywhere.
    guess = None if dm_zero is None else (dm_zero if dm_resp is None else dm_zero + np.einsum('x,xij->ij', F, dm_resp))
    t = time.time()
    converge(mf, dm0=guess, what=f'SCF at F = {F}')
    e = float(mf.e_tot - F @ nuc_dip)
    mu = dipole(mf.make_rdm1())
    records[key] = (e, mu)
    cycles[key] = getattr(mf, 'cycles', None)
    print(f'[nlo] F = ({F[0]:+.4f},{F[1]:+.4f},{F[2]:+.4f})  E = {e:.12f}  mu = {mu.round(7)}  '
          f'cycles {cycles[key]}  {time.time() - t:.0f} s', flush=True)
    return records[key]


cycles = {}
e0, mu0 = at_field(np.zeros(3))
dm_zero = mf.make_rdm1()
mo0 = (mf.mo_energy.copy(), mf.mo_coeff.copy(), mf.mo_occ.copy())
E = np.eye(3)

# ---- analytic alpha (CPKS) at zero field: independent check, and the first-order density for the guesses ----
# The static response in PCM includes the equilibrium solvent response, i.e. the same physics as the
# finite-field runs (where the reaction field relaxes with the perturbed density).
t = time.time()
mo_e, C, occ = mo0
o = occ > 0; orbo, orbv = C[:, o], C[:, ~o]
ws = getattr(mf, 'with_solvent', None)
with (lib.temporary_env(ws, equilibrium_solvation=True) if ws is not None else contextlib.nullcontext()):
    hvo = lib.einsum('xpq,pa,qi->xai', rint, orbv, orbo)
    vresp = mf.gen_response(hermi=1)

    def to_dm(x):
        dm1 = lib.einsum('xai,pa,qi->xpq', x.reshape(-1, *hvo.shape[1:]), orbv, orbo)
        return (dm1 + dm1.transpose(0, 2, 1)) * 2          # symmetric, x2 for double occupancy

    def fvind(x):                                      # the Krylov solver may pass a subset of the 3 columns
        return lib.einsum('xpq,pa,qi->xai', vresp(to_dm(x)), orbv, orbo).ravel()

    # residual 1e-8: alpha to ~1e-7 relative, far below the finite-field comparison; 1e-10 cost 4 SCFs' worth of
    # response builds for sulfanilamide (385 s on the VM) for no visible gain
    u = cphf.solve(fvind, mo_e, occ, hvo, max_cycle=100, tol=1e-8)[0]
e2 = lib.einsum('xai,yai->xy', hvo, u)
a_an = -(e2 + e2.T) * 2
dm_resp = to_dm(u)                                     # d(dm)/dF_j for h = h_core + F.r
solver = 'own CPKS in dft/nlo.py (algorithm of pyscf-properties prop/polarizability/rhf.py; pyscf.prop not installed)'
a_prop = None
if not args.no_analytic:
    try:                                               # second opinion if the pyscf-properties package is present
        from pyscf.prop.polarizability import rks as _polmod
        with (lib.temporary_env(ws, equilibrium_solvation=True) if ws is not None else contextlib.nullcontext()):
            with mol.with_common_orig(origin):
                a_prop = np.asarray(_polmod.Polarizability(mf).polarizability())
        solver = 'own CPKS in dft/nlo.py, cross-checked with pyscf.prop.polarizability.rks'
    except ImportError:
        pass
print(f'[nlo] CPKS alpha_iso = {np.trace(a_an) / 3:.5f} a.u. ({time.time() - t:.0f} s)', flush=True)


def derivs(h):
    """alpha(h), beta(h) by central differences of the dipole at field step h (no extrapolation)."""
    a = np.zeros((3, 3)); b = np.zeros((3, 3, 3)); ediag = np.zeros(3)
    for j in range(3):
        ep, mp = at_field(+h * E[j]); em, mm = at_field(-h * E[j])
        a[:, j] = (mp - mm) / (2 * h)
        b[:, j, j] = (mp + mm - 2 * mu0) / h ** 2
        ediag[j] = -(ep + em - 2 * e0) / h ** 2
    for j, k in ((0, 1), (0, 2), (1, 2)):
        _, mpp = at_field(h * (E[j] + E[k])); _, mpm = at_field(h * (E[j] - E[k]))
        _, mmp = at_field(h * (-E[j] + E[k])); _, mmm = at_field(-h * (E[j] + E[k]))
        b[:, j, k] = b[:, k, j] = (mpp - mpm - mmp + mmm) / (4 * h ** 2)
    return a, b, ediag


a1, b1, ed1 = derivs(F1)
a2, b2, ed2 = derivs(F2)
alpha_raw = (4 * a1 - a2) / 3
beta_raw = (4 * b1 - b2) / 3
alpha_e = (4 * ed1 - ed2) / 3
alpha = 0.5 * (alpha_raw + alpha_raw.T)
perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
beta = sum(beta_raw.transpose(p) for p in perms) / 6
del mf.get_hcore


def alpha_invariants(a):
    iso = np.trace(a) / 3
    aniso = np.sqrt(0.5 * ((a[0, 0] - a[1, 1]) ** 2 + (a[1, 1] - a[2, 2]) ** 2 + (a[2, 2] - a[0, 0]) ** 2
                           + 6 * (a[0, 1] ** 2 + a[0, 2] ** 2 + a[1, 2] ** 2)))
    return float(iso), float(aniso)


def beta_invariants(b, mu):
    bvec = np.array([sum(b[i, j, j] for j in range(3)) for i in range(3)])
    btot = float(np.linalg.norm(bvec))
    mun = np.linalg.norm(mu)
    along = float(bvec @ mu / mun) if mun > 1e-8 else None
    return bvec, btot, along


iso, aniso = alpha_invariants(alpha)
bvec, btot, bmu = beta_invariants(beta, mu0)
iso1, _ = alpha_invariants(0.5 * (a1 + a1.T))
_, btot1, _ = beta_invariants(sum(b1.transpose(p) for p in perms) / 6, mu0)
analytic = {'solver': solver, 'solvent_response': 'equilibrium (static)' if ws is not None else None,
            'alpha_tensor_au': a_an, 'alpha_iso_au': float(np.trace(a_an) / 3),
            'max_abs_diff_vs_finite_field_au': float(abs(a_an - alpha).max()),
            'rel_diff_iso': float((np.trace(a_an) / 3 - iso) / iso),
            'pyscf_prop_alpha_tensor_au': a_prop}
print(f"[nlo] analytic alpha_iso = {np.trace(a_an) / 3:.5f} a.u. vs finite field {iso:.5f} "
      f"(max |diff| {abs(a_an - alpha).max():.2e})", flush=True)

# ---- rotation onto the input (CIF-frame) structure ---------------------------------------------------------
frame = None
s_in, x_in = read_xyz(ROOT / case['xyz'])
s_g, x_g = read_xyz(xyz)
heavy = [i for i, s in enumerate(s_g) if s.upper() != 'H']
if len(heavy) >= 3:
    P = x_g[heavy] - x_g[heavy].mean(0); Q = x_in[heavy] - x_in[heavy].mean(0)
    U, S, Vt = np.linalg.svd(P.T @ Q)
    dsgn = np.sign(np.linalg.det(U @ Vt))
    Rm = (U @ np.diag([1, 1, dsgn]) @ Vt).T               # Q ~ P @ Rm.T, i.e. v_in = Rm @ v_geom
    rmsd = float(np.sqrt(((P @ Rm.T - Q) ** 2).sum(1).mean()))
    frame = {'rotation_geom_to_input': Rm, 'heavy_atom_rmsd_A': rmsd,
             'mu_au': Rm @ mu0, 'alpha_tensor_au': Rm @ alpha @ Rm.T,
             'beta_tensor_au': np.einsum('ia,jb,kc,abc->ijk', Rm, Rm, Rm, beta),
             'beta_vec_au': Rm @ bvec,
             'note': ('tensors rotated onto the input structure ' + case['xyz'] + ' by a Kabsch fit of the heavy '
                      'atoms; for the crystal-derived models this is the CIF Cartesian frame (gemmi orthogonalisation: '
                      'x || a, z || c*), for urea it is the RDKit frame')}

res = {
    'provenance': provenance(case, args, xyz, t_start, {'level': f"{(case.get('xc') or XC).upper()}"
                                                                   f"{'-' + case['disp'].upper() if case.get('disp') not in (None, 'none') else ''}"
                                                                   f"/{basis_description(case)}"}),
    'method': ('finite field on the SCF dipole, fields 0, +-F, +-2F per axis and (+-F,+-F), (+-2F,+-2F) per plane; '
               'one Richardson step [4X(F)-X(2F)]/3; Taylor convention E = E0 - mu F - 1/2 alpha F F - 1/6 beta F F F; '
               'PCM: field on the solute only, equilibrium (static) solvent response'),
    'field_au': [F1, F2], 'n_scf': len(records), 'gauge_origin_bohr': origin,
    'units': {'alpha_au_to_esu': AU_ALPHA_ESU, 'alpha_au_to_A3': AU_ALPHA_A3, 'beta_au_to_esu': AU_BETA_ESU,
              'mu_au_to_debye': AU_DEBYE},
    'E0_Ha': e0,
    'mu_au': mu0, 'mu_debye': mu0 * AU_DEBYE, 'mu_total_au': float(np.linalg.norm(mu0)),
    'mu_total_debye': float(np.linalg.norm(mu0) * AU_DEBYE),
    'alpha_tensor_au': alpha, 'alpha_iso_au': iso, 'alpha_aniso_au': aniso,
    'alpha_iso_esu': iso * AU_ALPHA_ESU, 'alpha_aniso_esu': aniso * AU_ALPHA_ESU, 'alpha_iso_A3': iso * AU_ALPHA_A3,
    'alpha_max_asymmetry_au': float(abs(alpha_raw - alpha_raw.T).max()),
    'alpha_diag_energy_route_au': alpha_e, 'alpha_diag_energy_minus_dipole_au': alpha_e - np.diag(alpha),
    'beta_tensor_au': beta, 'beta_tensor_raw_au': beta_raw,
    'beta_max_asymmetry_au': float(max(abs(beta_raw - beta_raw.transpose(p)).max() for p in perms)),
    'beta_vec_au': bvec, 'beta_tot_au': btot, 'beta_tot_esu': btot * AU_BETA_ESU,
    'beta_along_mu_au': bmu, 'beta_along_mu_esu': None if bmu is None else bmu * AU_BETA_ESU,
    'beta_par_3_5_au': None if bmu is None else 0.6 * bmu,
    'beta_par_3_5_esu': None if bmu is None else 0.6 * bmu * AU_BETA_ESU,
    # The other common convention ("B", perturbation series mu = mu0 + alpha F + beta F F) halves every beta.
    # Literature values are often quoted without saying which one they use (Willetts et al., J. Chem. Phys. 97
    # (1992) 7590).  The widely cited urea beta_tot = 0.3728e-30 esu at "B3LYP/6-311++G(d,p)" happens to lie within
    # 4% of HALF our T-convention urea value, but urea's beta also depends on its geometry (planar vs pyramidal NH2)
    # and the convention behind that number is unverified, so the closeness proves nothing about conventions: both
    # are written out and any comparison must name the convention.  [review 28 Sep 2026]
    'beta_convention': 'T (Taylor series: mu = mu0 + alpha F + 1/2 beta F F); all beta_* fields above',
    'beta_tot_B_convention_au': btot / 2, 'beta_tot_B_convention_esu': btot / 2 * AU_BETA_ESU,
    'richardson_effect': {'alpha_iso_F_only_au': iso1, 'beta_tot_F_only_au': btot1,
                          'alpha_iso_extrapolated_au': iso, 'beta_tot_extrapolated_au': btot},
    'analytic_alpha_check': analytic,
    'input_frame': frame,
    'fields': [{'F_au': list(k), 'E_Ha': v[0], 'mu_au': v[1], 'scf_cycles': cycles.get(k)} for k, v in records.items()],
    'note': ('static (zero-frequency) electronic properties; no vibrational contributions, no frequency dispersion. '
             + ('PCM: solute response in solution, no local-field factors. ' if case.get('solvent') else '')
             + ('Charged system: mu depends on the origin (centre of nuclear charge used); alpha and beta do not. '
                if case['charge'] else '')),
}
name = out_name('nlo', tag)
save_json(outdir / name, res)
rows = [('mu', f"{np.linalg.norm(mu0):.4f}", f"{np.linalg.norm(mu0) * AU_DEBYE:.3f} D"),
        ('alpha_iso', f'{iso:.3f}', f'{iso * AU_ALPHA_ESU:.4e} esu'),
        ('Delta alpha', f'{aniso:.3f}', f'{aniso * AU_ALPHA_ESU:.4e} esu'),
        ('beta_tot (T convention)', f'{btot:.2f}', f'{btot * AU_BETA_ESU:.4e} esu'),
        ('beta_tot (B convention = T/2)', f'{btot / 2:.2f}', f'{btot / 2 * AU_BETA_ESU:.4e} esu'),
        ('beta along mu', '-' if bmu is None else f'{bmu:.2f}', '-' if bmu is None else f'{bmu * AU_BETA_ESU:.4e} esu')]
(outdir / out_name('nlo', tag, '.md')).write_text(
    f"Static NLO, {case['name']} ({environment_label(case)}), geometry {res['provenance']['geometry']}\n\n"
    + md_table(rows, ['Property', 'a.u.', 'other units']))
if outdir == case['dir']:
    stage_note(case, f"NLO ({name}): mu {np.linalg.norm(mu0) * AU_DEBYE:.3f} D, alpha_iso {iso:.2f} a.u., "
                     f"beta_tot {btot:.1f} a.u. ({btot * AU_BETA_ESU:.3e} esu)")
print(f"[nlo] mu = {np.linalg.norm(mu0) * AU_DEBYE:.4f} D; alpha_iso = {iso:.4f} a.u. ({iso * AU_ALPHA_ESU:.4e} esu), "
      f"Delta alpha = {aniso:.4f} a.u.; beta_tot = {btot:.3f} a.u. ({btot * AU_BETA_ESU:.4e} esu); "
      f"beta along mu = {bmu if bmu is None else round(bmu, 3)} a.u.; energy-route alpha_ii - dipole route = "
      f"{(alpha_e - np.diag(alpha)).round(4)}; beta asymmetry {res['beta_max_asymmetry_au']:.2e} a.u.; "
      f"{len(records)} SCFs, {(time.time() - t_start) / 60:.1f} min -> {outdir / name}", flush=True)
