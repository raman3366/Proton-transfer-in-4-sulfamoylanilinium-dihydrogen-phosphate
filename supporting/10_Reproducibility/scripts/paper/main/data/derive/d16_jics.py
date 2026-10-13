"""d16: quantities added for the version of 5 Oct 2026 (Journal of the Indian Chemical Society), re-derived from the
raw result files with code of our own (derive/ convention; nothing of dft/ is imported).

1-D proton model   The model of dft/zpe_barrier.py rebuilt independently: lower envelope of the two relaxed scan chains
                   (scan/[AB]_*.json), PCHIP slopes set to zero at the two unconstrained minima, cubic Hermite
                   interpolation in d(N1-H1A), continuation beyond the scanned range by harmonic walls, sinc-DVR
                   (Colbert-Miller) on 801 points over the range plus 0.6 A on each side, N-H reduced mass.  The
                   reference run uses the wall constant stored in zpe_barrier.json (the N1-H1A local-mode force
                   constant of the salt) on both sides and must reproduce the stored level, its distance from the
                   barrier top and the O5-side probability.  Sensitivity: the O5-side wall (where the relaxed
                   coordinate is the soft N...O stretch, so the N-H constant has no physical basis there), the
                   N1-side wall, and the effective mass.
Zero-point energy  of the separated species: ZPE(sulfanilamide) + ZPE(H3PO4) - ZPE(cation) - ZPE(H2PO4-) from the
                   gas-phase Hessians (freq.json of each case), the correction that the continuum exchange energy
                   of Section 4.3 does not contain.
def2-TZVP          shift of the proton coordinate delta = d(N1-H1A) - d(H1A...O5) between the 6-311++G(d,p) and
                   def2-TZVP optima of the ion pair in IEF-PCM, and the P1-O4 length against the X-ray value.
TD-DFT anchor      S2 of neutral sulfanilamide in IEF-PCM (CAM-B3LYP, non-equilibrium solvation) against the
                   wavelength at which aqueous sulfanilamide is determined by its absorption in the literature (258 nm).
"""
import glob

from scipy.interpolate import CubicHermiteSpline, PchipInterpolator

from common import *

AMU_ME = 1822.888486209          # atomic mass unit in electron masses
MDYN_A_AU = 15.569141            # 1 Eh/bohr^2 in mdyn/A
EV_NM = 1239.841984              # h c in eV nm
LIT_UV_NM = 258.0                # aqueous sulfanilamide, wavelength of its UV determination (literature; see facts_jics.py)

out = {}

# ------------------------------------------------------------------------------------------------- 1-D proton model
SCAN = RES / 'sfa_pair_pcm_d3' / 'scan'
ZB = json.loads((SCAN / 'zpe_barrier.json').read_text())
pts = []
for f in sorted(glob.glob(str(SCAN / '[AB]_*.json'))):
    r = json.loads(Path(f).read_text())
    pts.append((r['metrics']['d_N1_H1A_A'], r['E_Ha'], r['kind'] == 'anchor'))
env = {}
for d_, e_, a_ in pts:                                     # lower envelope of the two chains
    k = round(d_, 3)
    if k not in env or e_ < env[k][1]:
        env[k] = (d_, e_, a_)
env = sorted(env.values())
d = np.array([e[0] for e in env])
E = np.array([e[1] for e in env])
anchor = np.array([e[2] for e in env])
V = E - E[0]
x = d / BOHR_A


def solve(mu_amu, k_left, k_right, n=801, pad=0.6):
    """(E0, E0 - barrier top, P on the O5 side) in kcal/mol; wall constants in mdyn/A."""
    kl, kr = k_left / MDYN_A_AU, k_right / MDYN_A_AU
    sl = PchipInterpolator(x, V).derivative()(x)
    sl[anchor] = 0.0
    spl = CubicHermiteSpline(x, V, sl)

    def pot(xx):
        xx = np.asarray(xx, float)
        v = spl(np.clip(xx, x[0], x[-1]))
        lo, hi = xx < x[0], xx > x[-1]
        v[lo] += sl[0] * (xx[lo] - x[0]) + 0.5 * kl * (xx[lo] - x[0]) ** 2
        v[hi] += sl[-1] * (xx[hi] - x[-1]) + 0.5 * kr * (xx[hi] - x[-1]) ** 2
        return v
    xg = np.linspace(x[0] - pad / BOHR_A, x[-1] + pad / BOHR_A, n)
    dx = xg[1] - xg[0]
    i = np.arange(n)
    dij = i[:, None] - i[None, :]
    with np.errstate(divide='ignore'):
        T = np.where(dij == 0, np.pi ** 2 / 3, 2.0 / np.where(dij == 0, 1, dij) ** 2) * (-1.0) ** dij
    T = T / (2 * mu_amu * AMU_ME * dx ** 2)
    ev, U = np.linalg.eigh(T + np.diag(pot(xg)))
    xf = np.linspace(x[0], x[-1], 4000)
    vf = pot(xf)
    it = int(np.argmax(vf))
    p = U[:, 0] ** 2
    return float(ev[0] * HA_KCAL), float((ev[0] - vf[it]) * HA_KCAL), float(p[xg > xf[it]].sum())


MU = ZB['effective_mass_amu']
KW = ZB['one_d_model']['wall_k_mdyn_per_A']
ref = solve(MU, KW, KW)
stored = (ZB['one_d_model']['E0_above_reference_kcal'], ZB['one_d_model']['E0_minus_barrier_top_kcal'],
          ZB['one_d_model']['levels'][0]['P_O_side'])
oned = dict(mass_amu=MU, wall_k_mdyn_per_A=KW, d_first_A=float(d[0]), d_last_A=float(d[-1]), n_points=len(d),
            grid_points=801, pad_A=0.6,
            reference=dict(E0=ref[0], E0_minus_top=ref[1], P_O=ref[2]),
            stored=dict(E0=stored[0], E0_minus_top=stored[1], P_O=stored[2]),
            max_abs_diff_vs_stored=max(abs(a - b) for a, b in zip(ref, stored)))
RIGHT = [0.5, 1.0, 5.0]                                    # O5-side wall, mdyn/A (softer walls need a wider grid)
oned['right_wall'] = {f'{k:g}': dict(zip(('E0', 'E0_minus_top', 'P_O'), solve(MU, KW, k, n=1201, pad=1.5 if k < 0.5 else 0.6)))
                      for k in RIGHT}
LEFT = [1.0, 5.0]
oned['left_wall'] = {f'{k:g}': dict(zip(('E0', 'E0_minus_top', 'P_O'), solve(MU, k, KW))) for k in LEFT}
MASS = [1.008, 1.5, 2.0]
oned['mass'] = {f'{m:g}': dict(zip(('E0', 'E0_minus_top', 'P_O'), solve(m, KW, KW))) for m in MASS}
allv = [oned['reference']] + list(oned['right_wall'].values()) + list(oned['left_wall'].values()) + \
    list(oned['mass'].values())
oned['E0_minus_top_min'] = min(v['E0_minus_top'] for v in allv)
oned['E0_minus_top_max'] = max(v['E0_minus_top'] for v in allv)
oned['P_O_min'] = min(v['P_O'] for v in allv)
oned['P_O_max'] = max(v['P_O'] for v in allv)
oned['n_variants'] = len(allv) - 1
out['oned'] = oned

# -------------------------------------------------------------------------- zero-point energy of the separated species


def zpe(case):
    return float(json.loads((RES / case / 'freq.json').read_text())['ZPE_Ha'])


Z = {c: zpe(c) for c in ('sulfanilamide_gas_d3', 'h3po4_gas_d3', 'cation_gas_d3', 'h2po4_gas_d3')}
dz = Z['sulfanilamide_gas_d3'] + Z['h3po4_gas_d3'] - Z['cation_gas_d3'] - Z['h2po4_gas_d3']
out['zpe_separated'] = dict(cases=Z, dZPE_Ha=dz, dZPE_kcal=dz * HA_KCAL)

# ---------------------------------------------------------------------------------------------- def2-TZVP geometry


def geom(case):
    syms, xyz, _ = read_xyz(RES / case / 'opt.xyz')
    dist = lambda a, b: float(np.linalg.norm(xyz[IDX[a]] - xyz[IDX[b]]))                   # noqa: E731
    return dict(d_NH=dist('N1', 'H1A'), d_HO=dist('H1A', 'O5'), d_NO=dist('N1', 'O5'), P1_O4=dist('P1', 'O4'),
                delta=dist('N1', 'H1A') - dist('H1A', 'O5'))


g0, g1 = geom('sfa_pair_pcm_d3'), geom('sfa_pair_pcm_d3_tz')
out['tz'] = dict(main=g0, tz=g1, delta_shift=g1['delta'] - g0['delta'])

# ------------------------------------------------------------------------------------------------- TD-DFT anchor
out['uv'] = dict(lit_nm=LIT_UV_NM, lit_eV=EV_NM / LIT_UV_NM)

dump('d16_jics.json', out)
print(f"1-D model: E0 {ref[0]:.3f}, E0 - top {ref[1]:.3f}, P(O5) {ref[2]:.3f} (stored {stored[0]:.3f}, {stored[1]:.3f}, "
      f"{stored[2]:.3f}; max diff {oned['max_abs_diff_vs_stored']:.1e}); variants: E0 - top "
      f"{oned['E0_minus_top_min']:.2f} to {oned['E0_minus_top_max']:.2f}, P(O5) {oned['P_O_min']:.2f} to {oned['P_O_max']:.2f}")
print(f"dZPE separated {dz * HA_KCAL:.3f} kcal/mol; delta shift def2-TZVP {out['tz']['delta_shift']:.4f} A; "
      f"P1-O4 TZ {g1['P1_O4']:.4f} A; 258 nm = {EV_NM / LIT_UV_NM:.4f} eV")
