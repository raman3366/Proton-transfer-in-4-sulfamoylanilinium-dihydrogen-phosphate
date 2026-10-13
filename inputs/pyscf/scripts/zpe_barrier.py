"""
Zero-point energy of the N1-H1A stretch against the proton-transfer barrier of the relaxed scan (1-D estimate).

Why: at delta-pKa ~ 0 a relaxed scan can show a shallow second well on the phosphate side, or a shallow barrier
in front of it.  A classical (electronic-energy) barrier smaller than the zero-point energy of the N-H stretch
does not localise the proton: the ground vibrational level then sits above the barrier top.  The referee
critique asks for exactly this comparison, so it is made here explicitly, with the approximation stated.

Usage: python dft/zpe_barrier.py <case> [--scan-dir DIR] [--hessian H.npy --xyz X.xyz] [--freq freq.json]
                                        [--mass nh|h]
Inputs (defaults): results/<case>/scan/ (dft/scan.py per-point JSONs), results/<case>/hessian.npy and opt.xyz
(dft/frequencies.py; for sfa_pair_pcm_d3 the finite-difference Hessian of the IEF-PCM minimum) and freq.json
(cross-check of the frequencies).  Outputs: <scan dir>/zpe_barrier.json and zpe_barrier.png.

Three numbers for the N-H stretch, from the same Hessian:
  1. the normal mode with the largest N1-H1A stretch character c = mu_NH (dr/dQ)^2 (1 for a pure N-H stretch,
     0 for none); its harmonic ZPE 1/2 h nu is the headline number (the local-mode value below replaces it when
     the stretch is spread over several modes, character < 0.5);
  2. the Konkoli-Cremer adiabatic local mode, k_a = 1 / (b^T H^+ b) with H^+ the pseudo-inverse of the
     vibrational Hessian: the N-H force constant with every other coordinate relaxed, i.e. exactly the
     curvature a relaxed scan along d(N1-H1A) has at the minimum -- so it can be checked against the scan;
  3. that scan curvature itself (cubic fit of the first scan points, anharmonic).
Then a 1-D quantum model: the lower envelope of the scan (both chains) interpolated with a monotone cubic
(PCHIP), extended beyond the scanned range by harmonic walls of force constant k_a, and solved with a
Colbert-Miller sinc DVR for the proton's vibrational levels.

APPROXIMATIONS (also written to the JSON):
  * one dimension: the proton coordinate is d(N1-H1A); every other coordinate follows the relaxed scan
    (adiabatic, minimum-energy path).  The ZPE of the other 3N-7 modes -- notably the two H1A bends, which
    soften as the proton is shared -- is not followed along the path, and no tunnelling beyond this 1-D
    model is considered;
  * effective mass = the N-H reduced mass (0.94 u; --mass h gives the bare proton, 1.008 u).  The heavy-atom
    relaxation along the path adds inertia that is ignored, so this is the light end: it gives the largest
    ZPE and the most delocalised proton, i.e. it is conservative for a claim that the proton stays on N1;
  * harmonic normal-mode ZPE (no anharmonicity; H-bonded N-H stretches are strongly anharmonic);
  * electronic energies of the relaxed scan (B3LYP-D3(BJ), equilibrium IEF-PCM at every point; B3LYP tends
    to underestimate proton-transfer barriers), and the barrier is only as fine as the 0.10 A grid.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from proton_tools import *     # noqa: F401,F403
import argparse
from pyscf.data import elements
from scipy.interpolate import PchipInterpolator, CubicHermiteSpline, CubicSpline

AU2CM = 219474.6313702
AMU2ME = 1822.888486209
KCAL_PER_CM = 1.0 / 349.7550882
MDYNA_PER_AU = 15.569141            # 1 Ha/bohr^2 in mdyn/A

# ------------------------------------------------------------------------------------- Hessian analysis ---
def normal_modes(syms, xyz_A, hess):
    """Frequencies (cm-1), Cartesian displacements per unit mass-weighted normal coordinate (3N x nvib, in
    bohr / sqrt(m_e)) and mass-weighted eigenvalues.  Same translation/rotation projection as vibrations.py."""
    n = len(syms)
    m = np.array([elements.MASSES[elements.charge(s)] for s in syms]) * AMU2ME
    x = np.asarray(xyz_A) / BOHR_A
    Hc = hess.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)
    w = np.repeat(m ** -0.5, 3)
    Hm = Hc * w[:, None] * w[None, :]; Hm = 0.5 * (Hm + Hm.T)
    r = x - (m[:, None] * x).sum(0) / m.sum()
    TR = np.zeros((6, 3 * n))
    for i in range(n):
        s = np.sqrt(m[i])
        TR[0, 3 * i] = TR[1, 3 * i + 1] = TR[2, 3 * i + 2] = s
        for k, ax in enumerate(np.eye(3)):
            TR[3 + k, 3 * i:3 * i + 3] = s * np.cross(ax, r[i])
    q, _ = np.linalg.qr(TR.T)
    P = np.eye(3 * n) - q @ q.T
    ev, vec = np.linalg.eigh(P @ Hm @ P)
    f = np.sign(ev) * np.sqrt(np.abs(ev)) * AU2CM
    keep = np.argsort(np.abs(f))[6:]
    keep = keep[np.argsort(f[keep])]
    return f[keep], vec[:, keep] * w[:, None], ev[keep], m


def nh_analysis(syms, xyz_A, hess, mu, N, H):
    f, C, ev, m = normal_modes(syms, xyz_A, hess)
    x = np.asarray(xyz_A) / BOHR_A
    e = (x[H] - x[N]) / np.linalg.norm(x[H] - x[N])
    b = np.zeros(3 * len(syms)); b[3 * H:3 * H + 3] = e; b[3 * N:3 * N + 3] = -e
    dr = b @ C                                   # bohr per unit mass-weighted normal coordinate
    char = mu * dr ** 2                          # in [0, 1]; = 1 for a pure N-H stretch
    order = np.argsort(-char)
    k_best = int(order[0])
    compliance = float(np.sum(dr ** 2 / ev))     # b^T H^+ b (bohr^2 / Ha)
    k_a = 1.0 / compliance
    nu_a = np.sqrt(k_a / mu) * AU2CM if k_a > 0 else float('nan')
    return {'n_imaginary': int((f < 0).sum()),
            'normal_mode': {'index': k_best, 'freq_cm-1': float(f[k_best]), 'nh_stretch_character': float(char[k_best]),
                            'zpe_half_h_nu_kcal': float(0.5 * f[k_best] * KCAL_PER_CM),
                            'top3': [{'freq_cm-1': float(f[i]), 'nh_stretch_character': float(char[i])}
                                     for i in order[:3]]},
            'local_mode': {'k_a_mdyn_per_A': float(k_a * MDYNA_PER_AU), 'freq_cm-1': float(nu_a),
                           'zpe_half_h_nu_kcal': float(0.5 * nu_a * KCAL_PER_CM),
                           'note': 'Konkoli-Cremer adiabatic local mode (compliance) of the N1-H1A bond'},
            'frequencies_cm-1': f}


# ---------------------------------------------------------------------------------------- 1-D model -----
def dvr_levels(V, a, b, mu, n=801, nlev=6):
    """Colbert-Miller sinc DVR on [a, b] (bohr); V(x) in Ha; mu in m_e.  Returns levels (Ha), grid, states."""
    xg = np.linspace(a, b, n); dx = xg[1] - xg[0]
    i = np.arange(n); dij = i[:, None] - i[None, :]
    with np.errstate(divide='ignore'):
        T = np.where(dij == 0, np.pi ** 2 / 3, 2.0 / np.where(dij == 0, 1, dij) ** 2) * (-1.0) ** dij
    T = T / (2 * mu * dx ** 2)
    E, U = np.linalg.eigh(T + np.diag(V(xg)))
    return E[:nlev], xg, U[:, :nlev] / np.sqrt(dx)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('case')
    ap.add_argument('--scan-dir'); ap.add_argument('--hessian'); ap.add_argument('--xyz'); ap.add_argument('--freq')
    ap.add_argument('--mass', choices=['nh', 'h'], default='nh')
    args = ap.parse_args()
    case = load_case(args.case)
    atoms = proton_atoms(case)
    N, H = atoms['N'], atoms['H']
    SD = pathlib.Path(args.scan_dir) if args.scan_dir else case['dir'] / 'scan'
    HESS = pathlib.Path(args.hessian) if args.hessian else case['dir'] / 'hessian.npy'
    XYZ = pathlib.Path(args.xyz) if args.xyz else case['dir'] / 'opt.xyz'
    FREQ = pathlib.Path(args.freq) if args.freq else case['dir'] / 'freq.json'
    pts = load_points(SD)
    if len(pts) < 3:
        raise SystemExit(f'[zpe] need at least 3 finished scan points in {SD}, found {len(pts)}')
    env = lower_envelope(pts)
    syms, xh = read_xyz(XYZ)
    hess = np.load(HESS)
    m_H = elements.MASSES[1]; m_N = elements.MASSES[7]
    mu_amu = m_H * m_N / (m_H + m_N) if args.mass == 'nh' else m_H
    mu = mu_amu * AMU2ME
    nh = nh_analysis(syms, xh, hess, mu, N, H)
    xcheck = None
    if FREQ.exists():
        ref = np.array(json.load(open(FREQ))['frequencies_cm-1'])
        if len(ref) == len(nh['frequencies_cm-1']):
            xcheck = float(np.abs(np.sort(ref) - nh['frequencies_cm-1']).max())

    d = np.array([r['metrics']['d_N1_H1A_A'] for r in env]); E = np.array([r['E_Ha'] for r in env])
    anchors = [r['kind'] == 'anchor' for r in env]
    side = [r['metrics']['proton_on'] for r in env]
    # N-protonated minimum = the lowest local minimum of the envelope with the proton on N (normally chain A's
    # unconstrained anchor).  None when the N side runs downhill into the O-protonated form: then there is no
    # barrier to compare a ZPE with, and that is reported instead of a number.
    loc = [i for i in range(len(E)) if (i == 0 or E[i] <= E[i - 1]) and (i == len(E) - 1 or E[i] <= E[i + 1])]
    nmins = [i for i in loc if side[i] == 'N']
    iN = min(nmins, key=lambda i: E[i]) if nmins else None
    E0 = E[iN] if iN is not None else E.min()
    V = (E - E0)                                     # Ha, relative to the N minimum (else the global minimum)
    # where the proton changes sides along the envelope (delta = 0), for the side populations
    dl = np.array([r['metrics']['delta_A'] for r in env])
    cross = next((d[i] + (d[i + 1] - d[i]) * (-dl[i]) / (dl[i + 1] - dl[i])
                  for i in range(len(d) - 1) if dl[i] < 0 <= dl[i + 1]), None)

    # monotone cubic through the lower envelope; zero slope at unconstrained anchors (they are stationary)
    dB = d / BOHR_A
    slope = PchipInterpolator(dB, V).derivative()(dB)
    slope[np.array(anchors)] = 0.0
    spl = CubicHermiteSpline(dB, V, slope)
    kA = nh['local_mode']['k_a_mdyn_per_A'] / MDYNA_PER_AU
    # the walls outside the scanned range need a stiff, positive curvature; if the Hessian is not a clean minimum
    # (negative compliance) fall back to a typical X-H stretch constant (5 mdyn/A) and say so
    wall_note = 'k_a of the N1-H1A local mode'
    if not kA > 1.0 / MDYNA_PER_AU:
        kA, wall_note = 5.0 / MDYNA_PER_AU, 'fallback 5 mdyn/A (local-mode k_a not usable)'

    def Vx(x):
        x = np.asarray(x, dtype=float); out = spl(np.clip(x, dB[0], dB[-1]))
        lo, hi = x < dB[0], x > dB[-1]
        out[lo] += slope[0] * (x[lo] - dB[0]) + 0.5 * kA * (x[lo] - dB[0]) ** 2
        out[hi] += slope[-1] * (x[hi] - dB[-1]) + 0.5 * kA * (x[hi] - dB[-1]) ** 2
        return out

    # barrier: highest point of the interpolated profile beyond the N minimum and before the O-side minimum
    top = None
    if iN is None:
        prof = {'N_minimum': False, 'barrier_kcal': None, 'E_ref_Ha': float(E0),
                'note': 'no N-protonated minimum on the profile: the N side runs downhill into the O-protonated '
                        'form, so the proton is not held on N1 at this level and there is no barrier to compare '
                        'a zero-point energy with; energies are relative to the global minimum',
                'd_global_min_A': float(d[int(np.argmin(E))]), 'proton_on_at_global_min': side[int(np.argmin(E))]}
    else:
        xf = np.linspace(dB[iN], dB[-1], 2000); vf = Vx(xf)
        imax = [i for i in range(1, len(xf) - 1) if vf[i] >= vf[i - 1] and vf[i] > vf[i + 1]]
        prof = {'N_minimum': True, 'E_N_min_Ha': float(E0), 'd_N_min_A': float(d[iN]),
                'N_min_is_unconstrained_minimum': bool(anchors[iN]),
                'N_min_at_scan_edge': bool(iN == 0 and not anchors[iN])}
        if imax:
            top = imax[0]
            after = xf[top:]; vaft = vf[top:]
            jmin = int(np.argmin(vaft))
            prof.update({'barrier_kcal': float(vf[top] * KCAL_PER_HA), 'd_barrier_A': float(xf[top] * BOHR_A),
                         'O_side_min_kcal': float(vaft[jmin] * KCAL_PER_HA),
                         'd_O_side_min_A': float(after[jmin] * BOHR_A),
                         'O_side_min_is_end_of_scan': bool(jmin == len(after) - 1),
                         # the last envelope point is a real minimum only if it is an unconstrained anchor
                         'O_side_min_is_open_end': bool(jmin == len(after) - 1 and not anchors[-1]),
                         'reverse_barrier_kcal': float((vf[top] - vaft[jmin]) * KCAL_PER_HA)})
            top = xf[top]
            # PCHIP never overshoots the data, so at a smooth maximum between two grid points it underestimates
            # the barrier; a not-a-knot cubic spline through the same points is the less biased estimate.  When
            # they differ noticeably the top is under-resolved and one more point there is worth its cost.
            cs = CubicSpline(dB, V)(xf)
            j = int(np.argmax(np.where(xf <= (top + 0.15 / BOHR_A), cs, -np.inf)))
            prof.update({'barrier_cubic_spline_kcal': float(cs[j] * KCAL_PER_HA),
                         'd_barrier_cubic_spline_A': float(xf[j] * BOHR_A)})
            if abs(cs[j] - vf[np.searchsorted(xf, top)]) * KCAL_PER_HA > 0.2:
                prof['refine_suggestion'] = (f"barrier top under-resolved: run one extra constrained point, e.g. "
                                             f"dft/scan.py {case['name']} --chain A --points "
                                             f"{round(xf[j] * BOHR_A, 3):.3f} (then --report)")
        else:
            prof.update({'barrier_kcal': None,
                         'note': 'single well: no maximum beyond the N-protonated minimum, the energy rises over the '
                                 'whole scanned range, so there is no O-protonated minimum on this profile',
                         'E_end_of_scan_kcal': float(V[-1] * KCAL_PER_HA), 'd_end_of_scan_A': float(d[-1])})

    # scan curvature at the N minimum: E = E0 + k/2 x^2 + c x^3 through the next two points on the N side
    curv = None
    right = [i for i in range(iN + 1, len(d)) if d[i] - d[iN] > 0.02][:2] if iN is not None else []
    if iN is not None and anchors[iN] and len(right) == 2:
        xs = (d[right] - d[iN]) / BOHR_A
        A = np.array([[0.5 * xs[0] ** 2, xs[0] ** 3], [0.5 * xs[1] ** 2, xs[1] ** 3]])
        k_s, _ = np.linalg.solve(A, V[right])
        curv = {'k_mdyn_per_A': float(k_s * MDYNA_PER_AU),
                'freq_cm-1': float(np.sqrt(k_s / mu) * AU2CM) if k_s > 0 else None,
                'note': 'cubic fit E0 + k/2 x^2 + c x^3 through the two scan points above d0; compare with '
                        'the local-mode k_a'}

    a, b = dB[0] - 0.6 / BOHR_A, dB[-1] + 0.6 / BOHR_A
    lev, xg, U = dvr_levels(Vx, a, b, mu)
    # dividing point between the N and O sides: the barrier top if there is one, else where delta changes sign
    xdiv = top if top is not None else (cross / BOHR_A if cross is not None else None)
    levels = []
    for k in range(len(lev)):
        p = U[:, k] ** 2 * (xg[1] - xg[0])
        levels.append({'E_kcal': float(lev[k] * KCAL_PER_HA), '<d>_A': float((p * xg).sum() * BOHR_A),
                       'P_O_side': (float(p[xg > xdiv].sum()) if xdiv is not None else None)})
    zpe1d = levels[0]['E_kcal']
    # Headline harmonic ZPE: the normal mode when it is mostly the N1-H1A stretch; when that stretch is spread over
    # several modes (e.g. mixed with the other two NH3+ stretches) no single mode is "the" N-H stretch and the
    # local-mode value is the meaningful one.
    if nh['normal_mode']['nh_stretch_character'] >= 0.5:
        zpe_h, zpe_src = nh['normal_mode']['zpe_half_h_nu_kcal'], 'normal mode'
    else:
        zpe_h, zpe_src = nh['local_mode']['zpe_half_h_nu_kcal'], 'local mode (normal modes mixed)'
    if iN is None:
        verdict = (f"No N-protonated minimum on this profile: the proton sits on O5 (global minimum at d(N1-H1A) = "
                   f"{prof['d_global_min_A']:.2f} A); ZPE vs barrier does not arise. 1-D ground level "
                   f"{zpe1d:.2f} kcal/mol above the global minimum, <d> = {levels[0]['<d>_A']:.2f} A.")
    elif prof.get('barrier_kcal') is None:
        verdict = (f"Single well on N1 (no O-protonated minimum up to d = {d[-1]:.2f} A, where E is "
                   f"{prof['E_end_of_scan_kcal']:.2f} kcal/mol above the N minimum). Harmonic N-H ZPE ({zpe_src}) "
                   f"{zpe_h:.2f} kcal/mol; 1-D ground level {zpe1d:.2f} kcal/mol above the N minimum.")
    else:
        Eb = prof['barrier_kcal']
        verdict = (f"Barrier {Eb:.2f} kcal/mol (cubic-spline estimate {prof['barrier_cubic_spline_kcal']:.2f}; O-side "
                   f"{'end of the scanned range (energy still falling there, not a minimum)' if prof['O_side_min_is_open_end'] else 'minimum'} "
                   f"{prof['O_side_min_kcal']:+.2f}) vs harmonic N-H ZPE "
                   f"({zpe_src}) {zpe_h:.2f} kcal/mol (ratio ZPE/barrier {zpe_h / Eb:.2f}); 1-D ground level "
                   f"{zpe1d:.2f} kcal/mol ({'ABOVE' if zpe1d > Eb else 'below'} the barrier top), "
                   f"P(O side) = {levels[0]['P_O_side']:.3f}.")
    out = {'case': case['name'],
           'inputs': {'scan_dir': str(SD), 'hessian': str(HESS), 'hessian_geometry': str(XYZ),
                      'freq_json': str(FREQ) if FREQ.exists() else None,
                      'max_freq_diff_vs_freq_json_cm-1': xcheck, 'n_scan_points': len(pts),
                      'n_envelope_points': len(env)},
           'effective_mass_amu': float(mu_amu), 'mass_model': args.mass,
           'nh_stretch': {k: v for k, v in nh.items() if k != 'frequencies_cm-1'},
           'headline_zpe_kcal': zpe_h, 'headline_zpe_source': zpe_src,
           'scan_curvature': curv, 'profile': prof,
           'one_d_model': {'levels': levels, 'energy_reference': ('N-protonated minimum' if iN is not None
                                                                  else 'global minimum (no N minimum)'),
                           'E0_above_reference_kcal': zpe1d,
                           'E0_minus_barrier_top_kcal': (zpe1d - prof['barrier_kcal']
                                                         if prof.get('barrier_kcal') is not None else None),
                           'interpolation': 'PCHIP slopes (zero at unconstrained minima), cubic Hermite; harmonic '
                                            'walls outside the scanned range',
                           'wall_k_mdyn_per_A': float(kA * MDYNA_PER_AU), 'wall_k_source': wall_note,
                           'dvr': 'Colbert-Miller sinc DVR, 801 points, range scanned +/- 0.6 A'},
           'verdict': verdict,
           'approximations': [
               '1-D: d(N1-H1A) is the only quantum coordinate; all others follow the relaxed scan (adiabatic); the '
               'ZPE of the other 3N-7 modes (notably the H1A bends) is not followed along the path; no '
               'multidimensional tunnelling.',
               f'Effective mass {mu_amu:.3f} u ({"N-H reduced mass" if args.mass == "nh" else "proton"}); heavy-atom '
               'relaxation adds inertia that is ignored (upper bound on ZPE and delocalisation).',
               'Harmonic ZPE (1/2 h nu) of the normal mode with the largest N1-H1A stretch character; no anharmonicity.',
               'Electronic energies from the relaxed scan (lower envelope of both chains, 0.10 A grid); equilibrium '
               'IEF-PCM at every point; B3LYP-D3(BJ) tends to underestimate proton-transfer barriers.']}
    atomic_json(SD / 'zpe_barrier.json', out)

    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    xx = np.linspace(a, b, 600)
    ax.plot(xx * BOHR_A, Vx(xx) * KCAL_PER_HA, '-', color='k', lw=1.2, label='1-D potential (interpolated + walls)')
    ax.plot(d, V * KCAL_PER_HA, 'o', mfc='none', color='C0', label='scan, lower envelope')
    for k, lv in enumerate(levels[:4]):
        ax.axhline(lv['E_kcal'], color='C3', lw=0.8, ls='-' if k == 0 else ':')
    ax.axhline(zpe_h, color='C2', lw=1.0, ls='--', label=f'harmonic N-H ZPE {zpe_h:.2f} kcal/mol')
    top_y = max(V.max() * KCAL_PER_HA, levels[min(3, len(levels) - 1)]['E_kcal']) + 2
    ax.set_ylim(-1, top_y); ax.set_xlim(d[0] - 0.25, d[-1] + 0.25)
    ax.set_xlabel('d(N1-H1A) (A)'); ax.grid(alpha=0.3)
    ax.set_ylabel('E - E(N min) (kcal/mol)' if iN is not None else 'E - E(global min) (kcal/mol)')
    ax.legend(fontsize=7); ax.set_title(f"{case['name']}: 1-D proton levels (red) vs profile", fontsize=9)
    fig.tight_layout(); fig.savefig(SD / 'zpe_barrier.png', dpi=150); plt.close(fig)
    print(f"[zpe] N-H normal mode {nh['normal_mode']['freq_cm-1']:.0f} cm-1 (character "
          f"{nh['normal_mode']['nh_stretch_character']:.2f}); local mode {nh['local_mode']['freq_cm-1']:.0f} cm-1, "
          f"k_a {nh['local_mode']['k_a_mdyn_per_A']:.2f} mdyn/A"
          + (f"; scan curvature {curv['k_mdyn_per_A']:.2f} mdyn/A" if curv else '')
          + (f"; freq.json max diff {xcheck:.2f} cm-1" if xcheck is not None else ''))
    print(f'[zpe] {verdict}')
    print(f"[zpe] -> {SD / 'zpe_barrier.json'}")


if __name__ == '__main__':
    main()
