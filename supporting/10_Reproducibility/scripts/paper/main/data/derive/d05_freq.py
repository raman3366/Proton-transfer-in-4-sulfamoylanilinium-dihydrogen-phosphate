"""d05: harmonic frequencies, IR intensities, key-band assignments and thermochemistry of the IEF-PCM ion pair.

Our own re-diagonalisation of the archived finite-difference Hessian (hessian.npy, Ha/Bohr^2, (27,27,3,3)) with
translations/rotations projected out, IR intensities from the archived dipole derivatives (dipole_derivatives.npy,
d mu / d x in a.u.), and an assignment of each mode from bond-stretch and angle-bend projections of its Cartesian
displacement (B-matrix rows of the internal coordinates) -- a heuristic character, not a full PED.
Thermochemistry: own rigid-rotor / harmonic-oscillator / ideal-gas sums at 298.15 K and 1 atm.
Scale factor: 0.9688 for fundamentals, B3-LYP/6-311+G(d,p), Merrick, Moran & Radom 2007 (lit key
merrick2007evaluation).  Checked 29 Sep 2026: Crossref record and abstract via lit.py; the full text is closed, so the
table row was read from two independent tabulations that agree (the liyuanhe211 compilation and Table 1 of
doi:10.1016/j.newar.2017.01.001).  6-311++G(d,p) is not in Merrick's table; 6-311+G(d,p) is the closest basis.
"""
import itertools
from common import *

H = np.load(ARES / 'sfa_pair_pcm_d3/hessian.npy')          # (N, N, 3, 3)
DD = np.load(ARES / 'sfa_pair_pcm_d3/dipole_derivatives.npy')   # (3N, 3)
syms, X, _ = read_xyz(ARES / 'sfa_pair_pcm_d3/opt.xyz')
fr = json.loads((RES / 'sfa_pair_pcm_d3/freq.json').read_text())
N = len(syms)
MASS_MAIN = {'H': 1.00782503207, 'C': 12.0, 'N': 14.0030740048, 'O': 15.99491461956, 'P': 30.97376163,
             'S': 31.97207100}
MASS_AVG = {'H': 1.00794, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'P': 30.973762, 'S': 32.065}
AU2CM = 5140.48767                                   # sqrt(Ha / (Bohr^2 amu)) -> cm-1
HH = H.transpose(0, 2, 1, 3).reshape(3 * N, 3 * N)
HH = 0.5 * (HH + HH.T)
Xb = X / BOHR_A


def tr_projector(m):
    """Projector onto the complement of translations and rotations in mass-weighted coordinates."""
    sm = np.sqrt(m)
    com = (m[:, None] * Xb).sum(0) / m.sum()
    r = Xb - com
    vecs = []
    for k in range(3):
        v = np.zeros((N, 3)); v[:, k] = sm; vecs.append(v.ravel())
    for k in range(3):
        e = np.zeros(3); e[k] = 1
        v = np.cross(e, r) * sm[:, None]; vecs.append(v.ravel())
    Q, _ = np.linalg.qr(np.array(vecs).T)
    return np.eye(3 * N) - Q @ Q.T


def analyse(massmap, project=True):
    m = np.array([massmap[s] for s in syms])
    mw = np.repeat(1 / np.sqrt(m), 3)
    Hm = HH * mw[:, None] * mw[None, :]
    if project:
        Pp = tr_projector(m)
        Hm = Pp @ Hm @ Pp
    w, L = np.linalg.eigh(Hm)
    order = np.argsort(w)
    w, L = w[order], L[:, order]
    nu = np.sign(w) * np.sqrt(np.abs(w)) * AU2CM
    return m, nu, L


res = {}
for tag, mm in [('main_isotope', MASS_MAIN), ('average', MASS_AVG)]:
    m, nu, L = analyse(mm)
    vib = nu[6:]
    ref = np.array(fr['frequencies_cm-1'])
    res[tag] = dict(max_abs_diff_vs_freq_json=float(np.abs(np.sort(vib) - np.sort(ref)).max()),
                    lowest=float(vib[0]), highest=float(vib[-1]), n_imag=int((nu[6:] < 0).sum()))
best = min(res, key=lambda k: res[k]['max_abs_diff_vs_freq_json'])
m, nu, L = analyse(MASS_MAIN if best == 'main_isotope' else MASS_AVG)
nu, L = nu[6:], L[:, 6:]
mw = np.repeat(1 / np.sqrt(m), 3)
# IR intensities: dmu/dQ = sum_i dmu/dx_i * L_i / sqrt(m_i)  (x in Bohr, Q in amu^1/2 Bohr), km/mol factor 974.8801
dmudQ = (DD * (L * mw[:, None])[:, :, None].transpose(1, 0, 2)[0][:, None] if False else None)
dmudQ = np.einsum('ik,ij->kj', L * mw[:, None], DD)             # (modes, 3)
I_km = 974.8801 * (dmudQ ** 2).sum(1)
ref_I = np.array([x['intensity_km_mol'] for x in fr['ir']])
ref_nu = np.array([x['freq_cm-1'] for x in fr['ir']])

# ---- internal coordinates for the assignment
LAB = PAIR_LABELS
RC = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'P': 1.07, 'S': 1.05}
bonds = [(i, j) for i, j in itertools.combinations(range(N), 2) if dist(X[i], X[j]) < 1.2 * (RC[syms[i]] + RC[syms[j]])]
nbr = {i: [] for i in range(N)}
for i, j in bonds:
    nbr[i].append(j); nbr[j].append(i)
bends = [(a, c, b) for c in range(N) for a, b in itertools.combinations(sorted(nbr[c]), 2)]


def bstretch(i, j):
    e = (X[j] - X[i]) / dist(X[i], X[j])
    row = np.zeros(3 * N); row[3 * i:3 * i + 3] = -e; row[3 * j:3 * j + 3] = e
    return row


def bbend(a, c, b):
    """Wilson B row of the angle a-c-b (radians per Angstrom)."""
    u, v = X[a] - X[c], X[b] - X[c]
    lu, lv = np.linalg.norm(u), np.linalg.norm(v)
    u, v = u / lu, v / lv
    cos = np.clip(u @ v, -1, 1); sin = np.sqrt(1 - cos * cos)
    ga = (cos * u - v) / (lu * sin)
    gb = (cos * v - u) / (lv * sin)
    row = np.zeros(3 * N); row[3 * a:3 * a + 3] = ga; row[3 * b:3 * b + 3] = gb; row[3 * c:3 * c + 3] = -(ga + gb)
    return row


coords = [(f'str {LAB[i]}-{LAB[j]}', bstretch(i, j), 1.0) for i, j in bonds]
# the H-bond stretch H1A...O5 and the N1...O5 distance
coords.append(('str H1A...O5', bstretch(IDX['H1A'], IDX['O5']), 1.0))
coords += [(f'bend {LAB[a]}-{LAB[c]}-{LAB[b]}', bbend(a, c, b), 1.0) for a, c, b in bends]
B = np.array([c[1] for c in coords])
names = [c[0] for c in coords]
# a length scale for bends so that stretches and bends are comparable: 1 rad ~ 1 A (common crude choice)
modes = []
for k in range(len(nu)):
    d = (L[:, k] * mw)                                      # Cartesian displacement per unit Q
    d = d / np.linalg.norm(d)
    q = B @ d
    w = q ** 2
    w = w / w.sum()
    top = np.argsort(-w)[:4]
    ke = (L[:, k] ** 2).reshape(N, 3).sum(1)                # kinetic-energy fraction per atom
    atoms = [(LAB[i], float(ke[i])) for i in np.argsort(-ke)[:3]]
    modes.append(dict(n=k + 1, nu=float(nu[k]), I_km=float(I_km[k]),
                      top_internal=[(names[i], float(w[i])) for i in top],
                      top_atoms=atoms,
                      w_str_N1H1A=float(w[names.index('str N1-H1A')]),
                      phase_S1O=float(np.sign(q[names.index('str O1-S1')] * q[names.index('str O2-S1')])),
                      phase_N2H=float(np.sign(q[names.index('str N2-H2A')] * q[names.index('str N2-H2B')]))))

# ---- thermochemistry (RRHO, ideal gas, 298.15 K, 101325 Pa), own implementation
kB, h, c, NA, R = 1.380649e-23, 6.62607015e-34, 2.99792458e10, 6.02214076e23, 8.314462618
T, p = 298.15, 101325.0
mtot = m.sum() * 1.66053906660e-27
Xc = X - (m[:, None] * X).sum(0) / m.sum()
I_t = sum(mi * (np.dot(r, r) * np.eye(3) - np.outer(r, r)) for mi, r in zip(m, Xc)) * 1.66053906660e-27 * 1e-20
Ip = np.linalg.eigvalsh(I_t)
th = np.array([h * h / (8 * np.pi ** 2 * I * kB) for I in Ip])
vt = h * c * nu / kB                                         # vibrational temperatures (K), real modes only
ZPE_J = 0.5 * h * c * nu.sum() * NA
E_vib = R * np.sum(vt * (0.5 + 1 / (np.exp(vt / T) - 1)))
E_tr, E_rot = 1.5 * R * T, 1.5 * R * T
H_corr = E_vib + E_tr + E_rot + R * T                        # J/mol
S_tr = R * (np.log((2 * np.pi * mtot * kB * T / h ** 2) ** 1.5 * kB * T / p) + 2.5)
sigma = 1
S_rot = R * (np.log(np.sqrt(np.pi) / sigma * np.sqrt(T ** 3 / np.prod(th))) + 1.5)
x = vt / T
S_vib = R * np.sum(x / (np.exp(x) - 1) - np.log(1 - np.exp(-x)))
S = S_tr + S_rot + S_vib
G_corr = H_corr - T * S
E0 = fr['E_Ha']
thermo = dict(T=T, p_Pa=p, ZPE_Ha=ZPE_J / (HA_KJ * 1000), ZPE_kcal=ZPE_J / 4184, ZPE_eV=ZPE_J / (EV_KJ * 1000),
              H_corr_Ha=H_corr / (HA_KJ * 1000), G_corr_Ha=G_corr / (HA_KJ * 1000),
              H298_Ha=E0 + H_corr / (HA_KJ * 1000), G298_Ha=E0 + G_corr / (HA_KJ * 1000),
              S_J_molK=S, S_tr=S_tr, S_rot=S_rot, S_vib=S_vib,
              freq_json=dict(ZPE_Ha=fr['ZPE_Ha'], H298_Ha=fr['H298_Ha'], G298_Ha=fr['G298_Ha'], E_Ha=E0),
              note='harmonic, unscaled, all 75 modes incl. the 20-52 cm-1 librations (no quasi-RRHO); PCM solute treated as an ideal-gas molecule')
out = dict(source_hessian=str(ARES / 'sfa_pair_pcm_d3/hessian.npy'),
           source_dipole_derivatives=str(ARES / 'sfa_pair_pcm_d3/dipole_derivatives.npy'),
           mass_choice=best, mass_tests=res, n_modes=len(nu), n_imag=int((nu < 0).sum()),
           max_abs_I_diff_vs_freq_json=float(np.abs(I_km - ref_I).max()),
           max_abs_nu_diff_vs_freq_json=float(np.abs(nu - ref_nu).max()),
           scale_factor=dict(value=0.9688, level='B3-LYP/6-311+G(d,p), fundamentals', lit_key='merrick2007evaluation',
                             provenance='Merrick 2007 row B3-LYP/6-311+G(d,p) (fundamentals), read from two independent tabulations that agree: '
                                        'github.com/liyuanhe211/Collection_of_Frequency_Scale_Factors and Table 1 of '
                                        'doi:10.1016/j.newar.2017.01.001; full text closed (Crossref record and abstract checked '
                                        'with lit.py on 29 Sep 2026); 6-311++G(d,p) is not in Merrick\'s table, 6-311+G(d,p) is '
                                        'the closest basis'),
           modes=modes, thermo=thermo)
dump('d05_freq.json', out)
