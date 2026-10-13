"""d13: the calculations run for the referees (29 Sep 2026, VM burst 3), re-derived from the raw outputs with our own
parsers and geometry code.  Reads only; writes derived/d13_referee_calcs.json.

(1) M2/M3  Periodic N1-H1A profile.  pw.x outputs of relax_prof{110,125,140,155,170} (d(N1-H1A) of both formula units
    held by CONSTRAINTS 'distance', damped dynamics), relax_i (salt reference), relax_iirfix (the earlier constrained
    configuration) and relax_iir1 (H1A moved in formula unit 1 only, unconstrained BFGS); PBE-D3(BJ)/SG15 and
    PBE0-D3(BJ)/SG15 single points (100/400 Ry) of every profile point and of relax_i.  From each output: the final
    energy, the final crystal coordinates and cell, the last raw force block (which carries the constraint force) and
    the last "Constrained forces" block, the final pressure; the held distances and atoms from the CONSTRAINTS card of
    the input.  Geometry with our own minimum-image code; atom labels from start_i_xray.extxyz (same atom order).
    Slope dE/dd = (F_N1 - F_H1A).u / 2 per formula unit, u the unit vector N1 -> H1A.
(2) M4  Hessian of the neutral (O5-protonated) IEF-PCM minimum (sfa_pair_pcm_d3_oform): our own diagonalisation with
    translations and rotations projected out, IR intensities from the dipole derivatives, the O5-H1A stretch found by
    projection on internal stretches, and our own RRHO thermochemistry (298.15 K, 1 atm) for both minima with the same
    code, so that the differences neutral - salt are like for like.
(3) M6  QTAIM of the nine crystal hydrogen bonds in the 135-atom clusters (sfa_xclu_salt, sfa_xclu_neutral): the raw
    Multiwfn CPprop.txt of each cluster is parsed; nuclei are the (3,-3) points; the bond critical point of each
    H...A contact is assigned by our own geometry (the (3,-1) point nearest to the H...A segment), not by Multiwfn's
    labels; cluster energies from scf.json; basis-function count re-counted from the 6-31+G(d,p) shell structure.
(4) M9  Relaxed C3-C4-S1-N2 scans of sulfanilamide and the sulfonamidate (B3LYP-D3(BJ)/6-311++G(d,p), IEF-PCM): energy
    of every constrained point and of the unconstrained reference from the point JSON files, tau_fold re-measured from
    each point's .xyz, strain = E - E(reference); our own cubic spline through the converged C3-branch points and the
    reference gives the strain at the crystal-ligand torsions, the docked representative poses, the job means and
    every docked seed.
"""
import re

from scipy.interpolate import CubicSpline

from common import *

RY_KJ = RY_EV * EV_KJ                    # kJ/mol per Ry
RY_MEV = RY_EV * 1000.0
RYB_EVA = RY_EV / BOHR_A                 # Ry/Bohr -> eV/A
RUNS = ARES / 'periodic_qe/runs'
LABELS = [l for l in read_extxyz(RES / 'periodic_qe/start_i_xray.extxyz')[2]]
LIX = {l: i for i, l in enumerate(LABELS)}

# ================================================================================================= (1) periodic


def qe_out(path):
    """Final energy, final crystal coordinates, cell (Angstrom rows), last raw and constrained force blocks (Ry/Bohr),
    final pressure, step counts, number of symmetry operations, from one pw.x output."""
    t = Path(path).read_text()
    lines = t.splitlines()
    out = {}
    m = re.findall(r'Final energy\s*=\s*(-?\d+\.\d+)\s*Ry', t)
    out['E_final'] = float(m[-1]) if m else None
    out['E_bang'] = [float(x) for x in re.findall(r'^!\s+total energy\s*=\s*(-?\d+\.\d+)\s*Ry', t, re.M)]
    out['E_bangbang'] = [float(x) for x in re.findall(r'^!!\s+total energy\s*=\s*(-?\d+\.\d+)\s*Ry', t, re.M)]
    out['exx_converged'] = 'EXX self-consistency reached' in t
    alat = float(re.search(r'lattice parameter \(alat\)\s*=\s*([\d.]+)\s*a\.u\.', t).group(1)) * BOHR_A
    ax = []
    for k in (1, 2, 3):
        mm = re.search(rf'a\({k}\) = \(\s*(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*\)', t)
        ax.append([float(mm.group(i)) * alat for i in (1, 2, 3)])
    out['cell'] = np.array(ax)
    # final coordinates: the last ATOMIC_POSITIONS (crystal) block
    starts = [i for i, l in enumerate(lines) if l.startswith('ATOMIC_POSITIONS')]
    if starts:
        s = starts[-1]
        assert 'crystal' in lines[s], lines[s]
        el, fr = [], []
        for l in lines[s + 1:]:
            p = l.split()
            if len(p) < 4 or not re.match(r'^[A-Z][a-z]?$', p[0]):
                break
            el.append(p[0])
            fr.append([float(x) for x in p[1:4]])
        out['elements'] = el
        out['frac'] = np.array(fr)
    for tag, key in (('Forces acting on atoms', 'F_raw'), ('Constrained forces', 'F_con')):
        idx = [i for i, l in enumerate(lines) if tag in l]
        if not idx:
            out[key] = None
            continue
        F = []
        for l in lines[idx[-1] + 1:]:
            mm = re.match(r'\s*atom\s+\d+\s+type\s+\d+\s+force =\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)', l)
            if mm:
                F.append([float(mm.group(i)) for i in (1, 2, 3)])
            elif F:
                break
        out[key] = np.array(F)
    pr = re.findall(r'P=\s*(-?\d+\.\d+)', t)
    out['pressure_kbar'] = float(pr[-1]) if pr else None
    mb = re.search(r'bfgs converged in\s+(\d+) scf cycles and\s+(\d+) bfgs steps', t)
    md = re.search(r'Damped Dynamics: convergence achieved in\s+(\d+) steps', t)
    out['steps'] = int(mb.group(2)) if mb else (int(md.group(1)) if md else None)
    out['scf_cycles'] = int(mb.group(1)) if mb else None
    out['converged'] = bool(mb or md) and out['E_final'] is not None
    ms = re.search(r'(\d+)\s+Sym\. Ops\.', t)
    out['n_symops'] = int(ms.group(1)) if ms else (1 if 'No symmetry found' in t else None)
    return out


def qe_in(path):
    """Crystal coordinates of the input and the CONSTRAINTS card (1-based atoms, target in Bohr)."""
    t = Path(path).read_text()
    lines = t.splitlines()
    s = [i for i, l in enumerate(lines) if l.startswith('ATOMIC_POSITIONS')][0]
    fr = []
    for l in lines[s + 1:]:
        p = l.split()
        if len(p) < 4 or not re.match(r'^[A-Z][a-z]?$', p[0]):
            break
        fr.append([float(x) for x in p[1:4]])
    cons = [(int(a), int(b), float(d)) for a, b, d in re.findall(r"'distance'\s+(\d+)\s+(\d+)\s+(-?[\d.]+)", t)]
    dyn = re.search(r"ion_dynamics\s*=\s*'(\w+)'", t)
    fct = re.search(r'forc_conv_thr\s*=\s*([\d.dDeE+-]+)', t)
    return dict(frac=np.array(fr), constraints=cons, ion_dynamics=dyn.group(1) if dyn else None,
                forc_conv_thr=float(fct.group(1).replace('d', 'e').replace('D', 'e')) if fct else None)


def mic(v, cell):
    f = np.asarray(v, float) @ np.linalg.inv(cell)
    f -= np.round(f)
    cands = [(f + np.array(s) - 1) @ cell for s in np.ndindex(3, 3, 3)]
    return min(cands, key=np.linalg.norm)


def fu_geom(X, cell, suf):
    """Proton metrics of one formula unit (suffix '' or '_2') from Cartesian X (Angstrom, cell rows)."""
    g = lambda l: X[LIX[l + suf]]
    v_nh = mic(g('H1A') - g('N1'), cell)
    v_ho = mic(g('O5') - g('H1A'), cell)
    # the O5 that H1A points to must be the nearest oxygen image of any formula unit
    d_o = {l: np.linalg.norm(mic(X[LIX[l]] - g('H1A'), cell)) for l in LABELS if l.startswith('O')}
    assert min(d_o, key=d_o.get) == 'O5' + suf, (suf, sorted(d_o.items(), key=lambda kv: kv[1])[:2])
    nh, ho = np.linalg.norm(v_nh), np.linalg.norm(v_ho)
    no = np.linalg.norm(v_nh + v_ho)
    po = np.linalg.norm(mic(g('O5') - g('P1'), cell))
    # O4-H4A...O5: the image of O5 (same formula unit label) nearest to H4A, O4 bonded to H4A
    h4 = g('H4A')
    v_h4o5 = mic(g('O5') - h4, cell)
    v_o4h4 = mic(h4 - g('O4'), cell)
    o4o5 = np.linalg.norm(v_o4h4 + v_h4o5)
    ang = np.degrees(np.arccos(np.clip((-v_nh) @ v_ho / nh / ho, -1, 1)))
    return dict(d_N1_H1A=float(nh), d_H1A_O5=float(ho), delta=float(nh - ho), d_N1_O5=float(no), d_P1_O5=float(po),
                d_O4_H4A=float(np.linalg.norm(v_o4h4)), d_H4A_O5=float(np.linalg.norm(v_h4o5)), d_O4_O5=float(o4o5),
                ang_N1_H1A_O5=float(ang), side='N' if nh < ho else 'O', u_NH=(v_nh / nh).tolist())


# the nine hydrogen bonds of the CIF (donor, H, acceptor, symmetry code of the acceptor): a code starting with 1 keeps
# the acceptor in the formula unit of its label, a code starting with 2 puts it in the other one (glide image)
HB9 = [('N1', 'H1A', 'O5', '1_554'), ('N1', 'H1B', 'O6', '2_554'), ('N1', 'H1C', 'O6', '2_654'), ('N2', 'H2A', 'O3', '2_564'),
       ('N2', 'H2B', 'O2', '1_655'), ('O3', 'H3A', 'O1', '1_455'), ('O4', 'H4A', 'O5', '1_455'), ('C3', 'H3', 'O6', '1_655'),
       ('C6', 'H6', 'O4', '1_554')]


def hbonds9(X, cell, suf=''):
    other = '_2' if suf == '' else ''
    out = {}
    for D, Hh, A, sym in HB9:
        a_lab = A + (suf if sym.startswith('1') else other)
        pH, pD = X[LIX[Hh + suf]], X[LIX[D + suf]]
        vHA = mic(X[LIX[a_lab]] - pH, cell)                  # acceptor image nearest to H
        vDH = mic(pH - pD, cell)
        out[f'{D}{Hh}_{A}'] = dict(DA=float(np.linalg.norm(vDH + vHA)), HA=float(np.linalg.norm(vHA)),
                                   DH=float(np.linalg.norm(vDH)))
    return out


def periodic_run(tag):
    o = qe_out(RUNS / tag / f'{tag}.out')
    i = qe_in(RUNS / tag / f'{tag}.in')
    # element check (labels like 'H1A' -> 'H', 'O5_2' -> 'O'; 'S1' -> 'S', 'P1' -> 'P', 'C1' -> 'C', 'N1' -> 'N')
    el_lab = [('H' if l.startswith('H') else l[0]) for l in LABELS]
    assert o['elements'] == el_lab, tag
    X = o['frac'] @ o['cell']
    # heavy-atom RMSD from the X-ray start (mean shift removed: Pc has a floating origin in the a-c plane)
    _, X0x, _, _, _ = read_extxyz(RES / 'periodic_qe/start_i_xray.extxyz')
    hv = [i_ for i_, l in enumerate(LABELS) if not l.startswith('H')]
    dd_ = np.array([mic(X[i_] - X0x[i_], o['cell']) for i_ in hv])
    rmsd_x = float(np.sqrt(np.mean(np.sum((dd_ - dd_.mean(0)) ** 2, axis=1))))
    res = dict(tag=tag, E_Ry=o['E_final'], rmsd_heavy_vs_xray=rmsd_x, converged=o['converged'], steps=o['steps'], scf_cycles=o['scf_cycles'],
               pressure_kbar=o['pressure_kbar'], n_symops=o['n_symops'], ion_dynamics=i['ion_dynamics'],
               forc_conv_thr=i['forc_conv_thr'], volume_A3=float(abs(np.linalg.det(o['cell']))),
               fu={suf or '_1': fu_geom(X, o['cell'], suf) for suf in ('', '_2')},
               hbonds={suf or '_1': hbonds9(X, o['cell'], suf) for suf in ('', '_2')})
    # start geometry from the input (same cell)
    X0 = i['frac'] @ o['cell']
    res['start'] = {suf or '_1': {k: v for k, v in fu_geom(X0, o['cell'], suf).items() if k in ('d_N1_H1A', 'd_H1A_O5', 'side')}
                    for suf in ('', '_2')}
    if i['constraints']:
        res['constraints'] = [dict(atoms=(LABELS[a - 1], LABELS[b - 1]), target_A=d * BOHR_A) for a, b, d in i['constraints']]
        # slope from the raw forces (they carry the constraint force); held distances re-measured
        Fr = o['F_raw']
        sl = {}
        for suf in ('', '_2'):
            u = np.array(res['fu'][suf or '_1']['u_NH'])
            fN, fH = Fr[LIX['N1' + suf]], Fr[LIX['H1A' + suf]]
            s_ryb = float((fN - fH) @ u / 2)
            sl[suf or '_1'] = dict(dEdd_Ry_bohr=s_ryb, dEdd_eV_A=s_ryb * RYB_EVA, dEdd_kJ_A_fu=s_ryb * RY_KJ / BOHR_A,
                                   F_H1A_perp_con_Ry_bohr=float(np.linalg.norm(o['F_con'][LIX['H1A' + suf]]
                                                                               - (o['F_con'][LIX['H1A' + suf]] @ u) * u)))
        res['slope'] = sl
        res['max_constrained_force_Ry_bohr'] = float(np.abs(o['F_con']).max())
        res['max_raw_force_Ry_bohr'] = float(np.abs(Fr).max())
    else:
        res['max_raw_force_Ry_bohr'] = float(np.abs(o['F_raw']).max())
    return res


PROF = ['prof110', 'prof125', 'prof140', 'prof155', 'prof170']
per = {}
ref = periodic_run('relax_i')
per['relax_i'] = ref
for p in PROF:
    per[p] = periodic_run(f'relax_{p}')
per['iirfix'] = periodic_run('relax_iirfix')
per['iir1'] = periodic_run('relax_iir1')
for k, r_ in per.items():
    r_['dE_meV_cell'] = (r_['E_Ry'] - ref['E_Ry']) * RY_MEV
    r_['dE_kJ_fu'] = (r_['E_Ry'] - ref['E_Ry']) * RY_KJ / 2
    r_['dE_kcal_fu'] = r_['dE_kJ_fu'] / 4.184

# hybrid single points (SG15, 100/400 Ry): PBE-D3(BJ) '!' and PBE0-D3(BJ) final '!!'; first '!' of the PBE0 run = PBE
# electrons with the pbe0 D3(BJ) parameters (the D3 part of the correction)
hy = {}
H0 = qe_out(RUNS / 'hybrid_i_sg15_e100/hybrid_i_sg15_e100.out')
P0 = qe_out(RUNS / 'pbe_i_sg15_e100/pbe_i_sg15_e100.out')
for p in PROF:
    H = qe_out(RUNS / f'hybrid_{p}_sg15_e100/hybrid_{p}_sg15_e100.out')
    P = qe_out(RUNS / f'pbe_{p}_sg15_e100/pbe_{p}_sg15_e100.out')
    d_pbe0 = (H['E_bangbang'][-1] - H0['E_bangbang'][-1]) * RY_MEV
    d_pbe = (P['E_bang'][-1] - P0['E_bang'][-1]) * RY_MEV
    d_first = (H['E_bang'][0] - H0['E_bang'][0]) * RY_MEV
    corr = d_pbe0 - d_pbe
    best = per[p]['dE_meV_cell'] + corr
    hy[p] = dict(PBE0_SG15_meV_cell=d_pbe0, PBE_SG15_meV_cell=d_pbe, corr_meV_cell=corr, corr_kJ_fu=corr / 1000 * EV_KJ / 2,
                 corr_D3_meV_cell=d_first - d_pbe, corr_exx_meV_cell=d_pbe0 - d_first,
                 best_meV_cell=best, best_kJ_fu=best / 1000 * EV_KJ / 2, best_kcal_fu=best / 1000 * EV_KJ / 2 / 4.184,
                 SG15_minus_SSSP_meV_cell=d_pbe - per[p]['dE_meV_cell'],
                 exx_converged=H['exx_converged'] and H0['exx_converged'])

# verdicts
dpbe = [per[p]['dE_kJ_fu'] for p in PROF]
dpbe0 = [hy[p]['best_kJ_fu'] for p in PROF]
slopes = [min(per[p]['slope'][s]['dEdd_kJ_A_fu'] for s in ('_1', '_2')) for p in PROF]
held = [max(abs(per[p]['fu'][s]['d_N1_H1A'] - c['target_A']) for s, c in zip(('_1', '_2'), per[p]['constraints'])) for p in PROF]
verdict = dict(monotonic_PBE=all(b > a for a, b in zip([0.0] + dpbe, dpbe)),
               monotonic_PBE0=all(b > a for a, b in zip([0.0] + dpbe0, dpbe0)),
               all_slopes_positive=all(s > 0 for s in slopes), min_slope_kJ_A_fu=min(slopes),
               held_distance_max_dev_A=max(held),
               all_converged=all(per[p]['converged'] for p in PROF),
               max_constrained_force_Ry_bohr=max(per[p]['max_constrained_force_Ry_bohr'] for p in PROF),
               forc_conv_thr=per['prof110']['forc_conv_thr'],
               both_fu_identical=max(abs(per[p]['fu']['_1'][k] - per[p]['fu']['_2'][k]) for p in PROF
                                     for k in ('d_N1_H1A', 'd_H1A_O5', 'd_N1_O5', 'd_O4_O5')),
               proton_side_at=[per[p]['fu']['_1']['side'] for p in PROF],
               first_O_side_point=next(p for p in PROF if per[p]['fu']['_1']['side'] == 'O'),
               min_N1_O5_A=min(per[p]['fu']['_1']['d_N1_O5'] for p in PROF),
               min_N1_O5_at=min(PROF, key=lambda p: per[p]['fu']['_1']['d_N1_O5']),
               iir1_start_fu1_side=per['iir1']['start']['_1']['side'], iir1_start_fu2_side=per['iir1']['start']['_2']['side'],
               iir1_final_sides=[per['iir1']['fu'][s]['side'] for s in ('_1', '_2')],
               iirfix_minus_prof170_kJ_fu=per['iirfix']['dE_kJ_fu'] - per['prof170']['dE_kJ_fu'])

# the earlier constrained configuration (O5-H1A held at 0.98 A) against the profile at the same d(N1-H1A): cubic Hermite
# interpolation between the two bracketing profile points with their energies and slopes (per formula unit)
d_fix = per['iirfix']['fu']['_1']['d_N1_H1A']
lo_, hi_ = [p for p in PROF if per[p]['fu']['_1']['d_N1_H1A'] <= d_fix][-1], [p for p in PROF if per[p]['fu']['_1']['d_N1_H1A'] > d_fix][0]
x0, x1 = per[lo_]['fu']['_1']['d_N1_H1A'], per[hi_]['fu']['_1']['d_N1_H1A']
y0, y1 = per[lo_]['dE_kJ_fu'], per[hi_]['dE_kJ_fu']
m0, m1 = per[lo_]['slope']['_1']['dEdd_kJ_A_fu'], per[hi_]['slope']['_1']['dEdd_kJ_A_fu']
hh = x1 - x0
tt = (d_fix - x0) / hh
E_herm = ((2 * tt ** 3 - 3 * tt ** 2 + 1) * y0 + (tt ** 3 - 2 * tt ** 2 + tt) * hh * m0 + (-2 * tt ** 3 + 3 * tt ** 2) * y1
          + (tt ** 3 - tt ** 2) * hh * m1)
verdict['iirfix_d_NH'] = d_fix
verdict['profile_at_iirfix_d_kJ_fu'] = E_herm
verdict['iirfix_excess_kJ_fu'] = per['iirfix']['dE_kJ_fu'] - E_herm
verdict['iirfix_bracket'] = [lo_, hi_]

# stress estimate for the cell-relaxation limitation: relaxation energy p^2 V / (2B) of each fixed cell, per f.u.,
# difference neutral end point (1.70 A) - salt, for B = 10 and 20 GPa
NA = 6.02214076e23


def relax_E_kJ_fu(p_kbar, V_A3, B_GPa):
    p = p_kbar * 1e8                               # Pa
    return p * p * V_A3 * 1e-30 / (2 * B_GPa * 1e9) * NA / 1000 / 2


stress = {f'B{B}': relax_E_kJ_fu(per['prof170']['pressure_kbar'], per['prof170']['volume_A3'], B)
          - relax_E_kJ_fu(ref['pressure_kbar'], ref['volume_A3'], B) for B in (10, 20)}

# cross-check against profile.json (the pipeline's own summary)
PJ = json.loads((RES / 'periodic_qe/profile.json').read_text())
chk = {}
for p in PROF:
    pj = PJ['points'][p]
    chk[p] = dict(E_Ry=abs(pj['energy_Ry'] - per[p]['E_Ry']),
                  dE_meV=abs(pj['E_minus_i_meV_cell'] - per[p]['dE_meV_cell']),
                  d_NO=abs(pj['formula_units'][0]['d_N1_O5_A'] - per[p]['fu']['_1']['d_N1_O5']),
                  slope_meV_A=abs(pj['slope']['dE_dd_cell_meV_A'] - 2 * per[p]['slope']['_1']['dEdd_eV_A'] * 1000),
                  best_meV=abs(PJ['hybrid']['sg15_e100'][p]['best_estimate_PBE0_D3BJ_meV_cell'] - hy[p]['best_meV_cell']))
periodic = dict(runs=per, hybrid=hy, verdict=verdict, stress_estimate_kJ_fu=stress, check_vs_profile_json=chk,
                max_check=dict((k, max(c[k] for c in chk.values())) for k in chk['prof110']))

# ================================================================================================= (2) M4 Hessian
# average atomic masses, as PySCF's thermo (d05_freq.py found them to reproduce freq.json to 0.1 cm-1)
MASS = {'H': 1.00794, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'P': 30.973762, 'S': 32.065}
AU2CM = 5140.48767


def normal_modes(case):
    Hs = np.load(ARES / case / 'hessian.npy')
    DD = np.load(ARES / case / 'dipole_derivatives.npy')
    syms, X, _ = read_xyz(ARES / case / 'opt.xyz')
    N = len(syms)
    HH = Hs.transpose(0, 2, 1, 3).reshape(3 * N, 3 * N)
    HH = 0.5 * (HH + HH.T)
    m = np.array([MASS[s] for s in syms])
    mw = np.repeat(1 / np.sqrt(m), 3)
    Hm = HH * mw[:, None] * mw[None, :]
    Xb = X / BOHR_A
    sm = np.sqrt(m)
    r = Xb - (m[:, None] * Xb).sum(0) / m.sum()
    vecs = []
    for k in range(3):
        v = np.zeros((N, 3)); v[:, k] = sm; vecs.append(v.ravel())
    for k in range(3):
        e = np.zeros(3); e[k] = 1
        vecs.append((np.cross(e, r) * sm[:, None]).ravel())
    Q, _ = np.linalg.qr(np.array(vecs).T)
    Pp = np.eye(3 * N) - Q @ Q.T
    w, L = np.linalg.eigh(Pp @ Hm @ Pp)
    o = np.argsort(w)
    w, L = w[o][6:], L[:, o][:, 6:]
    nu = np.sign(w) * np.sqrt(np.abs(w)) * AU2CM
    Lc = L * mw[:, None]
    I_km = 974.8801 * (np.einsum('ik,ij->kj', Lc, DD) ** 2).sum(1)
    return syms, X, m, nu, Lc, I_km


def stretch_character(X, Lc, i, j):
    e = (X[j] - X[i]) / np.linalg.norm(X[j] - X[i])
    out = []
    for k in range(Lc.shape[1]):
        d = Lc[:, k].reshape(-1, 3)
        d = d / np.linalg.norm(d)
        out.append(float((d[j] - d[i]) @ e))
    return np.array(out)


def rrho(syms, X, m, nu, E0, T=298.15, p=101325.0):
    kB, h, c, NA_, R = 1.380649e-23, 6.62607015e-34, 2.99792458e10, 6.02214076e23, 8.314462618
    mt = m.sum() * 1.66053906660e-27
    Xc = X - (m[:, None] * X).sum(0) / m.sum()
    I_t = sum(mi * (r_ @ r_ * np.eye(3) - np.outer(r_, r_)) for mi, r_ in zip(m, Xc)) * 1.66053906660e-27 * 1e-20
    th = [h * h / (8 * np.pi ** 2 * I * kB) for I in np.linalg.eigvalsh(I_t)]
    vt = h * c * nu / kB
    zpe = 0.5 * h * c * nu.sum() * NA_
    Ev = R * np.sum(vt * (0.5 + 1 / (np.exp(vt / T) - 1)))
    Hc = Ev + 3 * R * T + R * T
    S = (R * (np.log((2 * np.pi * mt * kB * T / h ** 2) ** 1.5 * kB * T / p) + 2.5)
         + R * (np.log(np.sqrt(np.pi) * np.sqrt(T ** 3 / np.prod(th))) + 1.5)
         + R * np.sum((vt / T) / (np.exp(vt / T) - 1) - np.log(1 - np.exp(-vt / T))))
    J2Ha = 1 / (HA_KJ * 1000)
    return dict(E=E0, ZPE_Ha=zpe * J2Ha, H_Ha=E0 + Hc * J2Ha, G_Ha=E0 + (Hc - T * S) * J2Ha, S_J_molK=S,
                Hcorr_Ha=Hc * J2Ha, Gcorr_Ha=(Hc - T * S) * J2Ha)


m4 = {}
for tag, case in (('salt', 'sfa_pair_pcm_d3'), ('neutral', 'sfa_pair_pcm_d3_oform')):
    syms, X, m, nu, Lc, I_km = normal_modes(case)
    fj = json.loads((ARES / case / 'freq.json').read_text())
    assert [s for s in syms] == [('H' if l.startswith('H') else l[0]) for l in PAIR_LABELS], case
    ch_NH = stretch_character(X, Lc, IDX['N1'], IDX['H1A'])
    ch_OH = stretch_character(X, Lc, IDX['O5'], IDX['H1A'])
    k_nh, k_oh = int(np.argmax(np.abs(ch_NH))), int(np.argmax(np.abs(ch_OH)))
    th = rrho(syms, X, m, nu, fj['E_Ha'])
    m4[tag] = dict(case=case, E_Ha=fj['E_Ha'], n_imag=int((nu < 0).sum()), lowest3=nu[:3].tolist(), n_modes=len(nu),
                   max_abs_nu_vs_freq_json=float(np.abs(nu - np.array(fj['frequencies_cm-1'])).max()),
                   ZPE_freq_json=fj['ZPE_Ha'], H298_freq_json=fj['H298_Ha'], G298_freq_json=fj['G298_Ha'],
                   thermo=th, d_N1_H1A=float(np.linalg.norm(X[IDX['N1']] - X[IDX['H1A']])),
                   d_H1A_O5=float(np.linalg.norm(X[IDX['O5']] - X[IDX['H1A']])),
                   proton_mode=dict(kind='N1-H1A' if tag == 'salt' else 'O5-H1A', n=(k_nh if tag == 'salt' else k_oh) + 1,
                                    nu=float(nu[k_nh if tag == 'salt' else k_oh]), I_km=float(I_km[k_nh if tag == 'salt' else k_oh])),
                   strongest_IR=dict(n=int(np.argmax(I_km)) + 1, nu=float(nu[int(np.argmax(I_km))]), I_km=float(I_km.max())),
                   OH_modes=[dict(n=int(k) + 1, nu=float(nu[k]), I_km=float(I_km[k])) for k in np.argsort(-np.abs(ch_OH))[:1]])
    # the two anion O-H stretches (O3-H3A, O4-H4A)
    m4[tag]['anion_OH'] = [dict(n=int(np.argmax(np.abs(stretch_character(X, Lc, IDX[o], IDX[hh])))) + 1,
                                nu=float(nu[int(np.argmax(np.abs(stretch_character(X, Lc, IDX[o], IDX[hh]))))]))
                           for o, hh in (('O3', 'H3A'), ('O4', 'H4A'))]
s_, n_ = m4['salt']['thermo'], m4['neutral']['thermo']
m4['diff_kcal'] = dict(dE=(n_['E'] - s_['E']) * HA_KCAL, dEZPE=((n_['E'] + n_['ZPE_Ha']) - (s_['E'] + s_['ZPE_Ha'])) * HA_KCAL,
                       dZPE=(n_['ZPE_Ha'] - s_['ZPE_Ha']) * HA_KCAL, dH=(n_['H_Ha'] - s_['H_Ha']) * HA_KCAL,
                       dG=(n_['G_Ha'] - s_['G_Ha']) * HA_KCAL, dS_J_molK=n_['S_J_molK'] - s_['S_J_molK'],
                       dH_freq_json=(m4['neutral']['H298_freq_json'] - m4['salt']['H298_freq_json']) * HA_KCAL,
                       dG_freq_json=(m4['neutral']['G298_freq_json'] - m4['salt']['G298_freq_json']) * HA_KCAL,
                       dEZPE_freq_json=((m4['neutral']['E_Ha'] + m4['neutral']['ZPE_freq_json'])
                                        - (m4['salt']['E_Ha'] + m4['salt']['ZPE_freq_json'])) * HA_KCAL)
m4['neutral_is_scan_B_anchor'] = float(np.abs(read_xyz(ARES / 'sfa_pair_pcm_d3_oform/opt.xyz')[1]
                                              - read_xyz(RES / 'sfa_pair_pcm_d3/scan/B_00_anchor.xyz')[1]).max())
m4['half_hnu_OH_kcal'] = 0.5 * m4['neutral']['proton_mode']['nu'] * CM1_KCAL

# ================================================================================================= (3) M6 cluster QTAIM


def cpprop(path):
    """Every critical point of a Multiwfn CPprop.txt: index, type, position (A), nucleus (for (3,-3)), rho, lap, G, V, H,
    ellipticity."""
    t = Path(path).read_text()
    out = []
    for blk in re.split(r'-{10,}\s+CP\s+', t)[1:]:
        mi = re.match(r'(\d+),\s+Type\s+\((\d),([+-]\d)\)', blk)
        g = lambda pat: float(re.search(pat, blk).group(1))
        d = dict(cp=int(mi.group(1)), type=f'({mi.group(2)},{mi.group(3)})',
                 pos=[float(x) for x in re.search(r'Position \(Angstrom\):\s+(\S+)\s+(\S+)\s+(\S+)', blk).groups()],
                 rho=g(r'Density of all electrons:\s+(\S+)'), lap=g(r'Laplacian of electron density:\s+(\S+)'),
                 G=g(r'Lagrangian kinetic energy G\(r\):\s+(\S+)'), V=g(r'Potential energy density V\(r\):\s+(\S+)'),
                 H=g(r'Energy density E\(r\) or H\(r\):\s+(\S+)'))
        me = re.search(r'Ellipticity of electron density:\s+(\S+)', blk)
        d['ellipticity'] = float(me.group(1)) if me else None
        mn = re.search(r'Corresponding nucleus:\s+(\d+)\((\w+)\s*\)', blk)
        if mn:
            d['nucleus'] = (int(mn.group(1)), mn.group(2))
        out.append(d)
    return out


def seg_dist(p, a, b):
    ab = b - a
    s = float(np.clip((p - a) @ ab / (ab @ ab), 0, 1))
    return float(np.linalg.norm(p - (a + s * ab))), s


BUILD = json.loads((ARES / 'crystal_cluster_qtaim/build.json').read_text())
QN = json.loads((ARES / 'crystal_cluster_qtaim/qtaim_nine_bonds.json').read_text())
clu = {}
for st in ('salt', 'neutral'):
    case = QN['cases'][st]
    cps = cpprop(ARES / case / 'mwfn_extra/CPprop.txt')
    nuc = {c['nucleus'][0] - 1: (c['nucleus'][1], np.array(c['pos'])) for c in cps if c['type'] == '(3,-3)'}
    bcps = [c for c in cps if c['type'] == '(3,-1)']
    counts = {t: sum(1 for c in cps if c['type'] == t) for t in ('(3,-3)', '(3,-1)', '(3,+1)', '(3,+3)')}
    inst = {}
    for h in BUILD['structures'][st]['hbond_instances']:
        D, Hh, A = h['atoms']
        assert nuc[Hh][0] == 'H', h['name']
        pD, pH, pA = nuc[D][1], nuc[Hh][1], nuc[A][1]
        dHA, dDH = float(np.linalg.norm(pA - pH)), float(np.linalg.norm(pH - pD))
        # the contact: H...A, or N1...H1A where the proton has moved to the acceptor
        rev = dDH > dHA
        a, b = (pD, pH) if rev else (pH, pA)
        best = min(bcps, key=lambda c: seg_dist(np.array(c['pos']), a, b)[0])
        dd, s = seg_dist(np.array(best['pos']), a, b)
        # heavy nuclei coincide with their (3,-3) points; the (3,-3) point of H sits a few hundredths of an angstrom
        # off the proton, so only D...A is compared with the geometry recorded in build.json
        inst[h['name']] = dict(type=h['type'], primary=h['primary'], reverse=rev,
                               D_A=float(np.linalg.norm(pA - pD)), DA_vs_build=abs(float(np.linalg.norm(pA - pD)) - h['D..A']),
                               cp=best['cp'], cp_off_segment_A=dd, cp_frac=s,
                               **{k: best[k] for k in ('rho', 'lap', 'G', 'V', 'H', 'ellipticity')})
    scf = json.loads((ARES / case / 'scf.json').read_text())
    elements = [nuc[k][0] for k in sorted(nuc)]
    nao = sum({'H': 5, 'C': 18, 'N': 18, 'O': 18, 'P': 22, 'S': 22}[e] for e in elements)   # 6-31+G(d,p), spherical d
    clu[st] = dict(case=case, n_atoms=len(nuc), counts=counts, poincare_hopf=counts['(3,-3)'] - counts['(3,-1)'] + counts['(3,+1)'] - counts['(3,+3)'],
                   instances=inst, E_Ha=scf['e_tot_Ha'], gap_eV=scf['gap_eV'], nao_counted=nao, nao_scf=scf['nao'],
                   n_electrons=scf['nelectron'], converged=scf['converged'], xyz_sha256=scf['xyz_sha256'],
                   basis=scf['level']['basis_print'], grid=scf['level']['grid_level'])
bonds = []
for b in QN['bonds']:
    s_i, n_i = clu['salt']['instances'][b['instance']], clu['neutral']['instances'][b['instance']]
    same_type = lambda st: [v['rho'] for v in clu[st]['instances'].values() if v['type'] == b['bond']]
    bonds.append(dict(bond=b['bond'], instance=b['instance'], salt=s_i, neutral=n_i,
                      d_rho=n_i['rho'] - s_i['rho'], d_rho_pct=100 * (n_i['rho'] - s_i['rho']) / s_i['rho'],
                      d_H=n_i['H'] - s_i['H'], d_lap=n_i['lap'] - s_i['lap'],
                      rho_range_salt=[min(same_type('salt')), max(same_type('salt'))],
                      rho_range_neutral=[min(same_type('neutral')), max(same_type('neutral'))],
                      n_instances=len(same_type('salt')),
                      check_rho_vs_json=max(abs(s_i['rho'] - b['salt']['bcp']['rho']), abs(n_i['rho'] - b['neutral']['bcp']['rho']))))
npairs = BUILD['structures']['salt']['n_pairs']
m6 = dict(bonds=bonds, clusters={k: {kk: vv for kk, vv in v.items() if kk != 'instances'} for k, v in clu.items()},
          n_pairs=npairs, dE_kJ_fu=(clu['neutral']['E_Ha'] - clu['salt']['E_Ha']) / npairs * HA_KJ,
          dE_kcal_fu=(clu['neutral']['E_Ha'] - clu['salt']['E_Ha']) / npairs * HA_KCAL,
          neutral_source=BUILD['structures']['neutral']['source'], salt_source=BUILD['structures']['salt']['source'],
          max_atom_shift=BUILD['structures']['neutral']['max_atom_shift_vs_salt_A'],
          dE_json_kJ=QN['cluster_energy']['neutral']['dE_per_fu_kJ'])

# ================================================================================================= (4) M9 torsion scans
F_ = json.loads((REPO / 'paper/main/data/facts.json').read_text())
PS = json.loads((DOCK / 'results/production_summary.json').read_text())
TS = json.loads((RES / 'torsion_scan/torsion_strain.json').read_text())


def tau_fold_xyz(path):
    syms, X, _ = read_xyz(path)
    # C3, C4, S1, N2 by position in the case xyz (constraint atoms 3, 4, 11, 8, 1-based), C5 = the other ring C on C4
    c3, c4, s1, n2 = 2, 3, 10, 7
    assert syms[c3] == 'C' and syms[c4] == 'C' and syms[s1] == 'S' and syms[n2] == 'N'
    ring = [i for i, s in enumerate(syms) if s == 'C' and i != c4 and np.linalg.norm(X[i] - X[c4]) < 1.6]
    c5 = [i for i in ring if i != c3][0]
    t3, t5 = dihedral(X[c3], X[c4], X[s1], X[n2]), dihedral(X[c5], X[c4], X[s1], X[n2])
    return max(abs(t3), abs(t5)), t3, t5, angle(X[c4], X[s1], X[n2])


m9 = {}
for sp, case in (('sfa_neutral', 'sulfanilamide_pcm_d3'), ('sfa_anion', 'sulfonamidate_pcm_d3')):
    d = RES / case / 'torsion_scan'
    rf = json.loads((d / 'ref.json').read_text())
    tref = tau_fold_xyz(d / 'ref.xyz')
    pts = []
    for p in sorted(d.glob('tau*.json')):
        j = json.loads(p.read_text())
        if j.get('kind') != 'constrained':
            continue
        tf, t3, t5, csn = tau_fold_xyz(p.with_suffix('.xyz'))
        pts.append(dict(file=p.name, target=j['target_tau_fold_deg'], branch=j['branch'], converged=j['converged'],
                        E_Ha=j['E_Ha'], strain=(j['E_Ha'] - rf['E_Ha']) * HA_KCAL, tau_fold=tf, tau_C3=t3, tau_C5=t5,
                        ang_CSN=csn, dev_target=abs(tf - j['target_tau_fold_deg']), labels=j.get('labels', [])))
    c3 = sorted([p for p in pts if p['branch'] == 'C3' and p['converged']], key=lambda p: p['tau_fold'])
    xs = [tref[0]] + [p['tau_fold'] for p in c3]
    ys = [0.0] + [p['strain'] for p in c3]
    o = np.argsort(xs)
    cs = CubicSpline(np.array(xs)[o], np.array(ys)[o])
    grid = [p for p in c3 if abs(p['target'] - round(p['target'] / 10) * 10) < 1e-6]
    evalf = lambda t: float(cs(t))
    tab = {}
    for r_ in TS['table']:
        if r_['case'] != case or r_['kind'] == 'grid':
            continue
        tab[r_['label']] = dict(kind=r_['kind'], tau=r_['tau_fold_deg'], strain_src=r_['strain_kcal'], how=r_['how'],
                                strain_spline=evalf(r_['tau_fold_deg']), mmff=r_['mmff_relaxed_kcal'])
    # docked seeds of every production job of this species: strain range over the seeds (own spline)
    jobs = {}
    for key, jb in PS['jobs'].items():
        if jb['ligand'] != sp or jb.get('diagnostic') or not key.endswith('__prod'):
            continue
        taus = [q['tau_fold'] for q in jb['torsion']['per_seed']]
        st = [evalf(t) for t in taus]
        jobs[key] = dict(n=len(taus), tau_min=min(taus), tau_max=max(taus), tau_mean=float(np.mean(taus)),
                         strain_min=min(st), strain_max=max(st), strain_at_mean=evalf(float(np.mean(taus))),
                         mmff_min=min(q['strain_mmff_relaxed'] for q in jb['torsion']['per_seed']),
                         mmff_max=max(q['strain_mmff_relaxed'] for q in jb['torsion']['per_seed']))
    m9[sp] = dict(case=case, E_ref=rf['E_Ha'], ref_converged=rf['converged'], tau_ref=tref[0], ang_CSN_ref=tref[3],
                  n_points=len(pts), n_converged=sum(p['converged'] for p in pts), points=pts,
                  grid=[dict(tau=p['target'], strain=p['strain']) for p in grid],
                  strain_min=min(p['strain'] for p in c3), strain_max=max(p['strain'] for p in c3),
                  tau_at_min=min(c3, key=lambda p: p['strain'])['tau_fold'],
                  tau_at_max=max(c3, key=lambda p: p['strain'])['tau_fold'],
                  max_dev_target=max(p['dev_target'] for p in pts),
                  mirror=[dict(tau=p['target'], C5=p['strain'], C3=next(q['strain'] for q in c3 if abs(q['target'] - p['target']) < 1e-6))
                          for p in pts if p['branch'] == 'C5'],
                  labelled=tab, jobs=jobs,
                  max_spline_vs_src=max(abs(v['strain_spline'] - v['strain_src']) for v in tab.values()),
                  src_ref_E=TS['per_case'][case]['reference']['E_Ha'],
                  monotonic_grid=all(b['strain'] > a['strain'] for a, b in zip(grid, grid[1:])))

dump('d13_referee_calcs.json', dict(periodic=periodic, m4=m4, m6=m6, m9=m9))
v = periodic['verdict']
print('profile PBE kJ/fu', [round(x, 3) for x in dpbe], 'PBE0', [round(x, 3) for x in dpbe0])
print('slopes kJ/A/fu', [round(x, 2) for x in slopes], 'verdict', {k: v[k] for k in ('monotonic_PBE', 'monotonic_PBE0', 'all_slopes_positive')})
print('checks vs profile.json', periodic['max_check'])
print('iir1', per['iir1']['dE_kJ_fu'], v['iir1_start_fu1_side'], v['iir1_start_fu2_side'], v['iir1_final_sides'])
print('M4 diff', {k: round(x, 3) for k, x in m4['diff_kcal'].items()}, 'OH', m4['neutral']['proton_mode'])
print('M6 dE kJ/fu', m6['dE_kJ_fu'], 'check', max(b['check_rho_vs_json'] for b in bonds))
for sp in m9:
    print('M9', sp, m9[sp]['n_converged'], m9[sp]['strain_min'], m9[sp]['strain_max'], 'spline vs src', m9[sp]['max_spline_vs_src'])
