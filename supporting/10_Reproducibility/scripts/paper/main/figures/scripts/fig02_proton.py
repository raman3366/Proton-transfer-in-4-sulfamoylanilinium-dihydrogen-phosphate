"""Where the proton sits: X-ray against optimised ion pairs, the relaxed proton scan, and the functional check
(DOUBLE width).

    cd /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate
    nice -n 10 env OMP_NUM_THREADS=2 /root/workspace/cqc-dft/tools/paper-env/bin/python \
        paper/main/figures/scripts/fig02_proton.py

(a)-(c) Stick overlays of the X-ray ion pair (black) with (a) the IEF-PCM (water) optimum, (b) the gas-phase optimum
        (H1A on O5) and (c) formula unit 1 of the relaxed periodic PBE-D3(BJ) cell.  Heavy atoms plus H1A and H6, the
        two hydrogens of the R2,2(8) ring; dashed: N1-H1A...O5 (O5-H1A...N1 in the gas phase) and C6-H6...O4.
        (a), (b): superposed on the cation heavy atoms (the numbers are the facts' separate cation and anion fits);
        (c): cell frame, no rotation, the mean heavy-atom shift of the whole cell removed (fact periodic.rmsd_heavy).
        One common view (normal of the benzene / R2,2(8) planes of the X-ray pair, S1 -> N1 horizontal) and one scale.
(d)     Relaxed scan of d(N1-H1A) in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p): chain A from the salt minimum,
        chain B from the neutral minimum; filled = converged, open = stopped at the step limit with its forces below
        the thresholds (A_02, the barrier top; SI Table S7);
        grey curve: cubic spline through the lower envelope (its maximum reproduces fact scan.barrier_spline);
        dashed: harmonic zero-point level 1/2 h nu of the N1-H1A stretch above the salt minimum (no ratio to the barrier).
(e)     Single points of four functionals and DF-MP2 at the B3LYP-D3(BJ) scan geometries: the d = 1.25 A point and the
        neutral end point, relative to the salt end point (sfa_pair_pcm_d3/fsp_scan.json, burst 2; facts section 7).
(f)     Constrained N1-H1A profile of the periodic cell (VM burst 3, facts section 32): PBE-D3(BJ) energies (with short
        strokes whose slope is dE/dd from the constraint force) and PBE0-D3(BJ) single points, per formula unit, relative
        to the relaxed salt (i); grey open markers: the earlier constrained single configuration (relax_iirfix, O5-H1A
        held at 0.98 A; facts sections 28 and 36), no longer the measure of the cost.  Shaded: H1A nearer O5 (delta > 0).
(g)     Donor-acceptor distances of the crystal along the profile: N1...O5 (the bond kept in the ion pair) and the three
        crystal-only bonds that respond (N1-H1B...O6, N1-H1C...O6, O4-H4A...O5), facts periodic.profile.*.hb.*.
Updated 29 Sep 2026 (VM burst 3): the "+52.8 kJ/mol" line under (c) replaced by panels (f) and (g).
Updated 29 Sep 2026 (prose audit, burst 3 round 2; re-review O9): the ZPE label of (d) no longer states the ratio to the
barrier ("4.3 × the barrier"); the abstract and text rest on the one-dimensional proton model, not on that ratio.
The withdrawn relax_iifix is not used.
Updated 11 Oct 2026 (last reading before submission): two labels moved clear of lines, "neutral pair" in (d), which the
grey zero line crossed, and "O5-H1A held at 0.98 A" in (f), whose second line lay on the top frame; no number changed.
All numbers from paper/main/data/facts.json; geometries re-measured and asserted against the facts.
Updated 29 Sep 2026 (burst 2): DF-MP2 added to (e) from the rerun of the functional check.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

import _saph as S  # noqa: E402
from figkit import style, geometry as geo  # noqa: E402
from figkit.draw2d import draw_overlay  # noqa: E402

STEM = 'fig02_proton'
W_MM = 175.0
H_MM = 185.0
PROF = ['p110', 'p125', 'p140', 'p155', 'p170']
XL = (1.03, 1.74)                 # common d(N1-H1A) range of the molecular scan (d) and the crystal profile (f, g)
C_PBE, C_PBE0 = '#009E73', '#6A3D9A'   # validated pair (dataviz validate_palette.js: CVD dE 22.4, normal 30.4, contrast >= 3:1)
HA_KCAL = 627.5095            # the conversion used by the facts derivations (derive/common.py)

SUB_LABELS = [l for l in S.PAIR_LABELS if not l.startswith('H')] + ['H1A', 'H6']
SUB = [S.IDX[l] for l in SUB_LABELS]
SUB_SYM = [S.SYMBOLS[i] for i in SUB]
SI = {l: k for k, l in enumerate(SUB_LABELS)}


# ------------------------------------------------------------------------------------------------ overlays
def view_frame(X):
    """Rotation R and centre c: z along the mean normal of the benzene ring and the R2,2(8) ring of the X-ray pair,
    x along S1 -> N1 projected on the view plane (cation left, anion right)."""
    def normal(labels):
        P = X[[S.IDX[l] for l in labels]]
        P = P - P.mean(0)
        return np.linalg.svd(P)[2][2]
    n1 = normal(['C1', 'C2', 'C3', 'C4', 'C5', 'C6'])
    n2 = normal(['H1A', 'N1', 'C1', 'C6', 'H6', 'O4', 'P1', 'O5'])
    if n1 @ n2 < 0:
        n2 = -n2
    z = n1 + n2
    z /= np.linalg.norm(z)
    v = X[S.IDX['N1']] - X[S.IDX['S1']]
    x = v - (v @ z) * z
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    R = np.vstack([x, y, z])
    c = X[[S.IDX[l] for l in S.PAIR_LABELS if not l.startswith('H')]].mean(0)
    return R, c


def overlay_geoms():
    X = S.xray_pair()
    cat_heavy = [S.IDX[l] for l in S.CATION if not l.startswith('H')]
    out = {}
    for m, case in (('pcm', 'sfa_pair_pcm_d3'), ('gas', 'sfa_pair_gas_d3')):
        Y = S.case_pair(case)
        Rk, t, rmsd = geo.kabsch(Y[cat_heavy], X[cat_heavy])
        out[m] = Y @ Rk.T + t
        assert abs(rmsd - S.F(f'{m}.geom.kabsch_cation')) < 1e-6, (m, rmsd)
        an = [S.IDX[l] for l in S.ANION if not l.startswith('H')]
        assert abs(geo.kabsch(Y[an], X[an])[2] - S.F(f'{m}.geom.kabsch_anion')) < 1e-6
    P, R0, rmsd = S.cell_shift_removed('relax_i')
    assert np.abs(R0 - X).max() < 1e-6
    assert abs(rmsd - S.F('periodic.rmsd_heavy')) < 5e-4, rmsd
    out['periodic'] = P
    return X, out


def draw_overlays(fig, boxes):
    """boxes: model -> (x, y, w, h) mm of the drawing area. Common view and scale."""
    X, G = overlay_geoms()
    R, c = view_frame(X)
    V = {'xray': (X - c) @ R.T}
    for m, Y in G.items():
        V[m] = (Y - c) @ R.T
    # common scale: largest extent over all three overlays (subset atoms) must fit every box
    allxy = np.vstack([V[m][SUB][:, :2] for m in V])
    lo, hi = allxy.min(0) - 0.55, allxy.max(0) + 0.55
    w_mm, h_mm = boxes['pcm'][2], boxes['pcm'][3]
    scale = min(w_mm / (hi[0] - lo[0]), h_mm / (hi[1] - lo[1]))        # mm per Angstrom
    centre = (lo + hi) / 2
    axes = {}
    hb_std = [(SI['N1'], SI['H1A'], SI['O5']), (SI['C6'], SI['H6'], SI['O4'])]
    hb_gas = [(SI['O5'], SI['H1A'], SI['N1']), (SI['C6'], SI['H6'], SI['O4'])]
    for m, (x, y, w, h) in boxes.items():
        ax = S.axes_mm(fig, x, y, w, h)
        ax.set_xlim(centre[0] - w / scale / 2, centre[0] + w / scale / 2)
        ax.set_ylim(centre[1] - h / scale / 2, centre[1] + h / scale / 2)
        ax.set_aspect('equal', adjustable='box')
        ax.set_autoscale_on(False)
        geoms = {'xray': V['xray'][SUB], m: V[m][SUB]}
        bonds = {k: [(i, j) for i, j in geo.detect_bonds(SUB_SYM, g)] for k, g in geoms.items()}
        # proton bookkeeping: H1A covalently bound to N1 except in the gas phase
        for k, g in geoms.items():
            nb = [j for i, j in bonds[k] if i == SI['H1A']] + [i for i, j in bonds[k] if j == SI['H1A']]
            assert nb == [SI['O5' if k == 'gas' else 'N1']], (k, nb)
        hbs = {'xray': hb_std, m: hb_gas if m == 'gas' else hb_std}
        # model at the bottom (thick), X-ray on top (thin black), so both stay visible where they coincide
        order = {m: geoms[m], 'xray': geoms['xray']}
        d = draw_overlay(ax, SUB_SYM, order, {'xray': S.MODEL_COLOURS['xray'], m: S.MODEL_COLOURS[m]},
                         bonds=bonds, include_h=True, lw={'xray': 0.8, m: 2.2}, labels=SUB_LABELS,
                         label_atoms=['N1', 'H1A', 'O5', 'O4', 'H6', 'P1', 'S1', 'N2'], legend=False, fit=False,
                         hbonds=hbs, hbond_obstacles=True, strict=True, fontsize=style.SMALL_PT, halo=False)
        # the proton: larger discs on H1A (model: colour, X-ray: black)
        from matplotlib.patches import Circle
        for k, r, fc, z in ((m, 0.17, S.MODEL_COLOURS[m], 30), ('xray', 0.11, S.MODEL_COLOURS['xray'], 31)):
            p = geoms[k][SI['H1A']]
            ax.add_patch(Circle((p[0], p[1]), r, facecolor=fc, edgecolor='white', lw=0.5, zorder=z))
        axes[m] = (ax, d)
    return axes, scale


# ------------------------------------------------------------------------------------------------ scan
SCAN = ['A00', 'A01', 'A02', 'A03', 'B00', 'B01', 'B02', 'B03']
SCAN_FILES = {'A00': 'A_00_anchor', 'A01': 'A_01_d1.150', 'A02': 'A_02_d1.250', 'A03': 'A_03_d1.350',
              'B00': 'B_00_anchor', 'B01': 'B_01_d1.550', 'B02': 'B_02_d1.450', 'B03': 'B_03_d1.350'}


def scan_data():
    E0 = S.F('scan.A00.E')
    pts = {}
    for k in SCAN:
        d, E, conv = S.F(f'scan.{k}.d_NH'), S.F(f'scan.{k}.E'), S.F(f'scan.{k}.converged')
        Erel = (E - E0) * HA_KCAL
        assert abs(Erel - S.F(f'scan.{k}.Erel')) < 6e-4, (k, Erel)
        X = S.scan_point(SCAN_FILES[k])
        assert abs(S.dist(X, 'N1', 'H1A') - d) < 6e-5, (k, S.dist(X, 'N1', 'H1A'), d)
        pts[k] = dict(d=d, E=Erel, conv=bool(conv), delta=S.F(f'scan.{k}.delta'))
    assert [k for k in SCAN if not pts[k]['conv']] == ['A02'] and S.F('scan.n_unconverged') == 1
    # lower envelope (at d = 1.35 A the chain-B point is the lower one, fact scan.hysteresis_1.35 > 0)
    assert S.F('scan.hysteresis_1.35') > 0 and pts['A03']['E'] > pts['B03']['E']
    env = sorted([pts[k] for k in ('A00', 'A01', 'A02', 'B03', 'B02', 'B01', 'B00')], key=lambda r: r['d'])
    from scipy.interpolate import CubicSpline
    cs = CubicSpline([r['d'] for r in env], [r['E'] for r in env])
    xx = np.linspace(env[0]['d'], env[-1]['d'], 2001)
    yy = cs(xx)
    top = (xx[1:-1][np.argmax(yy[1:-1])], yy[1:-1].max())
    # the spline drawn here is the one behind the fact (maximum and its position)
    assert abs(top[1] - S.F('scan.barrier_spline')) < 2e-3, top
    assert abs(pts['A02']['E'] - S.F('scan.barrier')) < 1e-6 and abs(pts['B00']['E'] - S.F('scan.Omin')) < 1e-6
    return pts, (xx, yy)


def draw_scan(ax, fig):
    pts, (xx, yy) = scan_data()
    cA, cB = S.SIDE_COLOURS['N'], S.SIDE_COLOURS['O']
    fs = style.SMALL_PT
    zpe = S.F('zpe.NH')
    ax.plot(xx, yy, color='0.6', lw=0.9, zorder=1)
    for chain, col, mkr in (('A', cA, 'o'), ('B', cB, 'D')):
        ks = [k for k in SCAN if k[0] == chain]
        ax.plot([pts[k]['d'] for k in ks], [pts[k]['E'] for k in ks], color=col, lw=0.7, ls=(0, (1.5, 1.2)),
                zorder=2)
        for k in ks:
            conv = pts[k]['conv']
            ax.plot([pts[k]['d']], [pts[k]['E']], ls='none', marker=mkr, ms=5.0 if mkr == 'o' else 4.6,
                    mfc=col if conv else 'white', mec=col if not conv else 'white', mew=1.0 if not conv else 0.5,
                    zorder=4)
    ax.axhline(zpe, color='0.25', lw=0.8, ls=(0, (4, 2)), zorder=1)
    ax.axhline(0.0, color='0.8', lw=0.5, zorder=0)
    ax.set_xlim(*XL)
    ax.set_ylim(-0.95, 3.75)
    ax.set_xlabel('$d$(N1–H1A) (Å)')
    ax.set_ylabel('$E$ − $E$(salt) (kcal mol$^{-1}$)')
    ax.xaxis.set_major_formatter(__import__('matplotlib').ticker.FuncFormatter(lambda v, p: S.num(v, 1)))
    ax.yaxis.set_major_formatter(__import__('matplotlib').ticker.FuncFormatter(lambda v, p: S.num(v, 0)))
    ax.set_xticks([1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7])
    ax.set_yticks([0, 1, 2, 3])
    # ZPE label (upper line) and the barrier / end-point values
    zlab = ('½$h\\nu$(N1–H1A) = ' + S.txt('zpe.NH') + ' kcal mol$^{-1}$ (' + S.txt('zpe.nu_NH') +
            ' cm$^{-1}$)')
    ax.text(1.69, zpe + 0.10, zlab, fontsize=fs, ha='right', va='bottom', color='0.2')
    a2 = pts['A02']
    ax.annotate(f'{S.txt("scan.barrier")}', (a2['d'], a2['E']), xytext=(0, 7), textcoords='offset points',
                ha='center', va='bottom', fontsize=fs, color='0.15')
    b0 = pts['B00']
    ax.annotate(f'{S.txt("scan.Omin")}', (b0['d'], b0['E']), xytext=(0, -6), textcoords='offset points',
                ha='center', va='top', fontsize=fs, color='0.15')
    a0 = pts['A00']
    ax.annotate('salt\nN1–H1A···O5', (a0['d'], a0['E']), xytext=(3, -7), textcoords='offset points', ha='left',
                va='top', fontsize=fs, color=cA, linespacing=1.0)
    ax.annotate('neutral pair\nO5–H1A···N1', (b0['d'], b0['E']), xytext=(-12, 19), textcoords='offset points',
                ha='right', va='bottom', fontsize=fs, color=cB, linespacing=1.0)
    handles = [Line2D([], [], ls='none', marker='o', ms=5.0, mfc=cA, mec='white', mew=0.5,
                      label='chain A, from the salt'),
               Line2D([], [], ls='none', marker='D', ms=4.6, mfc=cB, mec='white', mew=0.5,
                      label='chain B, from the neutral pair'),
               Line2D([], [], ls='none', marker='o', ms=5.0, mfc='white', mec=cA, mew=1.0,
                      label='step limit, forces below thresholds'),
               Line2D([], [], color='0.6', lw=0.9, label='spline, lower envelope')]
    ax.legend(handles=handles, loc='upper right', bbox_to_anchor=(1.0, 0.86), fontsize=style.TICK_PT,
              handletextpad=0.4, labelspacing=0.3, borderaxespad=0.3)
    for k in SCAN:
        for f in ('d_NH', 'E', 'Erel', 'converged'):
            S.USED[f'scan.{k}.{f}']['printed'].append('plotted')
    S.F('scan.barrier_spline')
    S.USED['scan.barrier_spline']['printed'].append('plotted (spline maximum)')
    return pts


# ------------------------------------------------------------------------------------------------ functionals
FUNCS = [('b3lyp', 'B3LYP-\nD3(BJ)'), ('wb97xd', 'ωB97X-D'), ('pbe0', 'PBE0-\nD3(BJ)'), ('m062x', 'M06-2X'), ('mp2', 'DF-MP2')]


def draw_fcheck(ax, fig):
    fs = style.SMALL_PT
    top = [S.F(f'fcheck.{k}.top_N') for k, _ in FUNCS]
    neu = [S.F(f'fcheck.{k}.dE_O_N') for k, _ in FUNCS]
    for k, _ in FUNCS:
        for f in ('top_N', 'dE_O_N', 'E.end_N', 'E.top1', 'E.end_O'):
            assert S.status(f'fcheck.{k}.{f}') == '', (k, f, S.status(f'fcheck.{k}.{f}'))
        # the differences follow from the three total energies of fsp_scan.json
        e = {f: S.F(f'fcheck.{k}.E.{f}') for f in ('end_N', 'top1', 'end_O')}
        assert abs((e['top1'] - e['end_N']) * HA_KCAL - S.F(f'fcheck.{k}.top_N')) < 1e-3
        assert abs((e['end_O'] - e['end_N']) * HA_KCAL - S.F(f'fcheck.{k}.dE_O_N')) < 1e-3
    # B3LYP-D3(BJ) single points reproduce the scan, the three geometry files are the ones recorded, and every
    # method puts the neutral end point below the salt
    assert S.F('fcheck.b3lyp_consistency') < 1e-4
    assert all(S.F(f'fcheck.geom.{g}.sha256_ok') for g in ('end_N', 'top1', 'end_O'))
    assert S.F('fcheck.n_methods') == len(FUNCS) == S.F('fcheck.n_neutral_lower')
    x = np.arange(len(FUNCS))
    w = 0.37
    cT, cN = '#9A9A9A', S.SIDE_COLOURS['O']
    ax.axhline(0.0, color='0.3', lw=0.6, zorder=1)
    b1 = ax.bar(x - w / 2 - 0.02, top, w, color=cT, edgecolor='white', lw=0.4, zorder=2)
    b2 = ax.bar(x + w / 2 + 0.02, neu, w, color=cN, edgecolor='white', lw=0.4, zorder=2)
    for (k, _), xi, v1, v2 in zip(FUNCS, x, top, neu):
        s1, s2 = S.txt(f'fcheck.{k}.top_N', sign=True), S.txt(f'fcheck.{k}.dE_O_N', sign=True)
        for xx, v, s in ((xi - w / 2 - 0.02, v1, s1), (xi + w / 2 + 0.02, v2, s2)):
            ax.text(xx, v + (0.07 if v >= 0 else -0.07), s, ha='center', va='bottom' if v >= 0 else 'top',
                    fontsize=fs, color='0.15', rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels([n for _, n in FUNCS], fontsize=style.TICK_PT - 1.0, linespacing=0.95, va='top')
    ax.tick_params(axis='x', length=0, pad=2)
    ax.set_xlim(-0.6, len(FUNCS) - 0.4)
    ax.set_ylim(-3.6, 2.3)
    ax.set_yticks([-3, -2, -1, 0, 1, 2])
    ax.yaxis.set_major_formatter(__import__('matplotlib').ticker.FuncFormatter(lambda v, p: S.num(v, 0)))
    ax.set_ylabel('$E$ − $E$(salt) (kcal mol$^{-1}$)')
    handles = [Patch(color=cT, label='$d$ = 1.25 Å'), Patch(color=cN, label='neutral end')]
    ax.legend(handles=handles, loc='lower left', fontsize=style.TICK_PT, handlelength=1.2, handletextpad=0.4,
              labelspacing=0.3, borderaxespad=0.3)
    S.F('scan.A02.d_NH')
    S.USED['scan.A02.d_NH']['printed'].append('1.25')
    return b1, b2


# ------------------------------------------------------------------------------------------------ crystal profile
def profile_data():
    d = [S.F('proton.periodic.d_N1_H1A')] + [S.F(f'proton.periodic.{k}.d_N1_H1A') for k in PROF]
    e = [0.0] + [S.F(f'periodic.profile.{k}.PBE_kJ') for k in PROF]
    e0 = [0.0] + [S.F(f'periodic.profile.{k}.PBE0_kJ') for k in PROF]
    sl = [S.F(f'periodic.profile.{k}.slope_kJ') for k in PROF]
    dl = [S.F('proton.periodic.delta')] + [S.F(f'proton.periodic.{k}.delta') for k in PROF]
    # statements made in the caption: rising at every point with both functionals, every constraint force towards N1
    assert all(b > a for a, b in zip(e, e[1:])) and all(b > a for a, b in zip(e0, e0[1:])) and all(x > 0 for x in sl)
    assert S.F('periodic.profile.no_neutral_min') is True and S.F('periodic.iir1.final_sides') == ['N', 'N']
    for k, dd in zip(PROF, S.F('periodic.profile.targets')):
        assert abs(S.F(f'proton.periodic.{k}.d_N1_H1A') - dd) < 1e-4
    return d, e, e0, sl, dl


def draw_profile(ax):
    d, e, e0, sl, dl = profile_data()
    fs = style.SMALL_PT
    # H1A nearer O5 than N1 (delta > 0): shade from the delta = 0 crossing (linear between the bracketing points)
    i = next(j for j in range(1, len(dl)) if dl[j] > 0)
    x0 = d[i - 1] + (d[i] - d[i - 1]) * (0 - dl[i - 1]) / (dl[i] - dl[i - 1])
    ax.axvline(x0, ymin=0.0, ymax=0.42, color='0.55', lw=0.7, ls=(0, (1.2, 1.4)), zorder=1)
    ax.text(x0 + 0.008, -2.5, 'H1A nearer O5 →', ha='left', va='center', fontsize=fs, color='0.35')
    ax.axhline(0.0, color='0.8', lw=0.5, zorder=0)
    ax.plot(d, e, color=C_PBE, lw=0.9, zorder=2)
    ax.plot(d, e0, color=C_PBE0, lw=0.9, ls=(0, (3.0, 1.6)), zorder=2)
    ax.plot(d, e, ls='none', marker='o', ms=4.6, mfc=C_PBE, mec='white', mew=0.5, zorder=4)
    ax.plot(d, e0, ls='none', marker='D', ms=4.0, mfc=C_PBE0, mec='white', mew=0.5, zorder=4)
    # constraint force: a stroke of slope dE/dd through each constrained PBE point
    for x, y, m in zip(d[1:], e[1:], sl):
        h = 0.04
        ax.plot([x - h, x + h], [y - m * h, y + m * h], color='0.1', lw=1.1, zorder=3, solid_capstyle='butt')
    # the earlier constrained single configuration (O5-H1A held at 0.98 A)
    xf = S.F('proton.periodic.iirfix.d_N1_H1A')
    yf, yf0 = S.F('periodic.iirfix.dE_kJ'), S.F('periodic.hybrid.e100.best_kJ')
    ax.plot([xf], [yf], ls='none', marker='o', ms=4.6, mfc='white', mec='0.45', mew=0.9, zorder=5)
    ax.plot([xf], [yf0], ls='none', marker='D', ms=4.0, mfc='white', mec='0.45', mew=0.9, zorder=5)
    ax.annotate('O5–H1A held\nat ' + S.txt('protocol.qe_start_OH') + ' Å', (xf, yf), xytext=(-6, 0), textcoords='offset points',
                ha='right', va='center', fontsize=fs, color='0.35', linespacing=1.0)
    ax.annotate('salt', (d[0], 0.0), xytext=(1, -4), textcoords='offset points', ha='left', va='top', fontsize=fs,
                color=S.SIDE_COLOURS['N'])
    ax.set_xlim(*XL)
    ax.set_ylim(-9, 62)
    ax.set_xticks([1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7])
    ax.set_yticks([0, 20, 40, 60])
    ax.xaxis.set_major_formatter(__import__('matplotlib').ticker.FuncFormatter(lambda v, p: S.num(v, 1)))
    ax.set_xlabel('$d$(N1–H1A) held (Å)')
    ax.set_ylabel('$E$ − $E$(salt)\n(kJ mol$^{-1}$ per f.u.)')
    handles = [Line2D([], [], color=C_PBE, lw=0.9, marker='o', ms=4.6, mfc=C_PBE, mec='white', mew=0.5, label='PBE-D3(BJ)'),
               Line2D([], [], color=C_PBE0, lw=0.9, ls=(0, (3.0, 1.6)), marker='D', ms=4.0, mfc=C_PBE0, mec='white', mew=0.5,
                      label='PBE0-D3(BJ)'),
               Line2D([], [], color='0.1', lw=1.1, label='slope from the constraint force')]
    ax.legend(handles=handles, loc='upper left', fontsize=style.TICK_PT, handletextpad=0.4, labelspacing=0.3,
              borderaxespad=0.3, frameon=True, framealpha=1.0, edgecolor='none')
    for k in PROF:
        for f_ in ('PBE_kJ', 'PBE0_kJ', 'slope_kJ'):
            S.USED[f'periodic.profile.{k}.{f_}']['printed'].append('plotted')
    return d


HB3 = [('N1H1A_O5', 'N1···O5', S.HB_PAIR, 'o'), ('N1H1B_O6', 'N1–H1B···O6', S.HB_CRYSTAL, 's'),
       ('N1H1C_O6', 'N1–H1C···O6', S.HB_CRYSTAL, '^'), ('O4H4A_O5', 'O4–H4A···O5', S.HB_CRYSTAL, 'v')]


def draw_network(ax):
    fs = style.SMALL_PT
    d = [S.F('proton.periodic.d_N1_H1A')] + [S.F(f'proton.periodic.{k}.d_N1_H1A') for k in PROF]
    for key, lab, col, mk in HB3:
        y = [S.F(f'periodic.hb.{key}.DA')] + [S.F(f'periodic.profile.{k}.hb.{key}.DA') for k in PROF]
        ax.plot(d, y, color=col, lw=0.8, zorder=2)
        ax.plot(d, y, ls='none', marker=mk, ms=3.8, mfc=col, mec='white', mew=0.4, zorder=3)
        # labels: the two N1-H...O6 bonds at the right end (above / below), N1...O5 and O4-H4A...O5 at the left end
        at_left = key in ('N1H1A_O5', 'O4H4A_O5')
        dy = {'N1H1A_O5': 3.0, 'N1H1B_O6': -4.0, 'N1H1C_O6': 3.0, 'O4H4A_O5': 6.5}[key]
        ax.annotate(lab, (d[0], y[0]) if at_left else (d[-1], y[-1]), xytext=(1 if at_left else -2, dy), textcoords='offset points',
                    ha='left' if at_left else 'right', va='top' if dy < 0 else 'bottom', fontsize=fs, color='0.15')
        for k in PROF:
            S.USED[f'periodic.profile.{k}.hb.{key}.DA']['printed'].append('plotted')
    assert S.F('periodic.profile.p170.hb.N1H1C_O6.DA') > S.F('periodic.hb.N1H1C_O6.DA')
    ax.set_xlim(*XL)
    ax.set_ylim(2.46, 3.02)
    ax.set_xticks([1.1, 1.3, 1.5, 1.7])
    ax.xaxis.set_major_formatter(__import__('matplotlib').ticker.FuncFormatter(lambda v, p: S.num(v, 1)))
    ax.yaxis.set_major_formatter(__import__('matplotlib').ticker.FuncFormatter(lambda v, p: S.num(v, 1)))
    ax.set_yticks([2.5, 2.6, 2.7, 2.8, 2.9, 3.0])
    ax.set_xlabel('$d$(N1–H1A) held (Å)')
    ax.set_ylabel('$d$(D···A) (Å)')


# ------------------------------------------------------------------------------------------------ figure
def main():
    plt = S.setup()
    fig = S.figure_mm(W_MM, H_MM)
    fs = style.SMALL_PT
    top = H_MM - 3.0
    pw, ph, gx = 55.0, 35.0, 3.0
    xs = [1.0, 1.0 + pw + gx, 1.0 + 2 * (pw + gx)]
    yd = top - 5.0 - ph
    boxes = {'pcm': (xs[0], yd, pw, ph), 'gas': (xs[1], yd, pw, ph), 'periodic': (xs[2], yd, pw, ph)}
    axes, scale = draw_overlays(fig, boxes)
    titles = {'pcm': 'IEF-PCM (water): H1A on N1', 'gas': 'Gas phase: H1A on O5',
              'periodic': 'Periodic PBE-D3(BJ): H1A on N1'}
    notes = {
        'pcm': [f'$\\delta$ = {S.txt("proton.pcm.delta")} Å, N1···O5 {S.txt("proton.pcm.d_N1_O5")} Å',
                f'RMSD cation {S.txt("pcm.geom.kabsch_cation")} Å, anion {S.txt("pcm.geom.kabsch_anion")} Å'],
        'gas': [f'$\\delta$ = {S.txt("proton.gas.delta", sign=True)} Å, N1···O5 {S.txt("proton.gas.d_N1_O5")} Å',
                f'RMSD cation {S.txt("gas.geom.kabsch_cation")} Å, anion {S.txt("gas.geom.kabsch_anion")} Å'],
        'periodic': [f'$\\delta$ = {S.txt("proton.periodic.delta")} Å, N1···O5 {S.txt("proton.periodic.d_N1_O5")} Å',
                     f'heavy-atom RMSD, whole cell {S.txt("periodic.rmsd_heavy")} Å'],
    }
    W, H = fig._saph_mm
    for (m, (x, y, w, h)), letter in zip(boxes.items(), 'abc'):
        S.panel_letter(fig, x + 0.2, top, letter, titles[m], fontsize=style.TICK_PT + 0.5)
        for k, line in enumerate(notes[m]):
            fig.text((x + 2.0) / W, (y - 1.5 - 3.1 * k) / H, line, fontsize=fs, ha='left', va='top', color='0.15')
    # shared overlay legend (row 1, right of the notes)
    lh = [Line2D([], [], color=S.MODEL_COLOURS['xray'], lw=0.8, label='X-ray'),
          Line2D([], [], color=S.MODEL_COLOURS['pcm'], lw=2.2, label='IEF-PCM (water)'),
          Line2D([], [], color=S.MODEL_COLOURS['gas'], lw=2.2, label='gas phase'),
          Line2D([], [], color=S.MODEL_COLOURS['periodic'], lw=2.2, label='periodic model'),
          Line2D([], [], color='0.3', lw=0.8, ls=(0, (2, 1.4)), label='hydrogen bond')]
    y_leg = yd - 12.4
    fig.legend(handles=lh, loc='upper center', bbox_to_anchor=(0.5, y_leg / H), ncol=5, fontsize=style.TICK_PT,
               handlelength=2.0, columnspacing=1.4, borderaxespad=0.0, frameon=False)
    # rows 2 and 3
    y3, h3 = 11.0, 38.0
    y2 = y3 + h3 + 15.0
    h2 = y_leg - 6.5 - 7.5 - y2
    axs = S.axes_mm(fig, 15.0, y2, 92.0, h2)
    draw_scan(axs, fig)
    S.panel_letter(fig, 0.4, y2 + h2 + 2.2, 'd', 'Relaxed proton scan, IEF-PCM (water)',
                   fontsize=style.TICK_PT + 0.5)
    axf = S.axes_mm(fig, 123.5, y2 + 4.5, 50.3, h2 - 12.5)
    draw_fcheck(axf, fig)
    S.panel_letter(fig, 112.0, y2 + h2 + 2.2, 'e', 'Functional check', fontsize=style.TICK_PT + 0.5)
    fig.text(112.0 / W, (y2 + h2 - 1.6) / H, 'single points at the scan geometries,', fontsize=fs,
             ha='left', va='top', color='0.25')
    fig.text(112.0 / W, (y2 + h2 - 4.4) / H, 'IEF-PCM (water), 6-311++G(d,p)', fontsize=fs,
             ha='left', va='top', color='0.25')
    axp = S.axes_mm(fig, 15.0, y3, 92.0, h3)
    draw_profile(axp)
    S.panel_letter(fig, 0.4, y3 + h3 + 2.2, 'f', 'Constrained N1–H1A profile, periodic cell', fontsize=style.TICK_PT + 0.5)
    axn = S.axes_mm(fig, 123.5, y3, 50.3, h3)
    draw_network(axn)
    S.panel_letter(fig, 112.0, y3 + h3 + 2.2, 'g', 'Hydrogen bonds along the profile', fontsize=style.TICK_PT + 0.5)
    S.save(fig, STEM, W_MM)


if __name__ == '__main__':
    main()
