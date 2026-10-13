"""Fig. 11 (fig11_torsion): the C3-C4-S1-N2 torsion from the crystal to the DFT minima, the crystal ligands and the
docked poses.

    nice -n 10 tools/paper-env/bin/python paper/main/figures/scripts/fig11_torsion.py

tau_fold = max(|C3-C4-S1-N2|, |C5-C4-S1-N2|) (ring-flip invariant; 90 deg: S1-N2 perpendicular to the ring plane,
180 deg: in the plane), as defined in dock/RESULTS.md section 6.
(a) Relaxed B3LYP-D3(BJ)/6-311++G(d,p) IEF-PCM torsion scans (VM burst 3, referee M9; facts section 35) of neutral
    sulfanilamide and the sulfonamidate: every computed point (markers) and the cubic spline through the points and the
    minimum (line), relative to the unconstrained minimum of the same species; thin dashed lines: the MMFF94s relaxed
    profiles of the same species (dock/results/production_summary.json), for comparison.  Black triangles: the torsions of
    the ligands of the deposited complexes 3TZF, 1AJ0 and 6RL9 (facts dock.cocrystal.*.tau) on the sulfanilamide curve;
    dashed vertical: the SAPH crystal.  The rigid DFT single points of burst 2 are retired (facts section 36).
(b) tau_fold of the top pose of every seed of every production job of the three species (DFT input geometries), one
    row per target and scoring function (the 280 poses of facts dock.tors.n_poses); black triangles: the crystal-ligand
    torsion of that receptor (open in the 3HS4 rows: the 6RL9 sulfanilamide carried over, the cross-docking reference).
Updated 29 Sep 2026 (VM burst 3): panel (a) now shows the relaxed DFT scans, the crystal ligands are marked in both panels.
"""
import json

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import PchipInterpolator

import _figb as B
from figkit import style

STEM = 'fig11_torsion'
SUMMARY = B.DOCK / 'results' / 'production_summary.json'
SPECIES = [('sfa_neutral', 'sulfanilamide', 'o', B.FRAG_COLOURS['sulfanilamide']),
           ('sfa_anion', 'sulfonamidate', 'D', B.FRAG_COLOURS['sulfonamidate']),
           ('sfa_cation', 'cation', 's', B.FRAG_COLOURS['cation'])]
ROWS = [('T1_3TZF', 'vina', '3TZF, Vina'), ('T1_3TZF', 'vinardo', '3TZF, Vinardo'),
        ('T5_7S2K', 'vina', '7S2K, Vina'), ('T5_7S2K', 'vinardo', '7S2K, Vinardo'),
        ('T3_6RL9', 'ad4zn', '6RL9, AD4Zn'), ('T3_6RL9', 'vina', '6RL9, Vina'), ('T3_6RL9', 'vinardo', '6RL9, Vinardo'),
        ('T4_3HS4', 'ad4zn', '3HS4, AD4Zn'), ('T4_3HS4', 'vina', '3HS4, Vina'), ('T4_3HS4', 'vinardo', '3HS4, Vinardo')]
XTAL = {'T1_3TZF': ('3TZF', True), 'T3_6RL9': ('6RL9', True), 'T4_3HS4': ('6RL9', False)}   # (source complex, own ligand)


def fold(t):
    t = np.abs(np.asarray(t, float))
    return np.maximum(t, 180.0 - t)


def scan(sp):
    """Computed points (tau_fold, strain) of the relaxed scan from facts, and a shape-preserving (PCHIP) curve through them
    and the minimum for drawing (a cubic spline overshoots below zero between the 90.0 grid point and the 90.4 minimum)."""
    F = B.facts()
    import re
    keys = [k for k in F if re.fullmatch(rf'tors\.scan\.{sp}\.t\d{{4}}', k)]
    pts = sorted((F[k + '.tau']['value'], F[k]['value']) for k in keys)
    for k in keys:
        B.USED.setdefault(k, 'plotted')
    x = [B.F(f'tors.ref.{sp}.tau')] + [p[0] for p in pts]
    y = [0.0] + [p[1] for p in pts]
    o = np.argsort(x)
    return np.array(pts), PchipInterpolator(np.array(x)[o], np.array(y)[o])


def main():
    style.apply()
    S = json.loads(SUMMARY.read_text())
    MM = style.MM
    W = style.ONE_HALF
    H = 132 * MM
    fig = plt.figure(figsize=(W, H))
    left, right = 22.0 * MM, 2.0 * MM
    ha_top, ha_h = 5.5 * MM, 40.0 * MM
    hb_h = 64.0 * MM
    bottom = 10.0 * MM
    axa = fig.add_axes([left / W, (H - ha_top - ha_h) / H, (W - left - right) / W, ha_h / H])
    axb = fig.add_axes([left / W, bottom / H, (W - left - right) / W, hb_h / H], sharex=axa)
    xlim = (85, 180)
    fs = style.SMALL_PT

    # ---------------- (a) relaxed DFT scans, MMFF94s for comparison
    splines = {}
    for sp, lab, mk, col in SPECIES[:2]:
        pts, cs = scan(sp)
        splines[sp] = cs
        xx = np.linspace(B.F(f'tors.ref.{sp}.tau') if sp == 'sfa_neutral' else 90.0, 180.0, 400)
        axa.plot(xx, cs(xx), color=col, lw=1.1, zorder=3)
        axa.plot(pts[:, 0], pts[:, 1], ls='none', marker=mk, ms=3.0, mfc=col, mec='white', mew=0.35, zorder=4)
        prof = np.array(S['mmff_torsion_profiles'][sp]['relaxed'], float)
        f = fold(prof[:, 0])
        xs = np.unique(np.round(f, 3))
        ys = np.array([prof[np.isclose(f, x, atol=1e-3), 1].min() for x in xs])
        axa.plot(xs, ys, color=col, lw=0.7, ls=(0, (2.4, 1.6)), zorder=2)
        # the scan minimum (the unconstrained re-optimisation)
        axa.plot(B.F(f'tors.ref.{sp}.tau'), 0.0, marker=mk, ms=3.8, mfc='white', mec=col, mew=0.9, ls='none', zorder=5,
                 clip_on=False)
    # every relaxed point converged; the spline reproduces the scan and the facts it prints
    assert B.F('tors.scan.sfa_neutral.n_converged') == B.F('tors.scan.sfa_neutral.n_points')
    assert B.F('tors.scan.sfa_anion.n_converged') == B.F('tors.scan.sfa_anion.n_points')
    tau_x = B.F('dock.tors.crystal')
    axa.axvline(tau_x, ymin=0.0, ymax=0.15, color='0.15', lw=0.7, ls=(0, (4, 2)), zorder=1)   # short: clear of labels and legend
    axb.axvline(tau_x, color='0.15', lw=0.7, ls=(0, (4, 2)), zorder=1)
    # crystal ligands of the deposited complexes on the sulfanilamide curve
    XT = [('3TZF', -4.0, 3.5, 'right'), ('1AJ0', -2.0, 9.0, 'right'), ('6RL9', -4.0, 3.5, 'right')]
    for code, dx, dy, ha in XT:
        t = B.F(f'dock.cocrystal.{code}.tau')
        e = B.F(f'tors.strain.sfa_neutral.cocrystal.{code}')
        assert abs(float(splines['sfa_neutral'](t)) - e) < 0.03, code
        axa.plot(t, e, marker='^', ms=4.6, mfc='0.1', mec='white', mew=0.4, ls='none', zorder=6)
        axa.annotate(code, (t, e), xytext=(dx, dy), textcoords='offset points', ha=ha, va='bottom', fontsize=fs, color='0.1')
    axa.set_ylim(-0.25, 2.25)
    axa.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
    axa.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: B.num(v, 1)))
    axa.set_ylabel('Relaxed strain\n(kcal mol$^{-1}$)')
    axa.axhline(0.0, color='0.8', lw=0.5, zorder=0)
    axa.text(tau_x + 1.2, -0.235, 'SAPH crystal ' + B.fmt('dock.tors.crystal', 1) + '°', ha='left', va='bottom', fontsize=fs)
    hl = [Line2D([], [], color=col, lw=1.1, marker=mk, ms=3.2, mfc=col, mec='white', mew=0.3, label=f'{lab}, DFT')
          for sp, lab, mk, col in SPECIES[:2]]
    hl += [Line2D([], [], color='0.4', lw=0.7, ls=(0, (2.4, 1.6)), label='MMFF94s, same species'),
           Line2D([], [], ls='none', marker='^', ms=4.4, mfc='0.1', mec='white', mew=0.4, label='crystal ligand of the complex')]
    axa.legend(handles=hl, loc='upper left', bbox_to_anchor=(0.0, 1.0), fontsize=style.TICK_PT - 0.5, handlelength=2.0,
               borderaxespad=0.3, labelspacing=0.25, ncol=1, frameon=True, framealpha=1.0, edgecolor='none')
    plt.setp(axa.get_xticklabels(), visible=False)
    B.panel_letter(fig, axa, 'a', 'Relaxed torsion scans, B3LYP-D3(BJ), IEF-PCM (water)', dx_mm=-left / MM + 0.6)

    # ---------------- (b) docked poses, with the crystal-ligand torsion of each receptor
    n = 0
    lo, hi = 999, -999
    ylab = []
    for r_i, (tg, scoring, lab) in enumerate(ROWS):
        y0 = len(ROWS) - 1 - r_i
        ylab.append((y0, lab))
        for s_i, (sp, slab, mk, col) in enumerate(SPECIES):
            key = f'{tg}/{sp}__{scoring}__prod'
            if key not in S['jobs']:
                continue
            ps = S['jobs'][key]['torsion']['per_seed']
            t = np.array([p['tau_fold'] for p in ps])
            n += len(t)
            lo, hi = min(lo, t.min()), max(hi, t.max())
            yy = y0 + (0.22, 0.0, -0.22)[s_i]
            axb.plot(t, np.full(len(t), yy), ls='none', marker=mk, ms=3.0, mfc=col, mec='white', mew=0.3,
                     alpha=0.9, zorder=3)
        if tg in XTAL:
            code, own = XTAL[tg]
            axb.plot(B.F(f'dock.cocrystal.{code}.tau'), y0, marker='^', ms=4.8, ls='none', zorder=5,
                     mfc='0.1' if own else 'white', mec='0.1' if own else '0.1', mew=0.4 if own else 0.8)
    assert n == B.F('dock.tors.n_poses'), n
    assert abs(lo - B.F('dock.tors.range_min')) < 0.05 and abs(hi - B.F('dock.tors.range_max')) < 0.05, (lo, hi)
    for sp, lab, mk, col in SPECIES:
        axb.axvline(B.F(f'dock.tors.dft.{sp}'), color=col, lw=0.6, ls=(0, (1.2, 1.4)), zorder=1)
    axb.set_yticks([y for y, _ in ylab])
    axb.set_yticklabels([l for _, l in ylab], fontsize=style.TICK_PT)
    axb.set_ylim(-0.6, len(ROWS) - 0.4)
    for y in (2.5, 5.5, 7.5):
        axb.axhline(y, color='0.8', lw=0.5, zorder=0)
    axb.set_xlim(*xlim)
    axb.set_xlabel(r'$\tau_\mathrm{fold}$ = max(|C3–C4–S1–N2|, |C5–C4–S1–N2|) (°)')
    axb.tick_params(axis='y', length=0, pad=2)
    axb.text(0.985, 0.985, f'{B.fmt("dock.tors.n_poses", 0)} top poses: '
             + B.fmt('dock.tors.range_min', 1) + '–' + B.fmt('dock.tors.range_max', 1) + '°',
             transform=axb.transAxes, ha='right', va='top', fontsize=fs,
             bbox=dict(boxstyle='square,pad=0.2', fc='white', ec='none'))
    hs = [Line2D([], [], ls='none', marker=mk, ms=3.4, mfc=col, mec='white', mew=0.3, label=slab)
          for sp, slab, mk, col in SPECIES]
    hs += [Line2D([], [], ls='none', marker='^', ms=4.4, mfc='0.1', mec='0.1', mew=0.4, label='crystal ligand')]
    axb.legend(handles=hs, loc='lower right', bbox_to_anchor=(1.0, 0.0), fontsize=style.TICK_PT - 0.5, handletextpad=0.3,
               borderaxespad=0.3, labelspacing=0.25, frameon=True, framealpha=1.0, edgecolor='none')
    for k in ('dock.cocrystal.3TZF.tau', 'dock.cocrystal.1AJ0.tau', 'dock.cocrystal.6RL9.tau'):
        B.USED[k] = B.fmt(k, 1)
    B.panel_letter(fig, axb, 'b', 'Top poses of the three species (DFT input geometries)', dx_mm=-left / MM + 0.6)
    B.save(fig, STEM, width_mm=W / MM, extra=dict(n_poses=n))


if __name__ == '__main__':
    main()
