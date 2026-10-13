"""Dielectric ladder: where the proton of N1-H1A...O5 sits as the continuum dielectric constant grows (SINGLE width).

    cd /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate
    nice -n 10 env OMP_NUM_THREADS=2 /root/workspace/cqc-dft/tools/paper-env/bin/python \
        paper/main/figures/scripts/fig03_dielectric.py

Three panels on a common logarithmic epsilon axis, ion pair at B3LYP-D3(BJ)/6-311++G(d,p):
    (a) proton metric delta = d(N1-H1A) - d(H1A...O5)
    (b) P1-O5 and (c) C1-N1, the two heavy-atom bonds that report the protonation state without relying on the
        position of H1A (P1-O5 lengthens when O5 carries the proton, C1-N1 shortens when N1 loses it)
epsilon = 1 is the gas-phase optimum (from the X-ray structure); epsilon = 4, 8, 20 started from the water optimum;
water (78.36) shows both minima: the salt (primary optimum) and the neutral O5-protonated pair (end point B of
the relaxed scan).  Reference lines: the X-ray structure and the periodic PBE-D3(BJ) model (relax i).
Every plotted value is a fact from paper/main/data/facts.json and is re-measured here from the geometry files
(asserted to 1e-4 A); the crossover bracket comes from facts eps.bracket.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FixedLocator, NullLocator, FixedFormatter  # noqa: E402

import _saph as S  # noqa: E402
from figkit import style  # noqa: E402

STEM = 'fig03_dielectric'
W_MM = 85.0

# environment -> (facts prefix, epsilon, geometry loader)
EPS_WATER = None


def data():
    eps_w = S.F('protocol.eps_water')
    assert S.F('eps.values') == [1, 4, 8, 20, eps_w]
    rows = [
        ('gas', 'proton.gas', 1.0, lambda: S.case_pair('sfa_pair_gas_d3')),
        ('eps4', 'proton.eps4', 4.0, lambda: S.eps_point(4)),
        ('eps8', 'proton.eps8', 8.0, lambda: S.eps_point(8)),
        ('eps20', 'proton.eps20', 20.0, lambda: S.eps_point(20)),
        ('pcm', 'proton.pcm', eps_w, lambda: S.case_pair('sfa_pair_pcm_d3')),
        ('scanO', 'proton.scanO', eps_w, lambda: S.scan_point('B_00_anchor')),
    ]
    out = []
    for name, pre, eps, load in rows:
        X = load()
        r = dict(name=name, eps=eps, side=S.F(pre + '.side'), delta=S.F(pre + '.delta'),
                 PO5=S.F(pre + '.P1O5'), CN=S.F(pre + '.C1N1'))
        # independent check of every plotted value against the geometry file
        assert abs(S.delta(X) - r['delta']) < 1e-4, (name, S.delta(X), r['delta'])
        assert abs(S.dist(X, 'P1', 'O5') - r['PO5']) < 1e-4, name
        assert abs(S.dist(X, 'C1', 'N1') - r['CN']) < 1e-4, name
        assert (r['delta'] < 0) == (r['side'] == 'N'), name
        out.append(r)
    ref = dict(xray=dict(delta=S.F('proton.xray.delta'), PO5=S.F('proton.xray.P1O5'), CN=S.F('proton.xray.C1N1')),
               periodic=dict(delta=S.F('proton.periodic.delta'), PO5=S.F('proton.periodic.P1O5'),
                             CN=S.F('proton.periodic.C1N1')))
    Xx = S.xray_pair()
    assert abs(S.delta(Xx) - ref['xray']['delta']) < 1e-6
    assert abs(S.dist(Xx, 'P1', 'O5') - ref['xray']['PO5']) < 1e-6
    Xp, _ = S.periodic_pair('relax_i')
    assert abs(S.delta(Xp) - ref['periodic']['delta']) < 1e-6
    assert abs(S.dist(Xp, 'C1', 'N1') - ref['periodic']['CN']) < 1e-6
    br = S.F('eps.bracket')
    assert br == '4 < ε_crossover < 8', br
    return out, ref, eps_w


def main():
    plt = S.setup()
    rows, ref, eps_w = data()
    fig = S.figure_mm(W_MM, 118.0)
    L, Wax = 16.5, 64.0
    hs = [40.0, 25.0, 25.0]
    gap = 3.0
    y0 = 10.5
    ys = [y0 + hs[2] + gap + hs[1] + gap, y0 + hs[2] + gap, y0]
    axes = [S.axes_mm(fig, L, y, Wax, h) for y, h in zip(ys, hs)]
    keys = [('delta', '$\\delta$ (Å)'), ('PO5', '$d$(P1–O5) (Å)'), ('CN', '$d$(C1–N1) (Å)')]
    xlim = (0.8, 110.0)
    ticks = [1, 4, 8, 20, eps_w]
    ticklabels = ['1', '4', '8', '20', S.num(eps_w, 1)]
    S.USED['protocol.eps_water']['printed'].append(ticklabels[-1])
    cN, cO = S.SIDE_COLOURS['N'], S.SIDE_COLOURS['O']
    mk = {'N': dict(marker='o', ms=5.0, mfc=cN, mec='white', mew=0.5, color=cN),
          'O': dict(marker='D', ms=4.6, mfc=cO, mec='white', mew=0.5, color=cO)}
    for k, (ax, (key, ylab)) in enumerate(zip(axes, keys)):
        ax.set_xscale('log')
        ax.set_xlim(*xlim)
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_minor_locator(NullLocator())
        if k == 2:
            ax.xaxis.set_major_formatter(FixedFormatter(ticklabels))
            ax.set_xlabel('Dielectric constant $\\varepsilon$ (dimensionless)')
        else:
            ax.tick_params(labelbottom=False)
        ax.set_ylabel(ylab)
        ax.axvspan(4.0, 8.0, color='0.90', lw=0, zorder=0)
        # reference lines
        ax.axhline(ref['xray'][key], color=S.MODEL_COLOURS['xray'], lw=0.8, ls=(0, (4, 2)), zorder=1)
        ax.axhline(ref['periodic'][key], color=S.MODEL_COLOURS['periodic'], lw=1.0, ls=(0, (1, 1.3)), zorder=1)
        # connectors within one protonation state
        for side in ('O', 'N'):
            pts = sorted([(r['eps'], r[key]) for r in rows if r['side'] == side and r['name'] != 'scanO'])
            ax.plot([p[0] for p in pts], [p[1] for p in pts], color=S.SIDE_COLOURS[side], lw=0.8, zorder=2)
        for r in rows:
            ax.plot([r['eps']], [r[key]], ls='none', zorder=4, **mk[r['side']])
        if key == 'delta':
            ax.axhline(0.0, color='0.55', lw=0.5, zorder=1)
            ax.set_ylim(-1.15, 1.0)
            ax.set_yticks([-1.0, -0.5, 0.0, 0.5])
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: S.num(v, 1)))
            continue
        elif key == 'PO5':
            ax.set_ylim(1.49, 1.595)
            ax.set_yticks([1.50, 1.52, 1.54, 1.56, 1.58])
        else:
            ax.set_ylim(1.40, 1.475)
            ax.set_yticks([1.40, 1.42, 1.44, 1.46])
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: S.num(v, 2)))
    axd, axp, axc = axes
    fs = style.SMALL_PT
    # region labels and direct labels (a)
    axd.text(1.9, 0.86, 'proton on O5', color=cO, fontsize=fs, ha='center', va='center')
    axd.text(1.9, -0.14, 'proton on N1', color=cN, fontsize=fs, ha='center', va='center')
    # the band marks where a salt minimum first appears; the text says it is no thermodynamic crossover (5 Oct 2026)
    axd.text(np.sqrt(32.0), 0.33, 'first salt', color='0.3', fontsize=fs, ha='center', va='center')
    axd.text(np.sqrt(32.0), 0.20, 'minimum', color='0.3', fontsize=fs, ha='center', va='center')
    sO = [r for r in rows if r['name'] == 'scanO'][0]
    sN = [r for r in rows if r['name'] == 'pcm'][0]
    axd.annotate('neutral\nminimum', (eps_w, sO['delta']), xytext=(-5, 0), textcoords='offset points',
                 ha='right', va='center', fontsize=fs, color='0.2', linespacing=0.95)
    axd.annotate('salt minimum', (eps_w, sN['delta']), xytext=(-4, -3.5), textcoords='offset points',
                 ha='right', va='top', fontsize=fs, color='0.2')
    axd.text(1.0, ref['xray']['delta'] + 0.04, 'X-ray (H1A as refined)', fontsize=fs, color=S.MODEL_COLOURS['xray'],
             ha='left', va='bottom', transform=axd.get_yaxis_transform() if False else axd.transData)
    axd.text(1.0, ref['periodic']['delta'] - 0.045, 'periodic model', fontsize=fs, color='#00664A',
             ha='left', va='top')
    # (b), (c): reference-line labels at the left, clear of the data
    axp.text(1.0, ref['xray']['PO5'] - 0.003, 'X-ray', fontsize=fs, color=S.MODEL_COLOURS['xray'], ha='left',
             va='top')
    axp.text(1.0, ref['periodic']['PO5'] + 0.003, 'periodic', fontsize=fs, color='#00664A', ha='left',
             va='bottom')
    axc.text(1.0, ref['xray']['CN'] + 0.003, 'X-ray', fontsize=fs, color=S.MODEL_COLOURS['xray'], ha='left',
             va='bottom')
    axc.text(1.62, ref['periodic']['CN'] - 0.003, 'periodic', fontsize=fs, color='#00664A', ha='left',
             va='top')
    for y, h, letter in zip(ys, hs, 'abc'):
        S.panel_letter(fig, 0.4, y + h - 2.6, letter)
    # legend above (a)
    handles = [Line2D([], [], ls='none', label='proton on N1 (salt)', **mk['N']),
               Line2D([], [], ls='none', label='proton on O5 (neutral pair)', **mk['O']),
               Line2D([], [], color=S.MODEL_COLOURS['xray'], lw=0.8, ls=(0, (4, 2)), label='X-ray'),
               Line2D([], [], color=S.MODEL_COLOURS['periodic'], lw=1.0, ls=(0, (1, 1.3)),
                      label='periodic PBE-D3(BJ)')]
    axd.legend(handles=handles, loc='lower left', bbox_to_anchor=(-0.02, 1.02), ncol=2, fontsize=style.TICK_PT,
               handlelength=1.8, columnspacing=1.0, handletextpad=0.5, borderaxespad=0.0, frameon=False)
    # facts printed through the plotted data (for the audit record)
    for pre in ('proton.gas', 'proton.eps4', 'proton.eps8', 'proton.eps20', 'proton.pcm', 'proton.scanO',
                'proton.xray', 'proton.periodic'):
        for f in ('.delta', '.P1O5', '.C1N1'):
            S.USED[pre + f]['printed'].append('plotted')
    S.save(fig, STEM, W_MM)


if __name__ == '__main__':
    main()
