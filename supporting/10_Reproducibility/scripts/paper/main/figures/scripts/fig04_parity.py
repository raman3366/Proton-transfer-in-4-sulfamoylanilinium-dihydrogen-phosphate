"""X-ray against optimised heavy-atom bond lengths, by bond class (DOUBLE width).

    cd /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate
    nice -n 10 env OMP_NUM_THREADS=2 /root/workspace/cqc-dft/tools/paper-env/bin/python \
        paper/main/figures/scripts/fig04_parity.py

(a) Parity plot of the 15 heavy-atom bonds of the ion pair: IEF-PCM (water) and gas-phase optima at
    B3LYP-D3(BJ)/6-311++G(d,p) and the periodic PBE-D3(BJ) model (relax i, formula unit 1) against the X-ray
    structure; the bond classes of Table 3 are labelled along the diagonal; y = x with a +-0.02 A band.
(b) Signed deviation (model - X-ray) of every heavy-atom bond, grouped by the same classes, same markers.
Bonds to hydrogen are left out (riding C-H, restrained O-H in the X-ray model).  Values: facts bulk.bond.<model>.<bond>
(re-measured here from the geometry files, asserted to 1e-4 A); legend MAE: facts <model>.geom.heavy.mae.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

import _saph as S  # noqa: E402
from figkit import style  # noqa: E402

STEM = 'fig04_parity'
W_MM = 175.0

# bond classes in the order and wording of Table 3 (tbl:geom); bond keys as in facts bulk.bond.*
CLASSES = [
    ('C–C (ring)', ['C1-C2', 'C2-C3', 'C3-C4', 'C4-C5', 'C5-C6', 'C1-C6']),
    ('C1–N1, C4–S1', ['C1-N1', 'C4-S1']),
    ('S1–O1, S1–O2', ['O1-S1', 'O2-S1']),
    ('S1–N2', ['N2-S1']),
    ('P1–O3, P1–O4 (P–OH)', ['P1-O3', 'P1-O4']),
    ('P1–O5, P1–O6', ['P1-O5', 'P1-O6']),
]
MODELS = [('pcm', 'pcm', 'sfa_pair_pcm_d3'), ('gas', 'gas', 'sfa_pair_gas_d3'), ('periodic', 'periodic_i', None)]
FACT_PREFIX = {'pcm': 'pcm', 'gas': 'gas', 'periodic': 'periodic'}


def pretty(b):
    """'O1-S1' -> 'S1–O1' (sulfur/phosphorus first, as in the text), 'C1-C2' -> 'C1–C2'."""
    a, c = b.split('-')
    if c[0] in 'SP' and a[0] in 'ON':
        a, c = c, a
    return f'{a}–{c}'


def data():
    geoms = {'xray': S.xray_pair(), 'pcm': S.case_pair('sfa_pair_pcm_d3'), 'gas': S.case_pair('sfa_pair_gas_d3'),
             'periodic': S.periodic_pair('relax_i')[0]}
    fkey = {'xray': 'xray', 'pcm': 'pcm', 'gas': 'gas', 'periodic': 'periodic_i'}
    tab = {}
    for cname, bonds in CLASSES:
        for b in bonds:
            a, c = b.split('-')
            row = {}
            for m, X in geoms.items():
                v = S.F(f'bulk.bond.{fkey[m]}.{b}')
                assert abs(S.dist(X, a, c) - v) < 1e-4, (m, b, S.dist(X, a, c), v)
                row[m] = v
            tab[b] = row
    assert len(tab) == 15
    # the class MAE values of Table 3 follow from these bonds (consistency check with facts)
    for m in ('pcm', 'gas', 'periodic'):
        mae = np.mean([abs(tab[b][m] - tab[b]['xray']) for b in tab])
        assert abs(mae - S.F(f'{FACT_PREFIX[m]}.geom.heavy.mae')) < 1e-6, (m, mae)
    return tab


def main():
    plt = S.setup()
    tab = data()
    fig = S.figure_mm(W_MM, 86.0)
    fs = style.SMALL_PT
    mk = {m: dict(marker=S.MODEL_MARKERS[m], color=S.MODEL_COLOURS[m]) for m in ('pcm', 'gas', 'periodic')}
    ms = {'pcm': 4.2, 'gas': 3.9, 'periodic': 4.4}
    style_of = {  # filled PCM, open gas and periodic so coinciding points stay visible
        'pcm': dict(mfc=S.MODEL_COLOURS['pcm'], mec='white', mew=0.4),
        'gas': dict(mfc='none', mec=S.MODEL_COLOURS['gas'], mew=0.85),
        'periodic': dict(mfc='none', mec=S.MODEL_COLOURS['periodic'], mew=0.85),
    }
    # ------------------------------------------------------------------ (a) parity
    ax = S.axes_mm(fig, 14.0, 11.0, 68.0, 68.0)
    lim = (1.36, 1.81)
    xx = np.array(lim)
    band = 0.02
    ax.fill_between(xx, xx - band, xx + band, color='0.90', lw=0, zorder=0)
    ax.plot(xx, xx, color='0.35', lw=0.7, zorder=1)
    for m in ('pcm', 'gas', 'periodic'):
        x = [tab[b]['xray'] for b in tab]
        y = [tab[b][m] for b in tab]
        ax.plot(x, y, ls='none', ms=ms[m], zorder={'pcm': 3, 'gas': 4, 'periodic': 5}[m], **mk[m], **style_of[m])
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('X-ray bond length (Å)')
    ax.set_ylabel('Optimised bond length (Å)')
    ticks = [1.40, 1.50, 1.60, 1.70, 1.80]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    fmt2 = plt.FuncFormatter(lambda v, p: S.num(v, 2))
    ax.xaxis.set_major_formatter(fmt2)
    ax.yaxis.set_major_formatter(fmt2)
    # class labels to the lower right of their cluster (the optimised bonds lie on or above the diagonal), each with
    # a short leader to the cluster; C4-S1 sits at the top right, so its label goes to the left
    clusters = [('C–C (ring)', ['C1-C2', 'C2-C3', 'C3-C4', 'C4-C5', 'C5-C6', 'C1-C6'], 'right'),
                ('S1–O1, S1–O2', ['O1-S1', 'O2-S1'], 'right'),
                ('C1–N1', ['C1-N1'], 'right'),
                ('P1–O5, P1–O6', ['P1-O5', 'P1-O6'], 'right'),
                ('P1–O3, P1–O4', ['P1-O3', 'P1-O4'], 'right'),
                ('S1–N2', ['N2-S1'], 'right'),
                ('C4–S1', ['C4-S1'], 'left')]
    for nm, bl, side in clusters:
        x0 = float(np.mean([tab[b]['xray'] for b in bl]))
        if side == 'right':
            tx, ty = x0 + 0.034, x0 - 0.006
            ax.plot([x0 + 0.006, tx - 0.003], [x0 - 0.002, ty], color='0.55', lw=0.45, zorder=2)
            ax.text(tx, ty, nm, fontsize=fs, ha='left', va='center', color='0.2')
        else:
            tx, ty = x0 - 0.03, x0 + 0.018
            ax.plot([x0 - 0.006, tx + 0.003], [x0 + 0.012, ty], color='0.55', lw=0.45, zorder=2)
            ax.text(tx, ty, nm, fontsize=fs, ha='right', va='center', color='0.2')
    maes = {m: S.txt(f'{FACT_PREFIX[m]}.geom.heavy.mae') for m in ('pcm', 'gas', 'periodic')}
    handles = [Line2D([], [], ls='none', ms=ms[m], label=f'{S.MODEL_LABELS[m]}, MAE {maes[m]} Å', **mk[m],
                      **style_of[m]) for m in ('pcm', 'gas', 'periodic')]
    handles.append(Patch(color='0.90', label='±0.02 Å'))
    ax.legend(handles=handles, loc='upper left', fontsize=style.TICK_PT, handletextpad=0.4, borderaxespad=0.5,
              labelspacing=0.35)
    S.panel_letter(fig, 0.4, 80.8, 'a', 'Parity of the 15 heavy-atom bonds')
    # ------------------------------------------------------------------ (b) deviations by class
    rows = []                     # (kind, text, bond)
    for cname, bonds in CLASSES:
        rows.append(('head', cname, None))
        for b in bonds:
            rows.append(('bond', pretty(b), b))
    n = len(rows)
    axb = S.axes_mm(fig, 118.0, 11.0, 55.0, 68.0)
    ypos = {k: n - 1 - k for k in range(n)}
    axb.set_ylim(-0.7, n - 0.3)
    axb.set_xlim(-0.065, 0.095)
    axb.axvspan(-band, band, color='0.90', lw=0, zorder=0)
    axb.axvline(0.0, color='0.35', lw=0.7, zorder=1)
    yt, ytl = [], []
    for k, (kind, text, b) in enumerate(rows):
        y = ypos[k]
        if kind == 'head':
            axb.text(-0.013, y, text, transform=axb.get_yaxis_transform(), ha='right', va='center', fontsize=fs,
                     color='0.2', fontweight='bold')
            if k:
                axb.axhline(y + 0.5, color='0.75', lw=0.4, zorder=0)
            continue
        yt.append(y)
        ytl.append(text)
        axb.axhline(y, color='0.93', lw=0.4, zorder=0)
        for m, off in (('pcm', 0.0), ('gas', 0.0), ('periodic', 0.0)):
            d = tab[b][m] - tab[b]['xray']
            axb.plot([d], [y + off], ls='none', ms=ms[m], zorder={'pcm': 3, 'gas': 4, 'periodic': 5}[m],
                     **mk[m], **style_of[m])
    axb.set_yticks(yt)
    axb.set_yticklabels(ytl, fontsize=fs)
    axb.tick_params(axis='y', length=0, pad=2.0)
    axb.set_xticks([-0.04, 0.0, 0.04, 0.08])
    axb.set_xticks([-0.06, -0.02, 0.02, 0.06], minor=True)
    axb.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: S.num(v, 2, sign=v > 1e-9)))
    axb.set_xlabel('Optimised − X-ray bond length (Å)')
    S.panel_letter(fig, 88.0, 80.8, 'b', 'Signed deviation by bond class')
    for b in tab:
        for m in ('xray', 'pcm', 'gas', 'periodic_i'):
            S.USED[f'bulk.bond.{m}.{b}']['printed'].append('plotted')
    S.save(fig, STEM, W_MM)


if __name__ == '__main__':
    main()
