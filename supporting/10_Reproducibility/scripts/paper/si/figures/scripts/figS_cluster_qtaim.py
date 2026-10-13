"""SI figure (figS_cluster_qtaim): QTAIM of the nine crystal hydrogen bonds in the five-pair cluster, salt against H1A on O5.

    nice -n 10 tools/paper-env/bin/python paper/si/figures/scripts/figS_cluster_qtaim.py

Data: facts section 34 (qtaim.xclu.<bond>.{salt,neutral}.{rho,H}, .d_rho_pct), from
molecules/sfa/results/crystal_cluster_qtaim/qtaim_nine_bonds.json, re-derived from the raw Multiwfn CPprop.txt by
paper/main/data/derive/d13_referee_calcs.py (VM burst 3, referee M6).  135-atom cluster of five ion pairs cut from the
relaxed periodic salt (i) and from the periodic cell with H1A held on O5 (O5-H1A 0.98 A), B3LYP-D3(BJ)/6-31+G(d,p),
gas phase, single points; values for the contact at the central ion pair.
(a) rho at the bond critical point of each contact, salt (blue circle) and H1A on O5 (orange diamond), with the change
    in per cent (the N1-H1A...O5 row: with H1A on O5 the point is that of N1...H1A).
(b) total energy density H at the same points; H < 0 marks a partly covalent contact.
Colours: the paper's "proton on N1 / proton on O5" pair (validated with dataviz validate_palette.js: CVD dE 28.0,
normal 36.2; the orange is below 3:1 against white, so both series also differ by marker and are labelled).
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'main' / 'figures' / 'scripts'))
import _figb as B  # noqa: E402
from figkit import style  # noqa: E402

STEM = 'figS_cluster_qtaim'
BONDS = [('N1H1A_O5', 'N1–H1A···O5'), ('N1H1B_O6', 'N1–H1B···O6'), ('N1H1C_O6', 'N1–H1C···O6'), ('N2H2A_O3', 'N2–H2A···O3'),
         ('N2H2B_O2', 'N2–H2B···O2'), ('O3H3A_O1', 'O3–H3A···O1'), ('O4H4A_O5', 'O4–H4A···O5'), ('C3H3_O6', 'C3–H3···O6'),
         ('C6H6_O4', 'C6–H6···O4')]
C_SALT, C_NEU = '#1F5AA6', '#E07B00'


def dumbbell(ax, key, fmt_x):
    n = len(BONDS)
    for i, (k, lab) in enumerate(BONDS):
        y = n - 1 - i
        a, b = B.F(f'qtaim.xclu.{k}.salt.{key}'), B.F(f'qtaim.xclu.{k}.neutral.{key}')
        ax.plot([a, b], [y, y], color='0.6', lw=1.0, zorder=1)
        ax.plot(a, y, ls='none', marker='o', ms=5.0, mfc=C_SALT, mec='white', mew=0.5, zorder=3)
        ax.plot(b, y, ls='none', marker='D', ms=4.4, mfc=C_NEU, mec='white', mew=0.5, zorder=3)
        B.USED[f'qtaim.xclu.{k}.salt.{key}'] = 'plotted'
        B.USED[f'qtaim.xclu.{k}.neutral.{key}'] = 'plotted'
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_yticks(range(n))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(fmt_x))
    ax.grid(axis='x', color='0.9', lw=0.5)
    ax.set_axisbelow(True)


def main():
    style.apply()
    MM = style.MM
    W, H = style.DOUBLE, 80 * MM
    fig = plt.figure(figsize=(W, H))
    left, gap, right, bottom, top = 28 * MM, 14 * MM, 3 * MM, 11 * MM, 8 * MM
    pw = (W - left - gap - right) / 2
    axa = fig.add_axes([left / W, bottom / H, pw / W, (H - bottom - top) / H])
    axb = fig.add_axes([(left + pw + gap) / W, bottom / H, pw / W, (H - bottom - top) / H], sharey=axa)
    fs = style.SMALL_PT
    # (a) rho with the change in per cent
    dumbbell(axa, 'rho', lambda v, p: B.num(v, 2))
    axa.set_yticklabels([lab for _, lab in BONDS][::-1], fontsize=style.TICK_PT)
    axa.set_xlim(0.0, 0.092)
    axa.set_xlabel(r'$\rho$ at the bond critical point (a.u.)')
    n = len(BONDS)
    for i, (k, lab) in enumerate(BONDS):
        y = n - 1 - i
        a, b = B.F(f'qtaim.xclu.{k}.salt.rho'), B.F(f'qtaim.xclu.{k}.neutral.rho')
        axa.text(max(a, b) + 0.003, y, B.fmt(f'qtaim.xclu.{k}.d_rho_pct', 0, signed=True) + ' %' + ('$^{a}$' if k == 'N1H1A_O5' else ''),
                 ha='left', va='center', fontsize=fs, color='0.15')
    # (b) total energy density
    dumbbell(axb, 'H', lambda v, p: B.num(v, 3))
    plt.setp(axb.get_yticklabels(), visible=False)
    axb.axvline(0.0, color='0.35', lw=0.7, zorder=0)
    axb.set_xlim(-0.0145, 0.0035)
    axb.set_xticks([-0.012, -0.008, -0.004, 0.0])
    axb.set_xlabel('$H$ at the bond critical point (a.u.)')
    axb.text(-0.0005, -0.42, '$H$ < 0: partly covalent', ha='right', va='center', fontsize=fs, color='0.35')
    # checks behind the caption
    assert B.F('qtaim.xclu.N1H1B_O6.d_rho_pct') < -30 and B.F('qtaim.xclu.O4H4A_O5.d_rho_pct') < -20
    assert B.F('qtaim.xclu.salt.poincare_hopf') and B.F('qtaim.xclu.neutral.poincare_hopf')
    hl = [Line2D([], [], ls='none', marker='o', ms=5.0, mfc=C_SALT, mec='white', mew=0.5, label='salt (relaxed cell)'),
          Line2D([], [], ls='none', marker='D', ms=4.4, mfc=C_NEU, mec='white', mew=0.5, label='H1A held on O5')]
    fig.legend(handles=hl, loc='upper center', bbox_to_anchor=(0.33, 1.0), ncol=2, fontsize=style.TICK_PT, frameon=False,
               handletextpad=0.4, columnspacing=1.6)
    B.panel_letter(fig, axa, 'a', None, dx_mm=-left / MM + 0.6)
    B.panel_letter(fig, axb, 'b', None, dx_mm=-3.0)
    B.save(fig, STEM, si=True, width_mm=W / MM)


if __name__ == '__main__':
    main()
