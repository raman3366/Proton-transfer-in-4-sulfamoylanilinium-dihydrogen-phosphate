"""SI figure (figS_hirshfeld_decomp): complete decomposition of the Hirshfeld surfaces of the two ions (X-ray
structure, neutron-normalised X-H; same data and conventions as Fig. 6).

    nice -n 10 tools/paper-env/bin/python paper/si/figures/scripts/figS_hirshfeld_decomp.py

(a) Share of each reciprocal contact type in the surface area of the cation and of H2PO4- (facts
    bulk.hs.<ion>.pct.<type>; labels for shares >= 0.1 %).
(b) Decomposed fingerprint plots of the cation, every contact type with a share >= PARAM hs_decomp_min;
(c) the same for the anion.  Full plot in grey behind each; colour = surface area per 0.02 x 0.02 A bin (power
    scale, the scale of the ion's full plot), as in Fig. 6.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3] / 'main' / 'figures' / 'scripts'))
import _figb as B  # noqa: E402
import fig06_hirshfeld as F6  # noqa: E402
from figkit import compose, style  # noqa: E402

STEM = 'figS_hirshfeld_decomp'


def contact_types(ion):
    fk = B.facts()
    pre = f'bulk.hs.{ion}.pct.'
    out = []
    for k, v in fk.items():
        if k.startswith(pre):
            t = k[len(pre):]
            parts = t.split('_')
            out.append((t, parts, v['value'], k))
    return sorted(out, key=lambda r: -r[2])


def pretty(t):
    return t.replace('...', '···').replace('_', '/')


def main():
    style.apply()
    MM = style.MM
    W = style.DOUBLE
    lim = B.PARAM('hs_fp_range')
    types = {ion: contact_types(ion) for ion in F6.IONS}
    thr = B.PARAM('hs_decomp_min')
    show = {ion: [r for r in types[ion] if r[2] >= thr] for ion in F6.IONS}
    ncol = 5
    size = (W - 9 * MM - (ncol - 1) * 2.2 * MM - 1 * MM) / ncol
    head = compose._header_mm() * MM
    rows_fp = {ion: int(np.ceil(len(show[ion]) / ncol)) for ion in F6.IONS}
    bar_h = 42 * MM
    fp_row_h = size + 8.0 * MM
    H = head + bar_h + 17 * MM + sum(head + rows_fp[i] * fp_row_h for i in F6.IONS) + 1 * MM
    fig = plt.figure(figsize=(W, H))
    # ---------------- (a) bars
    allt = []
    for ion in F6.IONS:
        for t, parts, v, k in types[ion]:
            if t not in allt:
                allt.append(t)
    y_top = H - head
    axa = fig.add_axes([12 * MM / W, (y_top - bar_h) / H, (W - 14 * MM) / W, bar_h / H])
    compose._header(fig, 0.6 * MM / W, (H - 0.2 * MM) / H, '(a)', 'Share of the Hirshfeld surface by contact type')
    xs = np.arange(len(allt))
    wbar = 0.38
    for j, (ion, lab) in enumerate((('cation', 'cation'), ('anion', 'H$_2$PO$_4^-$'))):
        vals = {t: (v, k) for t, parts, v, k in types[ion]}
        hv = [vals.get(t, (0.0, None))[0] for t in allt]
        axa.bar(xs + (j - 0.5) * wbar, hv, width=wbar, color=B.FRAG_COLOURS[ion], label=lab, zorder=2)
        for x_, t in zip(xs, allt):
            if t in vals and vals[t][0] >= 0.1:
                axa.text(x_ + (j - 0.5) * wbar, vals[t][0] + 1.0, B.fmt(vals[t][1], 1), ha='center', va='bottom',
                         fontsize=style.SMALL_PT - 0.5, rotation=90)
    axa.set_xticks(xs)
    axa.set_xticklabels([pretty(t) for t in allt], rotation=30, ha='right', fontsize=style.TICK_PT)
    axa.set_ylabel('Share of surface (%)')
    axa.set_ylim(0, 100)
    axa.set_xlim(-0.6, len(allt) - 0.4)
    axa.legend(loc='upper right', fontsize=style.TICK_PT)
    # ---------------- (b, c) fingerprints
    y = y_top - bar_h - 17 * MM
    for letter, ion in zip('bc', F6.IONS):
        fp = np.load(F6.HS / f'{ion}_fingerprint.npz')
        fp = {k: fp[k] for k in fp.files}
        compose._header(fig, 0.6 * MM / W, (y - 0.2 * MM) / H, f'({letter})',
                        'Decomposed fingerprint plots, ' + ('cation' if ion == 'cation' else 'H$_2$PO$_4^-$'))
        y -= head
        # the colour scale of the ion's full plot (as in Fig. 6)
        tot = F6.rebin(fp['total_A2'], int(B.PARAM('hs_fp_rebin')))
        vmax = np.nanpercentile(tot[tot > 0], 99.5)
        for n, (t, parts, v, k) in enumerate(show[ion]):
            r, cc = divmod(n, ncol)
            x0 = 9 * MM + cc * (size + 2.2 * MM)
            y0 = y - r * fp_row_h - size
            ax = fig.add_axes([x0 / W, y0 / H, size / W, size / H])
            a_, b_ = parts[0].split('...')
            rev = f'{b_}...{a_}'
            keys = [x for x in {parts[0], rev} if x in fp]
            assert keys, (ion, t)
            F6.fingerprint(ax, fp, keys, pretty(t), B.fmt(k, 1) + ' %', lim, vmax=vmax)
            if cc == 0:
                ax.set_ylabel(r'$d_\mathrm{e}$ (Å)', fontsize=style.TICK_PT, labelpad=1)
            else:
                ax.set_yticklabels([])
            ax.set_xlabel(r'$d_\mathrm{i}$ (Å)', fontsize=style.TICK_PT, labelpad=1)
        y -= rows_fp[ion] * fp_row_h
    B.save(fig, STEM, si=True, width_mm=W / MM)


if __name__ == '__main__':
    main()
