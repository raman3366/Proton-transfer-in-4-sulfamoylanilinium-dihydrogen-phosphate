"""SI figure (figS_ir): computed IR spectrum of the IEF-PCM (water) ion pair, B3LYP-D3(BJ)/6-311++G(d,p).

    nice -n 10 tools/paper-env/bin/python paper/si/figures/scripts/figS_ir.py

Data: molecules/sfa/results/sfa_pair_pcm_d3/freq.json (finite-difference Hessian, 75 modes, 0 imaginary; IR intensities
from the dipole derivatives).  Harmonic wavenumbers scaled by facts [freq.scale_factor] (Merrick 2007 value for
B3LYP/6-311+G(d,p), the closest basis in their table; checked 29 Sep 2026, facts section 9), sticks plus Lorentzian bands of FWHM
PARAM ir_fwhm.  The labelled bands and their assignments are facts freq.<band>.scaled (the scaled value printed,
0 decimals).  A harmonic single ion pair in a continuum: the 2122 cm-1 N1-H1A stretch of the very short N1-H1A...O5
bond and the free O-H stretches near 3680 cm-1 do not describe the crystal (facts caveats; caption).
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'main' / 'figures' / 'scripts'))
import _figb as B  # noqa: E402
from figkit import style  # noqa: E402

STEM = 'figS_ir'
BANDS = [('NH_Hbond', 'ν(N1–H1A)'), ('OH_2', 'ν(O–H)'), ('NH2_as', 'ν$_{as}$(NH$_2$)'), ('NH3_2', 'ν(NH$_3^+$)'),
         ('SO2_as', 'ν$_{as}$(SO$_2$)'), ('PO2_as', 'ν(P–O5/O6)'), ('CS', 'ν(C4–S1)'), ('SO2_s', 'ν$_s$(SO$_2$)'),
         ('POH_b', 'δ(P–O3–H3A)'), ('P_OH_s1', 'ν(P–OH)'), ('SN', 'ν(S1–N2)'), ('ring', 'ν(C=C)'),
         ('NH3_def1', 'δ(NH$_3^+$)'), ('CN', 'ν(C1–N1)')]


def main():
    style.apply()
    fr = json.loads((B.RES / 'sfa_pair_pcm_d3' / 'freq.json').read_text())
    assert fr['n_imaginary'] == B.F('freq.n_imag') == 0
    sf = B.F('freq.scale_factor')
    nu = np.array([m['freq_cm-1'] for m in fr['ir']]) * sf
    inten = np.array([m['intensity_km_mol'] for m in fr['ir']])
    fwhm = B.PARAM('ir_fwhm')
    x = np.linspace(400, 4000, 7201)
    g = fwhm / 2
    spec = (inten[None] * (g / np.pi) / ((x[:, None] - nu[None]) ** 2 + g ** 2)).sum(1)
    W = style.DOUBLE
    MM = style.MM
    H = 78 * MM
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([6 * MM / W, 10 * MM / H, (W - 19 * MM) / W, (H - 11 * MM) / H])
    ax.fill_between(x, spec, color='#9DB9D8', lw=0, zorder=1)
    ax.plot(x, spec, color='#1F4E8C', lw=0.8, zorder=2)
    ax2 = ax.twinx()
    ax2.vlines(nu, 0, inten, color='0.25', lw=0.6, zorder=3)
    ax.set_xlim(4000, 400)
    ymax = spec.max() * 1.62
    ax.set_ylim(0, ymax)
    ax2.set_ylim(0, inten.max() * 1.62)
    ax.set_xlabel('Wavenumber, scaled by ' + B.fmt('freq.scale_factor', 4) + ' (cm$^{-1}$)')
    ax.set_ylabel('Absorbance (Lorentzian, arb. units)')
    ax2.set_ylabel('IR intensity (km mol$^{-1}$)')
    ax.set_yticks([])
    # label x positions spread apart (min_sep cm-1) with short leaders down to their bands
    min_sep = 62.0
    rows = []
    for key, lab in BANDS:
        v = B.F(f'freq.{key}.scaled')
        k = int(np.argmin(np.abs(nu - v)))
        assert abs(nu[k] - v) < 0.01, key
        rows.append([v, v, key, lab])
    rows.sort(key=lambda r: r[0])
    xs = np.array([r[0] for r in rows])
    for _ in range(300):
        moved = False
        for i in range(len(xs) - 1):
            gap = xs[i + 1] - xs[i]
            if gap < min_sep:
                d = (min_sep - gap) / 2
                xs[i] -= d
                xs[i + 1] += d
                moved = True
        if not moved:
            break
    for (v, _, key, lab), xl in zip(rows, xs):
        top = spec[(x > min(v, xl) - 45) & (x < max(v, xl) + 45)].max()
        yl = top + 0.07 * ymax
        yp = spec[np.argmin(np.abs(x - v))] + 0.01 * ymax
        ax.plot([v, v, xl], [yp, yl - 0.03 * ymax, yl - 0.005 * ymax], color='0.4', lw=0.4, zorder=4)
        ax.text(xl, yl, B.fmt(f'freq.{key}.scaled', 0) + '  ' + lab, rotation=90, ha='center', va='bottom',
                fontsize=style.SMALL_PT, zorder=5)
    B.save(fig, STEM, si=True, width_mm=W / MM)


if __name__ == '__main__':
    main()
