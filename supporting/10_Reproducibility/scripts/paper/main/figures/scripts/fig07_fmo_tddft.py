"""Fig. 7 (fig07_fmo_tddft): frontier orbitals, density of states and excitations of the IEF-PCM ion pair.

    nice -n 10 tools/paper-env/bin/python paper/main/figures/scripts/fig07_fmo_tddft.py

(a, b) HOMO and LUMO of the IEF-PCM (water) ion pair, B3LYP-D3(BJ)/6-311++G(d,p), +-PARAM orbital_iso a.u.
       (archive homo.cube / lumo.cube), camera _figb.pair_view (shared with Figs. 5 and 8); energies in the
       titles from facts [fmo.pcm.homo, fmo.pcm.lumo].
(c)    Total DOS and fragment PDOS (Hirshfeld orbital composition, Multiwfn mwfn_extra/dos_curve.csv,
       FWHM facts [dos.fwhm]); HOMO / LUMO lines and the Kohn-Sham gap [fmo.pcm.gap].
(d)    CAM-B3LYP/6-311++G(d,p) absorption in IEF-PCM (water, non-equilibrium), 20 singlets, Gaussian FWHM
       [protocol.td_fwhm]: ion pair, cation and neutral sulfanilamide (spectrum_camb3lyp_noneq.csv); sticks of the
       ion-pair states (f, right axis).  Bands labelled with facts [td.pair.cam.bright.nm, td.pair.cam.max.nm,
       td.sfa.cam.bright.nm] (sfa = neutral sulfanilamide here).
(e)    Inter-ion charge-transfer number of every ion-pair state, CAM-B3LYP against B3LYP (tddft_*_noneq.json;
       each value checked against facts bulk.td.sfa.<xc>_noneq.Snn.CT); marker area proportional to f.
"""
import csv
import json

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import _figb as B
from figkit import compose, render3d as r3, style

STEM = 'fig07_fmo_tddft'
CASE = 'sfa_pair_pcm_d3'
RAW = B.RES_RAW
LABELS = ['N1', 'O5', 'P1', 'S1', 'N2', 'C4']
FC = B.FRAG_COLOURS


def read_csv(path):
    with open(path) as f:
        r = csv.reader(f)
        head = next(r)
        rows = np.array([[float(x) for x in row] for row in r])
    return [h.strip('"') for h in head], rows


def renders(mol, view):
    iso = B.PARAM('orbital_iso')
    d = RAW / CASE
    files = {k: d / f'{k}.cube' for k in ('homo', 'lumo')}
    memo = {}

    def scene(k):
        if k not in memo:
            memo[k] = r3.orbital_scene(files[k], mol, iso=iso, opacity=0.8)
        return memo[k]
    fitted = B.cached('fig07_fit', lambda: view.fit(scene('homo'), scene('lumo'), mol, pad=0.3),
                      files=list(files.values()), extra=(repr(view.as_dict()), iso))
    pxA = r3.px_per_angstrom(fitted, (175 - 2 * compose.OUTER_MM - 2.5) / 2)
    rs = [B.render_cached(f'fig07_{k}', lambda k=k: scene(k), fitted, pxA, files=[files[k]], extra=(iso, 'op0.8'),
                          crop=False) for k in ('homo', 'lumo')]
    return r3.crop_together(rs), fitted


def main():
    style.apply()
    mol = B.molecule(CASE)
    view, info = B.pair_view()
    (rh, rl), fitted = renders(mol, view)
    homo, lumo = B.fmt('fmo.pcm.homo', 2), B.fmt('fmo.pcm.lumo', 2)
    band_mm = 58.0
    fig, axes = compose.grid([rh, rl], ncols=2, width=style.DOUBLE,
                             titles=[f'HOMO, {homo} eV, ' + B.fmt('dos.pcm.HOMO.cation', 0) + ' % on the cation',
                                     f'LUMO, {lumo} eV, ' + B.fmt('dos.pcm.LUMO.cation', 0) + ' % on the cation'],
                             extra_bottom_mm=band_mm)
    for ax, name in zip(axes, ('HOMO', 'LUMO')):
        placed, dropped = compose.label_atoms(ax, LABELS, surface='soft')
        B.log(f'{name}: labels {[p.text for p in placed]} dropped {[d.text for d in dropped]}')

    W, H = fig.get_size_inches()
    MM = style.MM
    x0, y0, bw, bh = fig._figkit_grid['bottom_band']
    head = 4.6 * MM
    bottom = 9.0 * MM
    ax_h = bh - head - bottom
    gaps = [11.0 * MM, 13.0 * MM, 12.5 * MM]                      # left margin of each plot (tick labels)
    right_pad = [3.0 * MM, 9.5 * MM, 1.5 * MM]                   # (d) has a right axis
    widths_rel = np.array([1.0, 1.15, 0.95])
    usable = bw - sum(gaps) - sum(right_pad)
    ws = usable * widths_rel / widths_rel.sum()
    xs = []
    x = x0
    for g, w, rp in zip(gaps, ws, right_pad):
        xs.append(x + g)
        x += g + w + rp
    axc, axd, axe = [fig.add_axes([xx / W, (y0 + bottom) / H, w / W, ax_h / H]) for xx, w in zip(xs, ws)]
    for ax, letter, title, xx in zip((axc, axd, axe), 'cde',
                                     ('DOS / PDOS, B3LYP-D3(BJ)', 'CAM-B3LYP absorption', 'Inter-ion CT character'),
                                     xs):
        compose._header(fig, (xx - gaps[0] + 0.6 * MM) / W if letter == 'c' else (xx - 9.5 * MM) / W,
                        (y0 + bh - 0.2 * MM) / H, f'({letter})', title)

    # ---------------- (c) DOS / PDOS
    head_c, dos = read_csv(RAW / CASE / 'mwfn_extra' / 'dos_curve.csv')
    E = dos[:, 0]
    lo, hi = B.PARAM('dos_window')
    m = (E >= lo) & (E <= hi)
    tdos, pc, pa = dos[m, 1], dos[m, 2], dos[m, 3]
    assert 'cation' in head_c[2] and 'phosphate' in head_c[3], head_c
    scale = 1.0 / tdos.max()
    eh, el = B.F('fmo.pcm.homo'), B.F('fmo.pcm.lumo')
    axc.fill_between(E[m], pc * scale, color=FC['cation'], alpha=0.25, lw=0, zorder=1)
    axc.fill_between(E[m], pa * scale, color=FC['anion'], alpha=0.25, lw=0, zorder=1)
    axc.plot(E[m], tdos * scale, color='0.1', lw=0.9, label='total', zorder=3)
    axc.plot(E[m], pc * scale, color=FC['cation'], lw=0.8, label='cation', zorder=2)
    axc.plot(E[m], pa * scale, color=FC['anion'], lw=0.8, label=r'H$_2$PO$_4^-$', zorder=2)
    for e, name in ((eh, 'HOMO'), (el, 'LUMO')):
        axc.axvline(e, color='0.35', lw=0.6, ls=(0, (2.5, 1.5)), zorder=2)
    axc.set_xlim(lo, hi)
    axc.set_ylim(0, 1.32)
    yy = 0.93
    axc.annotate('', xy=(el, yy), xytext=(eh, yy), arrowprops=dict(arrowstyle='<->', lw=0.6, color='0.2',
                                                                  shrinkA=0, shrinkB=0))
    axc.text((eh + el) / 2, yy + 0.03, B.fmt('fmo.pcm.gap', 2) + ' eV', ha='center', va='bottom',
             fontsize=style.TICK_PT)
    axc.text(eh + 0.15, 1.27, 'HOMO', ha='left', va='top', fontsize=style.SMALL_PT, color='0.25')
    axc.text(el - 0.15, 1.27, 'LUMO', ha='right', va='top', fontsize=style.SMALL_PT, color='0.25')
    axc.set_xlabel('Orbital energy (eV)')
    axc.set_ylabel('DOS (normalised)')
    axc.set_yticks([0, 0.5, 1.0])
    axc.legend(loc='lower center', bbox_to_anchor=((0.5 * (eh + el) - lo) / (hi - lo), 0.12),
               fontsize=style.TICK_PT, handlelength=1.3, borderaxespad=0.0, labelspacing=0.25)

    # ---------------- (d) spectra
    specs = [('pair', CASE, 'ion pair', '-'), ('cation', 'cation_pcm_d3', 'cation', (0, (4, 1.6))),
             ('sulfanilamide', 'sulfanilamide_pcm_d3', 'sulfanilamide', (0, (1.3, 1.1)))]
    for key, case, lab, ls in specs:
        h, s = read_csv(RAW / case / 'spectrum_camb3lyp_noneq.csv')
        assert h[1] == 'wavelength_nm' and h[3].startswith('epsilon'), h
        axd.plot(s[:, 1], s[:, 3] / 1e4, color=FC[key], lw=1.0, ls=ls, label=lab)
    td = json.loads((RAW / CASE / 'tddft_camb3lyp_noneq.json').read_text())
    axd2 = axd.twinx()
    for st in td['states']:
        axd2.vlines(st['lambda_nm'], 0, st['f'], color='0.35', lw=0.7)
    xmin, xmax = 170, 290
    axd.set_xlim(xmin, xmax)
    axd.set_ylim(0, 13.5)
    axd2.set_ylim(0, 1.35)
    axd.set_xlabel('Wavelength (nm)')
    axd.set_ylabel(r'$\varepsilon$ (10$^4$ L mol$^{-1}$ cm$^{-1}$)')
    axd2.set_ylabel('Oscillator strength (ion pair)', fontsize=style.TICK_PT + 0.5)
    axd2.tick_params(labelsize=style.TICK_PT)
    # band labels (facts)
    b_pair = B.F('td.pair.cam.bright.nm')
    b_max = B.F('td.pair.cam.max.nm')
    b_sfa = B.F('td.sfa.cam.bright.nm')
    f_pair = B.F('td.pair.cam.bright.f')
    f_max = B.F('td.pair.cam.max.f')
    axd2.annotate(B.fmt('td.pair.cam.bright.nm', 0) + ' nm', xy=(b_pair, f_pair), xytext=(b_pair + 7, 0.52),
                  fontsize=style.SMALL_PT, ha='left', va='bottom',
                  arrowprops=dict(arrowstyle='-', lw=0.5, color='0.3', shrinkA=0.5, shrinkB=1.0))
    axd2.annotate(B.fmt('td.pair.cam.max.nm', 0) + ' nm', xy=(b_max, f_max), xytext=(b_max + 6, 1.18),
                  fontsize=style.SMALL_PT, ha='left', va='center',
                  arrowprops=dict(arrowstyle='-', lw=0.5, color='0.3', shrinkA=0.5, shrinkB=1.0))
    # sulfanilamide band maximum on its own curve
    h, s = read_csv(RAW / 'sulfanilamide_pcm_d3' / 'spectrum_camb3lyp_noneq.csv')
    k = np.argmin(np.abs(s[:, 1] - b_sfa))
    axd.annotate(B.fmt('td.sfa.cam.bright.nm', 0) + ' nm', xy=(s[k, 1], s[k, 3] / 1e4),
                 xytext=(s[k, 1] + 9, s[k, 3] / 1e4 + 2.3), fontsize=style.SMALL_PT, ha='left', va='bottom',
                 color=FC['sulfanilamide'],
                 arrowprops=dict(arrowstyle='-', lw=0.5, color=FC['sulfanilamide'], shrinkA=0.5, shrinkB=1.0))
    axd.legend(loc='upper right', fontsize=style.TICK_PT, handlelength=2.0, borderaxespad=0.3,
               labelspacing=0.25)

    # ---------------- (e) CT number, CAM-B3LYP vs B3LYP
    series = [('camb3lyp', 'CAM-B3LYP', 'o', '0.1', '0.1'), ('b3lyp', 'B3LYP', 'o', 'none', '#C0392B')]
    for xc, lab, mk, mfc, mec in series:
        dd = json.loads((RAW / CASE / f'tddft_{xc}_noneq.json').read_text())
        lam = np.array([st['lambda_nm'] for st in dd['states']])
        ct = np.array([st['ct']['CT_number'] for st in dd['states']])
        f = np.array([st['f'] for st in dd['states']])
        for st in dd['states']:                               # the plotted values are the dossier's values
            kk = f"bulk.td.sfa.{xc}_noneq.S{st['state']:02d}"
            assert abs(B.facts()[kk + '.CT']['value'] - st['ct']['CT_number']) < 1e-6, kk
            assert abs(B.facts()[kk + '.f']['value'] - st['f']) < 1e-6, kk
        size = 6 + 90 * f
        if mfc == 'none':
            axe.scatter(lam, ct, s=size, facecolors='none', edgecolors=mec, linewidths=0.8, label=lab, zorder=3)
        else:
            axe.scatter(lam, ct, s=size, color=mfc, edgecolors='white', linewidths=0.3, label=lab, zorder=4)
    axe.set_xlim(152, 245)
    axe.set_ylim(-0.04, 1.08)
    axe.set_xlabel('Wavelength (nm)')
    axe.set_ylabel('Inter-ion CT number')
    axe.axhline(0.5, color='0.6', lw=0.5, ls=(0, (2, 2)), zorder=1)
    axe.legend(loc='lower right', fontsize=style.TICK_PT, handletextpad=0.2, borderaxespad=0.3,
               labelspacing=0.3, markerscale=1.0, bbox_to_anchor=(1.0, 0.1))
    B.save(fig, STEM, extra=dict(camera=fitted.as_dict()))


if __name__ == '__main__':
    main()
