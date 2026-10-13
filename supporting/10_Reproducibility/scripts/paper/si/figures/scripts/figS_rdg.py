"""SI figure (figS_rdg): RDG against sign(lambda2)rho for the ion pair in IEF-PCM (water) and in the gas phase.

    nice -n 10 tools/paper-env/bin/python paper/si/figures/scripts/figS_rdg.py

Data: archive <case>/rdg_scatter.npz (Multiwfn grid points, RDG cut-off raised so that points with rho < 0.1 a.u. are
kept; facts nci.pcm.hb_spike) and nci_spikes.csv; spikes named by figkit rdg.identify_spikes from the grid-point
coordinates in rdg_scatter.txt (nearest atom pair; hydrogen bonds by name).  Colour: sign(lambda2)rho on the same
two-slope scale as the IGMH isosurface of Fig. 5 (PARAM igmh_range, green at 0).  In the gas phase the proton sits on
O5, so the strongest spike is O5-H1A...N1.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'main' / 'figures' / 'scripts'))
import _figb as B  # noqa: E402
from figkit import rdg, style  # noqa: E402

STEM = 'figS_rdg'
MODELS = [('sfa_pair_pcm_d3', 'IEF-PCM (water)'), ('sfa_pair_gas_d3', 'Gas phase')]


def green0_cmap(vmin, vmax):
    """render3d NCI colours with green exactly at 0 on a linear value axis (= TwoSlopeNorm + NCI_CMAP)."""
    g = (0.0 - vmin) / (vmax - vmin)
    return LinearSegmentedColormap.from_list('nci_green0', [(0.0, '#1E3CD2'), (g, '#23B43A'), (1.0, '#D7261E')],
                                             N=1024)


def main():
    style.apply()
    vmin, vmax = B.PARAM('igmh_range')
    xlim = tuple(B.PARAM('rdg_xlim'))
    ylim = (0.0, B.PARAM('rdg_ymax'))
    cmap = green0_cmap(vmin, vmax)
    MM = style.MM
    W = style.DOUBLE
    H = 74 * MM
    fig = plt.figure(figsize=(W, H))
    left, gap, right, bottom, top = 12.0 * MM, 6.0 * MM, 3.5 * MM, 21.0 * MM, 9.5 * MM
    pw = (W - left - right - gap) / 2
    B.F('nci.pcm.hb_spike')
    for k, (case, name) in enumerate(MODELS):
        ax = fig.add_axes([(left + k * (pw + gap)) / W, bottom / H, pw / W, (H - bottom - top) / H])
        d = B.RES_RAW / case
        data = rdg.load_rdg(d / 'rdg_scatter.npz', xlim=xlim, ylim=ylim)
        mol = B.molecule(case)
        spikes = rdg.read_spikes(d / 'nci_spikes.csv')
        spikes = B.cached(f'figS_rdg_spikes_{case}', lambda: rdg.identify_spikes(
            d / 'rdg_scatter.txt', spikes, mol.symbols, mol.coords, mol.labels, fragments=mol.fragments),
            files=[d / 'rdg_scatter.txt', d / 'nci_spikes.csv'])
        for sp in spikes:
            if sp['label'] and sp['label'].startswith('6-ring'):
                sp['label'] = 'benzene ring'
        B.log(f'{case}: {len(data["x"])} points; spikes ' + '; '.join(f'{s["value"]:+.4f} {s["label"]}' for s in spikes))
        shown = sorted([sp for sp in spikes if sp['label'] and xlim[0] <= sp['value'] <= xlim[1]],
                       key=lambda sp: sp['value'])
        num = {sp['label']: i + 1 for i, sp in enumerate(shown)}
        out = rdg.plot_rdg(ax, data['x'], data['y'], xlim=xlim, ylim=ylim, vrange=(vmin, vmax), cmap=cmap,
                           spikes=spikes, ylabel='RDG' if k == 0 else None, spike_style='numbered',
                           spike_key=False, region_bounds=rdg.REGION_BOUNDS, number_rings=7, numbering=num)
        assert not out['dropped'], out['dropped']
        B.log(f'{case}: numbering {num}')
        items = sorted(((n_, lab) for lab, n_ in num.items()), key=lambda t: t[0])
        ents = [f'{n_}  {lab}' for n_, lab in items]
        key = '\n'.join('     '.join(ents[i:i + 3]) for i in range(0, len(ents), 3))
        ax.text(0.0, -(10.8 * MM) / (H - bottom - top), key, transform=ax.transAxes, ha='left', va='top',
                fontsize=style.SMALL_PT, linespacing=1.3)
        if k:
            ax.tick_params(labelleft=False)
        ax.text(0.0, 1 + 5.2 * MM / (H - bottom - top), f'({"ab"[k]})', transform=ax.transAxes, ha='left',
                va='bottom', fontsize=style.PANEL_PT, fontweight='bold')
        ax.text(0.06, 1 + 5.2 * MM / (H - bottom - top), name, transform=ax.transAxes, ha='left', va='bottom',
                fontsize=style.BASE_PT)
    B.save(fig, STEM, si=True, width_mm=W / MM)


if __name__ == '__main__':
    main()
