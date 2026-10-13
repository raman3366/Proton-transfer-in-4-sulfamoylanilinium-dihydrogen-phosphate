"""SI figure (figS_elf_lol): ELF and LOL of the IEF-PCM (water) ion pair in three planes.

    nice -n 10 tools/paper-env/bin/python paper/si/figures/scripts/figS_elf_lol.py

Data: archive sfa_pair_pcm_d3/{elf,lol}_<plane>_plane.txt (Multiwfn 2026.9.20, 300 x 300 points, B3LYP-D3(BJ)/
6-311++G(d,p) density in IEF-PCM water).  Planes (Multiwfn input: three atoms): benzene ring (C1, C3, C5),
sulfonamide group (S1, N2, O1) and the R2,2(8) ring (N1, O5, O4).  Rows = planes, columns = ELF, LOL; one colour bar
(0-1).  Atoms within PARAM elf_atom_tol of the plane are drawn and labelled; the N1-H1A...O5 hydrogen bond is
dashed where its H and acceptor both lie in the plane.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'main' / 'figures' / 'scripts'))
import _figb as B  # noqa: E402
from figkit import maps, style  # noqa: E402

STEM = 'figS_elf_lol'
CASE = 'sfa_pair_pcm_d3'
PLANES = [('benzene_ring', 'benzene ring (C1, C3, C5)', ('C1', 'C3', 'C5')),
          ('sulfonamide_S1_N2_O1', 'sulfonamide (S1, N2, O1)', ('S1', 'N2', 'O1')),
          ('hbond_ring_R22_8', 'R$^2_2$(8) ring (N1, O5, O4)', ('N1', 'O5', 'O4'))]


def main():
    style.apply()
    mol = B.molecule(CASE)
    sym, X, lab = mol.symbols, mol.coords, mol.labels
    tol = B.PARAM('elf_atom_tol')
    d = B.RES_RAW / CASE
    MMm = 1.0                                    # layout in mm
    W = 175.0
    left, right, gap_x, cb_gap, cb_w, cb_lab = 10.5, 1.0, 3.0, 2.5, 3.0, 9.0
    top, bottom, gap_y = 5.5, 10.5, 8.0
    pw = (W - left - right - gap_x - cb_gap - cb_w - cb_lab) / 2
    geoms, grids = [], []
    for name, title, atoms in PLANES:
        e = maps.read_plane(d / f'elf_{name}_plane.txt')
        l_ = maps.read_plane(d / f'lol_{name}_plane.txt')
        assert np.abs(e.origin - l_.origin).max() < 1e-3 and np.abs(e.normal - l_.normal).max() < 1e-3, name
        uv, dist = e.project(X)
        for a in atoms:
            assert abs(dist[lab.index(a)]) < 0.01, (name, a, dist[lab.index(a)])
        g = maps.map_geometry(e, sym, X, atom_tol=tol, rotate='auto', crop='atoms', crop_margin=1.3)
        inp = [lab[i] for i in np.nonzero(g['inplane'])[0]]
        B.log(f'{name}: {len(inp)} atoms within {tol} A: {", ".join(inp)}')
        geoms.append(g)
        grids.append((e, l_))
    ph = [pw * (g['window'][3] - g['window'][2]) / (g['window'][1] - g['window'][0]) for g in geoms]
    H = top + bottom + sum(ph) + gap_y * (len(PLANES) - 1)
    fig = plt.figure(figsize=(W * style.MM, H * style.MM))
    y = H - top
    letters = iter('abcdef')
    mappable = None
    for (name, title, atoms), g, (e, l_), h in zip(PLANES, geoms, grids, ph):
        y0 = y - h
        for c_, (kind, grid) in enumerate((('ELF', e), ('LOL', l_))):
            x0 = left + c_ * (pw + gap_x)
            ax = fig.add_axes([x0 / W, y0 / H, pw / W, h / H])
            last = name == PLANES[-1][0]
            pm = maps.plot_plane_map(ax, grid, sym, X, labels=lab, kind=kind, geom=g, atom_tol=tol, colorbar=False,
                                     hbonds=mol.hbonds, strict=True, xlabel='x (Å)' if last else None,
                                     ylabel='y (Å)' if c_ == 0 else None, line_smooth=1.0, line_min_extent=0.15)
            if not last:
                ax.tick_params(labelbottom=False)
            if c_ == 1:
                ax.tick_params(labelleft=False)
            mappable = pm.mappable
            ax.text(0, 1 + 1.3 / h, f'({next(letters)})', transform=ax.transAxes, ha='left', va='baseline',
                    fontsize=style.PANEL_PT, fontweight='bold')
            ax.text(0.075, 1 + 1.3 / h, f'{kind}, {title}', transform=ax.transAxes, ha='left', va='baseline',
                    fontsize=style.BASE_PT)
        y = y0 - gap_y
    cax = fig.add_axes([(left + 2 * pw + gap_x + cb_gap) / W, bottom / H, cb_w / W, (H - top - bottom) / H])
    cb = fig.colorbar(mappable, cax=cax, ticks=np.linspace(0, 1, 6))
    cb.set_label('ELF, LOL')
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=style.TICK_PT)
    B.save(fig, STEM, si=True, width_mm=W)


if __name__ == '__main__':
    main()
