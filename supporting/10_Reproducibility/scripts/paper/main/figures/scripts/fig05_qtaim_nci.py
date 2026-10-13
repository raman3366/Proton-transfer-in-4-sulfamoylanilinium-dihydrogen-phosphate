"""Fig. 5 (fig05_qtaim_nci): inter-ion interactions of the IEF-PCM ion pair, B3LYP-D3(BJ)/6-311++G(d,p).

    nice -n 10 tools/paper-env/bin/python paper/main/figures/scripts/fig05_qtaim_nci.py

(a) IGMH: isosurface of the inter-fragment delta-g (cation | H2PO4-) at PARAM igmh_iso, coloured by
    sign(lambda2)rho over PARAM igmh_range with green at 0 (two slopes); Multiwfn cubes
    mwfn_extra/igmh_dg_inter.cube and igmh_sl2r.cube of the archive.
(b) QTAIM molecular graph in the same camera: atoms, covalent bonds, the two inter-ion bond paths (Multiwfn
    paths.txt, dotted), every (3,-1) bond critical point (small dark-red dots) and the (3,+1) ring critical points
    (open squares); the two inter-ion BCPs carry rho and the Laplacian (a.u.) from facts.json
    [qtaim.pcm.N1H1A_O5.*, qtaim.pcm.C6H6_O4.*].  CPs and paths are in Bohr in the Multiwfn files; they are
    checked against the nuclei of opt.xyz before drawing.
One camera for both panels (_figb.pair_view, shared with Figs. 7 and 8).
"""
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import _figb as B
from figkit import compose, draw2d, render3d as r3, style
from figkit.labels import LabelPlacer

STEM = 'fig05_qtaim_nci'
CASE = 'sfa_pair_pcm_d3'
BOHR = r3.BOHR
MWFN = B.RES_RAW / CASE / 'mwfn_extra'
LABELS_A = ['N1', 'H1A', 'O5', 'C6', 'H6', 'O4', 'P1', 'S1', 'C1']
PANEL_MM = 84.0


def read_cps(path):
    """CPs.txt: index, x, y, z (Bohr), type (1 = (3,-3), 2 = (3,-1), 3 = (3,+1), 4 = (3,+3)) -> dict."""
    lines = path.read_text().split('\n')
    n = int(lines[0].split()[0])
    out = {}
    for ln in lines[1:1 + n]:
        k, x, y, z, t = ln.split()
        out[int(k)] = (np.array([float(x), float(y), float(z)]) * BOHR, int(t))
    return out


def read_paths(path):
    """paths.txt: list of (N, 3) arrays in Angstrom (each bond path is written as two halves, BCP -> nucleus)."""
    txt = path.read_text()
    blocks = re.split(r'Path index:\s*\d+', txt)[1:]
    out = []
    for b in blocks:
        rows = [r.split() for r in b.strip().split('\n')]
        n = int(rows[0][0])
        P = np.array([[float(v) for v in r] for r in rows[1:1 + n]]) * BOHR
        out.append(P)
    return out


def qtaim_data(mol):
    cps = read_cps(MWFN / 'CPs.txt')
    paths = read_paths(MWFN / 'paths.txt')
    nuc = np.array([p for p, t in cps.values() if t == 1])
    # every nucleus of the Multiwfn run coincides with an atom of opt.xyz (same frame; CP order differs)
    assert len(nuc) == len(mol.coords)
    # (the (3,-3) CP of an H atom may sit a few hundredths of an Angstrom off its nucleus; heavy atoms < 2e-3)
    for q in nuc:
        d = np.linalg.norm(mol.coords - q, axis=1)
        k = int(np.argmin(d))
        assert d[k] < (0.05 if mol.symbols[k] == 'H' else 2e-3), (mol.labels[k], d[k])
    frag = np.zeros(len(mol.symbols), int)
    for fi, idx in enumerate(mol.fragments.values()):
        frag[idx] = fi
    inter = []
    for P in paths:
        a0 = int(np.argmin(np.linalg.norm(mol.coords - P[-1], axis=1)))
        inter.append(a0)
    # pair the two halves of each path by their common start (the BCP)
    halves = {}
    for P, a in zip(paths, inter):
        key = tuple(np.round(P[0], 3))
        halves.setdefault(key, []).append((P, a))
    bond_paths = []
    for key, hs in halves.items():
        assert len(hs) == 2, (key, len(hs))
        (P1, a1), (P2, a2) = hs
        full = np.vstack([P1[::-1], P2[1:]])
        bond_paths.append(dict(bcp=np.array(key), atoms=(a1, a2), xyz=full, inter=frag[a1] != frag[a2]))
    return cps, bond_paths


def panel_igmh(mol, view):
    iso = B.PARAM('igmh_iso')
    vmin, vmax = B.PARAM('igmh_range')
    files = [MWFN / 'igmh_dg_inter.cube', MWFN / 'igmh_sl2r.cube']
    memo = {}

    def scene():
        if 'sc' not in memo:
            sc = r3.nci_scene(*files, mol, iso=iso, vrange=(vmin, vmax))
            st = B.smooth_sl2r(sc, vmin, vmax)
            B.log(f'IGMH surface cleaned: {st}')
            sc.meta.update(sl2r_min=st['sl2r_min'], sl2r_max=st['sl2r_max'])
            memo['sc'] = sc
        return memo['sc']
    fitted = B.cached('fig05_fit', lambda: view.fit(scene(), mol, pad=0.35), files=files,
                      extra=(repr(view.as_dict()), iso, 'smooth'))
    pxA = r3.px_per_angstrom(fitted, PANEL_MM)
    r = B.render_cached('fig05_igmh', scene, fitted, pxA, files=files, extra=(iso, vmin, vmax, 'smooth3-taubin20'))
    return r, fitted, B.two_slope_scale(vmin, vmax)


def main():
    style.apply()
    mol = B.molecule(CASE)
    view, info = B.pair_view()
    B.log(f'camera: {info}')
    r, fitted, scale = panel_igmh(mol, view)
    B.log(f'IGMH render {r.image.shape}, sign(l2)rho on the surface {r.meta.get("sl2r_min"):.3f} .. '
          f'{r.meta.get("sl2r_max"):.3f} a.u.')

    # ---------------- layout: (a) render | (b) QTAIM graph, colour bar under (a)
    W = style.DOUBLE
    MM = style.MM
    h_img, w_img = r.image.shape[:2]
    a_w = PANEL_MM * MM
    a_h = a_w * h_img / w_img
    head = 4.6 * MM
    cb_band = 12.5 * MM
    H = head + a_h + cb_band + 0.8 * MM
    fig = plt.figure(figsize=(W, H))
    ax_a = fig.add_axes([0.6 * MM / W, (cb_band) / H, a_w / W, a_h / H])
    ax_a.imshow(r.image, extent=(0, w_img, h_img, 0), interpolation='antialiased')
    ax_a.set_xlim(0, w_img)
    ax_a.set_ylim(h_img, 0)
    ax_a.set_axis_off()
    ax_a._figkit_render = r
    b_x0 = (0.6 + PANEL_MM + 4.0) * MM
    b_w = W - b_x0 - 0.6 * MM
    ax_b = fig.add_axes([b_x0 / W, 0.8 * MM / H, b_w / W, (H - head - 0.8 * MM) / H])
    compose._header(fig, 0.6 * MM / W, (H - 0.3 * MM) / H, '(a)',
                    'IGMH, ' + r'$\delta g^{\mathrm{inter}}$ = ' + B.num(B.PARAM('igmh_iso'), 3) + ' a.u.')
    compose._header(fig, b_x0 / W, (H - 0.3 * MM) / H, '(b)', 'QTAIM molecular graph')
    # colour bar
    L = 0.8 * a_w
    cb = compose.colourbar(fig, (0.6 * MM + (a_w - L) / 2, 7.4 * MM, L, 2.2 * MM), scale, linear=True)
    fig._figkit_grid = dict(width=W, height=H)

    # atom labels on the render (never on an isosurface patch)
    placed, dropped = compose.label_atoms(ax_a, LABELS_A, surface='hard')
    B.log(f'(a) labels placed {[p.text for p in placed]}, dropped {[d.text for d in dropped]}')

    # ---------------- (b) QTAIM molecular graph in the same camera
    cps, bps = qtaim_data(mol)
    S = fitted.screen(mol.coords)                      # u right, v up, depth towards the camera (Angstrom)
    xyz2 = S.copy()
    labels = mol.labels
    lab_idx = [i for i, s in enumerate(mol.symbols) if s != 'H'] + [labels.index(x) for x in ('H1A', 'H6')]
    # inter-ion bond paths: dotted, drawn before the molecule (the placer treats them as obstacles)
    inter = [bp for bp in bps if bp['inter']]
    names = {}
    for bp in inter:
        a1, a2 = bp['atoms']
        names[frozenset((labels[a1], labels[a2]))] = bp
    hb = {'N1H1A_O5': ('H1A', 'O5'), 'C6H6_O4': ('H6', 'O4')}
    assert set(names) == {frozenset(v) for v in hb.values()}, names.keys()
    seg_obst = []
    for bp in inter:
        P = fitted.screen(bp['xyz'])
        ax_b.plot(P[:, 0], P[:, 1], color='#8B1A1A', lw=0.9, ls=(0, (1.2, 1.3)), zorder=1.5,
                  solid_capstyle='round', dash_capstyle='round')
        for k in range(0, len(P) - 1, 3):
            seg_obst.append((P[k, :2], P[min(k + 3, len(P) - 1), :2], 0.6))
    # BCP / RCP markers: skip those hidden behind an atom ball nearer to the camera
    rad = np.array([style.BALL_RADII.get(s, 0.35) for s in mol.symbols])

    def hidden(p):
        s = fitted.screen(p)[0]
        d = np.hypot(S[:, 0] - s[0], S[:, 1] - s[1])
        return bool(((d < rad * 0.95) & (S[:, 2] > s[2] + 0.05)).any())
    bcp_pts = [(p, t) for p, t in cps.values() if t in (2, 3)]
    n_bcp = sum(1 for _, t in bcp_pts if t == 2)
    n_rcp = sum(1 for _, t in bcp_pts if t == 3)
    B.log(f'CPs: {n_bcp} (3,-1), {n_rcp} (3,+1)')
    extra_discs = []
    inter_bcp = {}
    for key, (d, a) in hb.items():
        bp = names[frozenset((d, a))]
        inter_bcp[key] = bp['bcp']
    drawing_items = []
    # rho / Laplacian labels of the two inter-ion BCPs
    for key, (h, a) in hb.items():
        s = fitted.screen(inter_bcp[key])[0]
        rho = B.fmt(f'qtaim.pcm.{key}.rho', 3)
        lap = B.fmt(f'qtaim.pcm.{key}.lap', 3)
        txt = r'$\rho$ = ' + rho + '\n' + r'$\nabla^2\rho$ = ' + lap
        drawing_items.append(dict(text=txt, anchor=s[:2], radius_pt=2.0, priority=0, leader=True,
                                  kw=dict(fontsize=style.SMALL_PT, color='#8B1A1A', linespacing=1.05),
                                  key=('bcp', key), max_ring=7))
        extra_discs.append((s[:2], 0.0, 2.2))
    drawing = draw2d.draw_molecule(ax_b, mol.symbols, xyz2, labels=labels, bonds=mol.bonds, hbonds=None,
                                   label_atoms=lab_idx, label_h=True, hbond_distances=False, margin=0.5,
                                   extra_segments=seg_obst, extra_discs=extra_discs, extra_items=drawing_items,
                                   strict=True, placer_kw=dict(n_rings=8))
    for p, t in bcp_pts:
        if hidden(p):
            continue
        s = fitted.screen(p)[0]
        if t == 2:
            is_inter = any(np.linalg.norm(p - v) < 1e-3 for v in inter_bcp.values())
            ax_b.plot(s[0], s[1], 'o', ms=2.6 if is_inter else 1.9, mfc='#C0392B', mec='#5A0E0E', mew=0.35,
                      zorder=50)
        else:
            ax_b.plot(s[0], s[1], 's', ms=2.3, mfc='white', mec='#2E7D32', mew=0.6, zorder=50)
    # legend of the markers under the drawing
    hs = [Line2D([], [], ls='none', marker='o', ms=2.6, mfc='#C0392B', mec='#5A0E0E', mew=0.35,
                 label='bond critical point'),
          Line2D([], [], ls='none', marker='s', ms=2.3, mfc='white', mec='#2E7D32', mew=0.6,
                 label='ring critical point'),
          Line2D([], [], color='#8B1A1A', lw=0.9, ls=(0, (1.2, 1.3)), label='inter-ion bond path')]
    ax_b.legend(handles=hs, loc='lower center', bbox_to_anchor=(0.5, -0.01), ncol=3, fontsize=style.TICK_PT,
                handlelength=1.6, columnspacing=1.2, handletextpad=0.4, frameon=False)
    B.log(f'(b) labels dropped: {[d.text for d in drawing.dropped]}')
    B.save(fig, STEM, extra=dict(camera=fitted.as_dict(), n_bcp=n_bcp, n_rcp=n_rcp))


if __name__ == '__main__':
    main()
