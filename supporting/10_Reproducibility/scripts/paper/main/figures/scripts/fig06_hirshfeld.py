"""Fig. 6 (fig06_hirshfeld): Hirshfeld surfaces of the two ions in the crystal (X-ray structure, neutron-normalised
X-H) mapped with d_norm, and their fingerprint plots.

    nice -n 10 tools/paper-env/bin/python paper/main/figures/scripts/fig06_hirshfeld.py

Data: archive hirshfeld_xray/ (dft/hirshfeld.py: promolecule of spherical ROHF/aug-cc-pVQZ atoms, w = 0.5
surface, CrystalExplorer d_norm definition with Bondi radii; *_surface.npz, *_fingerprint.npz, *_red_spots.csv).
(a, b) d_norm on the cation / anion surface (PARAM hs_dnorm_range, red < 0 < blue, white at 0), surface at
       PARAM hs_opacity; camera chosen so that most strong red spots (d_norm < -0.4)
       face the viewer; visible spots named by their contact (facts bulk.hs.<ion>.spotNN).
(c, d) Fingerprint plots (d_i, d_e; area per bin of PARAM hs_fp_rebin x 0.01 A, i.e. 0.02 A), full and decomposed into the three largest
       reciprocal contact types of each ion (percentages: facts hs.<ion>.HO/HH/CH/OO); the full plot in grey
       behind each decomposed one.
"""
import csv

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

import _figb as B
from figkit import callouts, compose, geometry, render3d as r3, style

STEM = 'fig06_hirshfeld'
HS = B.RES_RAW / 'hirshfeld_xray'
IONS = ('cation', 'anion')
DECOMP = {'cation': [('H...O', 'O...H', 'HO', 'H···O/O···H'), ('H...H', None, 'HH', 'H···H'),
                     ('H...C', 'C...H', 'CH', 'H···C/C···H')],
          'anion': [('H...O', 'O...H', 'HO', 'H···O/O···H'), ('H...H', None, 'HH', 'H···H'),
                    ('O...O', None, 'OO', 'O···O')]}
DNORM_CMAP = mcolors.LinearSegmentedColormap.from_list('dnorm', [(0.0, '#C0152B'), (0.5, '#FFFFFF'),
                                                                 (1.0, '#2350A8')])
FP_CMAP = mcolors.LinearSegmentedColormap.from_list('fp', ['#2B3CC8', '#1FA3D6', '#2DB43C', '#F2C12E', '#D7261E'])
SPOT_MAX = -0.4            # red spots labelled when d_norm is below this (every one is a facts bulk.hs value)


def load_surface(ion):
    d = np.load(HS / f'{ion}_surface.npz')
    return {k: d[k] for k in d.files}


def spots(ion):
    rows = list(csv.DictReader(open(HS / f'{ion}_red_spots.csv')))
    fk = B.facts()
    out = []
    for k, r in enumerate(rows):
        key = f'bulk.hs.{ion}.spot{k:02d}'
        assert abs(fk[key + '.dnorm']['value'] - float(r['min_dnorm'])) < 5e-4, (key, r)
        a, b = geometry.clean_label(r['inside_atom']), geometry.clean_label(r['outside_atom'])
        out.append(dict(key=key, dnorm=float(r['min_dnorm']), pos=np.array([float(r[c]) for c in 'xyz']),
                        text=f'{a}···{b}', strong=float(r['min_dnorm']) < SPOT_MAX))
    return out


def ion_molecule(S):
    m = S['inside_mask']
    sym = [str(s) for s in S['atom_el'][m]]
    X = S['atom_xyz'][m]
    lab = geometry.clean_labels([str(s) for s in S['atom_labels'][m]])
    return r3.Molecule(symbols=sym, coords=X, labels=lab, raw_labels=lab, bonds=geometry.detect_bonds(sym, X))


def spot_view(S, sp, mol):
    """Direction (400-point Fibonacci sphere) that shows the most strong red spots: surface normal at the spot
    (nearest vertex) . d > 0.35; ties broken towards a view across the largest heavy-atom extent."""
    V, N = S['vertices'], S['normals']
    strong = [s for s in sp if s['strong']]
    n = []
    for s in strong:
        k = int(np.argmin(np.linalg.norm(V - s['pos'], axis=1)))
        n.append(N[k] / np.linalg.norm(N[k]))
    n = np.array(n)
    k = np.arange(400) + 0.5
    phi = np.arccos(1 - 2 * k / 400)
    th = np.pi * (1 + 5 ** 0.5) * k
    D = np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]
    c, ax, _ = mol.principal_axes()
    score = ((n @ D.T) > 0.35).sum(0) + 0.4 * np.abs(D @ ax[:, 2])
    d = D[int(np.argmax(score))]
    up = ax[:, 1] if abs(ax[:, 1] @ d) < 0.9 else ax[:, 0]
    v = r3._orient(r3.View(c, d, up), mol, None, True)
    return v, dict(n_strong=len(strong), n_facing=int(((n @ d) > 0.35).sum()))


def surface_scene(S, mol):
    import pyvista as pv
    lo, hi = B.PARAM('hs_dnorm_range')
    V, F = S['vertices'], S['faces']
    surf = pv.PolyData(V, np.c_[np.full(len(F), 3), F].ravel())
    # vertex normals recomputed with one consistent outward orientation (the stored field-gradient normals shade
    # in patches)
    surf = surf.compute_normals(cell_normals=False, point_normals=True, auto_orient_normals=True,
                                consistent_normals=True, split_vertices=False)
    norm = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=lo, vmax=hi)
    assert S['d_norm'].min() >= lo and S['d_norm'].max() <= hi, (S['d_norm'].min(), S['d_norm'].max())
    surf.point_data['t'] = np.clip(norm(S['d_norm']), 0, 1)
    op = B.PARAM('hs_opacity')
    sc = r3.Scene(mol, translucent=op < 1)
    if op < 1:                                   # the ion inside is only seen through a translucent surface
        for m, k in r3.molecule_meshes(mol, hbonds=False):
            sc.add(m, **k)
    sc.add(surf, scalars='t', cmap=DNORM_CMAP, clim=(0, 1), opacity=op, smooth_shading=True, specular=0.15,
           specular_power=15, diffuse=0.75, ambient=0.3, interpolate_before_map=True, show_scalar_bar=False,
           n_colors=512)
    sc.meta.update(kind='hirshfeld', surface_points=r3._subsample(V))
    return sc


def dnorm_scale():
    lo, hi = B.PARAM('hs_dnorm_range')
    norm = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=lo, vmax=hi)
    ticks = [-0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0]
    return r3.ColourScale(DNORM_CMAP, norm, r'$d_\mathrm{norm}$', [t for t in ticks if lo <= t <= hi], '{:.1f}')


def rebin(M, f):
    n = (M.shape[0] // f) * f
    return M[:n, :n].reshape(n // f, f, n // f, f).sum(axis=(1, 3))


def fingerprint(ax, fp, keys, title, pct, lim, full_bg=True, vmax=None):
    f = int(B.PARAM('hs_fp_rebin'))
    b = fp['bins_A'][::f]
    tot = rebin(fp['total_A2'], f)
    sel = np.zeros_like(tot)
    for k in keys:
        if k is not None:
            sel = sel + rebin(fp[k], f)
    ext = (b[0], b[-1], b[0], b[-1])
    # the npz matrices are indexed [i(d_i), j(d_e)]: transpose so that d_e runs up the y axis
    if full_bg:
        bg = np.where(tot.T > 0, 1.0, np.nan)
        ax.imshow(bg, origin='lower', extent=ext, cmap=mcolors.ListedColormap(['#D4D4D4']), vmin=0, vmax=1,
                  interpolation='nearest', aspect='equal')
    data = np.where(sel.T > 0, sel.T, np.nan)
    vmax = vmax or np.nanpercentile(tot[tot > 0], 99.5)
    ax.imshow(data, origin='lower', extent=ext, cmap=FP_CMAP, norm=mcolors.PowerNorm(0.5, vmin=0, vmax=vmax),
              interpolation='nearest', aspect='equal')
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ticks = np.arange(np.ceil(lim[0] / 0.4) * 0.4, lim[1] + 1e-9, 0.4)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.tick_params(labelsize=style.TICK_PT - 0.5, length=2, pad=1.5)
    ax.text(0.04, 0.96, title, transform=ax.transAxes, ha='left', va='top', fontsize=style.SMALL_PT)
    if pct is not None:
        ax.text(0.96, 0.04, pct, transform=ax.transAxes, ha='right', va='bottom', fontsize=style.SMALL_PT)
    return vmax


def main():
    style.apply()
    S = {ion: load_surface(ion) for ion in IONS}
    SP = {ion: spots(ion) for ion in IONS}
    mols = {ion: ion_molecule(S[ion]) for ion in IONS}
    views, fitted = {}, {}
    for ion in IONS:
        views[ion], info = spot_view(S[ion], SP[ion], mols[ion])
        B.log(f'{ion}: view {info}')
        fitted[ion] = views[ion].fit(S[ion]['vertices'], pad=0.1)
    pxA = r3.px_per_angstrom([fitted['cation']], 92.0)
    rr = {}
    for ion in IONS:
        r = B.render_cached(f'fig06_{ion}', lambda ion=ion: surface_scene(S[ion], mols[ion]), fitted[ion], pxA,
                            files=[HS / f'{ion}_surface.npz'], extra=(B.PARAM('hs_dnorm_range'),
                                                                     B.PARAM('hs_opacity'), 'v2-normals'), crop=False)
        rr[ion] = r3.crop_together(r, pad=4)
    # margins for callouts (mm at print scale)
    fig0, ax0 = compose.grid([rr['cation'], rr['anion']], ncols=2, width=style.DOUBLE)
    kpx = ax0[0].get_position().width * fig0.get_figwidth() / rr['cation'].image.shape[1]
    plt.close(fig0)
    px_mm = 1.0 / (kpx * 25.4)
    marg = dict(left=14.0, right=14.0, top=6.0, bottom=6.0)
    for ion in IONS:
        rr[ion] = callouts.expand_canvas(rr[ion], **{s: v * px_mm for s, v in marg.items()})
    band = 92.0
    fig, axes = compose.grid([rr['cation'], rr['anion']], ncols=2, width=style.DOUBLE,
                             titles=['Cation, ' + r'$d_\mathrm{norm}$', r'H$_2$PO$_4^-$, $d_\mathrm{norm}$'],
                             colourbar_scale=dnorm_scale(), colourbar_len=0.4, colourbar_linear=True,
                             extra_bottom_mm=band)
    for ion, ax in zip(IONS, axes):
        rd = ax._figkit_render
        items = []
        for s in SP[ion]:
            if not s['strong']:
                continue
            xy, _ = rd.project(np.atleast_2d(s['pos']))
            if not rd.visible(np.atleast_2d(s['pos']), towards_camera=0.15, window=3)[0]:
                continue
            B.fmt(s['key'] + '.dnorm', 3)
            items.append(dict(target=tuple(xy[0]), text=s['text'], radius=0.12 * rd.px_per_A,
                              kw=dict(fontsize=style.SMALL_PT)))
        callouts.place(ax, items, min_leader_mm=0.8, max_leader_mm=20.0, name=f'red spots {ion}', strict=True)
        B.log(f'{ion}: spots labelled {[it["text"] for it in items]}')

    # ---------------- fingerprint plots in the bottom band: two rows of four
    x0, y0, w, h = fig._figkit_grid['bottom_band']
    W, H = fig.get_size_inches()
    MM = style.MM
    lim = B.PARAM('hs_fp_range')
    head = compose._header_mm() * MM
    row_h = h / 2
    left_pad, gap = 9.0 * MM, 2.2 * MM
    size = min((w - left_pad - 3 * gap - 1.0 * MM) / 4, row_h - head - 7.5 * MM)
    for r_i, ion in enumerate(IONS):
        fp = np.load(HS / f'{ion}_fingerprint.npz')
        fp = {k: fp[k] for k in fp.files}
        ytop = y0 + h - r_i * row_h - head
        letter = 'cd'[r_i]
        compose._header(fig, (x0 + 0.6 * MM) / W, (ytop + head - 0.2 * MM) / H, f'({letter})',
                        f'Fingerprint plots, {"cation" if ion == "cation" else "H$_2$PO$_4^-$"}')
        vmax = None
        panels = [(None, None, None, 'full')] + DECOMP[ion]
        for c_i, (k1, k2, fkey, title) in enumerate(panels):
            ax = fig.add_axes([(x0 + left_pad + c_i * (size + gap)) / W, (ytop - size) / H, size / W, size / H])
            if k1 is None:
                vmax = fingerprint(ax, fp, ['total_A2'], 'full', None, lim, full_bg=False)
            else:
                pct = B.fmt(f'hs.{ion}.{fkey}', 1) + ' %'
                fingerprint(ax, fp, [k1, k2], title, pct, lim, vmax=vmax)
            if c_i == 0:
                ax.set_ylabel(r'$d_\mathrm{e}$ (Å)', fontsize=style.TICK_PT, labelpad=1)
            else:
                ax.set_yticklabels([])
            ax.set_xlabel(r'$d_\mathrm{i}$ (Å)', fontsize=style.TICK_PT, labelpad=1)
    B.save(fig, STEM, extra=dict(camera={ion: v.as_dict() for ion, v in fitted.items()}))


if __name__ == '__main__':
    main()
