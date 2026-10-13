"""Fig. 8 (fig08_mep): molecular electrostatic potential of the IEF-PCM ion pair and of neutral sulfanilamide.

    nice -n 10 tools/paper-env/bin/python paper/main/figures/scripts/fig08_mep.py

(a) Ion pair, (b) neutral sulfanilamide (the species docked in Figs. 9-10), both B3LYP-D3(BJ)/6-311++G(d,p) in
    IEF-PCM (water): ESP (esp.cube) mapped on the rho = PARAM mep_iso a.u. density surface (density.cube), one
    colour range PARAM mep_range for both (two slopes, green at 0).  Surface extrema from Multiwfn
    (mep_extrema.csv): minima cyan, maxima amber beads; the ion-pair extrema with |V| >= PARAM mep_label_min_abs
    and the global extrema of sulfanilamide are labelled with their values (facts bulk.mep.NN, mep.sfa.min/max)
    and the nearest atom.  Camera: _figb.pair_view for (a) (Figs. 5 and 7), the same cation orientation carried to
    sulfanilamide by a Kabsch fit on the 11 heavy atoms for (b).
(c) Surface area per ESP bin (Multiwfn vtx.txt values on a marching-cubes triangulation, figkit
    plots.esp_surface_areas; totals checked against mep_surface.json).
"""
import numpy as np

import _figb as B
from figkit import callouts, compose, plots, render3d as r3, style

STEM = 'fig08_mep'
CASES = {'pair': 'sfa_pair_pcm_d3', 'sfa': 'sulfanilamide_pcm_d3'}
TITLES = {'pair': 'Ion pair, IEF-PCM (water)', 'sfa': 'Sulfanilamide, IEF-PCM (water)'}
MARGIN_MM = dict(left=11.0, right=11.0, top=7.0, bottom=7.0)
BAND_MM = 40.0


def extrema(key, mol):
    """Rows of mep_extrema.csv with site labels; the printed value is the dossier's (checked)."""
    rows = r3.read_mep_extrema(B.RES_RAW / CASES[key] / 'mep_extrema.csv')
    out = []
    fk = B.facts()
    for r in rows:
        k = int(np.argmin(np.linalg.norm(mol.coords - r['pos'], axis=1)))
        r['site'] = mol.labels[k]
        if key == 'pair':
            match = [kk for kk in fk if kk.startswith('bulk.mep.') and kk.count('.') == 2
                     and abs(fk[kk]['value'] - r['kcal']) < 5e-4]
            assert len(match) == 1, (r, match)
            r['fact'] = match[0]
            r['label'] = abs(r['kcal']) >= B.PARAM('mep_label_min_abs')
        else:
            if r['is_global'] and r['type'] == 'minimum':
                r['fact'] = 'mep.sfa.min'
            elif r['is_global'] and r['type'] == 'maximum':
                r['fact'] = 'mep.sfa.max'
            r['label'] = 'fact' in r
            if 'fact' in r:
                assert abs(B.F(r['fact']) - r['kcal']) < 5e-3, r
        out.append(r)
    return out


def extrema_view(mol, rows, dens_path, ref_view):
    """Camera from which the labelled extrema face the viewer: over 400 directions on a Fibonacci sphere, the
    score is sum_i w_i [n_i . d > 0.3] (n_i = outward surface normal at extremum i from -grad rho, w = 3 for a
    global extremum, 1 otherwise) + 0.5 (d . d_ref) to stay close to the camera of Figs. 5 and 7; then landscape
    with the first fragment on the left (figkit _orient)."""
    dens = r3.read_cube(dens_path)
    P, N, Wt = [], [], []
    for r in rows:
        if not r.get('label'):
            continue
        q = r3._lift(dens, r['pos'], 0.3)
        n = q - r['pos']
        P.append(r['pos'])
        N.append(n / np.linalg.norm(n))
        Wt.append(3.0 if r['is_global'] else 1.0)
    N, Wt = np.array(N), np.array(Wt)
    k = np.arange(400) + 0.5
    phi = np.arccos(1 - 2 * k / 400)
    th = np.pi * (1 + 5 ** 0.5) * k
    D = np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]
    score = ((N @ D.T) > 0.3).astype(float).T @ Wt + 0.5 * (D @ ref_view.direction)
    d = D[int(np.argmax(score))]
    v = r3.View(ref_view.focal, d, ref_view.up if abs(ref_view.up @ d) < 0.95 else ref_view.right)
    v = r3._orient(v, mol, None, True)
    return v, dict(score=float(score.max()), n_labelled=int(len(Wt)), max_possible=float(Wt.sum() + 0.5))


def main():
    style.apply()
    mols = {k: B.molecule(c) for k, c in CASES.items()}
    view, _ = B.pair_view()
    iso = B.PARAM('mep_iso')
    vr = tuple(B.PARAM('mep_range'))
    ext = {k: extrema(k, mols[k]) for k in CASES}
    ref_sfa = r3.matched_view(view, mols['pair'], mols['sfa'], atoms=range(11))
    v_sfa, info = B.cached('fig08_view_sfa', lambda: extrema_view(mols['sfa'], ext['sfa'],
                                                                   B.RES_RAW / CASES['sfa'] / 'density.cube', ref_sfa),
                           files=[B.RES_RAW / CASES['sfa'] / 'mep_extrema.csv'], extra=(repr(ref_sfa.as_dict()), 'v1'))
    B.log(f'sfa: extrema camera {info}')
    panels = {'front': ('pair', view), 'back': ('pair', view.flipped()), 'sfa': ('sfa', v_sfa)}
    memo = {}

    def scene(k):
        if k not in memo:
            d = B.RES_RAW / CASES[k]
            sc = r3.mep_scene(d / 'density.cube', d / 'esp.cube', mols[k], iso=iso, extrema=ext[k])
            r3.recolour_mep(sc, vrange=vr, centre_zero=True)
            B.log(f'{k}: surface ESP {sc.meta["vmin_surface"]:.1f} .. {sc.meta["vmax_surface"]:.1f} kcal/mol')
            assert vr[0] <= sc.meta['vmin_surface'] and sc.meta['vmax_surface'] <= vr[1]
            memo[k] = sc
        return memo[k]
    files = {k: [B.RES_RAW / c / 'density.cube', B.RES_RAW / c / 'esp.cube'] for k, c in CASES.items()}
    fitted = {p: B.cached(f'fig08_fit_{p}', lambda p=p: panels[p][1].fit(scene(panels[p][0]), pad=0.15),
                          files=files[panels[p][0]], extra=(repr(panels[p][1].as_dict()), iso))
              for p in panels}
    pxA = r3.px_per_angstrom([fitted['front'], fitted['back']], 80.0)
    rr = {}
    for p, (k, _) in panels.items():
        r = B.render_cached(f'fig08_{p}', lambda k=k: scene(k), fitted[p], pxA, files=files[k],
                            extra=(iso, vr, 'v2'), crop=False)
        rr[p] = r3.crop_together(r, pad=4)
    order = ['front', 'back', 'sfa']
    # grid once without margins to learn the print scale, then expand the canvases by millimetres
    import matplotlib.pyplot as plt
    fig0, ax0 = compose.grid([rr[p] for p in order] + [None], ncols=2, width=style.DOUBLE)
    k_in_per_px = ax0[0].get_position().width * fig0.get_figwidth() / rr['front'].image.shape[1]
    plt.close(fig0)
    px_mm = 1.0 / (k_in_per_px * 25.4)
    for p in order:
        rr[p] = callouts.expand_canvas(rr[p], **{s_: v * px_mm for s_, v in MARGIN_MM.items()})
    titles = ['Ion pair, IEF-PCM (water), front view', 'Ion pair, back view',
              'Sulfanilamide, IEF-PCM (water)']
    fig, axes = compose.grid([rr[p] for p in order] + [None], ncols=2, width=style.DOUBLE, titles=titles,
                             colourbar_scale=mep_scale(vr), colourbar_len=0.5, colourbar_linear=True)
    done = set()
    for p, ax in zip(order, axes):
        k = panels[p][0]
        rd = ax._figkit_render
        items = []
        for i, r in enumerate(ext[k]):
            if not r.get('label') or (k, i) in done:
                continue
            mk = [m for m in rd.markers if abs(m['value'] - r['kcal']) < 1e-6]
            assert len(mk) == 1, (r, len(mk))
            pos = np.atleast_2d(mk[0]['pos'])                   # the bead drawn in the render (lifted off the surface)
            xy, _ = rd.project(pos)
            vis = rd.visible(pos, towards_camera=0.19, window=3)[0]
            if not vis:
                continue
            done.add((k, i))
            val = B.fmt(r['fact'], 1, signed=r['kcal'] > 0)
            items.append(dict(target=tuple(xy[0]), text=f'{val} ({r["site"]})',
                              radius=0.19 * rd.px_per_A, kw=dict(fontsize=style.SMALL_PT)))
        out = callouts.place(ax, items, min_leader_mm=0.8, max_leader_mm=20.0, name=f"MEP extrema {p}",
                             strict=True)
        B.log(f'{p}: labelled {[it["text"] for it in items]}')
    missing = [(k, r['site'], r['kcal']) for k in CASES for i, r in enumerate(ext[k]) if r.get('label')
               and (k, i) not in done]
    assert not missing, f'labelled extrema seen in no panel: {missing}'
    # ---------------- (d) ESP area histogram in the free cell
    W, H = fig.get_size_inches()
    MM = style.MM
    bb_b = axes[1].get_position()
    bb_c = axes[2].get_position()
    left = bb_b.x0 * W + 12.0 * MM
    right = bb_b.x1 * W - 3.2 * MM
    bottom = bb_c.y0 * H + 9.5 * MM
    top = bb_c.y1 * H - 1.0 * MM
    axd = fig.add_axes([left / W, bottom / H, (right - left) / W, (top - bottom) / H])
    compose._header(fig, bb_b.x0, (bb_c.y1 * H + compose._header_mm() * MM - 0.1 * MM) / H, '(d)',
                    'Surface area by electrostatic potential')
    data = {}
    for k, c in CASES.items():
        a = plots.esp_surface_areas(B.RES_RAW / c)
        ck = a['check']
        B.log(f'{k}: ESP area check {ck}')
        assert abs(ck['total_ratio'] - 1) < 0.03, ck
        data[k] = a
    bw = B.PARAM('mep_hist_bin')
    bins = np.arange(vr[0], vr[1] + 1e-9, bw)
    for k, col, lab in (('pair', B.FRAG_COLOURS['pair'], 'ion pair'),
                        ('sfa', B.FRAG_COLOURS['sulfanilamide'], 'sulfanilamide')):
        hst, _ = np.histogram(data[k]['esp_kcal'], bins=bins, weights=data[k]['area_A2'])
        axd.stairs(hst, bins, color=col, lw=1.0, label=lab)
    axd.axvline(0, color='0.5', lw=0.5, ls=(0, (2, 2)))
    axd.set_xlim(*vr)
    axd.set_xlabel('ESP on the surface (kcal mol$^{-1}$)')
    axd.set_ylabel('Area per bin (Å$^2$)')
    axd.legend(loc='upper right', fontsize=style.TICK_PT, handlelength=1.6)
    B.save(fig, STEM, extra=dict(camera={p: v.as_dict() for p, v in fitted.items()}))


def mep_scale(vr):
    s = r3.mep_colour_scale(None, vrange=vr, centre_zero=True)
    s.ticks = [t for t in range(-80, 101, 20)]
    s.fmt = '{:.0f}'
    s.label = 'Electrostatic potential (kcal mol$^{-1}$)'
    return s


if __name__ == '__main__':
    main()
