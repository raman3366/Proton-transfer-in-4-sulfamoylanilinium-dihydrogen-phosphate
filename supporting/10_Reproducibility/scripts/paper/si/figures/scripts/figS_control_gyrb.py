"""SI figure (figS_control_gyrb): the control receptor, the ATP site of E. coli DNA gyrase B (4DUH, chain A).

    nice -n 10 tools/paper-env/bin/python paper/si/figures/scripts/figS_control_gyrb.py

(a) Redocking of the crystal inhibitor with Vina, seeds of the median top-pose RMSD: the crystal pose (grey carbons),
    the top pose in the receptor as prepared (orange) and the top pose in the receptor with the ordered crystal water
    kept (purple), in the receptor frame (no superposition), with the crystal waters that the orange pose overlaps
    (red spheres; dock/control_gyrb.py, WATER_CLASH) and the pocket residues within PARAM dock_pocket_cut of any of
    the two docked poses.  Residues that PDBFixer built (109-116) are drawn with pink carbons.
(b) The representative poses of neutral sulfanilamide (DFT geometry, Vina) in the receptor as prepared (orange) and
    with the crystal water kept (purple), beside the crystal inhibitor (grey, thin).  The two poses lie in different
    sites, so the camera is restricted to directions that see them apart (view_apart); the crystal waters under the
    orange pose are drawn, and a dashed line joins the aniline N of the purple pose to Asp73 OD1.
(c) Mean top-pose Vina score against the number of heavy atoms: every species docked into 4DUH as prepared (filled
    circles) and with the crystal water kept (filled diamonds), and the same species in the four production
    structures (open, set side by side at each ligand size).  Facts ctrl.score.*, ctrl.wat.*, dock.prod.*.
Inputs: dock/control/control_summary.json, dock/C1_4DUH/{receptor.pdb, ref_RLI.sdf, runs/}, dock/C1_4DUH_wat/
(receptor.pdb for the water positions, runs/), dock/control/poses/.
"""
import json
import sys
from pathlib import Path

import numpy as np
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'main' / 'figures' / 'scripts'))
import _figb as B  # noqa: E402
import _dock3d as D  # noqa: E402
from figkit import compose, render3d as r3, style  # noqa: E402

STEM = 'figS_control_gyrb'
C1 = 'C1_4DUH'
CW = 'C1_4DUH_wat'
CS = json.loads((B.DOCK / 'control' / 'control_summary.json').read_text())
C_XTAL, C_POSE, C_NEAR, C_BUILT, C_WATER = '#D2D2D2', '#F28E2B', '#9C5FB5', '#E7A7B8', '#D62728'
BUILT = range(109, 117)
SPECIES = [('sfa_neutral', 'sulfanilamide'), ('sfa_cation', 'cation'), ('sfa_anion', 'sulfonamidate'),
           ('paba_anion', 'pABA$^-$'), ('azm_neutral', 'AZM$^0$'), ('azm_anion', 'AZM$^-$'),
           ('smx_neutral', 'SMX$^0$'), ('smx_anion', 'SMX$^-$'), ('cfx_zwit', 'ciprofloxacin'),
           ('rli_anion', '4DUH inhibitor')]
# one label per ligand size: (heavy atoms, text, offset in points, horizontal alignment)
GROUPS = [(10, 'pABA$^-$', (-9, 10), 'right'), (11, 'sulfanilamide,\ncation, sulfonamidate', (4, 23), 'left'),
          (13, 'AZM$^0$, AZM$^-$', (9, 12), 'left'), (17, 'SMX$^0$, SMX$^-$', (9, 12), 'left'),
          (24, 'ciprofloxacin', (9, -2), 'left'), (26, '4DUH inhibitor', (8, 9), 'left')]
# the four production structures: (PDB code, marker, legend, sideways shift in heavy atoms, drawing only)
TARGETS = [('3TZF', 's', '3TZF (DHPS)', -0.42), ('7S2K', '^', '7S2K (Sul2)', -0.24), ('6RL9', 'D', '6RL9 (hCA II)', 0.24),
           ('3HS4', 'v', '3HS4 (hCA II)', 0.42)]
SHIFT_PREP, SHIFT_WAT = -0.07, 0.07      # the two series of the control receptor, likewise


def waters():
    out = []
    for ln in (B.DOCK / 'C1_4DUH_wat' / 'receptor.pdb').read_text().splitlines():
        if ln.startswith('HETATM') and ln[17:20] == 'HOH':
            out.append(np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
    return np.array(out)


def scene(rec, keys, ligs, colours, radii, spheres=None, dashes=()):
    """As _dock3d.build_scene, with the residues PDBFixer built in their own carbon colour, optional spheres and
    optional dashed lines between pairs of points (drawn like the hydrogen bonds of figkit.render3d)."""
    sc = r3.Scene(None)
    atoms = D.pocket_atoms(rec, keys, ligs)
    for group, carbon in (([a for a in atoms if a['resnum'] not in BUILT], '#ABABAB'),
                          ([a for a in atoms if a['resnum'] in BUILT], C_BUILT)):
        if group:
            sym, X, bonds = D.residue_bonds(group)
            D.add_sticks(sc, sym, X, bonds, carbon, 0.065)
    for lig, col, rad in zip(ligs, colours, radii):
        D.add_sticks(sc, lig.symbols, lig.coords, lig.bonds, col, rad)
    if spheres is not None and len(spheres):
        sc.add(r3._sphere_mesh(list(spheres), [0.34] * len(spheres), [C_WATER] * len(spheres), 30),
               scalars='rgb', rgb=True, smooth_shading=True, specular=0.4)
    segs = []
    for p0, p1 in dashes:
        u = (p1 - p0) / np.linalg.norm(p1 - p0)
        start, end = p0 + 0.30 * u, p1 - 0.30 * u
        length, on, off = np.linalg.norm(end - start), 0.20, 0.13
        k = max(1, int((length + off) // (on + off)))
        lead = 0.5 * (length - (k * (on + off) - off))
        segs += [r3._cylinder(start + u * (lead + i * (on + off)), start + u * (lead + i * (on + off) + on), 0.05, 16)
                 for i in range(k)]
    if segs:
        import pyvista as pv
        sc.add(pv.merge(segs), color='#2A2A2A', smooth_shading=True, specular=0.1)
    return sc, atoms


def view_apart(ligs, context_xyz, a, b, n=900):
    """Camera chosen like _dock3d.ligand_view (fewest pose atoms hidden behind pocket atoms, the reference ligand
    face-on where possible), over the directions that see the points a and b apart: their separation is
    foreshortened by less than 7 % (|cos| < 0.35 between the view direction and a-b).  The picture is landscape
    along a-b, with a on the left."""
    L = np.vstack([l.coords for l in ligs[1:]])            # the two sulfanilamide poses; the inhibitor is drawn thin
    C = np.asarray(context_xyz)
    ab = (b - a) / np.linalg.norm(b - a)
    ref = ligs[0].coords
    _, V = np.linalg.eigh(np.cov((ref - ref.mean(0)).T))
    nrm = V[:, 0]
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n)
    th = np.pi * (1 + 5 ** 0.5) * k
    best = None
    for d in np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]:
        if abs(ab @ d) >= 0.35:
            continue
        v = r3.View(np.vstack([l.coords for l in ligs]).mean(0), d, V[:, 2] if abs(V[:, 2] @ d) < 0.95 else V[:, 1])
        sL, sC = v.screen(L), v.screen(C)
        occ = sum(int(((np.hypot(sC[:, 0] - q[0], sC[:, 1] - q[1]) < 0.9) & (sC[:, 2] > q[2] + 0.3)).any()) for q in sL)
        score = -occ + 0.8 * abs(nrm @ d)
        if best is None or score > best[0]:
            best = (score, v, occ)
    v = best[1]
    v = r3.View(v.focal, v.direction, np.cross(v.direction, ab))
    assert (v.screen(b)[0, 0] - v.screen(a)[0, 0]) > 0
    return v, dict(occluded_pose_atoms=best[2], face=round(float(abs(nrm @ v.direction)), 2),
                   separation_kept=round(float(np.sqrt(1 - (ab @ v.direction) ** 2)), 2))


def median_row(job):
    per = CS['water_diagnostic']['near_native_mode'][job]['per_seed']
    r = np.array([p['rmsd_top_A'] for p in per])
    return per[int(np.argmin(np.abs(r - np.median(r))))]


def panel_a():
    ref = D.read_sdf(B.DOCK / C1 / 'ref_RLI.sdf')
    row = median_row(f'{C1}/runs/rli_anion__vina')
    top = D.read_pdbqt(B.DOCK / C1 / 'runs' / 'rli_anion__vina' / f"seed{row['seed']:02d}" / 'out.pdbqt', model=1)
    roww = median_row(f'{CW}/runs/rli_anion__vina')
    topw = D.read_pdbqt(B.DOCK / CW / 'runs' / 'rli_anion__vina' / f"seed{roww['seed']:02d}" / 'out.pdbqt', model=1)
    W = waters()
    under = W[np.linalg.norm(W[:, None] - top.coords[None], axis=2).min(1) < B.F('ctrl.water_clash')]
    assert len(under) == row['waters_overlapped_top'], (len(under), row)
    B.log(f"(a) as prepared: seed {row['seed']}, {row['rmsd_top_A']} A (own {D.rmsd_same_order(ref, top):.2f}), "
          f"{len(under)} waters under the pose; with water: seed {roww['seed']}, {roww['rmsd_top_A']} A "
          f"(own {D.rmsd_same_order(ref, topw):.2f})")
    title = ('Inhibitor: ' + B.fmt('ctrl.val.primary.vina.rmsd', 2) + ' Å as prepared, '
             + B.fmt('ctrl.val.water.vina.rmsd', 2) + ' Å with crystal water')
    return dict(ligs=[ref, top, topw], colours=[C_XTAL, C_POSE, C_NEAR], radii=[0.16, 0.16, 0.11], spheres=under,
                dashes=(), apart=None, title=title, tag='a')


def panel_b():
    ref = D.read_sdf(B.DOCK / C1 / 'ref_RLI.sdf')
    ligs = [ref]
    for rec in (C1, CW):
        q = CS['jobs'][f'{rec}/runs/sfa_neutral__vina__prod']
        seed = q['representative']['seed']
        ligs.append(D.read_sdf(B.DOCK / 'control' / 'poses' / rec / 'sfa_neutral__vina__top_pose_per_seed.sdf',
                               prop=('seed', seed)))
        B.log(f"(b) {rec}: sulfanilamide seed {seed}, score {q['representative']['score']}")
    prep, wat = ligs[1], ligs[2]
    W = waters()
    under = W[np.linalg.norm(W[:, None] - prep.coords[None], axis=2).min(1) < B.F('ctrl.water_clash')]
    assert len(under) == B.F('ctrl.score.sfa_neutral.vina.waters_range')[0], len(under)
    # the contact of the pose with the crystal water: aniline N to Asp73 OD1 (fact ctrl.cont.water.sfa_neutral.vina.*)
    od1 = next(np.array(a['xyz']) for a in D.read_receptor(C1) if a['resnum'] == 73 and a['name'] == 'OD1')
    nit = [x for s_, x in zip(wat.symbols, wat.coords) if s_ == 'N']
    nan_ = min(nit, key=lambda x: np.linalg.norm(x - od1))
    key = next(k for k in sorted(B.facts()) if k.startswith('ctrl.cont.water.sfa_neutral.vina.')
               and 'ASP73 OD1···N(aniline)' in B.facts()[k]['meaning'])
    assert abs(np.linalg.norm(nan_ - od1) - B.F(key)) < 0.011, (np.linalg.norm(nan_ - od1), B.F(key))
    B.log(f"(b) {len(under)} crystal waters under the pose as prepared; N(aniline)...Asp73 OD1 "
          f"{np.linalg.norm(nan_ - od1):.2f} A; centroids {np.linalg.norm(prep.coords.mean(0) - wat.coords.mean(0)):.2f} A apart")
    title = ('Sulfanilamide: ' + B.fmt('ctrl.score.sfa_neutral.vina', 2) + ' and '
             + B.fmt('ctrl.wat.sfa_neutral.vina', 2) + ' kcal mol' + r'$^{-1}$')
    return dict(ligs=ligs, colours=[C_XTAL, C_POSE, C_NEAR], radii=[0.085, 0.17, 0.17], spheres=under,
                dashes=[(nan_, od1)], apart=(wat.coords.mean(0), prep.coords.mean(0)), title=title, tag='b')


def score_panel(fig, rect):
    ax = fig.add_axes(rect)
    for pdb, mk, lab, dx in TARGETS:
        xs, ys = [], []
        for sp, _ in SPECIES:
            key = f'dock.prod.{pdb}.{sp}.vina.score'
            if key in B.facts() and f'ctrl.score.{sp}.vina.heavy' in B.facts():
                xs.append(B.F(f'ctrl.score.{sp}.vina.heavy') + dx)
                ys.append(B.F(key))
        ax.plot(xs, ys, mk, ms=3.4, mfc='none', mec='0.35', mew=0.7, ls='none', label=lab, zorder=2)
    xs, ys, es, labs = [], [], [], []
    for sp, lab in SPECIES:
        if f'ctrl.score.{sp}.vina' not in B.facts():
            continue
        xs.append(B.F(f'ctrl.score.{sp}.vina.heavy'))
        ys.append(B.F(f'ctrl.score.{sp}.vina'))
        es.append(B.F(f'ctrl.score.{sp}.vina.sd'))
        labs.append(lab)
    ax.errorbar(np.array(xs) + SHIFT_PREP, ys, yerr=es, fmt='o', ms=4.0, mfc='#2B6CB0', mec='none', ecolor='#2B6CB0',
                elinewidth=0.8, capsize=1.5, label='4DUH (gyrase B) as prepared', zorder=3)
    xw = [B.F(f'ctrl.wat.{sp}.vina.heavy') for sp, _ in SPECIES if f'ctrl.wat.{sp}.vina' in B.facts()]
    yw = [B.F(f'ctrl.wat.{sp}.vina') for sp, _ in SPECIES if f'ctrl.wat.{sp}.vina' in B.facts()]
    ew = [B.F(f'ctrl.wat.{sp}.vina.sd') for sp, _ in SPECIES if f'ctrl.wat.{sp}.vina' in B.facts()]
    ax.errorbar(np.array(xw) + SHIFT_WAT, yw, yerr=ew, fmt='D', ms=3.2, mfc=C_NEAR, mec='none', ecolor=C_NEAR,
                elinewidth=0.8, capsize=1.5, label='4DUH with the crystal water', zorder=3)
    assert len(xs) == B.F('ctrl.size.vina.n'), (len(xs), B.F('ctrl.size.vina.n'))
    k, c = np.polyfit(xs, ys, 1)
    assert abs(np.corrcoef(xs, ys)[0, 1] - B.F('ctrl.size.vina.r')) < 2e-3
    xx = np.array([min(xs) - 1, max(xs) + 1])
    ax.plot(xx, k * xx + c, color='#2B6CB0', lw=0.7, ls=(0, (3, 2)), zorder=1)
    for ha_, lab, off, al in GROUPS:
        yy = [y for x, y in zip(xs, ys) if x == ha_]
        if yy:
            ax.annotate(lab, (ha_, max(yy)), xytext=off, textcoords='offset points', fontsize=style.SMALL_PT,
                        color='#1F4E79', ha=al, va='center', linespacing=1.05)
    ax.set_xlabel('Heavy atoms of the ligand', fontsize=style.BASE_PT)
    ax.set_ylabel('Mean top-pose Vina score (kcal mol$^{-1}$)', fontsize=style.BASE_PT)
    ax.set_xlim(min(xs) - 2.5, max(xs) + 7)
    allv = ys + yw + [B.F(f'dock.prod.{pdb}.{sp}.vina.score') for pdb, _, _, _ in TARGETS for sp, _ in SPECIES
                      if f'dock.prod.{pdb}.{sp}.vina.score' in B.facts()]
    ax.set_ylim(min(allv) - 0.45, max(allv) + 0.75)
    ax.grid(color='0.9', lw=0.5, zorder=0)
    ax.legend(loc='upper right', fontsize=style.TICK_PT, frameon=False, handletextpad=0.4, borderaxespad=0.2)
    ax.text(0.02, 0.05, '$r$ = ' + B.fmt('ctrl.size.vina.r', 2) + f" over the {B.F('ctrl.size.vina.n')} species in 4DUH "
            'as prepared', transform=ax.transAxes, fontsize=style.TICK_PT, color='#1F4E79', ha='left', va='bottom')
    return ax


def main():
    style.apply()
    MM = style.MM
    rec = D.read_receptor(C1)
    P = [panel_a(), panel_b()]
    views = []
    for p in P:
        keys, _ = D.pocket(rec, p['ligs'][1:], cutoff=B.PARAM('dock_pocket_cut'))      # residues near a docked pose
        sc, atoms = scene(rec, keys, p['ligs'], p['colours'], p['radii'], p['spheres'], p['dashes'])
        ctx = np.array([a['xyz'] for a in atoms])
        v, info = D.ligand_view(p['ligs'], ctx) if p['apart'] is None else view_apart(p['ligs'], ctx, *p['apart'])
        B.log(f"({p['tag']}) {len(keys)} pocket residues, built residues shown "
              f"{sorted({a['resnum'] for a in atoms if a['resnum'] in BUILT})}; view {info}")
        p.update(scene=sc, atoms=atoms, view=v.fit(sc, pad=0.4))
        views.append(p['view'])
    pxA = r3.px_per_angstrom(views, 84.0)
    rr = [B.render_cached(f"figS_control_{p['tag']}", lambda p=p: p['scene'], p['view'], pxA,
                          files=[B.DOCK / C1 / 'receptor.pdb'], extra=('v4', p['title'])) for p in P]
    plot_h = 62.0
    fig, axes = compose.grid(rr, ncols=2, width=style.DOUBLE, titles=[p['title'] for p in P],
                             extra_bottom_mm=plot_h + 14.0)
    for p, ax in zip(P, axes):
        # obstacles that are drawn but are no atoms of the label molecule: the dashed line and the water spheres
        extra = [D.Lig(symbols=['C'] * 7, coords=np.linspace(a, b, 7), bonds=[]) for a, b in p['dashes']]
        if p['spheres'] is not None and len(p['spheres']):
            extra.append(D.Lig(symbols=['O'] * len(p['spheres']), coords=np.array(p['spheres']), bonds=[]))
        mol, anchors = D.label_molecule(p['atoms'], [], p['ligs'] + extra)
        placed, dropped = compose.label_atoms(ax, anchors, mol=mol, surface=None, hydrogens_droppable=False,
                                              placer_kw=dict(n_rings=13, n_dirs=40, leader_from=1),
                                              name=f"residue labels {p['tag']}")
        B.log(f"({p['tag']}) labels {len(placed)}, dropped {[d.text for d in dropped]}")
    x0, y0, w, h = fig._figkit_grid['bottom_band']
    W, H = fig.get_size_inches()
    hs = [Line2D([], [], color=C_XTAL, lw=3.0, label='crystal inhibitor'),
          Line2D([], [], color=C_POSE, lw=3.0, label='pose as prepared'),
          Line2D([], [], color=C_NEAR, lw=3.0, label='pose with crystal water'),
          Line2D([], [], color=C_WATER, lw=0, marker='o', ms=4.5, label='crystal water under the pose'),
          Line2D([], [], color=C_BUILT, lw=2.0, label='built residues')]
    fig.legend(handles=hs, loc='center', bbox_to_anchor=((x0 + w / 2) / W, (y0 + h - 3.0 * MM) / H), ncol=5,
               fontsize=style.TICK_PT, handlelength=1.6, columnspacing=1.2, frameon=False)
    left = 17.0 * MM
    head_y = y0 + plot_h * MM + 2.2 * MM
    compose._header(fig, 0.6 * MM / W, head_y / H, '(c)',
                    'Vina score against ligand size: the control receptor and the four production structures')
    score_panel(fig, [left / W, (y0 + 9.0 * MM) / H, (W - left - 4.0 * MM) / W, (plot_h - 12.0) * MM / H])
    B.save(fig, STEM, si=True, extra=dict(cameras={p['tag']: p['view'].as_dict() for p in P}))


if __name__ == '__main__':
    main()
