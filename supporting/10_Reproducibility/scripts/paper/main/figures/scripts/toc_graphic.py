"""Graphical abstract (Elsevier, the default since 5 Oct 2026) or table-of-contents graphic (ACS): where the proton
H1A sits in the crystal and in water, and the species that reaches the target.

    TOC_FORMAT=elsevier (default)  graphical_abstract.png/.tif, 5.0 x 2.0 in (12.7 x 5.08 cm) at 600 dpi = 3000 x 1200 px:
                                   Elsevier asks for at least 1328 x 531 px (w x h) in the ratio 500:200, readable at
                                   13 x 5 cm, in TIFF, EPS, PDF or MS Office format.
    TOC_FORMAT=acs                 toc_graphic.png/.tif, the ACS box described below (8.25 x 4.45 cm).
The two formats share every drawing step; only the canvas and the positions of the three columns differ.

    cd /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate
    nice -n 10 env OMP_NUM_THREADS=2 /root/workspace/cqc-dft/tools/paper-env/bin/python \
        paper/main/figures/scripts/toc_graphic.py

ACS TOC specification (Guidelines for Table of Contents/Abstract Graphics): final size 3.25 x 1.75 in (8.25 x 4.45 cm),
colour TIF at >= 300 dpi (written here at 600 dpi, the house resolution), sans-serif text >= 6 pt (here 7 and 8 pt),
labels only and no specific results, new artwork (no panel of the manuscript is reused).  The canvas is saved at
exactly that size (no tight bounding box), so the files are 1950 x 1050 px at 600 dpi.

Left, top     Crystal: formula unit 1 of the relaxed periodic PBE-D3(BJ) cell (periodic_qe/relax_i.extxyz), H1A on
              N1 (facts proton.periodic.d_N1_H1A, proton.periodic.delta), ball and stick, among the ions of its
              molecular layer (pale thin sticks, those that reach into the tile): the lattice plane (hkl) whose normal
              is nearest to the pair's ring normal (h, k, l in -2..2, spacing >= 3 A; it is (1 2 -1), 6 deg), every
              image of the four ions of the cell with its centroid within half the plane spacing of the pair; seen
              along that normal.
Left, bottom  Water continuum: the IEF-PCM (water) neutral minimum of the relaxed proton scan
              (sfa_pair_pcm_d3/scan/B_00_anchor.xyz), H1A on O5 (facts scan.B00.d_NH, scan.B00.delta), seen along
              its own ring normal; the pale blue field with a white cavity is a schematic dielectric continuum (union
              of CAVITY_SCALE x Bondi radii, projected, soft edge), not a computed surface.
              Both pairs: one scale, S1 -> N1 pointing right and the anion below, so the proton moves left to right
              as along the energy axis; H1A enlarged on a soft spot of the tile colour; dashed N1-H1A...O5
              (O5-H1A...N1 in water).
Middle        Energy E against d(N1-H1A) (axis label "H+: N -> O"), no numbers or ticks: crystal = the constrained
              N1-H1A profile of the periodic cell, PBE0-D3(BJ) single points per formula unit relative to the relaxed
              salt (facts periodic.profile.p110..p170.PBE0_kJ at proton.periodic.p110..p170.d_N1_H1A, start at
              proton.periodic.d_N1_H1A with 0), rising at every step to periodic.profile.p170.PBE0_kJ = 50.0 kJ/mol;
              water = the IEF-PCM relaxed scan, lower envelope as in Fig. 2d (facts scan.<A00..B00>.d_NH and .Erel,
              kcal/mol x 4.184 -> kJ/mol, checked against scan.Omin_kJ and scan.barrier_kJ), whose neutral end lies
              below the salt (scan.hysteresis_1.35 > 0 picks chain B at 1.35 A).  One energy axis (kJ/mol), dotted
              line at the salt energy, monotone (PCHIP) curve through the crystal points and the cubic spline of
              Fig. 2d through the water points.
Right         Docking: neutral sulfanilamide, the species released at pH 7.4, in the hCA II site of 6RL9
              (AutoDock4Zn production run, cluster-1 representative pose of dock/results/production_contacts.json, as
              Fig. 9c), sulfonamide N on the zinc (distance checked against
              bulk.dockcont.T3_6RL9.sfa_neutral__ad4zn__prod.zn); the zinc ligands His94, His96, His119 and Thr199 in
              grey sticks, the other residues within 4 A pale; carbons in the sulfanilamide orange of Figs. 8 and 9.
Colours: crystal = periodic model green, water = IEF-PCM blue (figkit style.MODEL_COLOURS, the validated pair of the
paper's figures).  Printed text: labels only, no numbers.  The facts used (plotted or checked) are printed on stdout
and stored in the PNG text chunk "Source" and the TIFF ImageDescription, so the graphic needs no side file.
Writes only paper/main/figures/graphical_abstract.png and .tif (toc_graphic.png and .tif with TOC_FORMAT=acs; drafts:
TOC_OUT=<dir> redirects both, TOC_LAX=1 reports QA issues instead of stopping).
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True                  # leave no byte-code caches next to the scripts

import io  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
from pathlib import Path  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figb as B  # noqa: E402  (thread limits, niceness >= 10, docking paths)
import _saph as S  # noqa: E402  (facts with provenance, geometries)
import _dock3d as D  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import colors as mcolors  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

from figkit import geometry as geo, qa, render3d as r3, style  # noqa: E402

FORMAT = os.environ.get('TOC_FORMAT', 'elsevier')
assert FORMAT in ('elsevier', 'acs'), FORMAT
STEM = {'elsevier': 'graphical_abstract', 'acs': 'toc_graphic'}[FORMAT]
OUT = Path(os.environ.get('TOC_OUT', str(S.FIG)))
ACS_W_IN, ACS_H_IN = 3.25, 1.75                 # ACS TOC box (inches) = 8.25 x 4.45 cm; the layout below is given for it
W_IN, H_IN = {'elsevier': (5.0, 2.0), 'acs': (ACS_W_IN, ACS_H_IN)}[FORMAT]     # 5.0 x 2.0 in: the ratio 500:200
SX, SY = W_IN / ACS_W_IN, H_IN / ACS_H_IN       # stretch of the ACS layout to the chosen canvas
DPI = 600
W_MM, H_MM = W_IN * 25.4, H_IN * 25.4
PX_MM = DPI / 25.4                              # final pixels per mm

# colours: the paper's model colours (figkit style.MODEL_COLOURS; periodic / IEF-PCM validated for CVD)
C_XTAL = style.MODEL_COLOURS['periodic']        # '#009E73' crystal
C_WATER = style.MODEL_COLOURS['pcm']            # '#1F5AA6' water continuum
C_SFA = '#F28E2B'                               # sulfanilamide carbons, as Figures 8 and 9
C_RES = '#8A8A8A'                               # zinc-ligand residues
INK = '#222222'
FS_TAG, FS_LAB = 8.0, 7.0                       # pt at print size (ACS: preferably 8, never below 6)

# layout (mm from the lower-left corner of the canvas)
TILE_R = 1.4                                    # corner radius of the picture tiles
XT = (0.0, 36.0 * SX)                           # left column: crystal (top) and water (bottom) tiles
YT_GAP = 0.8
DOCK = (59.0 * SX, W_MM)                        # right column: docking tile
PLOT = (40.4 * SX, 56.6 * SX, 6.2 * SY, 41.0 * SY)   # x0, x1, y0, y1 of the energy axes

# rendering choices (not results)
CAVITY_SCALE = 1.0                              # x Bondi radius, schematic continuum cavity
CAVITY_EDGE = 1.0                               # A, soft edge of the cavity (white -> continuum tint)
BONDI = {'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80}
H1A_SCALE = 1.35                                # enlarged H1A sphere
GLOW_R = 1.25                                   # A, radius of the soft spot behind H1A
TAG_RES = 4.0                                   # mm kept free at the top of the pair tiles for the tag
DOCK_TURN = 330.0                               # deg about the ligand axis from the ring-face-on view (occlusion-free)
KCAL_KJ = 4.184

CAT = [S.IDX[l] for l in S.CATION]
AN = [S.IDX[l] for l in S.ANION]


def pale(c, f):
    """Colour c mixed with white (f = 0: c, f = 1: white)."""
    r = np.array(mcolors.to_rgb(c))
    return mcolors.to_hex(r + (1 - r) * f)


TINT_XTAL = pale(C_XTAL, 0.90)
TINT_WATER = pale(C_WATER, 0.80)
TINT_DOCK = '#F3F3F3'


# ------------------------------------------------------------------------------------------------ geometry
def ring_normal(X):
    """Mean normal of the benzene ring and the R2,2(8) ring of one ion pair (as Fig. 2a-c)."""
    def normal(labels):
        P = X[[S.IDX[l] for l in labels]]
        return np.linalg.svd(P - P.mean(0))[2][2]
    n1 = normal(['C1', 'C2', 'C3', 'C4', 'C5', 'C6'])
    n2 = normal(['H1A', 'N1', 'C1', 'C6', 'H6', 'O4', 'P1', 'O5'])
    if n1 @ n2 < 0:
        n2 = -n2
    return (n1 + n2) / np.linalg.norm(n1 + n2)


def pair_view(X, d0):
    """Camera along +-d0 with S1 -> N1 pointing right; the sign puts the anion below the cation."""
    sn = X[S.IDX['N1']] - X[S.IDX['S1']]
    for d in (d0, -d0):
        v = r3.View(X.mean(0), d, np.cross(d, sn))
        s = v.screen(X)
        if s[AN, 1].mean() < s[CAT, 1].mean():
            return v
    raise AssertionError('no view with the anion below the cation')


def layer_plane(X, lat):
    """Lattice plane (h k l), h, k, l in -2..2 with spacing >= 3 A, whose normal is nearest to the ring normal."""
    z = ring_normal(X)
    rec = np.linalg.inv(lat).T
    best = None
    for hkl in np.ndindex(5, 5, 5):
        h = np.array(hkl) - 2
        if not h.any():
            continue
        g = h @ rec
        dsp = 1.0 / np.linalg.norm(g)
        if dsp < 3.0:
            continue
        ang = np.degrees(np.arccos(min(1.0, abs(g @ z) / np.linalg.norm(g))))
        key = (round(ang, 6), tuple(-h) if tuple(h) < tuple(-h) else tuple(h))
        if best is None or key < best[0]:
            best = (key, h, g / np.linalg.norm(g), dsp)
    (ang, hkl), h, n, dsp = best
    return tuple(int(x) for x in hkl), (n if n @ z > 0 else -n), dsp, ang


def crystal_layer():
    """Central pair (formula unit 1 of relax_i, whole) and the ions of its molecular layer."""
    P1, lat = S.periodic_pair('relax_i', fu=1)
    P2, _ = S.periodic_pair('relax_i', fu=2)
    assert abs(S.dist(P1, 'N1', 'H1A') - S.F('proton.periodic.d_N1_H1A')) < 1e-4
    assert abs(S.delta(P1) - S.F('proton.periodic.delta')) < 1e-4 and S.delta(P1) < 0      # H1A on N1
    hkl, n, dsp, ang = layer_plane(P1, lat)
    c0 = P1.mean(0)
    ions = []
    for k, (Fr, idx) in enumerate(((P1[CAT], CAT), (P1[AN], AN), (P2[CAT], CAT), (P2[AN], AN))):
        sym = [S.SYMBOLS[i] for i in idx]
        for m in np.ndindex(9, 9, 9):
            t = (np.array(m) - 4) @ lat
            if k < 2 and not t.any():
                continue                                    # the central pair itself
            Y = Fr + t
            if abs((Y.mean(0) - c0) @ n) <= dsp / 2:
                ions.append((sym, Y))
    return P1, ions, dict(hkl=hkl, spacing=round(dsp, 3), angle_to_ring_normal=round(ang, 1)), n


def water_pair():
    """IEF-PCM neutral minimum of the proton scan (H1A on O5)."""
    X = S.scan_point('B_00_anchor')
    assert abs(S.dist(X, 'N1', 'H1A') - S.F('scan.B00.d_NH')) < 6e-5
    assert abs(S.delta(X) - S.F('scan.B00.delta')) < 1e-3 and S.delta(X) > 0        # H1A nearer O5
    return X


# ------------------------------------------------------------------------------------------------ 3D scenes
def ball_stick(sc, sym, X, carbon=None, atom_scale=0.95, bond_radius=0.09, big=(), res=28):
    """Ball and stick in the house element colours (optional carbon colour), `big` atoms enlarged by H1A_SCALE."""
    import pyvista as pv
    X = np.asarray(X, float)
    rad = np.array([style.BALL_RADII[s] * atom_scale * (H1A_SCALE if i in big else 1.0) for i, s in enumerate(sym)])
    col = [carbon if (s == 'C' and carbon) else style.ELEMENT_COLOURS[s] for s in sym]
    sc.add(r3._sphere_mesh(X, rad, col, res), scalars='rgb', rgb=True, smooth_shading=True, specular=0.35,
           specular_power=30, diffuse=0.8, ambient=0.18)
    cyl = []
    for i, j in geo.detect_bonds(sym, X):
        u = r3._unit(X[j] - X[i])
        L = np.linalg.norm(X[j] - X[i])
        m = X[i] + u * (rad[i] + 0.5 * (L - rad[i] - rad[j]))
        for a, b, c in ((X[i] + u * 0.5 * rad[i], m, col[i]), (m, X[j] - u * 0.5 * rad[j], col[j])):
            cy = r3._cylinder(a, b, bond_radius, 20)
            cy.point_data['rgb'] = np.tile((np.array(mcolors.to_rgb(c)) * 255).astype(np.uint8), (cy.n_points, 1))
            cyl.append(cy)
    sc.add(pv.merge(cyl), scalars='rgb', rgb=True, smooth_shading=True, specular=0.25, specular_power=20,
           diffuse=0.8, ambient=0.2)
    return rad


def thin_sticks(sc, sym, X, carbon, radius, fade):
    """Pale thin sticks (context: lattice neighbours, pocket residues)."""
    X = np.asarray(X, float)
    cols = [pale(carbon if s == 'C' else style.ELEMENT_COLOURS[s], fade) for s in sym]
    bonds = geo.detect_bonds(sym, X)
    if bonds:
        sc.add(D._stick_mesh(X, bonds, cols, radius), scalars='rgb', rgb=True, smooth_shading=True, specular=0.15,
               specular_power=20, diffuse=0.8, ambient=0.3)


def dashed(sc, a, b, gap_a, gap_b, colour='#2A2A2A', radius=0.04, dash=(0.17, 0.11)):
    """Dashed hydrogen bond a -> b, shortened by the ball radii."""
    import pyvista as pv
    u = r3._unit(b - a)
    s, e = a + u * gap_a, b - u * gap_b
    L = np.linalg.norm(e - s)
    per = dash[0] + dash[1]
    n = max(1, int((L + dash[1]) // per))
    off = 0.5 * (L - (n * per - dash[1]))
    parts = [r3._cylinder(s + u * (off + k * per), s + u * (off + k * per + dash[0]), radius, 16) for k in range(n)]
    sc.add(pv.merge(parts), color=colour, smooth_shading=True, specular=0.1)


def pair_scene(X, acceptor, context=()):
    """Ion pair in ball and stick (H1A enlarged), the H1A...acceptor hydrogen bond dashed, context ions pale."""
    sc = r3.Scene(None)
    for sym, Y in context:
        thin_sticks(sc, sym, Y, '#4A4A4A', 0.065, 0.42)
    h = S.IDX['H1A']
    rad = ball_stick(sc, list(S.SYMBOLS), X, big={h})
    a = S.IDX[acceptor]
    dashed(sc, X[h], X[a], rad[h] + 0.05, rad[a] + 0.05)
    return sc


def render_tile(sc, view, w_mm, h_mm, px_per_A):
    """Render exactly the tile: view width/height = tile size / scale, final pixel density PX_MM."""
    from dataclasses import replace
    v = replace(view, width=w_mm * PX_MM / px_per_A, height=h_mm * PX_MM / px_per_A)
    r = r3.render(sc, v, px_per_A=px_per_A, crop=False)
    assert r.image.shape[:2] == (int(round(v.height * px_per_A)), int(round(v.width * px_per_A))), r.image.shape
    return r


# ------------------------------------------------------------------------------------------------ energy profiles
def crystal_profile():
    prof = ['p110', 'p125', 'p140', 'p155', 'p170']
    d = [S.F('proton.periodic.d_N1_H1A')] + [S.F(f'proton.periodic.{k}.d_N1_H1A') for k in prof]
    e = [0.0] + [S.F(f'periodic.profile.{k}.PBE0_kJ') for k in prof]
    assert all(b > a for a, b in zip(d, d[1:]))
    assert all(b > a for a, b in zip(e, e[1:])), 'the crystal profile rises at every step'
    assert S.F('periodic.profile.no_neutral_min') is True
    assert S.num(e[-1], 1) == '50.0'
    for k in prof:
        S.USED[f'periodic.profile.{k}.PBE0_kJ']['printed'].append('plotted (TOC graphic, no number printed)')
    return np.array(d), np.array(e)


def water_profile():
    keys = ['A00', 'A01', 'A02', 'B03', 'B02', 'B01', 'B00']      # lower envelope, as Fig. 2d
    assert S.F('scan.hysteresis_1.35') > 0 and S.F('scan.A03.Erel') > S.F('scan.B03.Erel')
    d = np.array([S.F(f'scan.{k}.d_NH') for k in keys])
    e = np.array([S.F(f'scan.{k}.Erel') for k in keys]) * KCAL_KJ
    assert abs(e[-1] - S.F('scan.Omin_kJ')) < 1e-3 and e[-1] < 0, 'the neutral end lies below the salt'
    assert abs(S.F('scan.barrier') * KCAL_KJ - S.F('scan.barrier_kJ')) < 1e-3
    for k in keys:
        S.USED[f'scan.{k}.Erel']['printed'].append('plotted (TOC graphic, no number printed)')
    return d, e


# ------------------------------------------------------------------------------------------------ docking
def dock_pose():
    contacts = json.loads((B.DOCK / 'results' / 'production_contacts.json').read_text())
    rep = contacts['T3_6RL9/sfa_neutral__ad4zn__prod']['1']
    lig = D.read_sdf(B.DOCK / 'results' / 'poses' / 'T3_6RL9' / 'sfa_neutral__ad4zn__top_pose_per_seed.sdf',
                     prop=('seed', rep['seed']))
    rec = D.read_receptor('T3_6RL9')
    zn = [a['xyz'] for a in rec if a['resname'] == 'ZN']
    assert len(zn) == 1
    zn = zn[0]
    iN = [i for i, s in enumerate(lig.symbols) if s == 'N']
    nS = min(iN, key=lambda i: np.linalg.norm(lig.coords[i] - zn))      # sulfonamide N (on the zinc)
    nA = max(iN, key=lambda i: np.linalg.norm(lig.coords[i] - zn))      # aniline N
    dzn = float(np.linalg.norm(lig.coords[nS] - zn))
    assert abs(dzn - S.F('bulk.dockcont.T3_6RL9.sfa_neutral__ad4zn__prod.zn')) < 0.006, dzn
    assert rep['zinc']['nearest_ligand_atom'].startswith('N(sulfonamide)')
    return lig, rec, zn, nS, nA


def dock_view(lig, zn, nA):
    """Zinc at the bottom, the ligand's long axis (Zn -> aniline N) up; turned DOCK_TURN about that axis from the
    view with the benzene ring face-on, so that no residue atom covers the ligand or the zinc and the three His
    ligands of the zinc stand apart (checked in dock_scene)."""
    L = lig.coords
    up = r3._unit(L[nA] - zn)
    ring = [i for i, s in enumerate(lig.symbols) if s == 'C']
    nrm = np.linalg.svd(L[ring] - L[ring].mean(0))[2][2]
    d = r3._unit(nrm - (nrm @ up) * up)
    return r3.View(L.mean(0), r3._rot(up, DOCK_TURN) @ d, up)


def dock_scene(lig, rec, zn, nS, view):
    """Sulfanilamide (ball and stick, orange carbons), the zinc, its three His ligands and Thr199 (grey sticks),
    zinc-N contacts (dashed), and the other residues within 4 A (pale thin sticks; those drawn in front of the
    ligand and over it on screen, and those reaching above the ligand where the labels go, are left out)."""
    keys, _ = D.pocket(rec, [lig], cutoff=4.0, extra=('ZN',))
    atoms = [a for a in D.pocket_atoms(rec, keys, [lig]) if a['resname'] in D.RES3]
    # zinc ligands: His N within 2.3 A of the zinc (hCA II: His94, His96, His119)
    znN = [a for a in rec if a['resname'] == 'HIS' and a['el'] == 'N' and np.linalg.norm(a['xyz'] - zn) < 2.3]
    his = sorted({a['resnum'] for a in znN})
    assert his == [94, 96, 119], his
    main = {('HIS', k) for k in his} | {('THR', 199)}
    assert ('A', 199, 'THR') in keys
    hi = [a for a in rec if (a['resname'], a['resnum']) in main and a['name'] not in ('N', 'C', 'O')]
    sL = view.screen(lig.coords)
    lo_res = {}
    for a in atoms:
        if (a['resname'], a['resnum']) not in main:
            lo_res.setdefault((a['resname'], a['resnum']), []).append(a)
    lo = []
    for k, grp in lo_res.items():
        sG = view.screen(np.array([a['xyz'] for a in grp]))
        over = np.hypot(sG[:, None, 0] - sL[None, :, 0], sG[:, None, 1] - sL[None, :, 1]).min() < 1.2
        front = sG[:, 2].max() > sL[:, 2].min()
        above = sG[:, 1].max() > sL[:, 1].max() - 0.5          # would sit among the labels above the ligand
        if not ((over and front) or above):
            lo += grp
    sH, sZ = view.screen(np.array([a['xyz'] for a in hi])), view.screen(zn[None])[0]
    for p in np.vstack([sL, sZ[None]]):
        assert not ((np.hypot(sH[:, 0] - p[0], sH[:, 1] - p[1]) < 0.7) & (sH[:, 2] > p[2])).any(), 'occluded'
    assert abs(np.linalg.svd(lig.coords[[i for i, s in enumerate(lig.symbols) if s == 'C']]
                             - lig.coords[[i for i, s in enumerate(lig.symbols) if s == 'C']].mean(0))[2][2]
               @ view.direction) > 0.75, 'benzene ring not face-on enough'
    sc = r3.Scene(None)
    for grp, col, rad, fade in ((lo, '#9A9A9A', 0.06, 0.5), (hi, C_RES, 0.10, 0.0)):
        sym, X, bonds = D.residue_bonds(grp)
        cols = [pale(col if s == 'C' else style.ELEMENT_COLOURS[s], fade) for s in sym]
        sc.add(D._stick_mesh(X, bonds, cols, rad), scalars='rgb', rgb=True, smooth_shading=True, specular=0.2,
               specular_power=20, diffuse=0.8, ambient=0.25)
    rzn = 0.55
    sc.add(r3._sphere_mesh([zn], [rzn], [D.METAL_COLOURS['ZN']], 40), scalars='rgb', rgb=True, smooth_shading=True,
           specular=0.4, specular_power=30)
    rad = ball_stick(sc, list(lig.symbols), lig.coords, carbon=C_SFA, atom_scale=0.95, bond_radius=0.1)
    for a in znN:
        dashed(sc, zn, a['xyz'], rzn + 0.05, 0.12, colour='#5A5C8C', radius=0.045)
    dashed(sc, zn, lig.coords[nS], rzn + 0.05, rad[nS] + 0.05, colour='#5A5C8C', radius=0.045)
    return sc, np.array([a['xyz'] for a in hi])


# ------------------------------------------------------------------------------------------------ figure
def tile(fig, x0, y0, x1, y1, img):
    """Picture tile with rounded corners (mm); img is an (h, w, 3) uint8 array of exactly the tile size."""
    ax = S.axes_mm(fig, x0, y0, x1 - x0, y1 - y0)
    h, w = img.shape[:2]
    im = ax.imshow(img, extent=(0, w, h, 0), interpolation='nearest')
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.set_axis_off()
    r = TILE_R * PX_MM
    clip = FancyBboxPatch((0, 0), w, h, boxstyle=f'round,pad=0,rounding_size={r}', transform=ax.transData,
                          facecolor='none', edgecolor='none')
    ax.add_patch(clip)
    im.set_clip_path(clip)
    return ax


def multiply(img, bg):
    """Render (white background) laid over a background image: multiply blend."""
    return np.clip(np.round(img.astype(float) * bg.astype(float) / 255.0), 0, 255).astype(np.uint8)


def cavity_background(r, X, sym):
    """Pale blue continuum with a white cavity: union of CAVITY_SCALE x Bondi discs around the projected atoms,
    white inside, blending into the continuum tint over CAVITY_EDGE outside (smoothstep)."""
    h, w = r.image.shape[:2]
    xy, _ = r.project(X)
    yy, xx = np.mgrid[0:h, 0:w]
    xx = xx + 0.5
    yy = yy + 0.5
    dmin = np.full((h, w), np.inf)
    for (px, py), s in zip(xy, sym):
        dmin = np.minimum(dmin, np.hypot(xx - px, yy - py) - CAVITY_SCALE * BONDI[s] * r.px_per_A)
    t = np.clip(dmin / (CAVITY_EDGE * r.px_per_A), 0.0, 1.0)
    t = t * t * (3 - 2 * t)
    bg = np.array(mcolors.to_rgb(TINT_WATER)) * 255
    out = 255.0 * (1 - t[..., None]) + bg[None, None] * t[..., None]
    return np.clip(np.round(out), 0, 255).astype(np.uint8)


def flat(c, shape):
    return np.tile((np.array(mcolors.to_rgb(c)) * 255).round().astype(np.uint8), (shape[0], shape[1], 1))


def glow(bg, r, X, colour, radius_A=GLOW_R, strength=0.75):
    """Soft spot of `colour` behind H1A, painted into the background before the multiply blend."""
    h, w = bg.shape[:2]
    (px, py), = r.project(X[[S.IDX['H1A']]])[0]
    yy, xx = np.mgrid[0:h, 0:w]
    q = np.hypot(xx + 0.5 - px, yy + 0.5 - py) / (radius_A * r.px_per_A)
    a = strength * np.clip(1.0 - q, 0.0, 1.0) ** 0.8
    c = np.array(mcolors.to_rgb(colour)) * 255
    out = bg.astype(float) * (1 - a[..., None]) + c[None, None] * a[..., None]
    return np.clip(np.round(out), 0, 255).astype(np.uint8)


def tag(fig, x_mm, y_mm, text, colour, bg=None):
    """Tile name (bold, FS_TAG) with its top-left corner at (x, y) mm; optional background box."""
    return fig.text(x_mm / W_MM, y_mm / H_MM, text, fontsize=FS_TAG, fontweight='bold', color=colour, ha='left',
                    va='top',
                    bbox=None if bg is None else dict(boxstyle='round,pad=0.22', facecolor=bg, edgecolor='none'))


def draw_plot(fig):
    from scipy.interpolate import CubicSpline, PchipInterpolator
    x0, x1, y0, y1 = PLOT
    ax = S.axes_mm(fig, x0, y0, x1 - x0, y1 - y0)
    dc, ec = crystal_profile()
    dw, ew = water_profile()
    xl = (1.02, 1.76)
    yl = (-7.5, 55.0)
    ax.set_xlim(*xl)
    ax.set_ylim(*yl)
    ax.axhline(0.0, color='0.55', lw=0.6, ls=(0, (1.0, 1.6)), zorder=1)
    xc = np.linspace(dc[0], dc[-1], 400)
    ax.plot(xc, PchipInterpolator(dc, ec)(xc), color=C_XTAL, lw=1.5, zorder=3, solid_capstyle='round')
    xw = np.linspace(dw[0], dw[-1], 400)
    ax.plot(xw, CubicSpline(dw, ew)(xw), color=C_WATER, lw=1.5, zorder=3, solid_capstyle='round')
    # axes as two arrows, no ticks or numbers
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.patch.set_alpha(0.0)
    kw = dict(xycoords='axes fraction', textcoords='axes fraction',
              arrowprops=dict(arrowstyle='-|>,head_length=0.35,head_width=0.18', lw=0.8, color=INK,
                              shrinkA=0, shrinkB=0))
    ax.annotate('', xy=(1.0, 0.0), xytext=(0.0, 0.0), **kw)
    ax.annotate('', xy=(0.0, 1.0), xytext=(0.0, 0.0), **kw)
    ax.text(-0.05, 1.0, '$E$', transform=ax.transAxes, ha='right', va='top', fontsize=FS_TAG, color=INK)
    ax.text(0.5, -0.045, 'H$^{+}$: N → O', transform=ax.transAxes, ha='center', va='top', fontsize=FS_LAB, color=INK)
    ax.text(1.06, 47.0, 'crystal', ha='left', va='center', fontsize=FS_LAB, color=C_XTAL, fontweight='bold')
    ax.text(1.74, 5.0, 'water', ha='right', va='bottom', fontsize=FS_LAB, color=C_WATER, fontweight='bold')
    return ax


def place_label(fig, ax, img, bg_rgb, anchor, text, offsets, fontsize=FS_LAB, colour=INK):
    """Put `text` next to `anchor` (tile pixels) at the candidate offset (dx, dy px, ha, va) whose box covers the
    least ink of the tile image (pixels differing from the tile background); first candidate wins ties."""
    rend = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    bg = np.array(bg_rgb, float)
    ink = np.abs(img.astype(float) - bg[None, None]).max(axis=2) > 24
    best = None
    for dx, dy, ha, va in offsets:
        t = ax.text(anchor[0] + dx, anchor[1] + dy, text, ha=ha, va=va, fontsize=fontsize, color=colour)
        bb = t.get_window_extent(rend)
        (x0, y1), (x1, y0) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        h, w = ink.shape
        xa, xb = int(max(0, np.floor(x0))), int(min(w, np.ceil(x1)))
        ya, yb = int(max(0, np.floor(y0))), int(min(h, np.ceil(y1)))
        inside = xa >= 2 and ya >= 2 and xb <= w - 2 and yb <= h - 2
        cost = ink[ya:yb, xa:xb].mean() if (xb > xa and yb > ya) else 1.0
        cost += 0.0 if inside else 10.0
        t.remove()
        if best is None or cost < best[0] - 1e-9:
            best = (cost, dx, dy, ha, va)
    _, dx, dy, ha, va = best
    return ax.text(anchor[0] + dx, anchor[1] + dy, text, ha=ha, va=va, fontsize=fontsize, color=colour), best[0]


def main():
    from dataclasses import replace
    plt = S.setup()
    plt.rcParams['savefig.bbox'] = None             # the canvas is the ACS box: never crop
    plt.rcParams['savefig.pad_inches'] = 0.0
    fig = S.figure_mm(W_MM, H_MM)
    fig.patch.set_facecolor('white')

    # ---- crystal and water tiles: one scale for both pairs
    P, ions, layer, n = crystal_layer()
    Wp = water_pair()
    views = {'xtal': pair_view(P, n), 'water': pair_view(Wp, ring_normal(Wp))}
    geoms = {'xtal': P, 'water': Wp}
    tw = XT[1] - XT[0]
    th = (H_MM - YT_GAP) / 2
    rad = np.array([style.BALL_RADII[s] for s in S.SYMBOLS])
    ext = {}
    for k, v in views.items():
        s_ = v.screen(geoms[k])
        ext[k] = ((s_[:, 0] - rad).min(), (s_[:, 0] + rad).max(), (s_[:, 1] - rad).min(), (s_[:, 1] + rad).max())
    span_w = max(e[1] - e[0] for e in ext.values())
    span_h = max(e[3] - e[2] for e in ext.values())
    mm_per_A = min((tw - 2.0) / span_w, (th - TAG_RES - 0.6) / span_h)
    px_per_A = mm_per_A * PX_MM
    renders = {}
    for k, v in views.items():
        e = ext[k]
        # pair centred horizontally, a little below the tile centre (room for the tag)
        f = v.focal + 0.5 * (e[0] + e[1]) * v.right + (0.5 * (e[2] + e[3]) + 0.5 * TAG_RES / mm_per_A) * v.up
        v = replace(v, focal=f)
        ctx = []
        if k == 'xtal':                                    # lattice ions that show inside the tile
            hw, hh = tw / mm_per_A / 2 + 0.5, th / mm_per_A / 2 + 0.5
            for sym, Y in ions:
                sY = v.screen(Y)
                if ((np.abs(sY[:, 0]) < hw) & (np.abs(sY[:, 1]) < hh)).any():
                    ctx.append((sym, Y))
            B.log(f'crystal layer {layer}; {len(ctx)} lattice ions in the tile; pair scale {mm_per_A:.2f} mm/A')
        sc = pair_scene(geoms[k], 'O5', ctx) if k == 'xtal' else pair_scene(geoms[k], 'N1')
        renders[k] = render_tile(sc, v, tw, th, px_per_A)
    rx, rw = renders['xtal'], renders['water']
    img_x = multiply(rx.image, glow(flat(TINT_XTAL, rx.image.shape), rx, P, pale(C_XTAL, 0.25)))
    img_w = multiply(rw.image, glow(cavity_background(rw, Wp, S.SYMBOLS), rw, Wp, pale(C_WATER, 0.3)))
    tile(fig, XT[0], th + YT_GAP, XT[1], H_MM, img_x)
    tile(fig, XT[0], 0.0, XT[1], th, img_w)
    tag(fig, XT[0] + 1.5, H_MM - 1.2, 'crystal', C_XTAL, bg=TINT_XTAL)
    tag(fig, XT[0] + 1.5, th - 1.2, 'water continuum', C_WATER)

    # ---- energy plot
    draw_plot(fig)

    # ---- docking tile
    lig, rec, zn, nS, nA = dock_pose()
    vd = dock_view(lig, zn, nA)
    dsc, hi_xyz = dock_scene(lig, rec, zn, nS, vd)
    dw_mm = DOCK[1] - DOCK[0]
    s_ = vd.screen(np.vstack([lig.coords, [zn], hi_xyz]))
    # scale: ligand, zinc and zinc ligands fill the tile width (side-chain stems may meet the edge)
    mm_d = min((dw_mm + 1.0) / np.ptp(s_[:, 0]), (H_MM - 11.0) / np.ptp(s_[:, 1]))
    cu, cv = 0.5 * (s_[:, 0].min() + s_[:, 0].max()), 0.5 * (s_[:, 1].min() + s_[:, 1].max())
    vd = replace(vd, focal=vd.focal + cu * vd.right + (cv + 2.5 / mm_d) * vd.up)
    rd = render_tile(dsc, vd, dw_mm, H_MM, mm_d * PX_MM)
    img_d = multiply(rd.image, flat(TINT_DOCK, rd.image.shape))
    ax_d = tile(fig, DOCK[0], 0.0, DOCK[1], H_MM, img_d)
    B.log(f'docking scale {mm_d:.2f} mm/A')
    tag(fig, DOCK[0] + 1.5, H_MM - 1.2, 'hCA II', INK, bg=TINT_DOCK)
    fig.canvas.draw()
    bg_d = np.array(mcolors.to_rgb(TINT_DOCK)) * 255
    pz = rd.px_per_A
    (zx, zy), = rd.project(np.array([zn]))[0]
    g = 0.75 * pz
    place_label(fig, ax_d, img_d, bg_d, (zx, zy), 'Zn',
                [(g, 0, 'left', 'center'), (-g, 0, 'right', 'center'), (0, g, 'center', 'top'),
                 (g * 0.8, g * 0.8, 'left', 'top'), (-g * 0.8, g * 0.8, 'right', 'top')])
    (ax_, ay_), = rd.project(lig.coords[[nA]])[0]
    g = 0.55 * pz
    place_label(fig, ax_d, img_d, bg_d, (ax_, ay_), 'sulfanilamide',
                [(0, -g, 'center', 'bottom'), (g, -g * 0.3, 'left', 'bottom'), (-g, -g * 0.3, 'right', 'bottom')])

    # ---- checks and output
    lax = bool(os.environ.get('TOC_LAX'))              # drafts only: report QA issues instead of stopping
    issues = qa.check_figure(fig, min_pt=6.0, strict=not lax, name=STEM, glyphs=True)
    assert lax or not issues
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=DPI, facecolor='white')
    from PIL import Image, PngImagePlugin
    buf.seek(0)
    im = Image.open(buf).convert('RGB')
    assert im.size == (round(W_IN * DPI), round(H_IN * DPI)), im.size
    used = {k: v['value'] for k, v in sorted(S.USED.items())}
    src = ('paper/main/figures/scripts/toc_graphic.py; facts (paper/main/data/facts.json) plotted or checked, no '
           'number printed: ' + ', '.join(used))
    OUT.mkdir(parents=True, exist_ok=True)
    info = PngImagePlugin.PngInfo()
    info.add_text('Title', 'Graphical abstract' if FORMAT == 'elsevier' else 'Table of contents graphic')
    info.add_text('Source', src)
    png, tif = OUT / f'{STEM}.png', OUT / f'{STEM}.tif'
    im.save(png, dpi=(DPI, DPI), pnginfo=info, optimize=True)
    im.save(tif, dpi=(DPI, DPI), compression='tiff_lzw', tiffinfo={270: src})
    qa.check_png(png, min_dpi=DPI, max_border_mm=2.5, max_width_mm=W_MM, strict=not lax)
    with Image.open(tif) as t:
        assert t.size == im.size and t.mode == 'RGB' and round(t.info['dpi'][0]) == DPI
    print(f'{STEM}: {im.size[0]} x {im.size[1]} px at {DPI} dpi = {W_MM:.2f} x {H_MM:.2f} mm; wrote {png} and {tif}')
    print('facts used:', json.dumps(used, default=str))


if __name__ == '__main__':
    main()
