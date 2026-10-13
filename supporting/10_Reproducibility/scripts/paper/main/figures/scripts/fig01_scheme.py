"""Fig. 1: the ion pair, its R2,2(8) ring and the nine hydrogen bonds of the crystal (DOUBLE width).

    cd /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate
    nice -n 10 env OMP_NUM_THREADS=2 /root/workspace/cqc-dft/tools/paper-env/bin/python \
        paper/main/figures/scripts/fig01_scheme.py

(a) Scheme of 4-sulfamoylanilinium dihydrogen phosphate as the ion pair used in the calculations: the cation and
    the anion at x, y, z - 1 joined by N1-H1A...O5 (bond 1) and C6-H6...O4 (bond 9), which close the R2,2(8) ring.
    Chemistry from RDKit (atom-mapped SMILES, sanitised, formal charges asserted: cation +1, anion -1); the 2D layout
    is RDKit CoordGen computed on a copy in which the two hydrogen bonds are temporary ring bonds, so the R2,2(8) ring
    is drawn as a regular eight-membered ring fused to the benzene ring (not to scale). Drawn by figkit.scheme.
(b) Key: the nine hydrogen bonds of the deposited table in CIF order (numbers used in c and d and in Table 4), with
    D...A and D-H...A from facts xray.hb.*; which panel shows each one.
(c) One hydrogen-bonded layer (all ions related by lattice translations only), viewed along b: six ion pairs,
    bonds 1 and 9 (R2,2(8)) along c, 5 and 7 along a, 6 and 8 diagonal.
(d) Packing viewed along a: the layers of (c) and their glide images (x, -y, z + 1/2), joined by bonds 2, 3 and 4;
    bonds 5-8 run along a and are seen end-on here.
Every ion comes from the CIF by the space-group operations of Pc (gemmi); every hydrogen bond drawn is identified
by distance (D...A within 0.002 A of the deposited value, H...A within 0.02 A, angle > 120 deg) and the symmetry
code of each bond is checked against the CIF table (_crystal.check_codes). The highlighted ion pair is the one
modelled by DFT (atoms identical to molecules/sfa/input/sfa_pair.xyz, asserted).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Polygon, FancyArrowPatch, Patch  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import rdDepictor  # noqa: E402

import _saph as S  # noqa: E402
import _crystal as C  # noqa: E402
from figkit import style, scheme, geometry as geo  # noqa: E402
from figkit.labels import LabelPlacer  # noqa: E402

STEM = 'fig01_scheme'
W_MM = 175.0
FS = style.SMALL_PT
GREY = '0.40'
TEAL, BROWN = S.HB_PAIR, S.HB_CRYSTAL
HALO = '#F6E08A'              # the modelled ion pair
RING_FILL = '#0E7C7B'

# ------------------------------------------------------------------------------------------------ (a) scheme
MAP = {1: 'N1', 11: 'C1', 12: 'C2', 13: 'C3', 14: 'C4', 15: 'C5', 16: 'C6', 21: 'S1', 22: 'O1', 23: 'O2',
       24: 'N2', 10: 'P1', 3: 'O3', 4: 'O4', 5: 'O5', 6: 'O6', 91: 'H1A', 92: 'H1B', 93: 'H1C', 94: 'H4A',
       95: 'H3A', 96: 'H6', 97: 'H2A', 98: 'H2B'}
# layout copy: H1A and H6 are dummy atoms (*) with a temporary bond to their acceptors O5 and O4
SMI_LAYOUT = ('[N+:1]1([H:92])([H:93])[*:91][O:5][P:10](=[O:6])([O:3][H:95])[O:4]([H:94])[*:96]'
              '[c:16]2[cH:15][c:14]([S:21](=[O:22])(=[O:23])[N:24]([H:97])[H:98])[cH:13][cH:12][c:11]12')
SMI_CATION = ('[N+:1]([H:91])([H:92])([H:93])[c:11]1[cH:12][cH:13][c:14]([S:21](=[O:22])(=[O:23])[N:24]([H:97])'
              '[H:98])[cH:15][c:16]1[H:96]')
SMI_ANION = '[O-:5][P:10](=[O:6])([O:3][H:95])[O:4][H:94]'


def _mol(smi, sanitize=True):
    p = Chem.SmilesParserParams()
    p.removeHs = False
    p.sanitize = sanitize
    m = Chem.MolFromSmiles(smi, p)
    assert m is not None, smi
    return m


def scheme_structures():
    lay = _mol(SMI_LAYOUT, sanitize=False)
    lay.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(lay)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(lay)
    conf = lay.GetConformer()
    pos = {a.GetAtomMapNum(): np.array([conf.GetAtomPosition(a.GetIdx()).x, conf.GetAtomPosition(a.GetIdx()).y])
           for a in lay.GetAtoms()}
    L = np.median([np.linalg.norm(pos[MAP_INV[a]] - pos[MAP_INV[b]]) for a, b in
                   (('C1', 'C2'), ('C2', 'C3'), ('C3', 'C4'), ('C4', 'C5'), ('C5', 'C6'), ('C1', 'C6'))])
    out = {}
    for name, smi, q in (('cation', SMI_CATION, +1), ('anion', SMI_ANION, -1)):
        m = _mol(smi)
        assert Chem.GetFormalCharge(m) == q, name
        Chem.Kekulize(m, clearAromaticFlags=True)
        maps = {}
        xy = np.zeros((m.GetNumAtoms(), 2))
        for a in m.GetAtoms():
            maps[MAP[a.GetAtomMapNum()]] = a.GetIdx()
            xy[a.GetIdx()] = pos[a.GetAtomMapNum()] / L
            a.SetAtomMapNum(0)
        out[name] = scheme.Structure(m, xy, maps=maps)
    # chemistry cross-check against the CIF: which atoms carry hydrogens
    cat = out['cation'].mol
    nH = {lab: cat.GetAtomWithIdx(i).GetTotalNumHs(includeNeighbors=True) for lab, i in out['cation'].maps.items()}
    assert nH['N1'] == 3 and nH['N2'] == 2 and nH['C6'] == 1
    an = out['anion'].mol
    nHa = {lab: an.GetAtomWithIdx(i).GetTotalNumHs(includeNeighbors=True) for lab, i in out['anion'].maps.items()}
    assert nHa['O3'] == 1 and nHa['O4'] == 1 and nHa['O5'] == 0 and nHa['O6'] == 0
    return out


MAP_INV = {v: k for k, v in MAP.items()}


def draw_scheme(ax, origin, bond):
    """Draw (a) on the mm axes; returns dict label -> position (mm) and the Drawn objects."""
    st = scheme_structures()
    drawn = {}
    P = {}
    for name in ('cation', 'anion'):
        s = st[name]
        hcol = {s.maps['H1A']: TEAL} if name == 'cation' else {}
        d = scheme.draw_structure(ax, s, origin=origin, bond=bond, fontsize=8.0, lw=0.8, atom_colours=hcol)
        drawn[name] = (s, d)
        for lab, i in s.maps.items():
            P[lab] = d.pos[i]
    # R2,2(8) ring shading and the two hydrogen bonds (bond numbers as in the key)
    ring = ['H1A', 'N1', 'C1', 'C6', 'H6', 'O4', 'P1', 'O5']
    ax.add_patch(Polygon([P[l] for l in ring], closed=True, facecolor=RING_FILL, alpha=0.10, edgecolor='none',
                         zorder=0))
    cen = np.mean([P[l] for l in ring], 0)
    ax.text(cen[0], cen[1], 'R$_2^2$(8)', fontsize=7.0, ha='center', va='center', color=TEAL, zorder=5)
    segs = []
    for (dn, hn, an), n in ((('N1', 'H1A', 'O5'), 1), (('C6', 'H6', 'O4'), 9)):
        sh, dh = drawn['cation']
        sa, da = drawn['anion']
        p, q = P[hn], P[an]
        bh = dh.boxes.get(sh.maps[hn])
        ba = da.boxes.get(sa.maps[an])
        t0 = scheme._clip_t(p, q, bh) if bh else 0.0
        t1 = scheme._clip_t(q, p, ba) if ba else 0.0
        p0, p1 = p + t0 * (q - p), q - t1 * (q - p)
        p0 = p0 + 0.5 * (q - p) / np.linalg.norm(q - p)
        p1 = p1 - 0.5 * (q - p) / np.linalg.norm(q - p)
        ax.add_line(Line2D([p0[0], p1[0]], [p0[1], p1[1]], color=TEAL, lw=1.1, ls=(0, (2.0, 1.4)), zorder=4))
        segs.append((p0, p1, n))
    return P, drawn, segs


def label_scheme(ax, P, drawn, segs):
    pl = LabelPlacer(ax, name='scheme labels', n_rings=4, clearance_pt=1.0, bounds=None)
    for name in ('cation', 'anion'):
        scheme.add_obstacles(pl, drawn[name][1])
    for p0, p1, n in segs:
        pl.add_segment(p0, p1, hw_pt=1.0)
    for t in ax.texts:
        pl.add_text_obstacle(t)
    # bond numbers: outside the ring, next to the middle of each dashed line
    ring_c = np.mean([P[l] for l in ['H1A', 'N1', 'C1', 'C6', 'H6', 'O4', 'P1', 'O5']], 0)
    for p0, p1, n in segs:
        mid = (p0 + p1) / 2
        pl.add_item(str(n), mid, radius_pt=1.5, prefer=mid - ring_c, priority=0, leader=False, max_ring=2,
                    kw=dict(fontsize=7.0, color=TEAL, fontweight='bold'), key=('bond', n))
    labs = ['N1', 'H1A', 'O5', 'P1', 'O4', 'H6', 'S1', 'O1', 'O2', 'N2', 'O3', 'O6', 'H3A', 'H4A']
    for lab in labs:
        pl.add_item(lab, P[lab], radius_pt=4.4 if not lab.startswith(('C',)) else 2.2, priority=1, max_ring=3,
                    leader=False, kw=dict(fontsize=FS, color=GREY), key=lab)
    placed, dropped = pl.place()
    assert not dropped, [d.text for d in dropped]
    return placed


# ------------------------------------------------------------------------------------------------ crystal panels
def layer_ions():
    """(c): six ion pairs of one layer (translations only): cation (n, 0, m) + its partner anion (n, 0, m - 1)."""
    out = []
    for n in (-1, 0, 1):
        for m in (0, 1):
            out += [C.ion('cation', 1, (n, 0, m)), C.ion('anion', 1, (n, 0, m - 1))]
    return out


def packing_ions():
    """(d): every ion (both operations) with centroid in 0 <= x < 1 (one copy per stack along a) in the window."""
    return C.ions(((0.0, 1.0), (-0.30, 1.02), (-0.55, 1.00)))


REF = [('cation', 1, (0, 0, 0)), ('anion', 1, (0, 0, -1))]      # the modelled ion pair


def check_ref_pair():
    X = S.xray_pair()
    cat, an = C.ion(*REF[0]), C.ion(*REF[1])
    Y = np.vstack([cat.X, an.X])
    assert cat.labels + an.labels == S.PAIR_LABELS
    assert np.abs(Y - X).max() < 1e-3, np.abs(Y - X).max()


def ion_bonds(io):
    return geo.detect_bonds(io.symbols, io.X)


def draw_ions(ax, ions, R, lw=0.75, ref_keys=()):
    """Depth-sorted element-coloured wireframe (halves of each bond in the colour of its atom); the modelled pair
    gets a pale halo underneath."""
    prims = []
    for io in ions:
        X = io.X @ R.T
        for i, j in ion_bonds(io):
            m = (X[i] + X[j]) / 2
            for a, b, s in ((X[i], m, io.symbols[i]), (m, X[j], io.symbols[j])):
                col = '#A0A0A0' if s == 'H' else style.ELEMENT_COLOURS[s]
                prims.append(((a[2] + b[2]) / 2, a, b, col, (lw * 0.7 if s == 'H' else lw)))
            if io.key in ref_keys:
                ax.add_line(Line2D([X[i, 0], X[j, 0]], [X[i, 1], X[j, 1]], color=HALO, lw=4.6,
                                   solid_capstyle='round', zorder=0.5))
    prims.sort(key=lambda t: t[0])
    zs = np.linspace(2.0, 3.0, len(prims))
    for z, (_, a, b, c, w) in zip(zs, prims):
        ax.add_line(Line2D([a[0], b[0]], [a[1], b[1]], color=c, lw=w, solid_capstyle='round', zorder=z))


def ring_polygons(ax, ions, R):
    """Shade every complete R2,2(8) ring (cation and its N1-H1A...O5 partner anion) among the ions."""
    keys = {io.key: io for io in ions}
    n = 0
    for io in ions:
        if io.kind != 'cation':
            continue
        k = ('anion', io.op, (io.tr[0], io.tr[1], io.tr[2] - 1))
        if k not in keys:
            continue
        an = keys[k]
        pts = [io.at('H1A'), io.at('N1'), io.at('C1'), io.at('C6'), io.at('H6'), an.at('O4'), an.at('P1'),
               an.at('O5')]
        assert abs(np.linalg.norm(io.at('N1') - an.at('O5')) - S.F('xray.hb.N1_H1A_O5.DA')) < 0.002
        V = np.array(pts) @ R.T
        ax.add_patch(Polygon(V[:, :2], closed=True, facecolor=RING_FILL, alpha=0.13, edgecolor='none', zorder=1))
        n += 1
    return n


def draw_hbonds(ax, ions, hbs, R):
    segs = []
    for h in hbs:
        H = ions[h.h[0]].at(h.h[1]) @ R.T
        A = ions[h.a[0]].at(h.a[1]) @ R.T
        col = TEAL if C.HB[h.n - 1] in C.IN_PAIR else BROWN
        ax.add_line(Line2D([H[0], A[0]], [H[1], A[1]], color=col, lw=0.95, ls=(0, (2.0, 1.2)), zorder=3.5,
                           dash_capstyle='butt'))
        segs.append((H[:2], A[:2], h))
    return segs


def number_labels(ax, ions, R, segs, ref_keys, only, merge=None, extra=()):
    """One number per bond type, on the instance nearest the modelled pair; placed without overlaps."""
    pl = LabelPlacer(ax, name='bond numbers', n_rings=5, clearance_pt=1.2, bounds='axes')
    for io in ions:
        X = io.X @ R.T
        for i, j in ion_bonds(io):
            pl.add_segment(X[i, :2], X[j, :2], hw_pt=0.8)
        for k in range(len(X)):
            pl.add_disc(X[k, :2], r_pt=1.2)
    for p0, p1, h in segs:
        pl.add_segment(p0, p1, hw_pt=0.7)
    for p0, p1 in extra:
        pl.add_segment(p0, p1, hw_pt=1.2)
    for t in ax.texts:
        pl.add_text_obstacle(t)
    refc = np.mean([io.X.mean(0) for io in ions if io.key in ref_keys], 0) @ R.T
    chosen = {}
    for p0, p1, h in segs:
        if h.n not in only:
            continue
        mid = (p0 + p1) / 2
        onref = ions[h.d[0]].key in ref_keys or ions[h.a[0]].key in ref_keys
        score = (0 if onref else 1, float(np.linalg.norm(mid - refc[:2])))
        if h.n not in chosen or score < chosen[h.n][0]:
            chosen[h.n] = (score, mid, p1 - p0)
    merged = {}
    for m, into in (merge or {}).items():
        # bonds whose drawn segments share the acceptor position in this projection get one label ("2, 3")
        if m in chosen and into in chosen:
            ends = {h.n: p1 for p0, p1, h in segs if h.n in (m, into)}
            assert np.linalg.norm(ends[m] - ends[into]) < 0.05, (m, into)
            merged[into] = f'{into}, {m}'
            del chosen[m]
    for n, (_, mid, dv) in sorted(chosen.items()):
        perp = np.array([-dv[1], dv[0]])
        col = TEAL if C.HB[n - 1] in C.IN_PAIR else BROWN
        pl.add_item(merged.get(n, str(n)), mid, radius_pt=0.8, prefer=perp, priority=0, max_ring=4,
                    kw=dict(fontsize=7.0, color=col, fontweight='bold'), key=('bond', n))
    for into, txt in merged.items():
        chosen[int(txt.split(', ')[1])] = None
    placed, dropped = pl.place()
    assert not dropped, [d.text for d in dropped]
    for it in placed:
        it.artist.set_path_effects([pe.withStroke(linewidth=1.6, foreground='white')])
    return sorted(chosen)


def cell_outline(ax, R, corners):
    M = C.structure()['M']
    P = np.array(corners, float) @ M.T @ R.T
    ax.add_line(Line2D(P[:, 0], P[:, 1], color='0.45', lw=0.5, zorder=0.8))
    return [(P[k, :2], P[k + 1, :2]) for k in range(len(P) - 1)]


def axis_arrows(ax, R, origin_xy, names, length=2.2, fs=7.0):
    """Short arrows along the projected cell axes named in `names` (subset of 'abc'), from origin_xy (data)."""
    M = C.structure()['M']
    segs = []
    for k, nm in enumerate('abc'):
        if nm not in names:
            continue
        v = (M[:, k] / np.linalg.norm(M[:, k])) @ R.T
        v2 = v[:2] / np.linalg.norm(v[:2])
        p1 = np.array(origin_xy) + length * v2
        ax.add_patch(FancyArrowPatch(tuple(origin_xy), tuple(p1), arrowstyle='-|>,head_length=2.2,head_width=1.1',
                                     mutation_scale=1.0, lw=0.7, color='black', shrinkA=0, shrinkB=0, zorder=6))
        ax.text(p1[0] + 0.45 * v2[0], p1[1] + 0.45 * v2[1], nm, fontsize=fs, style='italic', ha='center',
                va='center', zorder=6)
        segs.append((np.array(origin_xy), p1 + 0.9 * v2))
    return segs


def extent(ions, R, pad=0.6):
    X = np.vstack([io.X for io in ions]) @ R.T
    return X[:, 0].min() - pad, X[:, 0].max() + pad, X[:, 1].min() - pad, X[:, 1].max() + pad


# ------------------------------------------------------------------------------------------------ figure
def main():
    plt = S.setup()
    C.check_codes()
    check_ref_pair()
    # crystal panels first: their common scale fixes the height of the figure
    Rc = np.array([[0, 0, 1.0], [1.0, 0, 0], [0, 1.0, 0]])       # right = c*, up = a, viewer on +b
    Rd = np.array([[0, 0, 1.0], [0, 1.0, 0], [-1.0, 0, 0]])      # right = c*, up = b, viewer on -a
    for Rm in (Rc, Rd):
        assert abs(np.linalg.det(Rm) - 1) < 1e-12            # proper rotations: no mirror images
    ic, idd = layer_ions(), packing_ions()
    ec, ed = extent(ic, Rc, pad=0.35), extent(idd, Rd, pad=0.35)
    wc_av = 84.0
    scale = min(wc_av / (ec[1] - ec[0]), wc_av / (ed[1] - ed[0]))      # mm per Angstrom
    hmax = max(ec[3] - ec[2], ed[3] - ed[2]) * scale
    TOP_ROW = 56.0
    H = TOP_ROW + 5.5 + hmax + 0.3
    assert H <= style.MAX_HEIGHT / S.MM + 1e-6, H
    W = W_MM
    fig = S.figure_mm(W, H)
    axm = scheme.mm_axes(fig)
    # ---------------- (a) scheme
    P, drawn, segs = draw_scheme(axm, origin=(48.0, H - 29.5), bond=6.3)
    y_names = H - 46.5
    label_scheme(axm, P, drawn, segs)
    axm.text(0.4, H - 3.4, '(a)', fontsize=style.PANEL_PT, fontweight='bold', va='baseline')
    axm.text(7.0, H - 3.4, 'The ion pair and its R$_2^2$(8) ring', fontsize=style.TICK_PT + 0.5, va='baseline')
    fig.canvas.draw()
    rr = fig.canvas.get_renderer()
    inv = axm.transData.inverted()
    for name, text in (('cation', '4-sulfamoylanilinium'), ('anion', 'dihydrogen phosphate')):
        s_, d_ = drawn[name]
        bbs = [a.get_window_extent(rr) for a in d_.texts + d_.lines]
        (x0_, y0_), (x1_, _) = inv.transform([[min(b.x0 for b in bbs), min(b.y0 for b in bbs)],
                                               [max(b.x1 for b in bbs), 0]])
        axm.text((x0_ + x1_) / 2, y_names, text, fontsize=FS + 0.5, ha='center', va='top', color='0.15')
    # ---------------- (b) key
    xk = 104.0
    axm.text(xk - 6.6, H - 3.4, '(b)', fontsize=style.PANEL_PT, fontweight='bold', va='baseline')
    axm.text(xk, H - 3.4, 'Hydrogen bonds of the crystal (X-ray)', fontsize=style.TICK_PT + 0.5, va='baseline')
    cols = [xk, xk + 5.0, xk + 44.0, xk + 60.0, xk + 63.5]      # no., bond, D...A (right), angle (right), panel
    y = H - 9.6
    heads = [('', 'left'), ('Bond', 'left'), ('D···A (Å)', 'right'), ('D–H···A (°)', 'right'), ('Panel', 'left')]
    for x, (t, ha) in zip(cols, heads):
        axm.text(x, y, t, fontsize=FS, ha=ha, va='baseline', color='0.25')
    axm.add_line(Line2D([xk, xk + 69.0], [y - 1.3, y - 1.3], color='0.6', lw=0.4))
    where = {1: 'c, d', 2: 'd', 3: 'd', 4: 'd', 5: 'c', 6: 'c', 7: 'c', 8: 'c', 9: 'c, d'}
    for n, (d_, h_, a_) in enumerate(C.HB, start=1):
        y -= 3.55
        key = f'xray.hb.{d_}_{h_}_{a_}'
        col = TEAL if (d_, h_, a_) in C.IN_PAIR else BROWN
        axm.text(cols[0] + 1.2, y, str(n), fontsize=FS + 0.5, ha='center', va='baseline', color=col,
                 fontweight='bold')
        axm.text(cols[1], y, f'{d_}–{h_}···{a_}', fontsize=FS, ha='left', va='baseline')
        axm.text(cols[2], y, S.txt(key + '.DA'), fontsize=FS, ha='right', va='baseline')
        axm.text(cols[3], y, S.txt(key + '.DHA'), fontsize=FS, ha='right', va='baseline')
        axm.text(cols[4], y, where[n], fontsize=FS, ha='left', va='baseline', color='0.25')
    assert S.F('xray.n_hbonds') == 9
    # legend of the crystal panels, under the key
    y -= 6.0
    items = [(Line2D([], [], color=TEAL, lw=1.1, ls=(0, (2.0, 1.4))), 'kept in the ion-pair model'),
             (Line2D([], [], color=BROWN, lw=1.1, ls=(0, (2.0, 1.4))), 'crystal only'),
             (Line2D([], [], color=HALO, lw=4.6, solid_capstyle='round'), 'ion pair modelled by DFT'),
             (Patch(facecolor=RING_FILL, alpha=0.13, edgecolor='none'), 'R$_2^2$(8) ring')]
    leg = axm.legend([h for h, _ in items], [t for _, t in items], loc='upper left',
                     bbox_to_anchor=(xk / W, (y + 2.5) / H), bbox_transform=fig.transFigure, ncol=2,
                     fontsize=FS, handlelength=2.0, columnspacing=1.2, handletextpad=0.5, borderaxespad=0.0,
                     frameon=False)
    # ---------------- (c) layer along b and (d) packing along a, one common scale
    top2 = H - TOP_ROW
    refk = {(k, o, t) for k, o, t in REF}
    panels = []
    for (ions, Rm, ex, x0, letter, title, only) in (
            (ic, Rc, ec, 1.0, 'c', 'One layer, viewed along $b$', {1, 5, 6, 7, 8, 9}),
            (idd, Rd, ed, 89.0, 'd', 'Packing viewed along $a$', {1, 2, 3, 4, 9})):
        w = (ex[1] - ex[0]) * scale
        h = (ex[3] - ex[2]) * scale
        ax = S.axes_mm(fig, x0 + (wc_av - w) / 2, top2 - 5.5 - h, w, h)
        ax.set_xlim(ex[0], ex[1])
        ax.set_ylim(ex[2], ex[3])
        ax.set_aspect('equal', adjustable='box')
        ax.axis('off')
        nring = ring_polygons(ax, ions, Rm)
        draw_ions(ax, ions, Rm, ref_keys=refk)
        if letter == 'c':
            hbs = C.hbonds(ions)
        else:
            # projection along a: the acceptor of bond 3 lies one translation along a from the drawn stack; take
            # acceptors from the neighbouring stacks too and draw them at their (identical) projected position
            pool = ions + [C.ion(io.kind, io.op, (io.tr[0] + s, io.tr[1], io.tr[2])) for io in ions for s in (-1, 1)]
            hbs = [h for h in C.hbonds(pool, donors=range(len(ions))) if h.n in only]
            ions = pool
        found = sorted({h.n for h in hbs})
        assert set(found) >= only, (letter, found)
        segs = draw_hbonds(ax, ions, [h for h in hbs if h.n in only], Rm)
        if letter == 'c':
            extra = axis_arrows(ax, Rm, (ex[0] + 0.6, ex[2] + 0.6), 'ac')
        else:
            extra = cell_outline(ax, Rm, [[0, 0, 0], [0, 1, 0], [0, 1, 1], [0, 0, 1], [0, 0, 0]])
            extra += axis_arrows(ax, Rm, (ex[0] + 0.6, ex[2] + 0.6), 'bc')
        numbered = number_labels(ax, ions, Rm, segs, refk, only, merge={3: 2} if letter == 'd' else None,
                                 extra=extra)
        assert set(numbered) == only, (letter, numbered)
        panels.append((ax, letter, title, x0))
        print(f'({letter}) {len(ions)} ions, bonds {found}, {nring} R22(8) rings, numbered {numbered}')
    for ax, letter, title, x0 in panels:
        axm.text(x0 - 0.6, top2, f'({letter})', fontsize=style.PANEL_PT, fontweight='bold', va='baseline')
        axm.text(x0 + 6.0, top2, title, fontsize=style.TICK_PT + 0.5, va='baseline')
    S.save(fig, STEM, W_MM, extra_used={'scale_mm_per_A': scale})


if __name__ == '__main__':
    main()
